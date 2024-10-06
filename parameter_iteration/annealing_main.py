import sys
import os
import math
import shutil
import numpy as np
from NMR import get_hsqc_frame,get_hsqc_exp,get_hsqc_mask,calc_contacts_site
from datafilemaker.datafile_lambda_2305 import write_datafile, new_datafile



#numpy.set_printoptions(threshold=numpy.nan)
# File to make the next alphas for a prbias simulation
# Inputs: alpha file from previous simulation (old_alpha_file), histogram of pr from previous
# simulation (f_sim.txt), pr data from experiment (pr_file), temperature (temp),
# tolerance for convergence (tol), factor to get rid of low frequency modes (cutEig),
# number of atoms used in the p(r) (natoms) number of previous iterations (num_iter)

totalframe=380
n=181
system_all=["WT","M1"]
system_all=["M2"]
system_all=["WT","M1","M2"]
temps=[278,288,298,303]

cutEig=176#int(sys.argv[6])


def get_next_lambda(hsqc_sim,hsqc_exp,last_lamda,dev,step):#To be modifed using max entropy
    new_lamda=last_lamda       
    for i in range (n):
        new_lamda[i]=last_lamda[i]+0.3*(hsqc_sim[i]-hsqc_exp[i])

        if new_lamda[i]<0.0:
            print("----WARNING: we did get a negative---")
            new_lamda[i]=0.0
            
    return new_lamda

def annealing(step):
    sumseries=[];nameseries=[];
    if step:
        os.chdir('../step_%d'%(step-1))
        #os.makedirs('../step_%d'%step, exist_ok=True)


    row=0
    for temp in temps:
        sn=0
        for system in system_all:
            print("=-=-=-=-=-=-=-=-=-=-=-=-\n",system, temp)

            
            if step!=0:
                # Get lambda from previous simulation
                last=np.loadtxt("%s_%d_lambdas.data"%(system,temp))
                print("total lambda(last):%.4f"%np.sum(abs(last)))

                # Check and Write in the optimization file
                previous=np.loadtxt("../annealing/lambdarecords_%s_%d.txt"%(system,temp))
                if step==1:
                    previous_errors=[[0],[0]]
                else:
                    previous_errors=np.loadtxt("../annealing/errors_%s_%d.txt"%(system,temp))
                print
                new=np.c_[previous,last]
                np.savetxt("../annealing/lambdarecords_%s_%d.txt"%(system,temp),new,fmt="%.4f")
                NEWT=new.T
                #print(np.shape(new[:][-1]),np.sum(NEWT[-1]),np.sum(NEWT[-2]),np.sum(last))

                # get HSQC from simulation
                hsqc_sim=get_hsqc_frame(fromfile=1,read_frame=totalframe,skip=5,r_cut=8.5,system="%s_%d"%(system,temp))
                hsqc_sim_ave=hsqc_sim#np.average(hsqc_sim,0)
                # get HSQC from experiment
                hsqc_exp=get_hsqc_exp(temp,system)
                
                #print("shape:",np.shape(hsqc_exp),"\n",np.shape(hsqc_sim_ave),"\n")
                # Calculate error
                diff=(hsqc_sim_ave-hsqc_exp)/hsqc_exp
                #np.savetxt('testing_diff.txt',diff)
                if not os.path.exists('../step_%d/file'%(step)):
                    #.mkdir('../step%d'%(step+1))
                    shutil.copytree("../file",'../step_%d'%(step)) 
                    
                # calculate the standard deviation of each element and overall error
                dev=[]
                for i in range(181):
                        dev.append((hsqc_sim_ave[i]-hsqc_exp[i])/hsqc_exp[i])
                print(dev[0:3],np.sqrt(np.mean(np.square(dev))))
                print("total error (stddev): ",np.sqrt(np.mean(np.square(dev))))
                newerror=np.c_[previous_errors,[step,np.sqrt(np.mean(np.square(dev)))]]
                np.savetxt("../annealing/errors_%s_%d.txt"%(system,temp),newerror,fmt="%.4f")   
                    
                    
                #write for the next step
                print("writing datafile....")
                new_lamdas=get_next_lambda(hsqc_sim_ave,hsqc_exp,last,dev,step)
                os.chdir("..")
                write_datafile(system,temp,new_lamdas,'step_%d'%(step))
                os.chdir('step_%d'%(step-1))
            else:
                os.chdir("..")
                if not os.path.exists('step_0/file'):
                    shutil.copytree("file",'step_0') 
                new_lamdas=new_datafile(system,temp,"Larsen",'step_%d'%(step))
                np.savetxt("annealing/lambdarecords_%s_%d.txt"%(system,temp),new_lamdas,fmt="%.4f")
                os.chdir('step_0')            
            print("total lambda(new):%.4f\n\n\n"%np.sum(new_lamdas))
            sumseries.append(np.sum(new_lamdas))
            nameseries.append("%s_%d"%(system,temp))
            sn+=1
        row+=1
    print(sumseries)	
    print(nameseries)

