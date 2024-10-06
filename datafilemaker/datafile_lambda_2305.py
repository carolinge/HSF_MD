import numpy as np
import pandas as pd
import re
import os
import time

bead_per_chain=181
temps=[298,278]
sys=["WT","M1","M2"]
files=["single","dense","trimer","dilute"]
def setsystem(con):
    if con=='trimer':
        total_chain=50
        skipline=81251#50: 81251,  20: 32531 ||
    elif con=='single':
        total_chain=1
        skipline=591
    else:
        total_chain=40#40
        skipline=21651#mono 5:L2751  20: 10851 40 21651
    return total_chain, skipline

def copy_oldfile(f_ori,f_new,skipline):
    for linec in range(0, skipline):
        line=f_ori.readline()                  
        f_new.write(line)
    print("Lambdas\n",file=f_new)
    f_ori.close()
    f_new.close()        
    
amino_one=["M","E","S","N","Q","G","A","I","R","Y",
       "V","P","H","D","F","T","W","K","C","L","Ca"]
HPS=pd.DataFrame( data=data_hps, index=amino_one, columns=parameters)

    
    
#-------------------------------------预设区请勿修改------------------------------------------

def write_datafile(system,temp,lambdas,newdir):
    for con in files:
        total_chain, skipline=setsystem(con)
        output_name="%s/file/%s_%d_%s.data"%(newdir,system,temp,con)
        f_ori = open("annealing/datafilemaker/%s/%s.data"%(con,system),'r')
        f_new = open(output_name,'w')    

    #------原来文件输出---
        copy_oldfile(f_ori,f_new,skipline)
    #--------------------
        with open (output_name, "a") as fout:
            atom_number=1
            total_interac=0.0
            for i in range (0,total_chain):
                if con=='trimer':
                    print("%d 0.0"%atom_number,file=fout)
                    atom_number=atom_number+1
                    for mono in range(0,3):
                        for bead in range(0,bead_per_chain):
                            res_lambda=lambdas[bead]
                            print("%d %.3f"%(atom_number,res_lambda),file=fout) 
                            atom_number=atom_number+1
                            
                else:
                    for bead in range(0,bead_per_chain):
                        res_lambda=lambdas[bead]
                        print("%d %.3f"%(atom_number,res_lambda),file=fout)
                        atom_number=atom_number+1
    np.savetxt("%s/%s_%d_lambdas.data"%(newdir,system,temp),lambdas)
    
    

def new_datafile(system,temp,std_lambda,newdir):
    dataNMR = pd.read_csv("annealing/datafilemaker/NMRinter.csv")
    final_lambdas=np.zeros(181) 
    for con in files:
        total_chain, skipline=setsystem(con)
        output_name="%s/file/%s_%d_%s.data"%(newdir,system,temp,con)
        f_ori = open("annealing/datafilemaker/%s/%s.data"%(con,system),'r')
        f_new = open(output_name,'w')    
        

    #------原来文件输出---
        copy_oldfile(f_ori,f_new,skipline)
    #--------------------
        with open (output_name, "a") as fout:
            atom_number=1
            total_interac=0.0
            for i in range (0,total_chain):
                if con=='trimer':
                    print("%d 0.0"%atom_number,file=fout)
                    atom_number=atom_number+1
                    for mono in range(0,3):
                        for bead in range(0,bead_per_chain):
                            restype=dataNMR[system][bead]
                            res_lambda=HPS[std_lambda][restype]
                            print("%d %.3f"%(atom_number,res_lambda),file=fout) 
                            atom_number=atom_number+1
                            
                else:
                    for bead in range(0,bead_per_chain):
                        restype=dataNMR[system][bead]
                        res_lambda=HPS[std_lambda][restype]
                        final_lambdas[bead]=res_lambda
                        print("%d %.3f"%(atom_number,res_lambda),file=fout)
                        atom_number=atom_number+1
                    
                        
    np.savetxt("%s/%s_%d_lambdas.data"%(newdir,system,temp),final_lambdas)
    return final_lambdas
    
