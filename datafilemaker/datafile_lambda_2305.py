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
parameters=["radius","charge","KR","urry","Larsen","FB","doolittle","Carter","capacity"]
data_hps= np.zeros((len(amino_one) ,len(parameters)))
data_hps[0]=[6.18,  0.0,    0.838,  0.676,  0.7256,  0.64648,  0.711, 0.000, 71.10]#MET  M
data_hps[1]=[5.92,  -1.0,    0.459,  0.0,  0.2217,  0.42621,  0.111,   0.000, -8.75]#GLU   E 酸
data_hps[2]=[5.18,  0.0,    0.595,  0.588,  0.3724,  0.11195,  0.411,    0.000, 10.43]#SER  S
data_hps[3]=[5.58,  0.0,    0.432,  0.588,  0.2648,   0.78447,  0.111,   0.000, 6.18
]#ASN  N 酰胺
data_hps[4]=[6.02,  0.0,    0.514,  0.559,  0.2546,  0.29516,  0.111,   0.000,  51.47]#GLN  Q 酰胺
data_hps[5]=[4.50,  0.0,    0.649,  0.574,  0.4573,  1.24153,  0.456,    0.000, 44.73]#GLY  G
data_hps[6]=[5.04,  0.0,    0.730,  0.603,  0.5332,  0.51507,  0.7, 0.000,  48.93]#ALA  A
data_hps[7]=[6.18,  0.0,    0.973,  0.706,  0.6717,  0.83907,  1.0, 0.000, 105.23
]#ILE I
data_hps[8]=[6.56,  1.0,    0.000,  0.588,  0.1583,  0.24025,  0.0,    0.000, 69.39
]#ARG  R
data_hps[9]=[6.46,  0.0,    0.865,  0.897,  0.6031,   1.04266,  0.356,   0.000, 71.09]#TYR  Y
data_hps[10]=[5.86,  0.0,    0.892,  0.664,  0.7797,  0.55645,  0.967,    0.000, 105.80]#VAL V 
data_hps[11]=[5.56,  0.0,    1.000,  0.759,  0.4372,   0.34128,  0.322,  0.000, 105.80]#PRO P 
data_hps[12]=[6.08,  0.0,    0.514,  0.765,  0.4902,   0.55537,  -3.2,  0.000, 38.01]#HIS H 
data_hps[13]=[5.58,  -1.0,    0.378,  0.294,  0.2075,  0.30525,  0.111,  0.000, -44.97]#ASP D 酸
data_hps[14]=[6.36,  0.0,    1.000,  0.824,  0.000,  1.17854,  0.811, 0.000, 102.24]#PHE F
data_hps[15]=[5.62,  0.0,    0.676,  0.588,  0.4255,  0.27538,  0.422,   0.000, 50.06]#THR T
data_hps[16]=[6.78,  0.0,    0.946,  1.000,  0.8206,  0.97588,  -0.4,   0.000,  108.10]#TRP W
data_hps[17]=[6.36,  1.0,    0.514,  0.382, 0.1851,   0.47106,  0.067,   0.000, 29.98]#LYS K
data_hps[18]=[5.48,  0.0,    0.595,  0.676,  0.7691,  0.46169,  0.778,    0.000, 43.61
]#CYS C
data_hps[19]=[6.18,  0.0,    0.973,  0.721,  0.8028,   0.51207,  0.922,   0.000,  109.38
]#LEU L
data_hps[20]=[6.00,  1.0,    0.00,  0.00,  0.00,   0.00,  0.00,   0.000, 0.000]#cation
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
    
