# -*- coding:utf-8 -*-
import re
import os
import numpy as np
from numba import jit
import pandas as pd
import time
from tools import read_file,read_file2





total_chain=40
bead_per_chain=181
total_atom=total_chain*bead_per_chain
reduce_method="notyet"
#---------画图参数
widbar = 0.01

#--------初始化
system_color=["tab:blue","#5CB8B3","#D78189"]
dataNMR = pd.read_csv("NMRinter.csv")


coordinates_dense = np.zeros((total_atom+1 ,5))
coordinates_dilute = np.zeros((total_atom+1 ,5))

#########################################函数区###########################################################


@jit(nopython=True)
def cal_contact(coordinates,box_len,r_cut,inter_matrix):
    for i in range (0, total_atom):
        for j in range (0, total_atom):  
            if i!=j :#and abs(i-j)>2:
                #get index
                chain_A=int(i/bead_per_chain)
                chain_B=int(j/bead_per_chain)     
                seq_A=int(i%bead_per_chain)
                seq_B=int(j%bead_per_chain)
                type_A=coordinates[i][0]
                type_B=coordinates[j][0]
                mod=0
                #get distance
                for k in range (0,3):
                    dist= abs(coordinates[i][k+1]-coordinates[j][k+1])     
                    if dist>0.5*box_len[k]:
                        dist=box_len[k]-dist                         
                    mod+=(dist*dist)
                mod=mod**0.5
                if mod < r_cut:
                    if chain_A!=chain_B:
                        inter_matrix[seq_A][seq_B]+=1.0#你似乎没有规定A和b谁大归谁
    return inter_matrix
                
@jit(nopython=True)
def find_neighbor(coordinates,box_len,r_cut,thehist ,inter_matrix,fi,frame):  
    count_i=0
    NMRchain = np.zeros((bead_per_chain,total_chain))
    for i in range (0, total_atom):
        for j in range (0, total_atom):  
            if i!=j :#and abs(i-j)>2:
                #get index
                chain_A=int(i/bead_per_chain)
                chain_B=int(j/bead_per_chain)     
                seq_A=int(i%bead_per_chain)
                seq_B=int(j%bead_per_chain)
                type_A=coordinates[i][0]
                type_B=coordinates[j][0]
                mod=0
                #get distance
                for k in range (0,3):
                    dist= abs(coordinates[i][k+1]-coordinates[j][k+1])     
                    if dist>0.5*box_len[frame][k]:
                        dist=box_len[frame][k]-dist                         
                    mod+=(dist*dist)
                #mod=mod**0.5
                if mod < r_cut*r_cut:
                    if chain_A!=chain_B:
                        inter_matrix[seq_A][seq_B]+=1.0
                        #thehist[1][seq_A]+=1
                        NMRchain[seq_A][chain_A]+=49/mod#认为r=6.5贴上去是1    
    for seq in range (0, bead_per_chain):
        '''for chain in range (0, total_chain):
            if NMRchain[seq][chain]>2:
                NMRchain[seq][chain]=1
            else:
                NMRchain[seq][chain]=0'''
        thehist[fi][seq]+=np.mean(NMRchain[seq])
        
        
                        
@jit(nopython=True)     
def reduce(matrix1,reduce_by):
    for i in range (0, bead_per_chain):
        for j in range (0, bead_per_chain):
            if matrix1[i][j]>0 and matrix1[i][j]<1:
                matrix1[i][j]/=reduce_by
            else:
                matrix1[i][j]=0.0
            

def normalize_NMR(NMRhist):
    summ=np.sum(NMRhist)
    NMRhist=1-NMRhist
    minn=np.max(NMRhist);maxx=np.min(NMRhist)
    #print(minn)
    NMRhist/=minn
    #for i in range (0, len(NMRhist)):
    
        #NMRhist[i]=NMRhist[i]/maxx
        #NMRhist[i]=NMRhist[i]/minn
        #NMRhist[i]=1-(NMRhist[i])*3
    print("max contact:",maxx,"min contact:",minn,"total NMR:",summ)
    #print(NMRhist)
    return NMRhist

      
def get_lambda_excel(temp,system,bead):
    NMR_int=1.0
    for trail in range(0,bead_per_chain):
        if dataNMR['Res(%s)'%(system)][trail]==bead+1:
            NMR_int=dataNMR['%dK_%s'%(temp, system)][trail]
            if NMR_int>1.0:
                NMR_int=1.0
            elif NMR_int<1.0 and NMR_int>0.0:
                pass
            else:
                NMR_int=1.0
            #print("excel",trail+2,"|""beadid",bead+1, system,NMR_int)
    return 1.0-NMR_int

def get_hsqc_exp(temperature,system):
    hsqc = np.zeros(bead_per_chain)
    for bead in range(bead_per_chain):
        hsqc[bead]=1-get_lambda_excel(temperature,system,bead)
    return hsqc

def get_hsqc_mask(temp,system):
    hsqc = np.zeros(bead_per_chain)
    mask = np.zeros(bead_per_chain)
    for bead in range(0,bead_per_chain):
        if  dataNvoid['Res(%s)'%(system)][bead]==bead+1:
            NMR_int=dataNvoid['%dK_%s'%(temp, system)][bead]
            if NMR_int>1.0:
                hsqc[bead]=1.0
            elif NMR_int<1.0 and NMR_int>0.0:
                hsqc[bead]=NMR_int
            else:
                #print("misssing site",bead)
                mask[bead]=1
        else:
            print("wrong at excel",bead,"|""beadid",bead+1, system)
    return hsqc, mask
    
        
############################################################################################################
#main program
def calc_contacts_site(read_frame=10,skip=5,r_cut=8.5,sitenum=1,system="no"):
    
    if not os.path.exists("processfile/sitecontact_%s_%.1fcut_%dframe"%(system,r_cut,read_frame)):
        os.mkdir("processfile/sitecontact_%s_%.1fcut_%dframe"%(system,r_cut,read_frame))
        
    contact_list_site=np.zeros((read_frame,181)) 
    dense_matrix = np.zeros((bead_per_chain ,bead_per_chain))
    if 0:#fromfile:
        contact_list_site=np.loadtxt("processfile/sitecontact_%s_%.1fcut_%dframe/%d.txt"%(system,r_cut,read_frame,sitenum),delimiter=' ')
    else:
        f_dense = open("./%s_dense/0/all.dump"%(system),"r")
        for frame in range(0, read_frame+skip):
            box_len=read_file2(coordinates_dense, f_dense, frame)
            if frame>=skip:
                contact_all=cal_contact(coordinates_dense,box_len,r_cut,dense_matrix)
                contact_list_site[frame-skip]=contact_all[sitenum]
        f_dense.close()
        np.savetxt("processfile/sitecontact_%s_%.1fcut_%dframe/%d.txt"%(system,r_cut,read_frame,sitenum),contact_list_site,fmt='%.3f')
    
    return contact_list_site
        
      

def get_hsqc_frame(fromfile=1,read_frame=50,skip=5,r_cut=8.5,system="no"):
    
    
    if not os.path.exists("processfile"):
        os.mkdir("processfile")
    box_len=np.zeros((read_frame+skip ,3))
    print(system,read_frame)
    NMRhist_dilute = np.zeros((read_frame ,bead_per_chain))
    NMRhist_dense = np.zeros((read_frame ,bead_per_chain))
    NMRhist = np.zeros((read_frame,bead_per_chain))
    dense_matrix = np.zeros((bead_per_chain ,bead_per_chain))
    dilute_matrix = np.zeros((bead_per_chain ,bead_per_chain))
    
    generate=1
    if fromfile:
        generate=0
        try:
            #NMRhist_dilute=np.loadtxt("processfile/histdilute_%s_matrix_c%.1f_%dframe.txt"%(system,r_cut,read_frame),delimiter=' ')
            NMRhist_dense=np.loadtxt("processfile/histdense_%s_matrix_c%.1f_%dframe.txt"%(system,r_cut,read_frame),delimiter=' ')      
        except:
            generate=1
    if generate:
        f_dense = open("./%s_dense/0/all.dump"%(system),"r")
        #f_dilute = open("./%s_dilute/0/all.dump"%(system),"r")
        fi=0
        for frame in range(0, read_frame+skip):
            read_file(coordinates_dense,box_len, f_dense, frame)
            #read_file(coordinates_dilute,box_len, f_dilute, frame)
            if frame>=skip:
                find_neighbor(coordinates_dense,box_len,r_cut,NMRhist_dense,dense_matrix,fi,frame)
                #find_neighbor(coordinates_dilute,box_len,r_cut,NMRhist_dilute,dilute_matrix,fi,frame)
                fi+=1
        f_dense.close()
        #f_dilute.close()
        #np.savetxt("processfile/histdilute_%s_matrix_c%.1f_%dframe.txt"%(system,r_cut,read_frame),NMRhist_dilute,fmt='%.3f')
        np.savetxt("processfile/histdense_%s_matrix_c%.1f_%dframe.txt"%(system,r_cut,read_frame),NMRhist_dense,fmt='%.3f')            
    diff_matrix=dilute_matrix
    reduce(diff_matrix,float(read_frame))
    #print(NMRhist_dense,NMRhist_dilute)
    ave_densehist=np.average(NMRhist_dense,0)
    ave_dilutehist=np.average(NMRhist_dilute,0)
    #print(ave_densehist,ave_dilutehist)
    
    NMRhist=ave_densehist#/np.sum(ave_dilutehist)

    NMRhist=normalize_NMR(NMRhist)
    '''plt.figure("NMR",figsize=(10,3))
    plt.bar(np.arange(0,181),np.average(NMRhist_dense,0),alpha=0.5)
    plt.bar(np.arange(0,181),np.average(NMRhist_dilute,0),alpha=0.6)
    plt.ylim(0.6,1.0)
    plt.show()'''
    

    #print(np.shape(NMRhist_dense),np.shape(NMRhist_dilute),np.shape(NMRhist))
    return NMRhist


