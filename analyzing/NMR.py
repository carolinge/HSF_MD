# -*- coding:utf-8 -*-
import re
import os
import numpy as np
from numba import jit
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time
from tools import read_file


temps=[278,298]
total_chain=40
bead_per_chain=181
total_atom=total_chain*bead_per_chain

#---------画图参数
widbar = 0.01
sns.set()
plt.rcParams['font.family']='SimHei'
plt.rcParams['axes.unicode_minus']=False
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['pdf.fonttype'] = 42

#--------初始化
system_all=["WT","M1","M2"]
system_color=["tab:blue","#5CB8B3","#D78189"]
colormap=["YlORBr","#5CB8B3","Blues"]
input_name="核磁整理数据_0822.xlsx"
dataNMR = pd.read_excel(input_name, skiprows=0, usecols="A:L",
                           index_col=None, sheet_name = "Sheet1", dtype = {"ID":str,"InStore":str})
coordinates = np.zeros((total_atom+1 ,5))
coordinates_dense = np.zeros((total_atom+1 ,5))
coordinates_dilute = np.zeros((total_atom+1 ,5))


NMRhist = np.zeros((3 ,bead_per_chain))
#########################################函数区###########################################################

                
@jit(nopython=True)
def find_neighbor(coordinates,box_len,r_cut,inter_matrix,intra_matrix,inter_intra_hist,frame):  
    count_i=0
    NMRchain = np.zeros((bead_per_chain,total_chain))
    for i in range (0, total_atom):
        for j in range (0, total_atom):  
            if i!=j :
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
                    addvalue=49/mod#认为r=6.5贴上去是1 
                    if chain_A!=chain_B:
                        inter_matrix[seq_A][seq_B]+=addvalue
                        inter_intra_hist[0][seq_A]+=addvalue
                    else:
                        intra_matrix[seq_A][seq_B]+=addvalue
                        inter_intra_hist[1][seq_A]+=addvalue                  
    #for seq in range (0, bead_per_chain):
        #thehist[fi][seq]+=np.mean(NMRchain[seq])
        
        
                        
@jit(nopython=True)     
def reduce(matrix1,reduce_by):
    for i in range (0, bead_per_chain):
        for j in range (0, bead_per_chain):
            if matrix1[i][j]>0:
                matrix1[i][j]/=reduce_by/5
            else:
                matrix1[i][j]=0.0

@jit(nopython=True)     
def reduce_method(matrix1,reduce_method):
    the_max=0.0
    total_contact=0.0
    for i in range (0, bead_per_chain):
        for j in range (0, bead_per_chain):
            total_contact+=matrix1[i][j]
            if reduce_method=="log":
                matrix1[i][j]=np.log(matrix1[i][j]+1.0)
                if  not (matrix1[i][j]>-0.1 and matrix1[i][j]<1000000):
                    matrix1[i][j]=0.0
            elif reduce_method=="root3":
                matrix1[i][j]=math.pow(matrix1[i][j], 0.3)
            if the_max<matrix1[i][j]:
                the_max=matrix1[i][j]
            pass                                
            
    return the_max,total_contact

############################################################################################################
#main program

def get_contacts(axnow,system, temperature,filename='',fromfile=1,frame=50,r_cut=8.5):
    #ratio.suptitle('dense/dilute')
    #ratio.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None) 
    read_frame=frame-20
    skip_frame=20
    box_len=np.zeros((frame ,3))
    reduce_method="notyet"
    if not os.path.exists("processfile"):
        os.mkdir("processfile")
    

    inter_matrix = np.zeros((bead_per_chain ,bead_per_chain))
    intra_matrix = np.zeros((bead_per_chain ,bead_per_chain))
    inter_intra_hist=np.zeros((2 ,bead_per_chain))
    generate=1
    if fromfile:
        generate=0
        try:
            inter_matrix=np.loadtxt("%s/processfile/intermatrix_%s_%s_c%.1f_%dframe.txt"%(filename,system,temperature,r_cut,read_frame),delimiter=' ')  
            intra_matrix=np.loadtxt("%s/processfile/intramatrix_%s_%s_c%.1f_%dframe.txt"%(filename,system,temperature,r_cut,read_frame),delimiter=' ')  
            inter_intra_hist=np.loadtxt("%s/processfile/interahist_%s_%s_c%.1f_%dframe.txt"%(filename,system,temperature,r_cut,read_frame),delimiter=' ')  

        except:
            print("could not find processed file on", os.getcwd())
            generate=1
    if generate:        
        f_dense = open("%s/%s_%d_dense/0/all.dump"%(filename,system,temperature),"r")
        for frame in range(0, read_frame+skip_frame):
            read_file(coordinates_dense,box_len, f_dense, frame)
            if frame>skip_frame:
                find_neighbor(coordinates_dense,box_len,r_cut,inter_matrix,intra_matrix,inter_intra_hist,frame)
                pass
        f_dense.close()

        np.savetxt("%s/processfile/intermatrix_%s_%s_c%.1f_%dframe.txt"%(filename,system,temperature,r_cut,read_frame),inter_matrix,fmt='%.3f')  
        np.savetxt("%s/processfile/intramatrix_%s_%s_c%.1f_%dframe.txt"%(filename,system,temperature,r_cut,read_frame),intra_matrix,fmt='%.3f')                 
        np.savetxt("%s/processfile/interahist_%s_%s_c%.1f_%dframe.txt"%(filename,system,temperature,r_cut,read_frame),inter_intra_hist,fmt='%.3f')                
        
                
    reduce(inter_matrix,float(read_frame))


    #Contact map
    vmin=0.1;vmax=0.5
    mappable=axnow.imshow(inter_matrix, 
       cmap="Spectral_r",   #Spectral_r https://matplotlib.org/stable/tutorials/colors/colormaps.html
       origin='lower',   # orign参数指定绘制热图时的方向，默认值为upper,  此时热图的右上角为(0, 0), 
       vmin=vmin, vmax=vmax  ,  # vmin和vmax参数用于限定数值的范围，只将vmin和vmax之间的值进行映射
       alpha=0.9
       )
    axnow.tick_params(axis='both', direction='in')#, length=6, width=1)  # 可以根据需要调整长度和宽度

    axnow.set_xticks(np.arange(0,180,10))
    axnow.set_xticklabels(np.arange(203,383,10))
    for i, tick in enumerate(axnow.xaxis.get_ticklabels()):
        if i % 5 != 1:
            tick.set_visible(False)
    axnow.set_yticks(np.arange(0,180,10))
    axnow.set_yticklabels(np.arange(203,383,10))
    for i, tick in enumerate(axnow.yaxis.get_ticklabels()):
        if i % 5 != 1:
            tick.set_visible(False)   
    axnow.grid(color='w', linestyle='-', linewidth=0)
    axnow.set_title('%s T=%dK'%(system,temperature))                                                                 
    if  1:
        plt.figure("contact 2D")
        plt.colorbar(mappable)
        #axnow.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)
        1
  
    return inter_matrix
    print("-------------------------------finished NMR")

    #plt.show()
    print("-------------------------------finished NMR")

    #plt.show()


