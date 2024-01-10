# -*- coding:utf-8 -*-
import re
import numpy as np
from numba import jit
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time
import math
import os
from tools import read_file


#---------画图参数
widbar = 0.01
sns.set()
plt.rcParams['font.family']='SimHei'
plt.rcParams['axes.unicode_minus']=False

#--------初始化
temps=[278,298]
total_chain=40
bead_per_chain=181
total_atom=total_chain*bead_per_chain
system_all=["WT","M1","M2"]
system_color=["tab:blue","#5CB8B3","#D78189"]
coordinates = np.zeros((total_atom+1 ,5))
coordinates_dense = np.zeros((total_atom+1 ,5))

dense_matrix = np.zeros((bead_per_chain ,bead_per_chain))
dense_3D= np.zeros((6,bead_per_chain ,bead_per_chain))
final_data=pd.DataFrame(np.zeros((bead_per_chain ,3)), index=range(0,bead_per_chain), columns=system_all)
#########################################函数区###########################################################

           
@jit(nopython=True)
def find_neighbor(coordinates,box_len,r_cut,thehist ,inter_matrix,frame):  
    count_i=0
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
                mod=mod**0.5
                if mod < r_cut:
                    if chain_A!=chain_B:
                        inter_matrix[seq_A][seq_B]+=1.0

@jit(nopython=True)     
def preduce(matrix,read_frame):
    for i in range (0, bead_per_chain):
        for j in range (0, bead_per_chain):
            matrix[i][j]/=read_frame            

@jit(nopython=True)     
def reduce(matrix1,reduce_method):
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

def contact_map(hist2D,fromfile=1,frame=50,r_cut=8.5,reduce_method="root3",vvmin=0.0,vvmax=0.5):    
    read_frame=frame-5
    skip_frame=5
    box_len=np.zeros((frame ,3))
    #hist2D=plt.figure("contact 2D",figsize=(10, 5.8))
    if not os.path.exists("processfile"):# os.makedirs 传入一个path路径，生成一个递归的文件夹；如果文件夹存在，就会报错,因此创建文件夹之前，需要使用os.path.exists(path)函数判断文件夹是否存在；
        os.makedirs("processfile")  # makedirs 创建文件时如果路径不存在会创建这个路径

    row=0
    countt=0
    for temperature in temps:
        sn=0
        for system in system_all:
            print(system,temperature,"rcut=", r_cut,reduce_method)
            NMRhist_dense = np.zeros((2 ,bead_per_chain))
            dense_matrix = np.zeros((bead_per_chain ,bead_per_chain))
            if fromfile:
                dense_matrix=np.loadtxt("processfile/dense_%s_%d_matrix_c%.1f_%dframe.txt"%(system,temperature,r_cut,read_frame),delimiter=' ')
            else:
                filename_dense="./%s_%d_dense/0/all.dump"%(system,temperature)
                f_dense = open(filename_dense,"r")
                for frame in range(0, read_frame+skip_frame):
                    read_file(coordinates_dense,box_len, f_dense, frame)
                    if frame>skip_frame:
                        find_neighbor(coordinates_dense,box_len,r_cut,NMRhist_dense ,dense_matrix,frame)
                print("total contact",preduce(dense_matrix,read_frame))
                f_dense.close()
                np.savetxt("processfile/dense_%s_%d_matrix_c%.1f_%dframe.txt"%(system,temperature,r_cut,read_frame),dense_matrix,fmt='%.3f')
            reduced_by,total_contact=reduce(dense_matrix,reduce_method)
            print("total before:%.2f after: %.2f |max: %.2f "%(total_contact,np.sum(dense_matrix),reduced_by))
        #------------------作图 
            ax2D=hist2D.add_subplot(231+sn+row*3)
            #ax=sns.heatmap(dense_matrix,vmin=vvmin, vmax=vvmax, cmap=sns.color_palette("Spectral_r", as_cmap=True),yticklabels=yticklabels)#vmax
            ax2D.imshow(dense_matrix, 
               cmap="Spectral_r",   # https://matplotlib.org/stable/tutorials/colors/colormaps.html
               origin='lower',   # orign参数指定绘制热图时的方向，默认值为upper,  此时热图的右上角为(0, 0), 
               vmin=0, vmax=0.5    # vmin和vmax参数用于限定数值的范围，只将vmin和vmax之间的值进行映射
               )
            ax2D.set_xticks([30,60,90,120,150])
            ax2D.set_yticks([30,60,90,120,150])
            ax2D.grid(color='w', linestyle='-', linewidth=0)
            dense_3D[countt]=inter/reduced_by
            ax2D.set_title('%s T=%dK'%(system,temperature))                                                                 
            sn+=1
            countt+=1
        row+=1
    plt.subplots_adjust(left=None, bottom=None,wspace=0.5, hspace=0.6) 
    plt.suptitle('Contact map, Dated 20%s.%s.%s,%s:%s \n rcut=%.1f, method=%s,vmax=%.2f, vmin=%.2f'%(
                            time.strftime("%y"),time.strftime("%m"),time.strftime("%d")
                            ,time.strftime("%H"),time.strftime("%M"),r_cut,reduce_method, vvmax, vvmin))    

    plt.savefig('contactmap_%s%s%s_%s%s.svg'%(time.strftime("%y"),time.strftime("%m"),time.strftime("%d"),time.strftime("%H"),time.strftime("%M")))
    print("-------------------------------finished contact")
    return plt
#contact_map()
