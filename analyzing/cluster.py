
# -*- coding:utf-8 -*-
import re
import numpy as np
import matplotlib.pyplot as plt
import time
import warnings
#warnings.filterwarnings("ignore")

start = time.clock()
###########################变量设置（建议检查每行）###########################

total_bead=544
system_all=["WT","M1","M2"]
temps=[278,298]

#slab参数
chunk_gap=10
chunk =300
totalframe=300
#cluster 参数
ensemble=1
frame=700

#约化系数
smooth_points=600

#---------画图参数
widbar = 0.01
system_color=["#F8991C","#ED797B","#21A6EC"]
hfont = {'fontname':'Arial'}

np.set_printoptions(suppress=True)


#-----数组
how_many_sys=len(system_all)*len(temps)
max_lines=( (chunk+1)*totalframe+3)
critical_temp=np.zeros((2 ,how_many_sys))
pattern_1= r'([0-9]+\s+-?[0-9]*\.?[0-9]+)'
data = np.zeros((2 ,chunk))
data[0]=range(0, chunk)
data_smoothed = np.zeros((2 ,chunk))
length_for_average= np.zeros((totalframe ,4))#存储高低密度相的长度、粒子
density_for_average= np.zeros((2 ,totalframe))#存储高低密度相的长度、粒子
smoothed_phase=np.zeros((2*how_many_sys ,2*smooth_points))
density_distribution=np.zeros((how_many_sys ,100))


###########################函数区###########################

def read_file(largest_data,second_largest_data,f,e):
    for wt in range(3):
        line = next(f)
    the_frame=0
    while True:
        line = f.readline()
        if line:
            match= re.findall(r'(-?[0-9\.]+)', line)
            if(len(match)==2):
                the_frame += 1
                maxsize=0
                secondmax=0
            else:
                if int(match[1])>maxsize:
                    secondmax=maxsize
                    maxsize=int(match[1])
                    largest_data[e][the_frame]=maxsize
                    second_largest_data[e][the_frame]=secondmax
                elif int(match[1])>secondmax:
                    secondmax=int(match[1])
                    second_largest_data[e][the_frame]=secondmax
        else:
            break
    return the_frame

###########################开始循环###########################


def get_max_cluster(axnow,temperature,frame=50,r_cut=8.5,filename=''):
    sn=0
    color_line=0
    for system in system_all:
        print(system)
        sn+=1
        largest_data = np.zeros((ensemble,frame))
        second_largest_data = np.zeros((ensemble,frame))
        averaged_data = np.zeros((how_many_sys+1,frame))
        averaged_data[0]=range(1, frame+1)
        for e in range(0,1):
            f=open("%s/%s_%d_trimer/%d/hsf.cluster"%(filename,system,temperature,e),"r")
            totalframe=read_file(largest_data,second_largest_data,f,e)
            for fr in range(0,totalframe):
                if largest_data[0][fr]>0:                
                    averaged_data[sn][fr]+=largest_data[e][fr]
                    averaged_data[sn][fr]+=second_largest_data[e][fr]
            for fr in range(0,totalframe):
                averaged_data[sn][fr]/=ensemble#*2
            axnow.plot(averaged_data[0][0:totalframe],averaged_data[sn][0:totalframe]/544,linewidth=2,label='%s %d'%(system,temperature),color= system_color[color_line])
            axnow.set_xlabel("time(*10^5 step)",fontdict={ 'size' : 15},**hfont)
            axnow.set_ylabel(r"Protein # in largest droplet",fontdict={ 'size' : 15},**hfont)
            axnow.set_title("T=%d"%temperature,fontdict={ 'size' : 15},**hfont)
            #axnow.set_yticks(fontsize=12,**hfont)
            #axnow.set_xticks(fontsize=12,**hfont)
            axnow.set_xlim([0,totalframe])
            axnow.legend()
    ###########################还是得靠拟合啊。。。###########################
        color_line+=1

    


#https://matplotlib.org/3.1.0/gallery/color/named_colors.html配色
