
# -*- coding:utf-8 -*-
import re
import numpy as np
import matplotlib.pyplot as plt
import time
import tools
import warnings


total_bead=544
system_all=["WT","M1","M2"]
temps=[278,298]

#slab参数
chunk_gap=10

totalframe=50

#约化系数
smooth_points=600

#---------画图参数
widbar = 0.01
system_color=["tab:grey","#ED797B","#21A6EC"]
hfont = {'fontname':'Arial'}
plt.rcParams['axes.unicode_minus']=False

np.set_printoptions(suppress=True)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.rcParams['ps.fonttype'] = 42
plt.rcParams['pdf.fonttype'] = 42
#-----数组
how_many_sys=len(system_all)*len(temps)

length_for_average= np.zeros((totalframe ,4))#存储高低密度相的长度、粒子
density_for_average= np.zeros((2 ,totalframe))#存储高低密度相的长度、粒子

###########################开始循环###########################

def draw_density(axdensity,system,temperature,color_line,datay):
    #print(datay[0])
    datax=np.arange(0, len(datay[0]))
    axdensity.plot(datax/10,datay[0],label='%s %d'%(system,temperature),color= system_color[color_line])
    plt.xlabel("distance along z-axis",fontdict={ 'size' : 15},**hfont)
    plt.ylabel(r"Reduced density $\rho ^* $",fontdict={ 'size' : 15},**hfont)
    plt.title("T=%d"%temperature,fontdict={ 'size' : 15},**hfont)
    plt.yticks(fontsize=12,**hfont)
    plt.xticks(fontsize=12,**hfont)
    plt.legend()
    #plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.15)
    
def draw_densitybox(axdensity,system,temperature,color_line,datay):
    #print(datay[0])
    datax=np.arange(0, len(datay[0]))
    axdensity.plot(datax[1:],datay[0][1:],label='%s %d'%(system,temperature),color= tools.color_dict[color_line])
    plt.xlabel("density",fontdict={ 'size' : 15},**hfont)
    plt.ylabel(r"PDF",fontdict={ 'size' : 15},**hfont)
    plt.title("T=%d"%temperature,fontdict={ 'size' : 15},**hfont)
    plt.yticks(fontsize=12,**hfont)
    plt.xticks(fontsize=12,**hfont)
    plt.legend()
    #plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.15)

    
def smooth_density(data,chunk):
    data_smoothed = np.zeros((2 ,chunk))
    for c in range (0,chunk_gap):
        data_smoothed[1][c]=0
        for p in range(0,2*chunk_gap):
            tmp_i=c-chunk_gap+p
            if tmp_i <0 :
                tmp_i=chunk+tmp_i
            data_smoothed[1][c]+= data[1][tmp_i]                    
    for c in range (chunk_gap,chunk-chunk_gap):
        data_smoothed[1][c]=0
        for p in range(0,2*chunk_gap):
            data_smoothed[1][c]+= data[1][c-chunk_gap+p]               
    for c in range (chunk-chunk_gap,chunk):
        data_smoothed[1][c]=0
        for p in range(0,2*chunk_gap):
            tmp_i=c-chunk_gap+p
            if tmp_i >=chunk :
                tmp_i=tmp_i-chunk
            data_smoothed[1][c]+= data[1][tmp_i]
    return data_smoothed

def center_density(data_smoothed,chunk,frame,smoothed_for_average):
    single_max=0
    center_chunk=[0]*totalframe
    for the_chunk in range (0, chunk):
        if data_smoothed[1][the_chunk]>single_max:
            single_max=data_smoothed[1][the_chunk]
            center_chunk[frame]=the_chunk
    
    for c in range (0,chunk):#进行居中操作
        new_chunk=c-center_chunk[frame]+int(0.5*chunk)
        if new_chunk<0:
            new_chunk+=chunk
        if new_chunk>chunk-1:
            new_chunk=new_chunk%chunk
            #print(frame,new_chunk)
        smoothed_for_average[frame][new_chunk]=data_smoothed[1][c]
        
def read_chunkframe(chunk,f,frame):
    data= np.zeros((2 ,chunk))
    data[0]=range(0, chunk)
    line = f.readline() #每一个frame的第一行，用于check和统计 
    for the_chunk in range (0, chunk):
        match= re.findall(r'(-?[0-9\.]+)', f.readline())
        data[1][int (match[0])-1]=int(match[1])
    return data

def attribute_density(density_box,data,slot):
    for den in data[1]:
        the_chunk=int(den/slot)
        if the_chunk<len(density_box[0]):
            density_box[0][the_chunk]+=1
        else:
            density_box[0][len(density_box[0])-1]+=1
  

def density_cube(f,chunk,slotsize,totalslots):
    smoothed_for_average = np.zeros((totalframe+1,chunk)) 
    density_box=np.zeros((1,totalslots))
    for i in range(0,3):
        f.readline()
    for frame in range(0, totalframe ):            
        data=read_chunkframe(chunk,f,frame)
        attribute_density(density_box,data,slotsize)
    return density_box
       
def density_distrib(chunk):
    thermo=plt.figure("density",figsize=(8, 6))
    start=time.clock()
    smoothed_for_average = np.zeros((totalframe+1,chunk))
    ti=0
    for temperature in temps:
        color_line=0
        axdensity=thermo.add_subplot(211+ti)
        thermo.subplots_adjust(left=None, bottom=None,wspace=None, hspace=0.4) 
        for system in system_all:
            print(system)
            smoothed_for_average = np.zeros((totalframe+1,chunk))
            f = open("./%s_%d_trimer/0/hsf.density"%(system,temperature),"r")
            for i in range(0,3):
                f.readline()
            for frame in range(0, totalframe ):            
                data=read_chunkframe(chunk,f,frame)
             
                data_smoothed=smooth_density(data,chunk)#近邻平均 ，存储在smooth  
                center_density(data_smoothed,chunk,frame,smoothed_for_average)
                                  
            for c in range (0,chunk):#平均   
                for frame in range(0, totalframe): 
                    smoothed_for_average[totalframe][c]+=smoothed_for_average[frame][c]
                smoothed_for_average[totalframe][c]/=totalframe
                
            draw_density(axdensity,system,temperature,color_line,smoothed_for_average)
            color_line+=1
        ti+=1
    print("-------------------------------finished  density Time used:%.1f"%(time.clock() - start))
    plt.savefig("density.pdf", format="pdf")
    
    return plt
#density_distrib()
