
# -*- coding:utf-8 -*-
import re
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline as spline
from scipy import interpolate
import math


np.set_printoptions(suppress=True)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
hfont = {'fontname':'Arial'}
plt.rcParams['font.family']='SimHei'
plt.rcParams['axes.unicode_minus']=False
system = [ "WT","M1","M2"]
temp=[278,298]
experiments=[49.8,49.8, 36.1,36.1,56.2,56.2]
epsilon=["0.11","0.115","0.12","0.13"]
epsilon=["0.115"]
how_many_sys=len(system)*len(temp)
color_dict=[ '#14517c', '#2f7fc1', '#c7dfda', '#96c37d', '#f3d266', '#d8383a', '#e7b1bd', '#f1e1e1', '#c497b2', '#a8b8c6']
system_color=['#14518c','#75a5cc',"#5aa8a3","#9CD8d3","#C78189","#F98189"]
#system_color=['#14517c',"#5CB8B3","#C78189",'#55717c',"#8CD8B3","#F98189"]
system_color=['#E48E40','#F3CBA7',"#C83232","#ED797B","#3E7BA4","#21A6EC"]
total_chain=1
bead_per_chain=181
total_atom=total_chain*bead_per_chain
coordinates_dense = np.zeros((total_atom+1 ,5))

gyration_average=np.zeros((4 ,how_many_sys))
#PDF
chunk_gap=10.0
total_chunk=30
PDF=np.zeros((how_many_sys+1,total_chunk))

for i in range (0,total_chunk):
    PDF[0][i]=i*chunk_gap+0.5*chunk_gap
    

def read_file(coordinates,f, frame):
    box_len=np.zeros((1 ,3))
    for i in range(0, 5):
        f.readline()
    for i in range(0, 3):   
        match= re.findall(r'(-?[0-9\.e+]+)', f.readline())
        box_len[0][i]=float (match[1])
    if frame%50==0:
        print("check frame: ",frame, box_len[0][0],box_len[0][1],box_len[0][2],total_atom)
    f.readline()  
    for i in range (0, total_atom):
        line = f.readline()
        match=line .split()
        if float(match[0])<=total_atom:
            coordinates[int(match[0])-1][0]=int(match[1])
            dist=[0, 0, 0]
            for x in range (0,3):                
                coordinates[int(match[0])-1][x+1]=float(match[x+2])*box_len[0][x]
    for i in range (1, total_atom):
        dist=[0, 0, 0]
        for x in range (0,3):
            dist[x]=coordinates[i][x+1]-coordinates[i-1][x+1]
            if  dist[x]>50:
                coordinates[i][x+1]-=box_len[0][x]
            elif dist[x]<-50:
                coordinates[i][x+1]+=box_len[0][x]
            dist[x]=coordinates[i][x+1]-coordinates[i-1][x+1]
        if dist[0]>15 or  dist[1]>15 or dist[2]>15:
            print(i,dist[0], dist[1], dist[2],coordinates[i][1],coordinates[i][2],coordinates[i][3])
            print(aaa)
                        

def radius_of_gyration(coordinates):
    '''
    Extracts data from population list to compute the weighted average of
    the radius of gyration Rg and its error.
    In:  population (list): grown population of polymers
                            Each element in the list is a list [polymer, pol_weight]
    Out: Rg (float): average square of the radius of gyration Rg
         Rg_err (float): error on Rg
    '''
    N=bead_per_chain
    center_of_mass = np.zeros(3)
    sum_rad_gy = 0
    for k in range (0,3):
        for i in range(0, N):
            center_of_mass[k] += coordinates_dense[i][k+1]
            #没有去除奇怪的点
        center_of_mass[k] /= bead_per_chain
    for k in range (0,3):    
        for i in range(0, N):
            diff = (center_of_mass[k] - coordinates_dense[i][k+1])
            sum_rad_gy += diff* diff
    Rg, Rg_err=math.sqrt(sum_rad_gy / N)+3.0 , sum_rad_gy / (2 * N * N)
    #Rg, Rg_err=math.sqrt(sum_rad_gy / N), sum_rad_gy / (2 * N * N)    
    return Rg, Rg_err


def draw_PDF(PDF,se,exper,t,read_frame):
    PDF_smoothed=np.zeros((1,total_chunk))
    for i in range(0,total_chunk):
        PDF[se+1][i]/=read_frame
    for i in range(1,total_chunk-1):
        PDF_smoothed[0][i]=(PDF[se+1][i-1]+ PDF[se+1][i]+PDF[se+1][i+1] )/3
        
    xnew = np.linspace(10,PDF[0].max(),300)
    tck=interpolate.splrep(PDF[0],PDF[se+1])
    y_spline=interpolate.splev(xnew,tck)
    if t==278:
        plt.plot(xnew,y_spline,label='T=%d'%t, color= system_color[se])
        #plt.plot(PDF[0],PDF[se+1],label='T=%d'%t, color= system_color[se])
        plt.fill_between(xnew, 0, y_spline, facecolor= system_color[se], alpha=0.35)
    else:
        plt.plot(xnew,y_spline,label='T=%d'%t, color= system_color[se],alpha=0.75)
        #plt.plot(PDF[0],PDF[se+1],label='T=%d'%t, color= system_color[se],alpha=0.75)
        plt.fill_between(xnew, 0, y_spline, facecolor= system_color[se], alpha=0.3)
    
    plt.xlabel("$R_g (\AA)$",fontdict={ 'size' : 15})
    plt.ylabel("PDF",fontdict={ 'size' : 15})
    plt.xlim(10, 70.0)
    plt.ylim(0, 0.6)
    if t==278:
        plt.axvline(x=gyration_average[0][se],ls="--", lw=1.2, c=system_color[se])#添加竖直直线
        plt.axvline(x=exper, lw=1,label='exp T=278k', c="black")#添加竖直直线
    else:
        plt.axvline(x=gyration_average[0][se],ls="--", lw=1.2, c=system_color[se], alpha=0.9)#添加竖直直线
        #plt.text(gyration_average[0][se]-6, 0.5, r'SAXS',fontsize=10)
    #plt.yticks(np.arange(0.1, 0.8, 0.1),fontsize=12)
    #plt.xticks(np.arange(10.0, 26.0, 3.0),fontsize=12)
    plt.legend(fontsize=10)
    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)

def gyration_single(skip_frame, read_frame):
    se=0
    PDFfig=plt.figure("PDF epsilon", figsize=(5,8))
    gyrations=np.zeros((how_many_sys ,read_frame+skip_frame))
    ti=0
    for j in system:
        axPDF=PDFfig.add_subplot(311+ti)
        axPDF.set_title('%s'%(j),fontsize=14)
        for t in temp:
            f = open('./%s_%d_single/0/all.dump'%(j,t),"r")
            for frame in range(0, read_frame+skip_frame):
                read_file(coordinates_dense, f, frame)
                gyrations[se][frame], Rg_err=radius_of_gyration(coordinates_dense)
                
            gyration_average[0][se]=np.mean(gyrations[se][skip_frame:read_frame+skip_frame])
            for frame in range(skip_frame, read_frame+skip_frame):
                the_chunk=int(gyrations[se][frame]/chunk_gap)
                PDF[se+1][the_chunk]+=1
            plt.figure("PDF epsilon")
            draw_PDF(PDF,se,experiments[se],t,read_frame)
            plt.figure("compare epsilon", figsize=(5,4))
            plt.annotate(j, xy = (experiments[se]+0.5,gyration_average[0][se]-0.1), fontsize=10 ,**hfont) # 这里xy是需要标记的坐标，xytext是对应的标签坐标
            plt.figure("boxplots epsilon", figsize=(5,4))
            plt.boxplot( gyrations[se], patch_artist=True,showfliers=False,widths=0.45,positions=[se],sym='o', whis=1.0,
                            boxprops=dict(facecolor=color_dict[int(se/2)],  color="black",alpha=0.6),
                            medianprops=dict(color="black"))
            plt.annotate(j, xy = (se,22), fontsize=10 ,**hfont) 
            plt.xticks([])
            se+=1
        ti+=1
    plt.figure("PDF epsilon")
    #plt.xticks(np.arange(10, 60, 5),fontsize=14)
    plt.subplots_adjust(left=None, bottom=None,wspace=None, hspace=0.5) 
    plt.savefig('PDF')
    plt.savefig("singleRg.svg", format="svg")
    for file in os.listdir("./"):
        if file.endswith(".txt") and not file.startswith("total"):
            filename=file
    #plt.savefig("../%s_s.png"%filename, format="png")
    
    #plt.close("PDF epsilon=%s"%para)
    
    plt.figure("compare epsilon", figsize=(5,4))
    plt.xlabel("experiments ",fontdict={ 'size' : 15})
    plt.ylabel("simulation",fontdict={ 'size' : 15})
    plt.title("$R_g (\AA)$",fontdict={ 'size' : 15})
    plt.yticks(np.arange(32, 50, 5),fontsize=12)
    plt.xticks(np.arange(32, 50, 5),fontsize=12)
    plt.xlim(36, 50)
    plt.ylim(36, 50)
    plt.plot([1,60],[1,60], color= 'tab:blue')
    plt.scatter(experiments,gyration_average[0], marker='o',s=60, color= 'tab:blue')
    #np.savetxt('gyration_dilute.txt',gyration_average,fmt='%.2f')
    #plt.savefig('Rg epsilon=%s.png'%para)
    plt.close("compare epsilon")
    print("-------------------------------finished gyration")
    plt.figure("boxplots epsilon", figsize=(5,4))
    plt.ylim(20, 70)
    plt.ylabel("$R_g(simulation) (\AA)$",fontdict={ 'size' : 15})    
    #plt.savefig('Rg boxplots epsilon=%s.png'%para)
    plt.close("boxplots epsilon")
    return plt

