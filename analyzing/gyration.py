
# -*- coding:utf-8 -*-
import re
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline as spline
from scipy import interpolate
import math


np.set_printoptions(suppress=True)

hfont = {'fontname':'Arial'}

system = [ "WT","M1","M2"]
temp=[278,298]
experiments=[49.8,49.8, 36.1,36.1,56.2,56.2]
sys_dict = {
    'WT': 49.8,
    'M1': 36.1,
    'M2': 56.2
}
color_dict = {
    'WT_278': '#E48E40','WT_288':'#F3CBA7','WT_298':'#F3CBA7','WT_303':'#F3CBA7',
    'M1_278': '#C83232','M1_288':'#ED797B','M1_298':'#ED797B','M1_303':'#ED797B',
    'M2_278': '#3E7BA4','M2_288':"#21A6EC",'M2_298':"#21A6EC",'M2_303':"#21A6EC",
}
color_dict = {
    #'WT_278': '#E48E40','WT_288':'#E6BA40','WT_298':'#E7D640','WT_303':'#E8E040',
    'WT_278': '#D07F61','WT_288':'#D0A362','WT_298':'#D0BA62','WT_303':'#D0C362',
    'M1_278': '#C83232','M1_288':'#ED797B','M1_298':'#ED797B','M1_303':'#ED797B',
    'M2_278': '#3E7BA4','M2_288':"#21A6EC",'M2_298':"#21A6EC",'M2_303':"#21A6EC",
}
dash_dict = {
    '278': 'solid',
    '288': 'dashed',
    '298': 'dashdot',
    '303': 'dotted',
}
# 查询不同温度对应的实验数据


epsilon=["0.11","0.115","0.12","0.13"]
epsilon=["0.115"]
how_many_sys=len(system)*len(temp)
total_chain=1
bead_per_chain=181
total_atom=total_chain*bead_per_chain
coordinates_dense = np.zeros((total_atom+1 ,5))

gyration_average=np.zeros((4 ,how_many_sys))
#PDF
chunk_gap=8.0
total_chunk=30
PDF=np.zeros((how_many_sys+1,total_chunk))

for i in range (0,total_chunk):
    PDF[0][i]=i*chunk_gap+0.5*chunk_gap+3.0
    

def read_file(coordinates,f, frame):
    box_len=np.zeros((1 ,3))
    for i in range(0, 5):
        f.readline()
    for i in range(0, 3):   
        match= re.findall(r'(-?[0-9\.e+]+)', f.readline())
        box_len[0][i]=float (match[1])
    if frame%1000==0:
        print(" | frame: ",frame, box_len[0][0],total_atom)
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


def draw_PDF(axnow,PDF,se,system,t,read_frame,showraw=1):
    PDF_smoothed=np.zeros((1,total_chunk))
    for i in range(0,total_chunk):
        PDF[se+1][i]/=read_frame
    for i in range(1,total_chunk-1):
        PDF_smoothed[0][i]=(PDF[se+1][i-1]+ PDF[se+1][i]+PDF[se+1][i+1] )/3
        
    xnew = np.linspace(10,PDF[0].max(),300)
    tck=interpolate.splrep(PDF[0],PDF[se+1])
    y_spline=interpolate.splev(xnew,tck)
    if t==278:
        axnow.plot(xnew,y_spline,label='T=%d'%t, color= color_dict["%s_%d"%(system,t)])
        if showraw:
            axnow.plot(PDF[0],PDF[se+1],label='T=%d'%t, color= color_dict["%s_%d"%(system,t)])
            axnow.fill_between(xnew, 0, y_spline, facecolor= color_dict["%s_%d"%(system,t)], alpha=0.35)
    else:
        axnow.plot(xnew,y_spline,label='T=%d'%t, color= color_dict["%s_%d"%(system,t)],alpha=1.0-(t-288)/40)
        if showraw:
            axnow.plot(PDF[0],PDF[se+1],label='T=%d'%t,color= color_dict["%s_%d"%(system,t)],alpha=0.75)
        axnow.fill_between(xnew, 0, y_spline, facecolor= color_dict["%s_%d"%(system,t)], alpha=0.3-(t-288)/100)

    axnow.set_xlabel("$R_g (\AA)$",fontdict={ 'size' : 15})
    axnow.set_ylabel("PDF",fontdict={ 'size' : 15})
    axnow.set_xlim(15, 70.0)
    axnow.set_ylim(0, 0.6)
    axnow.set_title(system)
    axnow.set_xticks(np.arange(15, 60, 10))
    axnow.tick_params(axis='x', labelsize=12) 
    axnow.legend()
    
    if t==278:
        axnow.axvline(x=gyration_average[0][se],ls="--", lw=1.2, c=color_dict["%s_%d"%(system,t)])#添加竖直直线
        #axnow.axvline(x= sys_dict.get(system), lw=1,label='exp T=278k', c="black")#添加竖直直线
    else:
        1
        axnow.axvline(x=gyration_average[0][se],ls="--", lw=1.2, c=color_dict["%s_%d"%(system,t)], linestyle=dash_dict['%d'%t], alpha=1.0-(t-288)/30)#添加竖直直线
        #axnow.set_text(gyration_average[0][se]-6, 0.5, r'SAXS',fontsize=10)
    #axnow.set_yticks(np.arange(0.1, 0.8, 0.1),fontsize=12)
    #axnow.set_xticks(np.arange(10.0, 26.0, 3.0),fontsize=12)


def gyration_single(axnow,skip_frame, read_frame,fromfile=1,filename="",showraw=1,filename2="",t=278,system=""):
    se=0
    gyrations=np.zeros((how_many_sys ,read_frame+skip_frame))
    ti=0

    
    try:
        gyrations=np.loadtxt("%s/processfile/gyraion_%s_%s_%d-%dframe.txt"%(filename2,system,t,skip_frame,read_frame),delimiter=' ')
        if not fromfile:
            ddd
    except:
        print("could not find processed file on", os.getcwd())
        f = open('./%s_single/0/all.dump'%(filename),"r")
        for frame in range(0, read_frame+skip_frame):
            read_file(coordinates_dense, f, frame)
            gyrations[se][frame], Rg_err=radius_of_gyration(coordinates_dense)
        np.savetxt("%s/processfile/gyraion_%s_%s_%d-%dframe.txt"%(filename2,system,t,skip_frame,read_frame),gyrations,fmt='%.3f')  
        f.close()
    gyration_average[0][se]=np.mean(gyrations[se][skip_frame:read_frame+skip_frame])
    for frame in range(skip_frame, read_frame+skip_frame):
        the_chunk=int(gyrations[se][frame]/chunk_gap)
        PDF[se+1][the_chunk]+=1
    draw_PDF(axnow,PDF,se,system,t,read_frame,showraw=showraw)
    plt.figure("trajectory")
    plt.plot(np.arange(0,len(gyrations[0])),gyrations[0])

    #plt.annotate(j, xy = (se,22), fontsize=10 ,**hfont) 
    '''plt.xticks([])
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
'''
