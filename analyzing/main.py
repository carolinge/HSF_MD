# -*- coding:utf-8 -*-
import os,time
import numpy as np
from gyration import gyration_single
from draw_evolutionfigure import draw_evolution,draw_diff
from density import density_distrib,density_cube,draw_phasediagram
#from contacts import get_contacts_single
from NMR import get_contacts
from cluster import get_max_cluster
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
np.set_printoptions(suppress=True)
os.chdir('..')
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['pdf.fonttype'] = 42
hfont = {'fontname':'Arial'}
plt.rcParams['axes.unicode_minus']=False
densities=[]
hist2D=plt.figure("contact 2D",figsize=(10, 10))
thermo=plt.figure("density",figsize=(13, 8))
single=plt.figure("gyration",figsize=(7, 12))
gyseperate=plt.figure("seperate gyration",figsize=(7, 12))
cluster=plt.figure("cluster",figsize=(13, 8))
phase=plt.figure("Phase diagram",figsize=(8, 6))
hist_list=[]
system_all=["WT","M1","M2"]
step_all=[7,7,7]
temps=[278, 288,298,303]
#temps=[278, 278,303,303]
row=0
draw_diff(step_all[0],totalframe=380)
draw_evolution()
for temp in temps:
    print("==========================",temp,"==========================")
    sn=0

    
    axcluster=cluster.add_subplot(2,2,1+row)
    get_max_cluster(axcluster,temp,frame=2000,filename='step_%d'%step_all[sn])
    for spine in axcluster.spines.values():
        spine.set_color('black')
    axcluster.patch.set_facecolor('white')
    print("cluster done......")
    axdensity=thermo.add_subplot(2,2,1+row)
    densities.append(density_distrib(300, axdensity,temp,totalframe=50,filename='step_%d'%step_all[sn]))
    for spine in axdensity.spines.values():
        spine.set_color('black')
    axdensity.patch.set_facecolor('white')
    print("density done......")
    for system in system_all:
        skip_frame=100; read_frame=500
        axsingle=single.add_subplot(3,1,1+sn)
        gyration_single(axsingle,skip_frame, read_frame,fromfile=1,filename='step_%d/%s_%d'%(step_all[sn],system,temp),t=temp,showraw=0,system=system,filename2='step_%d'%step_all[sn])#skip_frame, read_frame
        axsep=gyseperate.add_subplot(3,2,2*(1+sn)-1)
        if temp==278 or temp==303:
            axsep=gyseperate.add_subplot(3,2,2*(1+sn))
        gyration_single(axsep,skip_frame, read_frame,fromfile=1,filename='step_%d/%s_%d'%(step_all[sn],system,temp),t=temp,showraw=0,system=system,filename2='step_%d'%step_all[sn])#skip_frame, read_frame            

        axnow=hist2D.add_subplot(4,3,1+sn+row*3)
        
        hist_list.append(get_contacts(axnow,system, temp,filename='step_%d'%step_all[sn],frame=380,r_cut=8.5))
        sn+=1
    row+=1
    print("gyration done......")


plt.figure("contact 2D")
plt.legend()
plt.savefig('contactmap_step%d.pdf'%(step_all[0]),format="pdf")
plt.subplots_adjust(left=None, bottom=None,wspace=0.2, hspace=0.2)
plt.tight_layout()
plt.figure("Phase diagram")
axphase=phase.add_subplot(1,1,1)
axphase.patch.set_facecolor('white')
densities=np.array(densities).T
def c_2_phi(rho,c):
    """convert protein concentration in mg/mL into volume fraction
    1 mgmL is eq to 1 kg/m^3
    density is in SI units of kg/m^3  
    """
    return c / rho
for s in range(2):
    draw_phasediagram(axphase,densities[2*s:2*s+2],temps,system_all[s])
data1=[[ 5.,3.03,253.0 ], [10.  , 2.59 , 257.33], [15. ,   2.26,  264.33],
       [20.   , 2.17,  286.67], [25.    ,2.21,  290.33], [30.  ,  1.92 , 311.00], [35.   , 1.98,  321.67 ]]
data0=[[ 5.  ,   5.49,  5.49],[ 10.  ,   4.49,  66.53], [ 15.   ,  4.07 , 90.  ], [ 20.  ,   3.  , 119.73],
       [ 25.  ,   2.59 ,139.7], [ 30. ,    2.23 ,176.8 ], [ 35.  ,   1.93, 227.2 ]]
data=np.array(data0)
T = data[:,0]+273.1
phi1= c_2_phi(1000, data[:,1])
phi2= c_2_phi(1350, data[:,2])#density of protein
axphase.scatter(phi2,T,marker='o',s=60,color='gray',label='WT',alpha=0.6)#水0
axphase.scatter(phi1,T,marker='o',s=60, color='gray',alpha=0.6)

data = np.array(data1)
T = data[:,0]+273.15
phi1= c_2_phi(1000, data[:,1])
phi2= c_2_phi(1350, data[:,2])#density of protein
axphase.scatter(phi2,T,marker='o',s=60,color="#ED797B",label='M1',alpha=0.6)
axphase.scatter(phi1,T,marker='o',s=60,color="#ED797B",alpha=0.6)
for spine in axphase.spines.values():
    spine.set_color('black')

plt.legend()
plt.savefig('phasediagram_step%d.pdf'%step_all[0],format="pdf")
plt.subplots_adjust(left=None, bottom=None,wspace=0.2, hspace=0.2)
plt.tight_layout()

plt.figure("density")
plt.legend()
plt.savefig('density_step%d.pdf'%step_all[0],format="pdf")
plt.subplots_adjust(left=None, bottom=None,wspace=0.2, hspace=0.2)     
plt.tight_layout()

plt.figure("gyration")
plt.legend()
plt.savefig('gyration_step%d_%d.pdf'%(step_all[0],read_frame),format="pdf")
plt.legend(fontsize=10)
plt.subplots_adjust(left=None, bottom=None,wspace=None, hspace=0.5) 
#plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
plt.figure("seperate gyration")
plt.legend()
plt.savefig('gyrationseparate_step%d_%d.pdf'%(step_all[0],read_frame),format="pdf")
plt.legend(fontsize=10)
plt.subplots_adjust(left=None, bottom=None,wspace=None, hspace=0.5) 
#plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)

plt.figure("cluster")
plt.savefig("cluster_step%d.pdf"%step_all[0],format="pdf")
plt.subplots_adjust(left=None, bottom=None,wspace=0.2, hspace=0.2)     
plt.tight_layout()

plt.show()

