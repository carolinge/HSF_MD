# -*- coding:utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from NMR_figure import get_hsqc_exp,get_hsqc_mask
total_step=6

system_all=["WT","M1","M2"]
temps=[278,288,298,303]

os.chdir('../')


evolve=plt.figure("iteration error",figsize=(10,8))
lambdas=plt.figure("lambdachange",figsize=(15,6))
totforce=plt.figure("totalinteract",figsize=(8,4))
axtot=totforce.add_subplot(1,1,1)
totalforce=[]
totallabel=[]
totalNMR=[]
np.set_printoptions(suppress=True)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['pdf.fonttype'] = 42
hfont = {'fontname':'Arial'}
plt.rcParams['axes.unicode_minus']=False

sn=0
for system in system_all:
    row=0
    tempslambda=np.zeros((4,181))
    for temp in temps:
        errors=np.zeros(total_step)
        ax=evolve.add_subplot(3,4,1+row+sn*4)
        
        ax.set_title("%s %dK"%(system,temp),fontdict={ 'size' : 12},**hfont)
        ax.set_ylabel(f'error')
        ax.set_xlabel('iteration step')
        ax.set_xlim([-0.5,total_step+0.5])
        plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.15,hspace=0.2,wspace=0.2)
        #ax.set_ylim([0.0,0.17])
        #ax.set_yticks(size=12,**hfont)
        #ax.set_xticks(size=12,**hfont)
        ax.locator_params(nbins=5)
        
        read=np.loadtxt("annealing/errors_%s_%d.txt"%(system,temp))
        for i in range (len(read[0])):
            try:
                errors[int(read[0][i])-1]=read[1][i]
            except:
                pass
        ax.plot(np.arange(total_step),errors, marker='o') 
        
        
        #draw lambda change compared to larson
        before=np.loadtxt("step_0/%s_%d_lambdas.data"%(system,temp))
        after=np.loadtxt("step_%d/%s_%d_lambdas.data"%(total_step,system,temp))
        ax=lambdas.add_subplot(4,3,1+sn+row*3)
        ax.set_title("%s %dK"%(system,temp),fontdict={ 'size' : 12},**hfont)
        ax.set_ylabel(r'$\lambda_{NMR}-\lambda_0$')
        #ax.set_xlabel('iteration')
        #ax.set_xlim([-0.5,total_step+0.5])
        plt.subplots_adjust(left=0.15,bottom=0.15,hspace=0.4,wspace=0.4)
        #ax.set_ylim([0.0,0.17])
        plt.yticks(fontsize=12,**hfont)
        plt.xticks(fontsize=12,**hfont)
        ax.bar(np.arange(181),after-before,alpha=0.6)        
        #ax.locator_params(nbins=5)       
        tempslambda[row]=after
        totalforce.append(np.mean(after))
        totallabel.append("%s_%d"%(system,temp))
        intexp=get_hsqc_exp(temp,system)
        totalNMR.append(1-np.mean(intexp))
    
        
        row+=1 
    
        
    #effect of temperature
    '''ax=lambdas.add_subplot(4,3,1+sn+row*3)
    ax.set_title("%s %dK"%(system,temp),fontdict={ 'size' : 12},**hfont)
    ax.set_ylabel(r'$\lambda_{NMR}-\lambda_0$')
    #ax.set_xlabel('iteration')
    #ax.set_xlim([-0.5,total_step+0.5])
    plt.subplots_adjust(left=0.15,bottom=0.15,hspace=0.4,wspace=0.4)
    #ax.set_ylim([0.0,0.17])
    plt.yticks(fontsize=12,**hfont)
    plt.xticks(fontsize=12,**hfont)
    ax.bar(np.arange(181),after-before,alpha=0.6)        
    #ax.locator_params(nbins=5)       '''
        
    sn+=1
axtot.bar(totallabel,totalforce,label="simulated totalforce")
axtot.bar(totallabel,totalNMR,label="1-NMR")
axtot.legend()
for i, v in enumerate(totalforce):
    plt.text(totallabel[i], v + 0.02, "%.2f"%v, ha='center')
for i, v in enumerate(totalNMR):
    plt.text(totallabel[i], v + 0.02, "%.3f"%v, ha='center')

plt.tight_layout()
plt.savefig('iteration.pdf',format="pdf")

plt.show()


