# -*- coding:utf-8 -*-
import os
import numpy as np
from gyration import gyration_single
from density import density_distrib,density_cube
from contact import contact_map
from NMR import get_contacts
import matplotlib.pyplot as plt

os.chdir('../')

for ex in range(7,8):
#for ex in range(5,9):
    os.chdir('step%d/'%(ex))
    skip_frame=400; read_frame=1200
    #plt=gyration_single(skip_frame, read_frame)#skip_frame, read_fram


    
    #thermo=plt.figure("density",figsize=(8, 6))   
    #density_distrib(300,thermo)
    
    #gyration_single(skip_frame, read_frame)#skip_frame, read_frame
    #contact_map(hist2D,fromfile=0,frame=20,r_cut=8.5,reduce_method="root3",vvmin=0.0,vvmax=0.5)
    
    
    ratio=plt.figure("dense/dilute NMR",figsize=(13, 8))

    hist2D=plt.figure("contact 2D",figsize=(10, 5.8))
    get_contacts(ratio,hist2D,fromfile=1,frame=300,r_cut=8.5)
    

    #plt.savefig("gyration.pdf", format="pdf")
    os.chdir('../')
    plt.show()
    #plt.clf()#傻逼傻逼 找不到bug记得看这里

