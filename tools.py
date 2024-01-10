# -*- coding:utf-8 -*-
import re
import numpy as np
import pandas as pd
from numba import jit
import math
import os

read_frame=40#350
skip_frame=1
#--------初始化
temps=[278,298]
total_chain=40
bead_per_chain=181
total_atom=total_chain*bead_per_chain

box_len=np.zeros((read_frame+skip_frame ,3))

#########################################函数区###########################################################

color_dict=[ '#14517c', '#2f7fc1', '#c7dfda', '#96c37d', '#f3d266', '#d8383a', '#e7b1bd', '#f1e1e1', '#c497b2', '#a8b8c6']
color_dict=['#f3d266','#96c37d','#c1d1da', '#2f7fc1', '#24517c' ]

#########################################函数区###########################################################

def read_file(coordinates,box_len, f, frame):
    for i in range(0, 5):
        f.readline()
    for i in range(0, 3):   
        match= re.findall(r'(-?[0-9\.e+]+)', f.readline())
        box_len[frame][i]=float (match[1])
    if frame%20==0: 
        print("frame: ",frame, box_len[frame][0],box_len[frame][1],box_len[frame][2])
    f.readline()  
    for i in range (0, total_atom):
        line = f.readline()            
        match=line .split()
        if float(match[0])<total_atom:
            coordinates[int(match[0])-1][0]=int(match[1])
            for x in range (0,3):                
                coordinates[int(match[0])-1][x+1]=float(match[x+2])*box_len[frame][x]

@jit(nopython=True)                 
def reduce(matrix,value):
    for x in range (0,len(matrix)):  
        for y in range(0, b): 
            matrix[x][y]/=value


