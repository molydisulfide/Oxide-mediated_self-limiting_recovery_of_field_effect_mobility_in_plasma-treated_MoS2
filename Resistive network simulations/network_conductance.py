# -*- coding: utf-8 -*-
"""
Oxide-mediated self-limiting recovery of field effect mobility in plasma-treated MoS2
Figure 5 data generation
Simulate the lattice conductance. 
Authors Eamonn Weitz and Colin O'Callaghan
Last edited June 2017
Contact email: ocallaco@tcd.ie
"""

import numpy as np
import random as rand
#from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import time
import math as mt
from cvxopt.base import matrix as m
from symmlq import *
import itertools

##############################################################################
### Needed for symmlq linear equation solver 
###########################################################################  
def Gfun(x,y,trans='N'):
    '''Function that passes matrix A to the symmlq routine which solves Ax=B.'''
    gemv(Amatrix,x,y,trans)

tol=1e-10
show=False
maxit=None


Lead=30 #length of lattice
W=30    #Width of lattice
n=W*Lead*2 #number of site transformations
step=10 #The conductanc of the lattic will be measured after every 10 steps

Kirch_Dim=2*(Lead*W)-(Lead+W)   # dimesnion of the Kirchoff matrix
l_pairs=2*Lead-1 #number of links attached to one column of nodes on the lattice

a=1    #lattice spacing
C0=1.0  #Conductance of site in MoS2 phase
C1=10000.0  #Conductance of site in 2D-MoO3 phase
C2=0.01   #Conductance of site in MoO3 phase


Resistance=[0]*(n/step) # resistance list

for i in range(0,(n/step)):
    Resistance[i]=list()
Conductance=list()

p01 = 0.2 # probability of the MoS2 -> 2D-MoO3 transition
p12 = 0.8 # probability of the 2D-MoO3 ->MoO3 transition

# note that p01 and p12 are arbirary but want p01<p12

Sites = np.zeros(shape=(Lead,W)) #Lattice is redefined to be composed of only zeroes

data_file = open("Conductance_curve_p01_%s_p12_%s.txt"%(p01,p12),"w") # output file

for b in range(0,n,step):
    t=0
    coords =[] 
    pairs=[]
    N_Pairs=len(pairs)

	#############################################################################
	### This loop assigns to each node of the lattice a point in Euclidean     ##
	## space. The lower left corner of the lattice is the orgin.###########
	############################################################################
    for i in range(W):
        for j in range(Lead):
            coords.append([])
            coords[t].append(i*a)
            coords[t].append(j*a)
            t+=1
	#############################################################################
	## Here, each link in the lattice is identified by the two nodes which #####
	## it connects. These are stored in the pairs list.#########################
	#############################################################################       

    for i in range(W):
        sub_pairs=[]
        if i==W-1:
            for j in range(Lead-1):
                sub_pairs.append([])
            for j in range(Lead-1):
                sub_pairs[j].append(j+Lead*i)
                sub_pairs[j].append(j+1+Lead*i) 
            for j in range(Lead-1):
                pairs.append(sub_pairs[j])
    
        else:
            for j in range(l_pairs):
                sub_pairs.append([])         
            for j in range(Lead-1):
                sub_pairs[j].append(j+Lead*i)
                sub_pairs[j].append(j+1+Lead*i)
            for j in range(Lead):
                sub_pairs[j+Lead-1].append(j+Lead*i)
                sub_pairs[j+Lead-1].append(j+Lead+Lead*i)
            for j in range(l_pairs):
                pairs.append(sub_pairs[j])
    new_coords=[]
    new_pairs=[]
    list1=[]
    list2=[]
    N_Pairs=len(pairs)
    for i in range(len(coords)): # Construct these as nested lists
        
        new_coords.append([])
        list1.append([])        
        list2.append([])
	##############################################################################
	### A list entry in new_coords contains the link connected to a certain node##
	## as well as the overall number of links.################################
	##############################################################################
    for i in range(N_Pairs):
        for j in range(N_Pairs):
            if i in pairs[j]:
                new_coords[i].append(pairs[j]) 
            else:
                pass
    for i in range(N_Pairs):    #need to number pairs by node number for when dealing with Kirchoff matrix later
        pairs[i].append(i)

	##############################################################################
	#### For each node on the lattice, the ith entry of list1 tells us what are ###
	### number coordinates of the ith node's nearest neighbours. Note this coordinate is not##
	## the same as the one contained in the coords list. list2 contains the ######
	#### numbers of the links connected to the ith node.##########################
	##############################################################################
    for i in range(len(coords)):
        for j in range(len(new_coords[i])):
            for k in range(0,2):
                if new_coords[i][j][k]!=i:
                    list1[i].append(new_coords[i][j][k])
                else:
                    pass          
            list2[i].append(new_coords[i][j][2])
            
    for i in range(len(coords)):
        g=list()
        d=list(itertools.combinations(list1[i],2))
        f=list(itertools.combinations(list2[i],2))      #by finding all possible pairs in list1 and list2, have all links for new_coords
        g.append(d)
        g.append(f)
        new_pairs.append(g)
    #print new_pairs
    
    for q in range(0,step):
   

		# choose random row and column for site to change
        ran1=rand.randint(0,Lead-1) 
        ran2=rand.randint(0,W-1)
        transformation=False             
        while transformation==False:
        
            s=0
            if Sites[ran1][ran2]==0 and s==0: ########
                           #if site is in state 0
        
                r = rand.random()
                bp = (1 - p01)/4 
                if r<p01: # transform site
                    Sites[ran1][ran2]=Sites[ran1][ran2]+1 # change the site value
                    transformation=True
      			# otherwise choose a nearest neighbour at random
            
                if p01<r<p01+bp:
                    ran1=(ran1+1)%Lead
                
                if p01+bp<r<p01+2*bp:
                    ran1=(ran1-1)%Lead
                if p01+2*bp<r<p01+3*bp:
                    ran2=(ran2+1)%W
        
                if p01+3*bp<r<1:
                    ran2=(ran2-1)%W
        
                s=1
        


            if Sites[ran1][ran2]==1 and s==0: ########
                       #if box is in state 1

                r = rand.random()
                bp = (1 - p12)/4
                if r<p12:
                    Sites[ran1][ran2]=Sites[ran1][ran2]+1
                    transformation=True
      
            
                if p12<r<p12+bp:
                    ran1=(ran1+1)%Lead
                
                if p12+bp<r<p12+2*bp:
                    ran1=(ran1-1)%Lead
                if p12+2*bp<r<p12+3*bp:
                    ran2=(ran2+1)%W
        
                if p12+3*bp<r<1:
                    ran2=(ran2-1)%W
        
                s=1

    
            if Sites[ran1][ran2]==2 and s==0 :#####
                # box in final state
      
      
                r = rand.random()
      

                #choose new site		 
                if 0<r<0.25:
                    ran1=(ran1+1)%Lead
            
                if 0.25<r<0.5:
                    ran1=(ran1-1)%Lead
                if 0.5<r<0.75:
                    ran2=(ran2+1)%W
                if 0.75<r<1:
                    ran2=(ran2-1)%W
        
                s=1 
                
    for i in range(len(coords)):
        for j in range(0,2):
            for k in range(len(new_pairs[i][j])):
                new_pairs[i][j][k]=list(new_pairs[i][j][k]) #need new_pairs to be in list form 
    for i in range(len(coords)):
        coords[i].append(Sites[coords[i][0]][coords[i][1]])
        coords[i].append(i)
       #Gives links between new_coords conductance values based on phase
        if Sites[coords[i][0]][coords[i][1]]==0:
            for j in range(len(new_pairs[i][0])):
                new_pairs[i][0][j].append(C0)         # phase zero with conductance C0
        elif Sites[coords[i][0]][coords[i][1]]==1:
            for j in range(len(new_pairs[i][0])):
                new_pairs[i][0][j].append(C1)          # phase zero with conductance C1
        elif Sites[coords[i][0]][coords[i][1]]==2:
            for j in range(len(new_pairs[i][0])):        
                new_pairs[i][0][j].append(C2)          # phase zero with conductance C2





    Kirchoff_Mat = np.zeros((Kirch_Dim,Kirch_Dim)) # Kirchoff matrix



    for i in range(len(pairs)):
        for j in range(len(coords)):
            for k in range(len(new_pairs[j][1])):
                if i in new_pairs[j][1][k]:
                    x=new_pairs[j][1][k][0]                 #Builing up Kirchoff Matrix
                    y=new_pairs[j][1][k][1]
                    Kirchoff_Mat[x][y]=-new_pairs[j][0][k][2]
                    Kirchoff_Mat[y][x]=Kirchoff_Mat[x][y]
                else:
                    pass

    np.fill_diagonal((Kirchoff_Mat), abs(Kirchoff_Mat.sum(1)))
    I = 10.0 # total current being pumped through the system
    i0 = I/(Lead) # the current that will be put into each node on the left lead and is then extracted on the right lead
    ic = np.zeros(Kirch_Dim) # current vector
    for i in range(0,Lead-1):	
        # This loop fills in the current vector
        ic[i] = +i0 # an injection node
        ic[(Kirch_Dim)-i-1] = -i0 # an extraction node
    Kirchoff_Matrix = m(Kirchoff_Mat) # This just changes the format of Kirchoff_Mat so it can be solved by the symmlq routine
    Imatrix = m(ic) # puts ic in the format for symmlq
    # Solving Ax = B (Mr * V = I) with symmlq
    Amatrix = Kirchoff_Matrix # renamed Amatrix for Gfun function above
    elec_pot_mr = symmlq( Imatrix, Gfun, show=show, rtol=tol, maxit=maxit) # result from symmlq. elec_pot_mr[0] is the voltage vector. 
    V_left=list()
    V_right=list()

    for k in range(0,Lead-1):
        # adding total voltages on each lead
    
        V_left.append(elec_pot_mr[0][k]) 
        V_right.append(elec_pot_mr[0][(Kirch_Dim)-1-k]) 

    L=np.mean(V_left)
    R=np.mean(V_right)
    conductance = I/abs(L - R)
    if b==0:
        c_norm = conductance # normalise the conductance of network so it starts at 1
	#conductance = conductance/c_norm
    print conductance,conductance/c_norm
    data_file.write("%s\t%s\n"%(b,conductance/c_norm))


data_file.close()



