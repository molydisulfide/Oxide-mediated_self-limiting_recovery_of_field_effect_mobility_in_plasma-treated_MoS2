# -*- coding: utf-8 -*-
"""
Oxide-mediated self-limiting recovery of field effect mobility in plasma-treated MoS2
Figure 5 lattice visualisation
Authors Eamonn Weitz and Colin O'Callaghan
Last edited June 2017
Contact email: ocallaco@tcd.ie
"""

import numpy as np
import random as rand
import matplotlib.pyplot as plt

#############################################################################
# We begin with a lattice with sites in one of 3 possible states.
#############################################################################

l1=30   			#   length of matrix
l2=30   			#   width of matrix
n=2*l1*l2  			#   number of transistions
Sites = np.zeros(shape=(l1,l2)) # defines zero matrix of required size

step=10 			# calculate the conductance of netwo_drk every 10 iterations (arbitrary, just for efficiency)
mos2=list()  # number of sites in the initital phase, the MoS2 phase
two_d=list()          # number of sites in the second phase, 2D-MoO3 phase
moo3=list()        # number of sites in the third phase, MoO3 phase


z=l1*l2       # number of zero entries before any balls have been dropped
o=0           # number of one entries before any balls have been dropped
t=0           # number of two_d entries before any balls have been dropped
th=0          # number of moo3 entries before any balls have been dropped 

##############################################################################
### two_d random integers correspond to the location of transformation.      
## Another random number between zero and one is then generated. The transformation
## occurs with some probability or one of its nearest neighbouring site otherwise
## is chosen otherwise. When one of its neighbours is chosen another random number 
## between zero and one is generated and the process repeats until a Site accepts the ball. 
## At this point, transformation=True and the loop terminates.
## This loop is repeated "step" times after which the lattice is visualised.
##############################################################################

p01 = 0.2 # probability of the MoS2 -> 2D-MoO3 transition
p12 = 0.8 # probability of the 2D-MoO3 ->MoO3 transition
for j in range(0,n/step):  # conductance calculated every step number of iterations
                           
    for i in range(0,step):
		# choose random row and column for site to change
        ran1=rand.randint(0,l1-1) 
        ran2=rand.randint(0,l2-1)
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
                    ran1=(ran1+1)%l1
                
                if p01+bp<r<p01+2*bp:
                    ran1=(ran1-1)%l1
                if p01+2*bp<r<p01+3*bp:
                    ran2=(ran2+1)%l2
        
                if p01+3*bp<r<1:
                    ran2=(ran2-1)%l2
        
                s=1
        


            if Sites[ran1][ran2]==1 and s==0: ########
                       #if Site is in state 1

                r = rand.random()
                bp = (1 - p12)/4
                if r<p12:
                    Sites[ran1][ran2]=Sites[ran1][ran2]+1
                    transformation=True
      
            
                if p12<r<p12+bp:
                    ran1=(ran1+1)%l1
                
                if p12+bp<r<p12+2*bp:
                    ran1=(ran1-1)%l1
                if p12+2*bp<r<p12+3*bp:
                    ran2=(ran2+1)%l2
        
                if p12+3*bp<r<1:
                    ran2=(ran2-1)%l2
        
                s=1

    
            if Sites[ran1][ran2]==2 and s==0 :#####
                # Site in final state
      
      
                r = rand.random()
      

                #choose new site		 
                if 0<r<0.25:
                    ran1=(ran1+1)%l1
            
                if 0.25<r<0.5:
                    ran1=(ran1-1)%l1
                if 0.5<r<0.75:
                    ran2=(ran2+1)%l2
                if 0.75<r<1:
                    ran2=(ran2-1)%l2
        
                s=1 
        # After a ball is accepted into a Site, the number of
        # Sites with zero, one, etc. balls is adjusted accordingly. 
        
        if Sites[ran1][ran2]==1:
            z-=1; o+=1
        if Sites[ran1][ran2]==2:
            o-=1; t+=1
        if Sites[ran1][ran2]==3:
            t-=1; th+=1
    
        tot = 900
        mos2.append(float(z+o)/tot)
        two_d.append(float(t)/tot)
        moo3.append(float(th)/tot)

    print j
    plt.clf()
    Sitecop = np.array(Sites[:])
    for x in range(l1):
        for y in range(l2):
            if Sitecop[x][y]==1:
                Sitecop[x][y]=0
            if Sitecop[x][y]==2:
                Sitecop[x][y]=1.5

    img = plt.imshow(Sitecop, interpolation='nearest', cmap='viridis')
    cbar = plt.colorbar(ticks=[0, 1.5, 3], orientation='vertical')
    cbar.ax.set_yticklabels(['MoS$_2$', '2D MoO$_3$', 'MoO$_3$'], fontsize = 20) 
    #plt.colorbar() #Visualises matrix by assigning different colours to 
                    #Sites (entries of matrix) based on how many balls they contain

    plt.savefig("nw_snapshot_iteration_%s.png"%((j+1)*step))

#Want to see how the lists containing the number of 0 or 1's, 2's and 3's
#change as the number of dropped balls increases.
N=list()
for i in range(n):
    N.append(i)

f = open("2d_concentration.txt","w")
for i in range(len(N)):
    f.write("%s %s\n"%(N[i],two_d[i]))
f.close()


