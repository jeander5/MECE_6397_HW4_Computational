# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:37:40 2020

@author: johna
"""

#I used MATLAB last time so I will use python this time

# MECE 6397, SciComp, Problem 4, Computational
#Solve the helmholts equation for 1. Dirchlet. 2. Nuemann

#lets try and do this efficiently. and with correct formatting python style
#Note the Code Analyzer does NOT LIKE the lower case constant name or X_whatever,
#but im just gonna plow ahead


#imports
import numpy as np
#import math
import matplotlib.pyplot as plt
#import sympy as sym
#import scipy
#assigning shortcuts from imports
#from random import randint as randi
from math import sinh as sinh
#import time
#start_time=time.time()
#print('Start')
#print("--- %s seconds --" % (time.time()-start_time))

#Contstants given in problem statement, constant for both boundary conditions.
#Interval length, u(x=0), v is constant on hyperbolic sin function, A is exact solution of f(x).
# solve both problem for both values of k
L = 1
U_o = 1
v = 1
A = 1
K = [1, 10]

#Discretizing the interval length.
N = 10
h = L/(N+1)
x = np.linspace(0, L, N+2)

#Part 1, Dirchlet
#I will just start off with k=1
k = K[0]

#getting all u values, gonna do this all at once outide of any loop. 
#This Includes the boundaries
#u_exact = [((sinh(k*(L-x))+sinh(k*x))/sinh(k*L)-1)*A/k**2+U_o*sinh(k*(L-x))/sinh(k*L) for x in x]
#This does not
u_exact = [((sinh(k*(L-x))+sinh(k*x))/sinh(k*L)-1)*A/k**2+U_o*sinh(k*(L-x))/sinh(k*L)
for x in x[1:-1]]

#Pre Thomas algorith set up. for this problem these values are all constant
a =-(2+k**2*h**2)
b = 1
c = 1
#should I call this f a different varible? nah its ok
f = A*h**2

#should I use an array? nah,Im not actually gonna do any matrix stuff really
#lists should be fine, I will still create them ahead of time
alpha = [0]*N
g = [0]*N
u_appx = [0]*N
#I hope thats not inefficient

#Thomas algorith, folliwng the psuedo code
#zeroth element of this list corresponds to the first subscript in thomas algorith 
alpha[0] = a
g[0] = f-U_o

for j in range(1, N):
    alpha[j] = a-(b/alpha[j-1])*c
    g[j] = f-(b/alpha[j-1])*g[j-1]
    
u_appx[N-1] = g[N-1]/alpha[N-1]

for m in range(1,N):
    u_appx[-1-m] = (g[-1-m]-c*u_appx[-m])/alpha[-1-m]

#all right I am in business,
#small little mix up with the sign on the k and lamda term from the Helmhotz eq
#I knew I wrote that algorithm correctly
    
#lets do the grid convergence and the error stuff next
#will that have to be like in an outter loop?, the way im invisioning it is like an outter while loop
#I also need to submitt graphs and compare values so maybe not a while loop
 
# I think I will plot now, then start while loop after the initial calulation, do it again for 2*N, and compare???
# no then I will have the same lines of code just reapeated
# maybe make some functions.
#I gotta do k=10 as well
#I dont care how many physical lines it is I just want it to be computer efficient.
#If I made a thomas algoithm function I would want to feed in the a,b,c terms, cause what if they arent always the same?.
#im so indecisve
#plt.plot(x[1:-1],u_exact)
#I think making a Thomas algorithm function is the way to go, maybe not though.
#I should probably just make everything functions
# maybe not I cant decide.
# oh and we are supposed to use git.
#with functions I can easily change variables with out searching through the whole code
#will be nice for touble shooting    
    
    
