# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 14:07:44 2020

@author: johna
"""
#https://github.com/jeander5/MECE_6397_HW4_Computational
#Okay I am making a new file and will be working from this version from now on using 
#Github to document all changes.


# MECE 6397, SciComp, Problem 4, Computational
#Solve the helmholtzs equation for 1. Dirchlet. 2. Nuemann

#imports
import numpy as np

import math
#import matplotlib.pyplot as plt
from math import sinh as sinh
from math import log as log
from random import randint as randi

#Contstants given in problem statement, constant for both boundary conditions.
#Interval length, u(x=0), v is constant on hyperbolic sin function, A is exact solution of f(x).
# solve both problem for both values of k

L = 1
U_o = 1
v = 1
A = 1
K = [1, 10]

#The N value, this is gonna change
#N must be greater than 2
N=17

#helmotlz dirchlet part 1 problem

#discretize the interval function     
def DIF(L,N):
#Discretizing the interval length. This is the same for both problems
    h = L/(N+1)
    x = np.linspace(0, L, N+2)
    return(x[:],h)

#Helmholtz Thomas Algorith Function    
#I still need that A as an input for this problem 
def HTAF(N,h,lamda,U_o,A):   
#inputs are N, lamda, U_o=u(x=0), and for this problem A.
#Pre Thomas algorith set up. for this problem these values are all constant
# Note these values are now inside the function
    a =-(2-lamda*h**2)
    b = 1
    c = 1
    # the following line would be changed/removed or I would need a vector f as in input
    #depending on the problem
    f = A*h**2
    alpha = [0]*N
    g = [0]*N
    u_appx = [0]*N
#following the psuedo code
#zeroth element of this list corresponds to the first subscript in thomas algorith
    alpha[0] = a
    g[0] = f-U_o
    for j in range(1, N):
        alpha[j] = a-(b/alpha[j-1])*c
        g[j] = f-(b/alpha[j-1])*g[j-1]   
    u_appx[N-1] = g[N-1]/alpha[N-1]
    for m in range(1,N):
        u_appx[-1-m] = (g[-1-m]-c*u_appx[-m])/alpha[-1-m]
    return(u_appx)

#u exact function, for the helmhotlz dirchlet part 1 problem
def uEF(k,L,x,A,U_o):
    u_exact = [((sinh(k*(L-x))+sinh(k*x))/sinh(k*L)-1)*A/k**2+
U_o*sinh(k*(L-x))/sinh(k*L) for x in x[1:-1]]
    return(u_exact)

#Note:Functions are now defined. Moving along

 
#Calling the Discr. the Interval right here for now     
x,h=DIF(L,N)   
 
#im gonna need this eventually    
#lenK=len(K)
#for n in range(lenK)
#    k = K[n]

#just gonna start with k=1 for now
k = K[0]

#Note: lamda in the helmholtz eq defined here
lamda=-k**2

#Part 1, Dirchlet

#Grid Convergence
#im using the flag so I dont have to call the function before and inside the while statement
Flag = 0
m=1
#Lebron is a placeholder variable name
Lebron=1*10**-3
#Checking a random discretized point, not every value.
check_val=randi(0,N)
while Flag == 0:
#calling all my functions
    print(N)
    x, h = DIF(L,N)  
    u_appx=HTAF(N,h,lamda,U_o,A)
    u_appx_next = HTAF(2*N,L/(2*N+1),lamda,U_o,A)
    if abs(u_appx[check_val]-u_appx_next[2*check_val])<Lebron:
        Flag = 1
        print('%s Grid Points Needed' %(N))
        print('Doubling the Grid Points would result in less then %s differece between u values for the same Grid Point' %(Lebron))
    else:
        N=N+N


#I think this is now actually doing the job pretty good. If N is large enough it will return the original N value.
#I still think that Flag variable is gonna be inefficeint. 
#Still for not perfect tho. for example consider Lebron=1*10**-3 and N=767 N=768, N=769, produce very different results.
#so its kinda sensitive near key points
#just move on.        

       
#Note: here is the exact value function call
u_exact = uEF(k, L, x, A, U_o)
x1, h1 = DIF(L,N)   
x2, h2 = DIF(L,2*N)
u_exact_next = uEF(k, L, x2, A, U_o)     

#Next....Checking relative error, graphs and tables, then on to part two, 

#Ok relative error, sure function why not, nah I dont even need to store that. I just need that formal order equaition.
# I will just make them 
#formal order of accuracy.
#log 2 error term
#doing it like this so the computer doesnt have to calulate this value everytime.
log2=log(2)
#lets create the storage list ahead of time.fooa= formal order of accuracy... for now.
fooa=[0]*N
rel_err_1=[0]*N
rel_err_2=[0]*N
fooa3=[0]*N
for k in range(N):
    rel_err_1[k]=abs((u_appx[k]-u_exact[k])/u_exact[k])
    rel_err_2[k]=abs((u_appx_next[2*k]-u_exact_next[2*k])/u_exact_next[2*k])
    fooa[k]=(1/log2)*(log((abs((u_appx[k]-u_exact[k]))/u_exact[k]))
-log((abs((u_appx_next[2*k]-u_exact[k]))/u_exact[k])))
    fooa3[k]=log(rel_err_1[k]/rel_err_2[k])/log2

#im gonna keep these lines here.
#the stuff below here is just for copying and pasting into the command line for troubleshooting    
#call all functions
#def CAF(N,lamda,U_o,A,k, L):
#    x, h = DIF(L,N) 
#    u_appx=HTAF(N,h,lamda,U_o,A)
#    u_appx_next = HTAF(2*N,h/2,lamda,U_o,A)        
#    u_exact = uEF(k, L, x, A, U_o)
#    check_in_AF=randi(1,N-1)
#    print('%s Grid Points %(N)')
#    print('u appx=%s'%(u_appx[check_in_AF]))
#    print('u appx_next=%s'%(u_appx_next[2*check_in_AF]))
#    print('u exact=%s'%(u_exact[check_in_AF]))
#    return(x, h, u_appx, u_appx_next, u_exact)
#    
#N=196    
#x, h, uno,dos,tres=CAF(N,lamda,U_o,A,k, L) 
#u_appx[check_val], u_appx_next[2*check_val], u_exact[check_val]
    
#This is Now friday 1.
    
#it seems like doubling the gridpoints isnt getting my answer closer to th exact value
#yeah I think something is actually wrong here with HTAF or the u_appx_next call
# ahh i think it is the h/2 thing. that was actually wrong but im still not getting exactly what I was expecting

#rel_error 2 should be smaller than rel_error 1
#okay I kinda see it now. u_appx_next[0] is closer to U_o, then both u_appx, and u_exact. but it still isnt
#u_appx[0]=u_appx[x= 0+1*h1]
#and u_appx_next[0]=u_appx_next[x= 0+1*h2]
# h2!=h1/2
#so this is the same issue as the other day. Im not actually checking the same physical point. I thought i had solved this problem
#in the grid. conv. study. but no it is still there.
# lets just do this for now. call the DIF and make 2 seperate x lists and compare the values
#ok so in the grid conv study I cant actually check the same physical point.
#so I will need to call call the DIF again for 2*N to do the relative error part.
#yeah that seems about right.    
#and then that fooa will need to be changed because I will have lists of different lengths.
# I can check elements close to one another, or maybe I should just take the average of the whole list.
#this seems like a good time for a break. and to add to github.