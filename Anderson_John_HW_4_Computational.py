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

#import math
import matplotlib.pyplot as plt
from math import sinh as sinh
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
    u_appx_next = HTAF(2*N,h/2,lamda,U_o,A)
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

#Next....Checking relative error, graphs and tables, then on to part two, 

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
    
#This is wednesday 5