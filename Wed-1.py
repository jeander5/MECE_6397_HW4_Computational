# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:37:40 2020

@author: johna
"""
# MECE 6397, SciComp, Problem 4, Computational
#Solve the helmholts equation for 1. Dirchlet. 2. Nuemann

#imports
import numpy as np
#import math
import matplotlib.pyplot as plt
from math import sinh as sinh

#Contstants given in problem statement, constant for both boundary conditions.
#Interval length, u(x=0), v is constant on hyperbolic sin function, A is exact solution of f(x).
# solve both problem for both values of k
L = 1
U_o = 1
v = 1
A = 1
K = [1, 10]

#The N value, this is gonna change
N=10

#The code analyis is gonna hate these function names I bet
#Thomas Algorithm Function, note im not making a general thomas algorithm function, 
#just one for this homework assignments
#the a b c and f would be lists, cause they will not always be the same for every tri-di Matrix
#helmotlz dirchlet part 1 problem

#discretize the interval function     
def DIF(L,N):
#Discretizing the interval length. This is the same for both problems
    h = L/(N+1)
    x = np.linspace(0, L, N+2)
    return(x[:],h)

def TAF(N,a,b,c,f,U_o):
#inputs are N, the tridiagonal elements, the RHS, and U_o=u(x=0).
#create lists ahead of time
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

#u exact function, for the helmotlz dirchlet part 1 problem
def uEF(k,L,x,A,U_o):
    u_exact = [((sinh(k*(L-x))+sinh(k*x))/sinh(k*L)-1)*A/k**2+
U_o*sinh(k*(L-x))/sinh(k*L) for x in x[1:-1]]
    return(u_exact)
 
x,h=DIF(L,N)    

#im gonna need this eventually    
#lenK=len(K)
#for n in range(lenK)
#    k = K[n]
#just gonna start with k=1 for now
k = K[0]

#Part 1, Dirchlet
#Pre Thomas algorith set up. for this problem these values are all constant
a =-(2+k**2*h**2)
b = 1
c = 1
f = A*h**2

#grid convergence
#im using the flag so I dont have to call the function before and inside the while statement
Flag = 0
m=1
while Flag == 0:
#calling all my functions
    print(m)
    print(N)
    x, h = DIF(L,N)  
    u_appx=TAF(N, a, b, c, f, U_o)
    u_appx_next = TAF(2*N, a, b, c, f, U_o)
#Im only checking one value
#And I gotta have conditional statements somewhere     
    if u_appx[7] == u_appx_next[7]:
        Flag = 1
        print('%s Grid Points Needed' %(N))
    else:
        N=N+N
        m=m+1
        
##maybe something else is catastrophically wrong. yes that is why. ok I see it.
# my a and f values are changing but im not redifining them with new N.
#that stuff has gotta be in the function definition.
#Ill make the git hub later.

#This is Wednsday 3