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
N=13

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

#okay im gonna change this to be a Helmhotz TAF, with lamda.
#I still need that A as an input for this problem 
def HTAF(N,h,lamda,U_o,A):
#inputs are N, the tridiagonal elements, the RHS, and U_o=u(x=0).
#create lists ahead of time
    #Pre Thomas algorith set up. for this problem these values are all constant
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
#lamda in the helmoltz eq defined here
lamda=-k**2

#Part 1, Dirchlet

#grid convergence
#im using the flag so I dont have to call the function before and inside the while statement
#there is probably a better way to do this
#am i wasting memory using the flag? probably but its just one number
Flag = 0
m=1
#I need a variable name or something for this. I cant think of anything right now,
#closeness? diff, naw thats confusing, with differentiation. close enough
#wahtever 
Lebron=1*10**-3
#this is the point we will check. and what the hell if that ones good we will check another.  
check_val=randi(1,N-2)
while Flag == 0:
#calling all my functions
    print(m)
    print(N)
    x, h = DIF(L,N)  
    u_appx=HTAF(N,h,lamda,U_o,A)
    u_appx_next = HTAF(2*N,h/2,lamda,U_o,A)
    if abs(u_appx[check_val]-u_appx_next[2*check_val])<Lebron:
        print('checked once')
        check_val=randi(2,N-3)#I brought the interval in here just a little
        if abs(u_appx[check_val]-u_appx_next[2*check_val])<Lebron:
            print('checked twice')
        Flag = 1
        print('%s Grid Points Needed' %(N))
        print('Doubling the Grid Points would result in less then %s differece between u values for the same Grid Point' %(Lebron))
    else:
        N=N+N
        m=m+1
#here is the exact value function call
u_exact = uEF(k, L, x, A, U_o)     



#ok yes that was indeed "A" problem, but now the number of gridpoints is getting insane.
#lets use that "close enough"        
#ok good now the number of grid points needed depends on how accurate I want the numbers
#AND I really should probably be checking more than one number because im getting different results        
#im getting different values depending on which element I check.
#I really dont wanna check the whole list.
#hmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
# lets think about this?...?...?...?
#10 thousand grid points seems ok to me.
#but 100 thousand  seems insane.       
#okay do the earlier values need more than the middle? the middle more than the end.....?
#hmmmmmmmmmm        
#ok but hold up im not even checking the value for the same point here
#whatever I do I gotta be checking the same physical point and comparing the value.
#ah I see this now, very nice. I need to check u_appx[y]-u_appx_next[2*y])
#that was probably pretty simple
#oh well I figured it out pretty quick.       
## 
#im just gonna check a random value 
#       
#ok damn that still aint working quite right. something still isnt right
#ok I gotta divde the h by 2 as well, or call the DIF again. Yep.
#now the number of points needed is related to the differnce between the values needed/required/wanted/whatever
#1*10^-3----> appx 1000 grid points 1*10^-4 --->appx 10,000 grip points, 10^-5 appx 100,000 good, good   

#the stuff below here is just for copying and pasting into the command line for troubleshooting    
#N=100
#x, h = DIF(L,N)        
#        
#u_appx=HTAF(N,h,lamda,U_o,A)
#u_appx_next = HTAF(2*N,h/2,lamda,U_o,A)        
#u_exact = uEF(k, L, x, A, U_o)  

#call all functions,... CAF
#def AF(N,lamda,U_o,A,k, L):
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

#N=196    
#x, h, uno,dos,tres=AF(N,lamda,U_o,A,k, L) 
#u_appx[check_val], u_appx_next[2*check_val], u_exact[check_val]

#This is Wednsday 4    