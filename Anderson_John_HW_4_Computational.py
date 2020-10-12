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
import matplotlib.pyplot as plt
from math import sinh as sinh
from math import log2 as log2
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
N = 17
N_initial = N

#helmotlz dirchlet part 1 problem

#discretize the interval function  
#im gonna go ahead and change so that x does NOT include the endpoints
#this really doesnt need to be a function. I still llike having it tho   
def DIF(L, N):
#Discretizing the interval length. This is the same for both problems
    h = L/(N+1)
    x = np.linspace(h, L-h, N)
    return(x[:],h)

#Helmholtz Thomas Algorith Function    
#I still need that A as an input for this problem 
def HTAF(N, h, lamda, U_o, A):   
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
    
#changing this beacuse I changed the DIF    
def uEF(k, L, x, A, U_o):
    u_exact = [((sinh(k*(L-x))+sinh(k*x))/sinh(k*L)-1)*A/k**2+
U_o*sinh(k*(L-x))/sinh(k*L) for x in x]
    return(u_exact)

#Note now with changes to the above len(x) == len(u_exact)
#I think I should have just done this from the beginning
#Note:Functions are now defined. Moving along

 
#Calling the Discr. the Interval right here for now     
x,h = DIF(L,N)   
 
#im gonna need this eventually    
#lenK=len(K)
#for n in range(lenK)
#    k = K[n]

#just gonna start with k=1 for now
k = K[0]

#Note: lamda in the helmholtz eq defined here
lamda = -k**2

#Part 1, Dirchlet

#Grid Convergence

#This section I think is where the main problem is.

#"...If the results of the two calculations are different, it means
#that N nodes are too few.
#In this case you should compare the results for 2N and 4N and so forth until the
#results of two successive discretizations are about the same;..."

#maybe I should be checking the whole list idk

#the issue here is for large N, max_err_2>max_err_1.
#why is that the case. should that be the case? I definately dont think it should.
#as N goes up h goes down. and h^2 can be much smaller than the difference value betwwen 
#Is that really the issue? maybe
#it also seems like checking values earlier in the list gives different values for fooa 
#I should ceck the middle of the list probably,.
#still, if max_err_2>max_err_1 i get a negative result from that log equation
#maybe im just reaching the limits of the accuracy of the approximation, even with h getting smaller

#im using the flag so I dont have to call the function before and inside the while statement
Flag = 0
#Lebron is a placeholder variable name
Lebron = 1*10**-4
#Checking a random discretized point, not every value.
#check_val = randi(1,N-2)
while Flag == 0:
    check_val = round(N/2)
#calling all my functions
    print(N)
    # i shouldnt have this function in the loop
    #im calculating alot xs i dont use
    x, h = DIF(L,N)  
    u_appx=HTAF(N,h,lamda,U_o,A)
    #I bet I should define these so they are not in the function call.
    N2=2*N
    h2=L/(N2+1)
    u_appx_next = HTAF(N2,h2,lamda,U_o,A)
    # I still need to be comparing u values for the closest x points.  
    if abs(u_appx[check_val]-u_appx_next[2*check_val+1])<Lebron:
        Flag = 1
        print('%s Grid Points Needed' %(N))
        print('Doubling the Grid Points would result in less then %s differnce between u values for the closest Grid Point' %(Lebron))
    else:
        N=N+N   
print(check_val)
        

#Note: here is the exact value function calls
        
u_exact = uEF(k, L, x, A, U_o)
x2, h2= DIF(L, N2)
u_exact_next = uEF(k, L, x2, A, U_o)   

#Next....Checking absolute error, graphs and tables, then on to part two,
 
#To do that you have to calculate the maximum absolute errors
#in N and 2N nodes using the analytical solution. Then use following formula to find the formal order of
#accuracy.
#So I actually just need max absolute error for both N and 2N

#formal order of accuracy.
#lets create the storage list ahead of time.
abs_err_1 = [0]*N
abs_err_2 = [0]*N*2
for j in range(N):
    abs_err_1[j] = abs((u_appx[j]-u_exact[j]))
for j in range(2*N):
    abs_err_2[j] = abs((u_appx_next[j]-u_exact_next[j]))
max_err_1 = max(abs_err_1)
max_err_2 = max(abs_err_2)

#or should I do something like this, all in one line. 
#This kinds seems like I dont actually need to store the whole error
#BIG G is just generic place holder name   
G1 = max([abs((u_appx[j]-u_exact[j])) for j in range(N)]) 
G2 = max([abs((u_appx_next[j]-u_exact_next[j])) for j in range(2*N)])  

#OR should I check to see if the next element is greater than previous one 
#and only keep two values at a time?
#idk
stored_val = (u_appx[0]-u_exact[0])
for j in range(0,N):
    next_val = u_appx[j]-u_exact[j]
    if next_val>stored_val:
        stored_val = next_val
    
stored_val2 = (u_appx[0]-u_exact[0])
for j in range(0,N*2):
    next_val2 = u_appx_next[j]-u_exact_next[j]
    if next_val2>stored_val2:
        stored_val2=next_val2   
#all three should be the same number yep
#stored_val==G1==max_err_1
#stored_val2==G2==max_err_2
# I will just keep all three I like the middle cause its all in one line 
#but probably the 3rd  one is more computer efficient.
#I could time it real quick. Later        
#fooa= formal order of accuracy... for now.
fooa = log2(max_err_1/max_err_2)
print(fooa)

#ok so this kinda is working, still if my Initial N is too big Im not getting consistently 2.
#I think the grid convergence stuff still isnt right. Maybe I do need to be checking the whole list.
#okay sometimes Im getting a negative value, which means max_err_2>max_err_1
#which shouldnt be the case   
# maybe i should use relative error
#maybe i should actually shift the list so I am checking the same exact same x point for N and 2N somehow
#I think that will mess up the results from the algorithm
#i think im just stuck with checking closest valus


#I really dont need this havent used
#im gonna keep these lines here.
#the stuff below here is just for copying and pasting into the command line for troubleshooting    
#call all functions
#def CAF(N,lamda,U_o,A,k, L):
#    x, h = DIF(L,N) 
#    u_appx=HTAF(N,h,lamda,U_o,A)
#    u_appx_next = HTAF(2*N,h/2,lamda,U_o,A)        
#    u_exact = uEF(k, L, x, A, U_o)
#    x2, h2 = DIF(L,2*N)
#    u_exact_next = uEF(k, L, x2, A, U_o)     
#    check_in_AF=randi(1,N-1)
#    print('%s Grid Points %(N)')
#    print('u appx=%s'%(u_appx[check_in_AF]))
#    print('u appx_next=%s'%(u_appx_next[2*check_in_AF+1]))
#    print('u exact=%s'%(u_exact[check_in_AF]))
#    return(x, h, u_appx, u_appx_next, u_exact, u-exact_next)
#    
#N=196    
#x, h, uno,dos,tres,quatro=CAF(N,lamda,U_o,A,k, L) 
#u_appx[check_val], u_appx_next[2*check_val+1], u_exact[check_val], u_exact_next[check_val*2+1],    
#This is Now Mon 1.
