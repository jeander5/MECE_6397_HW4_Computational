# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 14:07:44 2020
@author: johna
"""
#https://github.com/jeander5/MECE_6397_HW4_Computational
#Okay I am making a new file and will be working from this version from now on using 
#Github to document all changes.

# MECE 6397, SciComp, Problem 4, Computational
#Solve the Helmholtzs equation for 1. Dirchlet. 2. Nuemann

#imports
import numpy as np
import matplotlib.pyplot as plt
from math import sinh as sinh
from math import cosh as cosh
from math import log2 as log2
import time

#Contstants given in problem statement, constant for both boundary conditions.
#Interval length, u(x=0) for Dirchlet, v is constant for the Neumann, A is exact solution of f(x).
L = 1
U_o = 1
v = 1
A = 1
K = [1, 10]

#The N_initial value must be greater than 2. 
#N is defined in the for lor and changes during the grid convergence study.
N_initial = 10

#Helmholtz Thomas Algorith Function, for Dirchlet part 1
def HTAF(N, h, lamda, U_o, A):   
#inputs are N, lamda, U_o=u(x=0), and for this problem A.
#Pre Thomas algorith set up. for this problem these values are all constant
    a = -(2-lamda*h**2)
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
    
#Helmholtz Thomas Algorith Function, for Neumann part 2
def NHTAF(N, h, lamda, v, A):   
#Pre Thomas algorith set up.
#I now need to make N one point larger to incorporate the ghost node method for
#u(x=0) whihch is unknown
#but I am still keeping h the same    
    N = N+1    
#these values are constant but the c's are not.    
    a = -(2-lamda*h**2)
    b = 1
#I now need c to be a list because they are now not all the same
#or I could use some conditinal statements but I wanna still follow the
#pseudo code for the algorithm closely    
    c = [1]*N
    c[0] = 2
    f= A*h**2    
#this line is added because of the ghost node method
#right hand side, initial equation    
    rhs_o = A*h**2+2*h*v
    alpha = [0]*N
    g = [0]*N
    u_appx = [0]*N
#following the psuedo code
#zeroth element of this list does infact correspond to  subscript zero in thomas algorith
#because of the ghost node method
    alpha[0] = a
    g[0] = rhs_o
    for j in range(1, N):
        alpha[j] = a-(b/alpha[j-1])*c[j-1]
        g[j] = f-(b/alpha[j-1])*g[j-1]   
    u_appx[N-1] = g[N-1]/alpha[N-1]
    for m in range(1,N):
        u_appx[-1-m] = (g[-1-m]-c[-1-m]*u_appx[-m])/alpha[-1-m]
    return(u_appx)    

#u exact function, for the Helmhotlz Dirchlet part 1 problem     
def uEF(k, L, x, A, U_o):
    u_exact = [((sinh(k*(L-x))+sinh(k*x))/sinh(k*L)-1)*A/k**2+
U_o*sinh(k*(L-x))/sinh(k*L) for x in x[1:-1]]
#x[1:-1] I dont need u(x=0) or u(x=L) beacause they are given
    return(u_exact)
  
#u exact function, for the Helmhotlz Neumann Part 2 problem
def uEF2(k, L, x, A, v):
    u2_exact = [((cosh(k*x)/cosh(k*L))-1)*(A/k**2)
-(v/k)*(sinh(k*(L-x))/cosh(k*L)) for x in x[0:-1]]
#x[0:-1]  I dont need  u(x=L) beacause it is given    
    return(u2_exact)    

#Pre-loop Plot formatting
#turning on grid, and setting up subplots, title, and labels.
#legends are defined inside the loop     
plt.rcParams['axes.grid'] = True
fig1, (ax1, ax2) = plt.subplots(2)
fig1.suptitle('Part 1 Dirchlet Boundary Conditions')
fig2, (ax3, ax4) = plt.subplots(2)  
fig2.suptitle('Part 2 Nuemann Boundary Conditions')
ax1.title.set_text('k = %s'%(K[0]))
ax1.set_xlabel('x')
ax1.set_ylabel('u(x)')
ax2.title.set_text('k = %s'%(K[1]))
ax2.set_xlabel('x')
ax2.set_ylabel('u(x)')
ax3.title.set_text('k = %s'%(K[0]))
ax3.set_xlabel('x')
ax3.set_ylabel('u(x)')
ax4.title.set_text('k = %s'%(K[1]))
ax4.set_xlabel('x')
ax4.set_ylabel('u(x)')

#changing this number right here changes my results!
#Note! Very Important!
#Difference between N and 2N variable
Diff_N2N = 2*10**-6
#1*10**-6 takes too long but 2*10**-6 works
#Outter for loop for the different k values
lenK = len(K)
for n in range(lenK):
    k = K[n]
    print('\nFor k = %s \n'%(k))
    
#using a different value here for k=10
    if n == 1:
        Diff_N2N = 1*10**-8
    
#Note: lamda in the Helmholtz eq defined here
    lamda = -k**2
    
#Part 1, Dirchlet
    print('Part 1. Dirchlet Boundary Conditions \n')
    
#Grid Convergence
    
#resetting the N value, im putting it here for now, needs to come before both grid convergenve studies
    N = N_initial    
#im using the flag so I dont have to call the function before and inside the while statement
    Flag = 0 
    while Flag == 0:
#comparing values near the middle of the interval now
        check_val = round(N/2)
#calling my functions
#changed this here so no I am only storing the xs I need to plot and call the exact function        
        h = L/(N+1) 
        u_appx = HTAF(N,h,lamda,U_o,A)
#I bet I should define these so they are not in the function call.
        N2 = 2*N
        h2 = L/(N2+1)
        u_appx_next = HTAF(N2,h2,lamda,U_o,A)
# I still need to be comparing u values for the closest x points.  
        if abs(u_appx[check_val]-u_appx_next[2*check_val+1])<Diff_N2N:
            Flag = 1
            print('%s Grid Points Needed' %(N))
            print('Doubling the Grid Points would result in less than %s differnce between u values for the closest Grid Point \n' %(Diff_N2N))
        else:
            N = N+N
            
#ls for legend string        
        ls1=('Approximate Value with %s gridpoints'%(N))
        ls2=('Apprixmate Value with %s gridpoints'%(N2))  
         
#Note: here are the exact value function calls
#I dont even really need these h's now    
    x=np.linspace(0, L, N+2)
    u_exact = uEF(k, L, x, A, U_o)
    x2=np.linspace(0, L, N2+2)
    u_exact_next = uEF(k, L, x2, A, U_o)   
    
    
#formal order of accuracy.
    
    
#Method One    
#just keeping this quick method for now. G is still just a generic varibale name
    start_time1=time.time()
    print('Start')
    G1 = max([abs((u_appx[j]-u_exact[j])) for j in range(N)]) 
    G2 = max([abs((u_appx_next[j]-u_exact_next[j])) for j in range(2*N)])  
    fooa = round(log2(G1/G2),2)
    print('log2(Error_N/Error_2N) = %s'%(fooa))
    print("--- %s seconds --" % (time.time()-start_time1))

#Method Two
    start_time2=time.time()
    print('Start')    
    stored_val = (u_appx[0]-u_exact[0])
    for j in range(1,N):
        next_val = u_appx[j]-u_exact[j]
        if next_val>stored_val:
            stored_val = next_val
    stored_val2 = (u_appx[0]-u_exact[0])
    for j in range(1,N*2):
        next_val2 = u_appx_next[j]-u_exact_next[j]
        if next_val2>stored_val2:
            stored_val2=next_val2
    print("--- %s seconds --" % (time.time()-start_time2))     
            
#Method 3
    start_time3=time.time()
    print('Start')
    print("--- %s seconds --" % (time.time()-start_time3))   
    abs_err_1 = [0]*N
    abs_err_2 = [0]*N*2
    for j in range(N):
        abs_err_1[j] = abs((u_appx[j]-u_exact[j]))
    for j in range(2*N):
        abs_err_2[j] = abs((u_appx_next[j]-u_exact_next[j]))
    max_err_1 = max(abs_err_1)
    max_err_2 = max(abs_err_2)            
    print("--- %s seconds --" % (time.time()-start_time3))
#part two, Nuemann
    print('\nPart 2. Neumann Boundary conditions \n')
    #grid convergenvce
    N = N_initial
    Flag = 0
    while Flag == 0:
#comparing values near the middle of the interval now
        check_val = round(N/2)
#calling my functions
        h3= L/(N+1)  
        u2_appx=NHTAF(N, h3, lamda, v, A)
        #I bet I should define these so they are not in the function call.
        N2 = 2*N
        h4 = L/(N2+1)
        u2_appx_next = NHTAF(N2, h4, lamda, v, A)
        # I still need to be comparing u values for the closest x points. 
        if abs(u2_appx[check_val]-u2_appx_next[2*check_val+1])<Diff_N2N:
            Flag = 1
            print('%s Grid Points Needed' %(N))
            print('Doubling the Grid Points would result in less than %s differnce between u values for the closest Grid Point \n' %(Diff_N2N))
        else:
            N=N+N   
    
#Note: here are the exact value function calls
    x3=np.linspace(0, L, N+2)        
    u2_exact = uEF2(k, L, x3, A, v)
    x4=np.linspace(0, L, N2+2)
    u2_exact_next = uEF2(k, L, x4, A, U_o)   
    
#formal order of accuracy.
#just keeping this quick method for now. G is still just a generic varibale name
    G3 = max([abs((u2_appx[j]-u2_exact[j])) for j in range(N+1)]) 
    G4 = max([abs((u2_appx_next[j]-u2_exact_next[j])) for j in range(2*N+1)])  
    fooa2 = round(log2(G3/G4),2)
    print('log2(Error_N/Error_2N) = %s'%(fooa2))

#plotting, inside the for loop still
    if n==1:
        ax1.plot(x2[1:-1], u_exact_next,'k')
        ax1.plot(x[1:-1], u_appx,'-.r')
        ax1.plot(x2[1:-1], u_appx_next,':b')
        ax1.legend(['Exact Value',ls1, ls2])
        ax3.plot(x4[0:-1], u2_exact_next,'k')
        ax3.plot(x3[0:-1], u2_appx,'-.r')
        ax3.plot(x4[0:-1], u2_appx_next,':b')
        ax3.legend(['Exact Value','Approximate Value with %s grid points'%(N),
'Apprixmate Value with %s grid points'%(N2)])
    else:
        ax2.plot(x2[1:-1], u_exact_next,'k')
        ax2.plot(x[1:-1], u_appx,'-.r')
        ax2.plot(x2[1:-1], u_appx_next,':b') 
        ax2.legend(['Exact Value',ls1,ls2])
        ax4.plot(x4[0:-1], u2_exact_next,'k')
        ax4.plot(x3[0:-1], u2_appx,'-.r')
        ax4.plot(x4[0:-1], u2_appx_next,':b')
        ax4.legend(['Exact Value','Approximate Value with %s grid points'%(N),
'Apprixmate Value with %s grid points'%(N2)])

#ok tables now
# I will just copy and paste what I need in excel
        
# ok nice it still takes a minute or so but I can get my uppx and my uappx next pretty dang close now.
#but the errorN<error2N is still happening
#you know it specifically says to check this AFTER the grid convergence so maybe this iant such a big deal.
#i still dont know why the 2N appx would be less acurate 
# also that searching the max error thing is taking awhile, maybe bring back those other methods 
#i mean its really not that long
#imagine subdividng a meter one micron at a time. crazy stuff       
#import time
#start_time=time.time()
#print('Start')
#print("--- %s seconds --" % (time.time()-start_time))
#I wanna time test those three different methods

