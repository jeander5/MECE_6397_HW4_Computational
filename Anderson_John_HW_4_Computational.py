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
from math import sinh as sinh
from math import cosh as cosh
from math import log2 as log2
import numpy as np
import matplotlib.pyplot as plt

#Contstants given in problem statement, constant for both boundary conditions.
#Interval length, u(x=0) for Dirchlet, v is constant for the Neumann, A is exact solution of f(x).
L = 1
U_0 = 1
v = 1
A = 1
K = [1, 10]

#The N_INITIAL value must be greater than 2.
#N is defined in the for lor and changes during the grid convergence study.
N_INITIAL = 10

#Helmholtz Thomas Algorithm Function, for Dirchlet part 1
def thomas_alg_one(N, h, lamda, U_0, A):
    """returns exact u values for the function from Part 1"""
#inputs are N, lamda, U_0=u(x=0), and for this problem A.
#Pre Thomas algorith set up. for this problem these values are all constant
    a = -(2-lamda*h**2)
    b = 1
    c = 1
#right hand side, the f*h**2 from the psuedocode
    rhs = A*h**2
    alpha = [0]*N
    g = [0]*N
    u_appx = [0]*N
#following the psuedo code
#zeroth element of this list corresponds to the first subscript in thomas algorith
    alpha[0] = a
    g[0] = rhs-U_0
    for j in range(1, N):
        alpha[j] = a-(b/alpha[j-1])*c
        g[j] = rhs-(b/alpha[j-1])*g[j-1]
    u_appx[N-1] = g[N-1]/alpha[N-1]
    for j in range(1, N):
        u_appx[-1-j] = (g[-1-j]-c*u_appx[-j])/alpha[-1-j]
    return u_appx

#Helmholtz Thomas Algorithm Function, for Neumann part 2
def thomas_alg_two(N, h, lamda, v, A):
    """returns approximate u values for the function from Part 2"""    
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
#this line is added because of the ghost node method
#right hand side, initial equation
    rhs = A*h**2
    alpha = [0]*N
    g = [0]*N
    u_appx = [0]*N
#following the psuedo code
#zeroth element of this list does infact correspond to  subscript zero in thomas algorith
#because of the ghost node method
    alpha[0] = a
    g[0] = rhs+2*h*v
    for j in range(1, N):
        alpha[j] = a-(b/alpha[j-1])*c[j-1]
        g[j] = rhs-(b/alpha[j-1])*g[j-1]
    u_appx[N-1] = g[N-1]/alpha[N-1]
    for j in range(1, N):
        u_appx[-1-j] = (g[-1-j]-c[-1-j]*u_appx[-j])/alpha[-1-j]
    return u_appx

#u exact function, for the Helmhotlz Dirchlet part 1 problem
def u_exact_func(k, L, x, A, U_0):
    """returns exact u values for the function from Part1"""
    func_vals = [((sinh(k*(L-x))+sinh(k*x))/sinh(k*L)-1)*A/k**2+
                 U_0*sinh(k*(L-x))/sinh(k*L) for x in x[1:-1]]
#x[1:-1] I dont need u(x=0) or u(x=L) beacause they are given
    return func_vals

#u exact function, for the Helmhotlz Neumann Part 2 problem
def u_exact_func2(k, L, x, A, v):
    """returns exact u values for the function from Part 2"""
    func_vals = [((cosh(k*x)/cosh(k*L))-1)*(A/k**2)
                 -(v/k)*(sinh(k*(L-x))/cosh(k*L)) for x in x[0:-1]]
#x[0:-1]  I dont need  u(x=L) beacause it is given
    return func_vals

#Pre-loop Plot formatting
#turning on grid, and setting up subplots, title, and labels.
#legends are defined inside the loop
plt.rcParams['axes.grid'] = True
fig1, (ax1, ax2) = plt.subplots(2)
fig1.suptitle('Part 1 Dirchlet Boundary Conditions')
fig2, (ax3, ax4) = plt.subplots(2)
fig2.suptitle('Part 2 Neumann Boundary Conditions')
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
CLOSE_ENOUGH = 6.33*10**-4
#This is really the nu,ber the controls the grid convergence
#as in "how close is close enough?"
#thats probably what should have named that variable
#CLOSE_ENOUGH
#1*10**-6 takes too long but 2*10**-6 works
#Outter for loop for the different k values
LEN_K = len(K)
for n in range(LEN_K):
    k = K[n]
    print('\nFor k = %s \n'%(k))

#using a different value here for k=10
    if n == 1:
        CLOSE_ENOUGH = 1.65*10**-6

#Note: lamda in the Helmholtz eq defined here
    lamda = -k**2

#Part 1, Dirchlet
    print('Part 1. Dirchlet Boundary Conditions \n')

#Grid Convergence

#resetting the N value, im putting it here for now
#needs to come before both grid convergenve studies
    N = N_INITIAL
#im using the Flag so I dont have to call the function before and inside the while statement
    Flag = 0
    while Flag == 0:
#comparing values near the middle of the interval now
        check_val = round(N/2)
#calling my functions
#changed this here so no I am only storing the xs I need to plot and call the exact function
        h = L/(N+1)
        u_appx = thomas_alg_one(N, h, lamda, U_0, A)
#I bet I should define these so they are not in the function call.
        N2 = 2*N
        h2 = L/(N2+1)
        u_appx_next = thomas_alg_one(N2, h2, lamda, U_0, A)
# I still need to be comparing u values for the closest x points.
        if abs(u_appx[check_val]-u_appx_next[2*check_val+1]) < CLOSE_ENOUGH:
            Flag = 1
            print('%s Grid Points Needed' %(N))
            print('Doubling the Grid Points would result in less than %s '
                  'differnce between u values for the closest grid points \n'
                  %(CLOSE_ENOUGH))
        else:
            N = N+N

#ls for legend string
        ls1 = ('Approximate Value with %s grid points'%(N))
        ls2 = ('Apprixmate Value with %s grid points'%(N2))

#Note: here are the exact value function calls
#I dont even really need these h's now
    x1 = np.linspace(0, L, N+2)
    u_exact = u_exact_func(k, L, x1, A, U_0)
    x2 = np.linspace(0, L, N2+2)
    u_exact_next = u_exact_func(k, L, x2, A, U_0)


#formal order of accuracy. fooa

#finding max absolute error
#keeping this method, after timing this proved to be the fastest
#for shorter lists tho the first method was faster     
    stored_val = (u_appx[0]-u_exact[0])
    for j in range(1, N):
        next_val = u_appx[j]-u_exact[j]
        if next_val > stored_val:
            stored_val = next_val
    stored_val2 = (u_appx[0]-u_exact[0])
    for j in range(1, N*2):
        next_val2 = u_appx_next[j]-u_exact_next[j]
        if next_val2 > stored_val2:
            stored_val2 = next_val2

    fooa1 = round(log2(stored_val/stored_val2), 2)
    print('Log2(Error_N/Error_2N) = %s'%(fooa1))

#part two, Nuemann
    print('\nPart 2. Neumann Boundary Conditions \n')
    #grid convergenvce
    N = N_INITIAL
    Flag = 0
    while Flag == 0:
#comparing values near the middle of the interval now
        check_val = round(N/2)
#calling my functions
        h3 = L/(N+1)
        u2_appx = thomas_alg_two(N, h3, lamda, v, A)
        #I bet I should define these so they are not in the function call.
        N2 = 2*N
        h4 = L/(N2+1)
        u2_appx_next = thomas_alg_two(N2, h4, lamda, v, A)
        # I still need to be comparing u values for the closest x points.
        if abs(u2_appx[check_val]-u2_appx_next[2*check_val+1]) < CLOSE_ENOUGH:
            Flag = 1
            print('%s Grid Points Needed' %(N))
            print('Doubling the Grid Points would result in less than %s '
                  'differnce between u values for the closest grid points \n'
                  %(CLOSE_ENOUGH))
        else:
            N = N+N

#Note: here are the exact value function calls
    x3 = np.linspace(0, L, N+2)
    u2_exact = u_exact_func2(k, L, x3, A, v)
    x4 = np.linspace(0, L, N2+2)
    u2_exact_next = u_exact_func2(k, L, x4, A, U_0)

#formal order of accuracy.
#just keeping this quick method for now. stored val and next val
# i could just redifine these numbers but its nice being to call them
    stored_val3 = (u2_appx[0]-u2_exact[0])
    for j in range(1, N):
        next_val3 = u2_appx[j]-u2_exact[j]
        if next_val3 > stored_val3:
            stored_val3 = next_val3
    stored_val4 = (u2_appx_next[0]-u2_exact_next[0])
    for j in range(1, N*2):
        next_val4 = u2_appx_next[j]-u2_exact_next[j]
        if next_val4 > stored_val4:
            stored_val4 = next_val4

    fooa2 = round(log2(stored_val3/stored_val4), 2)
    print('Log2(Error_N/Error_2N) = %s'%(fooa2))

#plotting, inside the for loop still
    if n == 1:
        ax1.plot(x2[1:-1], u_exact_next, 'k')
        ax1.plot(x1[1:-1], u_appx, '-.r')
        ax1.plot(x2[1:-1], u_appx_next, ':b')
        ax1.legend(['Exact Value', ls1, ls2])
        ax3.plot(x4[0:-1], u2_exact_next, 'k')
        ax3.plot(x3[0:-1], u2_appx, '-.r')
        ax3.plot(x4[0:-1], u2_appx_next, ':b')
        ax3.legend(['Exact Value', 'Approximate Value with %s grid points'%(N),
                    'Apprixmate Value with %s grid points'%(N2)])
    else:
        ax2.plot(x2[1:-1], u_exact_next, 'k')
        ax2.plot(x1[1:-1], u_appx, '-.r')
        ax2.plot(x2[1:-1], u_appx_next, ':b')
        ax2.legend(['Exact Value', ls1, ls2])
        ax4.plot(x4[0:-1], u2_exact_next, 'k')
        ax4.plot(x3[0:-1], u2_appx, '-.r')
        ax4.plot(x4[0:-1], u2_appx_next, ':b')
        ax4.legend(['Exact Value', 'Approximate Value with %s grid points'%(N),
                    'Apprixmate Value with %s grid points'%(N2)])

#ok tables now
# I will just copy and paste what I need in excel        
# ok nice it still takes a minute or so
#but I can get my uppx and my uappx next pretty dang close now.
#but the errorN<error2N is still happening
#you know it specifically says to check this AFTER
# the grid convergence so maybe this iant such a big deal.
#i still dont know why the 2N appx would be less acurate 
# also that searching the max error thing is taking awhile, maybe bring back those other methods 
#i mean its really not that long   
#import time
#start_time=time.time()
#print('Start')
#print("--- %s seconds --" % (time.time()-start_time))
#I wanna time test those three different methods
#looks like method two is the fastest. juts like I predicted nice
#???
#you know what? is it faster to make an 'a', 'b', and 'c' list outside the function in the for loop
#and have that be the input to the function then assign and call each a value.??
#I bet it woud be. It feel odd to me to make a whole a list with just same value repeated.
# but this cuts way down on the actual caluation time.
# no not really, no not all I think. It has to be calcaulate either before or
#inside the function, and it is outside of the for loops in the function.
#so nevermind
