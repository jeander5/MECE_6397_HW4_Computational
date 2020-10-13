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

#Contstants given in problem statement, constant for both boundary conditions
#Interval length, u(x=0) for Dirchlet, v is constant for the Neumann, A is exact solution of f(x)
L = 1
U_0 = 1
v = 1
A = 1
K = [1, 10]

#The N_INITIAL value must be greater than 2
#N is defined in the for loop and changes during the grid convergence study
N_INITIAL = 10

#Helmholtz Thomas Algorithm Function, for Dirchlet part 1
def thomas_alg_one(N, h, lamda, U_0, A):
    """returns exact u values for the function from Part 1"""
#inputs are N, lamda, U_0=u(x=0), and for this problem A.
#Pre Thomas algorithm set up. for this problem these values are all constant
    a = -(2-lamda*h**2)
    b = 1
    c = 1
#Right hand side, the f*h**2 from the psuedocode
    rhs = A*h**2
    alpha = [0]*N
    g = [0]*N
    u_appx = [0]*N
#Following the psuedocode
#Zeroth element of this list corresponds to the first subscript in thomas algorithm
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
#Pre Thomas algorithm set up.
#I now need to make N one point larger to incorporate the ghost node method for
#u(x=0) which is unknown. But I am still keeping h the same
    N = N+1
#These values are constant but the c's are not
    a = -(2-lamda*h**2)
    b = 1
#I now need c to be a list because they are now not all the same
#Or I could use some conditinal statements but I wanna still follow the
#Pseudocode for the algorithm closely
    c = [1]*N
    c[0] = 2
#This line is added because of the ghost node method
#Right hand side, initial equation
    rhs = A*h**2
    alpha = [0]*N
    g = [0]*N
    u_appx = [0]*N
#Following the psuedo code
#Zeroth element of this list does infact correspond to  subscript zero in thomas algorith
#Because of the ghost node method
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
#Turning on grid, and setting up subplots, titles, and labels.=
#Legends are defined inside the for loop
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

#Changing this number right here changes my results!
#Note! Very Important!
#Difference between N and 2N approximation
CLOSE_ENOUGH = 1*10**-3
#This is really the number the controls the grid convergence study
#As in "how close is close enough?"

#Outter for loop for the different k values
LEN_K = len(K)
for n in range(LEN_K):
    k = K[n]
    print('\nFor k = %s \n'%(k))

#Using a different value here for k=10
    if n == 1:
        CLOSE_ENOUGH = 1*10**-5

#Note: lamda in the Helmholtz equation is defined here
    lamda = -k**2

#Part 1, Dirchlet
    print('Part 1. Dirchlet Boundary Conditions \n')

#Grid Convergence Study

#Resetting the N value. Needs to come before both grid convergenve studies
    N = N_INITIAL
#I am using the Flag so I dont have to call the function before and inside the while statement
    Flag = 0
    while Flag == 0:
#Comparing values near the middle of the interval
        check_val = round(N/2)
#Calling the functions
        h = L/(N+1)
        u_appx = thomas_alg_one(N, h, lamda, U_0, A)
        N2 = 2*N
        h2 = L/(N2+1)
        u_appx_next = thomas_alg_one(N2, h2, lamda, U_0, A)
#I still need to be comparing u values for the closest x points, hence 2*check_val+1
        if abs(u_appx[check_val]-u_appx_next[2*check_val+1]) < CLOSE_ENOUGH:
            Flag = 1
            print('%s Grid Points Needed' %(N))
            print('Doubling the Grid Points would result in less than %s '
                  'differnce between u values for the closest grid points \n'
                  %(CLOSE_ENOUGH))
        else:
            N = N+N

#ls for legend string. This is here so I dont need to store this N and N2 value.
#I am storing the string instead
        ls1 = ('Approximate values with %s grid points'%(N))
        ls2 = ('Approximate values with %s grid points'%(N2))

#Note: here are the exact value function calls. x values are also defined here
    x1 = np.linspace(0, L, N+2)
    u_exact = u_exact_func(k, L, x1, A, U_0)
    x2 = np.linspace(0, L, N2+2)
    u_exact_next = u_exact_func(k, L, x2, A, U_0)

#formal order of accuracy, fooa

#Finding max absolute error, this method proved to be faster than others
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

#Part 2, Neumann
    print('\nPart 2. Neumann Boundary Conditions \n')

#Grid Convergence Study

    N = N_INITIAL
    Flag = 0
    while Flag == 0:
#Comparing values near the middle
        check_val = round(N/2)
#Calling the functions
        h3 = L/(N+1)
        u2_appx = thomas_alg_two(N, h3, lamda, v, A)
        N2 = 2*N
        h4 = L/(N2+1)
        u2_appx_next = thomas_alg_two(N2, h4, lamda, v, A)
# I still need to be comparing u values for the closest x points, hence 2*check_val+1
        if abs(u2_appx[check_val]-u2_appx_next[2*check_val+1]) < CLOSE_ENOUGH:
            Flag = 1
            print('%s Grid Points Needed' %(N))
            print('Doubling the Grid Points would result in less than %s '
                  'differnce between u values for the closest grid points \n'
                  %(CLOSE_ENOUGH))
        else:
            N = N+N

#Note: here are the exact value function calls. x values are also defined here
    x3 = np.linspace(0, L, N+2)
    u2_exact = u_exact_func2(k, L, x3, A, v)
    x4 = np.linspace(0, L, N2+2)
    u2_exact_next = u_exact_func2(k, L, x4, A, U_0)

#formal order of accuracy, fooa

#Finding max absolute error, this method proved to be faster than others

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

#Plotting, inside the for loop still
#k==1 plots
    if n == 1:
        ax1.plot(x2[1:-1], u_exact_next, 'k')
        ax1.plot(x1[1:-1], u_appx, '-.r')
        ax1.plot(x2[1:-1], u_appx_next, ':b')
        ax1.legend(['Exact values', ls1, ls2])
        ax3.plot(x4[0:-1], u2_exact_next, 'k')
        ax3.plot(x3[0:-1], u2_appx, '-.r')
        ax3.plot(x4[0:-1], u2_appx_next, ':b')
        ax3.legend(['Exact values', 'Approximate values with %s grid points'%(N),
                    'Approximate values with %s grid points'%(N2)])
#k==10 plots
    else:
        ax2.plot(x2[1:-1], u_exact_next, 'k')
        ax2.plot(x1[1:-1], u_appx, '-.r')
        ax2.plot(x2[1:-1], u_appx_next, ':b')
        ax2.legend(['Exact values', ls1, ls2])
        ax4.plot(x4[0:-1], u2_exact_next, 'k')
        ax4.plot(x3[0:-1], u2_appx, '-.r')
        ax4.plot(x4[0:-1], u2_appx_next, ':b')
        ax4.legend(['Exact values', 'Approximate values with %s grid points'%(N),
                    'Approximate values with %s grid points'%(N2)])
        