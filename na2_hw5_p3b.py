#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 5 Problem 3b

This program computes the error in position of a satellite
after 100 orbits
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

# Extract time of one orbit

x_0 = np.array([0, 6790*10**3]) # meausred in m
v_0 = np.array([(27578*10**3)/(60**2),0]) # measured in m/sec

def F_e(vec):
    v = vec[0:2]
    u = vec[2:4]
    vec_1 = ((-(6371000**2)*9.81)/((u[0]**2+u[1]**2)**(3/2)))*u
    vec_2 = v
    return np.concatenate([vec_1,vec_2])

def F_i(vec,current,delt_t):
    v = vec[0:2]
    u = vec[2:4]
    vec_1 = ((-(6371000**2)*9.81)/((u[0]**2+u[1]**2)**(3/2)))*u
    vec_2 = v
    tog = np.concatenate([vec_1,vec_2])
    func = vec - current - delt_t*tog
    return func

def Jac_i(vec,current,delt_t):
    c = (-(6371000**2)*9.81)/((vec[2]**2+vec[3]**2)**(3/2))
    d = (vec[2]**2+vec[3]**2)**(3/2)
    partial_13 = -delt_t*(c+vec[2]*(-(6371000)**2*9.81)*(-3/2)*d**(-5/2)*2*vec[2])
    partial_14 = -delt_t*(vec[2]*(-(6371000)**2*9.81)*(-3/2)*d**(-5/2)*2*vec[3])
    partial_23 = -delt_t*(vec[3]*(-(6371000)**2*9.81)*(-3/2)*d**(-5/2)*2*vec[2])
    partial_24 = -delt_t*(c+vec[3]*(-(6371000)**2*9.81)*(-3/2)*d**(-5/2)*2*vec[3])
    mat = np.array([[1,0,partial_13,partial_14],
                    [0,1,partial_23,partial_24],
                    [-delt_t,0,1,0],
                    [0,-delt_t,0,1]])
    mat_inv = inv(mat)
    return mat_inv

t = 0
delt_t = 2
u = x_0
u_last = 0
v = v_0
curr = np.concatenate([v,u])
while t <= 5570:
    k_1 = delt_t*F_e(curr)
    k_2 = delt_t*F_e(curr+(1/2)*k_1)
    k_3 = delt_t*F_e(curr+(1/2)*k_2)
    k_4 = delt_t*F_e(curr+k_3)
    update = curr + (1/6)*(k_1+2*k_2+2*k_3+k_4)
    curr = update
    if (t > 100 and np.abs(x_0[1]-curr[3]) < 1000):
        t_orbit = t
        break
    t = t + delt_t
    
step_vec = np.array([t_orbit/10,t_orbit/100,t_orbit/1000])
#step_vec = np.array([t_orbit/10,t_orbit/100,t_orbit/1000, t_orbit/100000, t_orbit/1000000])

# Explicit Euler
print('Explicit Euler')
plt.figure()
err_vec = np.zeros(len(step_vec))
for ii in range(0,len(step_vec)):
    delt_t = step_vec[ii]
    t = 0
    u = x_0
    v = v_0
    curr = np.concatenate([v,u])
    while t <= 100*t_orbit:
        update = curr + delt_t*F_e(curr)
        curr = update
        t = t + delt_t
    test = np.array([curr[2],curr[3]])
    err_vec[ii] = np.linalg.norm(test-x_0)
    txt_1 = str(test)
    txt_2 = str(err_vec[ii])
    print('Position after 100 orbits is ' + txt_1)
    print('The error in this approximation is ' + txt_2)
    
plt.loglog(step_vec,err_vec)
plt.title('Error in Explicit Euler')
plt.xlabel('Time Step')
plt.ylabel('Error in approximation after 100 orbits')

# Implicit Euler
print('Implicit Euler')
plt.figure()
err_vec = np.zeros(len(step_vec))
for ii in range(0,len(step_vec)):
    delt_t = step_vec[ii]
    t = 0
    u = x_0
    v = v_0
    curr = np.concatenate([v,u])
    while t <= 100*t_orbit:
        start = curr + delt_t*F_e(curr) # run explicit to extract start
        err_n = 10
        while err_n > 10**(-1):
            vec_up = np.matmul(Jac_i(start,curr,delt_t),np.reshape(F_i(start,curr,delt_t),(4,1)))
            vec_up = np.reshape(vec_up,(1,4))
            update = start - vec_up
            err_n = np.linalg.norm(update[0]-start)
            start = update[0]
        curr = start
        t = t + delt_t
    test = np.array([curr[2],curr[3]])
    err_vec[ii] = np.linalg.norm(test-x_0)
    txt_1 = str(test)
    txt_2 = str(err_vec[ii])
    print('Position after 100 orbits is ' + txt_1)
    print('The error in this approximation is ' + txt_2)
    
plt.loglog(step_vec,err_vec)
plt.title('Error in Implicit Euler')
plt.xlabel('Time Step')
plt.ylabel('Error in approximation after 100 orbits')

# Crank-Nicolson
print('Crank-Nicolson')
plt.figure()
err_vec = np.zeros(len(step_vec))
for ii in range(0,len(step_vec)):
    delt_t = step_vec[ii]
    t = 0
    u = x_0
    v = v_0
    curr = np.concatenate([v,u])
    while t <= 100*t_orbit:
        start = curr + delt_t*F_e(curr) # run explicit to extract start
        err_n = 10
        while err_n > 10**(-1):
            F = (1/2)*(F_i(start,curr,delt_t)+F_e(start))
            vec_up = np.matmul(Jac_i(start,curr,delt_t),np.reshape(F,(4,1)))
            vec_up = np.reshape(vec_up,(1,4))
            update = start - vec_up
            err_n = np.linalg.norm(update[0]-start)
            start = update[0]
        curr = start
        t = t + delt_t
    test = np.array([curr[2],curr[3]])
    err_vec[ii] = np.linalg.norm(test-x_0)
    txt_1 = str(test)
    txt_2 = str(err_vec[ii])
    print('Position after 100 orbits is ' + txt_1)
    print('The error in this approximation is ' + txt_2)
    
plt.loglog(step_vec,err_vec)
plt.title('Error in Crank-Nicolson')
plt.xlabel('Time Step')
plt.ylabel('Error in approximation after 100 orbits')

# RK-4
print('RK-4')
plt.figure()
err_vec = np.zeros(len(step_vec))
for ii in range(0,len(step_vec)):
    delt_t = step_vec[ii]
    t = 0
    u = x_0
    v = v_0
    curr = np.concatenate([v,u])
    while t <= 100*t_orbit:
        k_1 = delt_t*F_e(curr)
        k_2 = delt_t*F_e(curr+(1/2)*k_1)
        k_3 = delt_t*F_e(curr+(1/2)*k_2)
        k_4 = delt_t*F_e(curr+k_3)
        update = curr + (1/6)*(k_1+2*k_2+2*k_3+k_4)
        curr = update
        t = t + delt_t
    test = np.array([curr[2],curr[3]])
    err_vec[ii] = np.linalg.norm(test-x_0)
    txt_1 = str(test)
    txt_2 = str(err_vec[ii])
    print('Position after 100 orbits is ' + txt_1)
    print('The error in this approximation is ' + txt_2)
    
plt.loglog(step_vec,err_vec)
plt.title('Error in RK-4')
plt.xlabel('Time Step')
plt.ylabel('Error in approximation after 100 orbits')