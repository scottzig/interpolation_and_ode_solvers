#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 5 Problem 3a

This program computes the orbit of a satellite around
the Earth.
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

# preliminaries

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
    

# Explicit Euler
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 500
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (Explicit Euler) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
    update = curr + delt_t*F_e(curr)
    curr = update
    t = t + delt_t
    
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 80
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (Explicit Euler) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
    update = curr + delt_t*F_e(curr)
    curr = update
    t = t + delt_t
    
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 10
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (Explicit Euler) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
    update = curr + delt_t*F_e(curr)
    curr = update
    t = t + delt_t
    
# Implicit Euler
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 500
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (Implicit Euler) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
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
    
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 80
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (Implicit Euler) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
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
    
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 10
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (Implicit Euler) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
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
    
# Crank-Nicolson
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 500
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (Crank-Nicolson) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
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
    
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 80
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (Crank-Nicolson) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
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
    
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 10
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (Crank-Nicolson) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
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
    
# RK-4
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 500
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (RK-4) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
    k_1 = delt_t*F_e(curr)
    k_2 = delt_t*F_e(curr+(1/2)*k_1)
    k_3 = delt_t*F_e(curr+(1/2)*k_2)
    k_4 = delt_t*F_e(curr+k_3)
    update = curr + (1/6)*(k_1+2*k_2+2*k_3+k_4)
    curr = update
    t = t + delt_t
    
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 80
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (RK-4) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
    k_1 = delt_t*F_e(curr)
    k_2 = delt_t*F_e(curr+(1/2)*k_1)
    k_3 = delt_t*F_e(curr+(1/2)*k_2)
    k_4 = delt_t*F_e(curr+k_3)
    update = curr + (1/6)*(k_1+2*k_2+2*k_3+k_4)
    curr = update
    t = t + delt_t
    
theta = np.linspace(0,2*np.pi,100)
r = 6.371*10**6 # radius of Earth in meters
x_1 = r*np.cos(theta)
x_2 = r*np.sin(theta)
fig, ax = plt.subplots(1)
ax.plot(x_1,x_2)
t = 0
delt_t = 10
u = x_0
v = v_0
curr = np.concatenate([v,u])
while t <= 6000:
    plt.plot(curr[2],curr[3],'ko')
    txt = str(delt_t)
    plt.title('Position of Satellite (RK-4) with delt_t = ' +txt)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
    k_1 = delt_t*F_e(curr)
    k_2 = delt_t*F_e(curr+(1/2)*k_1)
    k_3 = delt_t*F_e(curr+(1/2)*k_2)
    k_4 = delt_t*F_e(curr+k_3)
    update = curr + (1/6)*(k_1+2*k_2+2*k_3+k_4)
    curr = update
    t = t + delt_t