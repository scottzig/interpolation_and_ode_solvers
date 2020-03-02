#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 5 Problem 2

This program computes solutions to an ODE with
a discontinuous right side using the explicit Euler
Method and RK-4.
"""
import numpy as np
import matplotlib.pyplot as plt

def x_true(t):
    if t<1/3:
        x_true = np.exp(1/3)
    else:
        x_true = (1+np.exp(-1.0/3))*np.exp(t)-1
    return x_true

# Explicit Euler
    
track = 0
err = np.zeros(10)
delt_t = np.zeros(10)
for ii in range(len(delt_t)):
    delt_t[ii] = float(1.0/2.0)**(ii+1)
for h, N in [(2**(-i), 2**i) for i in range(1,11)]:
    x_i = 1.0
    t_i = 0.0
    for i in range(N):
        #if t_i < (1.0 / 3.0):
        if i < (float(N) / 3):
            x_ip1 = ((1+h)*x_i)
            x_i = x_ip1
        else:
            x_ip1 = ((1+h)*x_i) + h
            x_i = x_ip1
        t_i += h
    err[track] = np.abs(x_true(1)-x_i)
    track = track + 1

plt.figure()
plt.loglog(delt_t,err,'-o', label = 'err')
plt.loglog(delt_t,delt_t,'-ro', label = 'h^1')
plt.title('Error in Explicit Euler Method')
plt.xlabel('h')
plt.ylabel('e(h)')
plt.legend()

##############################################################

track = 0
err = np.zeros(10)
delt_t = np.zeros(10)
for ii in range(len(delt_t)):
    delt_t[ii] = float(1.0/3.0)*float(1.0/2.0)**(ii)
for h, N in [((1.0/3.0)*2**(-i), 3*(2**i)) for i in range(0,10)]:
    x_i = 1.0
    t_i = 0.0
    for i in range(N):
        #if t_i < (1.0 / 3.0):
        if i < (N / 3):
            x_ip1 = ((1+h)*x_i)
            x_i = x_ip1
        else:
            x_ip1 = ((1+h)*x_i) + h
            x_i = x_ip1
        t_i += h
    err[track] = np.abs(x_true(1)-x_i)
    track = track + 1

plt.figure()
plt.loglog(delt_t,err,'-o', label = 'err')
plt.loglog(delt_t,delt_t,'-ro', label = 'h^1')
plt.title('Error in Explicit Euler Method (land on 1/3)')
plt.xlabel('h')
plt.ylabel('e(h)')
plt.legend()

# RK-4

track = 0
err = np.zeros(10)
delt_t = np.zeros(10)
for ii in range(len(delt_t)):
    delt_t[ii] = float(1.0/2.0)**(ii+1)
for h, N in [(2**(-i), 2**i) for i in range(1,11)]:
    x_i = 1.0
    t_i = 0.0
    for i in range(N):
        if i < (float(N) / 3):
            F_1 = (h * x_i)
            F_2 = (h * x_i) + ((h*F_1)/2.0)
            F_3 = (h * x_i) + ((h*F_2)/2.0)
            F_4 = (h * x_i) + (h*F_3)
            x_ip1 = x_i + ((1.0/6.0)*F_1) + ((2.0/6.0)*F_2) + ((2.0/6.0)*F_3) + ((1.0/6.0)*F_4)
            x_i = x_ip1
        else:
            F_1 = (h * x_i) + h
            F_2 = (h * x_i) + ((h*F_1)/2.0) + h
            F_3 = (h * x_i) + ((h*F_2)/2.0) + h
            F_4 = (h * x_i) + (h*F_3) + h
            x_ip1 = x_i + ((1.0/6.0)*F_1) + ((2.0/6.0)*F_2) + ((2.0/6.0)*F_3) + ((1.0/6.0)*F_4)
            x_i = x_ip1
        t_i += h
    err[track] = np.abs(x_true(1)-x_i)
    track = track + 1

plt.figure()
plt.loglog(delt_t,err,'-o', label = 'err')
plt.loglog(delt_t,delt_t,'-ro', label = 'h^1')
plt.title('Error in RK-4')
plt.xlabel('h')
plt.ylabel('e(h)')
plt.legend()

####################################################################

track = 0
err = np.zeros(10)
delt_t = np.zeros(10)
for ii in range(len(delt_t)):
    delt_t[ii] = float(1.0/3.0)*float(1.0/2.0)**(ii)
for h, N in [((1.0/3.0)*2**(-i), 3*(2**i)) for i in range(0,10)]:
    x_i = 1.0
    t_i = 0.0
    for i in range(N):
        if i < (N / 3):
            F_1 = (h * x_i)
            F_2 = (h * x_i) + ((h*F_1)/2.0)
            F_3 = (h * x_i) + ((h*F_2)/2.0)
            F_4 = (h * x_i) + (h*F_3)
            x_ip1 = x_i + ((1.0/6.0)*F_1) + ((2.0/6.0)*F_2) + ((2.0/6.0)*F_3) + ((1.0/6.0)*F_4)
            x_i = x_ip1
        else:
            F_1 = (h * x_i) + h
            F_2 = (h * x_i) + ((h*F_1)/2.0) + h
            F_3 = (h * x_i) + ((h*F_2)/2.0) + h
            F_4 = (h * x_i) + (h*F_3) + h
            x_ip1 = x_i + ((1.0/6.0)*F_1) + ((2.0/6.0)*F_2) + ((2.0/6.0)*F_3) + ((1.0/6.0)*F_4)
            x_i = x_ip1
        t_i += h
    err[track] = np.abs(x_true(1)-x_i)
    track = track + 1

test = np.zeros(10)
for ii in range(len(test)):
    test[ii] = delt_t[ii]**4
plt.figure()
plt.loglog(delt_t,err,'-o', label = 'err')
plt.loglog(delt_t,test,'-ro', label = 'h^4')
plt.title('Error in RK-4 (land on 1/3)')
plt.xlabel('h')
plt.ylabel('e(h)')
plt.legend()
    