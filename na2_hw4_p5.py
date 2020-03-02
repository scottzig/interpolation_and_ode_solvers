#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 4 Problem 5

This code computes the altitude of a rocket by solving
a coupled system of ODEs using the explicit Euler Method.
"""
import numpy as np
import matplotlib.pyplot as plt

def T(t):
    if t < 600:
        return 12
    else:
        return 0
    
def m(t):
    if t < 600:
        return float(1-(0.9*t)/600)
    else:
        return float(0.1)
    
t = 0
delt_t = 60
u = 6371000
u_last = 0
v = 0
while t <= 36000:
    plt.semilogy(t,u-6371000,'ko')
    txt = str(delt_t)
    plt.title('Altitude of Rocket with delt_t = ' +txt)
    plt.xlabel('Time (s)')
    plt.ylabel('Altitude (m)')
    v_new = v + delt_t*(T(t)/m(t)-(10*6371000**2)/float(u**2))
    u_new = u + delt_t*v
    u = u_new
    v = v_new
    t = t + delt_t

plt.figure()
t = 0
delt_t = 35
u = 6371000
u_last = 0
v = 0
while t <= 36000:
    plt.semilogy(t,u-6371000,'ko')
    txt = str(delt_t)
    plt.title('Altitude of Rocket with delt_t = ' +txt)
    plt.xlabel('Time (s)')
    plt.ylabel('Altitude (m)')
    v_new = v + delt_t*(T(t)/m(t)-(10*6371000**2)/float(u**2))
    u_new = u + delt_t*v
    u = u_new
    v = v_new
    t = t + delt_t
    
plt.figure()    
t = 0
delt_t = 20
u = 6371000
u_last = 0
v = 0
while t <= 36000:
    plt.semilogy(t,u-6371000,'ko')
    txt = str(delt_t)
    plt.title('Altitude of Rocket with delt_t = ' +txt)
    plt.xlabel('Time (s)')
    plt.ylabel('Altitude (m)')
    v_new = v + delt_t*(T(t)/m(t)-(10*6371000**2)/float(u**2))
    u_new = u + delt_t*v
    u = u_new
    v = v_new
    t = t + delt_t
        
#delt_t = 64
#u_last = 0
#err = 10000
#while err > 1000:
#    print(delt_t)
#    t = 0
#    u = 6371000
#    v = 0
#    while t <= 36000:
#        plt.semilogy(t,u-6371000,'ko')
#        txt = str(delt_t)
#        plt.title('Altitude of Rocket with delt_t = ' +txt)
#        plt.xlabel('Time (s)')
#        plt.ylabel('Altitude (m)')
#        v_new = v + delt_t*(T(t)/m(t)-(10*6371000**2)/float(u**2))
#        u_new = u + delt_t*v
#        u = u_new
#        v = v_new
#        t = t + delt_t
#    err = np.abs(u-u_last)
#    print(err)
#    u_last = u
#    delt_t = delt_t/2