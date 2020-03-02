#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical Analysis II Homework 2 Problem 3b

This code uses a cubic spline interpolation to approximate
f(x) = xe^x-1 and then uses that approximation in Newton's
Method to find the zero of f(x).
"""
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x*np.exp(x)-1

def der_f(x):
    return x*np.exp(x)+np.exp(x)

def fill_matrix(i, x, y, mat, vec):
    for j in range(0,i-1):
        x[j] = j/(i-1)
        y[j] = x[j]*np.exp(x[j])-1
    mat[0,0] = 2/(x[1]-x[0])
    mat[0,1] = 1/(x[1]-x[0])
    vec[0] = 3*(y[1]-y[0])/(x[1]-x[0])**2
    mat[i-2,i-3] = 1/(x[i-2]-x[i-3])
    mat[i-2,i-2] = 2/(x[i-2]-x[i-3])
    vec[i-2] = 3*(y[i-2]-y[i-3])/(x[i-2]-x[i-3])**2
    for j in range(1,i-2):
        mat[j,j-1] = 1/(x[j]-x[j-1])
        mat[j,j] = 2*(1/(x[j]-x[j-1]))+2*(1/(x[j+1]-x[j]))
        mat[j,j+1] = 1/(x[j+1]-x[j])
        vec[j] = 3*((y[j]-y[j-1])/(x[j]-x[j-1])**2)+3*((y[j+1]-y[j])/(x[j+1]-x[j])**2)
    return np.linalg.solve(mat,vec)
        
def q(x,x_1,x_2,k_1,k_2):
    y_1 = x_1*np.exp(x_1)-1
    y_2 = x_2*np.exp(x_2)-1
    t = (x-x_1)/(x_2-x_1)
    a = k_1*(x_2-x_1)-(y_2-y_1)
    b= -k_2*(x_2-x_1)+(y_2-y_1)
    q_fcn = (1-t)*y_1+t*y_2+t*(1-t)*((1-t)*a+t*b)
    return q_fcn

def der_q(x,x_1,x_2,k_1,k_2):
    y_1 = x_1*np.exp(x_1)-1
    y_2 = x_2*np.exp(x_2)-1
    t = (x-x_1)/(x_2-x_1)  
    a = k_1*(x_2-x_1)-(y_2-y_1)
    b= -k_2*(x_2-x_1)+(y_2-y_1)
    d_q = (y_2-y_1)/(x_2-x_1)+(1-2*t)*(a*(1-t)+b*t)/(x_2-x_1)+ \
            t*(1-t)*(b-a)/(x_2-x_1)
    return d_q
    
N = 101

#for i in range(3,N):
#    text = str(i-1)
#    plt.figure()
#    oo = np.linspace(0,1,100)
#    plt.plot(oo,f(oo))
#    plt.title('Number of Interpolation Points = ' + text)
#    mat = np.zeros([i-1,i-1])
#    vec = np.zeros(i-1)
#    x = np.zeros(i-1)
#    y = np.zeros(i-1)
#    k_vec = fill_matrix(i,x,y,mat,vec)
#    for j in range(0,i-2):
#        x_f = j/(i-2)
#        x_l = (j+1)/(i-2)
#        k_f = k_vec[j]
#        k_l = k_vec[j+1]
#        kk = np.linspace(x_f,x_l,100)
#        plt.plot(kk,q(kk,x_f,x_l,k_f,k_l))

x_0 = 0.5
err = 1
num = 4
while err > .000001:
    for i in range(3,num):
        mat = np.zeros([i-1,i-1])
        vec = np.zeros(i-1)
        x = np.zeros(i-1)
        y = np.zeros(i-1)
        k_vec = fill_matrix(i,x,y,mat,vec)
        for j in range(0,i-2):
            x_f = j/(i-2)
            x_l = (j+1)/(i-2)
            if x_f < x_0 < x_l:
                k_f = k_vec[j]
                k_l = k_vec[j+1]
                x_new = x_0 - q(x_0,x_f,x_l,k_f,k_l)/der_q(x_0,x_f,x_l,k_f,k_l)
                err = np.abs(x_new-x_0)
                x_0 = x_new
        num = num + 1

root = x_0    
print(root)

x_0 = 0.5
err = 1
err_vec = np.zeros(N)
for k in range(1,200):
    for i in range(3,N):
        mat = np.zeros([i-1,i-1])
        vec = np.zeros(i-1)
        x = np.zeros(i-1)
        y = np.zeros(i-1)
        k_vec = fill_matrix(i,x,y,mat,vec)
        for j in range(0,i-2):
            x_f = j/(i-2)
            x_l = (j+1)/(i-2)
            if x_f < x_0 < x_l:
                k_f = k_vec[j]
                k_l = k_vec[j+1]
                x_new = x_0 - q(x_0,x_f,x_l,k_f,k_l)/der_q(x_0,x_f,x_l,k_f,k_l)
                err = np.abs(x_new-root)
                x_0 = x_new
        err_vec[i] = err

plt.figure()
plt.plot(err_vec)
plt.title('Error x_N - x_* as a function of N')