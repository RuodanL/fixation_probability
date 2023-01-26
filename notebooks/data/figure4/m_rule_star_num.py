#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from sympy import *
import matplotlib.pyplot as plt

N=200
fix_major=[]
for r in np.arange(0.01,3.00,0.0001):
    a=[0]*N
    b=[0]*N
    m=[0]*N
    for i in range(2, N):
        a[i]=1+(N-i-1)*(N-1)/(r*(i-1)*(N+1))                    #assign value for alpha_i
        b[i]=r*(i-1)*(N-1)/(r*(i-1)*(N-1)+(N-i-1)*(N+1))        #assign value for beta_i
        m[i]=Matrix([[a[i], 1-a[i]], [a[i]*b[i], 1-a[i]*b[i]]]) #create the matrix A_i
    p=Matrix([[1,0],[0,1]])
    for i in range(N-1,1,-1):
        p=p*m[i]  #take the product between A_i's from i=N-1 to i=2
    fix_prob=2/(N*p[0,0])+(N-2)*a[2]*b[2]/(N*p[0,0]) #fixation probability of the initial configuration with two mutants
    fix_major.append(fix_prob)
np.savetxt('m2_star_N200_numerical.txt', fix_major)


# In[ ]:




