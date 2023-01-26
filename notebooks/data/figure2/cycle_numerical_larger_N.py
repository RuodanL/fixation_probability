#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pentapy as pp
import matplotlib.pyplot as plt
import scipy
from scipy import sparse

N = 200
u1=[0]* N
u2=[0]* (N-1)
d1=[0]* N
d2=[0]* (N-1)
m=[1]* (N+1)
b=[0]*N + [1]
V = np.array(b)

fix_prob=[]
for r in np.arange(0.01,3.00,0.0001):
    m[1]=(-2-r)/(r+N-1)
    m[N-1]=(-1-2*r)/(r*N-r+1)
    u1[N-1]=2*r/(r*N-r+1)
    u2[1]=r/(r+N-1)
    d1[0]=2/(r+N-1)
    d2[N-3]=1/(r*N-r+1)
    for i in range(1,N-2):
        d1[i]=4/(3*r*i+3*r+3*N-3*i-3) #the 1st lower diagonal
    for i in range(N-3):
        d2[i]=2/(3*r*i+6*r+3*N-3*i-6) #the 2nd lower diagonal
    for i in range(2,N-1):
        m[i]=(-2-2*r)/(r*i+N-i)       #the main diagonal
        u1[i]=4*r/(3*r*i+3*N-3*i)     #the 1st upper diagonal
        u2[i]=2*r/(3*r*i+3*N-3*i)     #the 2nd upper diagonal
    
    diagonals = [m,d1,u1,d2,u2]
    mat = scipy.sparse.diags(diagonals,[0,-1,1,-2,2]).toarray() #create the pentadiagonal matrix M
    X = pp.solve(mat, V, is_flat=False, solver=2)  #solve the linear system
    fix_prob.append(X[1])
np.savetxt('3-cyclic_N200_numerical.txt', fix_prob)


# In[ ]:


#r=np.arange(0.01,3.00,0.01)
#fix_moran=(1-1/r)/(1-1/(r**N))
#plt.xlim([0.01,3.0])
#plt.title('3-uniform cycle with N=10',fontsize=20)
#plt.tick_params(labelsize=15)
#plt.plot(r, fix_moran, label='moran process')
#plt.plot(r, fix_prob, label='3-uniform cycle')
#plt.xlabel('r',fontsize=18)
#plt.ylabel('fixation prob',fontsize=18)
#plt.legend(fontsize=15)
#plt.savefig('cycle_sim_N=10.pdf', bbox_inches="tight")
#plt.show()

