#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pentapy as pp
import matplotlib.pyplot as plt
import scipy
from scipy import sparse

N = 200
u1=[None]* N      
u2=[None]* (N-1)
d1=[None]* N
d2=[None]* (N-1)
m=[1]* (N+1)
b=[0]*N + [1]
V = np.array(b)

fix_prob=[]
for r in np.arange(0.01,3.00,0.0001):
    for i in range(N):
        u1[i]=(r*i/(r*i+N-i))*(2*(i-1)*(N-i)/((N-1)*(N-2))) #the 1st upper diagonal
        d1[i]=((N-i-1)/(r*i+r+N-i-1))*(2*(i+1)*(N-i-2)/((N-1)*(N-2))) #the 1st lower diagonal
    for i in range(N-1):
        u2[i]=(r*i/(r*i+N-i))*((N-i-1)*(N-i)/((N-1)*(N-2))) #the 2nd upper diag
        d2[i]=((N-i-2)/(r*i+2*r+N-i-2))*((i+2)*(i+1)/((N-1)*(N-2)))   #the 2nd lower diag
    for i in range(1,N):
        m[i]=((N-i)*(i-1)*(i+2*r*i)+(N-i)*(N-i-1)*(2*i+r*i))/((r*i+N-i)*(1-N)*(N-2)) #the main diag
    
    diagonals = [m,d1,u1,d2,u2]
    mat = scipy.sparse.diags(diagonals,[0,-1,1,-2,2]).toarray() #create the pentadiagonal matrix M
    X = pp.solve(mat, V, is_flat=False, solver=2)  #solve the linear system
    fix_prob.append(X[1])
np.savetxt('3-complete_N200_numerical.txt', fix_prob)


# In[ ]:


#r=np.arange(0.01,3.00,0.01)
#fix_moran=(1-1/r)/(1-1/(r**N))
#plt.xlim([0.01,3.0])
#plt.tick_params(labelsize=15)
#plt.plot(r, fix_moran, label='moran process')
#plt.plot(r, fix_prob, label='3-uniform complete')
#plt.xlabel('r',fontsize=18)
#plt.ylabel('fixation prob',fontsize=18)
#plt.legend(fontsize=13)
#plt.show()

