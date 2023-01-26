#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pentapy as pp
import matplotlib.pyplot as plt
import scipy
from scipy import sparse

N = 200
am = [1]* N
ad1 = [None]* (N-1)
ad2 = [None]* (N-2)

bm = [None]* N
bu1 = [None]* (N-1)

cm = [None]* N
cd1 = [None]* (N-1)

dm = [1]* N
du1 = [None]* (N-1)
du2 = [None]* (N-2)

b=[0]* (2*N-1) + [1]
V = np.array(b)

fix_prob=[]
for r in np.arange(0.01,3.00,0.0001):
    for i in range(N-1):
        ad1[i]=2*(i+1)*(N-i-2)/((r*(i+1)+N-i-1)*(N-1)*(N-2))+(N-i-2)*(i+1)/((r*(i+1)+N-i-1)*(N-2)) #the 1st lower diag of the 1st block
        du1[i]=2*r*i*(N-i-1)/((r*(i+1)+N-i-1)*(N-1)*(N-2))+r*i*(N-i-1)/((r*(i+1)+N-i-1)*(N-2))     #the 1st upper diag of the 4th block
    for i in range(N-2):
        ad2[i]=(i+2)*(i+1)/((r*(i+2)+N-i-2)*(N-1)*(N-2))       #the 2nd lower diag of the 1st block
        du2[i]=r*(N-i-1)*(N-i-2)/((r*(i+1)+N-i-1)*(N-1)*(N-2)) #the 2nd upper diag of the 4th block
    for i in range(1,N):
        am[i]=(r*i*(N-2)+i*(N-i-1))/((r*i+N-i)*(2-N))+(2*i*(N-i-1)+i*(i-1))/((r*i+N-i)*(N-1)*(2-N)) #the main diag of the 1st block
    for i in range(N-1):
        dm[i]=((N-i-1)*(N-2)+r*i*(N-i-1))/((r*(i+1)+N-i-1)*(2-N))+(2*r*i*(N-i-1)+r*(N-i-1)*(N-i-2))/((r*(i+1)+N-i-1)*(N-1)*(2-N))
                                                     #the main diag of the 4th block
        cd1[i]=(N-i-2)*(i+1)/((r*(i+2)+N-i-2)*(N-2)) #the 1st lower diag of the 3rd block 
        bu1[i]=r*i*(N-i-1)/((r*i+N-i)*(N-2))         #the 1st upper diag of the 2nd block
    for i in range(N):
        cm[i]=(N-i-1)*(N-i-2)/((r*(i+1)+N-i-1)*(N-2))#the main diag of the 3rd block
        bm[i]=r*i*(i-1)/((r*i+N-i)*(N-2))            #the main diag of the 2nd block
        
    Adiagonals = [am,ad1,ad2] #diagonals of the 1st block
    Bdiagonals = [bm,bu1]     #diagonals of the 2nd block
    Cdiagonals = [cm,cd1]     #diagonals of the 3rd block
    Ddiagonals = [dm,du1,du2] #diagonals of the 4th block
    A = scipy.sparse.diags(Adiagonals,[0,-1,-2]).toarray()
    B = scipy.sparse.diags(Bdiagonals,[0,1]).toarray()
    C = scipy.sparse.diags(Cdiagonals,[0,-1]).toarray()
    D = scipy.sparse.diags(Ddiagonals,[0,1,2]).toarray()
    P = np.block([
                  [A, B],
                  [C, D]
                ])   #creat the whole matrix with 4 blocks
    X = np.linalg.solve(P, V)  #solve the linear system
    f = X[N]/N +(1-1/N)*X[1]   #take the average fixation probability over all initial configuration with one single mutant
    fix_prob.append(f)
np.savetxt('3-star_N200_numerical.txt', fix_prob)


# In[ ]:


#r=np.arange(0.01,3.00,0.01)
#fix_moran=(1-1/r)/(1-1/(r**N))
#fix_star=(1-1/r**2)/(1-1/(r**(2*N)))
#plt.xlim([0.0,3.0])
#plt.ylim([0.0,1e-5])
#plt.tick_params(labelsize=15)
#plt.plot(r, fix_moran, label='moran process', color='royalblue')
#plt.plot(r, fix_prob, label='3-uniform star',color='magenta')
#plt.plot(r, fix_star, label='star')
#plt.xlabel('r',fontsize=18)
#plt.ylabel('fixation prob',fontsize=18)
#plt.legend(fontsize=15)
#plt.savefig('star_sim_N=50.pdf', bbox_inches="tight")
#plt.show()

