#!/usr/bin/env python
# coding: utf-8

# In[8]:


import numpy as np
import pentapy as pp
import matplotlib.pyplot as plt
import scipy
from scipy import sparse

N = 200
am = [1]* N
au1 = [0]* (N-1)
ad1 = [None]* (N-1)

bm = [0]* N

cm = [None]* N

dm = [1]* N
du1 = [None]* (N-1)
dd1 = [None]* (N-1)

b=[0]* (2*N-1) + [1]
V = np.array(b)

fix_prob=[]
for r in np.arange(0.01,3.00,0.0001):
    for i in range(1,N):
        am[i]=-(i/((r*i+N-i)*(N-1))+i*(N-i-1)/((r*i+N-i)*(2*N-4))+r*i*(N-i-1)/((r*i+N-i)*(2*N-4))+r*i/(2*r*i+2*N-2*i))
                                                                   #the main diagonal of the 1st block
        ad1[i-1]=i/((r*i+N-i)*(N-1))+(N-i-1)*i/((r*i+N-i)*(2*N-4)) #the 1st lower diag of the 1st block
        dm[i-1]=-((N-i)/(2*r*i+2*N-2*i)+(N-i)*(i-1)/((r*i+N-i)*(2*N-4))+r*(N-i)/((r*i+N-i)*(N-1))+r*(i-1)*(N-i)/((r*i+N-i)*(2*N-4)))
                                                                   #the main diag of the 4th block
        du1[i-1]=r*(N-i)/((r*i+N-i)*(N-1))+r*(i-1)*(N-i)/((r*i+N-i)*(2*N-4))
                                                                   #the 1st upper diag of the 4th block
        bm[i]=r*i/(2*(r*i+N-i))                                    #the main diag of the 2nd block
    for i in range(1,N+1):
        cm[i-1]=(N-i)/(2*(r*i+N-i))                                #the main diag of the 3rd block 
    for i in range(1,N-1):
        au1[i]=r*i*(N-i-1)/((r*i+N-i)*(2*N-4))                     #the 1st upper diag of the 1st block
    for i in range(2,N+1):
        dd1[i-2]=(N-i)*(i-1)/((r*i+N-i)*(2*N-4))                   #the 1st lower diag of the 4th block
    
    Adiagonals = [am,au1,ad1]
    Bdiagonals = [bm]
    Cdiagonals = [cm]
    Ddiagonals = [dm,du1,dd1]
    A = scipy.sparse.diags(Adiagonals,[0,1,-1]).toarray()
    B = scipy.sparse.diags(Bdiagonals,[0]).toarray()
    C = scipy.sparse.diags(Cdiagonals,[0]).toarray()
    D = scipy.sparse.diags(Ddiagonals,[0,1,-1]).toarray()
    P = np.block([
                  [A, B],
                  [C, D]
                ])  #creat the matrix with four blocks
    X = np.linalg.solve(P, V) #solve the linear system
    f = X[N]/N +(1-1/N)*X[1]  #take the average fixation probability over all initial configuration with one single mutant
    fix_prob.append(f)
np.savetxt('omp-star_N200_numerical.txt', fix_prob)


# In[ ]:


#r=np.arange(0.01,3.00,0.01)
#fix_moran=(1-1/r)/(1-1/(r**N))
#fix_star=(1-1/r**2)/(1-1/(r**(2*N)))
#fix_sim = np.loadtxt('wstar_fix_sim_N15fast.txt')
#fix_sim = np.loadtxt('one_mode_proj_6nodes.txt')
#fix_sim = np.loadtxt('one_mode_proj_davis_women.txt')
#plt.title('weighted star with N=15',fontsize=20)
#plt.xlim([0.0,3.0])
#plt.ylim([0.5,0.55])
#plt.tick_params(labelsize=15)
#plt.scatter(r, fix_sim, label='simulation', marker='.', s=10, color='dodgerblue')
#plt.plot(r, fix_moran, label='moran process', color='black')
#plt.plot(r, fix_prob, label='star with weight', color='magenta')
#plt.plot(r, fix_star, label='star')
#plt.xlabel('r',fontsize=18)
#plt.ylabel('fixation prob',fontsize=18)
#plt.legend(fontsize=15)
#plt.savefig('wstar_N=15ctm.pdf', bbox_inches="tight")
#plt.savefig('weighted_6nodes.pdf', bbox_inches="tight")
#plt.show()


# In[ ]:




