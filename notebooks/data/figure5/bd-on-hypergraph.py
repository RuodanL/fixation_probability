#!/usr/bin/env python
# coding: utf-8

# In[1]:


import networkx as nx
import random
import math
import numpy as np
import numpy.random as rnd
from collections import defaultdict
from networkx.algorithms import bipartite
from itertools import combinations
import matplotlib.pyplot as plt
import timeit


# In[2]:


#start = timeit.default_timer()
class _ListDict_(object):
     r'''
     
    The birth-death process will involve a step that randomly select an element from a set. 
    
    This class will allow me to select a random element uniformly, and then either add it or delete it from  
    a set.
    
     '''        
    def __init__(self):
        self.item_to_position = {}
        self.items = []
    
    def __len__(self):
        return len(self.items)
    
    def __contains__(self, item):
        return item in self.item_to_position
    
    def insert(self, item):
        if self.__contains__(item):
            self.remove(item)
            
    def update(self, item):
        if item in self:
            return
        self.items.append(item)
        self.item_to_position[item] = len(self.items)-1
    
    def remove(self, choice):
        position = self.item_to_position.pop(choice)
        last_item = self.items.pop()
        if position != len(self.items):
            self.items[position] = last_item
            self.item_to_position[last_item] = position
            
    def choose_random(self):
        return random.choice(self.items)
    
    def random_removal(self):
        choice = self.choose_random()
        self.remove(choice)
        return choice
    
    def total_weight(self):
        return len(self)


# In[3]:


def moran_bd(G, r, bottom_nodes, initial_mutants):
    r'''    
    
    Performs birth-death process simulations for bipartite networks.
    
    **G** networkx Graph
        The underlying bipartite network
    **r**
        the fitness of each mutant
    **bottom_nodes**
        nodes of the hypergraph 
    **initial_mutants**   
        one bottom node in G is initially infected
        
    '''  
    status = defaultdict(lambda : 'R')
    status[initial_mutants] = 'M'      #assign the status for each node, R for resident, M for mutant
        
    mutants = _ListDict_()    #set of mutants
    residents = _ListDict_()  #set of residents
    
    for node in bottom_nodes: 
        residents.update(node)  #add all bottom nodes to the set of residents
    
    mutants.update(initial_mutants)  #add the initial mutant node to the set of mutants
    residents.remove(initial_mutants)#remove the initial mutant from the set of residents
    
        
    total_mutants_fitness = r*mutants.total_weight() #r*mutants_sum
    total_residents_fitness = len(bottom_nodes)-mutants.total_weight() #residents_sum
    total_fitness = total_mutants_fitness + total_residents_fitness
        
    while mutants.total_weight() > 0 and mutants.total_weight() < len(bottom_nodes): #perform the birth-death process if neither mutants nor residents die out
        if random.random()<total_mutants_fitness/total_fitness: #the event in which a mutant reproduces
            reproducing_node = mutants.choose_random()
            hyperedge = random.choice(list(G.neighbors(reproducing_node)))
            for nbr in G.neighbors(hyperedge):
                if status[nbr] == 'R':
                    status[nbr] = 'M'
                    mutants.update(nbr)
                    residents.remove(nbr)
        else:
            reproducing_node = residents.choose_random() #the event in which a resident reproduces
            hyperedge = random.choice(list(G.neighbors(reproducing_node)))
            for nbr in G.neighbors(hyperedge):
                if status[nbr] == 'M':
                    status[nbr] = 'R'
                    mutants.remove(nbr)
                    residents.update(nbr)
    
        total_mutants_fitness = r*mutants.total_weight()
        total_residents_fitness = len(bottom_nodes)-mutants.total_weight()
        total_fitness = total_mutants_fitness + total_residents_fitness
                
    return mutants.total_weight()


# In[4]:


def fixation_prob(G, r, bottom_nodes, repetitions, initial_mutants):
    r'''    
    
    get the fraction of the runs in which mutants has fixated 
    
    '''  
    fixation_count = 0
    for i in range(repetitions):
        final_counts = moran_bd(G, r, bottom_nodes, initial_mutants)
        if final_counts > 0:
            fixation_count += 1
    
    return fixation_count / repetitions        


# In[5]:


repetitions = 3000 #number of runs for each initial mutant


# In[6]:


G=nx.read_edgelist("edgelist.txt")
bottom_nodes, top_nodes = nx.bipartite.sets(G)


# In[7]:


fix_avr=[]
for r in np.arange(1.00,3.00,0.01): #for each fitness r, take the average of the fixation probability of all bottom nodes 
    fix_sim=[]
    for initial_mutants in bottom_nodes:
        fix_prob=fixation_prob(G, r, bottom_nodes, repetitions, initial_mutants)
        fix_sim.append(fix_prob)
    fix_mean=np.mean(fix_sim)
    fix_avr.append(fix_mean)

#np.savetxt('club_membership.txt', fix_avr)

#stop = timeit.default_timer()
#print('Time: ', stop - start)


# In[ ]:




