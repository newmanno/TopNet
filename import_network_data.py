# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:01:46 2019

Written in Python v3.5.3

@author: Nolan K Newman <newmanno@oregonstate.edu>

This script takes a network file consisting of node correlations and outputs a pickled file of a networkx graph

Example command: python import_network_data.py <network_file>

See README file for file input format
"""

###################################################################################
 ##################### import network and make networkx object ###################
###################################################################################

# import the io and networkx module
import csv
import numpy as np 
import pickle
import sys
import networkx as nx

# corr_list = list that holds all correlations (in first column)

# Create a dictionary function to easily add keys and values to dictionary
class dictionary(dict):
    def __init__(self):
        self = dict()
        
    # Add key/value pair
    def add(self, key, value):
        self[key] = value

corr_dict = dictionary()

net_file = sys.argv[1]

# import specified file into python
with open(net_file) as csvfile:
    file = csv.reader(csvfile, delimiter = ',')
    for row in file:

        # Check if there are any NAs in the correlation column and if not then add the edge to the dictionary
        # The keys are nodes that make up the edge and values are a list of parameters (pval, comb pval, rho, etc.)        
        if row[12] != 'NA':
            nodes = row[1:3]
            list_to_tuple = tuple(nodes)
            corr_dict.add(list_to_tuple,row[3:len(row)])
        else:
            print("An NA was found for a correlation. Continuing anyways but omitting this correlation.")

    csvfile.close()

# This removes the header of the data frame essentially
del corr_dict['partner1', 'partner2']

# function to get unique values
def unique(list1):
    x = np.array(list1)
    return np.unique(x)

# add all nodes to a list
node_list = []
for key,value in corr_dict.items():
    node_list.append(key[0])
    node_list.append(key[1]) 
        
# find all nodes and store them in a new list
unique_nodes = list(unique(node_list))

G = nx.Graph() 

# Find the absolute value of each rho edge (to be added for node strength)
print("Node1\tnode2\tabs(Edge weight)")
for key,value in corr_dict.items():
    #print(corr_dict[e][1])
    G.add_edge(key[0],key[1],weight = abs(float(value[3])))
    print(key[0],"\t",key[1],"\t",abs(float(value[3])))        

print("Strength of each node (sum of rhos)")
print(G.degree(weight='weight'))

################################################################################
 ############################## Save to pickle ################################
################################################################################

# Pickles can be used to quickly save data in one program and be accessed in another (in this case
# calc_network_properties.py will access the pickled data
pickle_out = open(net_file+".pickle", "wb")
pickle.dump(G, pickle_out)
pickle_out.close()

print("Successfully saved to pickle. Now run calc_network_properties.py")























