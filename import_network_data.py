# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:01:46 2019

@author: Nolan
"""

# Script that takes as input the output of the create-network.R script and calculates network properties.

# Example: python import_network_data2.1.py <network_file>

###################################################################################
 ##################### import network and make networkx object ###################
###################################################################################

# import the io and networkx module
import csv
import re
import numpy as np 
import pickle
import sys
import os
import matplotlib
import networkx as nx
import matplotlib.pyplot as plt

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

# This removes the header of the data frame essentially
del corr_dict['partner1', 'partner2']

# function to get unique values, since python does not have a native 
# unique() function
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
     ############# Calculate ratio of +ve to -ve correlations #############
################################################################################
# Check to see if both nodes of each key of corr_dict are part of the giant component
# and if so, then check how many of those are positive correlations and how many are negative

g = max(nx.connected_component_subgraphs(G), key = len)

c = sorted(nx.connected_components(G), key = len, reverse = True)
largest_comp_size = len(c[0])

# Check if there is more than one giant component
try:
    second_largest_comp_size = len(c[1])
except:
    print("There is only one component")
    second_largest_comp_size = -1

with open(net_file + "_pos_and_neg_corrs.txt", "w") as file: 
 
    pos_corr = 0
    neg_corr = 0

    # If there is only one giant component (which should typically be the case)
    if (largest_comp_size != second_largest_comp_size): 

	# Loop through all key-value pairs in the corr_dict and for each pair...
        for key,value in corr_dict.items():
	   
           #print(key,value)		

           # Check the partner1 FoldChange direction (which should match partner 2) and count how many are negative correlations and how many are positive correlations. 
           if (key[0] in g.nodes()) & (key[1] in g.nodes()):
                if float(corr_dict[key[0],key[1]][9]) < 0:
                    neg_corr += 1 
                elif float(corr_dict[key[0],key[1]][9]) > 0:
                    pos_corr += 1
                else: 
                    print("ERROR: Correlation is neither positive or negative.")

	# If there are no negative correlations then calculate ratio of positive to negative edges, otherwise set to undefined (cannot divide by 0)
        if neg_corr != 0:
            pos_to_neg_ratio = pos_corr/neg_corr
        else:
            pos_to_neg_ratio = "undefined"

        file.write("Number_of_pos_corrs_in_GC\t " + str(pos_corr) + "\n")
        file.write("Number_of_neg_corrs_in_GC\t " + str(neg_corr) + "\n")
        file.write("\n")
        file.write("Pos_to_neg_ratio_in_GC\t " + str(pos_to_neg_ratio) + "\n")
        file.write("Giant_comp_size\t " + str(len(g)) + "\n")
        file.write("Num_of_total_nodes_in_net\t " + str(len(G)))
 
    else:
        for key,value in corr_dict.items():
            if float(corr_dict[key[0],key[1]][9]) < 0:
                neg_corr += 1
            elif float(corr_dict[key[0],key[1]][9]) > 0:
                pos_corr += 1
            else: 
                print("ERROR: Correlation is neither positive or negative.")

        file.write("There are multiple giant components. Each has " + str(len(g)) + " nodes. ")
        file.write("These have been calculated for the entire network, not the giant component.\n")
        file.write("Number_of_pos_corrs_in_GC\t " + str(pos_corr) + "\n")
        file.write("Number_of_neg_corrs_in_GC\t " + str(neg_corr) + "\n")
        file.write("Pos_to_neg_ratio_in_GC\t No_large_comp\n")
        file.write("Num_of_total_nodes_in_net\t " + str(len(G)))

    file.close()

########################################################
#Ignore this section, calculates pos. and neg. corr direction for whole graph, not giant component
########################################################
#with open(net_file + "_pos_and_neg_corrs.txt", "w") as file: 
#
#    #giant = max(nx.connected_component_subgraphs(G), key = len)  
#    #print(giant)
#
#
#    for key,value in corr_dict.items():
#        if float([value][0][0]) < 0:
#            neg_corr.append([value][0][0])
#        elif float([value][0][0]) > 0: 
#            pos_corr.append([value][0][0])
#        else:
#            print("Values were found to neither be positive or negative")
#

################################################################################
 ############################## Save to pickle ################################
################################################################################

# Pickles can be used to quickly save data in one program and be accessed in another (in this case
# calc_network_properties_vXX.py will access the pickled data
pickle_out = open(net_file+".pickle", "wb")
pickle.dump(G, pickle_out)
pickle_out.close()

print("Successfully saved to pickle. Now run calc_network_properties.py")























