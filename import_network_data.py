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

# import the io and networkx module for use later
import argparse
import csv
import numpy as np 
import pickle
import sys
import networkx as nx
from collections import Counter
from math import factorial

# corr_list = list that holds all correlations (in first column)

# Create a dictionary function to easily add keys and values to dictionary
class dictionary(dict):
    def __init__(self):
        self = dict()
        
    # Add key/value pair
    def add(self, key, value):
        self[key] = value

corr_dict = dictionary()
fc = {}

parser = argparse.ArgumentParser(description='Example: python dev_from_expected --input <network file> --num_groups <number of groups> \n\n See README.md for more info\n\n')
parser.add_argument("--input", help = 'Network file (see README.md for example)')
parser.add_argument("--num_groups", help = 'Number of groups correlations were initially performed in')

args = parser.parse_args()

if (args.input != '') and (args.num_groups != ''):
    net_file = args.input
    net_file_trimmed = net_file[:-4] # trim the ".csv" or ".txt" from the input file string   
    groups = int(args.num_groups)
else:
    raise Exception("Error: please provide both a network file and the number of groups correlations were performed in. See --help for more information.")

# Counters for puc calculation
puc_compliant = 0
puc_noncompliant = 0

row_count = 0

# import specified file into python
with open(net_file) as csvfile:
    file = csv.reader(csvfile, delimiter = ',')
    for row in file:

        # Check if there are any NAs in the correlation column and if not then add the edge to the dictionary
        # The keys are nodes that make up the edge and values are a list of parameters (pval, comb pval, rho, etc.)        
        #if row[12] != 'NA':
            
        # Take the index of the source and target node in the header of the file
        if row_count == 0: 
            p1 = int(row.index('partner1'))
            p2 = int(row.index('partner2')) + 1
        
        nodes = row[p1:p2]
        list_to_tuple = tuple(nodes)
        corr_dict.add(list_to_tuple,row[3:len(row)])
        
        fc_node1_column = 11 + groups
        fc_node2_column = 12 + groups
        
        # Find FC direction of each node
        fc[row[1]] = row[fc_node1_column].strip()
        fc[row[2]] = row[fc_node2_column].strip()
        
        # Is each edge PUC-compliant?
        puc_col = 14 + groups
        
        if row[puc_col].strip() == str(1):
            puc_compliant += 1
        elif row[puc_col].strip() == str(-1):
            puc_noncompliant += 1
            
        #else:
        #    print("An NA was found for a correlation. Continuing anyways but omitting this correlation.")

        row_count += 1

    csvfile.close()

# This removes the header of the data frame essentially
del corr_dict['partner1', 'partner2']
del fc['partner1']
del fc['partner2']

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

#print("Strength of each node (sum of rhos)")
#print(G.degree(weight='weight'))

# Save to picklke
pickle_out = open(net_file_trimmed + ".pickle", "wb")

################################################################################
 ##################### Calculate deviation from ideal #########################
################################################################################
# Get a dictionary of all correlation directions and count all positive and negative edges in observed network
rho_column = 7 + groups
pos_corr = 0 # counter for the number of positive edges
neg_corr = 0 # counter for the number of negative edges

for key,value in corr_dict.items():
    try:    
        if str(value[rho_column].strip()) == '1':
            pos_corr += 1
        elif str(value[rho_column].strip()) == '-1':
            neg_corr += 1
    except:
        print("ERROR: an incorrect value was supplied for correlation directions. Aborting.")         
    
nedges = pos_corr + neg_corr       

# Count the number of positive and negative nodes       
nodedir = Counter(fc.values())
pos_nodes = nodedir['1']
neg_nodes = nodedir['-1']         
total_nodes = int(pos_nodes) + int(neg_nodes)
obs_edge_node_ratio = nedges / total_nodes

obs_posneg_node_ratio = int(pos_nodes) / int(neg_nodes)
obs_negpos_node_ratio = int(neg_nodes) / int(pos_nodes)

# Find the ratio of positive:negative edges (and vice versa) in the observed graph
if int(neg_corr) != 0:
    obs_posneg_ratio = int(pos_corr)/int(neg_corr)
else:
    obs_posneg_ratio = 1.0
    
if int(pos_corr) != 0:
    obs_negpos_ratio = int(neg_corr)/int(pos_corr)
else:
    obs_negpos_ratio = 1.0

# Find the number of edges in a full graph      
expec_pos = int(factorial(pos_nodes)/(2 * factorial(pos_nodes - 2)) + factorial(neg_nodes)/(2 * factorial(neg_nodes - 2)))           
expec_neg = pos_nodes * neg_nodes
expec_total = expec_pos + expec_neg          
expec_edge_node_ratio = expec_total / total_nodes

# Find the ratio of positive:negative edges (and vice versa) in a full graph
ideal_ratio_posneg = expec_pos/expec_neg
ideal_ratio_negpos = expec_neg/expec_pos

# Find the ratio of positive:negative edges (and vice versa) in a complete graph
ideal_ratio_posneg = expec_pos/expec_neg
ideal_ratio_negpos = expec_neg/expec_pos

#Calculate the non-normalized deviation from the expected (full) graph
dev_posneg = obs_posneg_ratio/ideal_ratio_posneg
dev_negpos = obs_negpos_ratio/ideal_ratio_negpos
                      
#Calculate the normalized deviation from the expected (full) graph
dev_norm_posneg = (obs_posneg_ratio - ideal_ratio_posneg) / ideal_ratio_posneg
dev_norm_negpos = (obs_negpos_ratio - ideal_ratio_negpos) / ideal_ratio_negpos

# calculate the normalized deviation of the edge:node (density) from the full graph
dens_dev = (abs(obs_edge_node_ratio - expec_edge_node_ratio)) / expec_edge_node_ratio

# Calculate PUC (the proportion of edges that do not follow the expected direction)
puc = puc_noncompliant / nx.number_of_edges(G)

dev_dict = {}
dev_dict['OBSERVED_number_nodes_positive'] = pos_nodes
dev_dict['OBSERVED_number_nodes_negative'] = neg_nodes
dev_dict['OBSERVED_number_edges_positive'] = pos_corr
dev_dict['OBSERVED_number_edges_negative'] = neg_corr
dev_dict['OBSERVED_ratio_pos_to_neg_nodes'] = round(obs_posneg_node_ratio, 2)
dev_dict['OBSERVED_ratio_neg_to_pos_nodes'] = round(obs_negpos_node_ratio, 2)
dev_dict['OBSERVED_edge_node_ratio'] = round(obs_edge_node_ratio, 2)
dev_dict['OBSERVED_ratio_pos_to_neg_edges'] = round(obs_posneg_ratio, 2)
dev_dict['OBSERVED_ratio_neg_to_pos_edges'] = round(obs_negpos_ratio, 2)
dev_dict['IDEAL_total_number_edges_full_graph'] = expec_total
dev_dict['IDEAL_number_positive_edges_full_graph'] = expec_pos
dev_dict['IDEAL_number_negative_edges_full_graph'] = expec_neg
dev_dict['IDEAL_density_full_graph'] = round(expec_edge_node_ratio, 2)
dev_dict['IDEAL_ratio_pos_to_neg_edges'] = round(ideal_ratio_posneg, 2)
dev_dict['IDEAL_ratio_neg_to_pos_edges'] = round(ideal_ratio_negpos, 2)

dev_dict['DEVIATION_posneg_deviation_nonnormalized'] = round(dev_posneg, 2)
dev_dict['DEVIATION_negpos_deviation_nonnormalized'] = round(dev_negpos, 2)
dev_dict['DEVIATION_posneg_deviation_normalized'] = round(dev_norm_posneg, 2)
dev_dict['DEVIATION_negpos_deviation_normalized'] = round(dev_norm_negpos, 2)
dev_dict['DEVIATION_density_deviation'] = round(dens_dev, 2)
dev_dict['PUC'] = round(puc, 2)
dev_dict['PUC_noncompliant_edge_number'] = round(puc_noncompliant, 2)
dev_dict['PUC_compliant_edge_number'] = round(puc_compliant, 2)

#pickle.dump(G, pickle_out)
pickle_list = [G, dev_dict]
pickle.dump(pickle_list, pickle_out)
pickle_out.close()


print("Successfully saved to pickle. Now run calc_network_properties.py")























