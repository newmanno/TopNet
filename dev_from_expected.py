# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 11:47:38 2020

@author: Nolan K Newman <newmanno@oregonstate.edu>

This code takes as input the output of Richard's pipeline (see below) and will calculate the ratio of positive to negative nodes, edges/nodes, number of positive correlations, number of negative correlations, total number of correlations, and the deviation from expected raio. The deviation from expected is calculated by taking the observed positive:negative correlation ratio and dividing it by the ideal positive to negative ratio. The ideal ratio can be calculated by considering a complete graph with the same number of input nodes, then calculating how many edges would be negative and positive, assuming the complete graph is PUC-compliant (i.e. all pos-pos or neg-neg edges are a positive correlation and all pos-neg edges are a negative correlation.

NOTE: THESE VALUES ARE NOT NORMALIZED SO ONE NETWORK WITH A 50:50 RATIO VS ONE NETWORK THAT HAS A RATIO OF 1:99 ARE NOT YET COMPARABLE. WE WILL BE DOING THE NORMALIZATION AT A LATER. 

Python version: 3.5.1
Dependencies: pickle, networkx 

Example header for file (CSV)
pair,partner1,partner2,pval_E1,pval_E2,comb_pval,comb_rho,comb_FDR,partner1InFold,partner1_FoldChange,partner2InFold,partner2_FoldChange,corr_direction,partner1_FC_direction,partner2_FC_direction,IfFoldChangeDirectionMatch,PUC
    - pair: gene1 <==> gene2
    - partner1: gene1
    - partner2: gene2
    - pval_E1: correlation pvalue in experiment 1
    - pval_E2: correlation pvalue in Experiment 2
    - comb_pval: Fisher's combined pvalue across both experiments
    - comb_rho: combined rho coefficient across both experiments
    - comb_FDR: FDR calculated off the combined pvalue
    - partner1InFold: gene1
    - partner1_FoldChange: Fold change of gene1
    - partner2InFold: gene2
    - partner2_FoldChange: Fold change of gene2
    - corr_direction: correlation direction (either -1 or 1)
    - partner1_FC_direction: Fold change direction of gene1 (either -1 or 1)
    - partner2_FC_direction: Fold change direction of gene2 (either -1 or 1)    
    - IfFoldChangeDirectionMatch: Are the previous 2 values identical (1 if yes, -1 if no)
    - PUC: Is the correlation PUC-compliant (are neg-neg or pos-pos correlations +ve and pos-neg correlations -ve?)
"""

# import the io and networkx module
import csv
import sys
from collections import Counter
from math import factorial


# Create a dictionary function to easily add keys and values to dictionary
class dictionary(dict):
    def __init__(self):
        self = dict()
        
    # Add key/value pair
    def add(self, key, value):
        self[key] = value

corr_dict = dictionary()
fc = {}

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
            
            # Find FC direction of each node
            fc[row[1]] = row[13]
            fc[row[2]] = row[14]
        else:
            print("An NA was found for a correlation. Continuing anyways but omitting this correlation.")
            
del fc['partner2']
del fc['partner1']

# This removes the header of the data frame essentially
del corr_dict['partner1', 'partner2']

# Get a dictionary of all correlation directions and count all positive and negative edges in observed network
pos_corr = 0 # counter for the number of positive edges
neg_corr = 0 # counter for the number of positive edges
for key,value in corr_dict.items():
    try:    
        if value[9] == '1':
            pos_corr += 1
        elif value[9] == '-1':
            neg_corr += 1
    except:
        print("ERROR: an incorrect value was supplied for correlation directions. Aborting.")         
     
nedges = pos_corr + neg_corr        
        
print("There are " + str(pos_corr) + " positive edges and " + str(neg_corr) + " negative edges for a total of " + str(nedges) + " edges.")

# Count the number of positive and negative nodes       
nodedir = Counter(fc.values())
pos_nodes = nodedir['1']
neg_nodes = nodedir['-1']         
total_nodes = int(pos_nodes) + int(neg_nodes)
obs_edge_node_ratio = nedges / total_nodes

obs_posneg_node_ratio = int(pos_nodes) / int(neg_nodes)
obs_negpos_node_ratio = int(neg_nodes) / int(pos_nodes)

print("There are " + str(pos_nodes) + " positive nodes and " + str(neg_nodes) + " negative nodes.")
print("The observed negative:positive node ratio is " + str(obs_negpos_node_ratio))            
print("The observed positive:negative node ratio is " + str(obs_posneg_node_ratio))            


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
           
#Calculate the non-normalized deviation from the expected (full) graph
dev_posneg = obs_posneg_ratio/ideal_ratio_posneg
dev_negpos = obs_negpos_ratio/ideal_ratio_negpos
                      
#Calculate the normalized deviation from the expected (full) graph
dev_norm_posneg = (obs_posneg_ratio - ideal_ratio_posneg) / ideal_ratio_posneg
dev_norm_negpos = (obs_negpos_ratio - ideal_ratio_negpos) / ideal_ratio_negpos

# calculate the normalized deviation of the edge:node (density) from the full graph
dens_dev = (abs(obs_edge_node_ratio - expec_edge_node_ratio)) / expec_edge_node_ratio

with open("deviation_output.txt", "w") as file:
    file.write("Number of total nodes: " + str(total_nodes) + "\n")    
    file.write("Number of positive nodes: " + str(pos_nodes) + "\n")
    file.write("Number of negative nodes: " + str(neg_nodes) + "\n")
    file.write("Positive:negative node ratio: %.3f"  % obs_posneg_node_ratio + "\n")   
    file.write("Negative:positive node ratio: %.3f"  % obs_negpos_node_ratio + "\n\n") 
  
    file.write("### OBSERVED values ###\n") 
    file.write("Observed number of total edges: " + str(nedges) + "\n")        
    file.write("Observed number of positive edges: " + str(pos_corr) + "\n")  
    file.write("Observed number of negative edges: " + str(neg_corr) + "\n")
    file.write("Observed density (edge:node ratio): " + str(obs_edge_node_ratio) + "\n")     
    file.write("Observed positive:negative edge ratio: %.3f" % obs_posneg_ratio + "\n")   
    file.write("Observed negative:positive edge ratio: %.3f" % obs_negpos_ratio + "\n\n")
    
    file.write("### EXPECTED values ###\n")   
    file.write("Full graph number of total edges: " + str(expec_total) + "\n")   
    file.write("Full graph number of positive edges: " + str(expec_pos) + "\n")
    file.write("Full graph number of negative edges: " + str(expec_neg) + "\n")   
    file.write("Full graph density (edge:node ratio): " + str(expec_edge_node_ratio) + "\n")       
    file.write("Full graph positive:negative edge ratio: %.3f" % ideal_ratio_posneg + "\n")
    file.write("Full graph negative:positive edge ratio: %.3f" % ideal_ratio_negpos + "\n\n")   
    
    file.write("### Departures from full graph ###\n")       
    file.write('Non-normalized positive:negative deviation from full graph ratio (observed ratio / full graph ratio): %.2f' % dev_posneg + "\n")
    file.write('Non-normalized negative:positive deviation from full graph ratio (observed ratio / full graph ratio): %.2f' % dev_negpos + "\n")
    file.write('Normalized positive:negative deviation from full graph ((observed ratio - full graph ratio) / full graph ratio): %.2f' % dev_norm_posneg + "\n")      
    file.write('Normalized negative:positive deviation from full graph ((observed ratio - full graph ratio) / full graph ratio): %.2f' % dev_norm_negpos + "\n")          
    file.write('Normalized density (edge:node ratio) deviation from full graph ((observed ratio - full graph ratio) / full graph ratio): %.5f' % dens_dev + "\n")            
            
            
