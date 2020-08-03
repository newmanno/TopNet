# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 18:36:30 2020

@author: Nolan

Wrapper file to coordinate communication between the executable infomap program on the command line. Converts the NetworkX graph into a *.net file that infomap can accept as input, then runs infomap and returns a dictionary of partitions for each node

"""

import networkx as nx
import argparse
import pickle
import os

parser = argparse.ArgumentParser(description='Test')
parser.add_argument('pickle', help = 'The pickle file created with import_network_data.py')

args = parser.parse_args()

if __name__ == '__main__':

    def run_infomap(G):
          
        with open(network_name + "_infomap.net", "w") as im:
            im.write("*Vertices " + str(nx.number_of_nodes(G)) + "\n")
            
            # Convert nodes into integers for infomap input
            counter = 1 # Holds current node number in list
            infomap_int_name = {} # Holds numbers (ins str format) as keys and node names as values
            infomap_name_int = {} # Holds node names as keys and numbers (ins str format) as values 
        
            # Create pajek style input file for infomap
            for i in sorted(G.nodes):
                infomap_int_name[str(counter)] = i
                infomap_name_int[str(i)] = str(counter)
                counter += 1
              
            for key,value in sorted(infomap_int_name.items()):    
                im.write(" " + key + " \"" + value + "\"\n")
                
            im.write("*Edges " + str(nx.number_of_edges(G)) + "\n")
            
            for edge in sorted(list(G.edges)): 
                source = infomap_name_int[edge[0]]
                target = infomap_name_int[edge[1]]
                im.write(source + " " + target + "\n")
         
            ### ERROR: infomap does not recognize nodes when running from within python for some reason
            # os.system("echo infomap " + network_name + "_infomap.net" + " . -i pajek -vvv")
            # os.system("infomap " + network_name + "_infomap.net" + " . -i pajek -vvv")
          
        im.close()
    
        infomap_output = network_name + "_infomap.tree"
        
        # Read in the file created by infomap, which contains the partitions for each node        
        with open(infomap_output) as im_out:
            partition_lines = im_out.readlines()[8:] # Skip first 9 lines of file
            
            # Create and populate a dictionary where the keys are the node name and the values are the partition the node belongs to
            part_dict = {}
    
            for row in partition_lines:
                node_name = row.split()[2]
                node_name_strip = node_name.strip('\"') # Remove the quotes surrounding node name
                partition_num = row.split()[0][0]
                
                part_dict[node_name_strip] = partition_num
            
            # print("\nNew dictionary with node names and partitions from infomap\n\n")
            # for key,value in sorted(part_dict.items()):
            #     print("key: " + key + ", value: " + value)
                
        im_out.close()

        return(part_dict)
        
    # Unpack the pickle
    p = open(args.pickle, "rb")
    p = pickle.load(p)
    
    G = p[0] # Graph stored in first position of pickle
            
    network_name = args.pickle[:-7] # remove the pickle extension from file name        
      
    im_res = run_infomap(G)
        
    for key,value in sorted(im_res.items()):
        print("key: " + key + ", value: " + value)       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        