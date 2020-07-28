# network_analysis
 
 Author: Nolan K Newman <newmanno@oregonstate.edu>
 Date: 7/27/20
 Created with Python v3.5.3
 
 Purpose:
 Imports network and calculate analyses on node and network properties

 This repository consists of two programs: 
	- import_network_data.py: Imports a network file consisting of node correlations and outputs a pickled file of a networkx graph
	- calc_network_properties.py: Takes the output from import_network_data.py and performs analyses on the network, which are detailed below
	
######################################
import_network_data.py

Purpose:
This script takes a network file consisting of node correlations and outputs a pickled file of a networkx graph

Dependencies: pickle, networkx 

Required inputs:
	The input is a CSV file with the following details:

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
	
Arguments:
	This code does not currently accept any optional arguments

Example input:
	python import_network_data.py <network file>
	
Outputs:
	pickled file of the networkx graph, including each edge's abs(rho coefficient)
	

######################################
calc_network_properties.py

Purpose:
	Takes the output from import_network_data.py and performs analyses on the network

Description:
	This script takes the output of import_network_data.py, which is a pickled network file containing each node of the network and its associated rho coefficient (or rather the absolute value of its rho coefficient). It will calculate the following:
	
	Network properties:
		- Number of nodes in the network
		- Number of edges in the network
		- Ratio of positive to negative correlations
		- Mean degree (2 * # edges / # nodes)
		- Average clustering coefficient
		- Average path length
		- Giant component size (percent of network that is the giant component)
		- Number of connected components
		- Freeman centrality
		- Mean closeness centrality
		- Modularity
		- Median component size over the total number of nodes
		- Degree assortativity
		- Proportion of unexpected correlations (PUC): See here https://biologydirect.biomedcentral.com/articles/10.1186/s13062-016-0155-0
		- Network fragmentation

	Node properties:
		- Degree
		- Strength (sum of the absolute values of a node's edges)
		- Clustering coefficient
		- Closeness centrality
		- Eigenvector centrality
		- Betweenness centrality
		- Fragmentation centrality (A value of how the shortest path distances change as each node is removed, then replaced in the network)
		- Bipartite betweenness centrality (BiBC): See here https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4415676/
		
	Arguments:
		Required:
			- Pickled network file (output of import_network_data.py)
			
		Optional
			--PUC	-	

	Notes: 
		BiBC will not be calculated if the network size is either too small or if two separate groups are identified in the network









 
	
	
