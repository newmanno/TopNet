# TopNet: A network topology analysis program written using NetworkX
> This network analysis pipeline takes a user-submitted network and calculates topological properties of the network.

## Table of contents
* [General info](#general-info)
* [Dependencies](#dependencies)
* [Importing the network](#importing-the-network)
* [Calculating topological properties](#calculating-topological-properties)

## General info
* Author: [Nolan K Newman](http://blogs.oregonstate.edu/morgunshulzhenkolabs/members/nolan-newman/) - newmanno@oregonstate.edu
* Date: 8/2/2020
* Created with Python v3.5.3

## Dependencies
* NetworkX (v2.2)
* pickle (v4.0)
* Numpy (v1.16.4

## Importing the network
Purpose: This script takes a network file consisting of node correlations and outputs a pickled file of a NetworkX graph

Arguments:
> **Required**
> * Network file (see below for format)

> **Optional**
> * This code does not currently accept any optional arguments
	
Example input command:
* python import_network_data.py <network file>
	
Example input file format:
* The input is a **CSV** file with the following column format:

> * pair: gene1<==>gene2
> * partner1: gene1
> * partner2: gene2
> * pval_E1: correlation pvalue in experiment 1
> * pval_E2: correlation pvalue in Experiment 2
> * comb_pval: Fisher's combined pvalue across both experiments
> * comb_rho: combined rho coefficient across both experiments
> * comb_FDR: FDR calculated off the combined pvalue
> * partner1InFold: gene1
> * partner1_FoldChange: Fold change of gene1
> * partner2InFold: gene2
> * partner2_FoldChange: Fold change of gene2
> * corr_direction: correlation direction (either -1 or 1)
> * partner1_FC_direction: Fold change direction of gene1 (either -1 or 1)
> * partner2_FC_direction: Fold change direction of gene2 (either -1 or 1)    
> * IfFoldChangeDirectionMatch: Are the previous 2 values identical (1 if yes, -1 if no)
> * PUC: Is the correlation PUC-compliant (are neg-neg or pos-pos correlations +ve and pos-neg correlations -ve?)

Outputs:
* pickled file of the networkx graph, including each edge's abs(rho coefficient)
* a second pickled file that includes the values needed to calculate deviation from expected

## Calculating topological properties
Purpose: Takes the output from import_network_data.py and performs analyses on the network
	
Arguments:
> **Required**
> * Pickled network file (from import_network_data.py)

> **Optional**
> * --frag	-	Flag; Do you want to compute node fragmentation centrality? (Significantly increases run-time)
> * --bibc	-	Flag; Do you want to compute BiBC? (Significantly increases run-time)
> * --bibc_groups	-	choices are 'node_types' or 'modularity'; What to compute BiBC on, either distinct groups or on the two most modular regions of the network
> * --bibc_calc_type	-	choices are 'rbc' or 'bibc'; Would you like to normalize based on amount of nodes in each group (rbc) or not (bibc)?
> * --node_map	-	Required if node_types is specified for --bibc_groups. CSV of nodes and their types (i.e. otu, pheno, gene, etc.)		
> * --node_groups	-	Required if node_types is specified for --bibc_groups. Its the two groups of nodes to calculate BiBC/RBC on
	
Example input command:
> python calc_network_properties.py <pickled network file> --frag --bibc --bibc_groups node_types --bibc_calc_type bibc --node_map /path/to/mapping/file.csv --node_groups gene micro
	
Example node mapping file (only required if BiBC is being calculated based on pre-defined groups (i.e. genes and microbiota):
* The input is a **CSV** file with **no header**, in the following format:
> GAPDH,gene\
> EGFR,gene\
> TP53,gene\
> ASV1,micro\
> ASV2,micro\
> ASV5,micro

Notes: 
* BiBC will not be calculated if the network size is either too small or if two separate groups are identified in the network











