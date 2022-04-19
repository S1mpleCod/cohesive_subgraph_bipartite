# cohesive_subgraph_bipartite
This repository contains a maximal biclique enumeration algorithm ooMBEA. 
ooMBEA implements techniques proposed in Maximal Biclique Enumeration for Large Sparse Bipartite Graphs (PVLDB 22), 
which inlcude order optimisation, batch-pivots technique and array-based data strcutre. 


Utilities
Different formats of input are supproted such as adjlist (right and left), edge list or adjlist. Please see LCUtility.h/.cpp.
Please refine delimiters according to datasets. 

Data strctures
Adjacency hashset and the proposed data strcutre are implemented. Please load your graphs accrodingly.  Please see bigraph.h/.cpp.

Algorithms
Please see BiCore.h/.cpp for the unilateral core related algorithms
Please see ooMBE.h/.cpp for the main algorithm

