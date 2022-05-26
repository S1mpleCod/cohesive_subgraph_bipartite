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

How to reproduce major results
A) In main.cpp, please refine the parameter in bg.readEdgelist("---") and replace --- with an edge list of a dataset. 

  Notices 
  1) The default edge list supported should be: 1)vertx id is from 1, 2)directly start with an edge, 3) the delimiter shoud be and empty space
  between two node, and 3) each line contains one single edge.
  For instance:
  1 2
  2 1
  3 2
  2) The code will make ids of vertices staring from 0. 
  3) Change readEdgeList(-) in bgraph.cpp if you have different delimiters or ids.

B) Build the entire project using the provided cmake scripts using CMAKE and it will generate requried makefile and compile the code to a runable file. The compiling parameters are configured in CMkaeLists.txt.


C) Run the runable file and you will get results following the following format:
  call adv_mbeStart()
  switched!                                   //Which vertex set is enumerated 
  Running time in milliseconds: 14327         //total runing time excluding loading the data into the proposed data structure 
  Number of bicliques: 4899032                // number of maximal bicliques




