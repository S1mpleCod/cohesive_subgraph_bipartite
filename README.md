# Cohesive_Subgraph_Bipartite
This repository contains a maximal biclique enumeration algorithm ooMBEA. 
ooMBEA implements techniques proposed in Maximal Biclique Enumeration for Large Sparse Bipartite Graphs (PVLDB 22), 
which inlcude order optimisation, batch-pivots technique and array-based data strcutre. 


## Utilities
Different formats of input are supproted such as adjlist (right and left), edge list or adjlist. Please see LCUtility.h/.cpp.
Please refine delimiters according to datasets. 

## Data strctures
Adjacency hashset and the proposed data strcutre are implemented. Please load your graphs accrodingly.  Please see bigraph.h/.cpp.

## Algorithms
Please see BiCore.h/.cpp for the unilateral core related algorithms
Please see ooMBE.h/.cpp for the main algorithm

## How to reproduce major results

- In main.cpp, please refine the parameter in bg.readEdgelist("---") and replace --- with an edge list of a dataset.

- Build the entire project using the provided CMake scripts using CMake, and it will generate the requried makefile and compile the code to a runable file. The compiling parameters are configured in CMkaeLists.txt.

- Usage:
  
  ```shell
  cmake CMakeList.txt
  make
  ./mbbp <edgelist_path>
  ```

- Edge list format:
    - The ID should start from 1.
    - The left and right set's IDs cannot be duplicate.
    - One line (end with \t) represents one edge pair


- Run the runnable file, and you will get results following the following format: 
    
    call adv_mbeStart() 
    
    switched!                            //Which vertex set is enumerated 
    
    Running time in milliseconds: 14327  //total running time excluding loading the data into the proposed data structure 
    
    Number of bicliques: 4899032         // number of maximal bicliques
