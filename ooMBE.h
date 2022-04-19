//
// Created by Lu Chen on 5/07/21.
//

#ifndef MBBP_OOMBE_H
#define MBBP_OOMBE_H

#include "bgraph.h"
#include "BiCore.h"
#include <utility>
#include <functional>
#include <chrono>
using namespace std::chrono;
/**
 * Assume the vertex set to be enumerated has been decided
 * bgraph provides a function for swapping L and R
 */
struct PivotsNeighbour{ //local neighbour
    int vertex;
    vector<int> adj;
    vector<int> Q;
    vector<int> agg;
    int locAdj;
};
class ooMBE {
public:
    int recursive_tracker=0;
    int des_times=0;
    void iniVector(vector<vector<int>>&,vector<int> &,vector<int> & ,vector<int> &);
    void iniArrayValues(vector<vector<int>>& adjR, int *R, int *indexR, int *adjIndexR);
    bgraph &g;
    //bgraph g;
    int *counts= 0;
    bool maxPruning= false;
    vector<BFST> g_bfs;
    vector<bool> isPovit;
    ooMBE(bgraph &bg):g(bg){  nomb = 0;};
    int max_degree=0;

    //ooMBE(bgraph bg){  g=bg;nomb = 0;};
    //ooMBE(){nomb = 0;};
    //vector<int> R,indexR; // always works on R, a mapping shall be there in the MBE function
    //vector<int> L,indexL;
    //vector<int> adjIndexR,adjIndexL;
    int *R= 0, *indexR= 0, *L= 0, *indexL= 0, *adjIndexR= 0, *adjIndexL= 0, size_of_R=0, size_of_L=0;
    unsigned long int nomb;
    int *ln_for_vex;
    int *other_for_vex;
    int *two_hop;
    //pair<int,int> *pivots_arry;
    //this is the iMBEA with cached local subgraph and space optimisation
    void mbeStart();
    void mbeDSini(vector<vector<int>> &,vector<vector<int>> &);     //void mbeDSiniArray(vector<vector<int>> &,int*,int );
    void advIMBEA(int,int,int,int,vector<pair<int,int>>&,vector<pair<int,int>>&,vector<pair<int,int>>&);
    void advPIMBEA(int,int,int,int,vector<pair<int,int>>&,vector<pair<int,int>>&,vector<pair<int,int>>&);
    void advPIMBEA_local(int,int,int,int,vector<pair<int,int>>&,vector<pair<int,int>>&,vector<pair<int,int>>&);

    void adv_order_PIMBEA_local(int,int,int,int,vector<pair<int,int>>&,vector<pair<int,int>>&,vector<pair<int,int>>&);//under developing
    int max2Core =0;
    void adv_mbeStart();

    void adv_mbeStart_reuse_full();//under developing
    void adv_mbeStart_reuse(); //needs slightly more memory
    //void ooMBE::advIMBEA(int iA, int iQ, int iP, int iB, vector<pair<int,int>> &A,vector<pair<int,int>> &P,vector<pair<int,int>> &Q)
    void adjustment_after_moving_a_vertex(int u,
                                          int* L,
                                          int* indexL,
                                          int* adjIndexL,
                                          int size_of_L,
                                          int* R,
                                          int* indexR,
                                          int* adjIndexR,
                                          int size_of_R
                                          );
    void vertex_swap(int loc1,
                     int loc2,
                     int* array,
                     int* index
                     );


    void printVector(vector<int>);
    bool orderSort(int v, int v_prime) {
        return (g.adjR[v].size()>g.adjR[v_prime].size());
    };
    void orderBatchPivots(vector<int>&);
    virtual ~ooMBE();

    void orderBatchPivots_neighbours(vector<PivotsNeighbour> &);
};

struct Pivot{
    pair<int,int> p; //first vertex, second loc of local adjlist
    int loc=-1;
};

struct L_N{ //local neighbour
    int vertex;
    int loc_in_adj;
};



#endif //MBBP_OOMBE_H
