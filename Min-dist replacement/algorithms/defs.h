// defs.h
// Commonly used definitions
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 21/09/2012
//		Last Modified: 24/09/2012

#ifndef __DEFS_H
#define __DEFS_H

// dimensionality
#define DIMENSION 2

// files for the rtree data
#define CLIENT_RTREE "C.rtree"
#define POTENTIAL_LOCATION_RTREE "P.rtree"
#define FACILITY_RTREE "F.rtree"

#define CLIENT_RNN_TREE "C.rnntree"

#define CLIENT_MND_TREE "C_mnd.rtree"
#define CLIENT_RNN_TREE_WITH_DNN "C_dnn.rnntree"

#define CLIENT_BB_TREE "C.dnn_artree"

// file for the replacement query
#define CLIENT_MND_REPQUERY_TREE "C_2mnd.rtree"

#define LOCATION_BB_REPQUERY_TREE "P_is.rtree"
#define FACILITY_BB_REPQUERY_TREE "F_is.rtree"

// files for the sequential array data
#define CLIENT_ARRAY_FILE "C.array"
#define POTENTIAL_LOCATION_ARRAY_FILE "P.array"
#define FACILITY_ARRAY_FILE "F.array"

#define CLIENT_ARRAY_FILE_WITH_DNN "C_dnn.array"

// file for the replacement query
#define CLIENT_ARRAY_FILE_FOR_D2NN_COMPUTATION "C_d2nn_comp.array"
#define CLIENT_ARRAY_FILE_WITH_D2NN "C_d2nn.array"

// used in readInput()
#define MAXLEN 100

// number of page counters
#define PAGE_COUNTER_NUM 3

// number of cache blocks
#define CACHE_BLOCK_NUM 10

// data range
#define DATA_DOMAIN 1000

#endif