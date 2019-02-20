#ifndef VARIABLE_local_HEADER
#define VARIABLE_local_HEADER

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>

#include "graph.h"
#include "digraph.h"
#include "solution.h"
#include "extra_cuts.h"



using namespace std;



/////////////////////////////kHNDP DATA////////////////////////////////////////
extern string pathNameI;
extern string name;
extern unsigned k, L, demandsNumber;
extern graph G;
extern vector < vector < unsigned > > demands;
extern vector < digraph* > DGS;
extern vector < unsigned > linkX;
extern vector < vector < unsigned > > linkF;
extern extra_cuts myextra;
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////MODEL INFOS///////////////////////////////////////
extern unsigned NBCOLS;
extern unsigned nb_threads;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//BEST FITS
extern int UB;
extern double LB;
extern solution *sol;
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////MPI WORLD/////////////////////////////////////////
extern int world_rank;
extern int world_size;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//BEST FIT
//extern int binNumber;
//extern int *vector_item_bin;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//INPUT
extern double TIME_LIMIT;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//COMPUTATION TIME
extern clock_t time_start,time_finish,time_current;
extern double computation_time,current_time;
///////////////////////////////////////////////////////////////////////////////

#endif
