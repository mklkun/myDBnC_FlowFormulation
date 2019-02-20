

#include "global_variables.h"


////////////////////////////////GRAPHS/////////////////////////////////////////
string pathNameI;
string name;
unsigned k, L, demandsNumber;
graph G("");
vector < vector < unsigned > > demands;
vector < digraph* > DGS;
vector < unsigned > linkX;
vector < vector < unsigned > > linkF;
extra_cuts myextra;
///////////////////////////////////////////////////////////////////////////////

//////////////////////////////MODEL INFOS//////////////////////////////////////
unsigned NBCOLS(0);
unsigned nb_threads;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//BEST FITS
int UB;
double LB;
solution *sol;
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////MPI WORLD/////////////////////////////////////////
int world_rank;
int world_size;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//BEST FIT
//int binNumber;
//int *vector_item_bin;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//INPUT
double TIME_LIMIT;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//COMPUTATION TIME
clock_t time_start,time_finish,time_current;
double computation_time,current_time;
///////////////////////////////////////////////////////////////////////////////