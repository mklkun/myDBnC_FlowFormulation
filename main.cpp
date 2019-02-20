/*
 *		Main.cpp
 *		Created on: 31/08/2017
 *		Author: Mohamed Khalil Labidi
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <iomanip> // for precision

using namespace std;

#include "BC_FORM_EXT.h"
#include "global_functions.h"
#include "global_variables.h"
#include "control_panel.h"


int main(int argc, char** argv)
{
    ////////////////////////////////////////////
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    ////////////////////////////////////////////

    ////////////////////////////////////////////
    // Find out rank, size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    ////////////////////////////////////////////

    if ( argc == 6 )
    {
        ////////////////////////////////////////////
        pathNameI = "Instances/";

        name =  argv[1];
        TIME_LIMIT = atof(argv[2]);

        k = atof(argv[3]);;
        L = atof(argv[4]);;

        nb_threads = atof(argv[5]);
        ////////////////////////////////////////////

        ////////////////////////////////////////////
        MPI_MyBCinit();
        ////////////////////////////////////////////

        ////////////////////////////////////////////
        if (world_rank == 0)
        {
            cout <<         "---------------------------------------------" << endl ;
            //	cout << "Instance\t->" << argv[1] << endl;
            cout << "Instance\t-> " << name << endl;
            cout <<         "---------------------------------------------" << endl ;

            cout << "Time limit\t->" << TIME_LIMIT << endl;
        }
        ////////////////////////////////////////////

        ////////////////////////////////////////////
        kHNDP_instance_read(name.c_str());
        //	kHNDP_instance_read(argv[1]);
        ////////////////////////////////////////////

        ////////////////////////////////////////////
        master_kHNDP_initialize();
        ////////////////////////////////////////////

        ////////////////////////////////////////////
        if (world_rank == 0)
            master();
        else if (world_rank == 1)
            manager();
        else
            slave();
        ////////////////////////////////////////////

        ////////////////////////////////////////////
        master_kHNDP_free();
        ////////////////////////////////////////////

        ////////////////////////////////////////////
        if (world_rank == 0)
            master_kHNDP_Output();
        ////////////////////////////////////////////

        ////////////////////////////////////////////
        kHNDP_instance_free();
        ////////////////////////////////////////////

        ////////////////////////////////////////////
        MPI_MyBCfree();
        ////////////////////////////////////////////

    }
    else
    {
        if (world_rank == 0)
        {
            cout << "ERROR Invalid Number of Parameters " << endl;
            cout << "  usage: mpirun -np NB_PROCS ./FILENAME INSTANCE TIME_LIMIT K L NB_THREADS" << endl;
        }
    }

    ////////////////////////////////////////////
    //Finalizing the MPI environment
    MPI_Finalize();
    ////////////////////////////////////////////

    return 1;
}
