#ifndef EXTRA_CUTS_H
#define EXTRA_CUTS_H

#include <deque>

#include "I_graphe.h"
#include "I_digraphe.h"
#include "I_common.h"

//#include "CONSERVATION.h"
//#include "LIAISON.h"
#include "PARTITION.h"
#include "SP_PARTITION.h"
#include "DOUBLE_CUT.h"
#include "TRIPLE_PATH_CUT.h"

#include "COUPE_AGGREGEE.h" /// à finir de corriger dans la classe .cpp si on veut l'utiliser

#define EUC_2D 1
#define GEO 2
#define EXPLICIT_LOWER_ROW 3
#define EXPLICIT_UPPER_ROW 4
#define EXPLICIT_LOWER_DIAG_ROW 5
#define EXPLICIT_UPPER_DIAG_ROW 6
#define ATT 7

enum cut_types { DC, FP, PA, TPC, AG, SPP, P };
enum comparaison_signs { ET, GT, LT, GET, LET};

class extra_cuts
{
    public:
        extra_cuts();

        I_Graph G;
        vector<demande> vectDemande;
        I_digraph **grapheDemande;
        I_Graph G_dem;
        vector<int> source;
		vector<int> dest;
        int n,m;
        int k;

        long tsp_n_nodes,tsp_m_edges;
		double *x_coord;
		double *y_coord;

		bool isTSPGraph;
		int graphType;
		char File_Name[256];

        CUTINFO *c_info;

        virtual ~extra_cuts();

        void init(const char *filename,bool isTSPFile,int type,int k_val,bool demandeConnectee);
        void update(CUTINFO &Nc_info, IloNumArray &solX, IloArray<IloNumArray> &solF);

		int separation(deque < int * > &seperated_cons, int type, FILE *sortie);	/*Effectue la s�paration des contraintes*/
		int separation_globale(deque < int * > &seperated_cons, int type,FILE *sortie);

		simpleEdge *readTsp(char *file_name,int type);
		simpleEdge *readEUC2D(ifstream &fic);

		int AjouteFPartition(deque < int * > &seperated_cons, I_Graph *gr,int k,b_Node *frac_cycle,long cycle_sz,FILE *sortie);
		int AjoutePartition(deque < int * > &seperated_cons, I_Graph *gr,int k,b_Node *frac_cycle,long cycle_sz,FILE *sortie);
		int AjouteDoubleCut(deque < int * > &seperated_cons, I_Graph *gr,int k,FILE *sortie);
		int AjouteTriplePathCut(deque < int * > &seperated_cons, I_Graph *gr,int k,FILE *sortie);
		int AjouteCoupeAggregee(deque < int * > &seperated_cons, I_Graph *gr,int k);
		int AjouteSPPartition(deque < int * > &seperated_cons, I_Graph *gr,int k,b_Node *frac_cycle,long cycle_sz,FILE *sortie);
		int AjouteSPPartition2(deque < int * > &seperated_cons, I_Graph *gr,int k,b_Node *frac_cycle,long cycle_sz,FILE *sortie);

		static int nbIter;

		int nbCoupe;
		int nbDoubleCut;
		int nbTriplePathCut;
		int nbPartition;
		int nbSPPartition;
        int nbCoupeAggregee;

        double demandGraphConnected; /// where this takes its values?!?!?!?!

    protected:

    private:
};

#endif // EXTRA_CUTS_H
