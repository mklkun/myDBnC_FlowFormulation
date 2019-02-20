#ifndef I_GRAPH_H
#define I_GRAPH_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "I_common.h"
#include <iostream>
#include <fstream>
#include <vector>

#define CAP_VALUE_CAP 1
#define CAP_VALUE_X 2
#define CAP_VALUE_1 3
#define CAP_VALUE_0 4
#define CAP_VALUE_0_FOR_FRAC_EDGE 5
#define CAP_VALUE_M_FOR_FRAC_EDGE 6

#define NO_0_EDGES 1
#define ALL_EDGES 2
#define NO_CAP_GREATER_THAN 3

#define EPSILON 1e-4
#define INFINI 1000000000

#define NIL_E (I_Edge *)0
#define NIL_N (I_Node *)0

#define N_TAB 100000
#define M_TAB 2000000
#define N_SUR_2 50000

#define NODE_NOT_MARKED 0
#define NODE_MARKED 1
#define NODE_LIST_TERMINATED 2

#define MARKED 10

typedef struct demande
{
	int source;
	int dest;
	double val;
} demande;


typedef struct I_Node
{
	long 		id;
	struct I_Edge 	*first_edge;	/*Pointer to the head of the adjacency list*/
	bool   		terminal;

	struct I_Node 	*next_gr;	/*Suivant dans le graphe réduit*/
	long 		id_W;		/*Id de l'ensemble Wi auquel appartient le sommet*/
	struct I_Node	*next_W;	/*the following node of this node in the set of nodes to which they belongs to*/
	struct I_Node	*n_sh;		/*pointeur sur le représentant du sommet dans le graphe réduit*/
	struct I_Node	*next_sh;	/*suivant dans la liste des sommets qui sont contractés ensemble*/
	struct I_Edge	*first_edge_gr;	/*premier dans la liste d'incidence du sommet dans le graphe réduit*/

	/*Element annexe pour effectuer*/
	/*d'autres contractions*/
	struct I_Node 	*a_next_gr;	/*Suivant dans le graphe réduit*/
	struct I_Node	*a_n_sh;		/*pointeur sur le représentant du sommet dans le graphe réduit*/
	struct I_Node	*a_next_sh;	/*suivant dans la liste des sommets qui sont contractés ensemble*/

	/*Pour un deuxi�me niveau de sauvegarde*/
	struct I_Node 	*b_next_gr;	/*Suivant dans le graphe réduit*/
	struct I_Node	*b_n_sh;	/*pointeur sur le représentant du sommet dans le graphe réduit*/
	struct I_Node	*b_next_sh;	/*suivant dans la liste des sommets qui sont contractés ensemble*/
} I_Node;

typedef struct I_Edge
{
	struct I_Node	*adjac;	/*The destination node to which this edge is adjacent*/
	struct I_Edge	*next;	/*the following edge in the adjacency list to wich this edge belongs to*/
	struct I_Edge 	*back;	/*reverse edge*/
	double 		cap;	/*Capacity of the edge*/

	float 		X;	/*Value of X(e)*/
	long 		num;	/*number of this edge in the graph edge list*/

} I_Edge;


typedef struct I_Graph
{
	long n_Nodes;
	long m_Edges;
	I_Node *Nodes;
	I_Edge *Edges;

	Bool mark_tab[N_TAB];

} I_Graph;


bool AllocateGraph(long n,long m,I_Graph *g);	/*Allocates memory for the internal graph representation*/
bool DeleteGraph(I_Graph *g);			/*Deletes graph from memory*/
bool InitGraph_from_file(char *fich,I_Graph *g,int *k);	/*Reads graph from file input*/
bool InitGraph_from_list(simpleEdge *Te,long n,long m,I_Graph *g);	/*Reads graph from input list*/
bool InitGraph_from_list(vector<simpleEdge> &Te,long n,I_Graph *gr);
bool copy_Graph(I_Graph *gr,I_Graph *Ogr);
void PrintGraph(I_Graph *g);			/*Prints the graph*/
//void PrintReducedGraph(I_Graph *g,FILE *fout);	/*Print the shrunk graph*/
//void PrintGraph_Table(I_Graph *g);

void operationTheta(I_Graph *gr,int k);

/*Permet de sauvegarder la contraction courante*/
/*dans les champs annexes*/
void save_contraction_info(I_Graph *gr,int niveau);
void recall_contraction_info(I_Graph *gr,int niveau);

/*Retourne Vrai si l'entier z est pair et Faux sinon*/
Bool is_pair(long z);

long sup_part(long num,long denom);

//Verifie si l'ensemble de sommet W induit une coupe valide
//c-a-d qu'il existe une demande (s,t) avec s dans W et t dans V\W
//Pour cela, il suffit de verifier si |\delta_D(W)| >= 1
Bool admissible(I_Graph *gr_dem,b_Node *W);

//Verifie si une partition est admissible
Bool admissible(I_Graph *gr_dem,b_Node *partition[],int nbElem);

/*Gives the list of the edges of delta(W).*/
/*The method is similar to the method used in get_edge_list().*/
b_Edge *get_delta_w(I_Graph *g,b_Node *w,int edge_flag,int cap_flag);

//Renvoie les aretes de [u,v]
b_Edge *get_edge_u_v(I_Graph *g,long u,long v,int edge_flag,int cap_flag);

/*Contract the set of nodes contained in the linked list f_w*/
Bool contract_set_w(I_Graph *gr,b_Node *f_w);

void get_sp_partition_chaine_steiner(I_Graph *gr,b_Node **chaine,long *chaine_sz,long *nb_chaine);

Bool separation_sp_partition_2(I_Graph *gr,int k_ordre,b_Node *chaine[],long *chaine_sz,long nb_chaine,long
***sp_part_list,long *sp_p,long *nb_sp_part,long *rhs);

Bool separation_sp_partition_steiner(I_Graph *gr,int k_ordre,b_Node *chaine[],long *chaine_sz,
long nb_chaine,long ***sp_part_list,long *sp_p,long *nb_sp_part,long *rhs);

Bool separation_partition_steiner(I_Graph *gr,int k,b_Node *frac_cycle,long cycle_sz,long **partition,long *rhs,FILE *sortie_frac_gk);

/*Gives the list of the edges of delta(v).*/
/*The method is similar to the method used in get_delta_w().*/
b_Edge *get_delta_v(I_Graph *g,long v,int edge_flag,int cap_flag);

int trouverMinimum(double *L,int *T,int n);
void Dijkstra(I_Graph *gr,int s,double *&L,int *&P,int *&PredArete,int cap_flag);

Bool separation_double_cut(I_Graph *gr,b_Node *pi[],int k,int s,int t,int *&partition,int &F,int &rhs,FILE *sortie_frac_gk);
Bool separation_double_cut(I_Graph *gr,b_Node *pi[],int k,int s,int t,int *&partition,vector<int> &F,int &rhs,FILE *sortie_frac_gk);

Bool separation_triple_path_cut(I_Graph *gr,b_Node *pi[],int k,int s,int t1,int t2,int *&partition,vector<int> &F,int &rhs,FILE *sortie_frac_gk);

/*Merges two cyclic linked lists*/
/*The lists must have at least 1 node.*/
Bool cyclic_set_merge(I_Node *f1,I_Node *f2);

/*Checks if the cardinal of the set W is greater or equal to a.*/
/*You can use this function to check if |W| < b by checking if*/
/*the result of check_w_card(W,b) is false.*/
Bool check_w_card(b_Node *w,long a);

/*Checks if the cardinal of the set W is greater or equal to a.*/
/*You can use this function to check if |W| < b by checking if*/
/*the result of check_w_card(W,b) is false.*/
Bool check_delta_w_card(b_Edge *d_w,long a);

/*Calcule le plus court chemin entre*/
/*tous les couples de sommets dans le graphe G*/
/*On utilise l'algorithme de Floyd*/
long **Floyd(b_Edge *sh_eList,long n);

/*Renvoie la liste des ar�tes du graphe r�duit sans les ar�tes nulles*/
simpleEdge *GetReducedGraph_eList(I_Graph *g,long *m_edge);

bool InitGraph_from_vect(vector<demande> *vdem,long n,long m,I_Graph *gr);

void DrawGraphViz(I_Graph *g,char *nomfich,int parN,vector<simpleEdge> *label,int nbDem);

#endif




