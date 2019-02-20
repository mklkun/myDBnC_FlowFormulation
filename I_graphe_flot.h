#ifndef I_GRAPHE_FLOT_H
#define I_GRAPHE_FLOT_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <list>
#include <cstdlib>
#include "I_common.h"


using namespace std;

//#define DFS          1
#define BFS          2
#define MAX_GAIN     3
#define DINIC        4
#define DINIC_NEW    5
#define KARZANOV     6
#define GOLDBERG_1   7
#define GOLDBERG_2   8
#define GOLDBERG_3   9
#define GOLDBERG_4  10
#define GOLDBERG_5  11

#define INFINIT 10000000

#define FAILURE    0
#define SUCCESS    1
#define FALSE      0
#define TRUE       1


//#define MAX_N     20000
#define MAX_N 100

#define MAX_CAP   100000000

/* Dimacs problem types */
#define UNDEFINED        0
#define MINCOSTFLOW      1
#define MAXFLOW          2
#define ASSIGNMENT       3


typedef struct enode {
  struct enode *next;
  struct enode *mate;
  double c;
  double c_inf;
  double f;
  int h;
  int t;
  int flag;

  /************************************************/
  //int label;		//Label: ajouter pour connaitre les chaines augmentantes
  int numDiarc;		//numero du diarc correspondant

} I_Edge_Flot;


typedef struct {
  int head, tail, size;
  int *data;
} Queue;

typedef struct {
  int ptr, size;
  int *data;
} Stack;


class I_Graph_Flot
{
  public:

  	/*I_Edge_Flot *A[MAX_N];
  	int V[MAX_N];*/

	I_Edge_Flot **A;
  	int *V;
  	int size;
  	int max_v;
  	int edge_count;

	/*int nbLablel;	//Nombre de label utilises pour un flot
	int parLabel;	//Parametre utilise pour stocker les labels
			//on aura label(i) = label(i-1)*q + i;*/

  	/*Pour le flot*/
  	double *Excess;
  	int *Dist;
  	int RCount, UCount, SCount;
  	I_Edge_Flot **Current;

	/**************/
	~I_Graph_Flot();

	/**************/

	void OutputFlow(const char *fich,double s);
	void WriteVertex(int v,ofstream &f);
	void WriteVertex2(int v,ofstream &f);
	void WriteVertex3(int v,ofstream &f);
	I_Edge_Flot *EdgeLookup(int v1,int v2);
	void AddEdge(int v1,int v2,double a_inf,double a_sup,int corresp=-1);
	void AddVertex(int v);
	void InitGraph(int n);
	bool InputFlowGraph(char *fich,int &s,int &t);
	void InputFlowGraph(list<simpleEdge> &Te,int nodes,int &s,int &t);
	void GraphViz(char *fichier);
	void GraphViz(int s, int t,char *fichier);
	void DeleteGraph_Flot();

	//Fonction qui cherche un chemin augmentant dans G entre s et t
	//le parametre marquer donne les sommets qui sont marqués et ceux
	//qui ne le sont pas
	double AugmenteFlot(int s, int t,bool *&marquer,bool avecLabel=false);
	void CorrectionFlot();

	/**********************************************/
	/**********   Calcul de flot avec Goldberg   **/

	double FindFlow(int s,int t);
	double VertexFlow(int i);
	void ValidFlow(int s,int t);
	void MarkCut(int u,int C[],int *n);
	void PrintCut(int s);
	void InitFlow();
	void Goldberg(int s,int t);
	void Goldberg1(int s,int t);
	void InitGoldberg(int s);
	void SetLabels(int v1,int v2);
	void Discharge(int v,Queue *Q,int s,int t);
	void PushRelabel(int v,Queue *Q,int s,int t);
	void Relabel(int v);
	void Push(I_Edge_Flot *e);
	void EndGoldberg();

	double flotMaxAugmentation(int s,int t,bool *&marquer);
};

void Barf(const char *s);
char *Alloc(int n);
double Min(double x,double y);
double Max(double x,double y);
double Abs(double x);
void InitRandom (int seed);
Queue *MakeQueue(int n);
int Dequeue(Queue *Q);
void Enqueue(Queue *Q,int k);
bool QueueEmpty(Queue *Q);

#endif

