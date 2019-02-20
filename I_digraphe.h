#ifndef I_DIGRAPHE_H
#define I_DIGRAPHE_H

#include <iostream>
#include <vector>
#include <list>
#include <stdlib.h>
#include "I_graphe.h"
#include "I_graphe_flot.h"

#define ALL_ARCS 1
#define NO_0_ARCS 2

//Structures de donnee pour les graphes oriente de demandes

typedef struct I_dinode
{
	int 		num;		//numero du sommet dans le graphe (de 0 a |\tilde{V}|-1)
	int 		num_reel;	//numero reel du sommet dans le graphe non oriente
	struct I_diarc 	*premier;	//premier dans la liste d'adjacence du sommet
} I_dinode;

typedef struct I_diarc
{
	int 		num;		//numero d'ordre de l'arc
	I_dinode 		*source;	//pointeur sur noeud origine
	I_dinode 		*dest;		//pointeur sur le noeud destination
	float 		f;		    //valeur de flot
	float 		cap;		//Capacite
	I_Edge 		*pedge;		//pointeur sur l'arete correspondant dans le
                            //graphe non oriente

	struct I_diarc	*suivant;	//suivant dans la liste d'adjacence du sommet
}I_diarc;

#define NIL_A (I_diarc *) 0

bool operator==(I_dinode &n1,I_dinode &n2);
bool operator<(I_dinode &n1,I_dinode &n2);
bool operator>(I_dinode &n1,I_dinode &n2);
ostream &operator<<(ostream &os,I_dinode &n1);
ostream &operator<<(ostream &os,I_diarc &a);


class I_digraph
{
  public:
	int n;			//Nombre de sommet du graphe non oriente:
                    //On a : |\tilde{V}| = 2(n-2)
	int m;			//Nombre d'aretes de \tilde{G}

	I_dinode source;		//source de la demande
	I_dinode dest;		//destination de la demande
				//Dans le graphe oriente, la source et la dest
				//seront numero 1 et numero 2
				//On suppose aussi que s <= t

	I_dinode *V1;
	I_dinode *V2;
	I_diarc *larc;

	/*vector<I_dinode> V1;	//liste des sommets de N
	vector<I_dinode> V2;	//liste des sommets de N'

	vector<I_diarc> larc;	//liste des arcs du graphe
				//Remarque: on sait que tous les arcs qui sont
				//differents des arcs (s,u) ou (v',t) vont de
				//V1 vers V2. Donc, a partir d'un arc, on sait
				//que l'origine est dans V1 et la dest est dans V2
	*/
	/****************** Methodes *********************/

	//********************************
	//Constructeur a partir d'une liste d'une
	//liste d'arc fournie par simpleEdge

	I_digraph(int s,int t,I_Graph *g,int k=3);


	//********************************
	//Destructeur
	~I_digraph();


	//**************************************
	//Renvoie le nombre de sommets du graphe
	//Renvoie le nombre d'arcs du graphe

	int getNbNode();
	int getNbArc();




	//**************************************************
	//Ajoute l'arc a dans la liste des voisins de sommet

	bool lierArc(I_dinode &sommet,I_diarc *a);


	//********************************************************
	//Ajoute dans le graphe les arcs correpondant a l'arete e

	//bool addArc(I_dinode &sommet,I_diarc &a);


	//****************************************
	//Afiche les arcs du graphe oriente

	void print(ostream &os);
	void printVoisin(ostream &os,I_dinode &n1);
	void printTable(ostream &os);

	void graphViz(const char *nomfich,bool avec_val);
	void graphViz(const char *nomfich,simpleEdge *solution,bool avec_val);

	void getDelta(list<I_dinode *> &W,list<I_diarc *> &delta_w);
	void getDeltaPlusS(list<I_diarc *> &delta_s,int edge_flag);
	void getDeltaMoinsT(list<I_diarc *> &delta_t,int edge_flag);

    double makeFeasible(I_Graph *g,int k);
    void directInitGraphFlot(I_Graph_Flot *gr,int nodes);
    void tranfertFlot(I_Graph_Flot *gr);
    double calculFlot(list<I_dinode *> &W);
    bool arcLier2(I_diarc *a1,I_diarc *a2,int k,list<I_dinode *> &coupeW,list<I_diarc *> &lEntier);
};


/****************************************************/
//Recherche la position du sommet u dans la liste a
//entre les positions dep et arr
//Comme la liste est triee, on fera une recherche dichotomique

template<class T>
int rechercheDichotomique(T *L, T &u,int dep,int arr)
{
	int a,b,milieu;

	a = dep;
	b = arr;

	while(a <= b)
	{
		milieu = (a+b)/2;

		if(L[milieu] == u)
			return milieu;
		else if(u < L[milieu])
			b = milieu-1;
		else
			a = milieu+1;
	}

	return -1;
}

bool corresponde(I_diarc *a1,I_diarc *a2,int n);

#endif
