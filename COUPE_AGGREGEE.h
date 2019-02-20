#ifndef COUPE_AGGREGEE_H
#define COUPE_AGGREGEE_H

#include "I_graphe.h"
#include "I_digraphe.h"
#include "I_common.h"

#include "extra_cuts.h"

#include <list>
#include <vector>


class COUPE_AGGREGEE
{
	private:
		vector<list<I_dinode *> > LW;		    //Liste des ensembles W
		vector<list<I_diarc *> >  LEntier;	//Liste des ensembles d'arcs a 1
        vector<int>             LDemande;   //Liste des demandes correspondant à W et LEntier


		list<I_Edge*> Lx;	//Liste des variables x
		list<I_diarc*> Ly0;	//liste des variables y a 0
		list<I_diarc*> Ly1;	//liste des variables y a 1

		int n;			//nombre de sommets dans le graphe

	public:
		COUPE_AGGREGEE(int _n,int rhs_ca,vector<list<I_dinode*> > *_lw,
		vector<list<I_diarc*> > *_lentier, vector<int> *_ldemande,list<I_Edge*> *_lx,
		list<I_diarc*> *_ly0,list<I_diarc *> *_ly1);

		virtual ~COUPE_AGGREGEE();
		virtual double coeff(CUTINFO &c_info, int index);

		/*virtual int genRow(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *var,ABA_ROW &row);
  		virtual double slack(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *variables, double *x);*/
};

#endif
