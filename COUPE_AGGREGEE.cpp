#include "COUPE_AGGREGEE.h"
//#include "VAREDGE.h"
//#include "VARARC.h"

COUPE_AGGREGEE::COUPE_AGGREGEE(int _n,int
rhs_ca,vector<list<I_dinode*> > *_lw,vector<list<I_diarc*> >
*_lentier,vector<int> *_ldemande,list<I_Edge*> *_lx,list<I_diarc*> *_ly0,list<I_diarc *> *_ly1)
{
	LW.resize(_lw->size());
	LEntier.resize(_lentier->size());
	LDemande.resize(_ldemande->size());

	for(int i=0;i<(int)_lw->size();i++)
	{
		LW[i] = (*_lw)[i];
		LEntier[i] = (*_lentier)[i];
		LDemande[i] = (*_ldemande)[i];
	}

	Lx = (*_lx);
	Ly0 = (*_ly0);
	Ly1 = (*_ly1);

	n = _n;
}


COUPE_AGGREGEE::~COUPE_AGGREGEE()
{}


double COUPE_AGGREGEE::coeff(CUTINFO &c_info, int index)
{
	double coe = 0.0;

//	VAREDGE *ve = dynamic_cast<VAREDGE *>(v);
    lemonEdge e = c_info.G.edgeFromId(index);

	if(true)   ///(ve != NULL)	//S'il s'agit d'une variable arete...
	{
		//On regarde si l'arete n'est pas dans Lx
		list<I_Edge*>::iterator itx = Lx.begin();
		while(itx != Lx.end() && (*itx)->num != index+1 ) ///ve->getEptr()->num)
		{
		    //cout << "Lx == " << (*itx)->adjac->id << " " << (*itx)->back->adjac->id << endl;

            itx++;
		}

		if(itx != Lx.end())
			coe = 1.0;
		else
			coe = 0.0;
	}
	else
	{
	    /// Well, If we will use AG cuts we should correct this part too
//		VARARC *va = dynamic_cast<VARARC *>(v);
//
//		//Premierement, on regarde si l'arc est dans Ly1
//		list<diarc*>::iterator it1 = Ly1.begin();
//		while(it1 != Ly1.end() && *it1 != va->getAptr())
//			it1++;
//
//		if(it1 != Ly1.end())
//			coe = 1.0;
//		else	//sinon, on regarde ensuite dans Ly0
//		{
//			list<diarc*>::iterator it0 = Ly0.begin();
//			while(it0 != Ly0.end() && *it0 != va->getAptr())
//				it0++;
//
//			if(it0 != Ly0.end())
//				coe = 0.0;
//			else
//			{
//				//On regarde ensuite si l'arc n'est pas dans LEntier
//				bool arcEntierTrouve = false;
//				for(int i=0;i<(int)LEntier.size() && !arcEntierTrouve;i++)
//				{
//					list<diarc *>::iterator ite = LEntier[i].begin();
//
//					while(ite != LEntier[i].end() && *ite != va->getAptr())
//						ite++;
//
//					if(ite != LEntier[i].end())
//						arcEntierTrouve = true;
//				}
//
//				if(arcEntierTrouve)
//					coe = 0.0;
//				else
//				{
//					vector<bool> marquage(2*n);
//					list<dinode*>::iterator itLW;
//
//					int p = 0;
//					for(int i=0;i<(int)LW.size();i++)
//					{
//						if(va->getNumDemande() == LDemande[i])
//						{
//							for(int j=0;j<2*n;j++)
//								marquage[j] = false;
//
//							//On commence par marquer les sommets de la coupe
//							for(itLW=LW[i].begin();itLW!=LW[i].end();itLW++)
//							{
//								marquage[(*itLW)->num_reel-1] = true;
//							}
//
//							if(marquage[va->getAptr()->source->num_reel-1] && !marquage[va->getAptr()->dest->num_reel-1])
//							{
//								p++;
//							}
//							else
//							{
//								//cout << "NONONO" << endl;
//							}
//						}
//
//					}
//					coe = ceil((double)p/2.0);
//				}
//			}
//		}
	}

	//cout << "COE == " << coe << endl;

	return coe;
}


/*int COUPE_AGGREGEE::genRow(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *var, ABA_ROW &row)
{
	ABA_ARRAY<int> s(master_,Part_List.size());
  	ABA_ARRAY<double> c(master_,Part_List.size());
  	ABA_CSENSE sense(master_,'G');
  	int i;
  	EdgeSet_It it;

  	i=0;
  	it = it_begin;
  	while(it!=it_end)
    	{
      		s[i]=(*it -1);
		c[i]=1.0;
      		i++;
      		it++;
    	}

  	ABA_ROW row2(((ABA_GLOBAL*)master_),((int) Part_List.size()),s,c,sense,rhs_);

  	row.copy(row2);

  	return((int)Part_List.size());
}

double COUPE_AGGREGEE::slack(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *variables, double *x)
{
	double val=0;
  	EdgeSet_It it;

  	for(it=it_begin;it != it_end;it++)
		val+=x[(*it -1)];

  	return(rhs_-val);
}*/

