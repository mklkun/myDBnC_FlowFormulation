#include "PARTITION.h"
//#include "kHPPMASTER.h"
//#include "VAREDGE.h"

PARTITION::PARTITION(long *part,long sz,int rhs_p)
{
	pi = part;
	p_sz = sz;
}


PARTITION::~PARTITION()
{
	delete[] pi;
}


double PARTITION::coeff(CUTINFO &c_info, int index) //(ABA_VARIABLE *v)
{
	double coe = 0.0;

//	VAREDGE *ve = dynamic_cast<VAREDGE *>(v);
    lemonEdge e = c_info.G.edgeFromId(index);

//	if(ve != NULL)	//S'il s'agit d'une variable arete...
//	{
		long t = c_info.G.idFromNode( c_info.G.u(e) ) + 1;
		long h = c_info.G.idFromNode( c_info.G.v(e) ) + 1;

		if(pi[t-1] != pi[h-1])
		{
			coe = 1.0;
		}
		else
		{
			coe = 0.0;
		}
//	}
//	else
//	{
//		coe = 0.0;
//	}

	return coe;
}


/*int PARTITION::genRow(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *var, ABA_ROW &row)
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

double PARTITION::slack(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *variables, double *x)
{
	double val=0;
  	EdgeSet_It it;

  	for(it=it_begin;it != it_end;it++)
		val+=x[(*it -1)];

  	return(rhs_-val);
}*/

