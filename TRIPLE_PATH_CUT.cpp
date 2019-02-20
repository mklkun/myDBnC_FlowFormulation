#include "TRIPLE_PATH_CUT.h"
//#include "kHPPMASTER.h"
//#include "VAREDGE.h"

TRIPLE_PATH_CUT::TRIPLE_PATH_CUT(int *part,vector<int> &f_s,int rhs_p)
{
	pi = part;
	F = f_s;
}


TRIPLE_PATH_CUT::~TRIPLE_PATH_CUT()
{
	delete[] pi;
}


double TRIPLE_PATH_CUT::coeff(CUTINFO &c_info, int index) //(ABA_VARIABLE *v)
{
	double coe = 0.0;

//	VAREDGE *ve = dynamic_cast<VAREDGE *>(v);
    lemonEdge e = c_info.G.edgeFromId(index);

//	if(e != NULL)	//S'il s'agit d'une variable arete...
//	{
		int t = c_info.G.idFromNode( c_info.G.u(e) ) + 1;
		int h = c_info.G.idFromNode( c_info.G.v(e) ) + 1;

		//cout << "T == " << t << " H == " << h << endl;

		int j = 0;
		while(j < F.size() && F[j] != index + 1)
			j++;

		if(j < F.size())
			coe = 0;
		else if( (pi[t-1] == 3 && pi[h-1] == 4) || (pi[t-1] == 4 && pi[h-1] == 3) )
			coe = 1;
		else if( (pi[t-1] == 2 && pi[h-1] == 4) || (pi[t-1] == 4 && pi[h-1] == 2) )
			coe = 1;
		else if( (pi[t-1] == 3 && pi[h-1] == 5) || (pi[t-1] == 5 && pi[h-1] == 3) )
			coe = 1;
		else if( (pi[t-1] == 0 && pi[h-1] == 2) || (pi[t-1] == 2 && pi[h-1] == 0) )
			coe = 2;
		else if( (pi[t-1] == 0 && pi[h-1] == 3) || (pi[t-1] == 3 && pi[h-1] == 0) )
			coe = 2;
		else if( (pi[t-1] == 1 && pi[h-1] == 3) || (pi[t-1] == 3 && pi[h-1] == 1) )
			coe = 2;
		else if( (pi[t-1] == 0 && pi[h-1] == 5) || (pi[t-1] == 5 && pi[h-1] == 0))
			coe = 1;
		else if( (pi[t-1] == 0 && pi[h-1] == 4) || (pi[t-1] == 4 && pi[h-1] == 0))
			coe = 1;
		else if( (pi[t-1] == 1 && pi[h-1] == 4) || (pi[t-1] == 4 && pi[h-1] == 1))
			coe = 1;
		else if( (pi[t-1] == 1 && pi[h-1] == 5) || (pi[t-1] == 5 && pi[h-1] == 1) )
			coe = 1;
		else if( (pi[t-1] == 2 && pi[h-1] == 5) || (pi[t-1] == 5 && pi[h-1] == 2) )
			coe = 1;
		else if( (pi[t-1] == 4 && pi[h-1] == 5) || (pi[t-1] == 5 && pi[h-1] == 4) )
			coe = 1;
		else
			coe = 0;
//	}
//	else
//	{
//		coe = 0.0;
//	}

	return coe;
}


/*int TRIPLE_PATH_CUT::genRow(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *var, ABA_ROW &row)
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

double TRIPLE_PATH_CUT::slack(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *variables, double *x)
{
	double val=0;
  	EdgeSet_It it;

  	for(it=it_begin;it != it_end;it++)
		val+=x[(*it -1)];

  	return(rhs_-val);
}*/

