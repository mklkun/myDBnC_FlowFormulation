#include "SP_PARTITION.h"
//#include "kHPPMASTER.h"
//#include "VAREDGE.h"

SP_PARTITION::SP_PARTITION(long *sp_part_list,long sp_part_p,int rhs_sp)
{
	pi = sp_part_list;
	p = sp_part_p;
}


SP_PARTITION::~SP_PARTITION()
{
	delete[] pi;
}


double SP_PARTITION::coeff(CUTINFO &c_info, int index) //(ABA_VARIABLE *v)
{
//	VAREDGE *ve = dynamic_cast<VAREDGE *>(v);
    lemonEdge e = c_info.G.edgeFromId(index);

	long t,h;
	double coe = 0.0;

//	if(ve == NULL)
//		coe = 0.0;
//	else
//	{
		t = c_info.G.idFromNode( c_info.G.u(e) ) + 1;
		h = c_info.G.idFromNode( c_info.G.v(e) ) + 1;

		if(pi[t-1] != pi[h-1])
		{
			if((pi[t-1] == p) || (pi[h-1] == p))
			{
				coe = 1.0;
			}
			else
			{
				if(abs(pi[t-1] - pi[h-1]) <= 1)
				{
					coe = 1.0;
				}
				else
				{
					coe = 2.0;
				}
			}
		}
		else
		{
			coe = 0.0;
		}
//	}

	return coe;
}


/*int SP_PARTITION::genRow(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *var, ABA_ROW &row)
{
  	ABA_CSENSE sense(master_,'G');
  	int i;
  	EdgeSet_It it;
	MapNumCoef_It c_it;

	if(p == 3)
	{
		ABA_ARRAY<int> s(master_,SP_List.size());
  		ABA_ARRAY<double> c(master_,SP_List.size());

  		i = 0;
  		it = sp_begin;
  		while(it != sp_end)
    		{
      			s[i]=(*it -1);
			c[i]=1.0;
      			i++;
      			it++;
    		}

  		ABA_ROW row2(((ABA_GLOBAL*)master_),((int) SP_List.size()),s,c,sense,rhs_);

  		row.copy(row2);

		return((int)SP_List.size());
	}
	else
	{
		ABA_ARRAY<int> s(master_,SP_List_Couple.size());
  		ABA_ARRAY<double> c(master_,SP_List_Couple.size());

		i = 0;
		c_it = sp_couple_begin;
		while(c_it != sp_couple_end)
		{
			s[i] = (*c_it).first - 1;
			c[i] = (*c_it).second;

//			cout << " N1 == "
//<<((kECMASTER*)master_)->get_graph_ptr()->Edges[(*c_it).first-1].back->adjac->id << " N2 == " <<
//((kECMASTER*)master_)->get_graph_ptr()->Edges[(*c_it).first-1].adjac->id << " C == " <<
//(*c_it).first << " " << (*c_it).second << endl;

      			i++;
      			c_it++;
		}


		ABA_ROW row2(((ABA_GLOBAL*)master_),((int) SP_List_Couple.size()),s,c,sense,rhs_);

  		row.copy(row2);

		return((int)SP_List_Couple.size());
	}

	return 0;
}

double SP_PARTITION::slack(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *variables, double *x)
{
	double val=0;

	if(p == 3)
  	{
		EdgeSet_It it;

  		for(it=sp_begin;it != sp_end;it++)
			val+=x[(*it -1)];
	}
	else
	{
		MapNumCoef_It it;

		for(it=sp_couple_begin;it != sp_couple_end;it++)
			val = val + x[(*it).first-1]*(*it).second;
	}

	return(rhs_-val);
}*/

