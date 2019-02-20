#ifndef TRIPLE_PATH_CUT_H
#define TRIPLE_PATH_CUT_H

/*extern "C"
{
	#include "../Prog_kEC/graph.h"
	#include "../Prog_kEC/common.h"
}*/

#include "I_graphe.h"
#include "I_common.h"

#include "extra_cuts.h"

#include <list>
#include <vector>


class TRIPLE_PATH_CUT
{
	private:
		int *pi;
		vector<int> F;
	public:
		TRIPLE_PATH_CUT(int *part,vector<int> &f_s,int rhs_p);
		virtual ~TRIPLE_PATH_CUT();
		virtual double coeff(CUTINFO &c_info, int index);


		/*virtual int genRow(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *var,ABA_ROW &row);
  		virtual double slack(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *variables, double *x);*/
};

#endif
