#ifndef PARTITION_H
#define PARTITION_H

/*extern "C"
{
	#include "../Prog_kEC/graph.h"
	#include "../Prog_kEC/common.h"
}*/

#include "I_graphe.h"
#include "I_common.h"

#include <list>

class PARTITION
{
	private:
		long *pi;
		long p_sz;
	public:
		PARTITION(long *part,long sz,int rhs_p);
		virtual ~PARTITION();
		virtual double coeff(CUTINFO &c_info, int index);

		/*virtual int genRow(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *var,ABA_ROW &row);
  		virtual double slack(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *variables, double *x);*/
};

#endif
