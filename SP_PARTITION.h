#ifndef SP_PARTITION_H
#define SP_PARTITION_H

/*extern "C"
{
	#include "../Prog_kEC/graph.h"
	#include "../Prog_kEC/common.h"
}*/

#include "I_graphe.h"
#include "I_common.h"

#include "extra_cuts.h"

#include <list>
#include <map>


class SP_PARTITION
{
	private:
		long *pi;
		long p;
	public:
		SP_PARTITION(long *sp_part_list,long sp_part_p,int rhs_sp);
		virtual ~SP_PARTITION();
		virtual double coeff(CUTINFO &c_info, int index);

		//virtual int genRow(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *var,ABA_ROW &row);
  		//virtual double slack(ABA_ACTIVE<ABA_VARIABLE,ABA_CONSTRAINT> *variables, double *x);
};

#endif
