#include<stdio.h>
#include<math.h>
#include<nlopt.h>


float minimize(){
	double lb[2] = { -HUGE_VAL, 0 };
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LN_NELDERMEAD, 1); /* algorithm and dimensionality */
	nlopt_set_min_objective(opt, myfunc, NULL);

}
