/*
These codes have been developed by J. Olmedo (so far). You are free to use them for any
scientific or pedagogical purpose. However, the author is not responsible of any 
mistake or typo on these codes. Nevertheless, he will acknowledge any feedback
at this respect. Besides, if you plan to use these codes (or a modified version 
based on them) for scientific purposes (i.e. not pedagogical purposes), please,
he will very much appreciate if you cite any of the following references:
Date: 12/03/25.
Place: Universidad de Granada
*/


#include "functions.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

struct param_const;

/****************  Declaration of odes  *************************/

int func (double t, const double y[], double f[], void *params);

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

