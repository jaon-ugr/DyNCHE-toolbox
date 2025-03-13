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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

/****************  Declaration of functions  *************************/

struct param_const;

void Ltft_DC (const double t, void *params, double ft[]);

void Ltft_FD (const double t, void *params, double ft[]);

double norms (const double y[], void *params);

void Alpha (const double y1[], const double y2[], double alpha[], void *params);

void Beta (const double y1[], const double y2[], double beta[], void *params);

