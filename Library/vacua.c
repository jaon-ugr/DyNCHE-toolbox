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

#include "vacua.h"

/****************  Initial conditions for the field modes  *************************/

struct param_const {
  double L0;
  double v;
  int nmin;
  int nmax; 
  int N; //number of modes
};

void
id_Mink (double y[], void *params, int m1, double t)
{
//   (void)(t);
  //This initial data has been provided externally.
  struct param_const *my_params_pointer = params;

  double L = my_params_pointer->L0;
  int N = my_params_pointer->N;
  int nmin = my_params_pointer->nmin;
  int nmax = my_params_pointer->nmax;
  //ID for the inhomogeneities dof
  //Here y[9+4*i] = reu; y[9+4*i+1] -> reu'; y[9+4*i+2] -> imu; y[9+4*i+3] -> imu'; 
  double w; 
  int n, li = 0; 
  int n1;
  for(n = 0; n < N; n++)
  {
  n1 = n + nmin;
    if(n1==m1)
    {
      w = sqrt(((double)((n1) * (n1))) * (M_PI * M_PI));
      //y[n+li] = 1.0 / sqrt(w) * cos (w * t);
      //y[n+1+li] = 0.0; 
      //y[n+2+li] = 0.0;
      //y[n+3+li] = - sqrt(w);
      y[n+li] = 1.0 / sqrt(w) * cos (w * t);
      y[n+1+li] = sqrt(w) * sin (w * t); 
      y[n+2+li] = 1.0 / sqrt(w) * sin (w * t);
      y[n+3+li] = - sqrt(w) * cos (w * t);
    }
    else
    {
      y[n+li] = 0.0;
      y[n+1+li] = 0.0; 
      y[n+2+li] = 0.0;
      y[n+3+li] = 0.0;
    }
      li = li + 3;    
  }
  
}

