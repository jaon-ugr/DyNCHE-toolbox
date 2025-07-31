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

/****************  Some functions  *************************/

struct param_const {
/**
 *  Some global parameters
 */

  double L0;
  double v;
  int nmin;
  int nmax; 
  int N; //number of modes
  //additional parameters for some trajectories
  double t; 
  
};

void
Ltft_DC (const double t, void *params, double ft[])
{
/**
 * Here we compute the functions 1/L(t), 1/L * (dL/dt), 1/L * (df/dt).
 */
//   (void)(t);
  struct param_const *my_params_pointer = params;
  double L0 = my_params_pointer->L0;
  double phi = 1.0*M_PI, epsilon = 0.25;//for epsilon = 0.1 we get (dot f) and (dot g) order 1
  double exp2, sinqpi, sinphiqpi, sinphi;
  double qpi = 0.5 * M_PI;
  double tbump = 16.0;
  double s = 10.0; //steep of the bump function. 

  if(t < 0.01 * tbump)
  {
    ft[0] = 1.0 / L0;
    ft[1] = 0.0;
    ft[2] = 0.0;
  }
  else
  {
    if(t > 1.99 * tbump)
    {
      ft[0] = 1.0 / L0;
      ft[1] = 0.0;
      ft[2] = 0.0;
    }
    else
    {
      exp2 = exp(tbump * tbump / s / (2.0 * tbump  - t) / t) / exp(1.0 / s);
      sinqpi = sin((qpi*t)/L0);
      sinphiqpi = sin(phi + (qpi * t) / L0);
      sinphi = sin(phi);
      //if f(t) != 0
      //ft[0] = 1.0 / (L0 - 1.0 / exp2 * epsilon * (sinphi + sinqpi - sinphiqpi));
      //ft[1] = (4.0 * epsilon * sin(phi/2.0) * (2.0*L0*tbump*tbump * (tbump - t)*sin(qpi/2.0*t/L0)*sin(phi/2.0 + qpi / 2.0 *t/L0) + qpi / 2.0 * s * t*t*(-2.0*tbump + t)*(-2.0*tbump + t)*sin(phi/2.0 + (qpi*t)/L0))) / (L0*s*t*t*(-2.0*tbump + t)*(-2.0*tbump + t)*(-exp2*L0 + epsilon*(sinphi + sinqpi - sinphiqpi)));
      //ft[2] = (2.0 * epsilon*(-qpi / 2.0 *s*t*t*(-2.0*tbump + t)*(-2.0*tbump + t)*cos((qpi*t)/L0) + L0*tbump*tbump*(-tbump + t)*sinqpi)) / (L0*t*t*s*(-2.0*tbump + t)*(-2.0*tbump + t)*(-exp2*L0 + epsilon*(sinphi + sinqpi - sinphiqpi))); 
      //if f(t) = 0
      ft[0] = 1.0 / (L0 - 1.0 / exp2 * epsilon * (sinphi - sinphiqpi));
      ft[1] = (epsilon * (2.0*L0*tbump*tbump* (tbump - t) * (-sinphi + sinphiqpi) + qpi * s * t*t*(-2.0*tbump + t)*(-2.0*tbump + t)*cos(phi + (qpi*t)/L0))) / (L0*s*t*t*(-2.0*tbump + t)*(-2.0*tbump + t)*(exp2*L0 - epsilon*(sinphi - sinphiqpi)));
      ft[2] = 0.0; 
    }
  }
}



void
Ltft_FD (const double t, void *params, double ft[])
{
/**
 * Here we compute the functions 1/L(t), 1/L * (dL/dt), 1/L * (df/dt).
 */
//   (void)(t);
  struct param_const *my_params_pointer = params;
  double L0 = my_params_pointer->L0;
  double s = 12.5;
  double t1 = 1.0/8.0;
  double dl = 0.375;
  double k = s / dl;
  double kt = k * (t - t1);

  
  //function ft
  
  //1 or 2 plates moving right
  double L, dL;
  double plt = 1.0;
  L = L0 + dl / 2.0 / plt * (1.0 + (-log(cosh(s - kt)) + log(cosh(kt))) / s);
  dL = (dl * (k * tanh(s - kt) + k * tanh(kt)))/(2.0 * plt * s);
  
  ft[0] = 1.0 / L;
  ft[1] = dL * ft[0];
  ft[2] = (plt - 1.0) * dL * ft[0];
  

  //2 plates moving oposite directions
  //double L, dL;
  //L = L0 + dl * (1.0 + (-log(cosh(s - kt)) + log(cosh(kt))) / s);
  //dL = (dl * (k * tanh(s - kt) + k * tanh(kt)))/(2.0 * s);
  
  //ft[0] = 1.0 / L;
  //ft[1] = 2.0 * dL * ft[0];
  //ft[2] = - dL * ft[0];


  //2 plates moving rigidly
  /* double df = (dl * (k * tanh(s - kt) + k * tanh(kt)))/(2.0 * s); */
  
  /* ft[0] = 1.0 / L0; */
  /* ft[1] = 0.0; */
  /* ft[2] = df; */

}




double
PW_f (const double t, void *params)
{
double f;

//piecewise function f
  if(t <= 0.0)
    f = 0.0;
  else
    f = exp(-1.0/t);

return f;
}

//Norms of a solution

double
norms (const double y[], void *params)
{
//   (void)(t);
  struct param_const *my_params_pointer = params;
  int N = my_params_pointer->N;
  int n, li = 0;
  double norm = 0.0;
  for(n = 0; n < N; n++){
    norm += (y[1+n+li] * y[2+n+li] - y[n+li] * y[3+n+li]);
    li += 3;
  }
  return norm;
}



//bogo coefficients of a solution

void
Alpha (const double y1[], const double y2[], double alpha[], void *params)
{
//   (void)(t);
  struct param_const *my_params_pointer = params;
  int N = my_params_pointer->N;
  int n, li = 0;
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  for(n = 0; n < N; n++){
    alpha[0] += (-y1[3 + n + li] * y2[n + li] + y1[2 + n + li] * y2[1 + n + li] + y1[1 + n + li] * y2[2 + n + li] - y1[n + li] * y2[3 + n + li]);
    alpha[1] += (-y1[1 + n + li] * y2[n + li] + y1[n + li] * y2[1 + n + li] - y1[3 + n + li] * y2[2 + n + li] + y1[2 + n + li] * y2[3 + n + li]);
    li += 3;
  }
  alpha[0] = 0.5 * alpha[0];
  alpha[1] = 0.5 * alpha[1];
}


void
Beta (const double y1[], const double y2[], double beta[], void *params)
{
//   (void)(t);
  struct param_const *my_params_pointer = params;
  int N = my_params_pointer->N;
  int n, li = 0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  for(n = 0; n < N; n++){
    beta[0] += (-y1[3 + n + li] * y2[n + li] + y1[2 + n + li] * y2[1 + n + li] - y1[1 + n + li] * y2[2 + n + li] + y1[n + li] * y2[3 + n + li]);
    beta[1] += (y1[1 + n + li] * y2[n + li] - y1[n + li] * y2[1 + n + li] - y1[3 + n + li] * y2[2 + n + li] +  y1[2 + n + li] * y2[3 + n + li]);
    li += 3;
  }
  beta[0] = 0.5 * beta[0];
  beta[1] = 0.5 * beta[1];
}


