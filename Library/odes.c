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


#include "odes.h"




/****************  Definition of the odes  *************************/

struct param_const {
/**
 *  Some global parameters
 */

  double L0;
  double v;
  int nmin;
  int nmax; 
  int N; //number of modes
  
};

int
func (double t, const double y[], double f[],
      void *params)
{
/**
 */(void)(t); /* avoid unused parameter warning */
  
  struct param_const *my_params_pointer = params;
  double l0 = my_params_pointer->L0;
   
  //variables to store 1/L(t), 1/L * (dL/dt), 1/L * (df/dt). This improves efficiency ... ?
  double Ft[3]; 

  //Generation of the Ft[0] = 1/L(t), Ft[1] = 1/L * (dL/dt), Ft[2] = 1/L * (df/dt)
  //Ltft_DC(t, &*my_params_pointer, Ft);
  Ltft_FD(t, &*my_params_pointer, Ft);
  

  int N = my_params_pointer->N;
  int nmin = my_params_pointer->nmin;
  
  //These are some integers for the loops. 
  //n and m refer to modes running from 0 to N-1, 
  //li is an integer that accounts for 2 configs and 2 velocities. In total we have 4*N EOMS. 
  int n, m, li = 0, lj = 0;
  int n1, m1;
  
  //These are some doubles for the loops 
  double k2, coupln_phin_re, coupln_phin_im, coupln_pin_re, coupln_pin_im, aux1; 
  
  for(n = 0; n < N; n++){
    
    n1 = n + nmin;
  
    k2 = (((1.0*n1) * (1.0*n1))) * (M_PI * M_PI); //This is (n * pi)^2 
    
    //Here we compute the off diagonal couplings
    coupln_phin_re = 0.0;
    coupln_phin_im = 0.0;
    coupln_pin_re = 0.0;
    coupln_pin_im = 0.0;
    lj = 0; 
    for(m = 0; m < N; m++)
    {
      m1 = m + nmin;
      if(m1!=n1){
        aux1 = (Ft[2] * (-1.0 + pow(-1.0,n1 + m1)) + Ft[1] * pow(-1.0,n1 + m1)) * ((1.0*n1) * (1.0*m1)) / ((1.0*m1) * (1.0*m1) - (1.0*n1) * (1.0*n1));
        coupln_phin_re += aux1 * y[m+lj];
        coupln_phin_im += aux1 * y[m+2+lj];
        coupln_pin_re += aux1 * y[m+1+lj];
        coupln_pin_im += aux1 * y[m+3+lj];
      }
      lj = lj + 3;
    }
    
    //printf ("%.5e %.5e \n", t, y[4]*y[4]+y[4+2]*y[4+2]);
    
    f[n+li] = y[n+1+li] * Ft[0] - y[n+li] * Ft[1] / 2.0 + 2.0 * coupln_phin_re;
    f[n+1+li] = - y[n+li] * k2 * Ft[0] + y[n+1+li] * Ft[1] / 2.0 + 2.0 * coupln_pin_re;
    f[n+2+li] = y[n+3+li] * Ft[0] - y[n+2+li] * Ft[1] / 2.0 + 2.0 * coupln_phin_im;
    f[n+3+li] = - y[n+2+li] * k2 * Ft[0] + y[n+3+li] * Ft[1] / 2.0 + 2.0 * coupln_pin_im;
    
    li = li + 3;
    }
  
  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy, 
     double dfdt[], void *params)
{
/**
 *  Definition of the jacobian of odes of perturbations in Bianchi I spacetimes in GR as introduced in arXiv:0707.0736 and arxiv:0801.3596 necessary for implicit methods (in some cases we set it to zero since we focus on explicit methods).
 */
  (void)(t); /* avoid unused parameter warning */
  struct param_const *my_params_pointer = params;
  //double l = my_params_pointer->lambda;
  //double g = my_params_pointer->gamma;
  int n = my_params_pointer->N;
  int nn = 4*n;
  gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy, nn, nn);
  gsl_matrix * m = &dfdy_mat.matrix; 
  
  int i, j;
  for(i = 0; i < nn; i++){
    dfdt[i] = 0.0;
    for(j = 0; j < nn; j++){
       gsl_matrix_set (m, i, j, 0.0);
    }
  }
  
  return GSL_SUCCESS;
}


