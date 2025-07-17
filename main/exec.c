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

//exec.c

#include "../Library/functions.h"
#include "../Library/odes.h"
#include "../Library/vacua.h"



struct param_const {
  double L0;
  double v;
  int nmin;
  int nmax; 
  int N; //number of modes
};

int
main (void)
{
  //Here we initiallyze the values of some parameters
  
  //If variables initiated externally, uncomment and provide the corresponding data.
  char outfile1 [80];
  char outfile2 [80];
  char outfile3 [80];
  int nmin, nmax;
  int mode;//Initial state peaked on the mode = ...;
  scanf ("%79s %79s %79s %d %d %d %lf %lf", outfile1, outfile2, outfile3, &nmin, &nmax, &mode);
  double t0, tf;//Initial and final times for the simulation to be fixed below.
  //DC trajectory
  //t0 = -3.0;
  //tf = 3.0;
  //FD trajectory
  double dl = 0.25;
  t0 = -0.3+(0.125/dl)*(0.125/dl)/10.0;
  tf = 0.4 + dl;
  struct param_const my_params = {1.0, 1.0, nmin, nmax, nmax-nmin};

  //Here we initialize the odes
  int n = my_params.N;
  gsl_odeiv2_system sys = {func, jac, 4*n, &my_params};
  //These are the absolute and relative tolerances

  double hinit = 1.0e-10;  //For evolution backward use negative values!!!
  double abs_e = 1.0e-12;  //absolute error for the adaptative method
  double rel_e = 0.0;      //relative error for the adaptative method

  //Integrators of the EOMs are

  gsl_odeiv2_driver * d = 
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
				   hinit, abs_e, rel_e);   
  
  double y[4*n];
  id_Mink(y, &my_params, mode, 0.0); //Initial state
  
  FILE *fp1;
  fp1 = fopen(outfile1, "w+");
  FILE *fp2;
  fp2 = fopen(outfile2, "w+");
  FILE *fp3;
  fp3 = fopen(outfile3, "w+");
    
  //Header modes
  
  fprintf (fp1,"# nmin = %d nmax = %d \n", nmin, nmax);
  fprintf (fp1,"#   t        realpha           imalpha\n");
  fprintf (fp1,"#\n");
  fprintf (fp1,"#\n");
  
  fprintf (fp2,"# nmin = %d nmax = %d \n", nmin, nmax);
  fprintf (fp2,"#   t        rebeta             imbeta\n");
  fprintf (fp2,"#\n");
  fprintf (fp2,"#\n");
  
  fprintf (fp3,"# nmin = %d nmax = %d \n", nmin, nmax);
  fprintf (fp3,"#   t        norm y           norm y0\n");
  fprintf (fp3,"#\n");
  fprintf (fp3,"#\n");
    
  double step = 2.0e-4; //For evolution backward use negative values!!!
  double t = t0;
  double ti = t + step;
  double Norm;//, Norm0;
  double alpha[2], beta[2];
  double y0[4*n];//Initiating out state
  int n1, m;//Out states peaked on the mode n1
  int tn = 0;
  for (t = t0; t < tf; t=t)
    {
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
        {
          printf("error, return value=%d\n", status);
          break;
        }
      if (tn <= ((int)((t-t0)/(2.0e-2))))
	{
	  fprintf (fp1,"%.9e ", t);
	  fprintf (fp2,"%.9e ", t);
	  fprintf (fp3,"%.9e ", t);
	  Norm = norms(y, &my_params);//Norm of the solution
	  fprintf (fp3,"%.9e ", Norm);
	  for(m = 0; m < n; m++)
	    {
	      n1 = m + nmin;
	      id_Mink(y0, &my_params, n1, 0.0); //Out state 
	      alpha[0] = 0.0; alpha[1] = 0.0; beta[0] = 0.0; beta[1] = 0.0;
	      Alpha(y, y0, alpha, &my_params);//alpha coeficients between in and out states
	      Beta(y, y0, beta, &my_params);//beta coeficients between in and out states
	      fprintf (fp1,"%.9e %.9e ", alpha[0], alpha[1]);
	      fprintf (fp2,"%.9e %.9e ", beta[0], beta[1]);
	    }
      
	  fprintf (fp1,"\n");
	  fprintf (fp2,"\n");
	  fprintf (fp3,"\n");
	  tn += 1;
	}
      ti = ti + step;
      fflush(stdout);
    }
    printf ("evolution computed\n");
    fflush(stdout);
   
   

  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  gsl_odeiv2_driver_free (d);
  
  
  
  return 0;
}
