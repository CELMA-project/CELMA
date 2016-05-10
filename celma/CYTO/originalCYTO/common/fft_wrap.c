#undef DEBUG
#undef COMM_STATUS
#define COMM2(a)

#include "utilities.h"

#ifdef INTELFFT
#include <mkl_fft.h>
#endif

#ifdef FFTW
#include <fftw3.h>
#endif

#ifdef FFTW
/* 

   inplace 1D forward fft of  2D array 
 
   input  : f     2d array of dimension [0:ny][0:nx] or [-1:ny][-1;nx], slow index FFT'ed, note that two extra elements after ny are needed. 
   ny    dimension of the fft
   nx    number of fft operation

   output : f    
*/
void  Fft_1d_2d_f(double **f,int nx,int ny)
{
  static int i, j;
  static int FIRST=TRUE;
  static double *in;
  static fftw_complex *out;
  static double scale;
  static fftw_plan plan;
  static int nc = 0;

  BUGREPORT;
  
  if(FIRST)
    {
      nc    = (ny/2) +1;
      in    = (double *) fftw_malloc(sizeof(double) * ny);
      out   = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nc);
      plan  = fftw_plan_dft_r2c_1d(ny, in,out,  FFTW_ESTIMATE);
      scale = 1./(double)ny;
      FIRST = FALSE;
    } 

  BUGREPORT;
  for(j=0;j<nx;j++)
    {
        BUGREPORT;
        for(i=0;i<ny;i++) in[i] = scale*f[i][j];
        BUGREPORT;
        fftw_execute(plan);
        BUGREPORT;

        for(i=0;i<nc;i++) {
            f[2*i]  [j]   = out[i][0];
            f[2*i+1][j]   = out[i][1];
        }
    } 
}

void  Fft_1d_2d_b(double **f,int nx,int ny)
{
  static int i, j;
  static fftw_complex  *in;
  static double *out;
  static int FIRST=TRUE;
  static fftw_plan plan;
  static int nc = 0;


  BUGREPORT;
  
  if(FIRST)
    {
      nc    = (ny/2) + 1;
      out   = (double *) fftw_malloc(sizeof(double) * ny);
      in    = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nc);
      plan  =  fftw_plan_dft_c2r_1d(ny,in, out,FFTW_ESTIMATE);
      FIRST = FALSE;
    } 
 

  BUGREPORT;
  for(j=0;j<nx;j++)
    {
        for(i=0;i<nc;i++) 
        {
            in[i][0] = f[2*i][j];
            in[i][1] = f[2*i+1][j];
        }

       // in[0][1]=0.;/* Set phase to zero */
       // in[ny][1] = 0.;/* Set phase to zero for highest mode as well */
	   // in[ny+1][1] = 0.;/* Set phase to zero for highest mode as well */

        fftw_execute(plan);
        for(i=0;i<ny;i++) f[i][j] = out[i]; 
    }
}






#elif INTELFFT
/* 1D forward fft of an ordinary 2D array */ 
/* 
   input  : f     2d array of dimension [-1:ny][-1;nx], slow index FFT'ed 
            ny    dimension of the fft
	    nx    number of fft operation

   output : f     2d array of dimension [-1:ny][-1;nx]
*/
void  Fft_1d_2d_f(double **f,int lnx,int lny)
{
  static int i, j;
  static int FIRST=TRUE;
  static double *res;
  static double scale;


  if(FIRST)
    {
        res = Util_DVector(lny,2);
      FIRST = FALSE;
    } 
  
  for(j=0;j<lnx;j++)
    {
      for(i=0;i<lny;i++) res[i] = f[i][j];
      Fft_1d_f(res,lny);
      for(i=0;i<=lny;i++) f[i][j] = res[i];
    } 
}

void  Fft_1d_2d_b(double **f,int lnx,int lny)
{
  static int i, j;
  static double *res;
  static int FIRST=TRUE;


  if(FIRST)
    {
      res = Util_DVector(lny,2);
      FIRST = FALSE;
    } 
 

  for(j=0;j<lnx;j++)
    {
      for(i=0;i<=lny;i++) res[i] = f[i][j];
      Fft_1d_b(res,lny);
      for(i=0;i<lny;i++) f[i][j] = res[i]; 
    }
}

#endif


#ifdef INTELFFT

/***********************************************************************************************/
/* 1D backward (inverse) fft */ 
void  Fft_1d_b(double *f,int lny)
{
  static int ny_old =0;
  static double *wsave=NULL;
  int isign;


  if(ny_old != lny)
    {
      if(wsave != NULL) free(wsave);
      wsave =(double*) malloc((size_t)((2*lny+10)*sizeof(double)));
      isign = 0;
      zdfft1d(&f[0],&lny,&isign,&wsave[0]);
      ny_old = lny;
    }

  isign = 1;
  zdfft1d(&f[0],&lny,&isign,&wsave[0]);

}

/* 1D forward fft  */ 
void  Fft_1d_f(double *f,int lny)
{
  static int ny_old =0;
  static double *wsave=NULL;
  int isign;

  if(ny_old != lny)
    {
      if(wsave != NULL) free(wsave);
      wsave =(double*) malloc((size_t)((2*lny+10)*sizeof(double)));
      isign = 0;
      dzfft1d(&f[0],&lny,&isign,&wsave[0]);

      ny_old = lny; 
    }

  isign = 1;     
  dzfft1d(&f[0],&lny,&isign,&wsave[0]);

}
 

#include "fft_1D_transpose.c" 
#endif
 
