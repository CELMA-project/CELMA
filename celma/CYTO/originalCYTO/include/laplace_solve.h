#ifndef LAPLACE_H_
#define _LAPLACE_H_

/*! 
  Public Functions for solving Laplace Equation: Headers. 
*/
 
void Laplace_Solve3D(HDF_DS *data, PARA *p,double ***f, double ***w,
	       double **bdra,double **bdrb,double *hval,double *rcor,int mbdcnd,
	       double lamda,int scale);

void Laplace_Solve2D(HDF_DS *data, PARA *p,double **f, double **w,
               double *bda,double *bdb,double *hval,double *rcor,int mbdcnd,
               double lamda, int scale);

void Laplace_Solve2D_Par(HDF_DS *data, PARA *p,double **f, double **w,
                         double *ibdra,double *ibdrb,double *hval,double *rcor,
                         int mbdcnd, double lamda,int scale);

/* Local functions */
void LapSol_PeriodicSolver(HDF_DS *data, PARA *p,double **f, double **w,double lamda, int normalize);
int  LapSol_GaussElim(int n, double *dl, double *d, double *du, double *b);
int  LapSol_GaussElim_Complex(int n, double *dl, double *d, double *du,double *b);

#endif /* !_LAPLACE_H */ 
