/*! \file utilities.h */


#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923
#endif
#ifndef M_1_PI
#define M_1_PI      0.31830988618379067154
#endif
#ifndef M_PI_4
#define M_PI_4      0.78539816339744830962
#endif
#ifndef M_2_PI
#define M_2_PI      0.63661977236758134308
#endif
#ifndef M_2_SQRTPI
#define M_2_SQRTPI  1.12837916709551257390
#endif
#ifndef M_SQRT2
#define M_SQRT2     1.41421356237309504880
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2   0.70710678118654752440
#endif



#ifndef _UTILITIES_H_
#define _UTILITIES_H_



#define FORALL for(iz=0;iz<nz;iz++) for(ip=0;ip<ny;ip++) for(ir=0;ir<nx;ir++)
#define FORALL_BD for(iz=-d->offz;iz<(nz+d->offz);iz++) for(ip=-d->offy;ip<(ny+d->offy);ip++) for(ir=-d->offx;ir<(nx+d->offx);ir++)

#define FORYX    for(ip=0;ip<ny;ip++)    for(ir=0;ir<nx;ir++)
#define FORYX_BD for(ip=-1;ip<ny+1;ip++) for(ir=-1;ir<nx+1;ir++) 

#define FORZY for(iz=0;iz<nz;iz++) for(ip=0;ip<ny;ip++) 
#define FORZY_BD for(iz=-d->offz;iz<nz+d->offz;iz++)  for(ip=-d->offy;ip<(ny+d->offy);ip++)

#define FORZX  for(iz=0;iz<nz;iz++)  for(ir=0;ir<nx;ir++)
#define FORZX_BD  for(iz=-d->offz;iz<(nz+d->offz);iz++)  for(ir=-d->offx;ir<(nx+d->offx);ir++)

#define FORX    for(ir=0;ir<nx;ir++) 
#define FORX_BD for(ir=-d->offx;ir<(nx+d->offx);ir++) 

#define FORY    for(ip=0;ip<ny;ip++) 
#define FORY_BD for(ip=-d->offy;ip<(ny+d->offy);ip++) 

#define FORZ    for(iz=0;iz<nz;iz++) 
#define FORZ_BD for(iz=-d->offz;iz<(nz+d->offz);iz++) 

#define FORALL_BD_XY for(iz= 0;iz<nz;iz++) for(ip=-1;ip<(ny+1);ip++) for(ir=-1;ir<(nx+1);ir++)


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stddef.h>
#include <stdarg.h>
#include <ctype.h>
#include <float.h>
#include <signal.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h> 
#include <mpi.h>
#include <hdf.h>
#include <mfhdf.h>
#include <df.h>
#include <dfsd.h>
#include <stdint.h>

#include "definitions.h"
#include "fft_local.h"
#include "metric.h"
#include "util_alloc.h"
#include "file_utils.h"
#include "hdf_helper.h"
#include "ode_pack.h"
#include "laplace_solve.h"
#include "fft_local.h"
#include "metric.h"
#include "util_inifiles.h"

#include "diff.mac"

#define UT_COLORS 256

#ifdef CHECK
#define RGBFAK 256
/*65535*/
#define CMOFF 
#else
#define RGBFAK 256
#define CMOFF
#endif

#define PERIODIC 0
#define DIRDIR 1
#define DIRNEU 2
#define NEUNEU 3
#define NEUDIR 4
#define ZeroDIR 5
#define ZeroNEU 6


#define IS_ELECTRONVELOCITY 1
#define IS_IONVELOCITY 2
#define IS_DENSITY 3
#define IS_POTENTIAL 4
#define IS_TEMPERATURE 5
#define IS_VORTICITY 6
#define IS_CURRENT 7
#define IS_TE 8
#define IS_HEATFLUX_TE 9


#ifndef DATABASEPATH
#define DATABASEPATH "~/database"
#endif


#define SETRGB(x,y,z) { *r = (unsigned int) ((x) * (RGBFAK-1) CMOFF);\
			*g = (unsigned int) ((y) * (RGBFAK-1) CMOFF); \
			*b = (unsigned int) ((z) * (RGBFAK-1) CMOFF);}



#ifndef linux
#pragma ALLOCS_NEW_MEMORY  Util_DMatrix, Util_FMatrix, Util_DCube, Util_FCube, Util_DVector 
#endif 


#ifdef __cplusplus
extern "C"
{
#endif


    
/* Old ones */

    void Util_2DFullBdd_float(float **f,int nx,int ny,float *val1,float *val2);
    
    void Util_2DFullBdn_float(float **f,int nx,int ny,float *val1,float *val2,float dy);
    
    void Util_2DFullBdpf(float **f,int nx,int ny);
    
    
    
    float Util_DistF(float **v,float **b,int nx,int ny);
    
    double Util_Dist(double **v,double **b,int nx,int ny);
    
    
    
    
    double Util_Transport(HDF_DS *data,PARA *para, double **p, double **tmp);
    
    void Util_Noise(HDF_DS *data,double **f);
    
    void Util_CalcP(HDF_DS *data, PARA *para,double **fp1,double **f0, double **fm1, double **p,double *Sx,double *Sy);
    
    
/*****************************************/
/* Statistics                            */
    
    
    void Util_Moment(double **data,int nx, int ny,double *ave,double *adev,double *sdev,double *svar,double *skew,double *curt);
    
    double Util_Correlation(double **x,double ax,double **y,double ay,int nx,int ny);
    

    
/*****************************************/
/* Integration                           */
    
    
    double Util_1DSimpson(double *u,int n,double dx);
    float Util_SimpsonF(float *u,int n,float dx);
    double Util_2DSimpson(double **u,int nx, int ny, double dx, double dy);
    void Util_Integral(double **f,int order,double *norm,int np,int nr,double *result);
    void Util_3DIntegral(double ***f,int order, double *norm,int nx,int ny, int nz,double **help,double *result);
    
    double Util_CalcCirc(double *wa,double *wb,double *norm,int ny,int nx);
    
/*****************************************/
/* Differentiation                       */
    
    void Util_RectLaplace(double **out, double **in,int nx,int ny,double dx,double dy);
    void Util_RectDx(double **out, double **in,int nx,int ny,double dx,double dy);
    void Util_RectDy(double **out, double **in,int nx,int ny,double dx,double dy);
    void Util_RectStrain(double **out, double **in, int nx, int ny,double dx, double dy);
    void Util_RectWeiss(double **out, double **in, int nx, int ny,double dx, double dy);
    
    void Util_CylWeiss(double **out, double **vr, double **vp, double **w,
                       double *vval, double *hval, double *edr, int nx, int ny);
    
/*****************************************/
/* Plotting                              */
    
    
    void Util_CreateColors(char *kindof,unsigned int *, unsigned int *, unsigned int *);
    void  Util_WritePPM(FILE *file, int width, int height,unsigned char **in,char *color);
    void Util_MakeUCharField(double **in,unsigned char **out, int nx, int ny, double max, double min,int extra);
    void Util_MapRect(double **in, double  **out,int nxin, int nyin, int nx,int ny, double rmin, double rmax, double phimin, double phimax,int transpose,double *rcor, double *phicor,int remap);
    void Util_MapRect2(double  **in, double **out,int nxin, int nyin, int nx,int ny, double rmin, double rmax, double phimin, double phimax,int transpose,double *rcor, double *phicor,int remap);
    void Util_MapIn(double  **in, double **out,int nyin, int nxin, int ny, int nx, double xmin, double xmax,
                    double ymin, double ymax,int transpose,double *ycor, double *rcor,int newmapping);
    int Util_HsvToRgb(float h,float s,float v,unsigned int *r,unsigned int *g,unsigned int *b);
    
    
/*****************************************/
/* MaxMIn                                */
    
    
    double Util_AbsMaxField(double **in,int nx,int ny,int *xpos, int *ypos);
    double Util_MaxField(double **in,int nx,int ny,int *xpos, int *ypos);
    double Util_MinField(double **in,int nx,int ny,int *xpos, int *ypos);
    
    double Util_3DMaxField(double ***in,int nz, int ny,int nx,int *xpos, int *ypos,int *zpos);
    double Util_3DMinField(double ***in,int nz, int ny,int nx,int *xpos, int *ypos, int *zpos);
    
/*****************************************/
/* Upwind                                */
    
    void Util_VanLeerX(double **w_0,double **wxl,double **wxr,int nx,int ny);
    void Util_VanLeerY(double **w_0,double **wyl,double **wyr,int nx,int ny);
    
    void Util_ExtraX(double **w_0,double **wxl,double **wxr,int nx,int ny);
    void Util_ExtraY(double **w_0,double **wyl,double **wyr,int nx,int ny);
    
    void Util_UpWind(double **vrp,double **vrxl,double **vrxr,int nx,int ny);
    void Util_UpwindW(double **wxp,double **wxr,double **wxl,double **vrp,int nx,int ny);
    
    
/*****************************************/
/* Boundaries                            */
    

    void Util_BdPer(double **f,double *vala, double *valb,double *delta,int nx,int ny,int shifted);
    void bound_neumann(double **f,double *vala, double *valb,double *delta, int nx,int ny,int shifted);
    void Util_BdNeuNs(double **f,double *vala, double *valb,double *delta, int nx,int ny,int shifted);
    
    void Util_CylBdDir(double **f,int np,int nr,double *bdra,double *bdrb,double *dx);
    void Util_CylBdNeu(double **f,int np,int nr,double *bdra,double *bdrb,double *dx);

    
    void Util_2DFullBd(double **u,int ny,int nx,double *val1,double *val2,double *dx, int type);
    
    
    void Util_3DFullBd(double ***u,HDF_DS *d,int nx,int ny, int nz,
                       double **left,double **right,double **up, double **down,
                       double *dx,double dz, int typer, int typez);
    
    
    void Util_BdNeuA(double **f,double *val,double *delta, int nx,int ny,int shifted);  
    void Util_BdNeuB(double **f,double *val,double *delta, int nx,int ny,int shifted);
    
    void Util_BdDirA(double **f,double *val,double *delta, int nx,int ny,int shifted);
    void Util_BdDirB(double **f,double *val,double *delta, int nx,int ny,int shifted);
    
    void Util_BdPerX(double **f,double *val,double *delta, int nx,int ny,int shifted);
    void Util_BdCylA(double **f,double *val,double *delta, int nr,int np,int shifted);
    
    void Util_3DBdDirA(double ***f,double **val,HDF_DS *d,double delta, int nx,int ny,int nz);
    void Util_3DBdDirB(double ***f,double **val,HDF_DS *d,double delta, int nx,int ny,int nz);
    void Util_3DBdNeuB(double ***f,double **val,HDF_DS *d,double delta, int nx,int ny,int nz);
    void Util_3DBdNeuA(double ***f,double **val,HDF_DS *d,double delta, int nx,int ny,int nz);
    void Util_BdPerZ(double ***f,double **val,double hval, int nx,int ny,int nz);
    
    
    




/*****************************************/
/* Timestepping                          */
void Util_1DSsTimeStep(int iter,int nx, double dt,double mue,double *lamda,
                     double *w,double *w_0, double *w_1, double *w_2,
                       double *r_0, double *r_1, double *r_2);
    

void  Util_Ss1(double **w, double **w_0,double **r_0,int  ny, int nx, double dt);
void  Util_Ss2(double **w, double **w_0,double **w_1,double **r_0,double **r_1,int  ny, int nx, double dt);
void  Util_Ss3(double **w, double **w_0,double **w_1,double **w_2,double **r_0,double **r_1,double **r_2,int ny,int nx, double dt);

void  Util_1DSs1(double *w, double *w_0,double *r_0, int nx, double dt);
void  Util_1SSs2(double *w, double *w_0,double *w_1,double *r_0,double *r_1,int nx, double dt);
void  Util_1DSs3(double *w, double *w_0,double *w_1,double *w_2,double *r_0,double *r_1,double *r_2,int nx, double dt);


void Util_2DSsTimeStep(int iter,int nx, int ny, double dt,double mue,double *lamda,
	       double **w,double **w_0, double **w_1, double **w_2,
	       double **r_0, double **r_1, double **r_2);


void  ddd_Util_Ss1(double ***w, double ***w_0,double ***r_0,int nz,int  ny, int nx, double dt);
void  ddd_Util_Ss2(double ***w, double ***w_0,double ***w_1,double ***r_0,double ***r_1,int nz,int  ny, int nx, double dt);
void  ddd_Util_Ss3(double ***w, double ***w_0,double ***w_1,double ***w_2,double ***r_0,double ***r_1,double ***r_2,int nz,int ny,int nx, double dt);

void Util_3DSsTimeStep(int iter,int nx, int ny, int nz, double dt,double mue,double *lamda,
		     double ***w,double ***w_0, double ***w_1, double ***w_2,
		     double ***r_0, double ***r_1, double ***r_2);

/*****************************************/
/* Poisson Bracket                       */

void Util_ArakawaNl(double **res,double **a,double **b,double *fac,int nx,int ny);
void Util_Arakawa45(double **res,double **a,double **b,double *fac,int nx,int ny);
void Util_Arakawa46(double **res,double **a,double **b,double *fac,int nx,int ny);

void Util_3DArakawaNl(double ***res,double ***a,double ***b,double **fac,int nz,int nx,int ny);
void Util_3DSet2Zero(double ***res,double ***a,double ***b,double **fac,int nz,int nx,int ny);
    
/*****************************************/
/* Fourier                               */


void Util_ZeroPadPer(double **in,int cutx, int cuty, int nx, int ny);
void Util_ZeroPadDir(double **in,int cutx, int cuty, int nx, int ny);


/*****************************************/
/* Geometry                              */

int  Util_SetupGeom(double **rcor,double **edr, double **hval,double **vval,double **gval,
                    double **rs, double **edrs, double **hs, double **vs, double **gs,
                    double **zcor,double **zhval,double **zgval,
                    double **norm, HDF_DS *data,PARA *p,int rindex);


void Util_ForceNoSlip(HDF_DS *data,PARA *para,double **w,double **f,
				    double *hv,double *gv,double *edr,double *dyya,double *dyyb);




void Util_Glasser2Cart(double *x, double *y, double *z,double phi,double psi,double r,double q, double R0,double r0,double fr);



void Util_Glasser2Cyl(double *cylr, double *cylphi, double *cylz,double phi,double psi,double r,double q, double R0,double r0,double fr);


void Util_CalcNoSlip(double *bdra,double *bdrb,double **phi,
			     int mbdcnd, double *hval, double *hs, int nx,int ny);


int  Util_2DCheckCFL(double **f, double **fvr, double **fvp, HDF_DS *data,
	       PARA *para, double *hval, double *vval,double *cfl, double *cfr);


int  Util_3DCheckCFL(double ***f, double ***fvr, double ***fvp, 
		  HDF_DS *data,PARA *para, double *hval, double *vval,
		  double *cfr, double *cfp);


int  Util_CheckCFLMetric(double **f, double **fvr, double **fvp, HDF_DS *data,PARA *para, 
		      double *hval, double *vval,double *kai, double *grr, double *gff, double *grf,
		      double *norm_gr, double *norm_gp,double *cfr, double *cfp);

void Util_CalcParameters(PARA *p,HDF_DS *d);
    

/*********************************************/
/*MPI datatypes*/

MPI_Datatype hdfdata_type;
MPI_Datatype para_type;
MPI_Datatype ysend,xsend;



/*****************************************/
/* Fortran                               */

void hstcrt(double *rmin,double *rmax,int *nx,int *mbdcnd,double *bdra, double *bdrb,
	   double *thetamin, double *thetamax, int *ny,int *nbdcnd,double *bdrc, double *bdrd,
	    double *elmbda1,double *res,int *idimf,double *pertrb, int *ierror,double *workspace);
void hstplr(double *rmin,double *rmax,int *nx,int *mbdcnd,double *bdra, double *bdrb,
	   double *thetamin, double *thetamax, int *ny,int *nbdcnd,double *bdrc, double *bdrd,
	    double *elmbda1,double *res,int *idimf,double *pertrb, int *ierror,double *workspace);
void genbun(int *nperod,int *ny,int *mperod,int *nx,double *aphi,double *bphi_d,double *cphi,int *idimf,
           double *res,int *ierror,double *workspace);
void mgrid(int *key,double *local_field,double *local_x,double *local,int *m1,int *m2,int *m3,int *m4,
                double *expos,double *eypos, double *amplitude);
void dgtsvn(int *nx, double *dl,double *dc,double*du, double *res,int *stride, int *info);
void dgtsl(int *nx, double *dl,double *dc,double*du, double *res,int *stride, int *info);

void farakawa_45(double *res,double *a,double *b,double *fac,int *nx,int *ny);
void fcgtsvn(int *nx, double *dl,double *dc,double *du, double *res);
int multi_( double *res, double *trafo,double *tmp,int *nx,int *ny);
void Util_PrintStructures(HDF_DS *d,PARA *p,int line_no);
void Util_TrimString(char *str);

#ifdef __cplusplus
}
#endif




#endif /* _UTILITIES_H_ */
