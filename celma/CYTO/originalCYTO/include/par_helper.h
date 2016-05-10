#include <mpi.h>
#ifndef         PARHELPER_H_
#define         PARHELPER_H_


#define ROOT 0
#define ISROOT   (d->this_process == ROOT)

/* Include File for parallel Helpers */

/* Communication of structures */

MPI_Datatype create_mpi_hdfds(HDF_DS *d, int root, MPI_Comm comm,int rank);
MPI_Datatype create_mpi_para(PARA *d, int root, MPI_Comm comm,int rank);

#ifdef __cplusplus
extern "C" {
#endif
/* Integrals and the like */

double  PH_Integral(double **f, int order, double *norm,int ny,int nx);
void PH_3D_Geometry(HDF_DS *data);
void PH_3D_FS_Average(double ***f,double *result,HDF_DS *data);
void PH_3D_n0_Component(double ***f,double **result_n_0,HDF_DS *data);
void PH_3D_n0m0_Component(double ***f,double *result_n_0,HDF_DS *data);



int PH_CheckCFL3D(double ***f_0,double ***vr,double ***vp,HDF_DS *data,PARA *p,
			  double **hval,double **vval,double *cflr,double *cflp);

int PH_CheckCFL3D_Metric(double ***f_0,double ***vr,double ***vp,HDF_DS *data,PARA *p,
			   double **hval,double **vval,
			   double *kai, double **grr, double **gff, double **grf,
	  		   double **norm_gr, double **norm_gp,
			   double *cflr,double *cflp);

double PH_3DIntegral(double ***f, int order, double **norm, int nz,int ny,int nx);
double PH_3DIntegral_noreduce(double ***f, int order, double **norm, int nz,int ny,int nx);



/* Read Write HDF  */

int PHInt_ReadHDFbyNumber(int number, HDF_DS *data,PARA *para);
int PHInt_ReadHDFbyName(const char *name, HDF_DS *data,PARA *para);


int  PH_Read3DFieldbyName(double ***f_0,const char *name,HDF_DS *data,PARA *para);
int  PH_Read3DFieldbyNumber(double ***f_0,int number,HDF_DS *data,PARA *para);
int  PH_Read2DFieldbyName(double **f_0,const char* name,HDF_DS *data,PARA *para);
int  PH_Read2DFieldbyNumber(double **f_0,int number,HDF_DS *data,PARA *para);

void PH_3D_Write(double ***f_0,const char *fname,const char *name,int number, HDF_DS *data,PARA *p,int newfile);

/* Metric associated functions */
void PH_Untwist_Local(double **f,double **res,double *qprof,double dtheta,int dist, HDF_DS *data,int dir);

void PH_Update2dBoundaries(double ***f,int mbdcnd, double **a, double **b,
                           double **hval,HDF_DS *data);

void PH_UpdateZ(double ***f,HDF_DS *data,MPI_Comm comm);

void PH_UpdateZ_FourierField(double ***f,HDF_DS *data,MPI_Comm comm);


void PH_UpdateY(double ***f, HDF_DS *data);
void PH_UpdateX(double ***f, HDF_DS *data);
void PH_UpdateZ_fft(double ***f,HDF_DS *data,MPI_Comm comm);

void  PH_write_radial_profile(double *bar,char *name,HDF_DS *d);
void  PH_write_axial_profile(double *bar,char *name,HDF_DS *d);  
#ifdef __cplusplus
}
#endif

#endif       /*  !PARHELPER_H_ */
