#ifdef MPI
#include <mpi.h>
#endif 
#ifndef METRIC_H_ 
#define METRIC_H_

void Metric_IniMagneticField(double *qprof,double *qs,double *qss, double *kai,double *vcor,
			       double *rcor,double *edrcor,HDF_DS *data,PARA *para);



void Metric_ChangeGeometry(double *hval,double *vval, double *gm2vkm1,
			     double *kai,double *vcor, HDF_DS *data, PARA *para);


void Metric_CircularEquilibrium(double **grr, double **gff, double **gzz,
			double **grf, double **grz, double **gfz,
			double **gffm1, double **norm_gr, double **norm_gp,
			double **grrgfrmgrzgrf,double **grzgffmgrfgfz,
			double *theta,double *theta_shift,double *theta_trafo,double *norm,
			double *qprof,double *qs, double *kai, double *vcor,double *rcor, 
			HDF_DS *data,PARA *para);

void Metric_IniCurvature(double **kr, double **kf,double **kz, double *kai,
			  double *hval,double *vval, double *vcor, 
			  double *rcor,double *theta,double *theta_shift,
			  double **grr,double **gff,double **grf,
			  double **grz,double **gfz,
			  HDF_DS *data,PARA *para);

void Metric_CoeffCurvaturePassive(double **grr, double **gff, double **gzz,
                                  double **grf, double **grz, double **gfz,
                                  double **gffm1, double **norm_gr, double **norm_gp,
                                  double *theta,double *theta_sine,double *theta_trafo,double *norm,
                                  double *kai,double *vcor,double *rcor,double *edrcor, 
                                  double ***B_0, double ***Br,HDF_DS *data,PARA *p);


void Metric_Hamada2Glasser(double **res,double ***f,double *qprof,
		       double *theta_shift, HDF_DS *data,PARA *para, double dir);

void Metric_CoeffCurvatureITG(double *qprof,double *qs,double *qss, 
                                  double *kai,double *vcor,double *rcor,double *edrcor, 
                              double ***B_0, double ***Br,HDF_DS *data,PARA *p);
#endif       /*  !METRIC_H_ */
