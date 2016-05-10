/**************/
/* FFT.H      */
/* 15.09.1991 */
/**************/

#include <math.h>
#include <stddef.h>

/**************/
/* Prototypen */
/**************/

int bitrev(int v,int l);
void old_realft(double *data,int n,int isign);
void four1(double *data,int nn,int isign);
void realft(double *data,int n,int isign);
void four1c(double *DataRe, double *DataIm, int m, int isign);
void four2d(double **field, int nx, int ny);
void four2dre(double **field, int nx, int ny);

void powfour1(double **data,double **data_in,int nn,int ny,
			double *wrX, double *wiX,int max_xy,int isign);
void powfour1c(double *DataRe, double *DataIm, int m, int isign);
void powfour2d(double **field, int nx, int ny);
void powfour2dre(double **field, int nx, int ny);
int  powinitfft(int nx,int ny,double ***tempX,double **wrX,double **wiX,
        	  int **bitmask,double **tempr,double **tempi,int *max_xy);
void powrealft(double **data_out,double **data,int n,int ny,double *wrX,
                double *wiX,int max_xy,int isign);



void sinft(double *datar,int nx); 
void cosft(double *datar,int nx,int isign); 
void evenft2d(double **field, int nx, int ny);
void evenft2dre(double **field, int nx, int ny); 
void oddft2d(double **field, int nx, int ny);
void oddft2dre(double **field,int nx,int ny); 
void initeoft(int nx,int ny); 

void rft_2d_1d(double **f,int nx,int ny,int dir);




void vrffti(int *,double *); 
void vrfftf(int *nx,int *ny, double *f,double *r,int *MDIRM,double *wsave); 
void vrfftb(int *nx,int *ny, double *f,double *r,int *MDIRM,double *wsave); 




void vrft_2d_1d(double **r,double **f,int nx,int ny,int dir);
 void vrealft(double *r,double *f,int n,int dir);

