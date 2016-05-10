/******************************************/
/*   Allocating Memory                    */

#include <stdlib.h>

void *Util_AllocChk(size_t, size_t);
double *Util_DVector(int nx, int offx);
float *Util_FVector(int nx, int offx);
int *Util_IVector(int nx,int offx);
    
double **Util_DMatrix(int ny,int offy,int nx,int offx);
unsigned char **Util_UcMatrix(int ny,int offy,int nx,int offx);
float **Util_FMatrix(int ny,int offy,int nx,int offx);
int **Util_IMatrix(int ny,int offy,int nx,int offx);
    
double ***Util_DCube(int nx,int offx,int ny, int offy,int nz,int offz,int shift);
float ***Util_FCube(int nx,int offx,int ny, int offy,int nz,int offz);
int  ***Util_ICube(int nx,int offx,int ny,int offy,int nz,int offz, int shift);
    
void Util_FreeDCube(double ***field,int nx,int offx,int ny,int offy,int nz, int offz);
void Util_FreeDMatrix(double **field,int outside);

void Util_FreeFMatrix(float **field,int outside);
    
void Util_FreeIMatrix(int **field,int outside);
    
void Util_FreeUcMatrix(unsigned char **field,int outside);
    
void Util_3DAllocFields(double ****w_0,double ****w_1,double ****w_2,int sh1,
                        double ****dtw_0,double ****dtw_1,double ****dtw_2,int sh2,
                        double ***rbdra,double ***rbdrb,double ***zbdra,double ***zbdrb,
                        int nx,int ny,int nz,int offx,int offy,int offz);
    
void Util_2DAllocFields(double ***w_0,double ***w_1,double ***w_2,
                        double ***dtw_0,double ***dtw_1,double ***dtw_2,
                        double **rbdra,double **rbdrb,
                        int nx,int ny,int offx,int offy);
