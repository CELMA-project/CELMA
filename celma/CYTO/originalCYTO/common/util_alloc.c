#include <stdio.h>
#include <stdlib.h>
#include "util_alloc.h"

#define CALLOC Util_AllocChk

/*****************************************************************************/ 

void *Util_AllocChk(size_t n,size_t size)
{
    void *adr;
    
    if((adr=calloc(n,size))== 0 ) 
    {
        fprintf(stderr,"Allocation of %d elements of size %d failed\n",(int)n, (int)size );
        exit(0);
    }
    return adr;
}

/***************************************************************************/

double *Util_DVector(int n,int off)
{
    return ((double*) CALLOC((size_t)(n+2*off),sizeof(double) ))+off;
}
/***************************************************************************/

unsigned char **Util_UcMatrix(int ny,int offy,int nx,int offx)
{
/* y =sizey; x = sizex 		*/
/*
  Stellt ein Feld zur Verfeugung dessen letzter Index 
  nx-1 (bzw. nx+OFFSET-1) und dessen erster Index
  0  ( bzw. -OFFSET) ist.
  D.h. wir muessen (nx+1+2*OFFSET)*(ny+1+2*OFFSET)
  Elemente allokieren.
*/

    int		i;
    size_t     	sizey,sizex;
    unsigned char   **field;

    sizey = (size_t)(ny + 2*offy);
    sizex = (size_t)(nx + 2*offx);

    field 	= (unsigned char **)CALLOC(sizey,sizeof(unsigned char *))+offy;
    field[-offy]= (unsigned char *)CALLOC(sizex*sizey,sizeof(unsigned char))+offx;

    for(i = 1-offy;i<ny+offy;i++)   field[i]=field[i-1]+sizex;

    return field;

}

/***************************************************************************/

int **Util_IMatrix(int ny,int offy,int nx, int offx)
{

    int	i;
    size_t  sizey,sizex;
    int   	**field;

    sizey = (size_t)(ny + 2*offy);
    sizex = (size_t)(nx + 2*offx);

    field 	= (int **)CALLOC(sizey,sizeof(int *))+offy;
    field[-offy]= (int *)CALLOC(sizex*sizey,sizeof(int))+offx;

    for(i = 1-offy;i<ny+offy;i++)  field[i]=field[i-1]+sizex;

    return field;

}


/***************************************************************************/
double **Util_DMatrix(int ny,int offy,int nx,int offx)
{
/* y =sizey; x = sizex 		*/
/*
  Stellt ein Feld zur Verfeugung dessen letzter Index 
  nx-1 (bzw. nx+OFFSET-1) und dessen erster Index
  0  ( bzw. -OFFSET) ist.
  D.h. wir muessen (nx+2*OFFSET)*(ny+2*OFFSET)
  Elemente allokieren.
*/

    int		i;
    size_t     	sizey,sizex;
    double   	**field;

    sizey = (size_t)(ny  + 2*offy);
    sizex = (size_t)(nx  + 2*offx );

    field 	    = (double **)CALLOC(sizey,      sizeof(double *))+offy;
    field[-offy]= (double  *)CALLOC(sizex*sizey,sizeof(double)  )+offx;

    for(i = 1-offy;i<ny+offy;i++) field[i]=field[i-1]+sizex;

    return field;
}
/***************************************************************************/
void Util_FreeIMatrix(int **field,int off)
{
    free(field[-off]-off);
    free(field-off);
    field = (int **)NULL;	
}
/***************************************************************************/

float **Util_FMatrix(int ny,int offy,int nx,int offx)
{
/* y =sizey; x = sizex 		*/

    int		i;
    size_t     	sizey,sizex;
    float   	**field;

    sizey = (size_t)(ny + 2*offy);
    sizex = (size_t)(nx + 2*offx);

    field 	= (float **)CALLOC(sizey,sizeof(float *))+offy;
    field[-offy]= (float *)CALLOC(sizex*sizey,sizeof(float))+offx;

    for(i = 1-offy;i<ny+offy;i++)	field[i]=field[i-1]+sizex;
    return field;

}

/***************************************************************************/

double ***Util_DCube(int nx,int offx,int ny,int offy,int nz,int offz, int shift)
{
/* y =sizey; x = sizex 		*/
/*
  Allocates field with last index
  nx-1 (nx+OFFSET-1) 
  and first index
  0  (-OFFSET)
  This totals (nx+2*OFFSET)*(ny+2*OFFSET) elements.
*/

    int  		i,j;
    size_t     	sizey,sizex,sizez;
    double   	***field;

    sizez = (size_t)(nz  + 2*offz);
    sizey = (size_t)(ny  + 2*offy);
    sizex = (size_t)(nx  + 2*offx);

    field 	           = (double ***)CALLOC(sizez,sizeof(double **));
    field             += offz;
    
    field[-offz]       = (double **) CALLOC(sizez*sizey,sizeof(double *));
    field[-offz]      += offy;
    


    field[-offz][-offy]= (double *)  CALLOC((sizez*sizey*sizex+(size_t)shift),sizeof(double));
    field[-offz][-offy]+= offx+shift;
    

    for(i = 1-offz;i<nz+offz;i++)  field[i]=field[i-1]+sizey;

    for(i = 1-offz;i<nz+offz;i++)  field[i][-offy] = field[i-1][-offy]+sizey*sizex;

    for(i = -offz;i<nz+offz;i++) for(j = 1-offy;j<ny+offy;j++) field[i][j]=field[i][j-1]+sizex;


    return field;
}


int  ***Util_ICube(int nx,int offx,int ny,int offy,int nz,int offz, int shift)
{
/* y =sizey; x = sizex 		*/
/*
  Allocates field with last index
  nx-1 (nx+OFFSET-1) 
  and first index
  0  (-OFFSET)
  This totals (nx+2*OFFSET)*(ny+2*OFFSET) elements.
*/

    int  		i,j;
    size_t     	sizey,sizex,sizez;
    int   	***field;

    sizez = (size_t)(nz  + 2*offz);
    sizey = (size_t)(ny  + 2*offy);
    sizex = (size_t)(nx  + 2*offx);

    field 	= (int ***)CALLOC(sizez,sizeof(int **))+offz;
    field[-offz]= (int **)CALLOC(sizez*sizey,sizeof(int *))+offy;
    field[-offz][-offy]= (int *)CALLOC(sizez*sizey*sizex,sizeof(int))+offx+shift;

    for(i = 1-offz;i<nz+offz;i++) 
        field[i]=field[i-1]+sizey;


    for(i = 1-offz;i<nz+offz;i++) 
        field[i][-offy] = field[i-1][-offy]+sizey*sizex;
 

    for(i = -offz;i<nz+offz;i++)
        for(j = 1-offy;j<ny+offy;j++)
            field[i][j]=field[i][j-1]+sizex;


    return field;
}


/***************************************************************************/

void Util_FreeDCube(double ***field,int nx,int offx,int ny,int offy,int nz,int offz)
{
    free (field[-offz][-offy]-offx);
    free (field[-offz]-offy);
    free (field-offz);

    field = (double***)NULL;	
}


/******************************************************************************/
/* This function allocates a number of 3d fields and boundaries               */

void Util_3DAllocFields(double ****w_0,double ****w_1,double ****w_2,int sh1,
                        double ****dtw_0,double ****dtw_1,double ****dtw_2,int sh2,
                        double ***rbdra,double ***rbdrb,double ***zbdra,double ***zbdrb,
                        int nx,int ny,int nz,int offx,int offy,int offz)
{
    *w_0  = Util_DCube(nx,offx,ny,offy,nz,offz,sh1);
    *w_1  = Util_DCube(nx,offx,ny,offy,nz,offz,sh1);
    *w_2  = Util_DCube(nx,offx,ny,offy,nz,offz,sh1);
  
    *dtw_0  = Util_DCube(nx,offx,ny,offy,nz,offz,sh2);
    *dtw_1  = Util_DCube(nx,offx,ny,offy,nz,offz,sh2);
    *dtw_2  = Util_DCube(nx,offx,ny,offy,nz,offz,sh2);

    /* Allocate boundaries     */
  
    *rbdra = Util_DMatrix(nz,offz,ny,offy);
    *rbdrb = Util_DMatrix(nz,offz,ny,offy);

    *zbdra  = Util_DMatrix(ny,offy,nx,offx);
    *zbdrb  = Util_DMatrix(ny,offy,nx,offx);
}  


/***************************************************************************/

void Util_FreeUcMatrix(unsigned char **field,int off)
{
    if(field[-off] != (unsigned char*)NULL)
	{
        free((unsigned char *)field[-off]-off);
        if(field-off != (unsigned char**)NULL)
            free((unsigned char**)field-off);
	}
    field = (unsigned char**)NULL;	
}

/***************************************************************************/

void Util_FreeDMatrix(double **field,int off)
{
    if(field[-off]) free(field[-off]-off);
    if(field-off) free(field-off);
    field = (double**)NULL;	
}
/***************************************************************************/

void Util_FreeFMatrix(float **field,int off)
{
    if(field[-off]-off) free(field[-off]-off);
    if(field-off) free(field-off);
    field = (float**)NULL;	
}

/***************************************************************************/


void Util_2DAllocFields(double ***w_0,double ***w_1,double ***w_2,
                        double ***dtw_0,double ***dtw_1,double ***dtw_2,
                        double **rbdra,double **rbdrb,
                        int nx,int ny,int offx,int offy)
{
    *w_0  = Util_DMatrix(ny,offy,nx,offx);
    *w_1  = Util_DMatrix(ny,offy,nx,offx);
    *w_2  = Util_DMatrix(ny,offy,nx,offx);
  
    *dtw_0  = Util_DMatrix(ny,offy,nx,offy);
    *dtw_1  = Util_DMatrix(ny,offy,nx,offy);
    *dtw_2  = Util_DMatrix(ny,offy,nx,offx);
  
    /* Allocate boundaries     */
  
    *rbdra = Util_DVector(ny,offy);
    *rbdrb = Util_DVector(ny,offy);
}
