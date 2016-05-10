/**************************************************************/

void Util_BdCylA(double **f,double *val,double *hval, int nr,int np,int shifted)
{
    register int ip;
    int np2;

    np2 = np/2;

    for(ip=0;ip<np2;ip++) f[ip][-1] = f[ip+np2][0];
    for(ip=0;ip<np2;ip++) f[ip+np2][-1] = f[ip][0];
}
/**************************************************************/

void Util_CylBdDir(double **f,int np,int nr,double *bdra,double *bdrb,double *hval)
{
    Util_BdCylA(f,NULL,0,nr,np,FALSE);
    Util_BdDirB(f,bdrb,hval,nr,np,FALSE);
    Util_BdPer(f,NULL,NULL,0,nr,np,FALSE);
}


/***********************************************************************/

void Util_CylBdNeu(double **f,int np,int nr,double *bdra,double *bdrb,double *hval)
{
    register int ip;
    int np2;

    np2 = np/2;

    for(ip=-1;ip<np2;ip++) f[ip][-1]     = f[ip+np2][0];
    for(ip=0;ip<=np2;ip++) f[ip+np2][-1] = f[ip][0];

    for(ip=-1;ip<=np;ip++) f[ip][nr] = f[ip][nr-1];
    Util_BdPer(f,NULL,NULL,NULL,nr,np,FALSE);

}


/* 
   Sets up Ghost points for Neumann_boundaries on a non-staggerd grid at x=xmin or x=xmax
   Arguments are:
   f = field to add boundaries organized as f[y][x]
   vala = value of derivative at boundary xmin
   valb = value of derivative at boundary xmax
   delta = grid spacing
   nx = gridpints in x
   ny = gridpoints in y
   shifted = Boolean, indicates if grid is shifted. True for RIGTH approximated values 
*/

   
void Util_BdNeuNs(double **f,double *vala, double *valb,double *delta, int nx,int ny,int shifted)
{

    register int i;
  
    if(shifted)
    {   
        for(i=-1;i<=ny;i++) f[i][nx-1] = (1./3.) *(4.*f[i][nx-2] - f[i][nx-3] + 2.*delta[nx-1]*valb[i]);
    }
    else
    {
        for(i=-1;i<=ny;i++) f[i][-1] = (1./3.)   *(4.*f[i][0] - f[i][1] - 2.*delta[-1]*vala[i]);
    }
}


void Util_BdDirA(double **f,double *val,double *hval, int nx,int ny,int shifted)
{
    register int i;
  
    for(i=-1;i<=ny;i++) f[i][-1] = 2.*val[i] - f[i][0];

    if(shifted)  f[ny-1][nx]= f[-1][nx];   
    else 
    {
        f[-1][-1] = f[ny-1][-1];
        f[ny][-1] = f[0][-1];
    }
}

/**************************************************************/   
void Util_BdDirB(double **f,double *val,double *delta, int nx,int ny,int shifted)
{
    register int i; 
 
    for(i=-1;i<=ny;i++) f[i][nx] =  2.*val[i] - f[i][nx-1];

    if(shifted)  f[ny-1][nx]= f[-1][nx];
    else 
    {
        f[-1][nx] = f[ny-1][nx];
        f[ny][nx] = f[0][nx];
    }
   
}




/**************************************************************/

void Util_BdPerX(double **f,double *val,double *hval, int nx,int ny,int shifted)
{
    register int i; 
    for (i=-1; i<=ny; i++)
	{
        f[i][-1] = f[i][nx-1];
        f[i][nx] = f[i][0];
	}

}

/**************************************************************/

void Util_BdPerZ(double ***f,double **val,double hval, int nx,int ny,int nz)
{
    register int i,j;

    for (i=-1; i<=ny; i++)
        for (j=-1; j<=nx; j++)
        {
            f[-1][i][j] = f[nz-1][i][j];
            f[nz][i][j] = f[0][i][j];
        }

}

/* 
   Sets up Ghost points for Neumann_boundaries on a staggered grid at x=xmin or x=xmax
   Arguments are:
   f = field to add boundaries organized as f[y][x]
   vala = value of derivative at boundary xmin
   valb = value of derivative at boundary xmax
   delta = grid spacing
   nx = gridpints in x
   ny = gridpoints in y
   shifted = Boolean, indicates if grid is shifted. True for RIGTH approximated values 
*/

   
void Util_BdNeuA(double **f,double *val,double *hs, int nx,int ny,int shifted)
{
    register int i;
    for(i=-1;i<=ny;i++) f[i][-1] = -val[i]/hs[0]+f[i][0];
    if(shifted)  f[ny-1][nx]= f[-1][nx];  
    else 
    {
        f[-1][-1] = f[ny-1][-1];
        f[ny][-1] = f[0][-1];
    }  
}

   
void Util_BdNeuB(double **f,double *val,double *hs, int nx,int ny,int shifted)
{
    register int i; 
    for(i=-1;i<=ny;i++) f[i][nx] =  val[i]/hs[nx-1]+f[i][nx-1];
    if(shifted)  
        f[ny-1][nx]= f[-1][nx]; 
    else 
    {
        f[-1][nx] = f[ny-1][nx];
        f[ny][nx] = f[0][nx];
    }
}

/**************************************************************/


/*
  Set up ghost points for peridic boundaries a y=ymin and y = ymax

  Arguments are:
  f = field to add boundaries organized as f[y][x]
  vala = dummy value
  valb = dummy value
  delta = dummy value
  nx = gridpints in x
  ny = gridpoints in y
  shifted = Boolean, indicates if grid is shifted in case of Util_UpWinding 
*/

void Util_BdPer(double **f,double *vala, double *valb, double *delta,int nx,int ny,int shifted)
{
    /* Shifted is true for RIGHT approximated in y-direction */ 

    register int i;
    int offr=ny,offl=0;

    if(shifted) 
    {
        offr = 0; 
        offl = ny;
    }

    BUGREPORT;

    for(i=-1;i<=nx;i++) f[offl-1][i] = f[offr-1][i];
    for(i=-1;i<=nx;i++) f[ny][i]     = f[0][i];
    

}
void Util_2DFullBd(double **u,int ny,int nx,double *val1,double *val2,double *hval, int type)
{


    switch (type)
    {
        case 0:
            Util_BdPerX(u,val1,hval,nx,ny,FALSE);
            break;
        case 1:
            Util_BdDirA(u,val1,hval,nx,ny,FALSE);
            Util_BdDirB(u,val2,hval,nx,ny,FALSE);
            break;
        case 2:
            Util_BdDirA(u,val1,hval,nx,ny,FALSE);
            Util_BdNeuB(u,val2,hval,nx,ny,FALSE);
            break;
        case 3:
            Util_BdNeuA(u,val1,hval,nx,ny,FALSE);
            Util_BdNeuB(u,val2,hval,nx,ny,FALSE);
            break;
        case 4:
            Util_BdNeuA(u,val1,hval,nx,ny,FALSE);
            Util_BdDirB(u,val2,hval,nx,ny,FALSE);
            break;
        case 5:
            Util_BdCylA(u,val1,hval,nx,ny,FALSE);
            Util_BdDirB(u,val2,hval,nx,ny,FALSE);
            break;
        case 6:
            Util_BdCylA(u,val1,hval,nx,ny,FALSE);
            Util_BdNeuB(u,val2,hval,nx,ny,FALSE);   
            break;
        case 7:
            Util_BdDirA(u,val1,hval,nx,ny,FALSE);
            Util_BdDirB(u,val2,hval,nx,ny,FALSE);
            break;
    }

    Util_BdPer(u,NULL,NULL,NULL,nx,ny,FALSE);

}
/**********************************************************************************/



void Util_3DFullBd(double ***u,HDF_DS *d,int nx,int ny, int nz,
                   double **left,double **right,double **up, double **down,
                   double *dx,double dz, int typer, int typez)
{
    int i,j,k;

    double val;
 
    val = 1./(double)ny;
    

    switch (typez)
    {
        case 0:
            Util_BdPerZ(u,up,dz,nx,ny,nz);
            break;
        case 1:
            Util_3DBdDirA(u,down,d,dz,nx,ny,nz);
            Util_3DBdDirB(u,up,d,dz,nx,ny,nz);
            break;
        case 2:
            Util_3DBdDirA(u,down,d,dz,nx,ny,nz);
            Util_3DBdNeuB(u,up,d,dz,nx,ny,nz);
            break;
        case 3:
            Util_3DBdNeuA(u,down,d,dz,nx,ny,nz);
            Util_3DBdNeuB(u,up,d,dz,nx,ny,nz);
            break;
        case 4:
            Util_3DBdNeuA(u,down,d,dz,nx,ny,nz);
            Util_3DBdDirB(u,up,d,dz,nx,ny,nz);
            break;
        case 13:
            for(i=0;i<ny;i++) for(j=0;j<nx;j++)
                down[i][j]=up[i][j]=0.;
            for(i=0;i<ny;i++) for(j=0;j<nx;j++)
                down[0][j]+= u[0][i][j];
            for(i=0;i<ny;i++) for(j=0;j<nx;j++)
                up[0][j]+= u[nz-1][i][j];

        
            for(j=0;j<nx;j++) down[0][j]*=val;
            for(j=0;j<nx;j++) up[0][j]*=val;
        
            for(i=1;i<ny;i++) for(j=0;j<nx;j++)
                up[i][j] = up[0][j];    
            for(i=1;i<ny;i++) for(j=0;j<nx;j++)
                down[i][j] = down[0][j];    

            Util_3DBdDirA(u,down,d,dz,nx,ny,nz);
            Util_3DBdDirB(u,up,d,dz,nx,ny,nz);
            break;
        default:
            fprintf(stderr,"Z-BDRCND %d not implemented!!!!!\n",typez);
            break;
    }


    for(i = -d->offz; i< nz+d->offz;i++)
        Util_2DFullBd(u[i],ny,nx,left[i],right[i],dx,typer);

}
/**************************************************************/

void Util_3DBdDirA(double ***f,double **val,HDF_DS *d,double hval, int nx,int ny,int nz)
{
    register int ip,ir;
 
    for(ip=0;ip<ny;ip++) 
        for(ir=0;ir<nx;ir++)
            f[-1][ip][ir] = 2.*val[ip][ir]-f[0][ip][ir];
 
    if(d->offz == 2)
        for(ip=0;ip<ny;ip++) 
            for(ir=0;ir<nx;ir++)
                f[-2][ip][ir] =  2.*val[ip][ir] - f[1][ip][ir];


}
/**************************************************************/
void Util_3DBdDirB(double ***f,double **val,HDF_DS *d,double hval, int nx,int ny,int nz)
{
    register int ip,ir;
 
    for(ip=0;ip<ny;ip++) 
        for(ir=0;ir<nx;ir++)
            f[nz][ip][ir] = 2.*val[ip][ir]-f[nz-1][ip][ir];

    if(d->offz == 2)
        for(ip=0;ip<ny;ip++) 
            for(ir=0;ir<nx;ir++)
                f[nz+1][ip][ir] =  2.*val[ip][ir] - f[nz-2][ip][ir];
}

/**************************************************************/

void Util_3DBdNeuB(double ***f,double **val,HDF_DS *d,double hval, int nx,int ny,int nz)
{
    register int ip,ir;
 
    for(ip=0;ip<ny;ip++) 
        for(ir=0;ir<nx;ir++)
            f[nz][ip][ir] = val[ip][ir]/hval + f[nz-1][ip][ir];

    if(d->offz == 2)
        for(ip=0;ip<ny;ip++) 
            for(ir=0;ir<nx;ir++)
                f[nz+1][ip][ir] =   f[nz][ip][ir];


}

/****************************************************************/

void Util_3DBdNeuA(double ***f,double **val,HDF_DS *d,double hval, int nx,int ny,int nz)
{
    register int ip,ir;
 
  
   
    if(d->offz == 2)
        for(ip=0;ip<ny;ip++) 
            for(ir=0;ir<nx;ir++)
                f[-2][ip][ir] =   f[-1][ip][ir];
 
    for(ip=-1;ip<=ny;ip++) 
        for(ir=-1;ir<=nx;ir++)
            f[-1][ip][ir] =  - val[ip][ir]/hval + f[0][ip][ir];

}

