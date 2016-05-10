
void  Util_1DSs1(double *w, double *w_0,double *r_0,int  nx, double dt)
{
    register int ix;
    for(ix=0;ix<nx;ix++)
        w[ix] = w_0[ix] + dt*r_0[ix];
}

/***********************************************************************/

void  Util_1SSs2(double *w, double *w_0,double *w_1,double *r_0,double *r_1,int  nx, double dt)
{
    register int ix;
    for(ix=0;ix<nx;ix++)
        w[ix] = 2./3.*(2.*w_0[ix] -0.5*w_1[ix]  +dt*(2.* r_0[ix] -r_1[ix]));
}

/***********************************************************************/
void  Util_1DSs3(double *w, double *w_0,double *w_1,double *w_2,
                 double *r_0,double *r_1,double *r_2, int nx,double dt)

{
    register int ix;
    for(ix=0;ix<nx;ix++)
        w[ix] =  6./11.*
            (3.*w_0[ix] -1.5*w_1[ix] +1./3.*w_2[ix]+dt*(3.*r_0[ix] -3.*r_1[ix] +r_2[ix]) );  

}

/***********************************************************************/

void  Util_Ss1(double **w, double **w_0,double **r_0,int  ny,int  nx, double dt)
{
    register int ix,iy;
    for(iy=0;iy<ny;iy++) for(ix=0;ix<nx;ix++)
        w[iy][ix] = w_0[iy][ix] + dt*r_0[iy][ix];
}

/***********************************************************************/

void  Util_Ss2(double **w, double **w_0,double **w_1,double **r_0,double **r_1,int  ny,int  nx, double dt)
{
    register int ix,iy;
    for(iy=0;iy<ny;iy++) for(ix=0;ix<nx;ix++)
        w[iy][ix] = 2./3.*(2.*w_0[iy][ix] -0.5*w_1[iy][ix]  +dt*(2.* r_0[iy][ix] -r_1[iy][ix]));
}

/***********************************************************************/
void  Util_Ss3(double **w, double **w_0,double **w_1,double **w_2,
               double **r_0,double **r_1,double **r_2,int ny, int nx,double dt)

{
    register int ix,iy;
    for(iy=0;iy<ny;iy++) for(ix=0;ix<nx;ix++)
        w[iy][ix] =  3.*w_0[iy][ix] -1.5*w_1[iy][ix] +1./3.*w_2[iy][ix]+dt*(3.*r_0[iy][ix] -3.*r_1[iy][ix] +r_2[iy][ix]);  


    for(iy=0;iy<ny;iy++) for(ix=0;ix<nx;ix++)
        w[iy][ix] *=6./11.;

}
void Util_1DSsTimeStep(int iter,int nx, double dt,double mue,double *lamda,
                     double *w,double *w_0, double *w_1, double *w_2,
                     double *r_0, double *r_1, double *r_2)
{
  switch (iter) 
    {
    case 0:
      Util_1DSs1(w,w_0,r_0,nx,dt);
      break;
    case 1:
      Util_1SSs2(w,w_0,w_1,r_0,r_1,nx,dt);
      break;
    default:
      Util_1DSs3(w,w_0,w_1,w_2,r_0,r_1,r_2,nx,dt);
    }

if(mue != 0.0)
  {
  switch (iter) 
    {
    case 0:
      *lamda   = -1./(mue*dt);
      break;
    case 1:
      *lamda   = -3./2./(mue*dt);
      break;
    default:
      *lamda   =-11./6./(mue*dt);
     }
  }
else 
  *lamda = 0.;
}

void Util_2DSsTimeStep(int iter,int nx, int ny, double dt,double mue,double *lamda,
                       double **w,double **w_0, double **w_1, double **w_2,
                       double **r_0, double **r_1, double **r_2)
{
    switch (iter) 
    {
        case 0:
            Util_Ss1(w,w_0,r_0,ny,nx,dt);
            break;
        case 1:
            Util_Ss2(w,w_0,w_1,r_0,r_1,ny,nx,dt);
            break;
        default:
            Util_Ss3(w,w_0,w_1,w_2,r_0,r_1,r_2,ny,nx,dt);
    }

    if(mue > 0.0)
    {
        switch (iter) 
        {
            case 0:
                *lamda   = -1./(mue*dt);
                break;
            case 1:
                *lamda   = -3./2./(mue*dt);
                break;
            default:
                *lamda   =-11./6./(mue*dt);
        }
    }
    else 
        *lamda = 0.;

}



/*********************************************************/

/*
  does the Stiffly Stable Timestepping 
*/

void Util_3DSsTimeStep(int iter,int nx, int ny, int nz, double dt,double mue,double *lamda,
                       double ***w,double ***w_0, double ***w_1, double ***w_2,
                       double ***r_0, double ***r_1, double ***r_2)
{
    static int i;

    switch (iter) 
    {
        case 0:
            for(i=0;i<nz;i++)
                Util_Ss1(w[i],w_0[i],r_0[i],ny,nx,dt);
            break;
        case 1:
            for(i=0;i<nz;i++)
                Util_Ss2(w[i],w_0[i],w_1[i],r_0[i],r_1[i],ny,nx,dt);
            break;
        default:
            for(i=0;i<nz;i++)
                Util_Ss3(w[i],w_0[i],w_1[i],w_2[i],r_0[i],r_1[i],r_2[i],ny,nx,dt);
    }


    if(mue != 0.0)
    {
        switch (iter) 
        {
            case 0:
                *lamda   = -1./(mue*dt);
                break;
            case 1:
                *lamda   = -2./(3.*mue*dt);
                break;
            default:
                *lamda   =-11./(6.*mue*dt);
        }
    }
    else 
        *lamda = 0.;

}
