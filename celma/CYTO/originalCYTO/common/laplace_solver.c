#undef  DEBUG
#undef COMM_STATUS
#include "utilities.h"



#ifdef FFTW
#include <fftw3.h>
#endif

/***********************************************************************/
/*!

Call 2D solvers.

*/
void Laplace_Solve3D(HDF_DS *d, PARA *p,double ***f, double ***w,
               double **bdra,double **bdrb,double *hval,double *rcor,int mbdcnd,
               double lamda,int normalize)
{
int i;

 BUGREPORT;

 for(i=0;i<d->lnz;i++)
 {

     COMM(fprintf(stderr,"Proc %d %s line %d line %d of %d\n",d->this_process,__FILE__,__LINE__,i, d->lnz););
     if (d->N[2]>1)
         Laplace_Solve2D_Par(d,p,f[i],w[i], bdra[i],bdrb[i],hval,rcor,mbdcnd,lamda,normalize);
     else
         Laplace_Solve2D(d,p,f[i],w[i], bdra[i],bdrb[i],hval,rcor,mbdcnd,lamda,normalize);
 }

}
/***********************************************************************/



/* loeiten: Actually documented :D.
 * Helmholtz eq: https://en.wikipedia.org/wiki/Helmholtz_equation
 * lambda here is what is called k in wiki
 */
/*!
  Helmholtz solver on a disk/annulus with a staggered radial grid.
  Uses fouriertransform in poloidal and
  Gaussian elimination in radial direction
  to solve 2nd order finite difference approximation for each m-mode:

   w_m = 1/r dr f_m - m2/r2  f_m + drr f_m + lamda f_m


   To solve the implizit part of the SS3  w* = (1 + lamda^-1 lap ) w

   We here solve:

   lamda w* = (lamda + lap) w

   and renormalize w* with lamda!

   Thus for solving (1+rho^2 lap )f = phi
   lamda has to be 1/rho^2.

  von V. Naulin 30.10.1997


  Boundary conditions

  MBDCND:      R = A           R = B
     0                PERIODIC
     1          dir              dir
     2          dir              neu
     3          neu              neu
     4          neu              dir
     5          R = 0            dir
     6          R = 0            neu


  Used external functions :

        rft_2d_1d from libfft
        LapSol_GaussElim  modified dgtsvn from lapack to solve 1d equation

   Function gets on input:
        rmin, rmax, nx, mdcnd:       radial coordinate info
        bdra, bdrb:                  double fields wich contain boundary values
        rcor:                        double field contains postions of radial points
        edx:                         contains 1/dr
        thetamin,thetamax,ny:        poloidal coordinate info
        lamda:                       parameter lambda of helmholtz equation
        w_0:                         linke seite
        normalize:                   boolean to indicate if result should be normalized with lamda
        res:                         workspace of same dimension as w_0, fw


   on output:
        fw

 */
/***************************************************************************/
void Laplace_Solve2D(HDF_DS *data, PARA *p,double **f, double **w,
                     double *ibdra,double *ibdrb,double *hval,
                     double *rcor,int mbdcnd,
                     double lamda, int normalize)

{
    register int ir,ip;
    static int nx,ny,ny_cut;

    static double *dl,*dc,*du,*rsquare,*ky2;
    static double *bdra,*bdrb;
    static double **res;
    static int FIRST=TRUE, ispolar=FALSE;
    static double *dl_temp, *dc_temp,*du_temp;

    double ky,ehedx,r2d,kymin;
    double r,r2,mult;


    // Do doubly periodic solver in subroutine
    if((strncmp(data->coordsys,"pol",3) == 0) || strncmp(data->coordsys,"cyl",3) == 0) ispolar = TRUE;


    if(!ispolar && mbdcnd == 0)
    {
        LapSol_PeriodicSolver(data, p,f, w,lamda, normalize);
    }
    else if(ispolar && mbdcnd == 0)
    {
        fprintf(stderr,"Forbidden periodicity in polar mode!!!!!\n"); exit(-1);
    }
    else
    {
        /*loeiten: Enters here
         *         Set static variables
         *         http://stackoverflow.com/questions/572547/what-does-static-mean-in-a-c-program
         */
        if(FIRST)
        {


            nx = data->lnx;
            ny = data->lny;
            ny_cut=   5*(ny/6);
            ny_cut=   ny+2;


            bdra = Util_DVector(ny+2,1);
            bdrb = Util_DVector(ny+2,1);
            ky2= Util_DVector(ny+2,2);

            rsquare = Util_DVector(nx,2);
            dl = Util_DVector(nx,2);
            dc = Util_DVector(nx,2);
            du = Util_DVector(nx,2);

            BUGREPORT;
            dl_temp   = Util_DVector(nx,2);
            dc_temp   = Util_DVector(nx,2);
            du_temp   = Util_DVector(nx,2);
            res     = Util_DMatrix((ny+2),0,(nx+2),0);
            p->dky =  2.0*M_PI/(p->ymax-p->ymin);

            for(ip=0; ip< (ny+2);  ip++) ky2[ip]= (double)((ip)/2)*p->dky;
            for(ip=0; ip< (ny+2);  ip++) ky2[ip]*=ky2[ip];


            BUGREPORT;

            if(ispolar)
            {
                COMM(fprintf(stderr,"Setting up for POLAR coordinate system in Laplace solve.\n"););

                for(ir=0;ir<nx;ir++)
                {
                    r     = rcor[ir];
                    rsquare[ir] = r*r;
                    r2    = rsquare[ir];
                    r2d   = r2*hval[ir]*hval[ir];
                    ehedx = hval[ir]*r*0.5;

                    dl_temp[ir] = -ehedx +    r2d;
                    dc_temp[ir] = - 2.*r2d;
                    du_temp[ir] =  ehedx +    r2d;
                    }
            }
            else
            {
                COMM(fprintf(stderr,"Setting up for RECTANGULAR coordinate system in Laplace solve.\n"); );

                for(ir=0;ir<nx;ir++)
                {
                        r2d   = hval[ir]*hval[ir];
                        rsquare[ir] = 1.;
                        dl_temp[ir] =   r2d;
                        dc_temp[ir] = - 2.*r2d;
                        du_temp[ir] =   r2d;
                }
            }

            BUGREPORT;
            COMM(fprintf(stderr,"Proc %d %s line FIRST:  %d: nx  %d ny %d\n",data->this_process,__FILE__,__LINE__,nx,ny););
        }


        if(!normalize || lamda == 0.)
            /*loeiten: Enters here*/
            mult = 1.;
        else
            mult = lamda;

        COMM(fprintf(stderr,"Proc %d %s: line %d nx = %d,%d \n",data->this_process,__FILE__,__LINE__,nx,ny););

        for(ip=0;ip<ny;ip++) for(ir=0;ir<nx;ir++) res[ip][ir] = w[ip][ir]*rsquare[ir]*mult;


        COMM(fprintf(stderr,"Proc %d %s: line %d \n",data->this_process,__FILE__,__LINE__););

        /* Something to remember:

        m = 0 in f[0],f[1]
        m = 1 in f[2],f[3]
        m = 2 in f[4],f[5] ...etc.
        */

        // FFT boundary AND field
        BUGREPORT;
        for(ip=0;ip<ny;ip++) res[ip][nx]=ibdra[ip];
        for(ip=0;ip<ny;ip++) res[ip][nx+1]=ibdrb[ip];

        Fft_1d_2d_f(res,nx+2,ny);

        for(ip=0;ip<ny+2;ip++)  bdra[ip]=res[ip][nx];
        for(ip=0;ip<ny+2;ip++)  bdrb[ip]=res[ip][nx+1];

        // Set last two rows to zero
        for(ir=0;ir<nx;ir++) res[ny][ir] = res[ny+1][ir] = 0.;

        // solve 1-D problems
        COMM(fprintf(stderr,"Proc %d line %d Boundary Condition %d\n",data->this_process,__LINE__,mbdcnd););


        for(ip=0;ip<ny_cut;ip++)
        {
            /*setup vectors of tridi matrix and multiply rhs by r2 */
            for(ir=0;ir<nx;ir++)
            {
                dl[ir] =       dl_temp[ir];
                dc[ir] = -ky2[ip] + dc_temp[ir] + rsquare[ir] * lamda;
                du[ir] =       du_temp[ir];
            }

            // Enter boundary values

            switch (mbdcnd)
            {
                case PERIODIC:
                    if(ispolar)
                    {
                      fprintf(stderr,"%s: Boundary %d not allowed with poloidal geometry!\n",__func__,mbdcnd);
                        exit (-1);
                    }
                    break;
                case DIRDIR:
                case DIRNEU:
                case 11:
                    dc[0]         -=  dl[0];
                    res[ip][0]    -=  2.*dl[0]*bdra[ip];
                    break;
                case NEUNEU:
                case NEUDIR:
                    /*fprintf(stderr,"%s: Boundary %d selected!\n",__func__,mbdcnd);*/
                    dc[0]         +=  dl[0];
                    res[ip][0]    +=  dl[0]/hval[0]*bdra[ip];
                    break;
                case ZeroDIR:
                case ZeroNEU:
                    /*loeiten: Enters here*/
                                        // zero mode has zero derivative
                  if(ip  == 0)
                                                        dc[0]         +=  dl[0];
                     else
                                        // other  modes have zero value
                                                        dc[0]         -=  dl[0];

                     break;
                default:
                    fprintf(stderr,"%s: Boundary %d not allowed!\n",__func__,mbdcnd);
                    exit(-1);
            }

            ir = nx-1;

            switch (mbdcnd)
            {
                case PERIODIC:
                    if(ispolar )
                    {
                      fprintf(stderr,"%s: Boundary %d not allowed with poloidal geometry!\n",__func__,mbdcnd);
                        exit (-1);
                    }
                    break;
                case DIRDIR:
                case NEUDIR:
                case ZeroDIR:
                    dc[ir]      -=  du[ir];
                    res[ip][ir] -=  2.*du[ir]*bdrb[ip];
                    break;
                case DIRNEU:
                case NEUNEU:
                case ZeroNEU:
                    dc[ir]      +=  du[ir];
                    res[ip][ir] -=  du[ir]/hval[ir]*bdrb[ip];
                    break;
                case 11: /*dirchlet for n != 0, neuman for n == 0 */
                    dc[ir]      -=  du[ir];
                    res[ip][ir] -=  2.*du[ir]*bdrb[ip];
                    if (ip == 0)
                    {
                        dc[ir]      +=  du[ir];
                        res[ip][ir] -=  0.;
                    }
                    break;
                default:
                    fprintf(stderr,"%s: Boundary %d not allowed!\n",__func__,mbdcnd);
                    exit(-1);

            }
            COMM(fprintf(stderr,"Proc %d line %d \n",data->this_process,__LINE__););


            // Solve radial equation
            LapSol_GaussElim(nx, &dl[0]+1,dc,du,&res[ip][0]);
        } //End for loop ip

        for(ip=ny_cut;ip<ny+2;ip++) for(ir=0;ir<nx;ir++) res[ip][ir] = 0.0;
        // if polar remove high k from region around origin
        if(ispolar)for(ir=0;ir<nx;ir++) for(ip=2*(int)(ny/8)+(int)(10.*(double)ir/(double)nx);ip<ny+2;ip++) res[ip][ir] = 0.0;

        // Backfft
        Fft_1d_2d_b(res,nx,ny);
        for(ip=0;ip<ny;ip++) memcpy((void *)f[ip],(void *)res[ip],nx*sizeof(double));
    }
    FIRST = FALSE;
}

/***************************************************************************/
/*!
  Helmholtz solver
  to solve doubly periodic problem :

   \f[
   \omega_{m,n} = - m^2  f_{m,n} - n^2 f_{m,n} + \lambda f_{m,n}
   \f]

   To solve the implizit part of the SS3  \f[ \omega^* = (1 + \lambda^{-1} \nabla^2 ) \omega  \f]

   We here solve:

   \f[ \lambda \omega^* = (\lambda + \nabla^2) \omega  \f]

   and renormalize \f[ \omega^* \f]  with \f[ \lambda \f] !

 */

void LapSol_PeriodicSolver(HDF_DS *data, PARA *p,double **f, double **w,double lamda, int normalize)

{
    register int ir,ip;
    static int nx,ny;

    static double **k2,**k2m1;
    static double **fl,**flout;

    static int FIRST=TRUE;

#ifdef Intel_MPL
    static DFTI_DESCRIPTOR_HANDLE desc_handle;
    long dfti_status;
    static long dim[2];
    int mky,nkx;
#endif

#ifdef FFTW
  static fftw_complex *fftwcout;
  static fftw_plan plan_backward;
  static fftw_plan plan_forward;
  static int nxh;
#endif

    static double kx,ky,kymin;
    double mult;


    if(FIRST)
    {

        nx = data->lnx;
        ny = data->lny;
        fl     = Util_DMatrix(ny,0,nx,0);
        flout     = Util_DMatrix(ny,0,nx,0);

        k2   = Util_DMatrix(ny+2,data->offy,nx+2,data->offx);
        k2m1 = Util_DMatrix(ny+2,data->offy,nx+2,data->offx);

        p->dkx = 2.0*M_PI/(p->xmax-p->xmin);
        p->dky = 2.0*M_PI/(p->ymax-p->ymin);
        fprintf(stderr,"dky = %f\n",p->dky);

        kymin = -(double)(ny/2)*p->dky;

        for(ip=-1;ip<=ny+2;ip++) for(ir=-1;ir<=nx+2;ir++)
            k2[ip][ir] = k2m1[ip][ir] = 0.;


#ifdef Intel_MPL
        dim[0]=(long)ny;
        dim[1]=(long)nx;

        dfti_status = DftiCreateDescriptor(&desc_handle, DFTI_DOUBLE,DFTI_REAL,2,dim);
        dfti_status = DftiSetValue(desc_handle,DFTI_PACKED_FORMAT,DFTI_PERM_FORMAT);
        dfti_status = DftiSetValue(desc_handle,DFTI_BACKWARD_SCALE,1./((double)(ny*nx)));
        dfti_status =DftiCommitDescriptor(desc_handle);

        kymin = 0.;

        /* first 4 elements */

        ip = 0;ir = 0;
        ky = 0.;kx = 0.;
        k2[ip][ir]= kx*kx+ky*ky;


        ip = 1;
        ky = kymin+(double)(ny/2)*p->dky;
        k2[ip][ir]= kx*kx+ky*ky;


        ip = 0;ir = 1;
        ky = 0.;kx = ((double)(nx/2))*p->dkx;
        k2[ip][ir] = kx*kx+ky*ky;

        ip = 1;
        ky = kymin+(double)(ny/2)*p->dky;
        k2[ip][ir] = ky*ky+kx*kx;


        /* first row */


        for (ip=2,ky=p->dky;ip<ny;ip+=2,ky += p->dky)
        {
            ir = 0; kx = 0.;
            k2[ip+1][ir] = k2[ip][ir] =kx*kx+ky*ky;

            ir = 1;  kx = ((double)(nx/2))*p->dkx;
            k2[ip+1][ir] = k2[ip][ir] =kx*kx+ky*ky;
        }


        /* rest of matrix */

        for (ip=0;ip<=ny/2;ip++)
            for (ir=2;ir<nx;ir+=2)
            {
                ky = (double)(ip)*p->dky;
                kx = (double)((ir)/2)*p->dkx;

                k2[ny-ip][ir+1]= k2[ny-ip][ir]=k2[ip][ir+1]  = k2[ip][ir] =kx*kx+ky*ky;
            }
#endif

#ifdef FFTW

        nxh = ( nx / 2 ) + 1;
        fftwcout = (fftw_complex *) fftw_malloc ( sizeof ( fftw_complex ) * ny * nxh );

        plan_forward = fftw_plan_dft_r2c_2d(ny, nx,&fl[0][0], fftwcout,FFTW_ESTIMATE);
        plan_backward =  fftw_plan_dft_c2r_2d(ny, nx, fftwcout, &flout[0][0],FFTW_ESTIMATE );


        /* Positive kx,ky */

        for (ip=0,ky=0.0;ip<ny/2;ip++,ky += p->dky)
            for (ir=0,kx=0.;ir<(nx+2)/2;ir+=2,kx += p->dkx)
                {
                    k2[ip][ir+1] = k2[ip][ir] =(kx*kx+ky*ky);
                }
        /* Positive finite ky, negative kx */

        for (ip=ny-1,ky = p->dky;ip>=ny/2;ip--,ky += p->dky)
            for (ir=0,kx=0.;ir<(nx+2)/2;ir+=2,kx += p->dkx)
                {
                    k2[ip][ir+1] = k2[ip][ir] =(kx*kx+ky*ky);
                }


        BUGREPORT;
#else //OTHER
        for (ip=0,ky=kymin;ip<ny;ip++,ky += p->dky)
                for (ir=0,kx=0.;ir<nx;ir+=2,kx += p->dkx)
                {
                    k2[ip][ir+1] = k2[ip][ir] =(kx*kx+ky*ky);
                    COMM(fprintf(stderr,"x = %d,y = %d, k2 = %f\n",ir,ip,k2[ip][ir]););
                }
        BUGREPORT;
        for (ip=0,ky=kymin;ip<ny;ip++,ky += p->dky)  k2[ip][1] = (kx*kx+ky*ky);
        k2[ny/2][0]  = 0.;

#endif

        for (ip=0;ip<ny;ip++) for (ir=0;ir<nx+2;ir++)
            if (k2[ip][ir] != 0.)
                k2m1[ip][ir] = 1./k2[ip][ir];
            else
                k2m1[ip][ir] = 0.;
        BUGREPORT;


        COMM(fprintf(stderr,"dky = %g dkx = %g\n",p->dky,p->dkx););
        FIRST = FALSE;
        }


/* initialisation complete */

    if(!normalize || lamda == 0.)
        mult = 1.;
    else
        mult = lamda;

#ifdef FFTW

       mult /= (double)(nx*ny);

#endif

    COMM(fprintf(stderr,"Proc %d %s: line %d \n",data->this_process,__FILE__,__LINE__););

        for(ip=0;ip<ny;ip++) for(ir=0;ir<nx;ir++)  fl[ip][ir] = w[ip][ir]*mult;

#ifdef Intel_MPL
        DftiComputeForward(desc_handle,&fl[0][0]);

#endif // Intel

#ifdef  FFTW

     fftw_execute ( plan_forward );

/*
         for (ip=0;ip<ny;ip++) for (ir=0;ir<nxh;ir++)
             if ((fftwcout[ip*nxh+ir][0]*fftwcout[ip*nxh+ir][0]+ fftwcout[ip*nxh+ir][1] * fftwcout[ip*nxh+ir][1] ) > 1./(double)(nx*ny))
             {
                 fprintf(stderr, "-> (%d, %d ) \n", ir,ip);
             }
*/

     if(lamda == 0.)
         for (ip=0;ip<ny;ip++) for (ir=0;ir<nxh;ir++)
         {
             fftwcout[ip*nxh+ir][0] *= -k2m1[ip][2*ir];
             fftwcout[ip*nxh+ir][1] *= -k2m1[ip][2*ir+1];
         }
     else
         for (ip=0;ip<ny;ip++) for (ir=0;ir<nxh;ir++)
         {
             fftwcout[ip*nxh+ir][0]  /= (-k2[ip][2*ir]+lamda);
             fftwcout[ip*nxh+ir][1]  /= (-k2[ip][2*ir+1]+lamda);
         }

       for (ip=0;ip<ny;ip++)  fftwcout[ip*nxh][1] = 0.;
       /*   for (ip=0;ip<ny;ip++)  fftwcout[ip*nxh][0] *= 2.;*/

 /*for (ip=0;ip<ny;ip++)  fftwcout[ip*nxh][1] *= 2.;*/
/*       for (ip=0;ip<ny;ip++)  fftwcout[ip*nxh][0] = 0.;*/

        fftw_execute ( plan_backward );
        for(ip=0;ip<ny;ip++) for(ir=0;ir<nx;ir++) f[ip][ir] = flout[ip][ir];

#else // FFTW
        FHIN(fl,nx,ny);
#endif // OTHER

        BUGREPORT;

#ifndef FFTW
        if(lamda == 0.)
            for (ip=0;ip<ny;ip++) for (ir=0;ir<nx;ir++) fl[ip][ir] *= -k2m1[ip][ir];
        else
            for (ip=0;ip<ny;ip++) for (ir=0;ir<nx;ir++) fl[ip][ir] /= (-k2[ip][ir]+lamda);

#ifdef Intel_MPL
        DftiComputeBackward(desc_handle,&fl[0][0]);
 for(ip=0;ip<ny;ip++) for(ir=0;ir<nx;ir++) f[ip][ir] = fl[ip][ir];
#else
        FBACK(fl,nx,ny);
 for(ip=0;ip<ny;ip++) for(ir=0;ir<nx;ir++) f[ip][ir] = fl[ip][ir];
#endif

#endif //!FFTW





        COMM(fprintf(stderr,"Finished fft job...\n"););
        return;
}




/******************************************************************************************/
/*!
   Purpose
   =======

   LAPSOL_GAUSSELIM  solves the equation

      A*X = B,

   where A is an N-by-N tridiagonal matrix, by Gaussian elimination with
   partial pivoting.

   Note that the equation  A'*X = B  may be solved by interchanging the
   order of the arguments DU and DL.

   Arguments
   =========

   n       The order of the matrix A.  n >= 0.


   dl      dimension (N-1)
           On entry, DL must contain the (n-1) subdiagonal elements of
           A.
           On exit, DL is overwritten by the (n-2) elements of the
           second superdiagonal of the upper triangular matrix U from
           the LU factorization of A, in DL(1), ..., DL(n-2).

   d       dimension (N)
           On entry, D must contain the diagonal elements of A.
           On exit, D is overwritten by the n diagonal elements of U.

   du      dimension (N-1)
           On entry, DU must contain the (n-1) superdiagonal elements
           of A.
           On exit, DU is overwritten by the (n-1) elements of the first
           superdiagonal of U.

   b       dimension (N*STRIDE)
           On entry, the right hand side vector B.
           On exit, if function returns 0 , the solution  X.


   returns
           = 0:  successful exit
           < 0:  if INFO = -i, the i-th argument had an illegal value
           > 0:  if INFO = i, U(i,i) is exactly zero, and the solution
                 has not been computed.  The factorization has not been
                 completed unless i = N.

   =====================================================================
 */

int  LapSol_GaussElim(int n, double *dl, double *d, double *du, double *b)
{

 #ifdef _OPENMP
   static  int
    k,kp1;

  static double
    mult,temp;

#pragma omp threadprivate (k,kp1,mult,temp)
#else
  register int
    k,kp1;

  double
    mult,temp;
 #endif


  /*
  if( n <  0 ) return -1;
  if( n == 0 ) return -2;
  */

  for(k=0,kp1=1;k<n-1;k++,kp1++)
    {
      if(dl[k] == 0.0 )
        {
          /*     Subdiagonal is zero, no elimination is required.*/
          if( d[k] == 0.0)   return k;          /* Diagonal zero:return k; a unique
                                                   solution can not be found.*/
        }
/*       else if( fabs(d[k]) >= fabs(dl[k]) ) */
      else if( d[k]*d[k] >= dl[k]*dl[k] )
        {
          /*     No row interchange required */

          mult = dl[k] / d[k];

          d[kp1] -=  mult*du[k];
          b[kp1] -=  mult* b[k];

          if( k <  (n-2) ) dl[k] = 0.0;
        }
      else
        {
          /*Interchange rows K and K+1 */

          mult   = d[k] / dl[k];
          d[k]   =  dl[k];
          temp   = d[kp1];
          d[kp1] = du[k] - mult*temp;

          if(k < (n-2))
            {
              dl[k]   = du[kp1];
              du[kp1] = -mult*dl[k];
            }

          du[k] = temp;
          temp = b[k];
          b[k] = b[kp1];
          b[kp1] = temp - mult*b[kp1];
        }
    }
  if(d[n-1] == 0.0) return n;

  /*     Back solve with the matrix U from the factorization.*/

  b[n-1] = b[n-1]/d[n-1];

  if( n > 0 )
    b[n-2] = ( b[n-2]-du[n-2]*b[n-1] ) / d[n-2];

  for(k = (n-3);k >= 0;k--)
    {
     b[k]  = ( b[k] - du[k] * b[k+1] - dl[k]*b[k+2]) /d[k];
    }
  return 0;
}


/**************************************************************************/

/*!
 *  Purpose
 *  =======
 *
 *  COMPLEX LAPSOL_GAUSSELIM  solves the equation
 *
 *     A*X = B,
 *
 *  where A is an N-by-N tridiagonal matrix, by Gaussian elimination with
 *  partial pivoting.
 *
 *  Note that the equation  A'*X = B  may be solved by interchanging the
 *  order of the arguments DU and DL.
 *
 *  Arguments
 *  =========
 *
 *  n       The order of the matrix A.  n >= 0.
 *
 *
 *  dl      dimension 2*(N-1)
 *          On entry, DL must contain the (n-1) complex subdiagonal elements of
 *          A.
 *          On exit, DL is overwritten by the (n-2) elements of the
 *          second superdiagonal of the upper triangular matrix U from
 *          the LU factorization of A, in DL(1), ..., DL(n-2).
 *
 *  dr      dimension 2*(N)
 *          On entry, D must contain the complex diagonal elements of A.
 *          On exit, D is overwritten by the n diagonal elements of U.
 *
 *  du      dimension 2*(N-1)
 *          On entry, DU must contain the (n-1) superdiagonal elements
 *          of A.
 *          On exit, DU is overwritten by the (n-1) elements of the first
 *          superdiagonal of U.
 *
 *  b       dimension 2*(N)
 *          On entry, the right hand side vector B.
 *          On exit, if function returns 0 , the solution  X.
 *
 *
 *  Complex numbers are store real,imag in a double c-vector
 *
 *
 *  returns
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *          > 0:  if INFO = i, U(i,i) is exactly zero, and the solution
 *                has not been computed.  The factorization has not been
 *                completed unless i = N.
 *
 *  =====================================================================
 */
int  LapSol_GaussElim_Complex(int n, double *dl, double *d, double *du,double *b)
{



  register int
    kr,ki,krp1,kip1,klr,kli,klrm1,klim1;
  double
    mult_i,mult_r,
    tempi,tempr,tr2,ti2,
    betrag,betragl,betragc;


  if( n <  0 ) return -1;
  if( n == 0 ) return -2;




  for(kr=0,ki=1,krp1=2,kip1=3;kr<n-1;kr+=2,ki+=2,krp1+=2,kip1+=2)
  {
      betragl = dl[kr]*dl[kr]+dl[ki]*dl[ki];
      betragc = d[kr]*d[kr] + d[ki]*d[ki];
      if(betragl == 0.)
      {
          /*     Subdiagonal is zero, no elimination is required.*/
          if( betragc == 0. )
              return kr/2;                      /* Diagonal zero:return k; a unique
                                               solution can not be found.*/
      }
      else if(betragc  >= betragl )
      {
          /*     No row interchange required */
          betrag = 1./betragc;

          mult_r = betrag*(dl[kr]*d[kr]+dl[ki]*d[ki]);
          mult_i = betrag*(dl[ki]*d[kr]-dl[kr]*d[ki]);

          d[krp1] -=  mult_r*du[kr]-mult_i*du[ki];
          d[kip1] -=  mult_r*du[ki]+mult_i*du[kr];

          b[krp1] -=  mult_r*b[kr]-mult_i*b[ki];
          b[kip1] -=  mult_r*b[ki]+mult_i*b[kr];


          if( kr <  2*(n-2) ) {dl[kr] = 0.0; dl[ki] = 0.;}
        }
      else
        {
          /*Interchange rows K and K+1 */

          betrag = 1./betragl;

          mult_r = betrag*(d[kr]*dl[kr]+d[ki]*dl[ki]);
          mult_i = betrag*(d[ki]*dl[kr]-d[kr]*dl[ki]);

          d[kr]   =  dl[kr];      d[ki]   =  dl[ki];

          tempr   = d[krp1];      tempi   = d[kip1];


          d[krp1] = du[kr] - (mult_r*tempr-mult_i*tempi);
          d[kip1] = du[ki] - (mult_r*tempi+mult_i*tempr);

          if(kr < 2*n-4)
            {
              dl[kr]   = du[krp1];
              dl[ki]   = du[kip1];

              du[krp1] = -(mult_r*dl[kr]-mult_i*dl[ki]);
              du[kip1] = -(mult_r*dl[ki]+mult_i*dl[kr]);
            }



          du[kr]  = tempr;   du[ki] = tempi;

          tempr = b[kr];     tempi = b[ki];

          b[kr] = b[krp1];   b[ki] = b[kip1];


          tr2 = b[krp1];     ti2 = b[kip1];

          b[krp1] = tempr - (mult_r*tr2-mult_i*ti2);
          b[kip1] = tempi - (mult_i*tr2+mult_r*ti2);
        }
    }

  klr = 2*n-2;
  kli = 2*n-1;
  klrm1 = 2*n-4;
  klim1 = 2*n-3;

  if( d[klr]*d[klr]+d[kli]*d[kli] == 0.0) return n;


  /*     Back solve with the matrix U from the factorization.*/

  betrag = 1./(d[klr]*d[klr]+d[kli]*d[kli]);


  tr2 = b[klr];     ti2 = b[kli];

  b[klr] = betrag*(tr2*d[klr]+ti2*d[kli]);
  b[kli] = betrag*(ti2*d[klr]-tr2*d[kli]);


  if( n > 0 )
    {
    tempr = b[klrm1] - (du[klrm1]*b[klr] - du[klim1]*b[kli]);
    tempi = b[klim1] - (du[klrm1]*b[kli] + du[klim1]*b[klr]);

    betrag = 1./(d[klrm1]*d[klrm1]+d[klim1]*d[klim1]);



    b[klrm1]  =betrag*(tempr*d[klrm1]+ tempi*d[klim1]);
    b[klim1]  =betrag*(tempi*d[klrm1]- tempr*d[klim1]);

    }


  for(kr = 2*n-6,ki=2*n-5;kr >= 0;kr-=2,ki-=2)
    {
      tempr  = b[kr] - (du[kr]*b[kr+2]-du[ki]*b[ki+2]) - (dl[kr]*b[kr+4] - dl[ki]*b[ki+4]) ;
      tempi  = b[ki] - (du[kr]*b[ki+2]+du[ki]*b[kr+2]) - (dl[kr]*b[ki+4] + dl[ki]*b[kr+4]) ;


      betrag = 1./(d[kr]*d[kr]+d[ki]*d[ki]);


      b[kr]  = betrag*(tempr*d[kr]+tempi*d[ki]);
      b[ki]  = betrag*(tempi*d[kr]-tempr*d[ki]);
    }
  return 0;
}

/***********************************************************************/
/*!
  Helmholtz solver for a staggered radial grid.
   - Fouriertransform in poloidal
   - Gaussian elimination in radial direction
   - solves 2nd order finite difference approximation for each poloidal mode m:

   \f[
   \omega_m =  - m^2  g_{ff}  f_m + \partial_{rr} f_m + \lambda f_m
   \f]

   To solve the implizit part of the SS3  w* = (1 - lambda^-1 lap ) w we consider:

   lamda w* = (lamda -lap) w

   and renormalize w* with lamda, that is scale has to be set to true!


   - Boundary conditions

   MBDCND:      R = A           R = B
     0                PERIODIC
     1          dir              dir
     2          dir              neu
     3          neu              neu
     4          neu              dir
     5          R = 0            dir
     6          R = 0            neu

   -Function gets on input:

   data and para structures
   w:                           left side
   mdcnd:                      radial coordinate info
   ibdra, ibdrb:                double fields wich contain boundary values
   lamda:                       parameter lambda of helmholtz equation

   scale:                         boolean to indicate if result should be normalized with lamda (TRUE)
   res:                            workspace of same dimension as w_0, fw

   on output:
   f                 result,

*/
/***********************************************************************/
void Laplace_Solve2D_Par(HDF_DS *d, PARA *p,double **f, double **w,
                         double *ibdra,double *ibdrb,double *hval,
                         double *rcor,int mbdcnd,
                         double lamda,int scale)

{
    static int ir,ip,in,i,ipl;
    static int nx,ny,ny_cut,g_nx;
    static int yblock;

    static double *dl,*dc,*du;
    static double *dl_temp,*dc_temp,*du_temp;
    static double *bdra,*bdrb;
    static double *r2,*ky2;
    static double **res;
    static double **fft_bdr;
    static double *gather_res, *g_hval,*g_rcor;
    static double **g_res;

    static int FIRST=TRUE, ispolar=FALSE;
    static double ky,kymin,r,r2d,mult,ehedx;

    BUGREPORT;
    if((strncmp(d->coordsys,"pol",3) == 0) || strncmp(d->coordsys,"cyl",3) == 0)
        ispolar = TRUE;

    //  doubly periodic solver not implemented in parallel
    if(mbdcnd == 0)
    {
        fprintf(stderr,"Periodic radial domain and parallelisation in nx incompatible !!!!!\n");
        MPI_Finalize();
        exit(-1);
    }

    if(FIRST)
      {
        nx = d->lnx;
        ny = d->lny;
        ny_cut= 5*(int)(ny/6.);
        ny_cut=   ny+2;

                g_nx= d->dims[2];

                yblock=ny/d->N[2];
        // check if yblock is hole number
        if(ny != yblock*d->N[2])
        {
            fprintf(stderr,"Incompatible Number %d of processes. NY needs to be multiple of procs in radial direction %d !\n",
                    (int)d->num_procs,(int)d->N[2]);
            MPI_Finalize();
            exit(-1);
        }


        bdra  = Util_DVector(ny+2,1);
        bdrb  = Util_DVector(ny+2,1);
        ky2   = Util_DVector(ny+2,1);


        dl    = Util_DVector(g_nx,2);
        dc    = Util_DVector(g_nx,2);
        du    = Util_DVector(g_nx,2);

        dl_temp    = Util_DVector(g_nx,2);
        dc_temp    = Util_DVector(g_nx,2);
        du_temp    = Util_DVector(g_nx,2);

        r2         = Util_DVector(g_nx,2);
                g_hval     = Util_DVector(g_nx,2);
                g_rcor     = Util_DVector(g_nx,2);

        res        = Util_DMatrix(ny+2,0,nx,0);
        fft_bdr    = Util_DMatrix(ny+2,0,2,0);

                //gather_res is used as a buffer when res matrix is received in one MPI_Gather command
                gather_res = Util_DVector(yblock*g_nx,0);

                //g_res contains part of the res matrix from all processors in x-direction
                g_res = Util_DMatrix(yblock,0,g_nx,0);

        p->dky =  2.0*M_PI/(p->ymax-p->ymin);
        for(ip=0; ip < (ny+2);  ip++) ky2[ip] = (double)((ip)/2)*p->dky;
        for(ip=0; ip < (ny+2);  ip++) ky2[ip]*= ky2[ip];

        MPI_Allgather(&hval[0],nx,MPI_DOUBLE,&g_hval[0],nx,MPI_DOUBLE,d->xrow_comm);
        MPI_Allgather(&rcor[0],nx,MPI_DOUBLE,&g_rcor[0],nx,MPI_DOUBLE,d->xrow_comm);

        // generate temporary arrays needed in main loop

        if(ispolar)
            {
                COMM(fprintf(stderr,"Setting up for POLAR coordinate system in PAR Laplace solve.\n"););

                for(ir=0;ir<g_nx;ir++)
                {
                    r     = g_rcor[ir];
                    r2[ir]= r*r;

                    //fprintf(stderr,"%d %f %f.\n",ir,r, r2[ir]);
                    r2d   = r*r*g_hval[ir]*g_hval[ir];
                    ehedx = g_hval[ir]*r*0.5;

                    dl_temp[ir] = -ehedx +    r2d;
                    dc_temp[ir] = - 2.*r2d;
                    du_temp[ir] =  ehedx +    r2d;
                }
            }
            else
            {
                COMM(fprintf(stderr,"Setting up for RECTANGULAR coordinate system in PAR Laplace solve.\n"); );

                for(ir=0;ir<g_nx;ir++)
                {
                        r2d   = g_hval[ir]*g_hval[ir];
                        r2[ir]= 1.;
                        dl_temp[ir] =   r2d;
                        dc_temp[ir] = - 2.*r2d;
                        du_temp[ir] =   r2d;
                }
            }

            BUGREPORT;
            COMM(fprintf(stderr,"Proc %d %s line FIRST:  %d: nx  %d ny %d\n",d->this_process,__FILE__,__LINE__,nx,ny););
      }

    if(!scale || lamda == 0.)
        mult = 1.;
    else
        mult = lamda;

    BUGREPORT;
    for(ip=0;ip<ny;ip++) for(ir=0;ir<nx;ir++) res[ip][ir] = w[ip][ir]*r2[nx*d->grid_coords[2]+ir]*mult;

    //for(ir=0;ir<nx;ir++) fprintf(stderr,"%d %d %f\n",d->grid_coords[2],ir,r2[nx*d->grid_coords[2]+ir]);



    /* Something to remember:

    m = 0 in f[0],f[1]
    m = 1 in f[2],f[3]
    m = 2 in f[4],f[5] ...etc.
    */

    // FFT boundary AND field
    for(ip=0;ip<ny;ip++) fft_bdr[ip][0]=ibdra[ip];
    for(ip=0;ip<ny;ip++) fft_bdr[ip][1]=ibdrb[ip];

    // FFT boundary and field separately, otherwise communication in allgather fails.
    BUGREPORT;
    Fft_1d_2d_f(res,nx,ny);
    Fft_1d_2d_f(fft_bdr,2,ny);

    for(ip=0;ip<ny+2;ip++) bdra[ip] = fft_bdr[ip][0];
    for(ip=0;ip<ny+2;ip++) bdrb[ip] = fft_bdr[ip][1];

    // Set last two rows to zero
    for(ir=0;ir<nx;ir++) res[ny][ir] = res[ny+1][ir] = 0.;


    // Transpose data in sending
    BUGREPORT;
        for (i=0;i<d->N[2];i++) {
                MPI_Gather(&res[yblock*i][0],nx*yblock,MPI_DOUBLE,&gather_res[0],nx*yblock,MPI_DOUBLE,i,d->xrow_comm);
        }

    // resort to array
        for(ip=0;ip<yblock;ip++) for(i=0;i<d->N[2];i++) for(ir=0;ir<nx;ir++) {
                g_res[ip][ir+nx*i]=gather_res[ir+nx*yblock*i+nx*ip];
        }
    BUGREPORT;



    for(ipl= 0, ip=d->xrow_id*yblock;ip<(d->xrow_id+1)*yblock;ip++,ipl++)
      {

        /*setup vectors of tridi matrix and multiply rhs by r2 */
        for(ir=0;ir<g_nx;ir++)
        {
            dl[ir] =       dl_temp[ir];
            dc[ir] = -ky2[ip] + dc_temp[ir] +  r2[ir] * lamda;
            du[ir] =       du_temp[ir];
        }

        // Enter boundary values

        switch (mbdcnd)
          {
          case 0:
            fprintf(stderr,"Boundary %d not allowed with poloidal geometry!\n",mbdcnd);
            exit (-1);
            break;
          case 1:case 2:case 11:
            dc[0]         -=  dl[0];
            g_res[ipl][0]    -=  2.*dl[0]*bdra[ip];
            break;
          case 3:case 4:
            dc[0]         +=  dl[0];
            g_res[ipl][0]    +=  dl[0]/g_hval[0]*bdra[ip];
            break;
          case 5:case 6:
            break;
          default:
              fprintf(stderr,"%s: Boundary %d not allowed!\n",__func__,mbdcnd);
            exit(-1);
          }

        ir = g_nx-1;

        switch (mbdcnd)
          {
          case 0:
            if(ispolar )
              {
                  fprintf(stderr,"%s: Boundary %d not allowed with poloidal geometry!\n",__func__,mbdcnd);
                exit (-1);
              }
            break;
          case 1:case 4:case 5:
            dc[ir]      -=  du[ir];
            g_res[ipl][ir] -=  2.*du[ir]*bdrb[ip];
            break;
          case 2:case 3:case 6:
            dc[ir]      +=  du[ir];
            g_res[ipl][ir] -=  du[ir]/g_hval[ir]*bdrb[ip];
            break;
           case 11: //dirchlet for n != 0, neuman for n == 0
            dc[ir]      -=  du[ir];
            g_res[ipl][ir] -=  2.*du[ir]*bdrb[ip];
            if (ipl == 0)
              {
                dc[ir]      +=  du[ir];
                g_res[ipl][ir] -=  0.;
              }
            break;
          default:
            fprintf(stderr,"Subroutine nauplr: Boundary %d not allowed!\n",mbdcnd);
            exit(-1);
          }

        // Solve radial equation
        BUGREPORT;
        LapSol_GaussElim(g_nx, &dl[0]+1,dc,du,&g_res[ipl][0]);

      } //End loop ip

    BUGREPORT;
        //Transform data to vector for communication
        for(ip=0;ip<yblock;ip++) for(in=0;in<d->N[2];in++) for(ir=0;ir<nx;ir++) {
                gather_res[ir+nx*yblock*in+nx*ip]=g_res[ip][ir+nx*in];
        }

        //Scatter data back to each processor
        for (in=0;in<d->N[2];in++) {
                MPI_Scatter(&gather_res[0],nx*yblock,MPI_DOUBLE,&res[in*yblock][0],nx*yblock,MPI_DOUBLE,in,d->xrow_comm);
        }

    BUGREPORT;
    for(ip=ny_cut;ip<ny+2;ip++)  for(ir=0;ir<nx;ir++) res[ip][ir] = 0.0;
    // if polar remove high k from region around origin
    if(ispolar) for(ir=0;ir<nx;ir++) for(ip=2*(int)(ny/8)+(int)(10.*(double)(nx*d->grid_coords[2]+ir)/(double)g_nx);ip<ny+2;ip++) res[ip][ir] = 0.0;

    // Backfft
    Fft_1d_2d_b(res,nx,ny);
    for(ip=0;ip<ny;ip++) memcpy((void *)f[ip],(void *)res[ip],nx*sizeof(double));
    BUGREPORT;
    FIRST = FALSE;

}
