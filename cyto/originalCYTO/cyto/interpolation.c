/* This File contains some subroutines concerned with spectral approximations */

#include "utilities.h"
#include "interpolation.h"



/* Calculate velocity field at nop particle postions */


void ipol_vel_p_pos(double *u,double *v,double *xpos,double *ypos,int nop,double **phi,
					     int nx, int ny, double *Sx, double *Sy)
/* Subroutine gets as input:

   u,v:  empty fields of di

mension np to contain velocity values 

   xpos: double field with x positions of particles
   ypos: double field with y positions of particles
   nop:  number of particles
   phi:  double file containing electrostatic potential 
   nx,ny: Dimension of potential field
   Sx,Sy: Double filed containing k-values


on output:
   u,v: velocity at positions x,y


Restriction: nx < 2048
   */
{
int ix,iy,k;
double tcsx[2048];
double cy,sy;
double ukr,uki,vkr,vki;
double val,norm;    



/* u = u_x = -d_y phi 
   v = u_y =  d_x phi  */




#ifdef aix
norm = 1./((double)(nx*ny));
#else
norm = 1./(sqrt((double)ny));
#endif

for(k=0;k<nop;k++)
  {
    /* Precalculate kx values */
    for(ix=0;ix<nx;ix+=2)
      {
	val = -xpos[k]*Sx[ix];
	tcsx[ix]  = cos(val);
	tcsx[ix+1]= sin(val);
      }

    u[k] = 0.;
    v[k] = 0.;

    for(iy=0;iy<ny;iy++)
      {
	val = ypos[k]*Sy[iy];
	cy = cos(val);
	sy = sin(val);
	
	
      	ukr = -Sy[iy]*phi[iy][1];
	uki =  Sy[iy]*phi[iy][0];

	vkr = -Sx[0]*phi[iy][1];
	vki =  Sx[0]*phi[iy][0];	

	u[k]+= (tcsx[0]*cy-tcsx[1]*sy)*ukr
	  -(tcsx[1]*cy+tcsx[0]*sy)*uki;

	v[k]+= (tcsx[0]*cy-tcsx[1]*sy)*vkr
	          -(tcsx[1]*cy+tcsx[0]*sy)*vki;
		  

	for(ix=2;ix<nx;ix+=2)
	  {

	    ukr = -Sy[iy]*phi[iy][ix+1];
	    uki =  Sy[iy]*phi[iy][ix];

	    vkr = -Sx[ix]*phi[iy][ix+1];
	    vki =  Sx[ix]*phi[iy][ix];


	    u[k]+= 2.*((tcsx[ix  ]*cy-tcsx[ix+1]*sy)*ukr
	          -(tcsx[ix+1]*cy+tcsx[ix  ]*sy)*uki);

	    v[k]+= 2.*((tcsx[ix  ]*cy-tcsx[ix+1]*sy)*vkr
	          -(tcsx[ix+1]*cy+tcsx[ix  ]*sy)*vki);

	  }

      }
  }



for(k=0;k<nop;k++) u[k]*=-norm;
for(k=0;k<nop;k++) v[k]*=-norm;


#undef DEBUG
#ifdef DEBUG
for(k=0;k<nop;k++)
  fprintf(stderr,"%d\t%f\t%f\n",k,u[k],v[k]);

exit(0);
#endif

}


/***********************************************************************/

/* Interpolate value of field to position */

void spectral_interpol(double *ret_val,double *xpos,double *ypos,int nop,double **phi,int nx, int ny, double *Sx, double *Sy)
/* Subroutine gets as input:
   xpos: x positions
   ypos: y positions
   nop : number of positions
   phi:  fourier transformed field 
   nx,ny: Dimension of potential field
   Sx,Sy: Double filed containing k-values


on return:
   val value of field at positions x,y


Restriction: nx < 2048
   */

{
int i,ix,iy;
double tcsx[2048];
double cy,sy,cx,sx;
double val,norm=1.;    

#ifdef aix
norm = 1./((double)(nx*ny));
#else
norm = 1./sqrt((double)ny);
#endif


for(i=0;i<nop;i++)
  {
    /* printf("%f %f\n",xpos[i],ypos[i]);*/
    /* Precalculate kx values */
    for(ix=0;ix<nx;ix+=2)
      {
	val = -xpos[i]*Sx[ix];
	tcsx[ix]  = cos(val);
	tcsx[ix+1]= sin(val);
      }

    ret_val[i] = 0.;
    
    for(iy=0;iy<ny;iy++)
      {
	val = ypos[i]*Sy[iy];
	cy = cos(val);
	sy = sin(val);

	cx = tcsx[0];
	sx = tcsx[1];

	
       	ret_val[i]+= (cx*cy-sx*sy)*phi[iy][0] - (sx*cy+cx*sy)*phi[iy][1];
	
	for(ix=2;ix<nx;ix+=2)
	  {
	    cx = tcsx[ix];
	    sx = tcsx[ix+1];

	    ret_val[i]+= 2.*((cx*cy-sx*sy)*phi[iy][ix]-(sx*cy+cx*sy)*phi[iy][ix+1]); 

	    /*
	      if((phi[iy][ix] != 0.) || (phi[iy][ix+1] != 0.))
	      printf("%d %d (%f,%f)\n",iy,ix,phi[iy][ix],phi[iy][ix+1]);*/
	  }
	
      }
    ret_val[i]*=norm;
 
  }

}






