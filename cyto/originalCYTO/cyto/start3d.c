/*	HDF Startprogramm fuer 3D Verteilungen
	von
	V. Naulin initial revision 25.08.1997

*/
#undef DEBUG

#include "utilities.h"
#include "interpolation.h"
#ifndef RAND_MAX
#define RAND_MAX    32767
#endif

#define FHIN	powfour2d
#define FBACK 	powfour2dre
#define FHIND   oddft2d
#define FBACKD 	oddft2dre

#define USAGE	{fprintf(stderr,\
"This program produces initial 3D Util_Distributions.\n"\
"Parameters shall be listed in a file named Name.ini.\n"\
"\nusage : %s Name to produce a HDF file from Name.ini file\n"\
"\nusage : %s HDFfile -R to produce an .ini file from and existing HDF file \n"\
"\n or %s -help for help on ini-file.\n",argv[0],argv[0],argv[0]);exit(-1);}

void lamb(HDF_DS *, PARA *p);
void drift_Util_Noise(HDF_DS *data, PARA *p);
void random1(HDF_DS *data, PARA *p);
void map_Util_Noise(HDF_DS *data, PARA *p);
void peritest(HDF_DS *data, PARA *p);
void bdadjust(HDF_DS *, PARA *p);
void scott1(HDF_DS *data, PARA *p);
void scott2(HDF_DS *data, PARA *p);
void dent(HDF_DS *data, PARA *p);

int verteilung;

int main(int argc,char **argv)
{
HDF_DS  ds;
PARA   para;
static char rcsid[]="$Id: start3d.c,v 4.1 2002/06/26 08:27:03 vona Exp $";
char ininame[256],barename[256];
double val;

int i,j,k,l;

if(argv[1] == NULL) USAGE

if(! strcmp(argv[argc-1],"-help"))
	{
	fprintf(stderr,\
" The ini-file has the following syntax:\n\n"\
" The ini-file is a pure ascii file. "\
" Comment Lines start with #.\n"\
" In every line there may be one entry of the form:\n"\
"  [variable name] = [numerical value]\n"\
" there may be more than one occurrence of [variable name], \n"\
" but only the first one is evaluated.\n\n"\
" Giving an empty ini-file produces a template .ini file to stdout.\n"\
" Giving a -R filename as argument, reproduces an .ini file from a data-file\n\n");exit(0);
 	}



sprintf(ininame,"%s.ini",argv[1]);
sprintf(barename,"%s",argv[1]);

 para.desc = "empty";
FUtils_IniStructure(&ds,&para,argc,argv,rcsid);

if(! strcmp(argv[1],"-R") )
  {
    if(argv[2] == NULL) USAGE
    FUtils_ReadNetCDF(argv[2],-1,&ds,&para);   
    Util_ReadIniFile(&ds,&para,ininame,1);
    exit(0);
  } 


if((Util_ReadIniFile(&ds,&para,ininame,0)) == -1) 
  if ((Util_ReadIniFile(&ds,&para,barename,0)) == -1) exit(0);


ds.create = TRUE;

/* Erzeugen der Verteilung						 */

 if(ds.nx*ds.ny*ds.nz < 1) 
   {
     fprintf(stderr,"Number of points should be larger zero. Dims are (%d,%d,%d).\n Exiting\n",
	     ds.nx,ds.ny,ds.nz);
     exit(-1);
   }

ds.dddfelder[0] = Util_DCube(ds.nx,0,ds.ny,0,ds.nz,0);
ds.dddfelder[1] = Util_DCube(ds.nx,0,ds.ny,0,ds.nz,0);
ds.dddfelder[2] = Util_DCube(ds.nx,0,ds.ny,0,ds.nz,0);

/* Erzeugen des Koordinatensystems */

for(i=0,val = para.zmin+para.dz*.5;i<ds.nz;i++,val+=para.dz)  ds.coordinate[0][i] = val;
for(i=0,val = para.ymin+para.dy*.5;i<ds.ny;i++,val+=para.dy)  ds.coordinate[1][i] = val;
for(i=0,val = para.xmin+para.dx*.5;i<ds.nx;i++,val+=para.dx)  ds.coordinate[2][i] = val;

sprintf(ds.dim_label[0],"z");
sprintf(ds.dim_label[1],"phi");
sprintf(ds.dim_label[2],"r");

ds.anzahl = 3;


sprintf(ds.desc,"3d test verteilung.") ;

sprintf(ds.revision,"%s",rcsid);
sprintf(ds.integrator,"%s",__FILE__);
sprintf(ds.compiliert,"%s %s",__DATE__,__TIME__);
sprintf(ds.maschine,"%s","Start-Datei!");
sprintf(ds.jobid,"%s","Start-Datei!");
sprintf(ds.names[0],"%s","Potential");
sprintf(ds.names[1],"%s","Density");
sprintf(ds.names[2],"%s","Temperature");
 





switch(para.coordsys){
	case CARTESIAN:	
		strcpy(ds.coordsys,"cartesian, equidistant");
		strcpy(ds.dim_label[0],"x\0");
		strcpy(ds.dim_label[1],"y\0");
		strcpy(ds.dim_label[2],"z\0");
		break;
	case POLOIDAL:	
		strcpy(ds.coordsys,"polar, equidistant\0");
 		strcpy(ds.dim_label[0],"z\0");
		strcpy(ds.dim_label[1],"r\0");
		strcpy(ds.dim_label[1],"phi\0");
		break;
	case CYLINDRICAL:	
		strcpy(ds.coordsys,"cylindrical, equidistant\0");
		strcpy(ds.dim_label[0],"z\0");
		strcpy(ds.dim_label[1],"r\0");
		strcpy(ds.dim_label[1],"phi\0");
 		break;
	case SPHERICAL:	
		strcpy(ds.coordsys,"spherical, equidistant\0");
		strcpy(ds.dim_label[0],"theta\0");
		strcpy(ds.dim_label[1],"r\0");
		strcpy(ds.dim_label[1],"phi\0");
                break;
	case CART_TB_X:	
		strcpy(ds.coordsys,"cartesian, x-Tschebischeff\0");
		strcpy(ds.dim_label[0],"y\0");
		strcpy(ds.dim_label[1],"x\0");
		break;
	case POL_TB_R:	
		strcpy(ds.coordsys,"polar, r-Tschebischeff\0");
 		strcpy(ds.dim_label[0],"phi\0");
		strcpy(ds.dim_label[1],"r\0");
		break;
	case CYL_TB_R:	
		strcpy(ds.coordsys,"cylindrical, r-Tschebischeff\0");
		strcpy(ds.dim_label[0],"z\0");
		strcpy(ds.dim_label[1],"r\0");
		strcpy(ds.dim_label[1],"phi\0");
 		break;
	case SPHER_TB_R:	
		strcpy(ds.coordsys,"spherical, r-Tschebischeff\0");
		strcpy(ds.dim_label[0],"theta\0");
		strcpy(ds.dim_label[1],"r\0");
		strcpy(ds.dim_label[1],"phi\0");
                break;
 	}















ds.isfloat = -1;








switch(para.verteilung)
  {
  case 2:	
    random1(&ds, &para);
    sprintf(ds.desc,"Gaussian Noise.");
    break;
  case 7:	
    lamb(&ds, &para);
    sprintf(ds.desc,"Lamb Dipol.");
    break;	
  case 8:	
    peritest(&ds, &para);
    sprintf(ds.desc,"Test");
    break;	
  case 11:	
    drift_Util_Noise(&ds, &para);
    sprintf(ds.desc,"Drift Noise");
    break;	
  case 12:	
    map_Util_Noise(&ds, &para);
    sprintf(ds.desc,"Mapped Noise");
    break;	
  case 13:	
    scott1(&ds, &para);
    sprintf(ds.desc,"Alf1");
    break;
  case 14:	
    scott2(&ds, &para);
    sprintf(ds.desc,"Alf2");
    break;
  case 15:	
    dent(&ds, &para);
    sprintf(ds.desc,"dent");
    break;
  }




/* Adjust for boundary values */

for(l=0;l<3;l++) 
  for(i=0;i<ds.nz;i++)
    for (j=0;j<ds.ny;j++)
      for (k=0;k<ds.nx;k++)
	ds.dddfelder[l][i][j][k] += para.bdval[l][0] + (para.bdval[l][1]-para.bdval[l][0])*(ds.coordinate[2][k]-para.xmin)/(para.xmax-para.xmin);


if(para.coordsys != CARTESIAN)   bdadjust(&ds, &para);

printf(" Here we are, starting to write...\n");
ds.isfloat = -1;


if(para.nu != 0.)
  para.prandel = para.mue_w/para.nu;
else
  para.prandel =0.;
ds.rank = 3;
para.time=0.0;
para.energy = 0.0;
para.vorticity = 0.0;
para.nusselt  = 0.0;
para.reynolds  = 0.0;

for(i=0;i<ds.dims[2];i++)  
  para.shat[i] = para.qprof[i]  = 0.;

for(i=0;i<ds.rank;i++)  
  ds.elements[i] = ds.dims[i];

/* Remember to set ALL string constants.....otherwise no output!!*/
for(i=0;i<3;i++) 
  {
  FUtils_WriteNetCDF(argv[1], 0 ,&ds ,&para,ds.names[i],ds.dddfelder[i]);
  ds.create = FALSE;
  }
}


/***************************************************************************/

void lamb(HDF_DS *data, PARA *p)
{
int ir,ip,iz,i;
double r;
double gamma11 = 3.83170597020751;
double bes0,bes1,theta0,r0;
double amplitude,amp0;
double val;

bes0 = j0(gamma11);
amp0 = p->uvortex *2.0 * gamma11 /(p->radius*bes0);
theta0 = p->ymin+ 0.5*(p->ymax-p->ymin);
r0 =p->xmin+ 0.5*(p->xmax -p->xmin);

for(iz=0;iz<data->nz;iz++)
  {
  amplitude =amp0*sin(((double)iz+.5)/(double)data->nz*M_PI);
  /*printf("Amplitude = %g iz = %d\n",amplitude,iz);*/
  for(ip=0;ip<data->ny;ip++)
    for(ir=0;ir<data->nx;ir++)
      {
	r = sqrt( data->coordinate[2][ir]*data->coordinate[2][ir] 
		  + r0*r0-2.0*r0*data->coordinate[2][ir]
		  *cos(data->coordinate[1][ip]-theta0) );
	/*printf("ir %d \t |r| = %g\t r = %g  \n",ir,r,data->coordinate[2][ir] );*/
	if ((r > p->radius) || (r == 0.))
	  for(i=0;i<3;i++)  data->dddfelder[i][iz][ip][ir] = 0.;
	else 
	  {
	    bes1 = j1(gamma11*r/p->radius);
	    val = amplitude*bes1*data->coordinate[2][ir]*
	      (cos(data->coordinate[1][ip]-theta0-M_PI*.5))/r;
	    /*printf("Val = %g\n",val);*/
	    for(i=0;i<3;i++)
	      data->dddfelder[i][iz][ip][ir] = val; 
	 }
      }

  }
}
/***************************************************************************/

void peritest(HDF_DS *data, PARA *p)
{
int ir,ip,iz,i;
double mod;

for(i=0;i<3;i++) 
  for(iz=0;iz<data->nz;iz++)
    for(ip=0;ip<data->ny;ip++)
      for(ir=0;ir<data->nx;ir++)
	data->dddfelder[i][iz][ip][ir] = (double)sin(6.*M_PI/(double)(data->ny)*ip);
	
}

/***************************************************************************/

void bdadjust(HDF_DS *data, PARA *p)
{
int ir,ip,iz,i;
double dn,r,l;

for(i=0;i<3;i++)
  {
    for(iz=0;iz<data->nz;iz++)
      for(ip=0;ip<data->ny;ip++)
	  for(ir=0;ir<data->nx;ir++)
	    {
	      r   = ((double)(ir))/(double)(data->nx);
 
	      l    = tanh((1.- r)*5)*tanh(r*5.);
	      data->dddfelder[i][iz][ip][ir] *= l ;
	    }
	
  }
}
/***************************************************************************/

void drift_Util_Noise(HDF_DS *data, PARA *p)
{
int ir,ip,iz,i;
double mod,fac;


  for(iz=0;iz<data->nz;iz++)
    for(ip=0;ip<data->ny;ip++)
      for(ir=0;ir<data->nx;ir++)
	for(i=0;i<3;i++) 
	{
	data->dddfelder[i][iz][ip][ir] = p->amp_random[i]*((double)rand()/(double)RAND_MAX-.5);
	rand();rand();;rand();
      }	


}

/******************************************************************************/


void random1(HDF_DS *data, PARA *p)
{	
double 	
  norm,kz,kx,ky,
  r,r0,theta0 = 0.,mod,
  dkz,dkx,dky,
  kx2,ky2,khx,khy,kymin,
  fac,phase;
register int 	nzmodes,i,j,ix,iy,iz,l;
double **tfeld;

/* Erzeugen einer zuf"alligen Verteilung mit Gaussprofil im Fourierraum		
    Keine Vorbelegung von kx=0,ky=0  Werten !!!
*/

nzmodes = 2.;
norm = 1./(double)(nzmodes-1);



tfeld = Util_DMatrix(data->ny,0,data->nx,0);

dkz = 2.0*M_PI/(p->zmax - p->zmin);
p->dz = (p->zmax - p->zmin)/(double)data->nz;

printf("Vor erster FFT\n");

initeoft(data->nx,data->ny);	
 


	  
for(l=3;l>=0;l--)
  {
    i = l%3;
    for(iz=0;iz<data->nz;iz++)
      for (iy=0;iy<data->ny;iy++) 
	for (ix=0;ix<data->nx;ix++)  data->dddfelder[i][iz][iy][ix] = 0.;
    
    for (j=1;j<nzmodes;j++) 
    {
      kz = (double)dkz;

      if(p->boundary[i] == 1)
	{
	  dkx =     M_PI/(p->xmax - p->xmin);
	  dky = 2.0*M_PI/(p->ymax - p->ymin);
	  khx = p->width_random_x[i]*p->width_random_x[i]*dkx*dkx;
	  khy = p->width_random_y[i]*p->width_random_y[i]*dky*dky;
	 	  
	  for (iy=0,ky=0.;iy<data->ny;iy+=2,ky += dky) 
	    {
	      ky2 = ky*ky;
	      for (ix=0,kx=0.;ix<data->nx;ix++,kx += dkx) 
		{
		  kx2 = kx*kx;
		  phase= (double) rand()*2.*M_PI/(double)RAND_MAX;
		  fac = (double) exp(- (kx2*khx + ky2*khy));
		  tfeld[iy][ix]  = norm*sin(phase)*fac*p->amp_random[i]*.5*(double)(data->nx*data->ny/4);
		  tfeld[iy+1][ix]= norm*cos(phase)*fac*p->amp_random[i]*.5*(double)(data->nx*data->ny/4);
		  
		  if((iy == 0)|| (ix == 0) )
		    {
		      tfeld[iy][ix] = 0.;
		      tfeld[iy+1][ix] = 0.;
		    }
		}
	    }
	  BUGREPORT;
	  FBACKD(tfeld,data->nx,data->ny);


	  /* Move field half a point to the right ...*/

	  for (iy=0;iy<data->ny;iy++) 
	    {
	      ky=tfeld[iy][0];
	      for (ix=0;ix<data->nx-1;ix++) 
		tfeld[iy][ix] =0.5*(tfeld[iy][ix]+tfeld[iy][ix+1])  ;
	      tfeld[iy][data->nx-1] = 0.5*(tfeld[iy][data->nx-1]+ky) ;
	    }

	}
      else
	{
	  dkx = 2.0*M_PI/(p->xmax - p->xmin);
	  dky = 2.0*M_PI/(p->ymax - p->ymin);
	  kymin =-(double)(data->ny)*.5*dky;
	  
	  khx = 1./(p->width_random_x[i]*p->width_random_x[i]*dkx*dkx);
	  khy = 1./(p->width_random_y[i]*p->width_random_y[i]*dky*dky);
	  BUGREPORT;	
	  
	  for (iy=0,ky=kymin;iy<data->ny;iy++,ky += dky) 
	    {
	      ky2 = ky*ky;
	      for (ix=0,kx=dkx;ix<data->nx;ix+=2,kx += dkx) 
		{
		  kx2 = kx*kx;
		  phase= (double) rand()*2.*M_PI/(double)RAND_MAX;
		  fac = (double) exp(- (kx2*khx + ky2*khy));
		  tfeld[iy][ix] = norm*sin(phase)*fac*p->amp_random[i]*.5;
		  tfeld[iy][ix+1] = norm*cos(phase)*fac*p->amp_random[i]*.5;
		}
	    }
	  
	  for (iy=0;iy<data->ny;iy++)
	    {
	      tfeld[iy][0]=0.;
	      tfeld[iy][1]=0.;
	    } 
	  
	  for (ix=0;ix<data->nx;ix++)
	    tfeld[data->ny/2][ix] = 0.;
	  
	  FBACK(tfeld,data->nx,data->ny);
	}
    


      fac = (double) rand()*2.*M_PI/(double)RAND_MAX;
      mod = dkz*dkz/(kz*kz);

     for(iz=0;iz<data->nz;iz++)
	{
	  kx2 = sin(fac+dkz*iz*p->dz);
	  for (iy=0;iy<data->ny;iy++) 
	    for (ix=0;ix<data->nx;ix++) 
	      data->dddfelder[i][iz][iy][ix] +=/* mod* */ tfeld[iy][ix]*kx2;
	}






    /*
    if (p->radius > 0.)
      {
    	
	for(iz=0;iz<data->nz;iz++)  
	  for (iy=0;iy<data->ny;iy++) 
	    for (ix=0;ix<data->nx;ix++) 
		data->dddfelder[i][iz][iy][ix] *= exp(-(double)((data->nx-ix)*(data->nx-ix))*1./400.);


	if(i == 1 || i == 2)
	for(iz=0;iz<data->nz;iz++)  
	  for (iy=0;iy<data->ny;iy++) 
	    for (ix=0;ix<data->nx;ix++) 
		data->dddfelder[i][iz][iy][ix]+= p->radius*(1.+ exp(-(double)(ix*ix)*1./400.));
      }
    */		
    }
  }
}

/***************************************************************************/

/******************************************************************************/

 
void map_Util_Noise(HDF_DS *data, PARA *p)
{ 
register int ix,iy,iz,i,j,k;
double r,phi,x,y,dxp,dyp,cs,ss;
double **feld;
double mod;
int nor=64;


double 	kx,ky,r0,theta0 = 0.,
	dkx,dky,
	kx2,ky2,khx,khy,
	fac,phase,amp,norm;

 
 double *Sx,*Sy,*xpos,*ypos,*vals;




 Sx = (double *)calloc( (size_t)(nor+10),sizeof(double) );
 Sy = (double *)calloc( (size_t)(nor+10),sizeof(double) );

 xpos = (double *)calloc( (size_t)(data->nx+2),sizeof(double) );
 ypos = (double *)calloc( (size_t)(data->nx+2),sizeof(double) );
 vals = (double *)calloc( (size_t)(data->nx+2),sizeof(double) );


 feld = Util_DMatrix(nor,0,nor,0);
 dkx = M_PI;
 dky = M_PI;    

 for(i=0,kx=0.;                      i<nor;i+=2,kx+=dkx) Sx[i] = kx;

 for(i=0,ky=-.5*(double)nor*dky;i<nor;i++, ky+=dky) Sy[i] = ky;




/* Map rectangular Util_Noise to circle */
 
 for(i=0;i<3;i++)
   {
 
      
     khx = 1./(p->width_random_x[i]*p->width_random_x[i]*dkx*dkx);
     khy = 1./(p->width_random_y[i]*p->width_random_y[i]*dky*dky);
     
     for (iy=0;iy<nor;iy++) 
       {
	 ky = Sy[iy];
	 ky2 = ky*ky;
	 for (ix=2;ix<nor;ix+=2) 
	   {
	     kx = Sx[ix];
	     kx2 = kx*kx;
	     phase= (double) rand()*2.*M_PI/(double)RAND_MAX;
	     fac = (double) exp(- (kx2*khx + ky2*khy));
	     feld[iy][ix] = sin(phase)*fac*p->amp_random[i]*.5;
	     feld[iy][ix+1] = cos(phase)*fac*p->amp_random[i]*.5;
	   }
       }
     
     
     for (iy=0;iy<nor;iy++)
       feld[iy][0]= feld[iy][1] = 0.;
     
     for (iy=0;iy<data->ny;iy++)
	 {   
	   phi = data->coordinate[1][iy];
	   cs = cos(phi);
	   ss = sin(phi);

	   for (ix=0;ix<data->nx;ix++)
	     {
	       r   = data->coordinate[2][ix]/p->xmax;
	       xpos[ix] = r*cs;
	       ypos[ix] = r*ss;
	     }
	   
	       spectral_interpol(&data->dddfelder[i][0][iy][0],xpos,ypos,data->nx,feld,nor,nor,Sx,Sy);
	 }
 
     for (iy=0;iy<data->ny;iy++)
       for (ix=0;ix<data->nx;ix++)
	 {
	   r   = (data->coordinate[2][ix]-p->xmin)/(p->xmax-p->xmin); 
	   data->dddfelder[i][0][iy][ix] *= tanh((1.- r)*10.);
	 }


     if(p->xmin != 0.)
      for (iy=0;iy<data->ny;iy++)
       for (ix=0;ix<data->nx;ix++)
	 {
	   r   = (data->coordinate[2][ix]-p->xmin)/(p->xmax-p->xmin); 
	   data->dddfelder[i][0][iy][ix] *= tanh(r*10.);
	 }      
   

     
	for(iz=data->nz-1;iz>=0;iz--)
	  {
	    mod =  (1. + .1*((double)rand()/(double)RAND_MAX -0.5));
	    for (iy=0;iy<data->ny;iy++) 
	      for (ix=0;ix<data->nx;ix++)
		data->dddfelder[i][iz][iy][ix] =  mod*data->dddfelder[i][0][iy][ix];
	  }
     
   }


 /*
for(i=1;i<2;i++)
  for(iz=data->nz-1;iz>=0;iz--)
    for (iy=0;iy<data->ny;iy++) 
      for (ix=0;ix<data->nx;ix++)
	data->dddfelder[i][iz][iy][ix] +=  5.*(double)(data->nx-ix)/(double)(data->nx);
	
 */



}
/*************************************************************************************************/

void scott1(HDF_DS *data, PARA *p)
{
int ir,ip,iz,i;
 double x,y,z;
double mod,fac;

 for(i=0;i<3;i++) 
  for(iz=0;iz<data->nz;iz++)
    for(ip=0;ip<data->ny;ip++)
      for(ir=0;ir<data->nx;ir++)
	{
	  x =((double)ir+0.5)*p->dx+p->xmin;
	  y =((double)ip+0.5)*p->dy+p->ymin;
	  z =((double)iz+0.5)*p->dz+p->zmin; 

	  data->dddfelder[i][iz][ip][ir] = p->amp_random[i]*exp(-(x*x+y*y)/36.)*exp(-z*z);
	}
}

/******************************************************************************/


void dent(HDF_DS *data, PARA *p)
{
int ir,ip,iz,i;
double mod,fac;
 double x,y,z;

 for(i=0;i<3;i++) 
  for(iz=0;iz<data->nz;iz++)
    for(ip=0;ip<data->ny;ip++)
      for(ir=0;ir<data->nx;ir++)
	{
	  x =((double)ir+0.5-data->nx/2)/(double)data->nx;
	  y =((double)ip+0.5-data->ny/2)/(double)data->ny;
	  z =((double)iz+0.5)*p->dz+p->zmin; 

	  data->dddfelder[i][iz][ip][ir] = p->amp_random[i]*exp(-(x*x+y*y)*40.)*exp(-50.*z*z*z*z);

	}
}


/******************************************************************************/


void scott2(HDF_DS *data, PARA *p)
{
int ir,ip,iz,i;
double mod,fac;
 double x,y,z;

 for(i=0;i<3;i++) 
  for(iz=0;iz<data->nz;iz++)
    for(ip=0;ip<data->ny;ip++)
      for(ir=0;ir<data->nx;ir++)
	{
	  x =((double)ir+0.5)*p->dx+p->xmin;
	  y =((double)ip+0.5)*p->dy+p->ymin;
	  z =((double)iz+0.5)*p->dz+p->zmin; 

	  data->dddfelder[i][iz][ip][ir] = cos(M_PI/64.*x)*cos(M_PI/32.*y)*cos(z);
	}
}
