/*
  $Log: utilities.c,v $
*/
/* If the compiler supports SSE2, use that in the procedure Util_3DArakawaNl */
#ifdef USE_SSE2
#ifdef __SSE2__
#include <emmintrin.h>
#endif
#endif

#undef DEBUG
#undef COMM_STATUS
#include <utilities.h>

#define CALLOC Util_AllocChk

#include "util_arakawa.c"
#include "util_timestep.c"
#include "util_boundary.c"


/*********************************************************/
/* 
   Checks the cfl to second order and returns -1 on failure 
   and prints an error message 
     
   gets a streamfunction as input and assumes zero boundary velocity

   returns poloidal and radial velocity in vr and vp
   as well as the maximum value of cfr and cfl


   uses glasser coordinates

   / v_r \     /    v'*kai ( grf d_r \phi +  gff d_phi \phi)  \
   \vec u = |     |   = |                                              |
   \ v_p /     \  - v'*kai ( grr d_r  \phi + grf d_phi \phi)  /
*/


int  Util_CheckCFLMetric(double **f, double **fvr, double **fvp, HDF_DS *data,PARA *para, 
                         double *hval,double *vval, double *kai, double *grr, double *gff, double *grf,
                         double *norm_gr, double *norm_gp, double *cfr, double *cfp)
{
    int 
        ip,ir,np,result=0;
    double 
        maxvp=0,maxvr=0,
        vp,vr;

    COMM(fprintf(stderr,"Util_2DCheckCFL: Local Dimensions are (%ld,%ld), y- should be global...\n",data->lnx,data->lny););

    np = data->lny;

    for(ip=0;ip<data->lny;ip++)
        for(ir=0;ir<data->lnx;ir++)
        {
            vr = fvr[ip][ir] 
                =    -norm_gp[ir]*0.5*
                (grf[ir]*hval[ir]*(f[ip  ][ir+1]       -  f[ip  ][ir-1])
                 + gff[ir]*vval[ir]*(f[(ip+np+1)%np][ir] -  f[(ip-1+np)%np][ir]));
            maxvr = MAX(maxvr,vr*vr);
        }


    for(ip=0;ip<data->lny;ip++)
        for(ir=0;ir<data->lnx;ir++)
        {
            vp = fvp[ip][ir] 
                =      norm_gr[ir]*0.5*
                ( grr[ir]*hval[ir]*(f[ip  ][ir+1]        -  f[ip  ][ir-1])+
                  grf[ir]*vval[ir]*(f[(ip+np+1)%np][ir  ] - f[(ip-1+np)%np][ir  ]));
            maxvp = MAX(maxvp,vp*vp);
        }

    *cfr = sqrt(maxvr)*para->dt/para->adrhos;
    *cfp = sqrt(maxvp)*para->dt/para->adrhos;
 
  
    COMM(fprintf(stderr,"CFL_R = %f, CFL_P = %f\n",(*cfr),(*cfp)););

    if(*cfr > 2./3. || *cfp > 2./3.)
    {     
        fprintf(stderr," CFL VIOLATION: t= %f, CFL_P = %f, CFL_R = %f dt = %f\n",para->time,*cfr,*cfp,para->dt);
        result = -1;
    }
  
    return result;
}
/************************************************************************/
int Util_3DCheckCFL(double ***f, double ***fvr, double ***fvp, 
                     HDF_DS *data,PARA *para, double *hval, double *vval,
                     double *cfr, double *cfp)
{
    int i;
    int result=0;

    for(i=0;i<data->nz;i++)
    {
        if((result = Util_2DCheckCFL(f[i],fvr[i],fvp[i], data,para, hval, vval,cfr,cfp)) == -1) return result;
    }

    return result;   


}
/***********************************************************************/
/* Integrate value to the power of order                               */

void Util_Integral(double **f,int order, double *norm,int np,int nr,double *result)
{
    register int ir,ip,i=0,n;
    double val=0.;
    static double *help=NULL;
    static int oldn =0;
    
    n = nr*np;
    
    if( n > oldn) 
    {
        if(help != NULL) free(help);
        help =    Util_DVector(n,0);
        oldn = n;
        BUGREPORT;
    }

    for(i=0;i<n;i++)  help[i] = 1.;

    for(i=0;i<order;i++) for(ip=0;ip<np;ip++) for(ir=0;ir<nr;ir++) help[ip*nr+ir] *= f[ip][ir];

    /* Up to linear interpolation */
  
    for(ip=0;ip<np;ip++) for(ir=0;ir<nr;ir++) val += help[ip*nr+ir]*norm[ir];
  
    *result = val;
}



/***********************************************************************/
/* Integrate value to the power of order                               */

void Util_3DIntegral(double ***f,int order, double *norm,int nx,int ny, int nz,double **help,double *result)
{
    int i;
    double sum = 0;

    *result = 0;
    for(i=0;i<nz;i++)
    {
        Util_Integral(f[i],order,norm,ny,nx,&sum);
        *result += sum;
    }

}

/**********************************************************************/


int  Util_2DCheckCFL(double **f, double **fvr, double **fvp, HDF_DS *data,PARA *para, double *hval, double *vval,
                     double *cfr, double *cfp)
{
    int 
        ip,ir,np,result=0;
    double 
        maxvp=0,maxvr=0,
        vp,vr;

    COMM(fprintf(stderr,"Util_2DCheckCFL: Local Dimensions are (%ld,%ld), y- should be global...\n",data->lnx,data->lny););

    np = data->lny;

    for(ip=-1;ip<=data->lny;ip++)
        for(ir=0;ir<data->lnx;ir++)
        {
            vp = fvp[ip][ir] 
                =      0.5*hval[ir]*(f[ip  ][ir+1] -  f[ip  ][ir-1]);
            vp*=vval[ir];
            maxvp = MAX(maxvp,vp*vp);
        }



    for(ip=-1;ip<=data->lny;ip++)
        for(ir=0;ir<data->lnx;ir++)
        {
            vr = fvr[ip][ir] 
                =      -0.5*vval[ir]*(f[(ip+np+1)%np][ir  ] - f[(ip-1+np)%np][ir  ]);
            vr*=hval[ir];
            maxvr = MAX(maxvr,vr*vr);
        }

    *cfr = sqrt(maxvr)*0.5*para->dt*0.01;
    *cfp = sqrt(maxvp)*0.5*para->dt*0.01;
 
    if(*cfr > 1. || *cfp > 1.)
    {     
        /*fprintf(stderr," CFL VIOLATION: t= %f, CFL_P = %f, CFL_R = %f dt = %f\n",para->time,*cfr,*cfp,para->dt);*/
        result = -1;
    }
    return result;
}
/*************************************************************************/


void Util_CalcParameters(PARA *p,HDF_DS *d)
{
    double l_gamma = 1.;
    double Bgauss=0.;
    double n_cgs;
    double lamda_ei =0.;
    double boltzmann_k_si = 1.3807e-23;        // J/K
    double boltzmann_ev_per_kelvin =  8.617332478e-5; //eV/K
    double particles_per_pascal_ev = 6.24150934e18;
    
        
    double energy_per_ev = 1.1604e4;           // K
    double m_proton = 1.6726e-27;              // kg
    const char* sformat = "%-20s";
    double n_n=0.;                             // Neutral density from neutral pressure
   
#define ME (1./1.8362e3)

/* Formulas taken from NRL */ 
/* No check if all parameters are consistent */
   
    if(0. >=  p->Te ) return; // If we have no electron temperature we have no plasma, do not calculate
    

    Bgauss = p->B0*1.e4;
    n_cgs  = p->n0*1.e13;
    

    p->rho_s =  1.02 *sqrt(p->Mi)*sqrt(p->Te)/p->Z/Bgauss;             //     in [m]
    p->rho_i =  1.02 *sqrt(p->Mi)*sqrt(p->Ti)/p->Z/Bgauss;             //     in [m]
    p->rho_e =  2.38*sqrt(p->Te)/Bgauss/100.;                          //     in [m]

    p->omega_ci = 9.58 * 1.e3 * p->Z/p->Mi *Bgauss;                    //     in [rad/s]
    p->omega_ce = 1.76 * 1.e7 * Bgauss;                                //     in [rad/s]

    p->omega_pi = 1.32* 1e3 * p->Z/sqrt(p->Mi)*sqrt(n_cgs)   ;         //     in [rad/s]

    p->v_thi    = 9.79* 1.e5/sqrt(p->Mi)*sqrt(p->Ti)/100. ;            //     in [m/s]
    p->v_the    = 4.19* 1.e7* sqrt(p->Te)/100.;                        //     in [m/s]
    p->c_s      = 9.79* 1.e5 * sqrt(l_gamma*p->Z*p->Te /p->Mi)/100.;   //     in [m/s]
    p->v_alfven = 2.18* 1.e11 /sqrt(p->Mi) /sqrt(n_cgs) /100.*Bgauss;  //     in [m/s]


    fprintf(stderr,"Equation Parameters:\n");
    fprintf(stderr,"Lengths:\n");
   
    fprintf(stderr,"\t%-20s = %f [cm] \n",        "rho_s",p->rho_s*100.);
    fprintf(stderr,"\t%-20s = %f [cm] \n",        "rho_i",p->rho_i*100.);

    fprintf(stderr,"Frequencies:\n");
    fprintf(stderr,"\t%-20s = %f [10**6 rad/s (Mhz)] \n","omega_ci",p->omega_ci*1.e-6);
    fprintf(stderr,"\t%-20s = %f [10**9 rad/s (Ghz)] \n","omega_ce",p->omega_ce*1.e-9);
    fprintf(stderr,"\t%-20s = %f [10**9 rad/s (Ghz)] \n","omega_pi",p->omega_pi*1.e-9);

    fprintf(stderr,"Velocities:\n");
    fprintf(stderr,"\t%-20s = %f [km/s] \n","v_thi",p-> v_thi/1000.);
    fprintf(stderr,"\t%-20s = %f [km/s] \n","v_the",p->v_the/1000);
    fprintf(stderr,"\t%-20s = %f [km/s] \n","c_s",p->c_s/1000.);
    fprintf(stderr,"\t%-20s = %f [km/s] \n","rho_s*omega_ci",p->rho_s*p->omega_ci/1000);
    fprintf(stderr,"\t%-20s = %f [km/s] \n","v_alfven",p->v_alfven/1000.);
    
 
    /* Collisional parameters  BGK */

    lamda_ei = 24.- log(sqrt(n_cgs)/p->Te);

    // Use neutral pressure in pascal to calc neutral density at room temperature 300K
    p->n_n = n_n = p->p_n*particles_per_pascal_ev/(boltzmann_ev_per_kelvin*300.)*1.e-19; // This is per cubic m!! 

    
    // Print neutral density and plasma density
    //fprintf(stderr,"\t%-20s = %f %f  \n","Plasma density and neutral",p->n0,p->n_n);


    //NRL approximations 
    p->nu_ei = 2.91*1.e-6*n_cgs/(p->Te*sqrt(p->Te))*lamda_ei;
    p->nu_en = n_n* 5.e-15 *sqrt(boltzmann_k_si* energy_per_ev*p->Te/(ME*m_proton))*1.e19*1.e-4;
    p->nu_in = n_n* 5.e-15 *sqrt(boltzmann_k_si* energy_per_ev*p->Ti/(p->Mi*m_proton))*1.e19*1.e-4;
    
    //   fprintf(stderr,"\nn_cgs: %f 1e12 [cm^-3]",n_cgs/1e12);
    //   fprintf(stderr,"\nCoulomb Logarithm: %f\n",lamda_ei);
    

    fprintf(stderr,"\nFrequencies at n_n= %f [10**19*m-3]:\n",n_n);
    fprintf(stderr,"\t%-20s = %g [10**6 1/s] (MHz) \n","nu_ei",p->nu_ei*1e-6);
    fprintf(stderr,"\t%-20s = %g [10**6 1/s] (MHz) \n","nu_en",p->nu_en*1e-6);
    fprintf(stderr,"\t%-20s = %g [10**6 1/s] (MHz) \n\n","nu_in",p->nu_in*1.e-6);
    
    // Beta 
    p->beta = 0.5* 4.03*1e-11*n_cgs*(p->Te)/Bgauss/Bgauss; //This is the definition with 4 pi, NRL uses 8 pi .....

}

/*
  Print structures to file for debugging purposes
 */

void Util_PrintStructures(HDF_DS *d,PARA *p, int line_no)
{
  FILE * file;
  char file_name[DEFSTRLEN];
  int i,j;
  


  sprintf(file_name,"HDF_DS_p%d_l%d.dat",(int)d->this_process,(int)line_no);

  file=fopen(file_name,"w");
  for(i=0;i<4;i++)
      for(j=0;j<2;j++)
          fprintf(file,"Range[%d][%d]:  %f\n",i,j,d->range[i][j]);
  fprintf(file,"\n");

  fprintf(file,"OFF %d %d %d\n",(int)d->offx,(int)d->offy, (int)d->offz);
  fprintf(file,"NX %d %d %d\n",(int)d->nx,(int)d->ny,(int) d->nz);
  fprintf(file,"LNX %d %d %d\n",(int)d->lnx,(int)d->lny, (int)d->lnz);	  
  fprintf(file,"RANK %d %d\n",(int)d->isfloat,(int)d->rank);

  fprintf(file,"Dims %d %d %d %d\n",(int)d->dims[0],(int)d->dims[1],(int) d->dims[2],(int)d->dims[3]);
  fprintf(file,"Start %d %d %d %d\n",(int)d->start[0],(int)d->start[1],(int) d->start[2],(int)d->start[3]);
  fprintf(file,"End %d %d %d %d\n",(int)d->end[0],(int)d->end[1],(int) d->end[2],(int)d->end[3]);
  fprintf(file,"Edge %d %d %d %d\n",(int)d->edge[0],(int)d->edge[1],(int) d->edge[2],(int)d->edge[3]);
  fprintf(file,"N %d %d %d %d\n",(int)d->N[0],(int)d->N[1],(int) d->N[2],(int)d->N[3]);
  fprintf(file,"\n");

  fprintf(file,"Jobid: %s\n",d->jobid);
  fprintf(file,"CWD: %s\n",d->cwd);
  fprintf(file,"Desc: %s\n",d->desc);
  fprintf(file,"Filename: %s\n",d->filename); 
  fprintf(file,"Name: %s\n",d->name);
  fprintf(file,"Coordsys: %s\n",d->coordsys);
  fprintf(file,"Start_date: %s\n",d->start_date);
  fprintf(file,"Revision: %s\n",d->revision); 
  fprintf(file,"Integrator: %s\n",d->integrator);
  fprintf(file,"Compile_date: %s\n",d->compile_date);
  fprintf(file,"Maschine: %s\n",d->maschine);
  fprintf(file,"Write_date: %s\n",d->write_date); 
  fprintf(file,"Name_out: %s\n",d->name_out);
  fprintf(file,"Name_in: %s\n",d->name_in); 
  fprintf(file,"Erh_name: %s\n",d->erhname);

  fprintf(file,"\n");
  fprintf(file,"Grid_coords: %d %d %d\n",d->grid_coords[0],d->grid_coords[1],d->grid_coords[2]);

  fprintf(file,"\n\n");
	  
  fprintf(file,"P_time: %f %f %f %f\n",p->time, p->dt, p->end_time, p->out_time);
  fprintf(file,"P_xmax: %f %f %f %f\n",  p->xmax, p->xmin, p->dx, p->dkx);
  fprintf(file,"P_ymax: %f %f %f %f\n",  p->ymax, p->ymin, p->dy, p->dky);
  fprintf(file,"P_zmax: %f %f %f %f\n",  p->zmax, p->zmin, p->dz, p->dkz); fprintf(file,"Dims %d %d %d %d\n",(int)d->dims[0],(int)d->dims[1],(int) d->dims[2],(int)d->dims[3]);
  fprintf(file,"\n");

  fprintf(file,"P_kappan: %f %f %f %f\n",p->kappan, p->kappat, p->alpha, p->beta); 
  fprintf(file,"P_delta: %f %f %f %f\n",  p->delta, p->gamma, p->sigma, p->nu);
  fprintf(file,"P_r0: %f %f %f %f\n", p->r0, p->adrhos, p->betahat, p->muehat);
  fprintf(file,"P_mue_w: %f %f %f \n",  p->mue_w, p->mue_n, p->mue_t);
  fprintf(file,"\n");
 
  fprintf(file,"P_nprof: %f %f %f %f\n", p->nprof, p->bprof, p->tprof, p->phiprof); 
  fprintf(file,"P->source: %f %f %f %f %f\n", p->source, p->limiter, p->shat0, p->q0, p->k0);
  fprintf(file,"P->xbnd: %d %d %d %d %\n", p->xbnd,  p->ybnd,  p->zbnd,  p->otmult);
  fprintf(file,"P->roffset: %d %d %d %d\n", p->r_offset,  p->z_offset,  p->r_spacing,  p->z_spacing);
 
  fclose(file);

}
