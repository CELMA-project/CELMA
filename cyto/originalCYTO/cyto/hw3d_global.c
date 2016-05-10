/***********************************************************************/
/* This programm solves global HW-like equations for 
   density 

   dt N + {phi, N} =       -(V dzN +dzV) + mu_n (lap N + (nabla N)^2)

   vorticity
 
   dt w + {phi,w}  =     dzU -dzV + (U-V) dzN  + mu_w lap w -sigma w


   As a test to compare with the 2D ESEL code simple curvature Terms are build in. 

   Parallel Electron velocity

   nu V = dz phi - dzN

   Parallel Ion velocity

   dt u + {phi,u} = -dz phi -sigma u

   and with static relationship:


   w = nabla^2 phi 

   using the stiffy stable timestep  in 3D                   
   At r=0 a source injects plasma 
   while from R > 0  r_max we have limiter boundary conditions.
 
*/
/***********************************************************************/

#define HHSOLVER vnauplplr_3d
#define DRIVE
#undef VSTATIC
#undef DEBUG
#define PARALLEL_DAMPING


#define DX(F)  0.5*hval[ir]*( (F)[iz][ip][ir+1] - (F)[iz][ip][ir-1])
#define DY(F)  0.5*vval[ir]*( (F)[iz][ip+1][ir] - (F)[iz][ip-1][ir])

/* These Flags are for debugging */
#undef  COMM_STATUS

#include "../include/utilities.h"

 
/* Electron and Ion Mass */ 

#define ME 1.0
#define MI (4.*1800.)


#define IS_ELECTRONVELOCITY 1
#define IS_IONVELOCITY 2
#define IS_DENSITY 3
#define IS_POTENTIAL 4
#define IS_TEMPERATURE 5

/***********************************************************************/

#define FORYZ_BD for(iz=-data.offz;iz<nz+data.offz;iz++)  for(ip=-data.offy;ip<ny+data.offy;ip++) 
#define  FORALL_ZX         for(iz=0;iz<nz;iz++) for(ir=0;ir<nx;ir++)
#define  FORALL_ZX_BD         for(iz=-1;iz<=nz;iz++) for(ir=-1;ir<=nx;ir++)
#define  FORALL_ZY        for(iz=0;iz<nz;iz++) for(ip=0;ip<ny;ip++)
#define FORALL  for(iz=0;iz<nz;iz++) for(ip=0;ip<ny;ip++) for(ir=0;ir<nx;ir++)
#define FORALL_BD for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++)
#define  SOL_FORALL_XY_AND_BD  for(j=-1;j<=ny;j++) for(k=nsols;k<=nx;k++)
#define  SOL_FORALL   for(iz=0;iz<nz;iz++) for(ip=0;ip<ny;ip++) for(ir=nsols;ir<nx;ir++)

#define  FORALL_XY_BD  for(j=-1;j<=ny;j++) for(k=-1;k<=nx;k++)
#define FORALL_BD_XY for(i=0;i<nz;i++) for(j=-1;j<=ny;j++) for(k=-1;k<=nx;k++)

/***********************************************************************/
/* Local Functions */ 

void  HW3D_Background(double ***n_0,double **nbg,int nx,int ny,int nz);
void      HWG_DParallelSOL(int identity,double ***result,double ***result2,
                           double ***val,double ***f,double ***n,double ***V,double ***dzv,
                           HDF_DS *data, PARA *para);
void HWG_Laplace(double ***w,double ***f,double*hval,double*vval,
                 double *edrc,int nx,int ny,int nz,PARA *p,int ispolar);

/* void HWG_DSolveGlobVort(HDF_DS *data,PARA *p, */
/*                         double ***f,double **fa, double **fb, int bdf, */
/*                         double ***w,double **wa, double **wb, int bdw, */
/*                         double ***n,double **na, double **nb, int bdn, */
/*                         double ***res); */
void  HW3D_PolAverage(double ***s,double **t,int nx,int ny,int nz);


 
/***********************************************************************/
int main(int argc,char **argv)
{
    HDF_DS 
        data;
    PARA   
        para;
    FILE   
        *output;

    double
        ***res,***tmp,
        ***vr,***vp,
        ***res2;

    double ***omega;
 

    double 
        ***f_0,***fp,
        **fbdrup,**fbdrlow,
        **fbdra,**fbdrb; 

    double 
        ***V_0,***V_1,***V_2,
        ***dtV_0,***dtV_1,***dtV_2,
        **vbdrup,**vbdrlow,
        **vbdra,**vbdrb;


    double
        ***U_0,***U_1,***U_2,
        ***dtU_0,***dtU_1,***dtU_2,
        **ubdra,**ubdrb,
        **ubdrlow,**ubdrup;

    double
        ***w_0,***w_1,***w_2,
        ***dtw_0,***dtw_1,***dtw_2,
        **wbdra,**wbdrb,
        **ombdra,**ombdrb,
        **wbdrlow,**wbdrup;

    double
        ***n_0,***n_1,***n_2,
        ***dtn_0,***dtn_1,***dtn_2,
        **nbdra,**nbdrb,
        **enbdra,**enbdrb,
        **nbdrlow,**nbdrup;

    double
        ***t_0,***t_1,***t_2,
        ***dtt_0,***dtt_1,***dtt_2,
        **tbdra,**tbdrb,
        **tbdrlow,**tbdrup;

    double 
        ***vstatic,***exp_n_0=NULL,***exp_t_0=NULL,***b_0;


    double 
        ***dzU,***dzUN,***dzzU,
        ***dzV,***dzVN,***dzzV,
        ***dzF,***dzzF,
        ***dzN,***dzzN;


    double
        ***dxn,***dyn,
        ***dxf,***dyf,
        ***dtdxf,***dtdyf;

    double 
        ***nimp, ***drive,**nbg,***n_nobg;

    double *tprofile,*du,*dc,*dl;

    double
        *hval=NULL,*vval=NULL,*nlval=NULL,
        *gval=NULL,*rcor=NULL,*edrcor=NULL,
        *hs=NULL,*vs=NULL,*gs=NULL,*rs=NULL,*edrs=NULL,
        *norm=NULL,*dynbg=NULL;
 


    double  
        gamma_out,translim,fac,cflr, cflp,charge;

    double nu,nu_m,delta;
        
    int 	
        nx,ny,nz;

    register int  ir,ip,iz,i,j,k;

    int  ot=1,iter = 0,ispolar=FALSE,FIRST=TRUE;


    static char 
        rcsid[]="$Id: hw3d_global.c,v 4.1 2003/10/16 14:08:09 vona Exp $";

    double  
        lamda, gamma_n=0.,total_density = 0.;
    double 
        write_time;
    int   
        rbdcndw,rbdcndn,rbdcndnfluc,rbdcndf,rbdcndv,rbdcndu,rbdcndt;

    int   
        zbndcndw,zbndcndn,zbndcndf,zbndcndv,zbndcndu,zbndcndt;


    int 	restart,nsols;

/*****************************************************************************/

    FUtils_IniStructure(&data,&para,argc,argv,rcsid);
    data.anzahl = 1;

    para.codedesc =
        "# The HW_GLOB code solves Hasegawa-Wakatani like Equations in 3D:\n"\
        "#\n"\
        "# dt w + b[f,w] + [b,nT]                               =        + D_w*nabla^2 w \n"\
        "#\n"\
        "# dt N + b[f,N] + [b,T-f] + T[b,N]                     = + D_n*nabla^2 N + D_n((dN/dx)**2 + (dN/dy)**2)  \n"\
        "#\n"\
        "# dt T + b{f,T} - 2/3T[b,f] + 7/3T[b,T] + 2/3T*T[b,N]  = 5*c((T^+) - T) + S_T   - 5*lambda_T*(T-bg  ) + D_T*nabla^2 T  \n"\
        "#\n"\
        "# w = nabla^2 f \n"\
        "# b = 1/B: B is the magnetic field: \n"\
        "# N,T = ln(n,T) \n"\
        "# D_? is kinematic viscosity \n"\
        "#\n";


    para.desc = 
        "gamma: Artificial parallel damping;"\
        "sigma: ion neutral collisions;"\
        "nu: Coulomb collisions;"\
        "beta: -;"\
        "shear: -;"\
        "alpha: -;"\
        "delta: electron neutral collisions;"\
        "betahat: -;"\
        "adrhos: -;"\
        "mue_w: damping of vorticity;"\
        "mue_n: damping of density;"\
        "mue_t: damping of temperature;"\
        "kappan: -;"\
        "kappat: -;"\
        "r0: Length scale on magnetic field;"\
        "muehat: -;"\
        "bprof: magnetic field derivative;"\
        "tprof: -;"\
        "nprof: -;"\
        "phiprof: -;"\
        "hm_nl: -;"\
        "exb_ll: -;"\
        "limiter: SOLwidth;"\
        "source: -;"\
        "k0: -;"\
        "dt_pol: -;"\
        "uvortex: -;"\
        "radius: -;"; 

    para.codename = "hw_glob";
    data.rank = 3; 
    data.anzahl = 1;
    data.offz = 1;
    data.offx = 1;
    data.offy = 1; 

    restart = FUtils_ReadArguments(argc,argv,data.filename,&data.number,&data, &para);
    sprintf(data.erhname,"%s.erh",data.name_out);
    data.rank = 3; 
    data.anzahl = 1;
    data.offz = 1;
    data.offx = 1;
    data.offy = 1; 
 


    nx = data.nx;
    ny = data.ny;
    nz = data.nz;
    nsols = MIN(((-para.xmin)/para.dx),nx);
    if(nsols < 0) nsols = -1;
    

    /****************************************************/
    /*                                                  */
    /* WHAT ARE THE DEFINES                             */
    /*                                                  */
    /****************************************************/
      
    sprintf(data.desc,"SWITCHES:\n");

    strcat(data.desc,"Drive - ");
#ifdef DRIVE
    strcat(data.desc,"on.\n");
#else
    strcat(data.desc,"off.\n");
#endif
    fprintf(stderr,"Lz = %g,dz = %f\n", para.zmax,para.dz);
    fprintf(stderr,"%s",data.desc);


    /************************************************/
    /*                                              */
    /* setup of variables and fields                */
    /*                                              */
    /************************************************/
 
	
    /* Allocate Fields */

    Util_3DAllocFields(&w_0,&w_1,&w_2,&dtw_0,&dtw_1,&dtw_2,
                       &wbdra,&wbdrb,&wbdrlow,&wbdrup,
                       nx,ny,nz,data.offx,data.offy,data.offz);
   
    Util_3DAllocFields(&n_0,&n_1,&n_2,&dtn_0,&dtn_1,&dtn_2,
                       &nbdra,&nbdrb,&nbdrlow,&nbdrup,
                       nx,ny,nz,data.offx,data.offy,data.offz);
   
    Util_3DAllocFields(&t_0,&t_1,&t_2,&dtt_0,&dtt_1,&dtt_2,
                       &tbdra,&tbdrb,&tbdrlow,&tbdrup,
                       nx,ny,nz,data.offx,data.offy,data.offz);
   
   
    Util_3DAllocFields(&dzzF,&f_0,&dzN,&dzV,&dzF,&dzzN,
                       &fbdra,&fbdrb,&fbdrlow,&fbdrup,
                       nx,ny,nz,data.offx,data.offy,data.offz);
   
    Util_3DAllocFields(&V_0,&V_1,&V_2,&dtV_0,&dtV_1,&dtV_2,
                       &vbdra,&vbdrb,&vbdrlow,&vbdrup,
                       nx,ny,nz,data.offx,data.offy,data.offz);  

    Util_3DAllocFields(&U_0,&U_1,&U_2,&dtU_0,&dtU_1,&dtU_2,
                       &ubdra,&ubdrb,&ubdrlow,&ubdrup,
                       nx,ny,nz,data.offx,data.offy,data.offz);  
    
    dzVN     = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);   
    dzUN     = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);    
    dzU     = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);   
    dzzU     = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);   
    dzzV     = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);   
    drive     = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);   
    res         = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);     
    res2         = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);     
    vr         = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);     
    vp         = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);   


    dxn         = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);   
    dyn         = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);     
    dxf         = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);     
    dyf        = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);  

    dtdxf         = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);     
    dtdyf        = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);  

    exp_t_0         = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);     
    exp_n_0         = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);  
    vstatic         = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);    

    b_0        = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);    
    nimp        = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);  
    n_nobg        = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);  

    omega  = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz);     
   

    /* 2D Fields */
    nbg   =  Util_DMatrix(nz,data.offz,nx,data.offx); 
   
    enbdra = Util_DMatrix(nz,data.offz,ny,data.offy); 
    enbdrb = Util_DMatrix(nz,data.offz,ny,data.offy); 

    ombdra = Util_DMatrix(nz,data.offz,ny,data.offy);
    ombdrb = Util_DMatrix(nz,data.offz,ny,data.offy);

    /* Vectors */ 
   
    nlval           = Util_DVector(nx,data.offx);
    tprofile       = Util_DVector(nx,data.offx);

    du= Util_DVector(nz,data.offz);
    dc = Util_DVector(nz,data.offz);
    dl = Util_DVector(nz,data.offz);
    ispolar = Util_SetupGeom(&rcor,&edrcor,&hval,&vval,&gval,&rs,&edrs,&hs,&vs,&gs,&norm,&data,&para,1);
    
    fprintf(stderr,"The geometry is ");
    if(ispolar)fprintf(stderr,"polar.\n");
    else fprintf(stderr,"rectangular.\n");
    
/*
  
Some Constants

*/
    
    delta = para.delta;
    nu = para.nu;
    nu_m = para.nu*ME/MI;
    
    fprintf(stderr,"Electron neutral collisions        : %g\n",para.nu);
    fprintf(stderr,"Electron Ion collisions            : %g\n",para.delta);
    fprintf(stderr,"Ion neutral collisions             : %g\n",para.sigma);
    fprintf(stderr,"Viscous:                 : W = %g     N = %g\n",para.mue_w,para.mue_n);
    
    /* Value used for non-linearity */
    for(ir=-1;ir<=nx;ir++)     rcor[ir]  = para.xmin+((double)ir+0.5)*para.dx;
    for(ir=-1;ir<=nx;ir++)     nlval[ir] = hval[ir]*vval[ir];
    
    /* Drive */
    FORALL_BD  drive[iz][ip][ir] =  para.source* exp(- 0.1*(double)(ir-nx/2)*(ir-nx/2) -0.02*(double)(iz-nz/2)*(iz-nz/2) -0.02*(double)(ip-ny/2)*(ip-ny/2));


/* B-Field */
	FORALL_BD   b_0[iz][ip][ir] = 1.0 + para.bprof*(0.33 + sin((double)iz/(double)nz*2*M_PI)*para.r0*rcor[ir]);
   
    BUGREPORT;
    
    fprintf(stderr, "dx = %f\n",para.dx);
    for(ir=-1;ir<=nx;ir++) fprintf(stderr, "b0 = %f\n", b_0[0][0][ir]);

    /* Assumed Temperature profile, included to increase viscosity in the SOL */
  
    for(ir=-1;ir<=nx;ir++) tprofile[ir] = 1.; 


    /* 
       INITIALIZE BOUNDARY VALUES 
    */
  
    rbdcndf   = NEUDIR; 
    rbdcndw= NEUDIR;
    rbdcndn = NEUNEU;
    
    rbdcndnfluc   = DIRDIR; 
    rbdcndv   = NEUNEU; 
    rbdcndu   = NEUNEU;
    rbdcndt   = DIRDIR; 


    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) nbdra[iz][ip] = 0.0;
    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) nbdrb[iz][ip] = 0.0;

    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) tbdra[iz][ip] = 0.0;
    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) tbdrb[iz][ip] = 0.0;
 
    if(rbdcndn == DIRNEU){
        for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) enbdra[iz][ip]= exp(nbdra[iz][ip]);
        for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) enbdrb[iz][ip]= nbdrb[iz][ip];
	}	
    else if(rbdcndn == NEUDIR){
        for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) enbdra[iz][ip]= nbdra[iz][ip];
        for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) enbdrb[iz][ip]= exp(nbdrb[iz][ip]);
	}	
    else if(rbdcndn == DIRDIR){
        for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) enbdra[iz][ip]= exp(nbdra[iz][ip]);
        for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) enbdrb[iz][ip]= exp(nbdrb[iz][ip]);
	}	
    else {
        for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) enbdra[iz][ip]= nbdra[iz][ip];
        for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) enbdrb[iz][ip]= nbdrb[iz][ip];
	}	

    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) wbdra[iz][ip] = 0.;
    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) wbdrb[iz][ip] = 0.;

    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) vbdra[iz][ip] =  0.;
    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) vbdrb[iz][ip] =  0.;

    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) ubdra[iz][ip] = 0.;
    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) ubdrb[iz][ip] = 0.;

    /* Potential is relative to sheath potential */
    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) fbdra[iz][ip] =0.;
    for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) fbdrb[iz][ip] =0.;


    /*****************************************************************************/

    /* ZBDCND:   
    */

    zbndcndf = PERIODIC ; zbndcndw  = PERIODIC; zbndcndn = PERIODIC; 
    zbndcndt = PERIODIC; zbndcndv = PERIODIC;zbndcndu = PERIODIC;

    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) fbdrup[ip][ir] = 0.;
    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) fbdrlow[ip][ir] = 0.;

    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) vbdrup[ip][ir] = 0.;
    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) vbdrlow[ip][ir] = 0.;

    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) ubdrup[ip][ir] = 0.;
    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) ubdrlow[ip][ir] = 0.;

    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) nbdrup[ip][ir] = 0.;
    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) nbdrlow[ip][ir] = 0.;
 
    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) tbdrup[ip][ir] = 0.;
    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) tbdrlow[ip][ir] = 0.;
 
    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) wbdrup[ip][ir]  = 0.;
    for(ip=-1;ip<=ny;ip++) for(ir=-1;ir<=nx;ir++) wbdrlow[ip][ir] = 0.;

    /* 

    Copy initial value for n to n_0.
    We set initial potential and vorticity to zero

    */
 
    BUGREPORT;
 

    /* n gets initialized with the low value on the boundary at b plus a 3D Gaussian */
 
 
    if(restart)
    {
        BUGREPORT;
        FUtils_Read3DFieldbyName(w_0,"Vorticity",&data,&para);
        BUGREPORT;
        FUtils_Write3dSpace(w_0,"Vorticity",NULL,NULL,NULL,NULL,"RESTART",0,&data,&para,TRUE);
      
        FUtils_Read3DFieldbyName(n_0,"Density",&data,&para);
        FUtils_Write3dSpace(n_0,"Density",NULL,NULL,NULL,NULL,"RESTART",0,&data,&para,FALSE);
        FORALL_BD n_0[iz][ip][ir] = log(n_0[iz][ip][ir]);
     
        FUtils_Read3DFieldbyName(f_0,"Potential",&data,&para);
        FUtils_Write3dSpace(f_0,"Potential",NULL,NULL,NULL,NULL,"RESTART",0,&data,&para,FALSE);
      
        FUtils_Read3DFieldbyName(V_0,"V_0",&data,&para);
        FUtils_Write3dSpace(V_0,"V_0",NULL,NULL,NULL,NULL,"RESTART",0,&data,&para,FALSE);
                        
        FUtils_Read3DFieldbyName(U_0,"U_0",&data,&para); 
        FUtils_Write3dSpace(U_0,"U_0",NULL,NULL,NULL,NULL,"RESTART",0,&data,&para,FALSE);
  
        FUtils_Read3DFieldbyName(t_0,"T_e",&data,&para); 
        FUtils_Write3dSpace(t_0,"T_e",NULL,NULL,NULL,NULL,"RESTART",0,&data,&para,FALSE);
  

        FORALL_BD w_1[i][j][k] =   w_0[i][j][k];
        FORALL_BD n_1[i][j][k] =   n_0[i][j][k];
        FORALL_BD V_1[i][j][k] =   V_0[i][j][k];
        FORALL_BD U_1[i][j][k] =   U_0[i][j][k];
        FORALL_BD t_1[i][j][k] =    t_0[i][j][k];
 	
    }
    else
    { 
 	
        BUGREPORT;

        FORALL n_0[iz][ip][ir] = exp(-(0.05*(double)(ir-nx/2)*(double)(ir-nx/2)+0.01*(double)(iz-nz/2)*(double)(iz-nz/2)+0.1*(double)(ip-ny/2)*(double)(ip-ny/2)));

                   
        FORALL n_0[iz][ip][ir] 	= log( 1.0 +1.e-8* n_0[iz][ip][ir]);             
        
        FORALL_BD   w_0[iz][ip][ir] =0.0;
 
        /* Omega, w ,f, V, U  get initialized with zero values */
        FORALL_BD omega[iz][ip][ir]  = f_0[iz][ip][ir] = 0.;
        FORALL_BD  V_0[iz][ip][ir] = U_0[iz][ip][ir] = 0.;   
 
        /* Initialise temperature with 0 */
        FORALL_BD t_0[iz][ip][ir] = 0.;
 

/* Initialise V_0 and U_0 to be monotonous in the SOL region */

        SOL_FORALL   V_0[iz][ip][ir]= -0.0*cos(M_PI*sin(0.5*M_PI*((double)iz+0.5)/(double)nz));
        SOL_FORALL   U_0[iz][ip][ir]= 0.*tanh(10.*(  (  (double)iz+0.5)/(double)nz-0.5));
        SOL_FORALL   V_0[iz][ip][ir]= U_0[iz][ip][ir];
        
    }
    /*
      Put values into ghost points
    */
    Util_3DFullBd(n_0,&data,nx,ny,nz,nbdra,nbdrb,nbdrup,nbdrlow,hval,para.dz,rbdcndn,zbndcndn); 
    Util_3DFullBd(f_0,&data,nx,ny,nz,fbdra,fbdrb,fbdrup,fbdrlow,hval,para.dz,rbdcndf,zbndcndf);
    Util_3DFullBd(omega,&data,nx,ny,nz,wbdra,wbdrb,wbdrup,wbdrlow,hval,para.dz,rbdcndw,zbndcndw); 
    Util_3DFullBd(V_0,&data,nx,ny,nz,vbdra,vbdrb,vbdrup,vbdrlow,hval,para.dz,rbdcndv,zbndcndv);
    Util_3DFullBd(U_0,&data,nx,ny,nz,ubdra,ubdrb,ubdrup,ubdrlow,hval,para.dz,rbdcndu,zbndcndu); 
    Util_3DFullBd(t_0,&data,nx,ny,nz,tbdra,tbdrb,tbdrup,tbdrlow,hval,para.dz,rbdcndt,zbndcndt);
 
    /* Calculate omega = w  - d_r phi  d_r n*/
    FORALL  omega[iz][ip][ir]  = w_0[iz][ip][ir]-hval[ir]*hval[ir]*(f_0[iz][ip][ir+1]-f_0[iz][ip][ir-1])*(n_0[iz][ip][ir+1]-n_0[iz][ip][ir-1]);     
    Util_3DFullBd(omega,&data,nx,ny,nz,ombdra,ombdrb,wbdrup,wbdrlow,hval,para.dz,rbdcndw,zbndcndw); 
 

    FORALL_BD 	exp_t_0[iz][ip][ir] = exp(t_0[iz][ip][ir]) ;
    FORALL_BD 	exp_n_0[iz][ip][ir] = exp(n_0[iz][ip][ir]) ;

    /* calculate boundary for w  = d_r phi  d_r n */
    for(iz=0;iz<nz;iz++) 
        for(ip=-1;ip<ny+1;ip++) 
            wbdra[iz][ip] =   ombdra[iz][ip] +4.* hs[0]* hs[0]*(f_0 [iz][ip][0] - f_0 [iz][ip][-1])*(n_0 [iz][ip][0] - n_0 [iz][ip][-1]);
    BUGREPORT;
   /* calculate boundary for w  = d_r phi  d_r n */
    for(iz=0;iz<nz;iz++) 
        for(ip=-1;ip<ny+1;ip++) 
            wbdrb[iz][ip] =   ombdrb[iz][ip] +4.* hs[nx-1]* hs[nx-1]*(f_0 [iz][ip][nx] - f_0 [iz][ip][nx-1])*(n_0 [iz][ip][nx] - n_0 [iz][ip][nx-1]);
    BUGREPORT;
    Util_3DFullBd(w_0,&data,nx,ny,nz,wbdra,wbdrb,wbdrup,wbdrlow,hval,para.dz,rbdcndw,zbndcndw);

    /* Let the diffusion only work on the fluctuations */

    HW3D_Background(res,nbg,nx, ny, nz);
    FORALL_BD  n_nobg[iz][ip][ir]=res[iz][ip][ir]-nbg[iz][ir];


    /* Precalculate all parallel derivatives */
    FORALL_BD  dzV[iz][ip][ir]=0.0;

    HWG_DParallelSOL(IS_IONVELOCITY,dzU,dzzU,U_0,f_0,n_0,V_0,dzV,&data,&para);
    HWG_DParallelSOL(IS_DENSITY,dzN,dzzN,n_0,f_0,n_0,V_0,dzV,&data,&para);
    HWG_DParallelSOL(IS_POTENTIAL,dzF,dzzF,f_0,f_0,n_0,V_0,dzV,&data,&para);
#ifdef VSTATIC
    /*
      V       - Equation 
    */ 
    /* We neclect electron mass and em-effects a static relationship and
       not a dynamical equation,
                 
    */
    
    FORALL_BD 	res[iz][ip][ir] = para.nu*exp_n_0[iz][ip][ir];
    FORALL_BD	V_0[iz][ip][ir]  =  1./(para.delta+res[iz][ip][ir])
        *(MI*(dzF[iz][ip][ir]-dzN[iz][ip][ir]) +res[iz][ip][ir]*U_0[iz][ip][ir]);
#endif      
    HWG_DParallelSOL(IS_ELECTRONVELOCITY,dzV,dzzV,V_0,f_0,n_0,V_0,dzV,&data,&para);
    BUGREPORT;

    /*
      Write initial condition
    */
    BUGREPORT;
    FUtils_Write3dSpace(omega,"Vorticity",n_0,"Density",f_0,"Potential","ini",1,&data,&para,TRUE);
    FUtils_Write3dSpace(U_0,"U",drive,"Drive/Limiter",V_0,"V","ini",1,&data,&para,FALSE);
    FUtils_Write3dSpace(dzU,"CollProf",t_0,"T_e",b_0,"B-Field","ini",1,&data,&para,FALSE);
    FUtils_Write3dSpace(dzV,"dzV",dzzV,"dzzV",dzN,"dzN","ini",1,&data,&para,FALSE);
    FUtils_Write3dSpace(dzzN,"dzzN",dzF,"dzF",dzzF,"dzzF","ini",1,&data,&para,FALSE);    
    FUtils_Write3dSpace(dzzU,"dzzU",dzU,"dzU",NULL,NULL,"ini",1,&data,&para,FALSE);    
    BUGREPORT;   
/*****************************************************************/
/*                                                               */
/*                 End of Setup                                  */
/*                                                               */
/*****************************************************************/
/*          Start of Time Loop                                   */
/*****************************************************************/
    BUGREPORT;
 
    write_time     = para.time + para.out_time;
    while(para.time < para.end_time)
    {
        FORALL dzVN[iz][ip][ir] =  dzV[iz][ip][ir] +  V_0[iz][ip][ir]*dzN[iz][ip][ir]; 
        FORALL dzUN[iz][ip][ir] =  dzU[iz][ip][ir] +  U_0[iz][ip][ir]*dzN[iz][ip][ir]; 
   
        FORALL dxn[iz][ip][ir] = DX(n_0);
        FORALL dyn[iz][ip][ir] = DY(n_0);
     
        FORALL dxf[iz][ip][ir] = DX(f_0);
        FORALL dyf[iz][ip][ir] = DY(f_0);
     
        /* zero boundaries for dxf, dyf */
        Util_3DFullBd(dxf,&data,nx,ny,nz,wbdra,wbdrb,wbdrup,wbdrlow,hval,para.dz,NEUNEU,zbndcndw);
        Util_3DFullBd(dyf,&data,nx,ny,nz,wbdra,wbdrb,wbdrup,wbdrlow,hval,para.dz,DIRDIR,zbndcndw);
     
        BUGREPORT;   
        /*********************************************/
        /*       Density - Equation                  */ 
        /*********************************************/
        
            
        Util_3DArakawaNl(dtn_0,f_0,n_0,nlval,nz,nx,ny);		
        FORALL     dtn_0[iz][ip][ir] *= b_0[iz][ip][ir]; 
        BUGREPORT;     
  
        /* Curvature */
        
        FORALL_BD res[iz][ip][ir]  = -f_0[iz][ip][ir] + exp_t_0[iz][ip][ir];
        Util_3DArakawaNl(res2,b_0,res,nlval,nz,nx,ny);
        FORALL     dtn_0[iz][ip][ir] += res2[iz][ip][ir];
        BUGREPORT;       
        Util_3DArakawaNl(res,b_0,n_0,nlval,nz,nx,ny);
        FORALL     dtn_0[iz][ip][ir] += exp_t_0[iz][ip][ir]*res[iz][ip][ir];
        
        /* Parallel terms */
        BUGREPORT;
        FORALL dtn_0[iz][ip][ir] +=  -dzVN[iz][ip][ir];
     
        /* Source */
        BUGREPORT;
        FORALL  dtn_0[iz][ip][ir] += drive[iz][ip][ir]/exp_n_0[iz][ip][ir];
        
  
#ifdef PARALLEL_DAMPING
        FORALL       dtn_0[iz][ip][ir]+= para.gamma*(dzzN[iz][ip][ir]+dzN[iz][ip][ir]*dzN[iz][ip][ir]); 
#endif
        BUGREPORT;  

        Util_3DSsTimeStep(iter,nx,ny,nz,para.dt,para.mue_n,&lamda, 
                          n_2,n_0,n_1,n_2,dtn_0,dtn_1,dtn_2);

       /* Let the diffusion only work on the fluctuations */
        FORALL  nimp[iz][ip][ir] = exp(n_2[iz][ip][ir]);
        HW3D_Background(nimp,nbg,nx, ny, nz);
        FORALL   n_nobg[iz][ip][ir]=nimp[iz][ip][ir]-nbg[iz][ir];
        Util_3DFullBd(n_nobg,&data,nx,ny,nz,fbdra,fbdrb,nbdrup,fbdrlow,hval,para.dz,DIRDIR,zbndcndn);  

        vnauplplr_3d(&data,&para,nimp,n_nobg,fbdra,fbdrb,
                     hval,gval,rcor,DIRDIR,lamda,TRUE);
        FORALL n_2[iz][ip][ir] =log(nimp[iz][ip][ir]+nbg[iz][ir]);

       
        /*Calculate dt ln n .....with implicit term */
        /*
        HWG_Laplace(res,n_nobg,hval,vval,edrcor,nx,ny,nz,&para,ispolar);
        FORALL  nimp[iz][ip][ir] = para.mue_n/exp_n_0[iz][ip][ir]*res[iz][ip][ir];
     
        FORALL dtn_2[iz][ip][ir] = dtn_0[iz][ip][ir] + nimp[iz][ip][ir];*/
        BUGREPORT;
        /***********************************************************************/
        /* 
        
        Global Vorticity equation needs to be evaluated AFTER dt_n_0 is known
     
        */
        /***********************************************************************/
     
        Util_3DArakawaNl(dtw_0,f_0,omega,nlval,nz,nx,ny);
        /* there is a 1/b term on the convective term */
       FORALL     dtw_0[iz][ip][ir] *= b_0[iz][ip][ir];
        
        /* Curvature */
       FORALL_BD res2[iz][ip][ir]  = exp_n_0[iz][ip][ir]* exp_t_0[iz][ip][ir]; 
       Util_3DArakawaNl(res,b_0,res2,nlval,nz,nx,ny);
       FORALL     dtw_0[iz][ip][ir] += res[iz][ip][ir]/exp_n_0[iz][ip][ir];
        

        BUGREPORT;    
        /* Parallel Terms */
        FORALL dtw_0[iz][ip][ir] +=  dzUN[iz][ip][ir]-dzVN[iz][ip][ir];
        

        /* Collisional Term: Pedersen current */
/*        FORALL dtw_0[iz][ip][ir]+=  -para.sigma*omega[iz][ip][ir];*/
     
#ifdef PARALLEL_DAMPING
        FORALL       dtw_0[iz][ip][ir]-=para.gamma*dzzF[iz][ip][ir];
#endif

        /* \nabla d_t n \nabla \phi term */
        /* First zero boundaries for dt_n0 */
     
        Util_3DFullBd(dtn_2,&data,nx,ny,nz,fbdra,fbdrb,wbdrup,wbdrlow,hval,para.dz,NEUNEU,zbndcndw);
        FORALL dtw_0[iz][ip][ir] += (DX(dtn_2)*dxf[iz][ip][ir]+DY(dtn_2)*dyf[iz][ip][ir]);
     
        /* Nonlinear Term */ 
     
        Util_3DArakawaNl(dtdxf,f_0,dxf,nlval,nz,nx,ny);        
        Util_3DArakawaNl(dtdyf,f_0,dyf,nlval,nz,nx,ny);
        FORALL  dtdxf[iz][ip][ir] *=b_0[iz][ip][ir];	 
          FORALL  dtdyf[iz][ip][ir] *=b_0[iz][ip][ir];	   
     
        FORALL dtw_0[iz][ip][ir] += (dxn[iz][ip][ir]*dtdxf[iz][ip][ir]+dyn[iz][ip][ir]*dtdyf[iz][ip][ir]);
     
        /*
          Calculate nabla (sigma  phi - mue w) 
        */
        FORALL_BD res[iz][ip][ir] =   para.sigma*f_0[iz][ip][ir]- para.mue_w*omega[iz][ip][ir];
        FORALL dtw_0[iz][ip][ir] -= (dxn[iz][ip][ir]*DX(res) + dyn[iz][ip][ir]*DY(res));
        /* Time step */
     
        Util_3DSsTimeStep(iter,nx,ny,nz,para.dt,para.mue_w,&lamda, 
                          res,w_0,w_1,w_2,dtw_0,dtw_1,dtw_2);
          
        /* We treat viscosity implicitely */

        /* Damping only on vorticity */
        FORALL dtU_0[iz][ip][ir]  = dxn[iz][ip][ir]*dxf[iz][ip][ir]+ dyn[iz][ip][ir]*dyf[iz][ip][ir];
        FORALL  res[iz][ip][ir]  -=dtU_0[iz][ip][ir];
      
        vnauplplr_3d(&data,&para,w_2,res,wbdra,wbdrb,hval,gval,rcor,rbdcndw,lamda,TRUE);
        FORALL  w_2[iz][ip][ir]  +=dtU_0[iz][ip][ir];
     
        /*********************************************/
        /*      Ion Velocity - Equation                 */
        /*********************************************/
        BUGREPORT;
        Util_3DArakawaNl(dtU_0,f_0,U_0,nlval,nz,nx,ny);
        FORALL dtU_0[iz][ip][ir]*=b_0[iz][ip][ir];

        /* Curvature */
        
        FORALL_BD res[iz][ip][ir]  = -f_0[iz][ip][ir] + exp_t_0[iz][ip][ir];
        Util_3DArakawaNl(res2,b_0,res,nlval,nz,nx,ny);
        FORALL     dtU_0[iz][ip][ir] += res2[iz][ip][ir];

        FORALL dtU_0[iz][ip][ir]+=  -dzF[iz][ip][ir] -para.sigma*U_0[iz][ip][ir]
            +nu_m*exp_n_0[iz][ip][ir]*(V_0[iz][ip][ir]-U_0[iz][ip][ir]);
     
        /* SELFADVECTION  */
        FORALL dtU_0[iz][ip][ir]+=  -dzU[iz][ip][ir]*U_0[iz][ip][ir];
     
        BUGREPORT;
#ifdef PARALLEL_DAMPING
        FORALL       dtU_0[iz][ip][ir]+=para.gamma*dzzU[iz][ip][ir];
#endif
     
        BUGREPORT;
        /* Time step */
     
        Util_3DSsTimeStep(iter,nx,ny,nz,para.dt,para.mue_w,&lamda, 
                          U_2,U_0,U_1,U_2,dtU_0,dtU_1,dtU_2);
     
        /* We treat viscosity implicitely */
        BUGREPORT;	
        /*vnauplplr_3d(&data,&para,U_2,res,ubdra,ubdrb,
          hval,gval,rcor,rbdcndu,lamda,TRUE);*/
        BUGREPORT;
    
    
#ifndef VSTATIC
     
        /*********************************************/
        /*      Electron  Velocity - Equation                 */
        /*********************************************/
        BUGREPORT;
        Util_3DArakawaNl(dtV_0,f_0,V_0,nlval,nz,nx,ny);
        FORALL dtV_0[iz][ip][ir]*=b_0[iz][ip][ir];
        FORALL dtV_0[iz][ip][ir]-=  MI/ME*(dzN[iz][ip][ir]-dzF[iz][ip][ir]);

        /* selfadvection */
        FORALL dtV_0[iz][ip][ir]+= - 2.*dzV[iz][ip][ir]*V_0[iz][ip][ir];
  

     /* Curvature */
        
        FORALL     dtV_0[iz][ip][ir] += res2[iz][ip][ir];

      /* Neutral collision Term  */
        FORALL dtV_0[iz][ip][ir]+=  -para.delta*V_0[iz][ip][ir];   
     
        /* Coloumb collision Term */
     
        FORALL dtV_0[iz][ip][ir]+= -nu*tprofile[ir]*exp_n_0[iz][ip][ir]*(V_0[iz][ip][ir]-U_0[iz][ip][ir]) ;
     
        BUGREPORT;
#ifdef PARALLEL_DAMPING
        FORALL       dtV_0[iz][ip][ir]+=para.gamma*dzzV[iz][ip][ir];
#endif 
     
        BUGREPORT;
        /* Time step */
     
        Util_3DSsTimeStep(iter,nx,ny,nz,para.dt,0.1*para.mue_w,&lamda, 
                          V_2,V_0,V_1,V_2,dtV_0,dtV_1,dtV_2);
     
        /* We treat viscosity implicitely */
     
        /*       vnauplplr_3d(&data,&para,V_2,res,vbdra,vbdrb,
                 hval,gval,rcor,rbdcndv,lamda,TRUE);        */
#endif
    
#ifdef NEVER
        /********************************/
        /* Temperature - Equation       */ 
        /********************************/
        BUGREPORT;
         
        Util_3DArakawaNl(dtt_0,f_0,t_0,nlval,nz,nx,ny); 
        FORALL     dtt_0[iz][ip][ir] = b_0[iz][ip][ir]*dtt_0[iz][ip][ir];
         
        Util_3DArakawaNl(tres,b_0,f_0,nlval,nz,nx,ny);
        FORALL     dtt_0[iz][ip][ir] -= (2./3.)*tres[iz][ip][ir];
             
        Util_3DArakawaNl(tres,b_0,exp_t_0,nlval,nz,nx,ny);
        FORALL     dtt_0[iz][ip][ir] += (7./3.)*tres[iz][ip][ir];
         
        Util_3DArakawaNl(tres,b_0,n_0,nlval,nz,nx,ny);
        FORALL     dtt_0[iz][ip][ir] += (2./3.)*exp_t_0[iz][ip][ir]*tres[iz][ip][ir];
         
        /* Parallel Terms  */


        /*  (dx ln(T))**2 + (dy ln(T))**2 */
        FORALL 
        {
            at  = 0.5*hval[ir]*(t_0[iz][ip  ][ir+1] - t_0[iz][ip  ][ir-1]);
            bt  = 0.5*vval[ir]*(t_0[iz][ip+1][ir  ] - t_0[iz][ip-1][ir  ]);
            dtt_0[iz][ip][ir] += para.mue_t*(at*at + bt*bt);
        }
             
        Util_3DSsTimeStep(iter,nx,ny,para.dt,para.mue_t,&lamda,t,t_0,t_1,t_2,dtt_0,dtt_1,dtt_2); 
         
        vnauplplr_3d(&data,&para,t_2,t,tres2,tbdra,tbdrb,hval,gval,rcor,rbdcndt,lamda,TRUE); 
     
 
     
        tmp = t_2;     t_2 = t_1;     t_1 = t_0;     t_0  = tmp;
        tmp = dtt_2; dtt_2 = dtt_1; dtt_1 = dtt_0; dtt_0 = tmp;      
#endif 
        
    
        tmp = w_2;     w_2 = w_1;     w_1 = w_0;     w_0  = tmp;
        tmp = dtw_2; dtw_2 = dtw_1; dtw_1 = dtw_0; dtw_0 = tmp;
     
        tmp = n_2;     n_2 = n_1;     n_1 = n_0;     n_0  =tmp;
        tmp = dtn_2; dtn_2 = dtn_1; dtn_1 = dtn_0; dtn_0 = tmp;
     
        tmp = U_2;     U_2 = U_1;     U_1 = U_0;     U_0 = tmp;
        tmp = dtU_2; dtU_2 = dtU_1; dtU_1 = dtU_0; dtU_0 = tmp;
     
        tmp = V_2;     V_2 = V_1;     V_1 = V_0;     V_0 = tmp;
        tmp = dtV_2; dtV_2 = dtV_1; dtV_1 = dtV_0; dtV_0 = tmp;
              
     /*********************************************************************/
      /* Boundaries                                                        */       
      /*********************************************************************/
   
      BUGREPORT;
      FORYZ_BD  nbdrb[iz][ip] = n_0[iz][ip][nx-1];
      FORYZ_BD  nbdra[iz][ip] = n_0[iz][ip][0];
      Util_3DFullBd(n_0,&data,nx,ny,nz,nbdra,nbdrb,nbdrup,nbdrlow,hval,para.dz,rbdcndn,zbndcndn);
      FORALL_BD exp_n_0[iz][ip][ir] = exp(n_0[iz][ip][ir]) ;
      FORYZ_BD  enbdrb[iz][ip] = exp(nbdrb[iz][ip]);
      FORYZ_BD  enbdra[iz][ip] = exp(nbdra[iz][ip]);
     
      Util_3DFullBd(U_0,&data,nx,ny,nz,ubdra,ubdrb,ubdrup,ubdrlow,hval,para.dz,rbdcndu,zbndcndu);
  
      /*
        Calculate Potential and vorticity from 
        density n and GlobalVorticity  W = (Dxx + Dyy) phi + \nabla n \nabla phi 
        
        Iterative solution of equation 
        W = (Dxx + Dyy) phi + \nabla n nabla phi 
      */
      
      FORALL dxn[iz][ip][ir] = DX(n_0);
      FORALL dyn[iz][ip][ir] = DY(n_0);     
      
      for(i=0;i<3;i++)
      {
          FORALL  omega[iz][ip][ir] = w_0[iz][ip][ir]- dxn[iz][ip][ir]*dxf[iz][ip][ir]-dyn[iz][ip][ir]*dyf[iz][ip][ir];
          Util_3DFullBd(omega,&data,nx,ny,nz,ombdra,ombdrb,wbdrup,wbdrlow,hval,para.dz,rbdcndw,zbndcndw);
          HHSOLVER(&data,&para,f_0,omega,fbdra,fbdrb,hval,gval,rcor,rbdcndf,0.,FALSE);
          /*   for(iz=-1;iz<=nz;iz++) for(ip=-1;ip<=ny;ip++) fbdrb[iz][ip] = (1.-para.dt)*fbdrb[iz][ip]+para.dt*f_0[iz][ip][nx-1];   */
          Util_3DFullBd(f_0,&data,nx,ny,nz,fbdra,fbdrb,fbdrup,fbdrlow,hval,para.dz,rbdcndf,zbndcndf);
          FORALL dxf[iz][ip][ir] = DX(f_0);
          FORALL dyf[iz][ip][ir] = DY(f_0);
      }  
      /* End of iterative solve */ 
      
      /* calculate boundary for w  = om + d_r phi  d_r n */
      for(iz=-1;iz<nz+1;iz++) 
          for(ip=-1;ip<ny+1;ip++) 
              wbdra[iz][ip] =   ombdra[iz][ip] /*+ hs[0]* hs[0]*(f_0[iz][ip][0] - f_0[iz][ip][-1])*(n_0[iz][ip][0] - n_0[iz][ip][-1])*/ ;
       /* calculate boundary for w  = om + d_r phi  d_r n */
      for(iz=-1;iz<nz+1;iz++) 
          for(ip=-1;ip<ny+1;ip++) 
              wbdrb[iz][ip] =   ombdrb[iz][ip] /*+ hs[nx-1]* hs[nx-1]*(f_0[iz][ip][nx] - f_0[iz][ip][nx-1])*(n_0[iz][ip][nx] - n_0[iz][ip][nx-1])*/ ;
          
      Util_3DFullBd(w_0,&data,nx,ny,nz,wbdra,wbdrb,wbdrup,wbdrlow,hval,para.dz,rbdcndw,zbndcndw);                 
      /* Parallel derivatives */

      
    
      HWG_DParallelSOL(IS_POTENTIAL,dzF,dzzF,f_0,f_0,n_0,V_0,dzV,&data,&para);
      HWG_DParallelSOL(IS_DENSITY,dzN,dzzN,n_0,f_0,n_0,V_0,dzV,&data,&para);
      HWG_DParallelSOL(IS_IONVELOCITY,dzU,dzzU,U_0,f_0,n_0,V_0,dzV,&data,&para);
          
      /* Calculate Temperature  */
      FORALL_BD exp_t_0[iz][ip][ir]  = exp(t_0[iz][ip][ir]);
        

        /*
          V       - Equation 
        */ 
        /* We neclect electron mass and em-effects a static relationship and
           not a dynamical equation,
                 
        */
    
        FORALL_BD res[iz][ip][ir] = nu*exp_n_0[iz][ip][ir];
<<<<<<< .mine
        FORALL_BD	 V_0[iz][ip][ir]  =  1./(delta + res[iz][ip][ir])
            *(MI/ME*(dzF[iz][ip][ir]-dzN[iz][ip][ir]) +res[iz][ip][ir]*U_0[iz][ip][ir]);
=======
        FORALL_BD vstatic[iz][ip][ir]  =  1./(delta + res[iz][ip][ir]) *(MI/ME*(dzF[iz][ip][ir]-dzN[iz][ip][ir]) +res[iz][ip][ir]*U_0[iz][ip][ir]);
#ifdef VSTATIC
      FORALL_BD	V_0[iz][ip][ir]  = vstatic[iz][ip][ir];
>>>>>>> .r60
#endif      
              
        Util_3DFullBd(V_0,&data,nx,ny,nz,vbdra,vbdrb,vbdrup,vbdrlow,hval,para.dz,rbdcndv,zbndcndv); 
        HWG_DParallelSOL(IS_ELECTRONVELOCITY,dzV,dzzV,V_0,f_0,n_0,V_0,dzV,&data,&para);
    
        para.time +=para.dt;  
        if(iter < 200)  iter++;	
    
        /*
          Write each timestep 
        */
     
/*      FUtils_Write3dSpace(w_0,"Vorticity",n_0,"Density",f_0,"Potential","ptest",iter,&data,&para,TRUE);   */
/*      FUtils_Write3dSpace(V_0,"V",U_0,"U",dtw_1,"DTW","ptest",iter,&data,&para,FALSE);   */
/*      FUtils_Write3dSpace(dtn_1,"DTN",NULL,NULL,NULL,NULL,"ptest",iter,&data,&para,FALSE);   */

        BUGREPORT;
    
    
        if(para.time>= write_time-0.5*para.dt)
        {
            ot++;
            write_time+= para.out_time;
           
            /* 
               Check the CFL condition 
            */
           
            BUGREPORT;
           
            if(Util_3DCheckCFL(f_0,vr,vp,&data,&para,hval,vval,&cflr,&cflp) == -1)
            {
                BUGREPORT;
                FUtils_Write3dSpace(omega,"Vorticity",n_0,"Density",f_0,"Potential","cfl",0,&data,&para,TRUE);
                FUtils_Write3dSpace(f_0,"Potential",vp,"VP",vr,"VR","cfl",0,&data,&para,FALSE);
                fprintf(stderr,"CFL violation. Terminating.\n");
                exit(-1);
            }
           
            BUGREPORT;
           
   
           
            Util_3DIntegral(exp_n_0,1,norm,nx,ny,nz,res[1],&total_density);
           
            /* Calculate Density Fluctuations */
           
            HW3D_Background(exp_n_0,nbg,nx,ny,nz);
            FORALL dyn[iz][ip][ir] = exp_n_0[iz][ip][ir]-nbg[iz][ir];
           
            /* Calculate Potential Fluctuations */
           
            HW3D_Background(f_0,nbg,nx,ny,nz);
            FORALL dtn_0[iz][ip][ir] = f_0[iz][ip][ir]-nbg[iz][ir];	     
            BUGREPORT;
           
           
            /* Calculate Transport */
            FORALL dtU_0[iz][ip][ir] = exp_n_0[iz][ip][ir]*vr[iz][ip][ir];
            Util_3DIntegral(dtU_0,1,norm,nx,ny,nz,res[1],&gamma_n);	     
           
            /* and net transport */

            HW3D_PolAverage(dtU_0,nbg,nx,ny,nz);
           
            /* Calculate Energy */
           
            FORALL dtw_0[iz][ip][ir] = exp_n_0[iz][ip][ir]*(-f_0[iz][ip][ir]*omega[iz][ip][ir]
                                                        +1/MI*V_0[iz][ip][ir]*V_0[iz][ip][ir]
                                                        +U_0[iz][ip][ir]*U_0[iz][ip][ir]);
           


            Util_3DIntegral(omega,2,norm,nx,ny,nz,res[1],&para.vorticity);
            Util_3DIntegral(dtw_0,1,norm,nx,ny,nz,res[1],&para.energy);
            Util_3DIntegral(omega,1,norm,nx,ny,nz,res[1],&charge);
           
            /* Calculate Transport To Limiter Plates */
           
	 
            translim=0.;
            SOL_FORALL_XY_AND_BD  translim += -(para.dx*para.dy*para.dz)*(exp_n_0[0][j][k]*V_0[0][j][k]
                                                                          -exp_n_0[nz-1][j][k]*V_0[nz-1][j][k]);


            /* Calculate Transport into Limiter region */
           
           
            gamma_out = 0.;
            FORALL_ZY gamma_out +=	(para.dx*para.dy*para.dz)*dtU_0[iz][ip][nsols];

            if(FIRST)
            {
                output = fopen(data.erhname,"w");
                fprintf(output,para.codedesc);
                fprintf(output,"# nx %d \t ny %d\t nz %d\n",nx,ny,nz);
                fprintf(output,"# Background density level exp(-kappa_n):\t kappa_n = %g\n",para.kappan);
                fprintf(output,"# Viscosities:\t\t mue_w    = %g\t mue_n = %g\n",para.mue_w,para.mue_n);
                fprintf(output,"# Artificial parallel damping:                     gamma   = %g\n",para.gamma);
                fprintf(output,"# Ion-Neutral collisions:                             sigma     = %g\n",para.sigma);
                fprintf(output,"# Electron-Neutral collisions:                      delta      = %g\n",para.delta);
                fprintf(output,"# Electron-Ion collisions:                            nu          = %g\n",para.nu);
                fprintf(output,"# Density source:                                       source    = %g\n",para.source);
                fprintf(output,"# Ion to electron mass:                                         %g/%g \n",MI,ME);
                fprintf(output,"# \n");
                fprintf(output,"# 1\t 2\t 3\t 4\t 5\t 6\t 7\t 8\t 9\t 10\n");
                fprintf(output,"# t\t E\t W\t G(n)\t CFL_x\t CFL_y\t k_eff\t N_tot\t Charge\t Translim \t Gamm_out\n");
                fclose(output);
                FIRST = FALSE;
            }
	 
            output = fopen(data.erhname,"a");
            fprintf(output,"%g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n",
                    para.time,para.energy,para.vorticity,gamma_n,cflr,cflp,sqrt(para.vorticity/para.energy),total_density,charge,translim,gamma_out);
            fclose(output);
            BUGREPORT;
           

            if(ot >= para.otmult)
            {
                ot = 0;
                data.number++;
                BUGREPORT;
                FORALL dyn[iz][ip][ir]=nbg[iz][ir];
                FORALL res[iz][ip][ir]=V_0[iz][ip][ir]-U_0[iz][ip][ir];
               
                FUtils_Write3dSpace(omega,"Vorticity",n_0,"logN",f_0,"Potential",
                                    data.name_out,data.number,&data,&para,TRUE);  
                FUtils_Write3dSpace(V_0,"V_0",dtn_0,"f",NULL,NULL,
                                    data.name_out,data.number,&data,&para,FALSE);
                FUtils_Write3dSpace(U_0,"U_0",dtU_0,"Flux",exp_n_0,"Density",
                                    data.name_out,data.number,&data,&para,FALSE);
                FUtils_Write3dSpace(res,"Current",dzV,"dzV",dzN,"dzN",
                                    data.name_out,data.number,&data,&para,FALSE);
                FUtils_Write3dSpace(w_0,"w_0",dxn,"dxn",dxf,"dxf",
                                    data.name_out,data.number,&data,&para,FALSE);

            }
        } 
        
    }
    
    exit(0);
}

  


/*********************************************************************************/
/* Calculates background, as poloidal average!!  
 */
void  HW3D_Background(double ***n_0,double **nbg,int nx,int ny,int nz)
{
    register int i,j,k;
    double fac= 1./(double)(ny);
    double val = 0.;


    for(i=0;i<nz;i++) 
        for(k=-1;k<=nx;k++)
        {
            val = 0.;
            for(j=0;j<ny;j++) val += n_0[i][j][k];
            nbg[i][k] = val*fac;
        }
}
/*********************************************************************************/
/*
  Calculates poloidal average
*/

void  HW3D_PolAverage(double ***s,double **t,int nx,int ny,int nz)
{
    register int i,j,k;
    double fac= 1./(double)(ny);
    double val = 0.;


    for(i=0;i<nz;i++) 
        for(k=0;k<nx;k++)
        {
            val = 0.;
            for(j=0;j<ny;j++) val += s[i][j][k];
            t[i][k] = val*fac;
        }
}
/*************************************************************************************/
void      HWG_DParallelSOL(int identity,double ***result,double ***result2,
                           double ***val,double ***f,double ***n,double ***V,double ***dzv,
                           HDF_DS *data, PARA *para)
{
    /*   This routine calculates the parallel derivative: 
         and takes care of the parallel SOL boundary condition
       
         As parameters it gets: 
         int identity : kind of field to apply BD conditions right
         result: Stores the result on exit 
         val   : f 
         f     : Potential for sheath condition
         lf    : factor for the linear part of the derivative
    */



    static int FIRST= TRUE;
    static int nx,ny,nz,offx,offy,offz,nsols,nsole;
    static double lf_local,lf_local2,*sol_prof,*per_prof,fsheath;
    register int i,j,k;
    double bdval,tval;
    
    const double mass=1./sqrt(2.*M_PI)*sqrt(MI/ME); 

    if(FIRST)
    {
        nx = data->nx;
        ny = data->ny;
        nz = data->nz;
	
        offx = data->offx;
        offy = data->offy;
        offz = data->offz;
	
        /* Factor for derivative */
	
        lf_local = 1./(12. *para->dz);
        lf_local2 = 1./(12.* para->dz*para->dz);
	
        /* Profile of SOL */
        sol_prof = Util_DVector(nx,offx);
        per_prof = Util_DVector(nx,offx);
        for (k=-1;k<=nx;k++) sol_prof[k] = 1.;
        
        /* SOL starts at x=0 */
        nsols = -(int)(para->xmin/para->dx);

        /* There is a transition to full SOL within 5 grid points */ 
        nsole = nsols+5;
        
        if(nsols> nx) nsols = nx;
        if(nsole>nx) nsole = nx;
        
        if(nsols< 0) nsols = -1;
        if(nsole < 0 ) nsole = -1;
        
        fprintf("nsols = %d, nsole = %d\n",nsols,nsole);

        for (k=-1;k<nsols;k++) 
            sol_prof[k] = 0.;
        
        for (k=nsols;k<nsole;k++) 
            sol_prof[k] = (double)(k-nsols)/(double)(nsole -nsols);     
        
        for(k=nsole;k<=nx;k++) 
            sol_prof[k] = 1.;
        
        for (k=-1;k<=nx;k++) per_prof[k] = 1.-sol_prof[k];
        
        
        /* Potential is relative to sheath potential */
        fsheath = -log(1/mass);        
        FIRST = FALSE;
    }
    
    
    
    /* Initialisation finished */

    if(data->offz == 1)
        FORALL_BD_XY  
        {
            val[nz][j][k] = val[0][j][k];
            val[-1][j][k] = val[nz-1][j][k];
        }
    else
        FORALL_BD_XY  
        {
            val[nz][j][k] = val[0][j][k];
            val[nz+1][j][k] =  val[1][j][k];
            val[-1][j][k] = val[nz-1][j][k];
            val[-2][j][k] =  val[nz-2][j][k];
        }
    

    if(data->offz == 1)
    {
        switch(identity){
            case  IS_ELECTRONVELOCITY:
                SOL_FORALL_XY_AND_BD  
                /* Effective Electron speed is given by the sheath potential */
                {
                    /* Velocity is positive as electrons leave the plasma */
                    /* if f=0 the velocity is 1 */

                    bdval =    2.*exp(-f[nz-1][j][k]);
     
                    val[nz][j][k] = per_prof[k]*(val[0][j][k]) 
                        + sol_prof[k]*( bdval-val[nz-1][j][k]);
                        
                        
                    bdval =    -2.*exp(-f[0][j][k]);
                    val[-1][j][k] =  per_prof[k]*(val[nz-1][j][k]) 
                        + sol_prof[k]* (bdval-val[0][j][k]);
                }
                break;
            case  IS_IONVELOCITY:
                /* Ion speed is unity into the limiter */
                SOL_FORALL_XY_AND_BD
                {
                    bdval= 2.;
                    val[nz][j][k] = per_prof[k]*(val[0][j][k]) 
                        +sol_prof[k]* (bdval-val[nz-1][j][k]);
                                    
                    bdval= -2.;
                    val[-1][j][k] = per_prof[k]*(val[nz-1][j][k] )
                        +sol_prof[k]* (bdval-val[0][j][k]);
                }
                break;
            case IS_DENSITY:
                /* density flux remains constant at the sheath 
                   d(nv) = 0 => v dn = - n dv
                   dln  n = -  dv/v
                */
                SOL_FORALL_XY_AND_BD
                {
                    /*   val[nz][j][k] =  per_prof[k]*val[nz][j][nsols-1]
                         +sol_prof[k]* 0.02/(6.*lf_local)+val[nz-1][j][k];*/

                    /*val[-1][j][k] =  per_prof[k]*val[-1][j][nsols-1] 
                      +sol_prof[k]* -0.02/(6.*lf_local)+val[0][j][k]; */
             
                    val[nz][j][k] =  per_prof[k]*val[nz][j][nsols-1] 
                                       +sol_prof[k]*  -0.05*dzv[nz-1][j][k]/MAX(V[nz-1][j][k],0.1)/(6.*lf_local)+val[nz-1][j][k];

                     val[-1][j][k] =  -0.025*dzv[0][j][k]/MAX(-V[0][j][k],0.1)/(6.*lf_local)+val[0][j][k];
                }  
                break;
            case IS_POTENTIAL:
                /* No gradient in potential, sheeth works via parallel current */
                SOL_FORALL_XY_AND_BD
                {   
                    /*val[nz][j][k] =  per_prof[k]*val[nz][j][nsols-1]
                      +sol_prof[k]* 0.02/(6.*lf_local)+val[nz-1][j][k];*/               

                    /*val[-1][j][k] =  per_prof[k]*val[-1][j][nsols-1] 
                      +sol_prof[k]* -0.02/(6.*lf_local)+val[0][j][k];*/

                    val[nz][j][k] =  per_prof[k]*val[nz][j][nsols-1] 
                                       +sol_prof[k]*  -0.05*dzv[nz-1][j][k]/MAX(V[nz-1][j][k],0.1)/(6.*lf_local)+val[nz-1][j][k];

                     val[-1][j][k] =  -0.025*dzv[0][j][k]/MAX(-V[0][j][k],0.1)/(6.*lf_local)+val[0][j][k];          
                }
                break;
            case IS_TEMPERATURE:
                break;
            default:
                fprintf(stderr,"unidentified field object (UFO) in SOL boundary routine!\n");
        }

   
        /* dz */
        FORALL_BD_XY  
            result[i][j][k]  =  6.*lf_local*(val[i+1][j][k] - val[i-1][j][k]);
        
        /* dzz */
        FORALL_BD_XY  
            result2[i][j][k]  =  12.*lf_local2 *((val[i+1][j][k] + val[i-1][j][k]) -2.* val[i][j][k]);

    }
    else
    {
        switch(identity){
            case  IS_ELECTRONVELOCITY:
                SOL_FORALL_XY_AND_BD  
                /* Effective Electron speed is given by the sheath potential */
                {
                    bdval =    2.*exp(-f[nz-1][j][k]);
                    tval = 1.;
        
                    val[nz][j][k] = tval*per_prof[k]*val[nz][j][nsols] 
                        + sol_prof[k]*( bdval-val[nz-1][j][k]);
                    val[nz+1][j][k] =  tval*per_prof[k]*val[nz+1][j][nsols] 
                        +sol_prof[k]*  (bdval-val[nz-2][j][k]);
                    

                    bdval =    -2.*exp(-f[0][j][k]);

                    val[-1][j][k] =  tval*per_prof[k]*val[-1][j][nsols] 
                        + sol_prof[k]* (bdval-val[0][j][k]);
                    val[-2][j][k] =  tval*per_prof[k]*val[-2][j][nsols] 
                        +  sol_prof[k]* (bdval-val[1][j][k]);
                    

                    /* 2* the boundary value */
/*                    bdval =    2.*mass*exp(-(f[nz-1][j][k]+fsheath));
             
val[nz][j][k] = per_prof[k]* (tval-val[nz-1][j][k]) + sol_prof[k]*(bdval-val[nz-1][j][k]);
val[nz+1][j][k] =  per_prof[k]*(tval-val[nz-2][j][k]) +sol_prof[k]*(bdval-val[nz-2][j][k]);
                    
bdval =    -mass*exp(-(f[0][j][k]+fsheath));
val[-1][j][k] =  per_prof[k]*(tval-val[0][j][k]) + sol_prof[k]* (bdval-val[0][j][k]);
val[-2][j][k] =  per_prof[k]*(tval-val[-1][j][k]) +  sol_prof[k]*(bdval-val[1][j][k]);
*/

                }
                break;
            case  IS_IONVELOCITY:
                /* Ion speed is unity into the limiter */
                SOL_FORALL_XY_AND_BD
                {
                    bdval=2.;
                    tval =1.;
                    
                    val[nz][j][k] =  tval*per_prof[k]*val[nz][j][nsols] 
                        +sol_prof[k]* (bdval-val[nz-1][j][k]);
                    val[nz+1][j][k] =  tval*per_prof[k]*val[nz+1][j][nsols]
                        + sol_prof[k]* (2.*bdval-val[nz-2][j][k]);
                    
                    bdval=-2.;
                    val[-1][j][k] = tval*per_prof[k]*val[-1][j][nsols] 
                        +sol_prof[k]* (2.*bdval-val[0][j][k]);
                    val[-2][j][k] = tval*per_prof[k]*val[-2][j][nsols] 
                        + sol_prof[k]* (2.*bdval-val[1][j][k]);


                    /* 2* the boundary value */
                    /*                    bdval=2.; */
/*                     tval = 0.; */
/*                     val[nz][j][k] =  per_prof[k]*(tval-val[nz-1][j][k]) +sol_prof[k]* (bdval-val[nz-1][j][k]); */
/*                     val[nz+1][j][k] =  per_prof[k]*(tval-val[nz-2][j][k]) + sol_prof[k]* (bdval-val[nz-2][j][k]); */
                    
                    
                    /* 2* the boundary value */
                    /*     bdval=-2.; */
/*                     tval = 0.; */
/*                     val[-1][j][k] = per_prof[k]*(tval-val[0][j][k]) +sol_prof[k]* (bdval-val[0][j][k]); */
/*                     val[-2][j][k] = per_prof[k]*(tval-val[1][j][k]) + sol_prof[k]* (bdval-val[1][j][k]); */
                }
                break;
            case IS_DENSITY:
                /* density remains constant at the sheath */
                SOL_FORALL_XY_AND_BD
                {
                    val[nz][j][k] =  per_prof[k]*val[nz][j][nsols] 
                        +sol_prof[k]* val[nz-1][j][k];
                    val[nz+1][j][k] =  per_prof[k]*val[nz+1][j][nsols] 
                        +sol_prof[k]*  val[nz-2][j][k];
                    
                    val[-1][j][k] =  per_prof[k]*val[-1][j][nsols] 
                        +sol_prof[k]* val[0][j][k];
                    val[-2][j][k] =  per_prof[k]*val[-2][j][nsols] 
                        +sol_prof[k]* val[1][j][k];
                }           /* Density goes to background value */
                break;
            case IS_POTENTIAL:
                /* No gradient in potential, sheeth works via parallel current */
                SOL_FORALL_XY_AND_BD
                {
                    val[nz][j][k] =  per_prof[k]*val[nz][j][nsols-1] 
                        +sol_prof[k]*   (val[nz-1][j][k]);
                    val[nz+1][j][k] =  per_prof[k]*val[nz+1][j][nsols-1]
                        +sol_prof[k]*  (val[nz-2][j][k]);

                    val[-1][j][k] =  per_prof[k]*val[-1][j][nsols-1] 
                        +sol_prof[k]*  (val[0][j][k]);
                    val[-2][j][k] =  per_prof[k]*val[-2][j][nsols-1]
                        +sol_prof[k]* (val[1][j][k]);
                }
                break;
            case IS_TEMPERATURE:
                break;
            default:
                fprintf(stderr,"unidentified field object (UFO) in SOL boundary routine!\n");
        }

        /* dz */
        FORALL_BD_XY  
            result[i][j][k]  =  lf_local*(-val[i+2][j][k] + 8.*(val[i+1][j][k] - val[i-1][j][k]) + val[i-2][j][k]);
        
        /* dzz */
        FORALL_BD_XY  
            result2[i][j][k]  =  lf_local2 *(-val[i+2][j][k] + 16.*(val[i+1][j][k] + val[i-1][j][k]) -30.* val[i][j][k] - val[i-2][j][k]);
    }    
    
}


/****************************************************************************/
void HWG_Laplace(double ***w,double ***f,double*hval,double*vval,
                 double *edr,int nx,int ny,int nz,PARA *p,int ispolar)
{
    /* This function calculates the laplace of a scalar field, given the 
    
   
    w = Lap f =  d_rr f + 1/r d_r f + 1./r2 d_ff f 

    The laplacian is written out in field w all other fields are input

    w    : on output contains lap f
    f    : field of which laplce is to be calculated
    hval : 1/dx(x) according to numerical grid
    vval : 1/dy(x) according to numerical grid
    edr  : 1/r
    nx,ny,nz: Dimensions of fields
   
    */
    register  int i,j,k;
    double ky;
  
    /* Using FFT */
/* High accuracy 10e-11 */  
/* This is an exact inverse to the Poisson solver! */
  
    BUGREPORT;
    for(i=0;i<nz;i++) 
    {    
        for(j=-1;j<=ny;j++)  memcpy((void *)(w[i][j]-1),(void *)(f[i][j]-1),(nx+2)*sizeof(double));

        Fft_1d_2d_f(w[i],nx,ny);

        for(j=-1;j<=ny;j++)
        {
            ky = (j/2)*p->dky;  
            ky*=-ky;
            for(k=0;k<nx;k++) w[i][j][k] *= ky;
        }
     
        Fft_1d_2d_b(w[i],nx,ny);
        for(j=-1;j<=ny;j++) for(k=0;k<nx;k++)  w[i][j][k] *= edr[k]*edr[k];
        for(j=-1;j<=ny;j++) for(k=0;k<nx;k++)  w[i][j][k] += hval[k]*hval[k]*(f[i][j][k+1] -2.*f[i][j][k]  +f[i][j][k-1]);
        BUGREPORT;

        if(ispolar)   for(j=-1;j<=ny;j++) for(k=0;k<nx;k++)  w[i][j][k] += 0.5*edr[k]*hval[k]*(f[i][j][k+1] - f[i][j][k-1]);
        BUGREPORT;
    }
}




