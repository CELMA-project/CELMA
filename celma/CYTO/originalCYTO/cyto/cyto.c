/**********************************************************************/



// Defines
#ifndef TE
#undef TE // Use electron temperature equation
#endif

#define LAXF 2. //50.


#undef VSTATIC // Using Vstatic dereases necessary dt in first steps
#undef QSTATIC // Test simple electron heat flux
#define ME 0.01 //(1./1836.2)

#define DX(F)  0.5*hval[iz][ir]*( (F)[iz][ip][ir+1] - (F)[iz][ip][ir-1])
#define DY(F)  0.5*vval[iz][ir]*( (F)[iz][ip+1][ir] - (F)[iz][ip-1][ir])


#undef NONLINEAR
#ifdef NONLINEAR
#define CONVECT Util_3DArakawaNl
#else
#define CONVECT Util_3DSet2Zero
#endif

/* These Flags are for debugging */
#undef DEBUG
#undef CNTFILES
#undef COMM_STATUS
#define CHECK_STRUCTURE  Util_PrintStructures(d,p,__LINE__)
// #define CHECK_STRUCTURE


#define FORALL for(iz=0;iz<nz;iz++) for(ip=0;ip<ny;ip++) for(ir=0;ir<nx;ir++)
#define FORALL_BD for(iz=-d->offz;iz<(nz+d->offz);iz++) for(ip=-d->offy;ip<(ny+d->offy);ip++) for(ir=-d->offx;ir<(nx+d->offx);ir++)

#define FORYX    for(ip=0;ip<ny;ip++)    for(ir=0;ir<nx;ir++)
#define FORYX_BD for(ip=-1;ip<ny+1;ip++) for(ir=-1;ir<nx+1;ir++)

#define FORZY for(iz=0;iz<nz;iz++) for(ip=0;ip<ny;ip++)
#define FORZY_BD for(iz=-d->offz;iz<nz+d->offz;iz++)  for(ip=-d->offy;ip<(ny+d->offy);ip++)

#define FORZX  for(iz=0;iz<nz;iz++)  for(ir=0;ir<nx;ir++)
#define FORZX_BD  for(iz=-d->offz;iz<(nz+d->offz);iz++)  for(ir=-d->offx;ir<(nx+d->offx);ir++)

#define FORX    for(ir=0;ir<nx;ir++)
#define FORX_BD for(ir=-d->offx;ir<(nx+d->offx);ir++)

#define FORY    for(ip=0;ip<ny;ip++)
#define FORY_BD for(ip=-d->offy;ip<(ny+d->offy);ip++)

#define FORZ    for(iz=0;iz<nz;iz++)
#define FORZ_BD for(iz=-d->offz;iz<(nz+d->offz);iz++)

#define FORALL_BD_XY for(iz= 0;iz<nz;iz++) for(ip=-1;ip<(ny+1);ip++) for(ir=-1;ir<(nx+1);ir++)


/***********************************************************************/


#include <mpi.h>
#include <utilities.h>
#include <par_helper.h>

/***********************************************************************/
/*

  Local Functions

*/
void Cyto_Laplace(double ***w,double ***f,double **hval,double **vval,
                  double *edr,int nx,int ny,int nz,PARA *p,int ispolar);

void Cyto_DParallelSOL(int identity,double ***result,double ***result2,
                       double ***val,double ***f,double ***n,double ***exp_t,
                       double ***v,double ***dzv,
                       double *zhval, double *zgval,
                       HDF_DS *data, PARA *para);

void Cyto_CalcParameters(PARA *p, HDF_DS *d) ;


void Map_Noise(HDF_DS *data, PARA *p);

void spectral_interpol(double *ret_val,double *xpos,double *ypos,int nop,double **phi,int nx, int ny, double *Sx, double *Sy);



/***********************************************************************/


int main(int argc,char **argv)
{

  int
      W_UION   = TRUE,
      W_DZN    = TRUE,
      W_H0     = TRUE,
      W_F0     = TRUE,
      W_N0     = TRUE,
      W_lnN0   = FALSE,
      W_V0     = TRUE,
      W_W0     = TRUE,
      W_DZF    = TRUE,
      W_VSTATIC= FALSE,
      W_J0     = TRUE,
      W_PFLUX  = TRUE,
      W_DZVN   = TRUE,
      W_DZUN   = TRUE,
      W_NUVMU  = FALSE,
      W_TE0    = TRUE,
      W_lnT0   = TRUE,
      W_Q0     = TRUE,
      W_DZQ     = TRUE;

    const char code_desc[] ={
     "# The CYTO code solves the global, electrostatic fluid equations for a plasma in a cylinder:\n"\
     "#\n"                                                              \
     "#  N = ln n \n"                                                   \
     "#  H = omega + Nabla N nabla phi\n"                                 \
     "#  vorticity = omega = laplace phi\n"                             \
     "#\n"                                                              \
     "#  Dt N = - dz V -  V dz N  \n"                                   \
     "#\n"                                                              \
     "#  dt H = 1/Z*(dz (U-V) + (U-V) dz N  + nabla N ( mu nabla w - nu nabla f)" \
     "    +{phi, omega}  + nabla n { phi, nabla phi} + nabla dt n . nabla phi   " \
     "   + mu nabla^2 H -nu w \n"                                       \
     "#\n"                                                              \
     "#  V = U - tau  D_par (phi - N) + alpha\n"                                \
     "#\n"                                                              \
     "#  D_t U = D_par N - 2 U dz U \n"                                 \
         "#\n"};


    // variables with "-" are not read/written, change description if new variable is used

    const char para_desc[] = {
        "gamma: -;"                                                     \
        "sigma: ion-neutral collisions;"                                \
        "nu: electron-ion collisions;"                                  \
        "beta: -;"                                                      \
        "shear: background shear flow;"                                 \
        "alpha: background electron drift in cs;"
        "delta: electron-neutral collisions;"                           \
        "betahat: -;"                                                   \
        "adrhos: -;"                                                    \
        "mue_w: viscosity;"                                             \
        "mue_n: diffusion;"                                             \
        "mue_t: temperature;"                                                     \
        "kappan: width of gaussian density source (length);"            \
        "kappat: -;"                                                    \
        "r0: -;"                                                        \
        "muehat: -;"                                                    \
        "bprof: -;"                                                     \
        "tprof: temperature drive rate;"                                                     \
        "nprof: density source drive rate;"                             \
        "phiprof: -;"                                                   \
        "hm_nl: parallel advection (0/1);"                                                     \
        "exb_ll: -;"                                                    \
        "limiter: parallel damping (artificial);"                       \
        "source: background source to keep density from falling to zero;"   \
        "k0: -;"                                                        \
        "dt_pol: -;"                                                    \
        "uvortex: -;"                                                   \
        "radius: -;"};

    //string with primary geometry parameter descriptions

    const char geom_desc[] ={
        "R0: -;"                                        \
        "a: radius of plasma column      in [m];"       \
        "lpar: parallel length in [m];"                 \
        "ln: density gradient length scale  in [m];"    \
        "lTe: -;"                                       \
        "lTi: -;" };


    // These structures define most of the variables, they are defined in include/intergrator.h
    HDF_DS
        data = {0}, *d; // compiler will warn about missing braces, but this is c99 standard
    PARA
        para = {0}, *p;

    d = &data;
    p = &para;


    void (*ptrtoCalcParameters)(PARA *,HDF_DS *) = NULL;

    FILE
        *output;

    double
        ***res      = NULL,
        ***tmp      = NULL,
        ***f_0      = NULL,
        ***vr       = NULL,
        ***vp       = NULL,
        ***dsource  = NULL;

    double
        **fbdrup    = NULL,
        **fbdrlow   = NULL,
        **fbdr[2]   = {NULL,NULL};

    double
        ***w        = NULL,
        ***w_0      = NULL,
        ***w_1      = NULL,
        ***w_2      = NULL,
        ***dtw_0    = NULL,
        ***dtw_1    = NULL,
        ***dtw_2    = NULL,
        **wbdr[2]   = {NULL,NULL},
        **ombdr[2]  = {NULL,NULL},
        **wbdrlow   = NULL,
        **wbdrup    = NULL;

    double
        ***n        = NULL,
        ***n_0      = NULL,
        ***n_1      = NULL,
        ***n_2      = NULL,
        ***dtn_0    = NULL,
        ***dtn_1    = NULL,
        ***dtn_2    = NULL,
        ***exp_n_0  = NULL,
        ***exp_t_0  = NULL,
        ***target   = NULL,
        **nbdr[2]   = {NULL,NULL},
        **nbdrlow   = NULL,
        **nbdrup    = NULL;

    // Temperature
    double
        ***t_0             = NULL,
        ***t_1             = NULL,
        ***t_2             = NULL,
        ***dtt_0           = NULL,
        ***dtt_1           = NULL,
        ***dtt_2           = NULL,
        **tbdr[2]          = {NULL,NULL},
        **tbdz[2]          = {NULL,NULL};

        // Electron heatflux
    double
        ***q_0             = NULL,
        ***q_1             = NULL,
        ***q_2             = NULL,
        ***dtq_0           = NULL,
        ***dtq_1           = NULL,
        ***dtq_2           = NULL,
        **qbdr[2]          = {NULL,NULL},
        **qbdz[2]          = {NULL,NULL};


        // Electron velocity
    double
        ***V_0      = NULL,
        ***V_1      = NULL,
        ***V_2      = NULL,
        ***dtV_0    = NULL,
        ***dtV_1    = NULL,
        ***dtV_2    = NULL,
        **vbdr[2]   = {NULL,NULL},
        **vbdrlow   = NULL,
        **vbdrup    = NULL;

    // Ion velocity
    double
        ***U_0      = NULL,
        ***U_1      = NULL,
        ***U_2      = NULL,
        ***dtU_0    = NULL,
        ***dtU_1    = NULL,
        ***dtU_2    = NULL,
        **ubdr[2]   = {NULL,NULL},
        **ubdrlow   = NULL,
        **ubdrup    = NULL;

// LAX
        double
            ***laxt = NULL,
            ***laxw = NULL,
            ***laxn = NULL,
            ***laxq = NULL,
            ***laxv = NULL;


// Neutrals Slow species
        double
                ***H_slow_0        = NULL,
            ***H_slow_1      = NULL,
        ***H_slow_2      = NULL,
        ***dtH_slow_0    = NULL,
        ***dtH_slow_1    = NULL,
        ***dtH_slow_2    = NULL,
        **      H_slow_bdr[2]   = {NULL,NULL},
        **H_slow_bdrlow   = NULL,
        **H_slow_bdrup    = NULL;

// Neutrals Fast species
        double
                ***H_fast_0        = NULL,
            ***H_fast_1      = NULL,
        ***H_fast_2      = NULL,
        ***dtH_fast_0    = NULL,
        ***dtH_fast_1    = NULL,
        ***dtH_fast_2    = NULL,
        **      H_fast_bdr[2]   = {NULL,NULL},
        **H_fast_bdrlow   = NULL,
        **H_fast_bdrup    = NULL;


    double
        *pol_pot    = NULL,
        *pol_dens   = NULL,
        *u_par      = NULL;

    double
        ***omega    = NULL,
        ***vstatic  = NULL,
        ***dyf      = NULL,
        ***dxf      = NULL,
        ***dxn      = NULL,
        ***dyn      = NULL,
        ***dxt      = NULL,
        ***dyt      = NULL,
        ***dtdxf    = NULL,
        ***dtdyf    = NULL;

    double
        ***dzU      = NULL,
        ***dzzU     = NULL,
        ***dzV      = NULL,
        ***dzzV     = NULL;

    double
        ***VN       = NULL,
        ***dzVN     = NULL,
        ***dzzVN    = NULL;

    double
        ***UN       = NULL,
        ***dzUN     = NULL,
        ***dzzUN    = NULL;


    double
        ***dzF      = NULL,
        ***dzzF     = NULL,
        ***dzN      = NULL,
        ***dzzN     = NULL,
        ***dzt      = NULL,
        ***dzzt     = NULL,
        ***dzq      = NULL,
        ***dzzq     = NULL;

    double
        ***J_0      = NULL,
        ***dzJ      = NULL,
        ***dzzJ     = NULL,
        ***tfeld    = NULL,
        ***nu_density = NULL,
        ***nu_vmu   = NULL;

    double
        *ts_n       = NULL,
        *ts_f       = NULL,
        *ts_e       = NULL;



    double // radial profiles
        *ppar_prof  = NULL,
        *vpar_prof  = NULL,
        *upar_prof  = NULL,
        *npar_prof  = NULL,
        *fpar_prof  = NULL,
        *wpar_prof  = NULL;

   double // axial profiles
        *pax_prof  = NULL,
        *vax_prof  = NULL,
        *uax_prof  = NULL,
        *nax_prof  = NULL,
        *fax_prof  = NULL,
        *wax_prof  = NULL;

    double
        Kappa             = 1.6,
        Alpha             = 0.71,
        tf                = 2./3.,
        mum1              = 0.,
        aL                = 0.,
        fac_pot           = 0.,
        val         = 0.,
        norm2d      = 0.,
        nu_m        = 0.,
        nu          = 0.,
        delta       = 0.,
        padvection  = 0.; // 0 in GW version VN 15032010

    double
        dz   = 0.,
        edz  = 0.,
        zm   = 0.,
        g    = 0.,
        gm1  = 0.,
        rlen = 0.,
        xlen = 0.;


    double
        **hval      = NULL,
        **vval      = NULL,
        **nlval     = NULL,
        *rcor       = NULL,
        *edrcor     = NULL,
        **norm      = NULL,
        *dynbg      = NULL;


    double
        *zcor     = NULL,
        *zhval    = NULL,
        *zgval    = NULL;


    double
        *par_elec_mom_flux = NULL,
        *par_ion_mom_flux  = NULL;

    double
        totalf   = 0.,
        Jtotal   = 0.,
        Ipol     = 0.,
        G_par    = 0.,
        J_par    = 0.,
        N_in     = 0.,
        mu       = 0.,
        r        = 0.,
        z        = 0.,
        cflr     = 0.,
        cflp     = 0.;

    int
        nx            = 0,
        ny            = 0,
        nz            = 0,
        i             = 0,
        ir            = 0,
        ip            = 0,
        irg           = 0,
        iz            = 0,
        izg           = 0,
        first         = TRUE,
        counter       = 0,
        ot            = 1,
        iter          = 0,
        ispolar       = FALSE,
        WRITE         = TRUE,
        ErhFirstWrite = TRUE;

    char
        mom_name[256]    = "Empty",
        filename[256]    = "Empty",
        profilename[256] = "Empty";

    double
        lamda    = 0.,
        gamma_n  = 0.,
        phase    = 0.,
        min_nt   = 1.,
        lmin_nt  = 0.;
    double
        *par_rey_stress   = NULL,
        *flux_n           = NULL,
        *av_flux_n        = NULL;

    double
        density_source_integral = 0.,
        impulse           = 0.,
        dtimpulse         = 0.,
        nloss             = 0.,
        write_time        = 0.,
        Lz                = 0.;

    int
        rbdcndw  = 0 ,
        rbdcndv  = 0 ,
        rbdcndu  = 0 ,
        rbdcndn  = 0 ,
        rbdcndf  = 0 ,
        rbdcndt  = 0 ,
        rbdcndq  = 0 ,
        zbndcndv = 0 ,
        zbndcndw = 0 ,
        zbndcndu = 0 ,
        zbndcndn = 0 ,
        zbndcndf = 0 ,
        zbndcndt = 0 ,
        zbndcndq = 0 ;

    int
        off = 0,
        m = 0;
    // Variables for collecting communication to larger packages
    double
        reducing[25]      = {0},
        result_vector[25] = {0};

    /************************************************/
    /*                                              */
    /* start of body/read data and command line     */
    /*                                              */
    /************************************************/
    // Initialize data structures (necessary!)
    FUtils_IniStructure(d,p);
    CHECK_STRUCTURE;


    ptrtoCalcParameters = Cyto_CalcParameters;
    BUGREPORT;
    // Initialize MPI

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&d->this_process);
    COMM(fprintf(stderr,"Proc %d initializing...\n",d->this_process););
    MPI_Comm_size(MPI_COMM_WORLD,&d->num_procs);
    COMM(fprintf(stderr,"Procs %d\n",d->num_procs););

    // MPI_Datatypes
    hdfdata_type = create_mpi_hdfds(d,ROOT, MPI_COMM_WORLD,d->this_process);
    para_type = create_mpi_para(p, ROOT, MPI_COMM_WORLD,d->this_process);


    snprintf(p->codedesc,   DEFSTRLEN,"%s",code_desc);
    snprintf(p->desc,       DEFSTRLEN,"%s",para_desc);
    snprintf(p->Prim_Geom,  DEFSTRLEN,"%s",geom_desc );
    snprintf(p->codename,   DEFSTRLEN,"%s","cyto");


    /*
    Default values for parameters, these will be written into the sample.ini file if the code is called with the -H option
    Do not change these values, as they will be used to define a standard run!
    Data are typical values for Mirabelle.

    */

    d->nx                   = 40; /*Number of grid points in x*/
    d->ny                   = 64; /*Number of grid points in y*/
    d->nz                   = 96; /*Number of grid points in z*/
    p->coordsys             = CYLINDRICAL ;

    d->N[2]               = 1; /*! Number of procs in x*/
    d->N[1]               = 1; /*! Number of procs in y*/
    d->N[0]               = 0; /*! Number of procs in z*/


    p->dt                   = 1.e-4; //Delta T
    p->end_time             = 1.000000; //Simulation stops at this time
    p->out_time             = 5e-3; //Time between small outputs
    p->otmult               = 1; //Number of small outputs  before output of fields

    p->Ti                   = 0.030000; /*Central Ion Temperature [eV]*/
    p->Te                   = 2.000000; /*Central Electron Temperature [ev]*/
    p->B0                   = 0.070000; /*Magnetic field strength [T]*/
    p->n0                   = .005;/* central density [10^19 1/m**3]*/
    p->Mi                   = 40; /*Ion mass in [u]*/
    p->Z                    = 1.000000; /*Charge of main Ion species [e]*/
    p->p_n                  = 0.300000; /*Neutral pressure in [pa]*/

    p->a                    = 0.050000; /*radius of plasma column      in [m]*/
    p->lpar                 = 1.500000; /*parallel length in [m]*/
    p->ln                   = 0.0250000; /*density gradient length scale  in [m]*/

    p->part_source          = 0.0; // particle source in [10^19 m-3 s-1]
    p->temp_source          = 0.0; // energy source in [ eV 10^19 m-3 s-1]


    p->limiter              = 0.000000; // parallel damping (artificial)
    p->kappan               = 0.000000;  // width of gaussian density source
    p->nprof                = 0.250000;  // density source
    p->tprof                = 0.250000;  // density source
    p->alpha                = 0.000000;  // parallel electron drift

    d->rank   = 3; //dimensionality
    d->offz   = 1; // Number of ghostpoints
    d->offx   = 1;
    d->offy   = 1;
    d->number = 0;

    d->N[2]               = 1; /*! Number of procs in x*/
    d->N[1]               = 1; /*! Number of procs in y*/
    d->N[0]               = 0; /*! Number of procs in z*/


    MPI_Barrier(MPI_COMM_WORLD);

    /*****************************************************************************/
    /* First we only read information and not the actual fields....              */
    /*****************************************************************************/


    if(ISROOT)
    {
        d->read_data      = FALSE;
        d->ReadAttributes = TRUE;
        d->restart        = FUtils_ReadArguments(argc,argv,d, p,ptrtoCalcParameters);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Broadcast DATA  out to the other processors
    MPI_Bcast(d,1,hdfdata_type,ROOT,MPI_COMM_WORLD);
    CHECK_STRUCTURE;
    /***********************************/
    /* Geometry of processes           */
    /***********************************/
    d->zisperiodic = TRUE; // This is needed to set up periodic processor grid
    PH_3D_Geometry(d);

    COMM(if(ISROOT) {
        fprintf(stderr,"%ld %ld %ld \n",d->dims[2],d->dims[1],d->dims[0]);
        fprintf(stderr,"%ld %ld %ld \n",d->nx,d->ny,d->nz);
        fprintf(stderr,"%ld %ld %ld \n",d->offx,d->offy,d->offz);
         });

    BUGREPORT;

  /* Here the local sizes of fields are known, this is calculated on each proc */
   /* coordinates have not yet been defined or allocated */

    /************************************************/
    /*                                              */
    /* Allocation of variables and fields           */
    /*                                              */
    /************************************************/

    d->lnx = nx = d->elements[2];
    d->lny = ny = d->elements[1];
    d->lnz = nz = d->elements[0];

    //     Name of file to write flow energetics to
    sprintf(d->erhname,"%s.erh",d->name_out);
    sprintf(profilename,"%s.prof.dat",d->name_out);
    sprintf(mom_name,"%s.momentum.dat",d->name_out);

    // Allocation

    // vorticity
        Util_3DAllocFields(&w_0,&w_1,&w_2,0,&dtw_0,&dtw_1,&dtw_2,0,
                       &wbdr[0],&wbdr[1],&wbdrlow,&wbdrup,
                       nx,ny,nz,d->offx,d->offy,d->offz);

    // density
        Util_3DAllocFields(&n_0,&n_1,&n_2,0,&dtn_0,&dtn_1,&dtn_2,0,
                       &nbdr[0],&nbdr[1],&nbdrlow,&nbdrup,
                       nx,ny,nz,d->offx,d->offy,d->offz);

    //parallel electron speed
        Util_3DAllocFields(&V_0,&V_1,&V_2,0,&dtV_0,&dtV_1,&dtV_2,0,
                       &vbdr[0],&vbdr[1],&vbdrlow,&vbdrup,
                       nx,ny,nz,d->offx,d->offy,d->offz);

    // Parallel ion speed
        Util_3DAllocFields(&U_0,&U_1,&U_2,0,&dtU_0,&dtU_1,&dtU_2,0,
                       &ubdr[0],&ubdr[1],&ubdrlow,&ubdrup,
                       nx,ny,nz,d->offx,d->offy,d->offz);

    // electron temperature
        Util_3DAllocFields(&t_0,&t_1,&t_2,0,&dtt_0,&dtt_1,&dtt_2,0,
                       &tbdr[0],&tbdr[1],&tbdz[0],&tbdz[1],
                       nx,ny,nz,d->offx,d->offy,d->offz);

        // heatflux
        Util_3DAllocFields(&q_0,&q_1,&q_2,0,&dtq_0,&dtq_1,&dtq_2,0,
                       &qbdr[0],&qbdr[1],&qbdz[0],&qbdz[1],
                       nx,ny,nz,d->offx,d->offy,d->offz);


/* loeiten
  // Slow neutrals
        Util_3DAllocFields(&H_slow_0,&H_slow_1,&H_slow_2,0,&dtH_slow_0,&dtH_slow_1,&dtH_slow_2,0,
                       &H_slow_bdr[0],&H_slow_bdr[1],&H_slow_bdz[0],&H_slow_bdz[1],
                       nx,ny,nz,d->offx,d->offy,d->offz);

  // fast neutrals
        Util_3DAllocFields(&H_fast_0,&H_fast_1,&H_fast_2,0,&dtH_fast_0,&dtH_fast_1,&dtH_fast_2,0,
                       &H_fast_bdr[0],&H_fast_bdr[1],&H_fast_bdz[0],&H_fast_bdz[1],
                       nx,ny,nz,d->offx,d->offy,d->offz);
*/

        //  Some others
    Util_3DAllocFields(&f_0,&vr,&vp,0,&res,&n,&w,0,
                       &fbdr[0],&fbdr[1],&fbdrlow,&fbdrup,
                       nx,ny,nz,d->offx,d->offy,d->offz);



    BUGREPORT;

    dzF            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dzzF           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);

    dzN            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dzzN           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);

    dzV            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dzzV           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);

    VN             = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dzVN           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dzzVN          = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);

    UN             = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dzUN           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dzzUN          = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);


    dzU            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dzzU           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);

    dzt            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dzzt           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);

    dzq            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dzzq           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);

    dzJ            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dzzJ           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);

    dtdxf          = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dtdyf          = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dxn            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dyn            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dxt            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dyt            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dxf            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dyf            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    vstatic        = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    J_0            = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);

    tfeld          = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    dsource        = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    omega          = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    exp_n_0        = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    exp_t_0        = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    target         = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    nu_density     = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    nu_vmu         = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);

    laxn           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    laxt           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    laxq           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    laxv           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    laxw           = Util_DCube(nx,d->offx,ny,d->offy,nz,d->offz,0);
    BUGREPORT;

    ombdr[0]       = Util_DMatrix(nz,d->offz,ny,d->offy);
    ombdr[1]       = Util_DMatrix(nz,d->offz,ny,d->offy);
    nlval          = Util_DMatrix(nz,0,nx,0);// IMPORTANT: Tight field for Arakawa optimisation
    hval           = Util_DMatrix(nz,d->offz,nx,d->offx);
    vval           = Util_DMatrix(nz,d->offz,nx,d->offx);
    norm           = Util_DMatrix(nz,d->offz,nx,d->offx);

    // Vectors
    ppar_prof  = Util_DVector(nx,d->offx);
    vpar_prof  = Util_DVector(nx,d->offx);
    upar_prof  = Util_DVector(nx,d->offx);
    npar_prof  = Util_DVector(nx,d->offx);
    fpar_prof  = Util_DVector(nx,d->offx);
    wpar_prof  = Util_DVector(nx,d->offx);

    pax_prof  = Util_DVector(nz,d->offz);
    vax_prof  = Util_DVector(nz,d->offz);
    uax_prof  = Util_DVector(nz,d->offz);
    nax_prof  = Util_DVector(nz,d->offz);
    fax_prof  = Util_DVector(nz,d->offz);
    wax_prof  = Util_DVector(nz,d->offz);

    zcor   = Util_DVector(nz,d->offz);
    zhval  = Util_DVector(nz,d->offz);
    zgval  = Util_DVector(nz,d->offz);

    rcor              = Util_DVector(nx,d->offx);
    edrcor            = Util_DVector(nx,d->offx);

    par_rey_stress    = Util_DVector(nx,d->offx);
    flux_n            = Util_DVector(nx,d->offx);
    av_flux_n         = Util_DVector(nx,d->offx);
    dynbg             = Util_DVector(nx,d->offx);
    u_par             = Util_DVector(nx,d->offx);
    pol_pot           = Util_DVector(ny,d->offy);
    pol_dens          = Util_DVector(ny,d->offy);
    par_elec_mom_flux = Util_DVector(nx,d->offx);
    par_ion_mom_flux  = Util_DVector(nx,d->offx);
    ts_n              = Util_DVector((int)(p->out_time/p->dt)+100,1);
    ts_f              = Util_DVector((int)(p->out_time/p->dt)+100,1);
    ts_e              = Util_DVector((int)(p->out_time/p->dt)+100,1);
    BUGREPORT;




    /************************ INITIAL CONDITION **********************************************/

    /****************************************************/
    /*                                                  */
    /*   READ THE DATA FIELDS                           */
    /*                                                  */
    /****************************************************/

    BUGREPORT;
    MPI_Barrier(MPI_COMM_WORLD);
    d->ReadAttributes = FALSE;

    sprintf(d->integrator,"%s",argv[0]);
    sprintf(d->revision,"Empty");
    sprintf(d->compile_date,"%s %s",__DATE__,__TIME__);
    p->otmult = MAX(p->otmult,1);


    MPI_Bcast(p,1,para_type,ROOT,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    /*****************************************************************************/
    /*   Setup Geometry    */

    CHECK_STRUCTURE;
    if(p->coordsys == CYLINDRICAL) {
        ispolar = TRUE;
        sprintf(d->coordsys,"cylindrical");
        sprintf(d->dim_label[0],"z");
        sprintf(d->dim_label[2],"r");
        sprintf(d->dim_label[1],"phi");
        FORX_BD
        {
            r = p->xmin+(d->lnx*d->grid_coords[2]+ir+0.5)*p->dx;
            rcor[ir] = r;
            edrcor[ir] = 1./r;
            COMM(if(ISROOT)      fprintf(stderr,"%d r= %f\n",ir,r););
        }
    }
    else
    {
        ispolar = FALSE;
        sprintf(d->coordsys,"cartesian");
        sprintf(d->dim_label[0],"z");
        sprintf(d->dim_label[2],"x");
        sprintf(d->dim_label[1],"y");
        FORX_BD
        {
            rcor[ir] = 1;
            edrcor[ir] = 1.;
            COMM(if(ISROOT)      fprintf(stderr,"%d r= %f\n",ir,r););
        }
    }


    //if(ISROOT) {
        fprintf(stderr,"The geometry is ");
        if(ispolar)fprintf(stderr,"polar.\n");
        else fprintf(stderr,"rectangular.\n");

        COMM(fprintf(stderr,"%d coordsys: %s %d\t",__LINE__,d->coordsys, p->coordsys);
             fprintf(stderr,"with labels: %s %s %s \n",d->dim_label[2],d->dim_label[1],d->dim_label[0]););
        //fprintf(stderr,"%ld %ld %ld \n",d->dims[2],d->dims[1],d->dims[0]);
        //fprintf(stderr,"%ld %ld %ld \n",d->nx,d->ny,d->nz);
        //}




    BUGREPORT;

    // Check for z-Tschebischeff
    // This must take respect for domain decomposition
    // And if 0 is not in the z domain, we presume that the lower boundary has no sheath, ergo no denser points!


    COMM(fprintf(stderr,"Z spacing is %d\n", p->z_spacing););

    BUGREPORT;
    if (p->z_spacing == 1)
    {
        off = p->z_offset;

        COMM(fprintf(stderr,"%s: Setup z-Tschebicheff\n",__FILE__);    );

        if(p->zmin < -0.000001)
        {
            zm = 0.5*(p->zmax + p->zmin); // midpoint
            Lz = p->zmax - p->zmin; //Length of domain to be mapped
            dz = M_PI/(double)(d->dims[0]+2*off); // delta of cos distributed points
            xlen = cos( ((double)off+0.5) *dz)   -  cos( ((double)(off+d->dims[0]) + 0.5) *dz); // Length of Tschebischef domain, with offset

            COMM(fprintf(stderr,"xlen:  %.12f\n",xlen););
            g = xlen/Lz;
            gm1 = 1./g;

            BUGREPORT;

            // calculate locations of z planes
            // Note: These are LOCAL planes

            for(iz=-1;iz<=nz;iz++)
            {
                xlen = dz*( (double)(nz*d->grid_coords[0]+iz+off) +0.5);
                zcor[iz]  = -gm1*cos(xlen) + zm;
                COMM(fprintf(stderr,"%d xlen %.12f %.12f\n",iz,xlen, zcor[iz]););

            }
        }// now if zmin >= zero
        else
        {
            zm = p->zmin; // midpoint
            Lz = 2.*(p->zmax - p->zmin); //Length of domain to be mapped, but we only use half the points here
            dz = M_PI/(double)(2.*d->dims[0]+2*off); // delta of cos distributed points
            xlen = cos( ((double)off+0.5) *dz)   -  cos( ((double)(off+2*d->dims[0]) + 0.5) *dz); // Length of Tschebischef domain, with offset

            COMM(fprintf(stderr,"xlen:  %.12f\n",xlen););
            g = xlen/Lz;
            gm1 = 1./g;
            BUGREPORT;

            // calculate locations of z planes
            // Note: These are LOCAL planes

            for(iz=-1;iz<=nz;iz++)
            {
                xlen = dz*( (double)(d->dims[0]+ nz*d->grid_coords[0] + iz +off) +0.5);
                zcor[iz]  = -gm1*cos(xlen) + zm;
                COMM(fprintf(stderr,"%d xlen %.12f %.12f\n",iz,xlen, zcor[iz]););
            }

            BUGREPORT;

        }

        for(iz=-1;iz<=nz;iz++)
        {
            val = g*(zcor[iz]-zm);
            zhval[iz] = g/sqrt(1.-val*val)/dz;
            zgval[iz] = g*g*val/dz/((1.-val*val)*sqrt(1.-val*val));
        }
    }
    else
    {
        dz = (p->zmax-p->zmin)/(double)(d->dims[0]);

        for(iz=-1;iz<=nz;iz++)
        {
            zcor[iz]  = p->zmin+dz*((double)(nz*d->grid_coords[0]+iz)+0.5);
            COMM( fprintf(stderr,"%d %.12f\n",iz,zcor[iz] ););
        }

        for(iz=-1;iz<=nz;iz++)
        {
            zhval[iz] = 1./dz;
            zgval[iz] = 0.;
        }


    }

   BUGREPORT;


   COMM(MPI_Gather(&zhval[0],nz,MPI_DOUBLE,&d->coordinate[0][0],nz,MPI_DOUBLE,ROOT,d->zrow_comm);

   if(ISROOT) {
       output = fopen("zhval", "w");
       for(iz=0;iz<d->dims[0];iz++)  fprintf(output,"%d %.12f\n",iz,d->coordinate[0][iz] );
       fclose(output);
   }

   MPI_Gather(&zgval[0],nz,MPI_DOUBLE,&d->coordinate[0][0],nz,MPI_DOUBLE,ROOT,d->zrow_comm);

   if(ISROOT) {
       output = fopen("zgval", "w");
       for(iz=0;iz<d->dims[0];iz++)  fprintf(output,"%d %.12f\n",iz,d->coordinate[0][iz] );
       fclose(output);
   });





   // Gather results
   MPI_Gather(&zcor[0],nz,MPI_DOUBLE,&d->coordinate[0][0],nz,MPI_DOUBLE,ROOT,d->zrow_comm);
   BUGREPORT;

   MPI_Barrier;


   COMM(if(ISROOT) for(iz=0;iz<d->dims[0];iz++)  fprintf(stderr,"%d %.12f\n",iz,d->coordinate[0][iz] ););




    /*****************************************************************************/
    /* RBDCND:      X = A           X = B */
    BUGREPORT;
    /*****************************************************************************/
    if(p->xmin == 0. && ispolar)
    {
        rbdcndf   = ZeroDIR;
        rbdcndw   = ZeroDIR;
        rbdcndn   = ZeroDIR;
        rbdcndu   = ZeroNEU;
        rbdcndv   = ZeroNEU;
        rbdcndt   = ZeroNEU;
        rbdcndq   = ZeroNEU;
    }
    else
    {
        rbdcndf   = NEUDIR;
        rbdcndw   = NEUDIR;
        rbdcndt   = NEUDIR;
        rbdcndq   = NEUDIR;
        rbdcndn   = NEUDIR;
        rbdcndu   = NEUNEU;
        rbdcndv   = NEUNEU;
    }


    //Values on the boundaries
    BUGREPORT;

    for(i=0;i<2;i++)
    {
        FORZY_BD nbdr[i][iz][ip]  = 0.;
        FORZY_BD fbdr[i][iz][ip]  = 0.;
        FORZY_BD wbdr[i][iz][ip]  = 0.;
        FORZY_BD ombdr[i][iz][ip] = 0.;
        FORZY_BD ubdr[i][iz][ip]  = 0.;
        FORZY_BD vbdr[i][iz][ip]  = 0.;
        FORZY_BD tbdr[i][iz][ip]  = 0.;
        FORZY_BD qbdr[i][iz][ip]  = 0.;
/*loeiten
                 FORZY_BD H_slow_bdr[i][iz][ip]  = 0.;
        FORZY_BD H_fast_bdr[i][iz][ip]  = 0.;
*/
    }

    /*****************************************************************************/
    /* ZBDCND:      Z = down           Z = up*/


    zbndcndf = DIRNEU;
    zbndcndw = DIRNEU;
    zbndcndn = DIRNEU;
    zbndcndt = DIRNEU;
    zbndcndq = DIRNEU;
    zbndcndu = NEUNEU;
    zbndcndv = NEUNEU;

    BUGREPORT;
    FORYX_BD wbdrup[ip][ir]  = 0.;
    FORYX_BD wbdrlow[ip][ir] = 0.;
    FORYX_BD nbdrup[ip][ir]  = 0.;
    FORYX_BD nbdrlow[ip][ir] = 0.;
    FORYX_BD fbdrup[ip][ir]  = 0.;
    FORYX_BD fbdrlow[ip][ir] = 0.;
    FORYX_BD ubdrup[ip][ir]  = 0.;
    FORYX_BD ubdrlow[ip][ir] = 0.;
    FORYX_BD vbdrup[ip][ir]  = 0.;
    FORYX_BD vbdrlow[ip][ir] = 0.;
    FORYX_BD tbdz[0][ip][ir] = 0.;
    FORYX_BD tbdz[1][ip][ir] = 0.;
    FORYX_BD qbdz[0][ip][ir] = 0.;
    FORYX_BD qbdz[1][ip][ir] = 0.;

/*loeiten
        FORYX_BD H_slow_bdz[0][ip][ir] = 0.;
    FORYX_BD H_slow_bdz[1][ip][ir] = 0.;
        FORYX_BD H_slow_bdz[0][ip][ir] = 0.;
    FORYX_BD H_slow_bdz[1][ip][ir] = 0.;
*/
    BUGREPORT;






        // PUT RESTART CAPABILITIES IN SEPERATE SOURCE FILE




    if(d->restart == RESTART)
    {
     /******************RESTART FROM FILE MADE BY CYTO ****************************************/
     /* Read all the fields from a file, but allow overide of parameters from commandline,
        e.g. do not read parameters here ! */
     d->read_data = TRUE;
     d->ReadAttributes = FALSE;
     BUGREPORT;
     counter = 0;


     // Read fields which might appear under different names

     // Density
     if(PH_Read3DFieldbyName(n_0,"N",d,p) > 0 )
     {
         if(ISROOT) fprintf(stderr,"Restart success, read field N\n");
     }
     else if(PH_Read3DFieldbyName(n_0,"n",d,p)> 0  )
     {
         FORALL_BD  n_0[iz][ip][ir]  = log(n_0[nz][ip][ir]);
         if(ISROOT) fprintf(stderr,"Restart success, read field n \n");
     }
     else if(PH_Read3DFieldbyName(n_0,"Density",d,p)> 0  )
     {
         FORALL_BD  n_0[iz][ip][ir]  = log(n_0[nz][ip][ir]);
         if(ISROOT) fprintf(stderr,"Restart success, read field Density \n");
     }
     else
     {
         if(ISROOT) fprintf(stderr,"Restart failure, could not read field N\n");
         counter++;
     }
     COMM(PH_3D_Write( n_0 ,"N","RESTART",1,d,p,TRUE);     );



     // Field H
     if(PH_Read3DFieldbyName(w_0,"H",d,p) > 0  )
     {
         if(ISROOT)  fprintf(stderr,"Restart success, read field H\n");
         COMM(PH_3D_Write( w_0 ,"H","RESTART",2,d,p,TRUE);  );
     }
     else
     {
         if(ISROOT)  fprintf(stderr,"Restart failure, could not read field H\n");
         counter++;
     }


     // ION parallel velocity
     if(PH_Read3DFieldbyName(U_0,"UIon",d,p)> 0 )
     {
         if(ISROOT) fprintf(stderr,"Restart success, read field UIon\n");
         COMM(PH_3D_Write( U_0 ,"UIon","RESTART",3,d,p,TRUE););
     }
     else
     {
         if(ISROOT) fprintf(stderr,"Restart failure, could not read field UIon\n");
         counter++;
     }
     COMM(PH_3D_Write( U_0 ,"UIon","RESTART",3,d,p,TRUE););





     //Electron parallel velocity
     if(PH_Read3DFieldbyName(V_0,"VElectron",d,p)> 0  )
     {
         if(ISROOT) fprintf(stderr,"Restart success, read field VElectron\n");
         COMM(PH_3D_Write( V_0 ,"VElectron","RESTART",4,d,p,TRUE);     );
    }
     else if(PH_Read3DFieldbyName(V_0,"Velocity",d,p)> 0  )
     {
         if(ISROOT) fprintf(stderr,"Restart success, read field Velocity\n");
         COMM(PH_3D_Write( V_0 ,"V_0","RESTART",4,d,p,TRUE);     );
     }
     else
     {
        if(ISROOT) fprintf(stderr,"Restart failure, could not read field VElectron\n");
         counter++;
     }


#ifdef TE
     //Electron temperature
     if(PH_Read3DFieldbyName(t_0,"Te",d,p)> 0  )
     {
         if(ISROOT) fprintf(stderr,"Restart success, read field VElectron\n");
         COMM(PH_3D_Write( V_0 ,"VElectron","RESTART",4,d,p,TRUE);     );
     }
     else
     {
        if(ISROOT) fprintf(stderr,"Restart failure, could not read field VElectron\n");
         counter++;
     }
#endif //TE

     // Potential
     if(PH_Read3DFieldbyName(f_0,"Potential",d,p)> 0 )
     {
       if(ISROOT)  fprintf(stderr,"Restart success, read field Potential\n");
       COMM(PH_3D_Write( f_0 ,"Potential","RESTART",5,d,p,TRUE);     );
     }
     else
     {
       if(ISROOT)  fprintf(stderr,"Restart failure, could not read field Potential\n");
         counter++;
     }

     // Check if all fields were found
     if (counter > 0)
     {
       if(ISROOT)  fprintf(stderr,"Restart failure, could not read %d fields. Cyto cannot restart from this file.\n",counter);
         MPI_Finalize();
         exit(0);
     }
     else
     {
         fprintf(stderr,"Proc %d: Restart success. Can restart from this file.\n",d->this_process);
     }


   }
   else if (d->restart == START_FROM_INI || d->restart == DEFAULTSTART )
   {
       // printf("START_FROM_INI...\n");

       /******************START FROM INI FILE  ****************************************/

     Lz = (p->zmax - p->zmin);

       FORALL_BD
       {
           izg = nz*d->grid_coords[0]+iz;

           z = zcor[iz]/Lz;
           r = p->xmin+(d->lnx*d->grid_coords[2]+ir+0.5)*p->dx;
           n_0[iz][ip][ir] =  min_nt  + 0.1*exp(-(r*r*p->kappan*p->kappan + 10.*10.*z*z) );
       }


       for(m=2;m<ny/8;m++)
       {
           phase =  (double)rand()/(double)RAND_MAX ;

           FORALL_BD  {
               irg = ir + nx*d->grid_coords[2];
               n_0[iz][ip][ir]+=  n_0[iz][ip][ir]*0.001/(double)(m*m)
               *cos(2.*M_PI*m*((double)ip/(double)ny+phase))
                   *exp(-( (double)(irg-d->dims[2]/2)/(double)d->dims[2])*((irg-d->dims[2]/2)/(double)d->dims[2])*500.   );
           }

       }


       FORALL_BD  n_0[iz][ip][ir]  = log(n_0[iz][ip][ir]);
#ifdef TE
       FORALL_BD  t_0[iz][ip][ir]  = 0.;
#else
       FORALL_BD  t_0[iz][ip][ir]  = 0.;
#endif

       BUGREPORT;

       /* w ,f  get initialized with zero values */
       FORALL_BD omega[iz][ip][ir]  = 0.;


      FORALL_BD
       {
           r = p->xmin+(d->lnx*d->grid_coords[2]+ir+0.5)*p->dx;
           w_0[iz][ip][ir] =   0.1*exp(-(r*r*p->kappan*p->kappan) );
       }
       FORALL_BD f_0[iz][ip][ir]    = 0.;
       FORALL_BD q_0[iz][ip][ir]    = 0.;


       BUGREPORT;

       // U V are monotonous initially
       Lz = (p->zmax - p->zmin);
       if(p->zmin < 0)
       {
           if(p->zmax < 0)
           {
               fprintf(stderr,"ERROR:  zmax <=0 not allowed !\n");
               exit(-1);
           }

           FORALL
           {
               z = zcor[iz];
               if(z>0)
               {
                   V_0[iz][ip][ir] = z/p->zmax;
                   U_0[iz][ip][ir] = z/p->zmax;
               }
               else
               {
                   V_0[iz][ip][ir] = -z/p->zmin;
                   U_0[iz][ip][ir] = -z/p->zmin;
               }
           }
       }
       else
       {
           FORALL   V_0[iz][ip][ir] = zcor[iz]/Lz;
           FORALL   U_0[iz][ip][ir] = zcor[iz]/Lz;
       }

       BUGREPORT;
       // from zero to 1 with zero derivative at 0, note potential is relative to sheath potential!!!
       FORALL   V_0[iz][ip][ir]*=fabs(V_0[iz][ip][ir]); /*1./sqrt(2.*M_PI)*sqrt(p->Mi/ME)*/ ;
       FORALL   U_0[iz][ip][ir]*=fabs(U_0[iz][ip][ir]);

       BUGREPORT;
       for(m=2;m<ny/8;m++)
       {
           phase =  (double)rand()/(double)RAND_MAX ;

           /*   FORALL   V_0[iz][ip][ir]+=
               0.001/(double)(m*m)
               *cos(2.*M_PI*m*((double)ip/(double)ny+phase))
               *exp(-( (double)(ir-nx/2)/(double)nx)*((ir-nx/2)/(double)nx)*500.   );*/

           /*FORALL   U_0[iz][ip][ir]+=
               0.001/(double)(m*m)*sin(2.*M_PI*(double)iz*(1.+0.5*(double)rand()/(double)RAND_MAX)/(double)nz)
               *cos(2.*M_PI*m*(double)ip/(double)ny)*sin(M_PI*(double)ir/(double)nx)/(double)(2+iz);*/
       }
   }
   else if (d->restart == START_FROM_FILE )
   {
       printf("START_FROM_FILE...\n");

       /******************START FROM .000 FILE  ****************************************/
       fprintf(stderr, "Start from .000 file: nx = %d ny = %d nz = %d \n",nx,ny,nz);
       BUGREPORT;
       d->read_data = TRUE;
       d->ReadAttributes=FALSE;
       BUGREPORT;
       PH_Read3DFieldbyNumber(exp_n_0,1,d,p);


       for(ir=0;ir<nx;++ir) printf("%f ",exp_n_0[nz/2][ny/2][ir]); printf("\n");

       // Make sure that n is positive
       FORALL n_0[iz][ip][ir] =  log(fabs(exp_n_0[0][ip][ir]));


       // U V are monotonous initially
       FORALL   V_0[iz][ip][ir] = -1.+ 2./Lz*zcor[iz];
       FORALL   U_0[iz][ip][ir] = -1.+ 2./Lz*zcor[iz];
       FORALL   w_0[iz][ip][ir] = omega[iz][ip][ir] = f_0[iz][ip][ir]  =  0.;
       fprintf(stderr, "nx = %d ny = %d nz = %d \n",nx,ny,nz);
   }
   else
   {
       fprintf(stderr, "Restart undefined \n");
       MPI_Finalize();
       return 0;
   }
    MPI_Bcast(p,1,para_type,ROOT,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    BUGREPORT;
    /****************************************************/
    /*                                                  */
    /* Calculate  dimensionless parameters                             */
    /*                                                  */
    /****************************************************/

    FORZX_BD hval[iz][ir] = 1./(p->dx);
    FORZX_BD vval[iz][ir] = edrcor[ir]/(p->dy);
    FORZX  nlval[iz][ir]  = vval[iz][ir]*hval[iz][ir];
    FORZX  norm[iz][ir]       = 1./(zhval[iz]*hval[iz][ir]*vval[iz][ir]);


    COMM(FORX fprintf(stderr,"hval[%d] = %f \t vval[%d]=%f \t nlval = %f rcor = %f \n", ir,hval[0][ir],ir,vval[0][ir],nlval[0][ir],rcor[ir] );
         );


    mu     = p->Mi/ME;
    delta  = p->delta;
    nu     = p->nu;/* assume nu to be determined from plasma core while it should be plasma edge)*/
    nu_m   = nu/mu;
    mum1         = 1./mu;
    aL           = sqrt(mum1)/Lz;
    lmin_nt = log(min_nt);
    if(p->hm_nl > 0.0) padvection = 1.0;

    CHECK_STRUCTURE;
    if(ISROOT) fprintf(stderr,"kappan                             : %f\n",p->kappan);
    if(ISROOT) fprintf(stderr,"Electron parallel drift (alpha)    : %f\n",p->alpha);
    if(ISROOT) fprintf(stderr,"Electron neutral collisions (delta): %g\n",p->delta);
    if(ISROOT) fprintf(stderr,"Electron Ion collisions (nu)       : %g\n",p->nu);
    if(ISROOT) fprintf(stderr,"Ion neutral collisions (sigma)     : %g\n",p->sigma);
    if(ISROOT) fprintf(stderr,"Dissipations                       : W = %g     N = %g   T = %g\n",p->mue_w,p->mue_n,p->mue_t);
    if(ISROOT) fprintf(stderr,"Density source                     : %f\n",p->nprof);
    if(ISROOT) fprintf(stderr,"Heat source                        : %f\n",p->tprof);
    if(ISROOT) fprintf(stderr,"Avoid density/temperature to zero  : %f\n",p->source);
    if(ISROOT) fprintf(stderr,"Parallel diffusion                 : %f\n",p->limiter);
    if(ISROOT) fprintf(stderr,"Parallel advection                 : %d\n",(int)padvection);
    if(ISROOT) fprintf(stderr,"Zeff                               : %f\n",p->Z);
    if(ISROOT) fprintf(stderr,"Coordinates                        : \n");
    if(ISROOT) fprintf(stderr," cordsys                           : %s\n",d->coordsys);
    if(ISROOT) fprintf(stderr," X                                 : [%f %f] [%d]\n",p->xmin,p->xmax,p->r_spacing);
    if(ISROOT) fprintf(stderr," Y                                 : [%f %f]  [0]\n",p->ymin,p->ymax);
    if(ISROOT) fprintf(stderr," Z                                 : [%f %f] [%d]\n",p->zmin,p->zmax,p->z_spacing);
    if(ISROOT) fprintf(stderr," NX (nx,lnx)                       :  %d\t (%d,%d)\n",d->N[2],d->dims[2],nx);
    if(ISROOT) fprintf(stderr," NY (ny,lny)                       :  %d\t (%d,%d)\n",d->N[1],d->dims[1],ny);
    if(ISROOT) fprintf(stderr," NZ (nz,lnz)                       :  %d\t (%d,%d)\n",d->N[0],d->dims[0],nz);
    if(ISROOT) fprintf(stderr,"\n");

    if(ISROOT) fprintf(stderr,"Time Integration                   : \n");
    if(ISROOT) fprintf(stderr," dt                                : %f\n", p->dt);
    if(ISROOT) fprintf(stderr," Integral outputs                  : %f\n", p->out_time);
    if(ISROOT) fprintf(stderr," Field outputs                     : %f\n", p->out_time*p->otmult);
    if(ISROOT) fprintf(stderr," Stopping at                       : %f   (%ld timesteps)\n", p->end_time,(long)((p->end_time - p->time)/p->dt));
    BUGREPORT;


   //set up density source :Localisation of density source in first quarter of device

   FORALL_BD
   {
    r = p->xmin+(d->lnx*d->grid_coords[2]+ir+0.5)*p->dx;
    z = zcor[iz];

    dsource[iz][ip][ir] =  exp(-r*r*(p->kappan*p->kappan));  //    *exp(-((z/Lz)*(z/Lz)));
   }

   // set up a shape function to drive vorticity to zero in the vicinity of the radial wall

   FORALL {
       r = p->xmin+(d->lnx*d->grid_coords[2]+ir+0.5)*p->dx;
       target [iz][ip][ir] = exp(-4.*(p->xmax-r)*(p->xmax-r));
   }


   /****************************************************/
   /*                                                  */
   /* WHAT ARE THE DEFINES                             */
   /*                                                  */
   /****************************************************/

   sprintf(d->desc,"SWITCHES:\n");

   strcat(d->desc,"Order of Parallel Derivative - ");
   if(2 == d->offz)
     strcat(d->desc,"4th.\n");
   else
     strcat(d->desc,"2nd.\n");

   strcat(d->desc,"Electron parallel dynamics - ");
#ifdef VSTATIC
   strcat(d->desc,"static.\n");
#else
   strcat(d->desc,"dynamic.\n");
#endif

#ifdef TE
   strcat(d->desc,"Uses Electron Temperature equation.\n");
#else
   strcat(d->desc,"Electron Temperature constant.\n");
#endif //TE

   strcat(d->desc,"Electron parallel heat flux dynamics - ");
#ifdef QSTATIC
   strcat(d->desc,"static.\n");
#else
   strcat(d->desc,"dynamic.\n");
#endif


   if(ISROOT) fprintf(stderr,"%s\n",d->desc);



/*****************************************************************/
/*                                                               */
/*          Start of Time Loop                                   */
/*                                                               */
/*****************************************************************/
 MPI_Barrier(MPI_COMM_WORLD);
 write_time     = p->time;
 ot = p->otmult; // Make sure first file is written
 CHECK_STRUCTURE;
 while(p->time < p->end_time )
  {
      if(FALSE == first )
      {
          BUGREPORT;

           /***********************************************/
/*loeiten
                  / * Slow Neutral interaction * /

          C
*/





          /*********************************************/
          /*                   Electron Density        */
          /*********************************************/
          CONVECT(dtn_0,f_0,n_0,nlval,nz,nx,ny);
          FORALL dtn_0[iz][ip][ir]   +=  -dzVN[iz][ip][ir];

          // density source
          FORALL    dtn_0[iz][ip][ir] += p->nprof*dsource[iz][ip][ir]/exp_n_0[iz][ip][ir];

          // artificial "ionisation", keeps log density from falling below min_nt
          // FORALL    dtn_0[iz][ip][ir]-=  p->source* (n_0[iz][ip][ir]-lmin_nt -fabs(n_0[iz][ip][ir]-lmin_nt));

          // PARALLEL_DAMPING
          FORALL       dtn_0[iz][ip][ir]+=p->limiter*(dzzN[iz][ip][ir] +  dzN[iz][ip][ir] *dzN[iz][ip][ir]);//+ 0.1*laxn[iz][ip][ir];

          // Explicit part of damping
          FORALL dtn_0[iz][ip][ir] += p->mue_n*(dxn[iz][ip][ir]*dxn[iz][ip][ir] + dyn[iz][ip][ir]*dyn[iz][ip][ir]);

          Util_3DSsTimeStep(iter,nx,ny,nz,p->dt,p->mue_n,&lamda,
                            n,n_0,n_1,n_2,dtn_0,dtn_1,dtn_2);

          Laplace_Solve3D(d,p,n_2,n,nbdr[0],nbdr[1],hval[0],rcor,rbdcndn,lamda,TRUE);

          // Calculate dt ln n .....with implicit term
          Cyto_Laplace(res,n_0,hval,vval,edrcor,nx,ny,nz,p,ispolar);
          FORALL tfeld[iz][ip][ir] = dtn_0[iz][ip][ir] +  p->mue_n*res[iz][ip][ir];
          // UPDATE BOUNDARY ON dtn
          PH_Update2dBoundaries(tfeld,rbdcndf, qbdr[0], qbdr[1],hval,d);

          BUGREPORT;

          /*********************************************/
          /* Global Vorticity equation: Needs to be evaluated AFTER dt_n_0 is known  */
          /*********************************************/
          //  Note that w = omega + grad N grad phi
          CONVECT(dtw_0,f_0,omega,nlval,nz,nx,ny);

                  // density source, compensate w
                  FORALL_BD    tfeld[iz][ip][ir]= p->nprof*(dsource[iz][ip][ir]/exp_n_0[iz][ip][ir]);
          FORALL    dtw_0[iz][ip][ir] += DX(tfeld)*DX(f_0);


          // Parallel terms
          FORALL dtw_0[iz][ip][ir] +=  (dzUN[iz][ip][ir] - dzVN[iz][ip][ir]/p->Z);

          // Collisional Term: Pedersen current
          FORALL dtw_0[iz][ip][ir] +=  -p->sigma*omega[iz][ip][ir];

          // \nabla (d_t n) \nabla \phi term
          //FORALL dtw_0[iz][ip][ir] += DX(tfeld)*dxf[iz][ip][ir]+DY(tfeld)*dyf[iz][ip][ir];

          // Nonlinear Term: -grad N {phi, grad phi }, sign like in main nonlinear term!
          //CONVECT(dtdxf,f_0,dxf,nlval,nz,nx,ny);
          //CONVECT(dtdyf,f_0,dyf,nlval,nz,nx,ny);
          //FORALL dtw_0[iz][ip][ir] += dxn[iz][ip][ir]*dtdxf[iz][ip][ir]+dyn[iz][ip][ir]*dtdyf[iz][ip][ir];

          //     Calculate \nabla (sigma  phi - mue w)
          //FORALL_BD res[iz][ip][ir] =   p->mue_w*omega[iz][ip][ir] -p->sigma*f_0[iz][ip][ir];
          //FORALL dtw_0[iz][ip][ir] +=   dxn[iz][ip][ir]*DX(res) + dyn[iz][ip][ir]*DY(res);

          // damp to zero in target range
          //FORALL dtw_0[iz][ip][ir] -= w_0[iz][ip][ir]*target[iz][ip][ir];


          // PARALLEL_DAMPING
          //FORALL       dtw_0[iz][ip][ir]-=  p->limiter*dzzF[iz][ip][ir] ;//- laxw[iz][ip][ir];

          Util_3DSsTimeStep(iter,nx,ny,nz,p->dt,p->mue_w,&lamda,
                            w,w_0,w_1,w_2,dtw_0,dtw_1,dtw_2);


          // Damping only on vorticity, omega, thus substract grad n grad phi term ,
          FORALL tfeld[iz][ip][ir]  =  dxn[iz][ip][ir]*dxf[iz][ip][ir]+ dyn[iz][ip][ir]*dyf[iz][ip][ir];
          FORALL  w[iz][ip][ir]  -=  tfeld[iz][ip][ir];

          Laplace_Solve3D(d,p,w_2,w,ombdr[0],ombdr[1],hval[0],rcor,rbdcndw,lamda,TRUE);
          //add grad n grad phi term again
          FORALL  w_2[iz][ip][ir]  +=  tfeld[iz][ip][ir];


#ifdef TE
          /*********************************************/
          /*       log Electron temperature            */
          /*********************************************/
          CONVECT(dtt_0,f_0,t_0,nlval,nz,nx,ny);

          // source
          //FORALL  dtt_0[iz][ip][ir]  += p->tprof*dsource[iz][ip][ir]/(exp_t_0[iz][ip][ir]*exp_n_0[iz][ip][ir]);

          // parallel gradient
          FORALL  dtt_0[iz][ip][ir]  += -(V_0[iz][ip][ir]+p_>alpha)*dzt[iz][ip][ir] - 0.0*tf*dzV[iz][ip][ir];

          // artificial "heating", keeps log temperature from falling below t_0
          //FORALL    dtt_0[iz][ip][ir]-=  p->source* (t_0[iz][ip][ir] -fabs(t_0[iz][ip][ir]));

         // Curvature
         //FORALLBD   tpn1[iz][ip][ir] = 3.5* t_0[iz][ip][ir] + n_0[iz][ip][ir] - f_0[iz][ip][ir];
         //CONVECT(tfeld,tpn1,B_0,nlval,nz,nx,ny);
         //FORALL       dtt_0[iz][ip][ir] += tf*tfeld[iz][ip][ir];


         // Ohmic heating
         FORALL  dtt_0[iz][ip][ir]  += tf*nu_density[iz][ip][ir]*J_0[iz][ip][ir]/exp_n_0[iz][ip][ir]*J_0[iz][ip][ir]/exp_t_0[iz][ip][ir];

         // Parallel Divergence of heat flux
         FORALL  dtt_0[iz][ip][ir] +=  -tf*dzq[iz][ip][ir]/(exp_t_0[iz][ip][ir]*exp_n_0[iz][ip][ir]);

         // Profile screening
         //FORALL  dtt_0[iz][ip][ir] += -tf*screen[k]*tbg[i][k];

         // Thermal force
         FORALL  dtt_0[iz][ip][ir]  -= tf*Alpha*dzt[iz][ip][ir]*J_0[iz][ip][ir]/exp_n_0[iz][ip][ir];

         // PARALLEL_DAMPING
         FORALL  dtt_0[iz][ip][ir] += tf*p->limiter*( dzzt[iz][ip][ir]+  dzt[iz][ip][ir] *dzt[iz][ip][ir]) + laxt[iz][ip][ir];;

         // Explicit part of damping
         FORALL dtt_0[iz][ip][ir] += p->mue_t*(dxt[iz][ip][ir]*dxt[iz][ip][ir] + dyt[iz][ip][ir]*dyt[iz][ip][ir]);

         Util_3DSsTimeStep(iter,nx,ny,nz,p->dt,tf*p->mue_t,&lamda,tfeld,t_0,t_1,t_2,dtt_0,dtt_1,dtt_2);
         Laplace_Solve3D(d,p,t_2,tfeld,tbdr[0],tbdr[1],hval[0],rcor,rbdcndt,lamda,TRUE);


         /******************************************************************/
         /*  Parallel Electron Heat Flux  Equation                         */
         /******************************************************************/
         CONVECT(dtq_0,f_0,q_0,nlval,nz,nx,ny);

         // Parallel advection
         dtq_0[iz][ip][ir] -=  dzq[iz][ip][ir]*V_0[iz][ip][ir];

         // Temperature Gradient force
         FORALL  dtq_0[iz][ip][ir] += -5./2.*3.2*mum1*exp_n_0[iz][ip][ir]*exp_t_0[iz][ip][ir]*exp_t_0[iz][ip][ir]*dzt[iz][ip][ir];

         // Landau Damping
         FORALL  dtq_0[iz][ip][ir] += -aL*q_0[iz][ip][ir];

         // Parallel electron velocity
         //FORALL  dtq_0[iz][ip][ir] += -6./5.*q_0[iz][ip][ir]*dzV[iz][ip][ir];

         // collisional terms
         FORALL  dtq_0[iz][ip][ir] += -5./2.*nu_density[iz][ip][ir]/Kappa*(q_0[iz][ip][ir] + Alpha *exp_t_0[iz][ip][ir]*J_0[iz][ip][ir]);

         // PARALLEL_DAMPING
         FORALL  dtq_0[iz][ip][ir] += p->limiter*dzzq[iz][ip][ir] + laxq[iz][ip][ir]; ;

         Util_3DSsTimeStep(iter,nx,ny,nz,p->dt,p->mue_t,&lamda,tfeld,q_0,q_1,q_2,dtq_0,dtq_1,dtq_2);
         Laplace_Solve3D(d,p,q_2,tfeld,qbdr[0],qbdr[1],hval[0],rcor,rbdcndq,lamda,TRUE);

#endif // TE

          /*********************************************/
          /*       Ion velocity                        */
          /*********************************************/

          CONVECT(dtU_0,f_0,U_0,nlval,nz,nx,ny);

          // Parallel advection
          FORALL dtU_0[iz][ip][ir]  += -padvection*dzU[iz][ip][ir]*U_0[iz][ip][ir];
          BUGREPORT;

          FORALL dtU_0[iz][ip][ir]  +=  -p->Z * dzF[iz][ip][ir]; /* Z = 1 in GW version */

          // Neutral collision Term
          FORALL dtU_0[iz][ip][ir]  += -p->sigma*U_0[iz][ip][ir] ;

          // Resistivity
          FORALL dtU_0[iz][ip][ir]  += nu_vmu[iz][ip][ir];

          // Thermal force
          FORALL dtU_0[iz][ip][ir]  += -Alpha/mu*dzt[iz][ip][ir];

          //FORALL dtU_0[iz][ip][ir] -= U_0[iz][ip][ir]*target[iz][ip][ir];

          // PARALLEL_DAMPING
          FORALL       dtU_0[iz][ip][ir]+=  dzzU[iz][ip][ir]; //Factor is one in GWV0

          Util_3DSsTimeStep(iter,nx,ny,nz,p->dt,p->mue_n,&lamda,
                            tfeld,U_0,U_1,U_2,dtU_0,dtU_1,dtU_2);

          FORALL       U_2[iz][ip][ir] = tfeld[iz][ip][ir];
          // Laplace_Solve3D(d,p,U_2,tfeld,ubdr[0],ubdr[1],hval[0],rcor,rbdcndu,lamda,TRUE);


          /*********************************************/
          /*      Electron  Velocity                   */
          /*********************************************/
          BUGREPORT;
          CONVECT(dtV_0,f_0,V_0,nlval,nz,nx,ny);

          // Parallel Forces
          FORALL dtV_0[iz][ip][ir]   +=  mu*(dzF[iz][ip][ir] - exp_t_0[iz][ip][ir]*(dzN[iz][ip][ir]+dzt[iz][ip][ir]));

          // Parallel Advection
          FORALL dtV_0[iz][ip][ir] += -0.*padvection*dzV[iz][ip][ir]*(V_0[iz][ip][ir]+p->alpha);

          // Thermal force
          FORALL dtV_0[iz][ip][ir] += -0.*Alpha*dzt[iz][ip][ir];

          // Neutral collision Term
          FORALL dtV_0[iz][ip][ir]   += -0.*p->delta*(V_0[iz][ip][ir]+p->alpha);

          // Coloumb collision Term
          FORALL dtV_0[iz][ip][ir]   += -mu*nu_vmu[iz][ip][ir] ;

          // Parallel velocity towards zero at wall
          //FORALL dtV_0[iz][ip][ir] -= V_0[iz][ip][ir]*target[iz][ip][ir];

          //PARALLEL_DAMPING
          FORALL dtV_0[iz][ip][ir]   +=  10.*p->limiter*dzzV[iz][ip][ir]+ 100.*laxv[iz][ip][ir];

          BUGREPORT;

          Util_3DSsTimeStep(iter,nx,ny,nz,p->dt,p->mue_w,&lamda,
                            tfeld,V_0,V_1,V_2,dtV_0,dtV_1,dtV_2);

          FORALL       V_2[iz][ip][ir] = tfeld[iz][ip][ir];
          //Laplace_Solve3D(d,p,V_2,tfeld,vbdr[0],vbdr[1],hval[0],rcor,rbdcndv,lamda,TRUE);


          /*********************************************************************/
          /* Exchange pointers                                                 */
          /*********************************************************************/

          tmp = w_2;     w_2 = w_1;     w_1 = w_0;     w_0 = tmp;
          tmp = n_2;     n_2 = n_1;     n_1 = n_0;     n_0 = tmp;
          tmp = U_2;     U_2 = U_1;     U_1 = U_0;     U_0 = tmp;
          tmp = V_2;     V_2 = V_1;     V_1 = V_0;     V_0 = tmp;
          tmp = t_2;     t_2 = t_1;     t_1 = t_0;     t_0 = tmp;
          tmp = q_2;     q_2 = q_1;     q_1 = q_0;     q_0 = tmp;

          tmp = dtw_2; dtw_2 = dtw_1; dtw_1 = dtw_0; dtw_0 = tmp;
          tmp = dtn_2; dtn_2 = dtn_1; dtn_1 = dtn_0; dtn_0 = tmp;
          tmp = dtU_2; dtU_2 = dtU_1; dtU_1 = dtU_0; dtU_0 = tmp;
          tmp = dtV_2; dtV_2 = dtV_1; dtV_1 = dtV_0; dtV_0 = tmp;
          tmp = dtt_2; dtt_2 = dtt_1; dtt_1 = dtt_0; dtt_0 = tmp;
          tmp = dtq_2; dtq_2 = dtq_1; dtq_1 = dtq_0; dtq_0 = tmp;

          /*********************************************************************/
          /*   Advance time , iter                                             */
          /*********************************************************************/
          p->time +=p->dt;
          BUGREPORT;
          if(iter < 400)          iter++;
      } // End of if (FALSE == first)

      first = FALSE;

      /*********************************************************************/
      /* Static Relationships and Boundaries                               */
      /*********************************************************************/

      BUGREPORT;
      // Extrapolate value on boundary
      FORZY_BD  nbdr[1][iz][ip]   = 0.5*(3.*n_0[iz][ip][nx-1]-n_0[iz][ip][nx-2]);
      FORZY_BD  tbdr[1][iz][ip]   = 0.5*(3.*t_0[iz][ip][nx-1]-t_0[iz][ip][nx-2]);
      PH_Update2dBoundaries(n_0,rbdcndn, nbdr[0], nbdr[1],hval,d);
      PH_Update2dBoundaries(t_0,rbdcndt, tbdr[0], tbdr[1],hval,d);


      Cyto_DParallelSOL(IS_DENSITY,dzN,dzzN,n_0,f_0,exp_n_0,exp_t_0,U_0,dzU,zhval,zgval,d,p);
      FORALL laxn[iz][ip][ir] =   LAXF*dzzN[iz][ip][ir]; // LAXF*(n_0[iz+1][ip][ir] -2.*n_0[iz][ip][ir]+ n_0[iz-1][ip][ir]);

      Cyto_DParallelSOL(IS_TE,dzt,dzzt,t_0,f_0,exp_n_0,exp_t_0,U_0,dzU,zhval,zgval,d,p);
      FORALL laxt[iz][ip][ir] = LAXF*dzzt[iz][ip][ir];    //(t_0[iz+1][ip][ir] -2.*t_0[iz][ip][ir]+ t_0[iz-1][ip][ir]);

      FORALL_BD exp_n_0[iz][ip][ir] = exp(n_0[iz][ip][ir]) ;
      FORALL_BD exp_t_0[iz][ip][ir] = exp(t_0[iz][ip][ir]) ;
#ifndef TE
      FORALL_BD exp_t_0[iz][ip][ir] = 1.;
      FORALL_BD t_0[iz][ip][ir] =  dzt[iz][ip][ir] = 0.;
#endif



      FORALL  nu_density[iz][ip][ir] =  nu_m*exp_n_0[iz][ip][ir]*exp(-3.5*t_0[iz][ip][ir]);
      FORALL  nu_vmu[iz][ip][ir]     =  nu_density[iz][ip][ir] * (V_0[iz][ip][ir]+ p->alpha -0.* U_0[iz][ip][ir]);

      BUGREPORT;
      /*
        Calculate Potential and vorticity from
        density n and GlobalVorticity  W = (Dxx + Dyy) phi + \nabla n \nabla phi

        Iterative solution of equation
        W = (Dxx + Dyy) phi + \nabla n nabla phi
      */
      BUGREPORT;
      PH_Update2dBoundaries(omega,rbdcndw, ombdr[0], ombdr[1],hval,d);
      Laplace_Solve3D(d,p,f_0,omega,fbdr[0],fbdr[1],hval[0],rcor,rbdcndf,fac_pot,FALSE);
      PH_Update2dBoundaries(f_0,rbdcndf, fbdr[0], fbdr[1],hval,d);

      FORALL dxn[iz][ip][ir] = DX(n_0);
      FORALL dyn[iz][ip][ir] = DY(n_0);

      FORALL dxt[iz][ip][ir] = DX(t_0);
      FORALL dyt[iz][ip][ir] = DY(t_0);

      FORALL dxf[iz][ip][ir] = DX(f_0);
      FORALL dyf[iz][ip][ir] = DY(f_0);

      BUGREPORT;
      for(i=0;i<3;i++)
      {
          FORALL  omega[iz][ip][ir] = w_0[iz][ip][ir]- dxn[iz][ip][ir]*dxf[iz][ip][ir]-dyn[iz][ip][ir]*dyf[iz][ip][ir];
          PH_Update2dBoundaries(omega,rbdcndw, ombdr[0], ombdr[1],hval,d);
          Laplace_Solve3D(d,p,f_0,omega,fbdr[0],fbdr[1],hval[0],rcor,rbdcndf,fac_pot,FALSE);
          PH_Update2dBoundaries(f_0,rbdcndf, fbdr[0], fbdr[1],hval,d);
          FORALL dxf[iz][ip][ir] = DX(f_0);
          FORALL dyf[iz][ip][ir] = DY(f_0);

          /*PH_3D_Write(omega,"Vorticity","iter",i,d,p,TRUE);
          PH_3D_Write(w_0,"w_0","iter",i,d,p,FALSE);
          PH_3D_Write(dxn,"dxn","iter",i,d,p,FALSE);
          PH_3D_Write(dyn,"dyn","iter",i,d,p,FALSE);
          PH_3D_Write(dxf,"dxf","iter",i,d,p,FALSE);
          PH_3D_Write(dyf,"dyf","iter",i,d,p,FALSE);
          PH_3D_Write(f_0,"f_0","iter",i,d,p,FALSE);
          PH_3D_Write(n_0,"n_0","iter",i,d,p,FALSE);*/
      }

      // End of iterative solve
      // zero boundaries for dxf, dyf

      if(ispolar){
             PH_Update2dBoundaries(dxf,ZeroNEU, qbdr[0], qbdr[1],hval,d);
             PH_Update2dBoundaries(dyf,ZeroDIR, qbdr[0], qbdr[1],hval,d);
      }
      else{
             PH_Update2dBoundaries(dxf,NEUNEU, qbdr[0], qbdr[1],hval,d);
             PH_Update2dBoundaries(dyf,NEUDIR, qbdr[0], qbdr[1],hval,d);
      }

      BUGREPORT;

      // calculate boundary for w  = om + d_r phi  d_r n
      /* for(iz=0;iz<nz;iz++)
          for(ip=-1;ip<=ny;ip++)
              wbdr[1][iz][ip] =   ombdr[1][iz][ip] +
                  2.*hval[iz][nx-1]* 2.*hval[iz][nx-1]*(f_0[iz][ip][nx] - f_0[iz][ip][nx-1])
                  *(n_0[iz][ip][nx] - n_0[iz][ip][nx-1]) ;
                  */

      BUGREPORT;

      //  V-static       - Equation
      FORALL    vstatic[iz][ip][ir]    =
          (mu*(dzF[iz][ip][ir]- dzN[iz][ip][ir]) + nu_density[iz][ip][ir]*U_0[iz][ip][ir])/(delta +  nu_density[iz][ip][ir]);

#ifdef VSTATIC
      FORALL    V_0[iz][ip][ir]        = vstatic[iz][ip][ir];
#endif

#ifdef QSTATIC
      FORALL    q_0[iz][ip][ir]        = -Kappa/nu_density[iz][ip][ir]*mum1*exp_n_0[iz][ip][ir]*exp_t_0[iz][ip][ir]*exp_t_0[iz][ip][ir]*dzt[iz][ip][ir];
;
#endif

      PH_Update2dBoundaries(w_0,rbdcndw, wbdr[0], wbdr[1],hval,d);
      PH_Update2dBoundaries(U_0,rbdcndu, ubdr[0], ubdr[1],hval,d);
      PH_Update2dBoundaries(V_0,rbdcndv, vbdr[0], vbdr[1],hval,d);
      PH_Update2dBoundaries(t_0,rbdcndt, tbdr[0], tbdr[1],hval,d);
      PH_Update2dBoundaries(q_0,rbdcndq, qbdr[0], qbdr[1],hval,d);


      //  Cyto_DParallelSOL(2000,dzF,dzzF,w_0,f_0,exp_n_0,exp_t_0,U_0,dzU,zhval,zgval,d,p);
      // FORALL laxw[iz][ip][ir] = LAXF*dzzF[iz][ip][ir];//(w_0[iz+1][ip][ir] -2.*w_0[iz][ip][ir]+ w_0[iz-1][ip][ir]);

      // Parallel derivatives
      // Ion parallel derivative first
      Cyto_DParallelSOL(IS_IONVELOCITY,dzU,dzzU,U_0,f_0,exp_n_0,exp_t_0,U_0,dzU,zhval,zgval,d,p);

      Cyto_DParallelSOL(IS_POTENTIAL,dzF,dzzF,f_0,f_0,exp_n_0,exp_t_0,U_0,dzU,zhval,zgval,d,p);

      Cyto_DParallelSOL(IS_HEATFLUX_TE,dzq,dzzq,q_0,f_0,exp_n_0,exp_t_0,U_0,dzU,zhval,zgval,d,p);
      FORALL laxq[iz][ip][ir] = LAXF*dzzq[iz][ip][ir]; // (q_0[iz+1][ip][ir] -2.*q_0[iz][ip][ir]+ q_0[iz-1][ip][ir]);

      Cyto_DParallelSOL(IS_ELECTRONVELOCITY,dzV,dzzV,V_0,f_0,exp_n_0,exp_t_0,U_0,dzU,zhval,zgval,d,p);
      FORALL laxv[iz][ip][ir] = LAXF* dzzV[iz][ip][ir];//  (V_0[iz+1][ip][ir] -2.*V_0[iz][ip][ir]+ V_0[iz-1][ip][ir]);


      //  Calculate Parallel Fluxes
      FORALL_BD VN[iz][ip][ir] = exp_n_0[iz][ip][ir]*V_0[iz][ip][ir];
      FORALL_BD UN[iz][ip][ir] = exp_n_0[iz][ip][ir]*U_0[iz][ip][ir];


      //Cyto_DParallelSOL(IS_TE,dzVN,dzzVN,VN,f_0,exp_n_0,exp_t_0,U_0,dzU,zhval,zgval,d,p);
      //Cyto_DParallelSOL(IS_TE,dzUN,dzzUN,UN,f_0,exp_n_0,exp_t_0,U_0,dzU,zhval,zgval,d,p);
      //FORALL dzVN[iz][ip][ir] /= exp_n_0[iz][ip][ir];
      //FORALL dzUN[iz][ip][ir] /= exp_n_0[iz][ip][ir];



      // zero boundaries for dzzN, dzn
      if(ispolar) PH_Update2dBoundaries(dzzN,ZeroDIR, qbdr[0], qbdr[1],hval,d);
      else PH_Update2dBoundaries(dzzN,NEUDIR, qbdr[0], qbdr[1],hval,d);

      //  Calculate Current
      FORALL_BD J_0[iz][ip][ir] = -exp_n_0[iz][ip][ir]*(V_0[iz][ip][ir]+p->alpha - U_0[iz][ip][ir]);// contains electron bulk drift alpha


      FORALL dzVN[iz][ip][ir] =  dzV[iz][ip][ir] +  (V_0[iz][ip][ir]+p->alpha)*dzN[iz][ip][ir]; // contains electron bulk drift alpha
      FORALL dzUN[iz][ip][ir] =  dzU[iz][ip][ir] +  U_0[iz][ip][ir]*dzN[iz][ip][ir];


      COMM(fprintf(stderr,"t = %f, wt = %f, ot = %d, otm = %d\n",p->time,write_time,ot,p->otmult););
      BUGREPORT;

      /*******************************************/
      /*                                         */
      /*                                         */
      /*      Output                             */
      /*                                         */
      /*                                         */
      /*******************************************/

      //    Write initial condition
      if(iter < 10)
      {
          MPI_Barrier(MPI_COMM_WORLD);
          PH_3D_Write(omega,"Vorticity","ini",iter,d,p,TRUE);
          PH_3D_Write(n_0,"lnn","ini",iter,d,p,FALSE);
          PH_3D_Write(exp_n_0,"n","ini",iter,d,p,FALSE);
          PH_3D_Write(f_0,"Potential","ini",iter,d,p,FALSE);
          PH_3D_Write(U_0,"UIon","ini",iter,d,p,FALSE);
          PH_3D_Write(t_0,"lnTe","ini",iter,d,p,FALSE);
          PH_3D_Write(exp_t_0,"Te","ini",iter,d,p,FALSE);
          PH_3D_Write(V_0,"Velocity","ini",iter,d,p,FALSE);
          PH_3D_Write(J_0,"Current","ini",iter,d,p,FALSE);
          PH_3D_Write(target,"Target","ini",iter,d,p,FALSE);
          PH_3D_Write(dsource,"Source","ini",iter,d,p,FALSE);
          PH_3D_Write(dzU,"dzU","ini",iter,d,p,FALSE);
          PH_3D_Write(dzN,"dzN","ini",iter,d,p,FALSE);
          PH_3D_Write(dzF,"dzF","ini",iter,d,p,FALSE);
          PH_3D_Write(dzt,"dzt","ini",iter,d,p,FALSE);
          PH_3D_Write(w_0,"H","ini",iter,d,p,FALSE);
      }

      if(p->time>= write_time-0.5*p->dt)
      {
          BUGREPORT;
          ot++;
          write_time+= p->out_time;

          // Check the CFL condition
          if(  PH_CheckCFL3D(f_0,vr,vp,d,p,hval,vval,&cflr,&cflp) == -1)
          {
              PH_3D_Write(vp,"vp","VPR",0,d,p,TRUE);
              PH_3D_Write(vr,"vr","VPR",0,d,p,FALSE);
              PH_3D_Write(f_0,"Pot","VPR",0,d,p,FALSE);

              fprintf(stderr,"CFL violation detected! Aborting.\n");
              exit(0);
          }
          if(isnan(cflr*cflp))
          {
              fprintf(stderr,"NaNs detected! Aborting.\n");
              exit(0);
          }


          /*
             Integrals
          */
          /**************************************************************/

          // Energy
          FORALL dtw_0[iz][ip][ir] = -f_0[iz][ip][ir]*omega[iz][ip][ir];
          reducing[0] =  PH_3DIntegral_noreduce(dtw_0,1,norm,nz,ny,nx);   // p->energy
          reducing[1] =  PH_3DIntegral_noreduce(f_0,1,norm,nz,ny,nx);     // totalf
          reducing[2] =  PH_3DIntegral_noreduce(omega,2,norm,nz,ny,nx);   // p->vorticity
          reducing[3] =  PH_3DIntegral_noreduce(omega,1,norm,nz,ny,nx);   // p->circulation
          reducing[4] =  PH_3DIntegral_noreduce(exp_n_0,1,norm,nz,ny,nx); // gamma_n




          if(isnan(reducing[0]))
          {
              fprintf(stderr,"NaNs detected! Aborting.\n");
              exit(0);
          }


          BUGREPORT;
          // Calculate current
          reducing[5] =  PH_3DIntegral_noreduce(J_0,1,norm,nz,ny,nx);

          // Calculate parallel density transport to plates
          FORALL tfeld[iz][ip][ir] = exp_n_0[iz][ip][ir]*(V_0[iz][ip][ir]+p->alpha);
          if(d->grid_coords[0] == 0 ){
              Util_Integral(tfeld[0],1,norm[0],ny,nx,&reducing[6]);
              reducing[6] *= -1.;
          }
          else if(d->grid_coords[0] ==  (d->N[0]-1))
              Util_Integral(tfeld[nz-1],1,norm[nz-1],ny,nx,&reducing[6]);
          else reducing[6] = 0.;

          // current loss
          if(d->grid_coords[0] == 0 ) {
              Util_Integral(J_0[0],1, norm[0],ny,nx,&reducing[7]);
              reducing[7] *= -1.;
          }
          else if(d->grid_coords[0] ==  (d->N[0]-1))
              Util_Integral(J_0[nz-1],1, norm[nz-1],ny,nx,&reducing[7]);
          else reducing[7] = 0.;


          FORALL tfeld[iz][ip][ir] = p->nprof*dsource[iz][ip][ir];
          reducing[8] = PH_3DIntegral_noreduce(tfeld,1,norm,nz,ny,nx);

          BUGREPORT;
          MPI_Allreduce((void *)&reducing[0], (void *) &result_vector[0], 9, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

          p->energy               = result_vector[0];
          totalf                  = result_vector[1];
          p->vorticity            = result_vector[2];
          p->circulation          = result_vector[3];
          gamma_n                 = result_vector[4];
          Jtotal                  = result_vector[5];
          G_par                   = result_vector[6];
          J_par                   = result_vector[7];
          N_in                    = result_vector[8];

          BUGREPORT;

          if(ErhFirstWrite && ISROOT)
          {

              output = fopen(d->erhname,"a");
              fprintf(output,"# nx = %d ny = %d nz = %d\n",nx,ny,d->dims[0]);
              fprintf(output,"# Gamma %f; sigma %f; nu %f; alpha %f\n",p->gamma,p->sigma,nu,p->alpha);
              fprintf(output,"# mue_w %f; mue_n %f; kappan %f\n",p->mue_w,p->mue_n,p->kappan);
              fprintf(output,"# source %f; k0 %f\n#\n",p->source,p->k0);
              fprintf(output,"# 1:t\t 2:energy\t 3:vorticity\t 4:<n>"
                      "\t 5: <f>\t 6:circulation\t 7:cflr\t 8:cflp"
                      "\t 9:k_eff \t 10:J \t 11:G_par \t 12:J_par"
                      "\t 13: impulse \t 14:nloss \t 15:ipol \t 16:dtimpulse\n");
              fclose(output);


             output = fopen(profilename,"a");
             fprintf(output,"# nx = %d ny = %d nz = %d\n",nx,ny,d->dims[0]);
              fprintf(output,"# Gamma %f; sigma %f; nu %f; alpha %f\n",p->gamma,p->sigma,nu,p->alpha);
              fprintf(output,"# mue_w %f; mue_n %f; kappan %f\n",p->mue_w,p->mue_n,p->kappan);
              fprintf(output,"# source %f; k0 %f;\n#\n",p->source,p->k0);
              fprintf(output,"# 1: r\t 2: momentum \t3: vpar  \t 4: upar \t 5: n \t 6: phi \t 7: w\n");
              fclose(output);

              ErhFirstWrite=FALSE;
          }


          if(ISROOT)
          {
              output = fopen(d->erhname,"a");
              fprintf(output,"%g \t %g \t% g \t %g \t",
                      p->time,p->energy,p->vorticity,gamma_n);
              fprintf(output,"%g \t %g \t% g \t %g \t ",
                      totalf,p->circulation,cflr,cflp);
              fprintf(output,"%g \t %g \t %g \t %g\n",
                      sqrt(p->vorticity/p->energy),Jtotal,G_par, N_in);
              fclose(output);

          }
          //Momentum profile
          snprintf(profilename,DEFSTRLEN,  "%s.parallel_momentum.prof",d->name_out);
          FORALL dtw_0[iz][ip][ir] = exp_n_0[iz][ip][ir]*(ME*V_0[iz][ip][ir]+p->Mi*U_0[iz][ip][ir]);
          PH_3D_n0m0_Component(dtw_0,ppar_prof,d);
          PH_write_radial_profile(ppar_prof,profilename,d);

          //V-electron profile
          snprintf(profilename,DEFSTRLEN,  "%s.uelectron.prof",d->name_out);
          PH_3D_n0m0_Component(V_0,vpar_prof,d);
          PH_write_radial_profile(vpar_prof,profilename,d);

          //U-ion Profile
          snprintf(profilename,DEFSTRLEN,  "%s.uion.prof",d->name_out);
          PH_3D_n0m0_Component(U_0,upar_prof,d);
          PH_write_radial_profile(upar_prof,profilename,d);

          //density Profile
          snprintf(profilename,DEFSTRLEN,  "%s.density.prof",d->name_out);
          PH_3D_n0m0_Component(n_0,npar_prof,d);
          PH_write_radial_profile(npar_prof,profilename,d);

          //Potential Profile
          snprintf(profilename,DEFSTRLEN,  "%s.potential.prof",d->name_out);
          PH_3D_n0m0_Component(f_0,fpar_prof,d);
          PH_write_radial_profile(fpar_prof,profilename,d);

          //Vorticity Profile
          snprintf(profilename,DEFSTRLEN,  "%s.vorticity.prof",d->name_out);
          PH_3D_n0m0_Component(w_0,wpar_prof,d);
          PH_write_radial_profile(wpar_prof,profilename,d);
              BUGREPORT;
              MPI_Barrier(MPI_COMM_WORLD);


          // Axial profiles

          //density Profile
          snprintf(profilename,DEFSTRLEN,  "%s.n_z",d->name_out);
          FORZ nax_prof[iz] = n_0[iz][0][0];
          PH_write_axial_profile(nax_prof,profilename,d);


          //velocity Profile
          snprintf(profilename,DEFSTRLEN,  "%s.ve_z",d->name_out);
          FORZ vax_prof[iz] = V_0[iz][0][0];
          PH_write_axial_profile(vax_prof,profilename,d);


          //velocity Profile
          snprintf(profilename,DEFSTRLEN,  "%s.vi_z",d->name_out);
          FORZ uax_prof[iz] = U_0[iz][0][0];
          PH_write_axial_profile(uax_prof,profilename,d);


          /**************************************************************/
          /*
             Output 3D
          */

          /**************************************************************/
         /**************************************************************/
          if((ot >= p->otmult) /*&& p->time < 2000*/)
          {
              ot = 0;

              BUGREPORT;

              PH_3D_Write(exp_n_0,"n",d->name_out,d->number,d,p,TRUE);
              if(W_UION)PH_3D_Write(U_0,"UIon",d->name_out,d->number,d,p,FALSE);
              if(W_DZN) PH_3D_Write(dzN,"dzn",d->name_out,d->number,d,p,FALSE);
              if(W_DZN) PH_3D_Write(dzzN,"dzzn",d->name_out,d->number,d,p,FALSE);
              if(W_H0)  PH_3D_Write(w_0,"H",d->name_out,d->number,d,p,FALSE);
              if(W_F0)  PH_3D_Write(f_0,"Potential",d->name_out,d->number,d,p,FALSE);
              if(W_lnN0)  PH_3D_Write(n_0,"N",d->name_out,d->number,d,p,FALSE);
              if(W_V0)  PH_3D_Write(V_0,"VElectron",d->name_out,d->number,d,p,FALSE);
              if(W_W0)  PH_3D_Write(omega,"Vorticity",d->name_out,d->number,d,p,FALSE);
              if(W_DZF) PH_3D_Write(dzF ,"dzf",d->name_out,d->number,d,p,FALSE);
              if(W_DZF) PH_3D_Write(dzzF ,"dzzf",d->name_out,d->number,d,p,FALSE);
              if(W_VSTATIC) PH_3D_Write(vstatic,"Vstatic",d->name_out,d->number,d,p,FALSE);
              if(W_J0) PH_3D_Write(J_0,"Current",d->name_out,d->number,d,p,FALSE);

              //FORALL res[iz][ip][ir] = dzF[iz][ip][ir]-dzN[iz][ip][ir];
              //PH_3D_Write(res,"ParForce",d->name_out,d->number,d,p,FALSE);

              FORALL tfeld[iz][ip][ir] =  exp_n_0[iz][ip][ir]*V_0[iz][ip][ir];
              if(W_PFLUX) PH_3D_Write(tfeld,"ParticleFlux",d->name_out,d->number,d,p,FALSE);

              PH_3D_Write(dzVN,"dzVN",d->name_out,d->number,d,p,FALSE);

              PH_3D_Write(VN,"VN",d->name_out,d->number,d,p,FALSE);

              if(W_DZVN) PH_3D_Write(dzV,"dzV",d->name_out,d->number,d,p,FALSE);
              if(W_DZUN)  PH_3D_Write(dzU,"dzU",d->name_out,d->number,d,p,FALSE);
              if(W_NUVMU) PH_3D_Write(nu_vmu,"nu_vmu",d->name_out,d->number,d,p,FALSE);

              //FORALL tfeld[iz][ip][ir] =   dxn[iz][ip][ir]*dxf[iz][ip][ir]+ dyn[iz][ip][ir]*dyf[iz][ip][ir];
              //PH_3D_Write(tfeld,"GradN_GradPhi",d->name_out,d->number,d,p,FALSE);


              //FORALL tfeld[iz][ip][ir] = exp_n_0[iz][ip][ir]*(ME*V_0[iz][ip][ir]+p->Mi*U_0[iz][ip][ir]);
              //PH_3D_Write(tfeld,"ParMomentum",d->name_out,d->number,d,p,FALSE);
#ifdef TE
              if(W_TE0) PH_3D_Write(exp_t_0,"Te",d->name_out,d->number,d,p,FALSE);
              if (W_lnT0) PH_3D_Write(t_0,"lnTe",d->name_out,d->number,d,p,FALSE);
              if(W_Q0) PH_3D_Write(q_0,"Heatflux",d->name_out,d->number,d,p,FALSE);
              if(W_DZQ) PH_3D_Write(dzq,"dzq",d->name_out,d->number,d,p,FALSE);
#endif //Te


              if(ISROOT) fprintf(stderr,"\rWrote %s.%03d  J = %f Jtotal %f G = %f  N_in = %f         ",d->name_out,d->number,J_par,Jtotal,G_par, N_in);
              d->number++;
              BUGREPORT;
          }
      }
      BUGREPORT;
  }
 BUGREPORT;

 exit(0);
}

/*************************************************************************************/
void Cyto_DParallelSOL(int identity,double ***result,double ***result2,
                       double ***val,double ***f,double ***exp_n,double ***exp_t,
                       double ***v,double ***dzv,double *zhval, double *zgval,
                       HDF_DS *data, PARA *para)
{
  /*!   This routine calculates the parallel derivative:
       and takes care of the parallel SOL boundary condition

       As parameters it gets:
       int identity : kind of field to apply BD conditions right
       result: Stores the result on exit
       val   : f
       f     : Potential for sheath condition
       lf    : factor for the linear part of the derivative

       v     : parallel ion velocity
       dzv   : parallel gradient of ion velocity

       zhval : geometry for derivative
       zgval : curvature of geometry

       at lower end all parallel boundaries are neumann with zero derivative
  */



    static int FIRST= TRUE;
    static int nx,ny,nz,offz;
    static double lf_local,lf_local2;
    register int iz,ip,ir;
    double bdval,tval;

    if(FIRST)
      {
          nx = data->lnx;
          ny = data->lny;
          nz = data->lnz;

          offz = data->offz;
          COMM(fprintf(stderr,"%d(%ld) %d(%ld) %d(%ld) \n",nx,data->offx,ny,data->offy,nz,offz);     );
          FIRST = FALSE;
      }

    BUGREPORT;

    // Communicate Ghostpoints
    PH_UpdateZ(val,data,MPI_COMM_WORLD);



    BUGREPORT;
    if(offz == 1)
    {
        if(data->grid_coords[0] == 0)
        {
            if(para->zmin < 0.0)
            {
              // a wall for zmin  negative
                BUGREPORT;
                switch(identity){
                    case IS_CURRENT:
                        FORYX_BD
                        {
                            //Velocity is negative as electrons leave the plasma
                            // if f=0 the current is zero
                            bdval =    -2.*sqrt(exp_t[0][ip][ir])*( exp_n[0][ip][ir]*(1.-exp(-f[0][ip][ir]) ));
                            val[-1][ip][ir] =  bdval-val[0][ip][ir];
                        }
                    break;
                    case  IS_ELECTRONVELOCITY:
                        FORYX_BD
                        // Effective Electron speed is given by the sheath potential
                        {
                            // Velocity is negative as electrons leave the plasma
                            // if f=0 the velocity is - unity
                            bdval =    -sqrt(exp_t[0][ip][ir])* exp(-f[0][ip][ir]/exp_t[0][ip][ir]);
                            val[-1][ip][ir] =  2.*bdval-val[0][ip][ir];
                        }
                        break;
                    case  IS_IONVELOCITY:
                        // Ion speed is c_s ~ sqrt(Te) into the limiter
                        FORYX_BD
                        {
                            // bdval=-2.*sqrt(exp_t[0][ip][ir]);
                            val[-1][ip][ir] =1.;// bdval-val[0][ip][ir];
                        }
                        break;
                    case IS_DENSITY:
                        // density flux remains constant at the sheath that is dz N = -1/V dz V
                        FORYX_BD
                        {
                            val[-1][ip][ir] =   val[0][ip][ir];
                        }
                        break;
                    case IS_POTENTIAL:
                        // No gradient in potential, sheeth works via parallel current
                        FORYX_BD val[-1][ip][ir] =  val[0][ip][ir];

                        // continous gradient
                        // FORYX_BD  val[-1][ip][ir] =  -val[1][ip][ir]+2.*val[0][ip][ir];
                        break;
                    case IS_TEMPERATURE:
                    case IS_TE:  // !!Here use sheath transparancy formula
                        FORYX_BD
                        {
                            tval = 0.;
                            //val[-1][ip][ir] =   tval + val[0][ip][ir];
                            val[-1][ip][ir] =  -val[1][ip][ir]+(1.+zhval[0]/zhval[1])*val[0][ip][ir];
                        }
                        break;
                    case IS_HEATFLUX_TE:
                        FORYX_BD
                        {
                            bdval =  exp_n[-1][ip][ir]*(2.*exp_t[-1][ip][ir]+f[-1][ip][ir])*
                                ( -sqrt(exp_t[-1][ip][ir])* exp(-f[0][ip][ir]/exp_t[-1][ip][ir]));
                               //tval = 0.;
                            val[-1][ip][ir] =   2.*bdval - val[0][ip][ir];
                            //val[-1][ip][ir] =  -val[1][ip][ir]+(1.+zhval[0]/zhval[1] ) *val[0][ip][ir];
                        }
                        break;

                    default:
                        // No gradient
                        FORYX_BD val[-1][ip][ir] =  val[0][ip][ir];
                        break;

                }
            }
            else // zmin >= 0 , case when we do not have boundary at lower end
            {
                // Lower boundary, start of plasma, a nonconducting wall
                BUGREPORT;
                switch(identity){
                    case IS_CURRENT: // Zero current
                        FORYX_BD val[-1][ip][ir] = -val[0][ip][ir];
                        break;
                    case  IS_ELECTRONVELOCITY: // Zero velocities
                        FORYX_BD                  val[-1][ip][ir] = -val[0][ip][ir];
                        break;
                    case  IS_IONVELOCITY: // Zero velocities
                        FORYX_BD                  val[-1][ip][ir] = -val[0][ip][ir];
                        break;
                    case IS_DENSITY: // no gradient in density
                        FORYX_BD                  val[-1][ip][ir] =   val[0][ip][ir];
                        break;
                    case IS_POTENTIAL: //no gradient in potential
                        FORYX_BD                  val[-1][ip][ir] =  val[0][ip][ir];
                        break;
                    case IS_TEMPERATURE: // no gradient in temperature
                    case IS_TE:
                        FORYX_BD                  val[-1][ip][ir] =  val[0][ip][ir];
                        break;
                    case IS_HEATFLUX_TE: // zero heatflux
                    default:
                        FORYX_BD                  val[-1][ip][ir] =  -val[0][ip][ir];
                        break;

                }
            }

        } // gridcoords = 0

        BUGREPORT;



        // Upper z boundary, end of plasma, the sheath
        if(data->grid_coords[0] == (data->N[0]-1) ) // Last proc
        {
            BUGREPORT;
            switch(identity){
                case IS_CURRENT:
                    FORYX_BD
                    {
                        // Velocity is positive as electrons leave the plasma
                        // if f=0 the current is zero
                        bdval =    2.*( exp_n[nz-1][ip][ir]*(1.-exp(-f[nz-1][ip][ir]) ));
                        val[nz][ip][ir] =  bdval-val[nz-1][ip][ir];
                    }
                    break;
                case  IS_ELECTRONVELOCITY:
                    // Effective Electron speed is given by the sheath potential
                    // Velocity is positive as electrons leave the plasma
                    // if f=0 the velocity is unity, a positive potential deccellerates electrons
                    tval = 2./(1.+ 4.*zgval[nz-1]/zhval[nz-1]/zhval[nz-1]);
                    FORYX_BD
                    {
                        val[nz][ip][ir] =  tval*f[nz-1][ip][ir]  + (1.-tval)*f[nz-2][ip][ir];//extrapolate phi

                        val[nz][ip][ir] =  sqrt(exp_t[nz-1][ip][ir])*exp(-val[nz][ip][ir]/exp_t[nz-1][ip][ir]); // set velocity variable directly


                        // bdval =  (-f[nz-2][ip][ir]+(1.+zhval[nz]/zhval[nz-1] )*f[nz-1][ip][ir]+f[nz-1][ip][ir])*0.5;
                        // bdval =   sqrt(exp_t[nz-1][ip][ir])*exp(-bdval/exp_t[nz-1][ip][ir]);
                        // bdval =   sqrt(exp_t[nz-1][ip][ir])*exp(-f[nz-1][ip][ir]/exp_t[nz-1][ip][ir]);
                        // val[nz][ip][ir] =  2.*bdval-val[nz-1][ip][ir];

                    }
                    break;
                case IS_IONVELOCITY:
                    // Ion speed is cs ~ sqrt(Te) into the limiter
                    FORYX_BD
                    {
                        bdval           =    sqrt(exp_t[nz-1][ip][ir]);
                        val[nz][ip][ir] = 2.* bdval - 0.5*(val[nz-1][ip][ir]+val[nz-2][ip][ir]);
                    }
                    break;
                case IS_DENSITY:
                    // density flux remains constant at the sheath that is dz N = -dz U
                    FORYX_BD   val[nz][ip][ir] = -0.01*(30.*dzv[nz-1][ip][ir]+30.*dzv[nz-1][ip][ir]+20.*dzv[nz-1][ip][ir]+20.*dzv[nz-1][ip][ir])/(0.5*zhval[nz-1]) + val[nz-2][ip][ir];
                    break;
                case IS_POTENTIAL:
                    // No gradient in potential, sheeth works via parallel current
                    // FORYX_BD  val[nz][ip][ir] =  val[nz-1][ip][ir];

                    // continous gradient
                    FORYX_BD   val[nz][ip][ir] = -0.1*(4.*dzv[nz-1][ip][ir]+3.*dzv[nz-1][ip][ir]+2.*dzv[nz-1][ip][ir]+dzv[nz-1][ip][ir])/(0.5*zhval[nz-1]) + val[nz-2][ip][ir];
                    break;
                case IS_TEMPERATURE:
                case IS_TE:
                    FORYX_BD
                    {
                        tval = 0.;
                        //val[nz][ip][ir] =   tval + val[nz-1][ip][ir];
                        val[nz][ip][ir] =  -val[nz-2][ip][ir]+(1.+zhval[nz]/zhval[nz-1] )*val[nz-1][ip][ir];
                    }
                    break;
                case IS_HEATFLUX_TE:
                    FORYX_BD
                    {
                        //bdval =  exp_n[nz][ip][ir]*(exp_t[nz][ip][ir]+f[nz][ip][ir]/exp_t[nz][ip][ir])*         ( sqrt(exp_t[nz][ip][ir])* exp(-f[nz][ip][ir]/exp_t[nz][ip][ir]));
                        //bdval =0.*5.5*exp_n[nz-1][ip][ir]*exp_t[nz-1][ip][ir];
                        //val[nz][ip][ir] =   2.*bdval - val[nz-1][ip][ir];
                        //val[nz][ip][ir] =  -val[nz-2][ip][ir]+(1.+zhval[nz]/zhval[nz-1] )*val[nz-1][ip][ir];
                        val[nz][ip][ir] =  val[nz-1][ip][ir];
                    }
                    break;
               default:
                   // No gradient
                    FORYX_BD              val[nz][ip][ir] =  val[nz-1][ip][ir];
                    break;
            }
        } // Lastproc

        BUGREPORT;
        // dz
        FORALL_BD_XY result[iz][ip][ir]  =  0.5*zhval[iz]*(val[iz+1][ip][ir] - val[iz-1][ip][ir]);

        BUGREPORT;

        // dzz
        FORALL_BD_XY result2[iz][ip][ir]  =  0.25*zhval[iz]*zhval[iz] *((val[iz+1][ip][ir] + val[iz-1][ip][ir]) -2.* val[iz][ip][ir])
            +zgval[iz]*(val[iz+1][ip][ir] - val[iz-1][ip][ir]);
        BUGREPORT;
    }
    else /* Two boundary points, higher accuracy, but not implemented in parallel, boundary exchange is neither */
    {

        switch(identity){
            case IS_CURRENT:
                /* Velocity is positive as electrons leave the plasma */
                /* if f=0 the velocity is 1 */
                FORYX_BD
                {
                    bdval =    2.*( exp_n[nz][ip][ir]*(1.-exp(-f[nz][ip][ir])));
                    val[nz][ip][ir] =  bdval-val[nz-1][ip][ir];
                    val[nz+1][ip][ir] =  bdval-val[nz-2][ip][ir];

                    val[-1][ip][ir] = -val[0][ip][ir];
                    val[-2][ip][ir] = -val[1][ip][ir];
                }
                break;
            case  IS_ELECTRONVELOCITY:
                FORYX_BD
                    /* Effective Electron speed is given by the sheath potential */
                    {
                        bdval = 2.*exp(-f[nz][ip][ir]);
                        val[nz][ip][ir] = bdval-val[nz-1][ip][ir];
                        val[nz+1][ip][ir] =  bdval-val[nz-2][ip][ir];

                        //boundary value is zero at upper vall, no flux

                        val[-1][ip][ir] =  -val[0][ip][ir];
                        val[-2][ip][ir] = -val[1][ip][ir];

                  }
                break;
            case  IS_IONVELOCITY:
                /* Ion speed is unity into the limiter */
                FORYX_BD
                    {
                        bdval=2.;
                        val[nz][ip][ir] =   bdval-val[nz-1][ip][ir];
                        val[nz+1][ip][ir] =  bdval-val[nz-2][ip][ir];


                        val[-1][ip][ir] = - val[0][ip][ir];
                        val[-2][ip][ir] = - val[1][ip][ir];

                    }
                break;
            case IS_DENSITY:
                /* density flux remains constant at the sheath dz N = */
                FORYX_BD
                {
/*                    tval = (exp(-f[nz][ip][ir])-v[nz-1][ip][ir])/para->dz;
                    tval *=   -0.5/MAX(exp(-f[nz][ip][ir]),0.4)/(lf_local);



                    val[nz][ip][ir] =   tval + val[nz-1][ip][ir];
                    val[nz+1][ip][ir] =  -6.*tval +  8.*val[nz][ip][ir] -8.*val[nz-1][ip][ir] + val[nz-2][ip][ir];
*/
                    val[nz][ip][ir] =  val[nz-1][ip][ir];
                    val[nz+1][ip][ir] =  val[nz-2][ip][ir];


                    val[-1][ip][ir] =  val[0][ip][ir];
                    val[-2][ip][ir] =  val[1][ip][ir];
                }
                break;
            case IS_POTENTIAL:
                /* No gradient in potential, sheeth works via parallel current */
                FORYX_BD
                    {
                        val[nz][ip][ir] =  val[nz-1][ip][ir];
                        val[nz+1][ip][ir] =  val[nz-2][ip][ir];

                        val[-1][ip][ir] =  val[0][ip][ir];
                        val[-2][ip][ir] =  val[1][ip][ir];
                    }
            break;
            case IS_TEMPERATURE:
                break;
            case IS_TE:
            case IS_HEATFLUX_TE:
            default:
                fprintf(stderr,"unidentified field object (UFO) in  boundary routine!\n");
        }


        // dz
        FORALL_BD_XY
            result[iz][ip][ir]  =  1./6.*zhval[iz]*(-val[iz+2][ip][ir] + 8.*(val[iz+1][ip][ir] - val[iz-1][ip][ir]) + val[iz-2][ip][ir]);

        // dzz
        FORALL_BD_XY
            result2[iz][ip][ir]  =  1./3.*zhval[iz]*zhval[iz] *(-val[iz+2][ip][ir] + 16.*(val[iz+1][ip][ir] + val[iz-1][ip][ir]) -30.* val[iz][ip][ir] - val[iz-2][ip][ir]);
    }

    BUGREPORT;
}


/****************************************************************************/
void Cyto_Laplace(double ***w,double ***f,double **hval,double **vval,
                       double *edr,int nx,int ny,int nz,PARA *p,int ispolar)
{
/*! This function calculates the laplace of a scalar field, given the


     w = Lap f =  d_rr f + 1/r d_r f + 1./r2 d_ff f

     The laplacian is written out in field w,
     all other fields are input

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
        /*
          for(j=0;j<ny;j++) for(k=0;k<nx;k++) w[iz][ip][ir]= f[iz][ip][ir];
          //BUGREPORT;
          Fft_1d_2d_f(w[i],nx,ny);
          //BUGREPORT;
          for(j=0;j<ny;j++)
          {
          ky = (j/2)*p->dky;
          ky*=-ky;
          for(k=0;k<nx;k++) w[i][j][k] *= ky*edr[k]*edr[k];
          }
          //BUGREPORT;
          Fft_1d_2d_b(w[i],nx,ny);
          //BUGREPORT;

        */

        // If complete FD: Uncomment below

        for(j=0;j<ny;j++) for(k=0;k<nx;k++)   w[i][j][k] =  vval[i][k]*vval[i][k]*(f[i][j+1][k] -2.*f[i][j][k]  +f[i][j-1][k]);
        for(j=0;j<ny;j++) for(k=0;k<nx;k++)  w[i][j][k] += hval[i][k]*hval[i][k]*(f[i][j][k+1] -2.*f[i][j][k]  +f[i][j][k-1]);
        if(ispolar)
        {

            for(j=0;j<ny;j++) for(k=0;k<nx;k++)  w[i][j][k] += 0.5*edr[k]*hval[i][k]*(f[i][j][k+1] - f[i][j][k-1]);
        }

  }
}

/*! This routine calculates secondary parameters from primary ones */
void Cyto_CalcParameters(PARA *p,HDF_DS *d)
{
        double rho_s_cm;
    double val;
    double lgamma = 2.;
    double Bgauss=0.;
    double lamda_ei =0.;
    double n_cgs;
    int i;


    Util_CalcParameters(p,d);

    /* Calculate equation parameters */
    /* returns SI units, i.e. rho_s is in [m] */

    p->unit_time = p->rho_s/p->c_s;
    p->unit_length = p->rho_s;
        rho_s_cm = p->rho_s *100.;

    fprintf(stderr,"Parameters:\n");
    fprintf(stderr,"\t%-20s = %f [mu_s]\n","Unit of time",p->unit_time*1e6);
    fprintf(stderr,"\t%-20s = %f [cm]\n\n","Unit of length",p->unit_length*100);


    if(p->nu < DBL_MIN) p->nu = p->nu_ei/p->omega_ci;
    if(DBL_MIN > p->mue_w) p->mue_w = 0.3*(p->Ti/p->Te)*p->nu_in/p->omega_ci;
    if(DBL_MIN > p->mue_n) p->mue_n = p->mue_w ;


    if(0. == p->kappan)
    {
        if (p->ln > 0.) p->kappan = p->rho_s/p->ln;
    }
    if(0. == p->sigma) p->sigma = p->nu_in/p->omega_ci;
    if(0. == p->delta) p->delta = p->nu_en/p->omega_ci;


    // sources
    // Calculate normalised source strength
    // part_source is in particles per second per cubic meter


    p->nprof = p->part_source/p->n0*p->unit_time;
    p->tprof = p->temp_source/p->n0/p->Te*p->unit_time;

    fprintf(stderr,"\t%-20s = %f \n","nu",p->nu);
    fprintf(stderr,"\t%-20s = %f \n","mue_w",p->mue_w);
    fprintf(stderr,"\t%-20s = %f \n","mue_n",p->mue_n);
    fprintf(stderr,"\t%-20s = %f \n","delta(nu_en)",p->delta);
    fprintf(stderr,"\t%-20s = %f \n","kappan",p->kappan);
    fprintf(stderr,"\t%-20s = %f \n","sigma(nu_in)",p->sigma);
    fprintf(stderr,"\t%-20s = %f \n","particle source",p->nprof);
    fprintf(stderr,"\t%-20s = %f \n","energy source",p->tprof);
    /* Change length scales */

    p->xmin = 0.;
    p->xmax = p->a/p->rho_s;

    p->ymin = 0.;
    p->ymax = 2.*M_PI;

    p->zmin = 0.;
    p->zmax = p->lpar/p->rho_s;



    fprintf(stderr,"Domain:\n");
    fprintf(stderr,"\t [ %f: %f ] x [ %f : %f ] x [ %f : %f ] with  (%d,%d,%d)\n\n",
            p->xmin,p->xmax,p->ymin,p->ymax,p->zmin,p->zmax, (int)d->nx,(int)d->ny,(int)d->nz);

}

/*!Function maps spectral noise onto poloidal domain */
void Map_Noise(HDF_DS *data, PARA *p)
{
int
    ix,
    iy,
    iz,
    i,
    j,
    k,
    nor = 64;

double
    r,
    phi,
    x,
    y,
    dxp,
    dyp,
    cs,
    ss,
    mod,
        kx,
    ky,
    r0,
    theta0 = 0.,
        dkx,
    dky,
        kx2,
    ky2,
    khx,
    khy,
        fac,
    phase,
    amp,
    norm;


 double
     *Sx,
     *Sy,
     *xpos,
     *ypos,
     *vals;

 double
     **feld;


 Sx = (double *)calloc( (size_t)(nor+10),sizeof(double) );
 Sy = (double *)calloc( (size_t)(nor+10),sizeof(double) );

 xpos = (double *)calloc( (size_t)(data->dims[2]+2),sizeof(double) );
 ypos = (double *)calloc( (size_t)(data->dims[2]+2),sizeof(double) );
 vals = (double *)calloc( (size_t)(data->dims[2]+2),sizeof(double) );


 feld = Util_DMatrix(nor,0,nor,0);
 dkx = M_PI;
 dky = M_PI;

 for(i=0,kx=0.;i<nor;i+=2,kx+=dkx) Sx[i] = kx;
 for(i=0,ky=-.5*(double)nor*dky;i<nor;i++, ky+=dky) Sy[i] = ky;

 /* Map rectangle  to circle */


 khx = 1./(p->width_random_x[i]*p->width_random_x[i]*dkx*dkx);
 khy = 1./(p->width_random_y[i]*p->width_random_y[i]*dky*dky);


 /*
    Define field in fourier space

 */
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

 /*
    No m = 0 component
 */
 for (iy=0;iy<nor;iy++)
     feld[iy][0]= feld[iy][1] = 0.;



 /*
    Interpolate to field
 */
 for (iy=0;iy<data->dims[1];iy++)
 {
     phi = data->coordinate[1][iy];
     cs = cos(phi);
     ss = sin(phi);

     for (ix=0;ix<data->dims[2];ix++)
     {
         r   = data->coordinate[2][ix]/p->xmax;
         xpos[ix] = r*cs;
         ypos[ix] = r*ss;
     }

     spectral_interpol(data->datarw,xpos,ypos,data->nx,feld,nor,nor,Sx,Sy);

 }


 /*
    smooth to outer boundary
 */

 for (iy=0;iy<data->dims[1];iy++)
         for (ix=0;ix<data->dims[2];ix++)
         {
             r   = (data->coordinate[2][ix]-p->xmin)/(p->xmax-p->xmin);
             //   data->dddfelder[0][iy][ix] *= tanh((1.- r)*10.);
         }


 /*
   smooth to inner boundar if annulus only
 */

 if(p->xmin != 0.)
     for (iy=0;iy<data->ny;iy++)
         for (ix=0;ix<data->nx;ix++)
         {
             r   = (data->coordinate[2][ix]-p->xmin)/(p->xmax-p->xmin);
             // data->dddfelder[0][iy][ix] *= tanh(r*10.);
         }


 /*
   do z -dependence
 */

 for(iz=data->dims[0]-1;iz>=0;iz--)
 {
     mod =  (1. + .1*((double)rand()/(double)RAND_MAX -0.5));
     //   for (iy=0;iy<data->dims[1];iy++)
     //  for (ix=0;ix<data->dims[2];ix++)
             //  data->dddfelder[iz][iy][ix] =  mod*data->dddfelder[0][iy][ix];
 }

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
double val;

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

        }
    }
}

}
