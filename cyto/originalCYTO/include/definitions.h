/*! \file definitions.h */

#ifndef         DEFINITIONS_H_
#define         DEFINITIONS_H_


#define true  1
#define false 0

#define NEW_HDF 1
#define APPEND_HDF 0



#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923
#endif
#ifndef M_1_PI
#define M_1_PI      0.31830988618379067154
#endif
#ifndef M_PI_4
#define M_PI_4      0.78539816339744830962
#endif
#ifndef M_2_PI
#define M_2_PI      0.63661977236758134308
#endif
#ifndef M_2_SQRTPI
#define M_2_SQRTPI  1.12837916709551257390
#endif
#ifndef M_SQRT2
#define M_SQRT2     1.41421356237309504880
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2   0.70710678118654752440
#endif

#define OFFSET 1

#ifdef DEBUG
#define BUGREPORT fprintf(stderr,"%s %s compiled %s %s, Line %i\n",__FILE__,__func__,__DATE__,__TIME__,__LINE__);
#define BUGREPORT_PAR  if (ISROOT) fprintf(stderr,"ROOT: %s %s compiled %s %s, Line %i\n",__FILE__,__func__,__DATE__,__TIME__,__LINE__);
#else
#define BUGREPORT
#define BUGREPORT_PAR 
#endif


#ifdef COMM_VERBOSE
#define VERB(a) a
#else
#define VERB(a) 
#endif

#ifdef COMM_STATUS
#define COMM(a) a
#else
#define COMM(a) 
#endif

#ifndef MAX
#define MAX(x,y)                (((x)>(y))?(x):(y))
#endif

#ifndef MIN
#define MIN(x,y)                (((x)<(y))?(x):(y))
#endif

#ifndef HEAVYSIDE
#define HEAVYSIDE(x)                (((x)<=(0.))?(0.):(1.))
#endif

#ifndef SIGN
#define SIGN(x)                 (((x)>=0.)?(((x)>0.)?1.:0.):-1.)
#endif


#define PERIODIC 0
#define DIRDIR 1
#define DIRNEU 2
#define NEUNEU 3
#define NEUDIR 4
#define ZeroDIR 5
#define ZeroNEU 6


#define START_FROM_FILE FALSE
#define START_FROM_INI 10
#define RESTART TRUE
#define DEFAULTSTART 20


#define CARTESIAN   0
#define POLOIDAL    1
#define CYLINDRICAL 2
#define SPHERICAL   3
#define CART_TB_X   4
#define POL_TB_R    5
#define CYL_TB_R    6
#define SPHER_TB_R  7

#ifndef DEFSTRLEN 
#define DEFSTRLEN 1024
#endif 


typedef struct hdf_ds HDF_DS;
typedef struct para PARA;


/*! \struct para

Contains all parameters necessary for the code. Part of them are only used in some codes, so thisis a reservoir of predefined varaiables available. 
Should the need arise to add more parameters, please be aware that in this case you also have to change in par_helper.c
   the broadcast function for para, as it will otherwise break parallel programms!!!!!!
*/


struct para {
  double time, dt, end_time, out_time,
    xmax, xmin, dx, dkx,
    ymax, ymin, dy, dky,
    zmax, zmin, dz, dkz; // 16
  
  double kappan, kappat, alpha, beta, 
    delta, gamma, sigma, nu,
    r0, adrhos, betahat, muehat,
    mue_w, mue_n, mue_t; // 16+15 = 31

  double uvortex, radius, epsilon, v_pol,
    nprof, bprof, tprof, phiprof, 
    source, limiter, k0,
    dt_pol, nlf, nlp, 
    hm_nl, exb_ll, // 31+17 = 47
    q0,            // q in middle of domain
    shat0;         // shear in middle of domain  +2 = 49


  // Integral quantities

  double energy, dte, vorticity, Transport,
    heatflux, circulation, reynolds, prandel,
    nusselt; // 49 + 9 = 58

  double 
	qprof[1024]; //1024

  double 
	bdval[5][2]; // 10 

  double  
    amp[3], pos_x[3], pos_y[3], 
    width_x[3], width_y[3], width_random_x[3],
    width_random_y[3], amp_random[3]; // 8*3 = 24 

    // Boundaries

    int32 xbnd, ybnd, zbnd, otmult, r_offset, z_offset, r_spacing, z_spacing;
    int32 boundary[3], coordsys, verteilung;  // 13 int 

    /*
     * 
     *  END Transmit info 
     * 
     */

    double shat[1024],qss[1024],shift[1024]; 
    char desc[2*DEFSTRLEN],      // holds integrator specific descriptions. Not transmitted  24.03.00 
        codename[2*DEFSTRLEN],
        codedesc[2*DEFSTRLEN];   // description of the code in ini-comment format 
    char Prim_Phys[2*DEFSTRLEN]; //string for primary physics parameters 
    char Sec_Phys[2*DEFSTRLEN];  //string with secondary parameter descriptions
    char Prim_Geom[2*DEFSTRLEN]; // string with secondary parameter descriptions 

      /* 
       Primary units in [eV],[T],[1/m**3],[u], these are read in the .ini file 
       and potentially changed on the command line 
    */

    double 
        Ti,  /* in [eV] */
        Te,            /* in [eV] */
        B0,            /* in [T] */
        n0,            /* in [m**-3] */
        Mi,            /* in [u] atomic units */
        Z,               /* in electron charges [e] */
        p_n;          /* Neutral pressure in [pa] */ 

    double 
        part_source,  /* in [10 19 m-3 s-1] */
        temp_source;            /* in [eV 10 19 m-3 s-1] */   // 9 double

    /*
     * 
     *  Secondary, derived  units 
     * 
     */
    
    double rho_s,                            /* 1.02 *sqrt(Mi) sqrt(Te) / B                     in [m] */
        rho_i,                                       /* 1.02 *sqrt(Mi) sqrt(Ti) / B                      in [m] */ 
        rho_e,                                      
        omega_ci,                              /* 9.58 * 10**3 * Z/Mi  B                               in [rad/s] */ 
        omega_ce,                             /* 1.76 * 10**7 *  B                                          in [rad/s] */ 
        omega_pi,                             /* 1.32*  Z/sqrt(Mi) sqrt(n0)                        in [rad/s] */ 
        v_thi,                                       /* 9.79* 10**3* /sqrt(Mi)*sqrt(Ti)              in [m/s] */
        v_the,                                       /* 4.19* 10**5* sqrt(Te)                               in [m/s] */
        c_s,                                            /* 9.79*10**3 * sqrt(gamma  Z Te / Mi) in [m/s] */
        v_alfven,                               /* 2.18* 10**12 /sqrt(Mi) /sqrt(n_i)         in [m/s] */
        nu_ei,                                       /* in [1/s] */
        nu_in,                                 /* in [1/s] */
        nu_en,                                  /* in [1/s] */
        n_n,                                    /* Neutral density per cubic meter */
        unit_time,
        unit_length; // 16 double
    




    /* 
     * 
     * Primary dimensional geometry parameters 
     * 
     * 
     */
     
    
    double R0,            // major radius in [m]
        a,                // minor radius in [m]
        lpar,             // parallel length in [m]
        ln,               // density gradient scale length im [m]
        lTe,              // electron temperature gradient scale length in [m]
        lTi;              // ion temperature gradient scale length in [m] // 6 double 


    // Not used in paralle communication after here 
    /*
     * 
     *  BDS definitions for comparison added 30.11.2009 
     * 
     * 
     */
     
    double bds_epsilonhat,
        bds_betahat,
        bds_muehat,
        bds_nu,
        bds_wb,
        bds_A,
        bds_K,
        bds_C,
        bds_s;   

    double
        cflr, cflp;
    
}; // end para structure


/*! \struct HDF_DS

HDF data structure, contains basic information for field structure, such a dimensions, coordinates etc. 

   If you change this definition, don't forget to change the braoadcast in par_helper.c!! 
 */





struct hdf_ds 
{
    double  range[4][2];    /* 8 (marker)*/

    int32 offx,offy,offz,        /* 3 (marker) */ /*Number of ghostpoints in each direction*/
        nx, ny, nz, 			/*number of points Globally in each direction */
        lnx,lny,lnz,             /* 6 *//*number of points Locally in each direction */
        isfloat, rank,               /* 2 *//*rank denotes dimension not MPI rank (in )*/
        dims[4], 				/*Same as nx,ny,nz */
        start[4], 
        end[4], 
        edge[4],
        N[4];        /* Contains number of processes in each direction */ 
  /* 16   Sum= 11+16+4 = 31 */

    int32  zisperiodic,  /* marker */
        run_no, 
        read_data, 
        number,
        restart; /* 5 */

    char jobid[DEFSTRLEN], /* marker */
        cwd[DEFSTRLEN],
        desc[DEFSTRLEN],
        filename[DEFSTRLEN], /* 4*DEFSTRLEN */
        name[DEFSTRLEN],            /*  1* DEFSTLEN (marker)*/ 
        coordsys[DEFSTRLEN], 
        start_date[DEFSTRLEN], 
        revision[DEFSTRLEN], 
        integrator[DEFSTRLEN], 
        compile_date[DEFSTRLEN],        /* 5 * DEFSTRLEN */
        maschine[DEFSTRLEN], 
        write_date[DEFSTRLEN], 
        name_out[DEFSTRLEN], 
        name_in[DEFSTRLEN], 
        erhname[DEFSTRLEN];              /* 5 * DEFSTRLEN */


    // Not to be transmitted ......

    int ReadAttributes, 
        number_out, 
        this_process, // MPI rank
        num_procs,
        create;

 
    int POS[4];      // Own Position in Grid of Processors 
 
    int32 data_sets, file_attrs;
    int32 elements[4];
    char dim_label[4][DEFSTRLEN];
    int periodic[3]; /*! \var Pediodicity of axis*/
    int grid_id; /*Processor rank in cartesian communicator*/
    int xrow_id; /* Processor id in x row*/
    int zrow_id; /* Processor id in z row*/
    int grid_coords[3]; /*Coordinates of processor in cartesian grid of processors*/
    int neighbour[3][2];/*Contains ranks of neighbours for each processors: 
                          first number is axis, 0=Z,1=Y,2=X, 
                          second number determine the direction on axis, 
                          0 = decrease coordinate, 1= increase coordinate*/

    double *coordinate[4]; /* 4 */
    void  *datarw;
    int datarw_size;
    

    MPI_Comm cart_comm;
    MPI_Comm xrow_comm;
    MPI_Comm zrow_comm;
};

#endif       /*  !DEFINITIONS_H_ */
