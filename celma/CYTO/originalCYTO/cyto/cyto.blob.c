/**********************************************************************/

#ifdef FILEIO
#include <field.h>
#include <fileio.h>

int verbosity_level=0;
#endif

/* Defines */

#undef mathias_debug

#define PARALLEL_DAMPING

// Stefan
/*
#undef VSTATIC // Using Vstatic dereases necessary dt in first steps

#define ME (1./1836.2)
*/
#define VSTATIC
#define ME (1./1800)


#define DX(F)  0.5*hval[iz][ir]*( (F)[iz][ip][ir+1] - (F)[iz][ip][ir-1])
#define DY(F)  0.5*vval[iz][ir]*( (F)[iz][ip+1][ir] - (F)[iz][ip-1][ir])


/* These Flags are for debugging */
#undef DEBUG
#undef CNTFILES
#undef COMM_STATUS

#define HHSOLVER Laplace_Solve3D 

#define IS_ELECTRONVELOCITY 1
#define IS_IONVELOCITY 2
#define IS_DENSITY 3
#define IS_POTENTIAL 4
#define IS_TEMPERATURE 5
#define IS_VORTICITY 6
#define IS_CURRENT 7


#define FORALL for(iz=0;iz<nz;iz++) for(ip=0;ip<ny;ip++) for(ir=0;ir<nx;ir++)
#define FORALL_BD for(iz=-data.offz;iz<(nz+data.offz);iz++) for(ip=-data.offy;ip<(ny+data.offy);ip++) for(ir=-data.offx;ir<(nx+data.offx);ir++)

#define FORYX    for(ip=0;ip<ny;ip++)    for(ir=0;ir<nx;ir++)
#define FORYX_BD for(ip=-1;ip<ny+1;ip++) for(ir=-1;ir<nx+1;ir++) 

#define FORZY for(iz=0;iz<nz;iz++) for(ip=0;ip<ny;ip++) 
#define FORZY_BD for(iz=-data.offz;iz<nz+data.offz;iz++)  for(ip=-data.offy;ip<(ny+data.offy);ip++)

#define FORZX  for(iz=0;iz<nz;iz++)  for(ir=0;ir<nx;ir++)
#define FORZX_BD  for(iz=-data.offz;iz<(nz+data.offz);iz++)  for(ir=-data.offx;ir<(nx+data.offx);ir++)

#define FORX    for(ir=0;ir<nx;ir++) 
#define FORX_BD for(ir=-data.offx;ir<(nx+data.offx);ir++) 

#define FORY    for(ip=0;ip<ny;ip++) 
#define FORY_BD for(ip=-data.offy;ip<(ny+data.offy);ip++) 

#define FORZ    for(iz=0;iz<nz;iz++) 
#define FORZ_BD for(iz=-data.offz;iz<(nz+data.offz);iz++) 

#define FORALL_BD_XY for(iz= 0;iz<nz;iz++) for(ip=-1;ip<(ny+1);ip++) for(ir=-1;ir<(nx+1);ir++)


/***********************************************************************/

 
#include <mpi.h>
#include <utilities.h>
#include <par_helper.h>





#define blob_amplitude 16
#define bc_sheath      false
#define ic_velocity_amplitude 0


void debug_mean(double ***A, char *name, int time,
		int nz, int ny, int nx,
		int proc_id, int num_procs)
{
  double *total = (double*) malloc (nz* sizeof(double));  
  int iz, ip, ir;
  FORZ total[iz] = 0.0;
  FORALL{
    total[iz] += A[iz][ir][ip];
  }
  
  char filename[256];
  sprintf(filename, "Totalmass_time%d_%d_of_%d.dat", time, proc_id, num_procs);
  
  FILE *file; file=fopen(filename, "w");
  FORZ fprintf(file, "%lg \n", total[iz]);
  fclose(file);
  
  free(total);
}


/***********************************************************************/
/*

  Local Functions

*/
void Cyto_Laplace(double ***w,double ***f,double **hval,double **vval,
                  double *edr,int nx,int ny,int nz,PARA *p,int ispolar);

void Cyto_DParallelSOL(int identity,double ***result,double ***result2,
                       double ***val,double ***f,double ***n,double ***v,double ***dzv,
                       HDF_DS *data, PARA *para);

void Cyto_CalcParameters(PARA *p, HDF_DS *d) ;


void Map_Noise(HDF_DS *data, PARA *p);

void spectral_interpol(double *ret_val,double *xpos,double *ypos,int nop,double **phi,int nx, int ny, double *Sx, double *Sy);

double cyto_blob_distribution_z(int iz, int nz, HDF_DS *data);


/***********************************************************************/


int main(int argc,char **argv)
{ 
    // Stefan: {0} initializes *all* elements as 0, 0. or NULL
    HDF_DS 
        data = {0};
    PARA   
        para = {0};
    
    void (*ptrtoCalcParameters)(PARA *,HDF_DS *) = NULL;



    FILE   
        *output;

    double
        ***res      = NULL,
        ***tmp      = NULL,
        ***f_0      = NULL,
        ***vr       = NULL,
        ***vp       = NULL,
        ***tf       = NULL;
    
    double 
        **fbdrup    = NULL,
        **fbdrlow   = NULL,
        **fbdra     = NULL,
        **fbdrb     = NULL;
    
    double
        ***w        = NULL,
        ***w_0      = NULL,
        ***w_1      = NULL,
        ***w_2      = NULL,
        ***dtw_0    = NULL,
        ***dtw_1    = NULL,
        ***dtw_2    = NULL,
        **wbdra     = NULL,
        **wbdrb     = NULL,
        **ombdra    = NULL,
        **ombdrb    = NULL,
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
        ***target_density     = NULL,
        **nbdra     = NULL,
        **nbdrb     = NULL,
        **enbdra    = NULL,
        **enbdrb    = NULL,
        **nbdrlow   = NULL,
        **nbdrup    = NULL;
 
    double
        ***V_0      = NULL,
        ***V_1      = NULL,
        ***V_2      = NULL,
        ***dtV_0    = NULL,
        ***dtV_1    = NULL,
        ***dtV_2    = NULL,
        **vbdra     = NULL,
        **vbdrb     = NULL,
        **vbdrlow   = NULL,
        **vbdrup    = NULL;
    
    double
        ***U_0      = NULL,
        ***U_1      = NULL,
        ***U_2      = NULL,
        ***dtU_0    = NULL,
        ***dtU_1    = NULL,
        ***dtU_2    = NULL,
        **ubdra     = NULL,
        **ubdrb     = NULL,
        **ubdrlow   = NULL,
        **ubdrup    = NULL;
    
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
        ***dtdxf    = NULL,
        ***dtdyf    = NULL;
    
    double 
        ***dzU      = NULL,
        ***dzzU     = NULL,
        ***dzUN     = NULL,
        ***dzV      = NULL,
        ***dzzV     = NULL,
        ***dzVN     = NULL,
        ***dzF      = NULL,
        ***dzzF     = NULL,
        ***dzN      = NULL,
        ***dzzN     = NULL;
    
    double 
        ***J        = NULL,
        ***dzJ      = NULL,
        ***dzzJ     = NULL,
        ***nu_density = NULL;
    
    
    double
        *ts_n       = NULL,
        *ts_f       = NULL,
        *ts_e       = NULL;
    
    double 
        val,
        norm2d      = 0.0,
        nu_m,
        nu,
        delta,
        //  padvection  = 1.;
        padvection  = 0.; // GW version VN 15032010
    double
        **hval      = NULL,
        **vval      = NULL,
        **nlval     = NULL, 
        *rcor       = NULL,
        *edrcor     = NULL,
        *norm       = NULL,
        *dynbg      = NULL;
    
    double
        *par_elec_mom_flux = NULL,
        *par_ion_mom_flux  = NULL;
    
    double  
        totalf,
        Jtotal,
        Ipol,
        G_par, 
        J_par,
        mu, 
        r,
        cflr, 
        cflp;
    
    int 	
        nx,
        ny,
        nz,
        i,
        ir,
        ip,
        iz,
        counter,
        ot            = 1,
        iter          = 0,
        ispolar       = FALSE,
        WRITE         = TRUE,
        ErhFirstWrite = TRUE;
  
    char
        mom_name[256] = "Empty",
        filename[256] = "Empty",
        profilename[256] = "Empty";
    
    static char 
        rcsid[]="$Id: cyto.c,v 4.8 2003/10/16 14:08:09 vona Exp $";
    
    double  
        lamda    = 0., 
        gamma_n  = 0.,
        phase    = 0.;
    
    double 
        *par_rey_stress   = NULL,
        *flux_n           = NULL,
        *av_flux_n        = NULL;
  
    double
        density_source_integral,
        impulse           = 0.,
        dtimpulse         = 0.,
        nloss             = 0.,
        write_time,
        Lz                = 0.;
    
    int   
        rbdcndw, 
        rbdcndv,
        rbdcndu,
        rbdcndn,
        rbdcndf,
        zbndcndv,
        zbndcndw,
        zbndcndu,
        zbndcndn,
        zbndcndf;
    
    int
        m; 
    /* Variables for collecting communication to larger packages */

    double 
      reducing[25],
      result_vector[25];

#ifdef FILEIO
    char fname[256];

    /* field representations of para */
    field_dim D1k={1024}, D10={10}, D03={3};

    field
        S_time           = {"time"          , &para.time          , &FT_d, NULL, NULL, NULL},
        S_dt             = {"dt"            , &para.dt            , &FT_d, NULL, NULL, NULL},
        S_end_time       = {"end_time"      , &para.end_time      , &FT_d, NULL, NULL, NULL},
        S_out_time       = {"out_time"      , &para.out_time      , &FT_d, NULL, NULL, NULL},
        S_xmax           = {"xmax"          , &para.xmax          , &FT_d, NULL, NULL, NULL},
        S_xmin           = {"xmin"          , &para.xmin          , &FT_d, NULL, NULL, NULL},
        S_dx             = {"dx"            , &para.dx            , &FT_d, NULL, NULL, NULL},
        S_dkx            = {"dkx"           , &para.dkx           , &FT_d, NULL, NULL, NULL},
        S_ymax           = {"ymax"          , &para.ymax          , &FT_d, NULL, NULL, NULL},
        S_ymin           = {"ymin"          , &para.ymin          , &FT_d, NULL, NULL, NULL},
        S_dy             = {"dy"            , &para.dy            , &FT_d, NULL, NULL, NULL},
        S_dky            = {"dky"           , &para.dky           , &FT_d, NULL, NULL, NULL},
        S_zmax           = {"zmax"          , &para.zmax          , &FT_d, NULL, NULL, NULL},
        S_zmin           = {"zmin"          , &para.zmin          , &FT_d, NULL, NULL, NULL},
        S_dz             = {"dz"            , &para.dz            , &FT_d, NULL, NULL, NULL},
        S_dkz            = {"dkz"           , &para.dkz           , &FT_d, NULL, NULL, NULL}; /* 16 */

    field
        S_kappan         = {"kappan"        , &para.kappan        , &FT_d, NULL, NULL, NULL},
        S_kappat         = {"kappat"        , &para.kappat        , &FT_d, NULL, NULL, NULL},
        S_alpha          = {"alpha"         , &para.alpha         , &FT_d, NULL, NULL, NULL},
        S_beta           = {"beta"          , &para.beta          , &FT_d, NULL, NULL, NULL},
        S_delta          = {"delta"         , &para.delta         , &FT_d, NULL, NULL, NULL},
        S_gamma          = {"gamma"         , &para.gamma         , &FT_d, NULL, NULL, NULL},
        S_sigma          = {"sigma"         , &para.sigma         , &FT_d, NULL, NULL, NULL},
        S_nu             = {"nu"            , &para.nu            , &FT_d, NULL, NULL, NULL},
        S_r0             = {"r0"            , &para.r0            , &FT_d, NULL, NULL, NULL},
        S_adrhos         = {"adrhos"        , &para.adrhos        , &FT_d, NULL, NULL, NULL},
        S_betahat        = {"betahat"       , &para.betahat       , &FT_d, NULL, NULL, NULL},
        S_muehat         = {"muehat"        , &para.muehat        , &FT_d, NULL, NULL, NULL},
        S_mue_w          = {"mue_w"         , &para.mue_w         , &FT_d, NULL, NULL, NULL},
        S_mue_n          = {"mue_n"         , &para.mue_n         , &FT_d, NULL, NULL, NULL},
        S_mue_t          = {"mue_t"         , &para.mue_t         , &FT_d, NULL, NULL, NULL}; /* 15 */

    field
        S_uvortex        = {"uvortex"       , &para.uvortex       , &FT_d, NULL, NULL, NULL},
        S_radius         = {"radius"        , &para.radius        , &FT_d, NULL, NULL, NULL},
        S_epsilon        = {"epsilon"       , &para.epsilon       , &FT_d, NULL, NULL, NULL},
        S_v_pol          = {"v_pol"         , &para.v_pol         , &FT_d, NULL, NULL, NULL},
        S_nprof          = {"nprof"         , &para.nprof         , &FT_d, NULL, NULL, NULL},
        S_bprof          = {"bprof"         , &para.bprof         , &FT_d, NULL, NULL, NULL},
        S_tprof          = {"tprof"         , &para.tprof         , &FT_d, NULL, NULL, NULL},
        S_phiprof        = {"phiprof"       , &para.phiprof       , &FT_d, NULL, NULL, NULL},
        S_source         = {"source"        , &para.source        , &FT_d, NULL, NULL, NULL},
        S_limiter        = {"limiter"       , &para.limiter       , &FT_d, NULL, NULL, NULL},
        S_shear          = {"shear"         , &para.shear         , &FT_d, NULL, NULL, NULL},
        S_k0             = {"k0"            , &para.k0            , &FT_d, NULL, NULL, NULL},
        S_dt_pol         = {"dt_pol"        , &para.dt_pol        , &FT_d, NULL, NULL, NULL},
        S_nlf            = {"nlf"           , &para.nlf           , &FT_d, NULL, NULL, NULL},
        S_nlp            = {"nfp"           , &para.nlp           , &FT_d, NULL, NULL, NULL},
        S_hm_nl          = {"hm_nl"         , &para.hm_nl         , &FT_d, NULL, NULL, NULL},
        S_exb_ll         = {"exb_ll"        , &para.exb_ll        , &FT_d, NULL, NULL, NULL}; /* 17 */

    field
        S_energy         = {"energy"        , &para.energy        , &FT_d, NULL, NULL, NULL},
        S_dte            = {"dte"           , &para.dte           , &FT_d, NULL, NULL, NULL},
        S_vorticity      = {"vorticity"     , &para.vorticity     , &FT_d, NULL, NULL, NULL},
        S_Transport      = {"Transport"     , &para.Transport     , &FT_d, NULL, NULL, NULL},
        S_heatflux       = {"heatflux"      , &para.heatflux      , &FT_d, NULL, NULL, NULL},
        S_circulation    = {"circulation"   , &para.circulation   , &FT_d, NULL, NULL, NULL},
        S_reynolds       = {"reynolds"      , &para.reynolds      , &FT_d, NULL, NULL, NULL},
        S_prandel        = {"prandel"       , &para.prandel       , &FT_d, NULL, NULL, NULL},
/* ! */ S_nusselt        = {"Nusselt"       , &para.nusselt       , &FT_d, NULL, NULL, NULL}; /* 9 */

    field
        S_qprof          = {"qprof"         , &para.qprof         , &FT_d, &D1k, NULL, NULL},
        S_shat           = {"shat"          , &para.shat          , &FT_d, &D1k, NULL, NULL},
        S_qss            = {"qss"           , &para.qss           , &FT_d, &D1k, NULL, NULL},
        S_shift          = {"shift"         , &para.shift         , &FT_d, &D1k, NULL, NULL}; /* 1024*4 */

    field
        S_bdval          = {"bdval"         , &para.bdval         , &FT_d, &D10, NULL, NULL}, /* 10 */
        S_amp            = {"amp"           , &para.amp           , &FT_d, &D03, NULL, NULL},
        S_pos_x          = {"pos_x"         , &para.pos_x         , &FT_d, &D03, NULL, NULL},
        S_pos_y          = {"pos_y"         , &para.pos_y         , &FT_d, &D03, NULL, NULL},
        S_width_x        = {"width_x"       , &para.width_x       , &FT_d, &D03, NULL, NULL},
        S_width_y        = {"width_y"       , &para.width_y       , &FT_d, &D03, NULL, NULL},
        S_width_random_x = {"width_random_x", &para.width_random_x, &FT_d, &D03, NULL, NULL},
        S_width_random_y = {"width_random_y", &para.width_random_y, &FT_d, &D03, NULL, NULL},
        S_amp_random     = {"amp_random"    , &para.amp_random    , &FT_d, &D03, NULL, NULL}; /* 8*3 */

    field
        S_xbnd           = {"xbnd"          , &para.xbnd          , &FT_i, NULL, NULL, NULL},
        S_ybnd           = {"ybnd"          , &para.ybnd          , &FT_i, NULL, NULL, NULL},
        S_zbnd           = {"zbnd"          , &para.zbnd          , &FT_i, NULL, NULL, NULL},
/* ! */ S_otmult         = {"ot_mult"       , &para.otmult        , &FT_i, NULL, NULL, NULL},
        S_offset         = {"offset"        , &para.offset        , &FT_i, NULL, NULL, NULL},
/* ! */ S_boundary       = {"bndcnd"        , &para.boundary      , &FT_i, &D03, NULL, NULL},
        S_coordsys       = {"coordsys"      , &para.coordsys      , &FT_i, NULL, NULL, NULL},
        S_verteilung     = {"verteilung"    , &para.verteilung    , &FT_i, NULL, NULL, NULL}; /* 10*int */

    field
        S_Ti             = {"Ti"            , &para.Ti            , &FT_d, NULL, NULL, NULL},
        S_Te             = {"Te"            , &para.Te            , &FT_d, NULL, NULL, NULL},
        S_B0             = {"B0"            , &para.B0            , &FT_d, NULL, NULL, NULL},
        S_n0             = {"n0"            , &para.n0            , &FT_d, NULL, NULL, NULL},
        S_Mi             = {"Mi"            , &para.Mi            , &FT_d, NULL, NULL, NULL},
        S_Z              = {"Z"             , &para.Z             , &FT_d, NULL, NULL, NULL},
        S_p_n            = {"p_n"           , &para.p_n           , &FT_d, NULL, NULL, NULL}; /* 7 */

#define PARA_GROUP_1 \
        &S_time, &S_dt, &S_end_time, &S_out_time, \
        &S_xmax, &S_xmin, &S_dx, &S_dkx, \
        &S_ymax, &S_ymin, &S_dy, &S_dky, \
        &S_zmax, &S_zmin, &S_dz, &S_dkz

#define PARA_GROUP_2 \
        &S_kappan, &S_kappat, &S_alpha, &S_beta, \
        &S_delta, &S_gamma, &S_sigma, &S_nu, \
        &S_r0, &S_adrhos, &S_betahat, &S_muehat, \
        &S_mue_w, &S_mue_n, &S_mue_t

#define PARA_GROUP_3 \
        &S_uvortex, &S_radius, &S_epsilon, &S_v_pol, \
        &S_nprof, &S_bprof, &S_tprof, &S_phiprof, \
        &S_source, &S_limiter, &S_shear, &S_k0, \
        &S_dt_pol, &S_nlf, &S_nlp, &S_hm_nl, &S_exb_ll

#define PARA_GROUP_4 \
        &S_energy, &S_dte, &S_vorticity, &S_Transport, \
        &S_heatflux, &S_circulation, &S_reynolds, &S_prandel, &S_nusselt

#define PARA_GROUP_5 \
        &S_qprof, &S_shat, &S_qss, &S_shift

#define PARA_GROUP_6 \
        &S_bdval, &S_amp, &S_pos_x, &S_pos_y, \
        &S_width_x, &S_width_y, &S_width_random_x, &S_width_random_y, \
        &S_amp_random

#define PARA_GROUP_7 \
        &S_xbnd, &S_ybnd, &S_zbnd, &S_otmult, \
        &S_offset, &S_boundary, &S_coordsys, &S_verteilung

#define PARA_GROUP_8 \
        &S_Ti, &S_Te, &S_B0, &S_n0, \
        &S_Mi, &S_Z, &S_p_n

    field 
        *SP_para1[] = {PARA_GROUP_1, NULL},
        *SP_para2[] = {PARA_GROUP_2, NULL},
        *SP_para3[] = {PARA_GROUP_3, NULL},
        *SP_para4[] = {PARA_GROUP_4, NULL},
        *SP_para5[] = {PARA_GROUP_5, NULL},
        *SP_para6[] = {PARA_GROUP_6, NULL},
        *SP_para7[] = {PARA_GROUP_7, NULL},
        *SP_para8[] = {PARA_GROUP_8, NULL};

    field *SP_para[]  = {
        PARA_GROUP_1, PARA_GROUP_2, PARA_GROUP_3, PARA_GROUP_4, 
        PARA_GROUP_5, PARA_GROUP_6, PARA_GROUP_7, PARA_GROUP_8, NULL
    };

    para.Z = 1.;
    para.SP = SP_para;

    /* field representations of actual fields */
    float ***f_hel;
    field_dim D, D_hel;

    /* helper field is float and described by field_dim object D_hel */
    field F_hel = {NULL, &f_hel, &FT_f, &D_hel};

    field
        F_exp_n_0 = {"n"          , &exp_n_0, &FT_d, &D, &F_hel, "Density"              },
        F_omega   = {"Vorticity"  , &omega  , &FT_d, &D, &F_hel, "Vorticity"            },
        F_n_0_ini = {"Density"    , &n_0    , &FT_d, &D, &F_hel, "Log density"          },
        F_n_0     = {"N"          , &n_0    , &FT_d, &D, &F_hel, "Log density"          },
        F_f_0     = {"Potential"  , &f_0    , &FT_d, &D, &F_hel, "Potential"            },
        F_U_0     = {"UIon"       , &U_0    , &FT_d, &D, &F_hel, "Ion velocity"         },
        F_V_0_ini = {"Velocity"   , &V_0    , &FT_d, &D, &F_hel, "Electron velocity"    },
        F_V_0     = {"VElectron"  , &V_0    , &FT_d, &D, &F_hel, "Electron velocity"    },
        F_tf      = {"Source"     , &tf     , &FT_d, &D, &F_hel, "Source"               },
        F_dzU     = {"dzU"        , &dzU    , &FT_d, &D, &F_hel, "dzU"                  },
        F_dzV     = {"dzV"        , &dzV    , &FT_d, &D, &F_hel, "dzV"                  },
        F_dzN_cnt = {"dzN"        , &dzN    , &FT_d, &D, &F_hel, "dzN"                  },
        F_dzN     = {"dzn"        , &dzN    , &FT_d, &D, &F_hel, "dzN"                  },
        F_dzF_cnt = {"dzF"        , &dzF    , &FT_d, &D, &F_hel, "dzF"                  },
        F_dzF     = {"dzf"        , &dzF    , &FT_d, &D, &F_hel, "dzF"                  },
        F_w_0     = {"H"          , &w_0    , &FT_d, &D, &F_hel, "Generalized Vorticity"},
        F_dtw_0   = {"DiffHW"     , &dtw_0  , &FT_d, &D, &F_hel, "DiffHW"               },
        F_vst_cnt = {"VEstatic"   , &vstatic, &FT_d, &D, &F_hel, "Static EVelocity"     },
        F_vstatic = {"Vstatic"    , &vstatic, &FT_d, &D, &F_hel, "Static EVelocity"     },
        F_J       = {"Current"    , &J      , &FT_d, &D, &F_hel, "Current"              },
        F_dyn     = {"ParForce"   , &dyn    , &FT_d, &D, &F_hel, "Parallel force"       },
        F_dxn     = {"ParMomentum", &dxn    , &FT_d, &D, &F_hel, "Parallel momentum"    },
        F_vr      = {"vr"         , &vr     , &FT_d, &D, &F_hel, "vr"                   },
        F_vp      = {"vp"         , &vp     , &FT_d, &D, &F_hel, "vp"                   };

    const field *FP_ini[] = {
        &F_omega, &F_n_0_ini, &F_f_0, &F_U_0, &F_V_0_ini, &F_tf, &F_dzU, &F_dzN, &F_dzF, &F_w_0, NULL
    };

    const field *FP_cnt[] = {
        &F_exp_n_0, &F_U_0, &F_dzV, &F_w_0, &F_f_0, &F_n_0, &F_V_0, &F_omega, &F_dtw_0, &F_vst_cnt,
        &F_dzN_cnt, &F_dzF_cnt, &F_J, NULL
    };

    const field *FP_write[] = {
        &F_exp_n_0, &F_U_0, &F_dzN, &F_w_0, &F_f_0, &F_n_0, &F_V_0, &F_omega, &F_dzF, &F_vstatic,
        &F_J, &F_dyn, &F_dxn, NULL
    };

    const field *FP_VPR[] = {
        &F_vp, &F_vr, &F_f_0, NULL
    };

#endif


    
   /************************************************/
    /*                                              */
    /* start of body/read data and command line     */
    /*                                              */
    /************************************************/  

    BUGREPORT;
    /* Initialize MPI */

    MPI_Init(&argc,&argv); 
    MPI_Comm_rank(MPI_COMM_WORLD,&data.this_process);
    COMM(fprintf(stderr,"Proc %d initializing...\n",data.this_process););
    MPI_Comm_size(MPI_COMM_WORLD,&data.num_procs);
    COMM(fprintf(stderr,"Procs %d\n",data.num_procs););


    BUGREPORT;

    /* MPI_Datatypes */

    hdfdata_type = create_mpi_hdfds(&data,ROOT, MPI_COMM_WORLD,data.this_process);
    para_type = create_mpi_para(&para, ROOT, MPI_COMM_WORLD,data.this_process);


    sprintf(para.codedesc,
     "# The CYTO code solves the global, electrostatic fluid equations for a plasma in a cylinder:\n"\
     "#\n"                                                              \
     "#  N = ln n \n   and H = omega + Nabla N nabla phi"               \
     "#  vorticity = omega = laplace phi\n"                             \
     "#\n"                                                              \
     "#  Dt N = - dz V -  V dz N  \n"                                   \
     "#\n"                                                              \
     "#  dt H = 1/Z*(dz (U-V) + (U-V) dz N  + nabla N ( mu nabla w - nu nabla f)" \
     "    +{phi, omega}  + nabla n { phi, nabla phi} + nabla dt n . nabla phi   " \
     "   + mu nabla^2 H -nu w \n"                                       \
     "#\n"                                                              \
     "#  V = U - tau  D_par (phi - N)\n"                                \
     "#\n"                                                              \
     "#  D_t U = D_par N - 2 U dz U \n"                                 \
	 "#\n");
 
 
    sprintf(para.desc, 
         "gamma: background current;"                                   \
         "sigma: ion-neutral collisions;"                               \
         "nu: electron-ion collisions;"                                 \
         "beta: -;"                                                     \
         "shear: background shear flow;"                                \
         "alpha: -;"                                                    \
         "delta: electron-neutral collisions;"                          \
         "betahat: -;"                                                  \
         "adrhos: -;"                                                   \
         "mue_w: viscosity;"                                            \
         "mue_n: diffusion;"                                            \
         "mue_t: -;"                                                    \
         "kappan: width of gaussian density source (length);"           \
         "kappat: -;"                                                   \
         "r0: -;"                                                       \
         "muehat: -;"                                                   \
         "bprof: -;"                                                    \
         "tprof: -;"                                                    \
         "nprof: density source drive rate;"                            \
         "phiprof: -;"                                                  \
         "hm_nl: -;"                                                    \
         "exb_ll: -;"                                                   \
         "limiter: parallel damping (artificial);"                      \
         "source: background source to keep density from falling to zero (GW = 1);"                                                   \
         "k0: -;"                                                       \
         "dt_pol: -;"                                                   \
         "uvortex: -;"                                                  \
         "radius: -;");



 sprintf(para.Prim_Phys,     /* holds string for primary physics parameters */
         "Ti: Central Ion Temperature [eV];"     \
         "Te: Central Electron Temperature [ev];"\
         "B0: Magnetic field strength [T];"      \
         "n0: central density [1e19 m**-3];"     \
         "Mi: Ion mass in [u];"                  \
         "Z: Charge of main Ion species [e];"    \
         "p_n: Neutral pressure in [pa];"        \
         );     



 sprintf(para.Sec_Phys, /* holds string with secondary parameter describtions  */
         "rho_s: 1.02 *sqrt(Mi)*sqrt(Te) / B    in [m];"       \
         "rho_i: 1.02 *sqrt(Mi)*sqrt(Ti) / B      in [m];"     \
         "omega_ci: 9.58 * 10**3 * Z/Mi  B       in [rad/s];"  \
         "omega_ce: 1.76 * 10**7 *  B     in [rad/s] ;"        \
         "omega_pi: 1.32*  Z/sqrt(Mi)*sqrt(n0)  in [rad/s] ;"  \
         "v_thi:  9.79* 10**3* /sqrt(Mi)*sqrt(Ti)   in [m/s];" \
         "v_the:  4.19* 10**5* sqrt(Te)  in [m/s];"            \
         "c_s:  9.79*10**3 * sqrt(gamma  Z Te / Mi) in [m/s];" \
         "v_alfven:  2.18* 10**12 /sqrt(Mi) /sqrt(n_i)   in [m/s];");

 sprintf(para.Prim_Geom, /* holds string with primary geometry parameter describtions  */
         "R0: -;"                                       \
         "a: radius of plasma column      in [m];"      \
         "lpar: parallel length in [m];"                \
         "ln: density gradient length scale  in [m];"   \
         "lTe: -;"                                      \
         "lTi: -;" );



 strcpy(para.codename,"cyto");
 sprintf(para.codename,"%s_blob_%d", para.codename, blob_amplitude);

 if( bc_sheath ) sprintf(para.codename,"%s_sheath", para.codename);
 else sprintf(para.codename,"%s_no_sheath", para.codename);

 sprintf(para.codename,"%s_ic_velocity_%d", para.codename, ic_velocity_amplitude);


#ifndef FILEIO 
 /* Initialize data structures (necessary!) */
 FUtils_IniStructure(&data,&para,argv);
#endif

 ptrtoCalcParameters = Cyto_CalcParameters;
 
 
 /* 
    Default values for parameters, these will be written into the sample.ini file if the code is called with the -H option 
    Do not change these values, as they will be used to define a standard run!

*/

 data.nx                   = 40; /*Number of grid points in x*/
 data.ny                   = 64; /*Number of grid points in y*/
 data.nz                   = 70; /*Number of grid points in z*/
 para.coordsys             = CYLINDRICAL ;


 para.dt                   = 1.e-6; /*Delta T*/
 para.end_time             = 1.000000; /*Simulation stops at this time*/
 para.out_time             = 5e-3; /*Time between small outputs*/
 para.otmult               = 1; /*Number of small outputs  before output of fields*/

 para.Ti                   = 0.030000; /*Central Ion Temperature [eV]*/
 para.Te                   = 2.000000; /*Central Electron Temperature [ev]*/
 para.B0                   = 0.070000; /*Magnetic field strength [T]*/
 para.n0                   = .005;/* central density [10^19 1/m**3]*/
 para.Mi                   = 40; /*Ion mass in [u]*/
 para.Z                    = 1.000000; /*Charge of main Ion species [e]*/
 para.p_n                  = 0.300000; /*Neutral pressure in [pa]*/

 para.a                    = 0.050000; /*radius of plasma column      in [m]*/
 para.lpar                 = 1.500000; /*parallel length in [m]*/
 para.ln                   = 0.0250000; /*density gradient length scale  in [m]*/


 para.limiter              = 10.000000; /*parallel damping (artificial)*/
 para.kappan               = 0.000000; /*width of gaussian density source*/
 para.nprof                = 0.250000; /*density source*/


 /*





 */

 data.rank   = 3;
 data.offz   = 1;
 data.offx   = 1;
 data.offy   = 1; 
 data.number = 0;


 MPI_Barrier(MPI_COMM_WORLD);

 /*****************************************************************************/
 /* First we only read information and not the actual fields....              */
 /*****************************************************************************/

  
 if(ISROOT)
 {
     data.read_data      = FALSE;
     data.ReadAttributes = TRUE;
     data.restart        = FUtils_ReadArguments(argc,argv,data.filename,&data.number,&data, &para,ptrtoCalcParameters);
 }
 MPI_Barrier(MPI_COMM_WORLD);
 
 /* Broadcast DATA  out to the other processors */
 MPI_Bcast(&data,1,hdfdata_type,ROOT,MPI_COMM_WORLD);
 
   
 sprintf(data.erhname,"%s.erh",data.name_out);
 sprintf(profilename,"%s.prof.dat",data.name_out);
 sprintf(mom_name,"%s.momentum.dat",data.name_out);
    
 /***********************************/
 /* Geometry of processes           */
 /***********************************/
 data.zisperiodic = TRUE; /* This is needed to set up periodic processor grid  */
 PH_3D_Geometry(&data);
 
 data.lnx = nx = data.elements[2]; 
 data.lny = ny = data.elements[1];  
 data.lnz = nz = data.elements[0];    
 if(ISROOT) {
     fprintf(stderr,"%ld %ld %ld \n",data.dims[2],data.dims[1],data.dims[0]);
     fprintf(stderr,"%ld %ld %ld \n",data.nx,data.ny,data.nz);
     fprintf(stderr,"%ld %ld %ld \n",data.offx,data.offy,data.offz);
 }

 
 BUGREPORT;
 /************************************************/
 /*                                              */
 /* Allocation of variables and fields                */
 /*                                              */
 /************************************************/


 
// Stefan: shifts were 3,2
 Util_3DAllocFields(&w_0,&w_1,&w_2,0,&dtw_0,&dtw_1,&dtw_2,0,
                    &wbdra,&wbdrb,&wbdrlow,&wbdrup,
                    nx,ny,nz,data.offx,data.offy,data.offz);

 Util_3DAllocFields(&n_0,&n_1,&n_2,0,&dtn_0,&dtn_1,&dtn_2,0,
                    &nbdra,&nbdrb,&nbdrlow,&nbdrup,
                    nx,ny,nz,data.offx,data.offy,data.offz);  
 
 Util_3DAllocFields(&V_0,&V_1,&V_2,0,&dtV_0,&dtV_1,&dtV_2,0,
                    &vbdra,&vbdrb,&vbdrlow,&vbdrup,
                    nx,ny,nz,data.offx,data.offy,data.offz);  
 
 Util_3DAllocFields(&U_0,&U_1,&U_2,0,&dtU_0,&dtU_1,&dtU_2,0,
                    &ubdra,&ubdrb,&ubdrlow,&ubdrup,
                    nx,ny,nz,data.offx,data.offy,data.offz);  

// Stefan: shifts were 1,1
 Util_3DAllocFields(&f_0,&vr,&vp,0,&res,&n,&w,0,
                    &fbdra,&fbdrb,&fbdrlow,&fbdrup,
                    nx,ny,nz,data.offx,data.offy,data.offz);

 BUGREPORT;
/* Cubes */
// Stefan, shifts were 1,1,1,1,2,2,5,2,3,3,4,2,2,2,...
 dzF            = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);   
 dzzF           = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);
 dzN            = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);   
 dzzN           = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);
 dzV            = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);   
 dzzV           = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);
 dzVN           = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);
 dzU            = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);   
 dzzU           = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);   
 dzUN           = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);
 dtdxf          = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0); 
 dtdyf          = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);
 dxn            = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);   
 dyn            = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);     
 dxf            = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);     
 dyf            = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);
 vstatic        = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);  
 J              = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);  
 dzJ            = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);  
 dzzJ           = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);
 tf             = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0); 
 omega          = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);     
 exp_n_0        = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);     
 target_density = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);
 nu_density     = Util_DCube(nx,data.offx,ny,data.offy,nz,data.offz,0);     

 BUGREPORT;
 /* 2d fields */
 enbdra         = Util_DMatrix(nz,data.offz,ny,data.offy); 
 enbdrb         = Util_DMatrix(nz,data.offz,ny,data.offy);
 ombdra         = Util_DMatrix(nz,data.offz,ny,data.offy);
 ombdrb         = Util_DMatrix(nz,data.offz,ny,data.offy);
 nlval          = Util_DMatrix(nz,0,nx,0);/* Tight field for Arakawa optimisation */
 hval           = Util_DMatrix(nz,data.offz,nx,data.offx);
 vval           = Util_DMatrix(nz,data.offz,nx,data.offx);
 
   
 /* Vectors */
 rcor              = Util_DVector(nx,data.offx);
 edrcor            = Util_DVector(nx,data.offx);
 norm              = Util_DVector(nx,data.offx);
 par_rey_stress    = Util_DVector(nx,data.offx);
 flux_n            = Util_DVector(nx,data.offx);
 av_flux_n         = Util_DVector(nx,data.offx);
 dynbg             = Util_DVector(nx,data.offx);
 u_par             = Util_DVector(nx,data.offx);
 pol_pot           = Util_DVector(ny,data.offy);
 pol_dens          = Util_DVector(ny,data.offy);
 par_elec_mom_flux = Util_DVector(nx,data.offx);
 par_ion_mom_flux  = Util_DVector(nx,data.offx);
 ts_n              = Util_DVector((int)(para.out_time/para.dt)+100,1);
 ts_f              = Util_DVector((int)(para.out_time/para.dt)+100,1);
 ts_e              = Util_DVector((int)(para.out_time/para.dt)+100,1);
 BUGREPORT;


 /*****************************************************************************/  
 /*****************************************************************************/  

 /*   Setup Geometry    */

 ispolar = TRUE; /* Always True, maybe change later to allow rectangular geometry again */
 BUGREPORT;
 sprintf(data.coordsys,"polar, equidistant");
 sprintf(data.dim_label[0],"z");
 sprintf(data.dim_label[2],"r");
 sprintf(data.dim_label[1],"phi");

 if(ISROOT) {
   fprintf(stderr,"The geometry is ");
   if(ispolar)fprintf(stderr,"polar.\n");
   else fprintf(stderr,"rectangular.\n");

   fprintf(stderr,"%d cordsys %s\n",__LINE__,data.coordsys);
   fprintf(stderr,"%s %s %s \n",data.dim_label[2],data.dim_label[1],data.dim_label[0]);
   fprintf(stderr,"%ld %ld %ld \n",data.dims[2],data.dims[1],data.dims[0]);
   fprintf(stderr,"%ld %ld %ld \n",data.nx,data.ny,data.nz);
 }


 /*****************************************************************************/  
 /* RBDCND:      X = A           X = B */
 BUGREPORT;
 /*****************************************************************************/  
 if(para.xmin == 0. && ispolar) 
 {
     rbdcndf   = ZeroDIR; 
     rbdcndw   = ZeroDIR; 
     rbdcndn   = ZeroDIR;
     rbdcndu   = ZeroNEU;
     rbdcndv   = ZeroNEU;  
 } 
 else
 {
     rbdcndf   = DIRDIR; 
     rbdcndw   = DIRDIR; 
     rbdcndn   = DIRNEU;
     rbdcndu   = NEUNEU;    
     rbdcndv   = NEUNEU;
 }

 
/* Values on the boundaries */
   BUGREPORT;
   
 FORZY_BD nbdra[iz][ip]  = 0.0;
 FORZY_BD nbdrb[iz][ip]  = 0.0;
 FORZY_BD enbdra[iz][ip] = exp(nbdra[iz][ip]);
 FORZY_BD enbdrb[iz][ip] = exp(nbdrb[iz][ip]);
 FORZY_BD fbdra[iz][ip]  = 0.;
 FORZY_BD fbdrb[iz][ip]  = 0.;
 FORZY_BD wbdra[iz][ip]  = 0.;
 FORZY_BD wbdrb[iz][ip]  = 0.;
 FORZY_BD ombdra[iz][ip] = 0.;
 FORZY_BD ombdrb[iz][ip] = 0.;
 FORZY_BD ubdra[iz][ip]  = 0.;
 FORZY_BD ubdrb[iz][ip]  = 0.;
 FORZY_BD vbdra[iz][ip]  = 0.;
 FORZY_BD vbdrb[iz][ip]  = 0.;

/*****************************************************************************/
/* ZBDCND:      Z = down           Z = up*/

 if( false ){ // mathias: can't use for the blob study

   zbndcndf = DIRNEU; 
   zbndcndw = DIRNEU; 
   zbndcndn = DIRNEU; 
   zbndcndu = NEUNEU;
   zbndcndv = NEUNEU;

 } else { // mathias : necessary for the blob study
   
   zbndcndf = DIRDIR; 
   zbndcndw = DIRDIR; 
   zbndcndn = DIRDIR; 
   zbndcndu = DIRDIR;
   zbndcndv = DIRDIR;
   
 }  

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

 BUGREPORT;

/************************ INITIAL CONDITION **********************************************/
 
 /****************************************************/
 /*                                                  */
 /*   READ THE DATA FIELDS                           */
 /*                                                  */
 /****************************************************/
 
 BUGREPORT;
 
 data.ReadAttributes = FALSE;
 sprintf(data.integrator,"%s",argv[0]);
 sprintf(data.revision,"Empty");  
 sprintf(data.compile_date,"%s %s",__DATE__,__TIME__);


 if(ISROOT) fprintf(stderr,"%d cordsys %s\n",__LINE__,data.coordsys);

 para.otmult = MAX(para.otmult,1);
 MPI_Bcast(&para,1,para_type,ROOT,MPI_COMM_WORLD);
 MPI_Barrier(MPI_COMM_WORLD);

 BUGREPORT;
 if(ISROOT) fprintf(stderr,"%d cordsys %s\n",__LINE__,data.coordsys);
 if(data.restart == RESTART)
   {
     printf("RESTART...\n");

     /******************RESTART FROM FILE MADE BY CYTO ****************************************/
     /* Read all the fields from a file, but allow overide of parameters from commandline, 
        e.g. do not read parameters here ! */
     data.read_data = TRUE;
     data.ReadAttributes = FALSE;
     BUGREPORT;     
     counter = 0;
     
#ifdef FILEIO
     fprintf(stderr,"Please implement FileIO for RESTART or configure with --disable-fileio\n");
     exit(0);
#else
#ifdef HDF4 
     /* 
        Read fields which might appear under different names
     */

     if(PH_Read3DFieldbyName(n_0,"N",&data,&para) > 0 && (ISROOT) )
     {
         fprintf(stderr,"Restart success, read field N\n");
         
     }
     else if(PH_Read3DFieldbyName(n_0,"n",&data,&para)> 0 && (ISROOT) ) 
     {    
         FORALL_BD  n_0[iz][ip][ir]  = log(n_0[nz][ip][ir]);
         fprintf(stderr,"Restart success, read field n \n");
     }
     else if(PH_Read3DFieldbyName(n_0,"Density",&data,&para)> 0 && (ISROOT) )
     {    
         FORALL_BD  n_0[iz][ip][ir]  = log(n_0[nz][ip][ir]);
         fprintf(stderr,"Restart success, read field Density \n");
     }
     else if (ISROOT) 
     {
         fprintf(stderr,"Restart failure, could not read field N\n");
         counter++;
     }
     PH_3D_Write( n_0 ,"N","RESTART",1,&data,&para,TRUE);     



     if(PH_Read3DFieldbyName(w_0,"H",&data,&para) > 0 && (ISROOT) ) 
     {
         fprintf(stderr,"Restart success, read field H\n");
     }
     else if (ISROOT)
     {
         fprintf(stderr,"Restart failure, could not read field H\n");
         counter++;
     }
     PH_3D_Write( w_0 ,"H","RESTART",2,&data,&para,TRUE);  
   
  

     if(PH_Read3DFieldbyName(U_0,"UIon",&data,&para)> 0 && (ISROOT) ) 
     {
         fprintf(stderr,"Restart success, read field UIon\n");
     }
     else if (ISROOT)
     {
         fprintf(stderr,"Restart failure, could not read field UIon\n");
         counter++;    
     }
     PH_3D_Write( U_0 ,"UIon","RESTART",3,&data,&para,TRUE);     


     
     if(PH_Read3DFieldbyName(V_0,"VElectron",&data,&para)> 0 && (ISROOT) ) 
     {
         fprintf(stderr,"Restart success, read field VElectron\n");
     }
     else if(PH_Read3DFieldbyName(V_0,"Velocity",&data,&para)> 0 && (ISROOT) ) 
     {
         fprintf(stderr,"Restart success, read field Velocity\n");
     }
     else if (ISROOT)
     {
         fprintf(stderr,"Restart failure, could not read field VElectron\n");
         counter++;
     }
     PH_3D_Write( V_0 ,"V_0","RESTART",4,&data,&para,TRUE);     



     if(PH_Read3DFieldbyName(f_0,"Potential",&data,&para)> 0 && (ISROOT) ) 
     {
         fprintf(stderr,"Restart success, read field Potential\n");
     }
     else if (ISROOT)
     {
         fprintf(stderr,"Restart failure, could not read field Potential\n");
         counter++;
     }
     PH_3D_Write( f_0 ,"Potential","RESTART",5,&data,&para,TRUE);     
#elif HDF5
     printf("HDF5 and Restart is not avaible yet \n");
#endif /* !HDF4/5 */
#endif /*FILEIO*/

     if (counter > 0)
     {
         fprintf(stderr,"Restart failure, could not read %d fields. Cyto cannot restart from this file.\n",counter);
         MPI_Finalize();
         exit(0);
         
     }

   }
   else if (data.restart == START_FROM_INI ) 
   {
	printf("START_FROM_INI...\n");

     /******************START FROM INI FILE  ****************************************/

	if(data.this_process > data.nz){
	  printf("Error: the initialisation of n is not proper! \n");
	  printf("Number of CPU > data.nz\n");
	  MPI_Finalize();
	  exit(0);
	}

       FORALL_BD  n_0[iz][ip][ir] 	= 1.  /*+ exp(-rcor[ir]*rcor[ir]*(para.kappan*para.kappan))
                                            +(1.e-5*(double)rand()/(double)RAND_MAX)*exp(-( (double)(ir-nx/2)/(double)nx)*((ir-nx/2)/(double)nx)*100.)
                                            *(1.0 - 0.5*(double)iz/(double)nz)*/;
       if(false)
       for(m=2;m<ny/8;m++)
       {
           phase =  (double)rand()/(double)RAND_MAX ;
           
           FORALL_BD   n_0[iz][ip][ir]+= 
	     0.001 / (double)(m*m);
	     //  *cos(2.*M_PI*m*((double)ip/(double)ny+phase));
             //  *exp(-( (double)(ir-nx/2)/(double)nx)*((ir-nx/2)/(double)nx)*500.   );

	   
       }
       
       int 
	 x_c = 1, y_c = ny/2 , z_c = ((double)data.dims[0]-1.0) / 2.0;

       double loc_x, loc_y, loc_z;

       if(false)
       for(iz=-data.offz;iz<(nz+data.offz);iz++) {
	 printf("proc = %d \t", data.this_process );
	 printf("z_c = %d \t", z_c );
	 printf("iz  = %d \n", (nz* data.grid_coords[0] + iz));
       }
       MPI_Barrier(MPI_COMM_WORLD);

       int w = 10;
       double 
	 w_x  = 1./para.dx* (double)w, w_z = 1./para.dz* (double)w;

       FORALL_BD{
	 

	 n_0[iz][ip][ir] = 1.001 + (double)blob_amplitude*
	   exp(- (x_c - ir)*(x_c - ir)/ w_x )
	   //	 *exp(- (y_c - ip)*(y_c - ip)/ w_y )
	   *exp(- (z_c - (nz* data.grid_coords[0] + iz))*(z_c - (nz* data.grid_coords[0] + iz))/ w_z );
       }
       //       MPI_Barrier(MPI_COMM_WORLD);
       //       PH_3D_Write(n_0,"Density","test",0,&data,&para,TRUE);
       /* FORZ */
       /* 	 printf("(%d) %lg \t %lg \n", iz, n_0[iz][0][0], n_0[iz][10][10]); */
       /* MPI_Barrier(MPI_COMM_WORLD); */
       /* exit(0); */

       //cyto_blob_distribution_z(iz, nz, &data)
       FORALL_BD  n_0[iz][ip][ir]  = log(abs(n_0[iz][ip][ir]));
       //       FORALL_BD  n_0[iz][ip][ir]  = - n_0[iz][ip][ir];
       // PH_3D_Write(n_0,"Density","test2",0,&data,&para,TRUE);

       /* MPI_Barrier(MPI_COMM_WORLD); */
       /* MPI_Finalize(); */
       
#ifdef FILEIO
       if( ISROOT ){

       int loc_rank    = 3;
       int loc_size[]  = {data.nz, data.ny, data.nx};
       int loc_offs[]  = {data.offz,data.offy,data.offx};
       int loc_zeros[] = {0,0,0};
       
       /* Construct and setup D */
       FieldDim_Construct(&D, loc_rank, loc_size, loc_offs);

       FieldDim_CopyCoordsys(&D, data.coordsys);
       for (i=0; i<D.rank; ++i) {
           FieldDim_CopyName(&D, data.dim_label[i], i);
           D.xm[i] = data.range[i][0];
           D.xM[i] = data.range[i][1];
           memcpy(D.x[i], data.coordinate[i], D.n[i]*sizeof(double));
       }

       /* Construct and setup D_hel (same as D, but without ghost points) */
       FieldDim_Construct(&D_hel, loc_rank, loc_size, loc_zeros);
       FieldDim_CopyDimsCore(&D_hel, &D);

       /* Construct F_hel */
       Field_Construct(&F_hel);
       } // end if ISROOT

       /* Read initial values from file (NOT OK for -I startup) */
#ifdef HDF4
       sprintf(fname,"SET009_h4.000");
#endif
#ifdef HDF5
       sprintf(fname,"SET009_h5.000");
#endif
       /* Use a brand new field for reading so we don't interfere with anything else */
       /* field *F_new = Field_New(); */
       /* Field_CopyName(F_new, "n"); */
       /* FileIO_GetInfo(fname, F_new); */
       /* Field_Construct(F_new); */
       /* FileIO_Read(fname, F_new); */

       /* /\* Copy field read into F_new over to exp_n_0 *\/ */
       /* float *tmp = Field_First(F_new); */
       /* FORALL exp_n_0[iz][ip][ir] = *tmp++; */

       /* Field_Delete(F_new); */

       /* /\* Set initial value for n_0 *\/ */
       /* FORALL n_0[iz][ip][ir] = log(fabs(exp_n_0[0][ip][ir])); */

#endif /* !FILEIO */
       BUGREPORT;

       /* w ,f  get initialized with zero values */   
       FORALL_BD omega[iz][ip][ir]  = 0.;
       FORALL_BD w_0[iz][ip][ir]    = 0.;
       FORALL_BD f_0[iz][ip][ir]    = 0.;
       BUGREPORT;    

       /* U V are monotonous initially, V with noise */
       /* FORALL   V_0[iz][ip][ir] = ((double)(iz+ data.lnz*data.grid_coords[0]) + 0.5)/(double)data.dims[0]; */
       /* FORALL   U_0[iz][ip][ir] = ((double)(iz+ data.lnz*data.grid_coords[0]) + 0.5)/(double)data.dims[0]; */

       FORALL{
	 //	 V_0[iz][ip][ir] = cyto_blob_distribution_z(iz, nz, &data) / (double)data.dims[0];
	 //	 U_0[iz][ip][ir] = V_0[iz][ip][ir];
	 V_0[iz][ip][ir] = (double)ic_velocity_amplitude*
	   exp(- (x_c - ir)*(x_c - ir)/ w_x )
	   *exp(- (z_c - (nz* data.grid_coords[0] + iz))*(z_c - (nz* data.grid_coords[0] + iz))/ w_z );
	 U_0[iz][ip][ir] = V_0[iz][ip][ir];
       }

       BUGREPORT;
       /* from zero to 1 with zero derivative at 0, note potential is relative to sheath potential!!! */
       FORALL   V_0[iz][ip][ir]*=V_0[iz][ip][ir] /*1./sqrt(2.*M_PI)*sqrt(para.Mi/ME)*/ ;
       FORALL   U_0[iz][ip][ir]*=U_0[iz][ip][ir];

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
   else if (data.restart == START_FROM_FILE ) 
   {
	printf("START_FROM_FILE...\n");

       /******************START FROM .000 FILE  ****************************************/
       fprintf(stderr, "Start from .000 file: nx = %d ny = %d nz = %d \n",nx,ny,nz);
       BUGREPORT;
       data.read_data = TRUE;
       data.ReadAttributes=FALSE;
       BUGREPORT;     

#ifdef FILEIO /* Stefan: Handle HDF4 and HDF5 (complex) */
       int size[]={nz,ny,nx}, offs[]={data.offz,data.offy,data.offx};
       FieldDim_Construct(&D,3,size,offs);

       sprintf(fname,"%s.%03d",data.name_in,data.number);

       FileIO_Read(fname,&F_exp_n_0);

#elif HDF4 /* Mathias: simple HDF4/HDF5 way */
       PH_Read3DFieldbyNumber(exp_n_0,1,&data,&para);
#elif HDF5
       printf("Missing Read Field Function\n");
       exit(0);
#else
       printf("Not avaible yet!\n");
       exit(0);
#endif
       if(ISROOT) fprintf(stderr,"%d cordsys %s\n",__LINE__,data.coordsys);

       if(false)
       for(ir=0;ir<nx;++ir) printf("%f ",exp_n_0[nz/2][ny/2][ir]); printf("\n");

       /* Make sure that n is positive */
       
       FORALL n_0[iz][ip][ir] =  log(fabs(exp_n_0[0][ip][ir])); 
       /* Setup axial homogeneous plasma to start with */

       /* U V are monotonous initially */

       FORALL   V_0[iz][ip][ir] = -1.+ 2./(double)data.dims[0]*((double)(iz+ nz*data.grid_coords[0])+0.5);
       FORALL   U_0[iz][ip][ir] = -1.+ 2./(double)data.dims[0]*((double)(iz+ nz*data.grid_coords[0])+0.5);
       FORALL   w_0[iz][ip][ir] = omega[iz][ip][ir] = f_0[iz][ip][ir]  =  0.;
       fprintf(stderr, "nx = %d ny = %d nz = %d \n",nx,ny,nz);
   }
   else
   {
       fprintf(stderr, "Restart undefined \n");
       MPI_Finalize();
       return 0;
   }

 if(ISROOT) fprintf(stderr,"%d cordsys %s\n",__LINE__,data.coordsys);

 MPI_Bcast(&para,1,para_type,ROOT,MPI_COMM_WORLD);
 MPI_Barrier(MPI_COMM_WORLD);
if(ISROOT) fprintf(stderr,"%d cordsys %s\n",__LINE__,data.coordsys);
 BUGREPORT;
 /****************************************************/
 /*                                                  */
 /* Calculate  dimensionless parameters                             */
 /*                                                  */
 /****************************************************/
 
 FORZX_BD
 {
     r = para.xmin+(data.lnx*data.grid_coords[2]+ir+0.5)*para.dx;
     rcor[ir] = r;
     edrcor[ir] = 1./r;
 }
 
 FORZX_BD hval[iz][ir] = 1./(para.dx); 
 FORZX_BD vval[iz][ir] = edrcor[ir]/(para.dy);
 FORZX  nlval[iz][ir]  = vval[iz][ir]*hval[iz][ir];
 FORYX  norm[ir]       = 1./(hval[0][ir]*vval[0][ir]);

 COMM(if(ISROOT) {
          FORX fprintf(stderr,"hval[%d] = %f \t vval[%d]=%f \t nlval = %f \n", ir,hval[0][ir],ir,vval[0][ir],nlval[0][ir] );
      });

 // Mathias, to avoid that memory couting 

 MPI_Bcast(&para.Mi,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

 MPI_Bcast(&para.rho_s,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
 MPI_Bcast(&para.rho_i,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

 MPI_Bcast(&para.omega_ci,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
 MPI_Bcast(&para.omega_ce,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
 MPI_Bcast(&para.omega_pi,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
 MPI_Bcast(&para.v_thi,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
 MPI_Bcast(&para.v_the,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
 MPI_Bcast(&para.c_s,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
 MPI_Bcast(&para.v_alfven,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
 MPI_Bcast(&para.nu_ei,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
 MPI_Bcast(&para.nu_in,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
 MPI_Bcast(&para.nu_en,1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

 mu     = para.Mi/ME;
 delta  = para.delta;
 nu     = para.nu;/* assume nu to be determined from plasma core while it should be plasma edge)*/
 nu_m   = nu/mu;


 if(ISROOT) fprintf(stderr,"cordsys                            : %s\n",data.coordsys);
 if(ISROOT) fprintf(stderr,"kappan                             : %f\n",para.kappan);
 if(ISROOT) fprintf(stderr,"Electron neutral collisions        : %g\n",para.delta);
 if(ISROOT) fprintf(stderr,"Electron Ion collisions            : %g\n",para.nu);
 if(ISROOT) fprintf(stderr,"Ion neutral collisions             : %g\n",para.sigma);
 if(ISROOT) fprintf(stderr,"Viscous                            : W = %g     N = %g\n",para.mue_w,para.mue_n);
 if(ISROOT) fprintf(stderr,"z coord                            : [%g:%g]  dz =  %g\n",para.zmin,para.zmax,para.dz);
 if(ISROOT) fprintf(stderr,"nprof                              : %f\n",para.nprof);
 if(ISROOT) fprintf(stderr,"source                             : %f\n",para.source);
 if(ISROOT) fprintf(stderr,"limiter                            : %f\n",para.limiter);
 if(ISROOT) fprintf(stderr,"Zeff                               : %f\n",para.Z);
 if(ISROOT) fprintf(stderr,"X:  %f %f  [%f %f]\n",para.xmin,para.xmax,data.range[2][0],data.range[2][1]);
 if(ISROOT) fprintf(stderr,"Y:  %f %f  [%f %f]\n",para.ymin,para.ymax,data.range[1][0],data.range[1][1]);
 if(ISROOT) fprintf(stderr,"Z:  %f %f  [%f %f]\n",para.zmin,para.zmax,data.range[0][0],data.range[0][1]);
     
 BUGREPORT;
 if(ISROOT) FORX fprintf(stderr,"xmin,xmax = [%f,%f], ir  %d r=%f\n",para.xmin,para.xmax,ir,data.coordinate[2][ir]);
 FORALL target_density[iz][ip][ir] = 1.+10.*exp(-rcor[ir]*rcor[ir]*(para.kappan*para.kappan));
 
 BUGREPORT; //VN
   FORALL_BD	exp_n_0[iz][ip][ir] = exp(n_0[iz][ip][ir]) ;
 FORALL_BD	nu_density[iz][ip][ir] = nu*exp_n_0[iz][ip][ir]; //VN
 BUGREPORT; //VN
 if(ISROOT) fprintf(stderr,"%d cordsys %s\n",__LINE__,data.coordsys); 

//VN Code is in evil state already here 
/*
  Code breaks on restart from old file, is not related to read of fields but happens earlier.
  It is also not related to MPI routines
  Code btreaks on accessing negative indices

*/




 

/*
     Put values into xy ghost points, parallel derivative takes care of z boundaries!
*/


   BUGREPORT;   
   PH_Update2dBoundaries(n_0,rbdcndf, nbdra, nbdrb,hval,&data,&para);
   BUGREPORT;

   FORALL_BD	exp_n_0[iz][ip][ir] = exp(n_0[iz][ip][ir]) ;
   BUGREPORT;
   PH_Update2dBoundaries(f_0,rbdcndf, fbdra, fbdrb,hval,&data,&para);
   PH_Update2dBoundaries(U_0,rbdcndu, ubdra, ubdrb,hval,&data,&para);
   PH_Update2dBoundaries(V_0,rbdcndv, vbdra, vbdrb,hval,&data,&para);
   BUGREPORT;
   /* Calculate omega = w  - DX phi  DX n - DY phi DY n*/

   FORALL dxn[iz][ip][ir] = DX(n_0);
   FORALL dyn[iz][ip][ir] = DY(n_0);
   
   FORALL dxf[iz][ip][ir] = DX(f_0);
   FORALL dyf[iz][ip][ir] = DY(f_0);
   BUGREPORT;
   FORALL  omega[iz][ip][ir]  = w_0[iz][ip][ir]-dxn[iz][ip][ir]*dxf[iz][ip][ir]-dyn[iz][ip][ir]*dyf[iz][ip][ir];   
   PH_Update2dBoundaries(omega,rbdcndw, ombdra, ombdrb,hval,&data,&para);
 
   BUGREPORT;

   /* calculate boundary for w  = 0 + d_r phi  d_r n */
   for(iz=0;iz<nz;iz++) 
       for(ip=-1;ip<ny+1;ip++) 
           wbdrb[iz][ip] =   ombdrb[iz][ip] +
               4.* hval[iz][nx-1]* hval[iz][nx-1]*(f_0 [iz][ip][nx] - f_0 [iz][ip][nx-1])*(n_0 [iz][ip][nx] - n_0 [iz][ip][nx-1]);
   BUGREPORT;  
   PH_Update2dBoundaries(w_0,rbdcndw,wbdra,wbdrb,hval,&data,&para);

   BUGREPORT;

   /* Precalculate all parallel derivatives */
  
   if(ISROOT) fprintf(stderr,"%d cordsys %s\n",__LINE__,data.coordsys);
   Cyto_DParallelSOL(IS_IONVELOCITY,dzU,dzzU,U_0,f_0,exp_n_0,V_0,dzV,&data,&para);
   BUGREPORT;

   Cyto_DParallelSOL(IS_POTENTIAL,dzF,dzzF,f_0,f_0,exp_n_0,V_0,dzV,&data,&para);
   BUGREPORT;

    /*
      V-static       - Equation , assume dzN == 0 
    */ 
   FORALL	nu_density[iz][ip][ir] = nu*exp_n_0[iz][ip][ir];
   BUGREPORT;

   FORALL	vstatic[iz][ip][ir]  =  1./(delta +nu_density[iz][ip][ir])
       *(mu*(dzF[iz][ip][ir] /*-dzN[iz][ip][ir]*/) +nu_density[iz][ip][ir]*U_0[iz][ip][ir]);

#ifdef VSTATIC
// Stefan   
//   FORALL	V_0[iz][ip][ir] = vstatic[iz][ip][ir] ;
#endif 

   Cyto_DParallelSOL(IS_ELECTRONVELOCITY,dzV,dzzV,V_0,f_0,exp_n_0,V_0,dzV,&data,&para);
   /* Remember density parallel derivative depends on electron velocity */
   Cyto_DParallelSOL(IS_DENSITY,dzN,dzzN,n_0,f_0,exp_n_0,V_0,dzV,&data,&para);

   BUGREPORT;
   /*********************************************************/
   /* 
      
   set up density source 
   
   */
   
   /* Localisation of density source in first quarter of device */
   // Stefan
   /*
   Lz = 0.25*(para.zmax-para.zmin);
   FORALL_BD tf[iz][ip][ir] = para.nprof*(exp(-rcor[ir]*rcor[ir]/(para.kappan*para.kappan))*
                                          exp(-((double)(nz*data.POS[0]+iz)+0.5)*para.dz/Lz 
                                              *((double)(nz*data.POS[0]+iz)+0.5)*para.dz/Lz));
   */
#ifdef mathias_debug

   double *myvec; myvec=(double*)calloc(nz+2, sizeof(double));
   double *myvec2; myvec2=(double*)calloc(nz+2, sizeof(double));

   FORZ_BD myvec[iz+1]  = exp(-((double)(iz)-0.5) *((double)(iz)-0.5)*para.dz*para.dx/20./20.);
   FORZ_BD myvec2[iz+1] = exp(-((double)(nz*data.grid_coords[0]+iz)-0.5)*para.dz/20
			    *((double)(nz*data.grid_coords[0]+iz)-0.5)*para.dx/20);
#endif
   // Mathias: This was the orig. batch_job definition, but it's not proper set for >1 CPUs
   //   FORALL_BD tf[iz][ip][ir] = para.nprof*(exp(-rcor[ir]*rcor[ir]*(para.kappan*para.kappan))*
   //					  exp(-((double)(iz)-0.5) *((double)(iz)-0.5)*para.dz*para.dx/20./20.));

   double val_z, val_x;

   FORALL_BD{

     val_x = rcor[ir]*para.kappan;
     val_z = cyto_blob_distribution_z(iz, nz, &data)* para.dz / 10.0;
     
     tf[iz][ip][ir] = para.nprof*(  exp(- val_x* val_x )
				  * exp(- val_z* val_z ) );

     tf[iz][ip][ir] = 0.0;

   }

   //   FORZ printf("%lg \n", cyto_blob_distribution_z(iz, nz, &data));
   BUGREPORT;
   
   /*
     Write initial condition
   */
   MPI_Barrier(MPI_COMM_WORLD);
   if(ISROOT) fprintf(stderr,"%d cordsys %s\n",__LINE__,data.coordsys);

#ifdef FILEIO
   FileIO_Write_FP("ini.000", TRUE,FP_ini);
   FileIO_Write_FP("ini.000", FALSE, (const field**) SP_para); // Save Structure Values
#else
   /*
   PH_3D_Write(omega,"Vorticity","ini",0,&data,&para,TRUE);
   PH_3D_Write(n_0,"Density","ini",0,&data,&para,FALSE);
   PH_3D_Write(f_0,"Potential","ini",0,&data,&para,FALSE);
   PH_3D_Write(U_0,"UIon","ini",0,&data,&para,FALSE);
   PH_3D_Write(V_0,"Velocity","ini",0,&data,&para,FALSE);
   PH_3D_Write(tf,"Source","ini",0,&data,&para,FALSE);
   PH_3D_Write(dzU,"dzU","ini",0,&data,&para,FALSE);
   PH_3D_Write(dzN,"dzN","ini",0,&data,&para,FALSE);
   PH_3D_Write(dzF,"dzF","ini",0,&data,&para,FALSE);
   PH_3D_Write(w_0,"H","ini",0,&data,&para,FALSE);
   */
   /* debug_mean(n_0, "N", 0, nz, ny, nx, */
   /* 	      data.this_process, data.num_procs); */

   PH_3D_Write(exp_n_0,"n",data.name_out,0,&data,&para,TRUE);
   PH_3D_Write(U_0,"UIon",data.name_out,0,&data,&para,FALSE);
   PH_3D_Write(dzN,"dzn",data.name_out,0,&data,&para,FALSE);
   PH_3D_Write(w_0,"H",data.name_out,0,&data,&para,FALSE);
   PH_3D_Write(f_0,"Potential",data.name_out,0,&data,&para,FALSE);
   PH_3D_Write(n_0,"N",data.name_out,0,&data,&para,FALSE);
   PH_3D_Write(V_0,"VElectron",data.name_out,0,&data,&para,FALSE);
   PH_3D_Write(omega,"Vorticity",data.name_out,0,&data,&para,FALSE);
   PH_3D_Write(dzF ,"dzf",data.name_out,0,&data,&para,FALSE);
   PH_3D_Write(vstatic,"Vstatic",data.name_out,0,&data,&para,FALSE);
   PH_3D_Write(J,"Current",data.name_out,0,&data,&para,FALSE);
   PH_3D_Write(dyn,"ParForce",data.name_out,0,&data,&para,FALSE);
   PH_3D_Write(dxn,"ParMomentum",data.name_out,0,&data,&para,FALSE);

#endif
   MPI_Barrier(MPI_COMM_WORLD);
   if(ISROOT) fprintf(stderr,"%d cordsys %s\n",__LINE__,data.coordsys);
   
   /* Test of laplace solver */
    /*
   for(m=1;m<ny/2*2/3;m++)
   {
       FORALL dtw_0[iz][ip][ir] = cos(2.*M_PI*m*((double)ip/(double)ny+(double)iz/(double)nz))*exp(-( (double)((ir-nx/2)*(ir-nx/2))/(double)(nx*nx)*100.));
  
       PH_Update2dBoundaries(dtw_0,ZeroDIR, fbdra, fbdrb,hval,&data,&para);
     
       Cyto_Laplace(dtw_1,dtw_0,hval,vval,edrcor,nx,ny,nz,&para,ispolar);


       FORALL dtn_1[iz][ip][ir] = 
	 vval[iz][ir]*vval[iz][ir]*(dtw_0[iz][ip+1][ir]-2.*dtw_0[iz][ip][ir]+dtw_0[iz][ip-1][ir]) 
	 + hval[iz][ir]*hval[iz][ir] * (dtw_0[iz][ip][ir+1] - 2.*dtw_0[iz][ip][ir] +dtw_0[iz][ip][ir-1]) 
	 + edrcor[ir]*0.5*hval[iz][ir]*(dtw_0[iz][ip][ir+1]-dtw_0[iz][ip][ir-1]);



       HHSOLVER(&data,&para,dtw_2,dtw_1,fbdra,fbdrb,hval[0],rcor,rbdcndf,0.,FALSE);
       PH_3D_Write(dtw_0,"PHI_orig",dtw_1,"Vorticity",dtw_2,"RecPhi","laptest",m,&data,&para,TRUE);

       
       FORALL dtn_2[iz][ip][ir] = dtw_1[iz][ip][ir] - dtn_1[iz][ip][ir];
       FORALL dtn_0[iz][ip][ir] = dtw_0[iz][ip][ir] - dtw_2[iz][ip][ir];


       
       FUtils_Write3dSpace(dtn_0,"Diff",dtn_1,"LapFD",dtn_2,"FFT-FD","laptest",m,&data,&para,FALSE);
       
   }
   exit(0);
    */
   /****************************************************/
   /*                                                  */
   /* WHAT ARE THE DEFINES                             */
   /*                                                  */
   /****************************************************/
 
   sprintf(data.desc,"SWITCHES:\n");
   
   strcat(data.desc,"Order of Parallel Derivative - ");
   if(2 == data.offz) 
     strcat(data.desc,"4th.\n");
   else
     strcat(data.desc,"2nd.\n");
   
   strcat(data.desc,"Electron parallel dynamics - ");
#ifdef VSTATIC
   strcat(data.desc,"static\n");
#else
   strcat(data.desc,"dynamic.\n");
#endif
   
   if(ISROOT) fprintf(stderr,"%s",data.desc);
         


/*****************************************************************/
/*                                                               */
/*          Start of Time Loop                                   */
/*                                                               */
/*****************************************************************/
 MPI_Barrier(MPI_COMM_WORLD);
 write_time     = para.time + para.out_time;
 while(para.time < para.end_time )
  {
    if( isnan(w_0[nz/2][ny/2][nx/2]) || isnan(f_0[nz/2][ny/2][nx/2])
	|| isnan(n_0[nz/2][ny/2][nx/2]) || isnan(omega[nz/2][ny/2][nx/2]) ) {
      printf("Error nan occurs at time:%lg \n", para.time);
      exit(0);
    }



      BUGREPORT;
	
      FORALL dzVN[iz][ip][ir] =  dzV[iz][ip][ir] +  V_0[iz][ip][ir]*dzN[iz][ip][ir]; 
      FORALL dzUN[iz][ip][ir] =  dzU[iz][ip][ir] +  U_0[iz][ip][ir]*dzN[iz][ip][ir]; 

      FORALL dxn[iz][ip][ir]  = DX(n_0);
      FORALL dyn[iz][ip][ir]  = DY(n_0);
     
      FORALL dxf[iz][ip][ir]  = DX(f_0);
      FORALL dyf[iz][ip][ir]  = DY(f_0);

      /* zero boundaries for dxf, dyf */
      PH_Update2dBoundaries(dxf,ZeroNEU, fbdra, fbdrb,hval,&data,&para);
      PH_Update2dBoundaries(dyf,ZeroDIR, fbdra, fbdrb,hval,&data,&para);

BUGREPORT;

      /*********************************************/
      /*       
               Electron Density - Equation                  
      */ 
      /*********************************************/
      
      Util_3DArakawaNl(dtn_0,f_0,n_0,nlval,nz,nx,ny); 
      FORALL dtn_0[iz][ip][ir]   +=  -dzVN[iz][ip][ir];

      /* density source */
      FORALL    dtn_0[iz][ip][ir]+= tf[iz][ip][ir]/exp_n_0[iz][ip][ir];

      /* This localised source is absent in GW code */
      /*
        for(iz=0;iz<nz/16;iz++) { 
          FORYX dtn_0[iz][ip][ir]+= para.nprof*(1.-(double)(16*16*iz*iz)/(double)(nz*nz))*(target_density[iz][ip][ir]/exp_n_0[iz][ip][ir]-1.);
      }
      */


      /* ionisation, keeps log density from falling below yero  */

      // Stefan
      para.source = 1.;

      FORALL    dtn_0[iz][ip][ir]-=  para.source* (n_0[iz][ip][ir]-fabs(n_0[iz][ip][ir])); 
      /* VN introduced parameter source 15032010*/

#ifdef PARALLEL_DAMPING
      FORALL       dtn_0[iz][ip][ir]+=para.limiter*(dzzN[iz][ip][ir] +  dzN[iz][ip][ir] *dzN[iz][ip][ir]);
#endif
      
      /* Explicit part of damping */
      FORALL dtn_0[iz][ip][ir] += para.mue_n*(dxn[iz][ip][ir]*dxn[iz][ip][ir] + dyn[iz][ip][ir]*dyn[iz][ip][ir]);

      Util_3DSsTimeStep(iter,nx,ny,nz,para.dt,para.mue_n,&lamda, 
                        n,n_0,n_1,n_2,dtn_0,dtn_1,dtn_2);
      
      HHSOLVER(&data,&para,n_2,n,nbdra,nbdrb,hval[0],rcor,rbdcndn,lamda,TRUE);
            
      /*Calculate dt ln n .....with implicit term */
      
      Cyto_Laplace(dtw_0,n_0,hval,vval,edrcor,nx,ny,nz,&para,ispolar);
      FORALL  dtw_0[iz][ip][ir]*= para.mue_n;
      FORALL dzzN[iz][ip][ir] = dtn_0[iz][ip][ir] + dtw_0[iz][ip][ir];
      
      BUGREPORT;
      
      /*********************************************/
      /*
        
      Global Vorticity equation needs to be evaluated AFTER dt_n_0 is known
      
      */
      /*********************************************/
      /*
        Note that w = omega + grad N grad phi 
      */
      Util_3DArakawaNl(dtw_0,f_0,omega,nlval,nz,nx,ny);

      /* Parallel terms, Z = 1 in GW code */
      FORALL dtw_0[iz][ip][ir] +=  (dzUN[iz][ip][ir] - dzVN[iz][ip][ir])/para.Z;

      
      /* Collisional Term: Pedersen current */
      FORALL dtw_0[iz][ip][ir] +=  -para.sigma*omega[iz][ip][ir];
      
#ifdef PARALLEL_DAMPING
      FORALL dtw_0[iz][ip][ir] +=  -para.limiter*dzzF[iz][ip][ir];
#endif
      
      /* \nabla (d_t n) \nabla \phi term */
      /* Boundaries for dt_n0 */

      /* Was zeroNEU, but zeroDIR in GW version */
      PH_Update2dBoundaries(dzzN,ZeroDIR, fbdra, fbdrb,hval,&data,&para);


      FORALL dtw_0[iz][ip][ir] += DX(dzzN)*dxf[iz][ip][ir]+DY(dzzN)*dyf[iz][ip][ir];
      
      /* Nonlinear Term: -grad N {phi, grad phi }, sign like in main nonlinear term!*/ 
      
      Util_3DArakawaNl(dtdxf,f_0,dxf,nlval,nz,nx,ny);        
      Util_3DArakawaNl(dtdyf,f_0,dyf,nlval,nz,nx,ny);
      FORALL dtw_0[iz][ip][ir] += dxn[iz][ip][ir]*dtdxf[iz][ip][ir]+dyn[iz][ip][ir]*dtdyf[iz][ip][ir];
      
      /*   
                   Calculate \nabla (sigma  phi - mue w) 
      */
      BUGREPORT;




      FORALL res[iz][ip][ir] =   para.mue_w*omega[iz][ip][ir] -para.sigma*f_0[iz][ip][ir];
      
      
      FORALL dtw_0[iz][ip][ir] +=   dxn[iz][ip][ir]*DX(res) + dyn[iz][ip][ir]*DY(res);
              
      /* Step */         
      Util_3DSsTimeStep(iter,nx,ny,nz,para.dt,para.mue_w,&lamda, 
                        w,w_0,w_1,w_2,dtw_0,dtw_1,dtw_2);

      BUGREPORT;
/*
      // Damping only on vorticity, omega, thus substract grad n grad phi term , 
       //   not done in GW verison 
      FORALL dtU_0[iz][ip][ir]  = dxn[iz][ip][ir]*dxf[iz][ip][ir]+ dyn[iz][ip][ir]*dyf[iz][ip][ir];

      FORALL  w[iz][ip][ir]    += -dtU_0[iz][ip][ir];

       HHSOLVER(&data,&para,w_2,w,ombdra,ombdrb,hval[0],rcor,rbdcndw,lamda,TRUE);   
       //add grad n grad phi term again  

      FORALL  w_2[iz][ip][ir]  +=  dtU_0[iz][ip][ir]; 
*/

      // VN 15032010 GWV0 Version 
      // calculate boundary for w  = om + d_r phi  d_r n 
      for(iz=-1;iz<nz+1;iz++) 
              for(ip=-1;ip<ny+1;ip++) 
                  wbdrb[iz][ip] =   ombdrb[iz][ip] + 2.*hval[0][nx-1]* 2.*hval[0][nx-1]*(f_0 [iz][ip][nx] - f_0 [iz][ip][nx-1])*(n_0 [iz][ip][nx] - n_0 [iz][ip][nx-1]) ;
      HHSOLVER(&data,&para,w_2,w,wbdra,wbdrb,hval[0],rcor,rbdcndw,lamda,TRUE);  
   
            
      BUGREPORT;
      
      /*********************************************/
      /*       Ion velocity - Equation             */ 
      /*********************************************/
      
      Util_3DArakawaNl(dtU_0,f_0,U_0,nlval,nz,nx,ny);

      /*Nonlinear parallel advection */ 
      FORALL dtU_0[iz][ip][ir]  += - padvection*dzU[iz][ip][ir]*U_0[iz][ip][ir]; 
      BUGREPORT;

      FORALL dtU_0[iz][ip][ir]  +=  -para.Z * dzF[iz][ip][ir]; /* Z = 1 in GW version */
      
      /* Neutral collision Term  */
      FORALL dtU_0[iz][ip][ir]  += -para.sigma*U_0[iz][ip][ir] ;
      
      /* Resistivity  */

//      FORALL dtU_0[iz][ip][ir]  += nu_density[iz][ip][ir]/mu*(V_0[iz][ip][ir]-U_0[iz][ip][ir]);
// VN 15032010: Not in GWV0 code   
              
#ifdef PARALLEL_DAMPING
//      FORALL       dtU_0[iz][ip][ir]+=  0.01*para.limiter*dzzU[iz][ip][ir];
// VN 15032010: one GWV0 code
   
      FORALL       dtU_0[iz][ip][ir]+=  para.limiter*dzzU[iz][ip][ir];  
  
#endif
    
      // Util_3DSsTimeStep(iter,nx,ny,nz,para.dt,para.mue_n,&lamda, 
      //                         res,U_0,U_1,U_2,dtU_0,dtU_1,dtU_2);

        
      // HHSOLVER(&data,&para,U_2,res,ubdra,ubdrb,hval[0],rcor,rbdcndu,lamda,TRUE); 
      // VN 15032010: one GWV0 code

      Util_3DSsTimeStep(iter,nx,ny,nz,para.dt,para.mue_n,&lamda, 
                                U_2,U_0,U_1,U_2,dtU_0,dtU_1,dtU_2);
    
      /*********************************************/
      /*      Electron  Velocity - Equation                 */
      /*********************************************/
      BUGREPORT;
      Util_3DArakawaNl(dtV_0,f_0,V_0,nlval,nz,nx,ny);
      FORALL dtV_0[iz][ip][ir]   +=  mu*(dzF[iz][ip][ir]-dzN[iz][ip][ir]);
      
      /*FORALL dtV_0[iz][ip][ir] += - padvection*dzV[iz][ip][ir]*V_0[iz][ip][ir];*/
      
      /* Neutral collision Term  */
      FORALL dtV_0[iz][ip][ir]   += -para.delta*V_0[iz][ip][ir];   
      
      /* Coloumb collision Term */
      
      FORALL dtV_0[iz][ip][ir]   += -nu_density[iz][ip][ir]*(V_0[iz][ip][ir]-U_0[iz][ip][ir]) ;
     
      BUGREPORT;
#ifdef PARALLEL_DAMPING
      FORALL dtV_0[iz][ip][ir]   += 10.*para.limiter*dzzV[iz][ip][ir];
#endif

      BUGREPORT;
      Util_3DSsTimeStep(iter,nx,ny,nz,para.dt,para.mue_w,&lamda, 
                        res,V_0,V_1,V_2,dtV_0,dtV_1,dtV_2);
            
      HHSOLVER(&data,&para,V_2,res,vbdra,vbdrb,hval[0],rcor,rbdcndv,lamda,TRUE); 

    
      /*********************************************************************/
      /* Exchange pointers                                                 */       
      /*********************************************************************/
      
      tmp = w_2;     w_2 = w_1;     w_1 = w_0;     w_0 = tmp;
      tmp = n_2;     n_2 = n_1;     n_1 = n_0;     n_0 = tmp;
      tmp = U_2;     U_2 = U_1;     U_1 = U_0;     U_0 = tmp;   
      tmp = V_2;     V_2 = V_1;     V_1 = V_0;     V_0 = tmp;   
      
      tmp = dtw_2; dtw_2 = dtw_1; dtw_1 = dtw_0; dtw_0 = tmp;
      tmp = dtn_2; dtn_2 = dtn_1; dtn_1 = dtn_0; dtn_0 = tmp;
      tmp = dtU_2; dtU_2 = dtU_1; dtU_1 = dtU_0; dtU_0 = tmp;
      tmp = dtV_2; dtV_2 = dtV_1; dtV_1 = dtV_0; dtV_0 = tmp;
      

      MPI_Barrier(MPI_COMM_WORLD);
      
  
      /*********************************************************************/
      /* Boundaries                                                        */       
      /*********************************************************************/
  
      BUGREPORT;
      // FORZY_BD  nbdrb[iz][ip] = n_0[iz][ip][nx-1];
      // VN 15032010: Not in GWV0 code   

      BUGREPORT; 
      PH_Update2dBoundaries(n_0,rbdcndn, nbdra, nbdrb,hval,&data,&para);

      BUGREPORT; 
      FORALL_BD exp_n_0[iz][ip][ir] = exp(n_0[iz][ip][ir]) ;
      BUGREPORT; 
      FORZY_BD  enbdrb[iz][ip] = exp(nbdrb[iz][ip]);

      PH_Update2dBoundaries(U_0,rbdcndu, ubdra, ubdrb,hval,&data,&para);
      BUGREPORT;    
      /* calculate boundary for w  = om + d_r phi  d_r n */
/*      for(iz=0;iz<nz;iz++) 
          for(ip=-1;ip<=ny;ip++) 
              wbdrb[iz][ip] =   ombdrb[iz][ip] 
                  + 2.*hval[iz][nx-1]* 2.*hval[iz][nx-1]*(f_0[iz][ip][nx] - f_0[iz][ip][nx-1])
                  *(n_0[iz][ip][nx] - n_0[iz][ip][nx-1]) ; 
                  
*/
      for(iz=0;iz<nz;iz++) for(ip=-1;ip<=ny;ip++) 
          wbdrb[iz][ip] = 0.;
      // VN 15032010: Boundary incorrect in GWV0 code 

     

      PH_Update2dBoundaries(w_0,rbdcndw, wbdra, wbdrb,hval,&data,&para);
   
      BUGREPORT;   
      /*
        Calculate Potential and vorticity from 
        density n and GlobalVorticity  W = (Dxx + Dyy) phi + \nabla n \nabla phi 
        
        Iterative solution of equation 
        W = (Dxx + Dyy) phi + \nabla n nabla phi 
      */
      
      FORALL dxn[iz][ip][ir] = DX(n_0);
      FORALL dyn[iz][ip][ir] = DY(n_0);     
       BUGREPORT; 
      for(i=0;i<3;i++)
      {
          FORALL  omega[iz][ip][ir] = w_0[iz][ip][ir]- dxn[iz][ip][ir]*dxf[iz][ip][ir]-dyn[iz][ip][ir]*dyf[iz][ip][ir];
          PH_Update2dBoundaries(omega,rbdcndw, ombdra, ombdrb,hval,&data,&para);
          HHSOLVER(&data,&para,f_0,omega,fbdra,fbdrb,hval[0],rcor,rbdcndf,0.,FALSE);  
          PH_Update2dBoundaries(f_0,rbdcndf, fbdra, fbdrb,hval,&data,&para);
          FORALL dxf[iz][ip][ir] = DX(f_0);
          FORALL dyf[iz][ip][ir] = DY(f_0);
      }  
      /* End of iterative solve */ 
      BUGREPORT;  
  
      /* calculate boundary for w  = om + d_r phi  d_r n */
  
      /* for(iz=0;iz<nz;iz++) 
          for(ip=-1;ip<=ny;ip++) 
          wbdrb[iz][ip] =   ombdrb[iz][ip] + 
          2.*hval[iz][nx-1]* 2.*hval[iz][nx-1]*(f_0[iz][ip][nx] - f_0[iz][ip][nx-1])
          *(n_0[iz][ip][nx] - n_0[iz][ip][nx-1]) ;
          */
      for(iz=0;iz<nz;iz++) for(ip=-1;ip<=ny;ip++) 
          wbdrb[iz][ip] = 0.;
      // VN 15032010: Boundary incorrect in GWV0 code 

     


      PH_Update2dBoundaries(w_0,rbdcndw, wbdra, wbdrb,hval,&data,&para);
      BUGREPORT;
        
      /* Parallel derivatives */   
      Cyto_DParallelSOL(IS_POTENTIAL,dzF,dzzF,f_0,f_0,exp_n_0,V_0,dzV,&data,&para);
      Cyto_DParallelSOL(IS_DENSITY,dzN,dzzN,n_0,f_0,exp_n_0,V_0,dzV,&data,&para);
      Cyto_DParallelSOL(IS_IONVELOCITY,dzU,dzzU,U_0,f_0,exp_n_0,V_0,dzV,&data,&para);
      BUGREPORT;
      /*
        V-static       - Equation 
      */
      FORALL nu_density[iz][ip][ir] = nu*exp_n_0[iz][ip][ir];
      FORALL	vstatic[iz][ip][ir]    = 
	(mu*(dzF[iz][ip][ir]- dzN[iz][ip][ir]) + nu_density[iz][ip][ir]*U_0[iz][ip][ir])/(delta +  nu_density[iz][ip][ir]);

#ifdef VSTATIC
      FORALL	V_0[iz][ip][ir]        = vstatic[iz][ip][ir];
#endif

      PH_Update2dBoundaries(V_0,rbdcndv, vbdra, vbdrb,hval,&data,&para);
      Cyto_DParallelSOL(IS_ELECTRONVELOCITY,dzV,dzzV,V_0,f_0,n_0,V_0,dzV,&data,&para);
      BUGREPORT;
      Cyto_DParallelSOL(IS_DENSITY,dzN,dzzN,n_0,f_0,exp_n_0,V_0,dzV,&data,&para);

      BUGREPORT;
      /* Calculate Current */
      FORALL J[iz][ip][ir] = exp_n_0[iz][ip][ir]*(V_0[iz][ip][ir]-U_0[iz][ip][ir]);
    
      /*********************************************************************/
      /*   Advance time , iter                                             */       
      /*********************************************************************/
      
      para.time +=para.dt;  
      BUGREPORT;
      
      MPI_Barrier(MPI_COMM_WORLD);   
      if(iter < 1000)
      {
#ifdef CNTFILES        
          fprintf(stderr,"ITER=%d\n",iter);
          
          FORALL dtw_0[iz][ip][ir] =  dxn[iz][ip][ir]*dxf[iz][ip][ir]+dyn[iz][ip][ir]*dyf[iz][ip][ir];
          BUGREPORT;

#ifdef FILEIO
          sprintf(fname,"%s.%03d","cnt",iter);
          FileIO_Write_FP(fname,TRUE,FP_cnt);
#else 
          PH_3D_Write(exp_n_0,"n","cnt",iter,&data,&para,TRUE);
          PH_3D_Write(U_0,"UIon","cnt",iter,&data,&para,FALSE);
          PH_3D_Write(dzV,"dzV","cnt",iter,&data,&para,FALSE);
          PH_3D_Write(w_0,"H","cnt",iter,&data,&para,FALSE);
          PH_3D_Write(f_0,"Potential","cnt",iter,&data,&para,FALSE);
          PH_3D_Write(n_0,"N","cnt",iter,&data,&para,FALSE);
          PH_3D_Write(V_0,"VElectron","cnt",iter,&data,&para,FALSE);
          PH_3D_Write(omega,"Vorticity","cnt",iter,&data,&para,FALSE);
          PH_3D_Write(dtw_0,"DiffHW","cnt",iter,&data,&para,FALSE);
          PH_3D_Write(vstatic,"VEstatic","cnt",iter,&data,&para,FALSE);
          PH_3D_Write(dzN,"dzN","cnt",iter,&data,&para,FALSE);
          PH_3D_Write(dzF,"dzF","cnt",iter,&data,&para,FALSE);
          PH_3D_Write( J ,"Current","cnt",iter,&data,&para,FALSE);		   
#endif
#endif 
          iter++;	  
      }
      else
      {
          COMM( exit(0););
      }
      

      COMM(fprintf(stderr,"t = %f, wt = %f, ot = %d, otm = %d\n",para.time,write_time,ot,para.otmult););
      BUGREPORT;
      
      /*******************************************/
      /*                                         */
      /*                                         */
      /*      Output                             */
      /*                                         */
      /*                                         */
      /*******************************************/
      
      
      /* Only root writes ascii files */ 
      
      if(para.time>= write_time-0.5*para.dt)
      {
	BUGREPORT;
          ot++;
          write_time+= para.out_time; 
          
          /* 
             Check the CFL condition 
          */
          if(  PH_CheckCFL3D(f_0,vr,vp,&data,&para,hval,vval,&cflr,&cflp) == -1) 
          {
#ifdef FILEIO
              FileIO_Write_FP("VPR.000",TRUE,FP_VPR);
#else
              PH_3D_Write(vp,"vp","VPR",0,&data,&para,TRUE);
              PH_3D_Write(vr,"vr","VPR",0,&data,&para,FALSE);
	      PH_3D_Write(f_0,"Pot","VPR",0,&data,&para,FALSE);
#endif
              fprintf(stderr,"CFL violation detected! Aborting.\n");
              exit(0);
          }
          if(isnan(cflr*cflp))
          {   
              fprintf(stderr,"NaNs detected! Aborting.\n");
              exit(0);
          }
          
          /* Density fluctuations scaled to n_background */
          
          FORALL_BD res[iz][ip][ir] = 0.;
          FORALL res[iz][0][ir] += exp_n_0[iz][ip][ir]/(double)ny;
          FORALL dxn[iz][ip][ir] = (exp_n_0[iz][ip][ir]-res[iz][0][ir])/res[iz][0][ir];
          
          
          
          /**************************************************************/
          /* 
             Output 3D 
          */
          /**************************************************************/
          if((ot >= para.otmult) /*&& para.time < 2000*/)
          {
              ot = 0;
              data.number++;

              FORALL dxn[iz][ip][ir] = exp_n_0[iz][ip][ir]*(ME*V_0[iz][ip][ir]+para.Mi*U_0[iz][ip][ir]);
              FORALL dyn[iz][ip][ir] = dzF[iz][ip][ir]-dzN[iz][ip][ir];

              BUGREPORT;

#ifdef FILEIO
              sprintf(fname,"%s.%03d",data.name_out,data.number);
              FileIO_Write_FP(fname,TRUE,FP_write);
#else
	      /* debug_mean(n_0, "N", data.number, nz, ny, nx, */
	      /* 		 data.this_process, data.num_procs); */

              PH_3D_Write(exp_n_0,"n",data.name_out,data.number,&data,&para,TRUE);
              PH_3D_Write(U_0,"UIon",data.name_out,data.number,&data,&para,FALSE);
              PH_3D_Write(dzN,"dzn",data.name_out,data.number,&data,&para,FALSE);
              PH_3D_Write(w_0,"H",data.name_out,data.number,&data,&para,FALSE);
              PH_3D_Write(f_0,"Potential",data.name_out,data.number,&data,&para,FALSE);
              PH_3D_Write(n_0,"N",data.name_out,data.number,&data,&para,FALSE);
              PH_3D_Write(V_0,"VElectron",data.name_out,data.number,&data,&para,FALSE);
              PH_3D_Write(omega,"Vorticity",data.name_out,data.number,&data,&para,FALSE);
              PH_3D_Write(dzF ,"dzf",data.name_out,data.number,&data,&para,FALSE);
              PH_3D_Write(vstatic,"Vstatic",data.name_out,data.number,&data,&para,FALSE);
	      PH_3D_Write(J,"Current",data.name_out,data.number,&data,&para,FALSE);
              PH_3D_Write(dyn,"ParForce",data.name_out,data.number,&data,&para,FALSE);
              PH_3D_Write(dxn,"ParMomentum",data.name_out,data.number,&data,&para,FALSE);
#endif
              if(ISROOT) fprintf(stderr,"Wrote %s.%03d\n",data.name_out,data.number);
              
              BUGREPORT;
          }
          
          /**************************************************************/
          /* 
             Integrals
          */
          /**************************************************************/         
          
          /* Energy */
          FORALL dtw_0[iz][ip][ir] = -f_0[iz][ip][ir]*omega[iz][ip][ir];
          reducing[0] =  PH_3DIntegral_noreduce(dtw_0,1,norm,nz,ny,nx); /* para.energy */ 
          reducing[1] =  PH_3DIntegral_noreduce(f_0,1,norm,nz,ny,nx); /* totalf */ 
          reducing[2] =  PH_3DIntegral_noreduce(omega,2,norm,nz,ny,nx); /* para.vorticity */ 
          reducing[3] =  PH_3DIntegral_noreduce(omega,1,norm,nz,ny,nx); /* para.circulation */ 
          reducing[4] =  PH_3DIntegral_noreduce(exp_n_0,1,norm,nz,ny,nx); /* gamma_n */ 
        


          BUGREPORT;
          /* Calculate current through plates */
          reducing[5] =  PH_3DIntegral_noreduce(J,1,norm,nz,ny,nx);

          /* Calculate parallel density transport to plates */
          FORALL dtw_0[iz][ip][ir] = exp_n_0[iz][ip][ir]*V_0[iz][ip][ir];
          if(data.grid_coords[0] == 0 ) Util_Integral(dtw_0[0],1, norm,ny,nx,&reducing[6]);
          else if(data.grid_coords[0] ==  (data.N[0]-1))  Util_Integral(dtw_0[nz-1],1, norm,ny,nx,&reducing[6]);
          else reducing[6] = 0.;
       
          /* current loss */
          if(data.grid_coords[0] == 0 ) Util_Integral(J[0],1, norm,ny,nx,&reducing[7]);
          else if(data.grid_coords[0] ==  (data.N[0]-1))  Util_Integral(J[nz-1],1, norm,ny,nx,&reducing[7]);
          else reducing[7] = 0.;

          BUGREPORT;
          MPI_Allreduce((void *)&reducing[0], (void *) &result_vector[0], 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

          para.energy             = result_vector[0];
          totalf                  = result_vector[1];
          para.vorticity          = result_vector[2];
          para.circulation        = result_vector[3];
          gamma_n                 = result_vector[4];
          Jtotal                  = result_vector[5];
          G_par                   = result_vector[6];
          J_par                   = result_vector[7];
       

          BUGREPORT;   
          
          if(ErhFirstWrite && ISROOT)
          {
              
              output = fopen(data.erhname,"a");
              fprintf(output,"# nx = %d ny = %d nz = %ld\n",nx,ny,data.dims[0]);
              fprintf(output,"# Gamma %f; sigma %f; nu %f; alpha %f\n",para.gamma,para.sigma,para.nu,para.alpha);
              fprintf(output,"# mue_w %f; mue_n %f; kappan %f\n",para.mue_w,para.mue_n,para.kappan);
              fprintf(output,"# source %f; k0 %f; shear %f\n#\n",para.source,para.k0,para.shear);
              fprintf(output,"# 1:t\t 2:energy\t 3:vorticity\t 4:<n>"
                      "\t 5: <f>\t 6:circulation\t 7:cflr\t 8:cflp"
                      "\t 9:k_eff \t 10:J \t 11:G_par \t 12:J_par"
                      "\t 13: impulse \t 14:nloss \t 15:ipol \t 16:dtimpulse\n");
              fclose(output);


             output = fopen(profilename,"a");
	     fprintf(output,"# nx = %d ny = %d nz = %ld\n",nx,ny,data.dims[0]);
              fprintf(output,"# Gamma %f; sigma %f; nu %f; alpha %f\n",para.gamma,para.sigma,para.nu,para.alpha);
              fprintf(output,"# mue_w %f; mue_n %f; kappan %f\n",para.mue_w,para.mue_n,para.kappan);
              fprintf(output,"# source %f; k0 %f; shear %f\n#\n",para.source,para.k0,para.shear);
              fprintf(output,"# 1: r\t 2: imp_par \t3: <w>  \t 5: vpar \t 6:Upar\n");
              fclose(output);

              ErhFirstWrite=FALSE;
          }
          
          
          if(ISROOT)
          {
              output = fopen(data.erhname,"a");
              fprintf(output,"%g \t %g \t% g \t %g \t",
                      para.time,para.energy,para.vorticity,gamma_n);
              fprintf(output,"%g \t %g \t% g \t %g \t ",
                      totalf,para.circulation,cflr,cflp);
              fprintf(output,"%g \t %g \n",
                      sqrt(para.vorticity/para.energy),Jtotal);
              fclose(output);
              
          }

#undef PARALLELWORKDONE              
#ifdef PARALLELWORKDONE       
        /* profiles do not work in parallel ! */ 
              /*Momentum profile */
              FORALL dtw_0[iz][ip][ir] = exp_n_0[iz][ip][ir]*(ME*V_0[iz][ip][ir]+para.Mi*U_0[iz][ip][ir]);
              FORX  dtw_0[-1][0][ir]=0.;
              FORALL dtw_0[-1][0][ir]+=dtw_0[iz][ip][ir];
              FORX  dtw_0[-1][0][ir]/=(double)(ny*nz);
              
	      
              /* V-electron profile */
              
              FORX  dtV_0[-1][0][ir]=0.;
              FORALL dtV_0[-1][0][ir]+=V_0[iz][ip][ir];
              FORX  dtV_0[-1][0][ir]/=(double)(ny*nz); 
	      
              /* U-ion Profile */ 
              
              FORX  dtU_0[-1][0][ir]=0.;
              FORALL dtU_0[-1][0][ir]+=U_0[iz][ip][ir];
              FORX  dtU_0[-1][0][ir]/=(double)(ny*nz); 
              
              
              /* Parallel Reynoldstress */
              FORX  par_rey_stress [ir]=0.;
              FORALL par_rey_stress [ir]+=V_0[iz][ip][ir]*vr[iz][ip][ir];
              
              output = fopen(profilename,"a");
              FORX {fprintf(output,"%g \t %g \t %g  \t %g \t %g \t %g \n",
                            rcor[ir], dtw_0[-1][0][ir], 
                            dtn_0[-1][0][ir],dtV_0[-1][0][ir],dtU_0[-1][0][ir], par_rey_stress [ir]);}
              fprintf(output,"\n");
              fclose(output);
     
#endif
	      BUGREPORT;
	      MPI_Barrier(MPI_COMM_WORLD);
      }
      BUGREPORT;
  }
 BUGREPORT;
 
 exit(0);
}

void Cyto_no_sheath_condition(int identity,double ***result,double ***result2,
			   double ***val,double ***f,double ***exp_n,double ***v,double ***dzv,
			   HDF_DS *data, PARA *para)
{
  static int FIRST= TRUE;
  static int nx,ny,nz,offz;
  register int iz,ip,ir;
  double bdval,tval;
  
  if(FIRST){
    nx = data->lnx;
    ny = data->lny;
    nz = data->lnz;
    
    offz = data->offz;
    
    FIRST = FALSE;
  }

  /* Lower Boundary - A */
  if(data->grid_coords[0] == 0) // First Processor
    switch(identity)
      {
	
      case IS_CURRENT:
	FORYX_BD            val[-1][ip][ir] = -val[0][ip][ir];
	break;
	
      case  IS_ELECTRONVELOCITY: /* Zero */
	FORYX_BD 		  val[-1][ip][ir] = val[0][ip][ir];
	break;
      
      case  IS_IONVELOCITY: /* Zero */
	FORYX_BD		  val[-1][ip][ir] = val[0][ip][ir];
	break;
	
      case IS_DENSITY:
	FORYX_BD		  val[-1][ip][ir] = val[0][ip][ir];
	break;
	
      case IS_POTENTIAL:
	FORYX_BD		  val[-1][ip][ir] =  val[0][ip][ir];
	break;
	
      case IS_TEMPERATURE:
	break;
	
      default:
	FORYX_BD		  val[-1][ip][ir] =  val[0][ip][ir];
      }

  /* Upper Boundary - B */
  if(data->grid_coords[0] == (data->N[0]-1)) // Last Processor
    switch(identity)
      {
	
      case IS_CURRENT:
	FORYX_BD	          val[nz][ip][ir] = -val[nz-1][ip][ir];
	break;
	
      case  IS_ELECTRONVELOCITY: /* Zero */
	FORYX_BD 		  val[nz][ip][ir] = val[nz-1][ip][ir];
	break;
	
      case  IS_IONVELOCITY: /* Zero */
	FORYX_BD		  val[nz][ip][ir] = val[nz-1][ip][ir];
	break;
	
      case IS_DENSITY:
	FORYX_BD		  val[nz][ip][ir] =  val[nz-1][ip][ir];
	break;
	
      case IS_POTENTIAL:
	FORYX_BD		  val[nz][ip][ir] =  val[nz-1][ip][ir];
	break;
	
      case IS_TEMPERATURE:
	break;
	
      default:
	FORYX_BD		  val[nz][ip][ir] =  val[nz-1][ip][ir];
      }

}

void Cyto_sheath_condition(int identity,double ***result,double ***result2,
			   double ***val,double ***f,double ***exp_n,double ***v,double ***dzv,
			   HDF_DS *data, PARA *para)
{
  static int FIRST= TRUE;
  static int nx,ny,nz,offz;
  register int iz,ip,ir;
  double bdval,tval;
  
  if(FIRST){
    nx = data->lnx;
    ny = data->lny;
    nz = data->lnz;
    
    offz = data->offz;
    
    FIRST = FALSE;
  }


  BUGREPORT;
  /* Lower Boundary - A */
  if(data->grid_coords[0] == 0) // First Processor
    switch(identity)
      {

      case IS_CURRENT:
	/* Velocity is positive as electrons leave the plasma */
	/* if f=0 the velocity is 1 */                 
  	FORYX_BD {
	  bdval =    2.*( exp_n[0][ip][ir]*(1. - exp(-f[0][ip][ir]) ));
	  val[-1][ip][ir] =  bdval-val[0][ip][ir];
	}
	break;
	
      case  IS_ELECTRONVELOCITY:
	/* Effective Electron speed is given by the sheath potential */
	/* Velocity is positive as electrons leave the plasma */
	/* if f=0 the velocity is unity */                 
	
	FORYX_BD{ 
	  bdval =    2. * exp(-f[0][ip][ir]);
	  val[-1][ip][ir] =  bdval-val[0][ip][ir];
	}
	break;
    
      case  IS_IONVELOCITY:
	/* Ion speed is unity into the limiter */
	FORYX_BD{	
	  bdval=2.;
	  val[-1][ip][ir] = bdval-val[0][ip][ir];
	}
	break;
      
      case IS_DENSITY:
	/* density flux remains constant at the sheath that is dz N = -1/V dz V*/
	FORYX_BD{
	  tval = 0.;
	  val[-1][ip][ir] = tval + val[0][ip][ir];
	} 
	break;
        
      case IS_POTENTIAL:
	/* No gradient in potential, sheeth works via parallel current */
	FORYX_BD  val[-1][ip][ir] = val[0][ip][ir];
	break;
      
      case IS_TEMPERATURE:
	break;
      
      default:
	/* No gradient  */
	FORYX_BD  val[-1][ip][ir] = val[0][ip][ir];
	
      }


  /* Upper Boundary - B */
  if(data->grid_coords[0] == (data->N[0]-1)) // Last Processor
    switch(identity)
      {
      
      case IS_CURRENT:
	/* Velocity is positive as electrons leave the plasma */
	/* if f=0 the velocity is 1 */                 
  	FORYX_BD {
	  bdval =    2.*( exp_n[nz-1][ip][ir]*(1.-exp(-f[nz-1][ip][ir]) ));
	  val[nz][ip][ir] =  bdval-val[nz-1][ip][ir];
	}
	break;
	
      case  IS_ELECTRONVELOCITY:
	/* Effective Electron speed is given by the sheath potential */
	/* Velocity is positive as electrons leave the plasma */
	/* if f=0 the velocity is unity */                 
	FORYX_BD {
	  bdval =    2. * exp(-f[nz-1][ip][ir]);
	  val[nz][ip][ir] =  bdval-val[nz-1][ip][ir];
	}
	break;
    
      case  IS_IONVELOCITY:
	/* Ion speed is unity into the limiter */
	FORYX_BD{
	  bdval=2.;
	  val[nz][ip][ir] = bdval-val[nz-1][ip][ir];
	}
	break;
      
      case IS_DENSITY:
	/* density flux remains constant at the sheath that is dz N = -1/V dz V*/
	FORYX_BD{
	  tval = 0.;
	  val[nz][ip][ir] = tval + val[nz-1][ip][ir];
	}  
	break;
        
      case IS_POTENTIAL:
	/* No gradient in potential, sheeth works via parallel current */
	FORYX_BD  val[nz][ip][ir] = val[nz-1][ip][ir];
	break;
      
      case IS_TEMPERATURE:
	break;
      
      default:
	/* No gradient  */
	FORYX_BD		  val[nz][ip][ir] = val[nz-1][ip][ir];
	
      }
  
}

/*************************************************************************************/
void Cyto_DParallelSOL(int identity,double ***result,double ***result2,
                       double ***val,double ***f,double ***exp_n,double ***v,double ***dzv,
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
          /* Factor for derivative */
          
          lf_local  = 1./(12.* para->dz);
          lf_local2 = 1./(12.* para->dz*para->dz);        
          /* Potential is relative to sheath potential */

          FIRST = FALSE;
      }
    
    BUGREPORT;
    
    /* Communicate Ghostpoints */

    PH_UpdateZ(val,data,MPI_COMM_WORLD);

    BUGREPORT;
    if(offz == 1)   
    {
        BUGREPORT;
	if( bc_sheath ) Cyto_no_sheath_condition(identity, result, result2, val, f, exp_n, v, dzv, data, para);
	else     Cyto_sheath_condition(identity, result, result2, val, f, exp_n, v, dzv, data, para);

        BUGREPORT;  

        /* dz*/ 
        FORALL_BD_XY  result[iz][ip][ir] = 6.*lf_local*   (  val[iz+1][ip][ir] - val[iz-1][ip][ir] );
	

        /* dzz */
        FORALL_BD_XY result2[iz][ip][ir] = 12.*lf_local2* ( (val[iz+1][ip][ir] + val[iz-1][ip][ir]) 
							    - 2.* val[iz][ip][ir] );

        BUGREPORT;
    }
    else /* Two boundary points, higher accuracy, but not implemented in parallel, boundarz exchange is neither */
    {
        // Stefan
	printf("Never get here!!!\n");

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
  
			// Stefan
                        val[-1][ip][ir] = val[0][ip][ir];
                        val[-2][ip][ir] = val[1][ip][ir];
			/*
                        val[-1][ip][ir] =  -val[0][ip][ir];
                        val[-2][ip][ir] = -val[1][ip][ir];
			*/
                  }
                break;
            case  IS_IONVELOCITY:
                /* Ion speed is unity into the limiter */
                FORYX_BD
                    {
                        bdval=2.;
                        val[nz][ip][ir] =   bdval-val[nz-1][ip][ir];
                        val[nz+1][ip][ir] =  bdval-val[nz-2][ip][ir];
     

                        val[-1][ip][ir] =  val[0][ip][ir];
                        val[-2][ip][ir] = val[1][ip][ir];

                    }
                break;
            case IS_DENSITY:
                /* density flux remains constant at the sheath dz N = */
                FORYX_BD
                {
/*                    tval = (exp(-f[nz][ip][ir])-v[nz-1][ip][ir])/para->dz;
                    tval *=   -0.5/MAX(exp(-f[nz][ip][ir]),0.4)/(6.*lf_local);

                    
                  
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
            default:
                fprintf(stderr,"unidentified field object (UFO) in  boundary routine!\n");
        }

        
        /* dz */
        FORALL_BD_XY  
            result[iz][ip][ir]  =  lf_local*(-val[iz+2][ip][ir] + 8.*(val[iz+1][ip][ir] - val[iz-1][ip][ir]) + val[iz-2][ip][ir]);
        
        /* dzz */
        FORALL_BD_XY  
            result2[iz][ip][ir]  =  lf_local2 *(-val[iz+2][ip][ir] + 16.*(val[iz+1][ip][ir] + val[iz-1][ip][ir]) -30.* val[iz][ip][ir] - val[iz-2][ip][ir]);
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
          for(j=0;j<ny;j++) for(k=0;k<nx;k++) w[i][j][k]= f[i][j][k];
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
      /* If complete FD: Uncomment here */
      
        for(j=0;j<ny;j++) for(k=0;k<nx;k++)   w[i][j][k] =  vval[i][k]*vval[i][k]*(f[i][j+1][k] -2.*f[i][j][k]  +f[i][j-1][k]);
      

      for(j=0;j<ny;j++) for(k=0;k<nx;k++)  w[i][j][k] += hval[i][k]*hval[i][k]*(f[i][j][k+1] -2.*f[i][j][k]  +f[i][j][k-1]);
     

      /*if(ispolar)*/
      for(j=0;j<ny;j++) for(k=0;k<nx;k++)  w[i][j][k] += 0.5*edr[k]*hval[i][k]*(f[i][j][k+1] - f[i][j][k-1]);
      
  }
}






/* This routine calculates secondary parameters from primary ones */
void Cyto_CalcParameters(PARA *p,HDF_DS *d)
{
    double val;
    double lgamma = 2.;
    double Bgauss=0.;
    double   lamda_ei =0.;
    double n_cgs;
    double boltzmann_k_si = 1.3807e-23; /* J/K */
    double energy_per_ev = 1.1604e4; /* K */
    double m_proton = 1.6726e-27; /*kg*/
    int i;
    
        
    Util_CalcParameters(p,d);

    /* Calculate equation parameters */


    p->unit_time = p->rho_s/p->c_s;
    p->unit_length = p->rho_s;

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
    
    fprintf(stderr,"\t%-20s = %f \n","nu",p->nu);
    fprintf(stderr,"\t%-20s = %f \n","mue_w",p->mue_w);
    fprintf(stderr,"\t%-20s = %f \n","mue_n",p->mue_n);
    fprintf(stderr,"\t%-20s = %f \n","delta(nu_en)",p->delta);
    fprintf(stderr,"\t%-20s = %f \n","kappan",p->kappan);
    fprintf(stderr,"\t%-20s = %f \n","sigma(nu_in)",p->sigma);


    /* Change length scales */

    p->xmin = 0.;
    p->xmax = p->a/p->rho_s;
    
    p->ymin = 0.;
    p->ymax = 2.*M_PI;
    
    p->zmin = 0.;
    p->zmax = p->lpar/p->rho_s;
    
    p->dx = (p->xmax-p->xmin)/((double)d->nx) ;
    p->dkx = 2.0*M_PI/(p->xmax-p->xmin);
    
    p->dy = (p->ymax-p->ymin)/((double)d->ny) ;
    p->dky = 2.0*M_PI/(p->ymax-p->ymin);
     
    p->dz = (p->zmax-p->zmin)/((double)d->nz) ;
    p->dkz = 2.0*M_PI/(p->zmax-p->zmin);
    
    
    
    /* Make sure coodinates are consistent */
    
    for(i=0;i<3;i++)
        if(d->coordinate[i] != NULL)  free(d->coordinate[i]);
    for(i=0;i<3;i++)
        d->coordinate[i] = (double *)malloc(sizeof(double)*d->dims[i]);
    
    if(d->rank ==3)
    {
        
        for(i=0,val = p->zmin+p->dz*.5;i<d->nz;i++,val+=p->dz)  d->coordinate[0][i] = val;
        for(i=0,val = p->ymin+p->dy*.5;i<d->ny;i++,val+=p->dy)  d->coordinate[1][i] = val;
        for(i=0,val = p->xmin+p->dx*.5;i<d->nx;i++,val+=p->dx)  d->coordinate[2][i] = val;
    }
    else if (d->rank == 2)
    {
        for(i=0,val = p->ymin+p->dy*.5;i<d->ny;i++,val+=p->dy)  d->coordinate[0][i] = val;
        for(i=0,val = p->xmin+p->dx*.5;i<d->nx;i++,val+=p->dx)  d->coordinate[1][i] = val;
    }
    
    

    fprintf(stderr,"Domain:\n");
    fprintf(stderr,"\t [ %f: %f ] x [ %f : %f ] x [ %f : %f ] \n\t (%f,%f,%f) \n\t (%d,%d,%d)\n\n",
            p->xmin,p->xmax,p->ymin,p->ymax,p->zmin,p->zmax,
            p->dx,p->dy,p->dz,
            (int)d->nx,(int)d->ny,(int)d->nz);

}

/*
  Function maps spectral noise onto poloidal domain */

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



double cyto_blob_distribution_z(int iz, int nz, HDF_DS *data)
{
  double middle = ((double)data->dims[0]-1.0) / 2.0;
  int real_iz   = nz* data->grid_coords[0] + iz;

  double result = (abs(middle-real_iz)-0.5);

  if(result < 0) result=0.0;

  return result;

}
