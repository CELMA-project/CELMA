/* Solves drift wave dispersion in a cylinder with variable collision terms and background rotation */


#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <malloc.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#define FALSE	0
#define TRUE	1
#define OUTPUT_PHASE /* Write n-phi phase relationship, otherwise write raw components of solution */

/* Minimum density  */
#define NMIN 1.e-18
#define MINERROR 1.e-11 /* Minimal error */
#define MINOM 1.e-8/* Accuracy Omega minimal */
#define MAXOM 1.e-4 /* Do not refine range below */
#define ARGON 40.
#define PROTONMASS 1800.
#define NDIM 4

void  rQ(double r,double *Qr,double *Qi,double *kappar,double *kappai);
int func (double t, const double y[], double f[],
          void *params);
int jac (double t, const double y[], double *dfdy, 
     double dfdt[], void *params);



#define ADAPTIVE
#define NSAMPLE 31
#define LASTM 4
#define N 256
#ifndef MAX
#define MAX(x,y)                (((x)>(y))?(x):(y))
#endif
#define INTERP(a,b) glob_delta*a[(b)] + (1.-glob_delta)*a[(b)+1]

/*********************/

/* Global Variables */

double m,
    kp,                /* k parallel in 1/rho_s*/
    v_par,        /* Electron parallel drift in multiples of C_s */
    omr,          /* Real part of eigenvalue */
    omi;         /* Imaginary  part of eigenvalue */
   
int glob_index;
double glob_delta;


double density[N+10],
    dx_density[N+10],
    potential[N+10],
    velocity[N+10],
    S_p[N+10],
    nei[N+10],
    nen[N+10],
    nin[N+10];

double rcor[N+10];


int main(int argc,char **argv)
{
    extern double omr,omi;
    extern double m,v_par;

    double error,emin;
    double omcr=0.,omci=0.;
    double     omrs,omre,omis,omie,domr,domi,LOM=3.0;    
    double y[N+100][6], ys[N+100][6];
    double result[LASTM+1][N+100][6];
    double omega[LASTM+1][4]; /* contains omr, omi, found, lom */
    double lp,Te=0.0,Vp=0.,norm;
    double Qr,Qi;
    double DR,DI,ZR,ZI,NORM;
    double nfr,nfi,OR,OI,tomr,tomi,wstar;
    double kappa,mvp,P_r,P_i,om1;
    double max_nu_in;
    

    double kappar,kappai;
    double r,dr,par_start=0.,par_end=0.,dpar=0.,par;
    double A=1.e-12,B; /* Bounds of integration */
    int i,j,iter= 0;
    int k,l,im;
    int in,ind,seekn=0,scan_steps;
    char name[512];
    char lineout[512];
    char profile_file[512];
    char tmpc[2050];
    int FIRST = TRUE;
    int WriteEM=FALSE;
    int mode_file = 1;
    int FOUND = FALSE;
    
    /* Temporary variable to read profile file */
    double rt,nt,tt,tnei,tnin,ten;
    double te,nnorm,rho_s,wci;
    
    
    
    FILE *output;
    FILE *input;
    

    /* Stuff for the gnu solver */

    int status;
    
    const int Pdim=NDIM;

/* Define the stepping algorithm */
/*
  Step Type: gsl_odeiv_step_rk2 
  Embedded 2nd order Runge-Kutta with 3rd order error estimate. 
  
  Step Type: gsl_odeiv_step_rk4 
  4th order (classical) Runge-Kutta. 
  
  Step Type: gsl_odeiv_step_rkf45 
  Embedded 4th order Runge-Kutta-Fehlberg method with 5th order error estimate. This method is a good general-purpose integrator. 
  
  Step Type: gsl_odeiv_step_rkck 
  Embedded 4th order Runge-Kutta Cash-Karp method with 5th order error estimate. 
  
  Step Type: gsl_odeiv_step_rk8pd 
  Embedded 8th order Runge-Kutta Prince-Dormand method with 9th order error estimate. 
  
  Step Type: gsl_odeiv_step_rk2imp 
  Implicit 2nd order Runge-Kutta at Gaussian points 
  
  Step Type: gsl_odeiv_step_rk4imp 
  Implicit 4th order Runge-Kutta at Gaussian points 
  
  Step Type: gsl_odeiv_step_bsimp 
  Implicit Bulirsch-Stoer method of Bader and Deuflhard. This algorithm requires the Jacobian. 
  
  Step Type: gsl_odeiv_step_gear1 
  M=1 implicit Gear method 

  Step Type: gsl_odeiv_step_gear2 
  M=2 implicit Gear method 
*/
    const gsl_odeiv_step_type *T= gsl_odeiv_step_rk8pd;
    gsl_odeiv_step  *s= gsl_odeiv_step_alloc (T, Pdim);

    double dfdy[16];
#ifdef ADAPTIVE
   gsl_odeiv_control * c 
        = gsl_odeiv_control_y_new (1e-8, 0);
#endif
    gsl_odeiv_evolve * e 
        = gsl_odeiv_evolve_alloc (Pdim);
    gsl_odeiv_system   sys = {func, jac, Pdim, NULL};
    double t,t1,h=2.e-3;
    int mode,scan_kind;
    char yesno;
    double ln=1.;
    
  

    /******************* End variable declaration *********************************/


    /* Some values for variables */
    kp = 0.01;
    v_par = 35.;
    lp = M_PI/(2.*kp);

    


    /* Read user defined variables */


    printf("*************************************************************************************\n");
    printf("This programm solves the dispersion relation for Driftwaves in a cylinder.\n");
    printf("It uses a shooting method to find the Eigenvalue Omega\n");
    printf("and the Eigenfunction PHI**2/(n_0*n_0) for poloidal modes m=1 to m=5.\n");
    printf("It might also scan the eigenvakues given a parameter range.\n");
    printf("The Eigenfuction is choosen to have less than n zeros in the Intervall [0:B].\n\n");
    printf("*************************************************************************************\n");
    printf("The user must provide the following parameters:\n");
    printf("\tThe lenght L -- in multiples of rho_s --of the discharge, we asume that the wave-length is L_par = 4 L\n");
    printf("\tThe parallel electron drift velocity v_par in multiples of ion-sound speed\n");
    printf("\tThe collision frequencies of ions-electrons, and with neutrals, each given as \n");
    printf("\t\tmultiples of the ion-cyclotron frequency.\n");
    printf("*************************************************************************************\n");
    printf("On output the file sol.dat contains the mode structure normalized to one.\n");
    printf("The Eigenvalue and the number of radial zeros is printed on the screen.\n");
    printf("*************************************************************************************\n\n");
    
#ifdef ADAPTIVE
    fprintf(stderr,"Adaptive control method is '%s'\n",gsl_odeiv_control_name(c));
#endif

/********************************************************************/
    fflush(stdin);
    printf("Input:  Read profiles from file [0] or assume Gaussian [1]? \n");
    scanf("%d",&mode_file);
    fflush(stdin);


#define MODE_READFROMFILE 0
#define MODE_ASSUMEGAUSS 1    
    if(mode_file != MODE_READFROMFILE) mode_file = MODE_ASSUMEGAUSS;
    
    if(mode_file == MODE_READFROMFILE)
    {
        printf("Input: Name of profile file?\n");
        scanf("%s",&profile_file);
        fflush(stdin);
        printf("Read profiles from file >%s<\n",profile_file);
    }
    else
        printf("Gaussian Profile is used.\n");




/***********************************************************************************/


#define MODE_EIGENFUNCTION 0
#define MODE_PARAMETERSCAN 1

    fflush(stdin);
    printf("Input:  Choose Eigenfunction calculation [%d] or parameter scan [%d]? \n",MODE_EIGENFUNCTION,MODE_PARAMETERSCAN);
    scanf("%d",&mode);
    fflush(stdin);

  
    if(mode != MODE_EIGENFUNCTION) mode = MODE_PARAMETERSCAN;

    if(mode == MODE_PARAMETERSCAN)
    {
        printf("Input: Write Eigenmodes to file? [y] or [n]?\n");
        scanf("%1s",&yesno);
        fflush(stdin);
        if (yesno == 'Y' || yesno == 'y') 
        {
            WriteEM = TRUE;
            printf("Write also Eigenmodes.\n");
        }
        else 
        {
            WriteEM = FALSE;
            printf("Not writing Eigenmodes.\n");
        }
        
    }


/******************************************** Read in common stuff ********************************/

    printf("Input:  Parallel Length in Units of rho_s:\n");
    if( scanf("%lf",&lp) != 1) lp = 100.;
    fflush(stdin);
    kp=M_PI/(2.*lp);
    
    printf("Input: parallel electron drift velocity v_par (c_s):\n");
    scanf("%lf",&v_par); 
    fflush(stdin);

    printf("Input: Radial order of solution:\n");
    scanf("%d",&seekn);
    fflush(stdin);


/********************************** Read from user ******************************************************/

    if(mode_file != MODE_READFROMFILE)
        {
            printf("Input:  Gradient  Length (rho_s):\n");
            scanf("%lf",&ln);
    
            printf("Input:  Rotation frequency (Omega_ci):\n");
            scanf("%lf",&Vp);

            printf("Input: Radius of machine (rho_s):\n");
            if( scanf("%lf",&B) != 1) B=3.5;
            fflush(stdin);
   
            printf("Input: collision frequencies: nu_ei, nu_en, nu_in (Omega_ci):\n");
            scanf("%lf",&nei[0]);    fflush(stdin);
            scanf("%lf",&nen[0]);      fflush(stdin);
            scanf("%lf",&nin[0]);    fflush(stdin);
     
            printf("Starting calculation with following values:\n");
            printf("\t kp = %f \t lp = %f\t Vp 0 %f\n",kp,lp,Vp);
            printf("\t v_par = %f\n",v_par);   
            printf("\t nu_ei = %f\n",nei[0]);   
            printf("\t nu_en = %f\n",nen[0]);   
            printf("\t nu_in = %f\n",nin[0]);   
            printf("\t n = %d\n",seekn);  
            printf("\t Radial interval [%g:%g]\n",A,B);
        }
    
/************************Parameterscan *******************************************************/


    if(mode == MODE_PARAMETERSCAN)
    {
        printf("Scanning of some parameter is choosen.\n Please select parameter to vary:\n");   
        scan_kind = 7;

        /* We can't scan on collisions and gradient length, when we read these from a file */

        if(mode_file != MODE_READFROMFILE)
        {
            printf("\t 0: L parallel\n");   
            printf("\t 1: v parallel\n");
            printf("\t 2: nu_ei\n");
            printf("\t 3: nu_en\n");
            printf("\t 4: nu_in\n");
            printf("\t 5: ln\n");
            printf("\t 6: Vp\n");    
            
            while (scan_kind >6 || scan_kind < 0 ) 
            { 
                printf("Please select  valid option!\n");
                scanf("%d",&scan_kind);
                fflush(stdin);
           }
            
        }
        else
        {
            printf("\t 0: L parallel\n");   
            printf("\t 1: v parallel\n");
            printf("\t 6: V_p\n");
            while (scan_kind !=1 && scan_kind != 0 && scan_kind != 6 ) 
            { 
                printf("Please select  valid option!\n");
                scanf("%d",&scan_kind);
                fflush(stdin);
            }
           
        } 
       
        /*********************** read in scan range *******************************/
        switch(scan_kind)
        {
            case 0:
                printf("Scanning L parallel from %g to ? \n",lp);
                par_start = lp;
                break;
            case 1:
                printf("Scanning v parallel from %g to ? \n",v_par);
                par_start = v_par;
                break;     
            case 2:
                printf("Scanning nu_ei from %g to ? \n",nei[0]);
                par_start = nei[0];
                break;     
            case 3:
                printf("Scanning nu_en from %g to ? \n",nen[0]);
                par_start = nen[0];
                break;     
            case 4:
                printf("Scanning nu_in from %g to ? \n",nin[0]);
                par_start = nin[0];
                break;     
            case 5:
                printf("Scanning ln from %g to ? \n",ln);
                par_start = ln;
                break;     
            case 6:
                printf("Scanning Vp from %g to ? \n",Vp);
                par_start = Vp;
                break;

        }
      
      
        while( scanf("%lf",&par_end) != 1);
        fflush(stdin);
        
        printf("Input: Number of steps ? \n");
        scanf("%d",&scan_steps);
        fflush(stdin);

        printf("End scanning at %g after %d steps\n",par_end,scan_steps);
        dpar =(par_end-par_start)/(double)(scan_steps-1);
        fflush(stdin);
    }

    /*********************** density profile and derivative ***********************************/


    if(mode_file == MODE_READFROMFILE)
    {
        /* Read Profile from file */
        if((input = fopen(profile_file,"r")) == NULL)
        {
            fprintf(stderr,"Error, couldn't open  profile file %s.\n;",profile_file);
            exit(1);
        }
        /* Read normalisation values rhos and wci*/
        
        rho_s = wci = 0.;
        
        while ( fgets(tmpc,2047,input) != NULL)
        {
            if(strncmp(tmpc,"%#",2)== 0){
                sscanf(tmpc+2,"%lf %lf %lf",&wci,&rho_s,&Te);
                fprintf(stderr,"Read rho_s= <%f>, w_ci= <%f>, and T_e= <%f>.\n",rho_s,wci,Te);
            }
            
        }
        

        /* Check if we have normalization in code */

        if (0. == (rho_s*wci*Te) ) 
        { 
            fprintf(stderr,"Error, could not read rho_s or w_ci or T_e from %s.\n;",profile_file);
            exit(1);
        }
            
        rewind(input);


        /* Check how many entries there are in the file */
        
        i=0;
        while ( fgets(tmpc,2047,input) != NULL)
            if(strchr(tmpc,'#') == NULL && strchr(tmpc,'%') == NULL&& strchr(tmpc,'/')  == NULL) i++;
        
        if(i != N) 
        { 
            fprintf(stderr,"Error, file %s contains only %d positions.\n;",profile_file,i);
            exit(1);
        }

        rewind(input);
        
        
        /******************************* Read the data *******************************************/
        
        
        i=0; 
        while ( fgets(tmpc,2047,input) != NULL)
            if(strchr(tmpc,'#') == NULL && strchr(tmpc,'%') == NULL&& strchr(tmpc,'/')  == NULL)
            {
                sscanf(tmpc,"%lf %lf %lf %lf %lf %lf",&rt,&nt,&tt,&tnei,&tnin,&ten);
                /*fprintf(stderr,"%f %f %f %f %f %f\n",rt,nt,tt,tnei,tnin,ten);*/

               
                /************* Normalize to values in the center *******************/

                if(i == 0)
                {

                    if(rt <  A) 
                    {
                        rt = A*rho_s;
                        
                    }
                    else
                    {
                        
                        A = rt/rho_s;
                    }
                    
                    nnorm = MAX(nt,NMIN);
                }
                rcor[i] = rt/rho_s;
                density[i] = nt/nnorm;
                potential[i] = tt/Te;
                nei[i] = tnei/wci;
                nin[i] = tnin/wci;
                nen[i] = ten/wci;
                i++;
            }
        fclose(input);


        /* Ready reading file */
    
        B = rt/rho_s;
        dr = (B-A)/(double)N;


    }
    else
    {
        /********************** Set up Gaussian profile *************************************/
        /* The gridpoints */

        dr = (B-A)/(double)N;
        for(i=0;i<N;i++) rcor[i] = A + ((double)i+0.5)*dr;
        for(i=0;i<N;i++)  density[i]= exp(-(rcor[i]*rcor[i])/(ln*ln));
   
        /* Potential is never used */
        for(i=0;i<N;i++)  potential[i]= 0.;
        for(i=0;i<N;i++)  velocity[i]= Vp;

        /* Shear is zero for the solid state rotation */
        for(i=0;i<N;i++) S_p[i]= 0.;

        for(i=1;i<N;i++) nei[i] = nei[0]* density[i];
        for(i=1;i<N;i++) nin[i] = nin[0];
        for(i=1;i<N;i++) nen[i] = nen[0];

    }
    
    for(i=0;i<N;i++) density[i] = MAX(density[i],NMIN);
    

    /* Calculate dx_densisty */

    for(i=1;i<N-1;i++) dx_density[i] = 0.5/dr*(log(density[i+1])-log(density[i-1]));
    dx_density[0] = 0 ;
    dx_density[N-1]= 2.*dx_density[N-2]-dx_density[N-3];

    dr = (B-A)/(double)N;

    /* Calculate V_0 = 1/r*dx_potential */
    /* Only if we did read potential from file! */

    if(mode_file == MODE_READFROMFILE)
    {
        
        for(i=1;i<N-1;i++) velocity[i] =  0.5/dr/rcor[i]*(potential[i+1]-potential[i-1]);
        velocity[0] = velocity[1] ;
        velocity[N-1]= 2.*velocity[N-2]-velocity[N-3];
        

        /* Calculate S_p */

        for(i=1;i<N-1;i++) 
            S_p[i]  = 3.0*0.5/dr*(velocity[i+1]-velocity[i-1]);
        
        for(i=1;i<N-1;i++) 
               S_p[i]  += rcor[i] *(velocity[i+1]-2.*velocity[i]+velocity[i-1])/dr/dr;
        
        S_p[0] = 2.*S_p[1]-S_p[2] ;
        S_p[N-1]= 2.*S_p[N-2]-S_p[N-3];
    
        /* for(i=1;i<N-1;i++) fprintf(stderr,"r %f dr %f S_p %f\n",rcor[i],dr,S_p[i]);*/
        
    
    }
    

    
    /* Write data to file */
    
    output = fopen("profile","w");
    fprintf(output,"%% A %f; B %f, dr %f Te=%f\n",A,B,dr,Te);
    fprintf(output,"%% i \t r     \t \t   n(r) \t dx_n \t nei \t nin \t nen \t phi \t Vp \t Sp\n");

    for(i=0;i<N;i++) fprintf(output,"%d \t %.8g \t %.8g \t %.8g \t %.8g \t %.8g \t %.8g \t %.8g \t %.8g \t %.8g\n",
                             i,rcor[i],density[i],dx_density[i],
                             nei[i],nin[i],nen[i],potential[i],
                             velocity[i],S_p[i]);
    fclose(output);
    
    


    /* Ready making profiles */

    /*********************************************************************************************/

    /* Write header for output file */

    sprintf(name,"eigenvalue_n=%d.dat",seekn);
    output = fopen(name,"w");
    fprintf(output,"%%\tL_par\t= %g\n",M_PI/(2.*kp));
    fprintf(output,"%%\tV_par\t= %g\n",v_par);
    fprintf(output,"%%\tnu_ei\t= %g\n",nei[0]); 
    fprintf(output,"%%\tnu_en\t= %g\n",nen[0]); 
    fprintf(output,"%%\tnu_in\t= %g\n %%\n",nin[0]);
    fprintf(output,"%%\tln\t= %g\n%%\n",ln);   
    fprintf(output,"%%\tVp\t= %g\n%%\n",Vp);   
    fprintf(output,"%%\tRadial Range [%g:%g]\t= %g\n%%\n",A,B,dr);   
    fprintf(output,"%% m n omr omi eps\n");
    fclose(output);

    /* Initialize Eigenvalues */


    for(im=1;im<=LASTM;im++) omega[im-1][0]= 4.5;
    for(im=1;im<=LASTM;im++) omega[im-1][1] = 0.;
    for(im=1;im<=LASTM;im++) omega[im-1][2] = FALSE;
    for(im=1;im<=LASTM;im++) omega[im-1][3] = 9.0;  /*LOM */ 


    for(im=1;im<=LASTM;im++) result[im-1][0][5] = -1;
    

    /* Loop over parameters */

    if(mode != MODE_PARAMETERSCAN) {
        par_end = par_start;   
        dpar = 100.;
        
    }
    
    

    for(par=par_start;dpar*par <= dpar*par_end;par+=dpar)
    {
        if(mode != MODE_PARAMETERSCAN) par_end = par_start;    
        else if (mode == MODE_PARAMETERSCAN)
        {
            /* Set up parameter scan */
            switch(scan_kind)
            {
                case 0:
                    lp=par;
                    kp=M_PI/(2.*lp);
                    break;
                case 1:
                    v_par=par;
                    break;     
                case 2:
                    if(mode_file == MODE_READFROMFILE)
                    {
                        fprintf(stderr,
                                "Error: You can't change collisions if the profile is read from a file\n");
                        exit(-1);
                    }        
                    for(i=0;i<N;i++) nei[i]=par;
                    break;     
                case 3:
                    if(mode_file == MODE_READFROMFILE)
                    {
                        fprintf(stderr,
                                "Error: You can't change collisions if the profile is read from a file\n");
                        exit(-1);
                    }        
                    for(i=0;i<N;i++) nen[i]=par;
                    break;     
                case 4:
                    if(mode_file == MODE_READFROMFILE)
                    {
                        fprintf(stderr,
                                "Error: You can't change collisions if the profile is read from a file\n");
                        exit(-1);
                    }        
                    for(i=0;i<N;i++) nin[i]=par;
                    break;     
                case 5:
                    if(mode_file == MODE_READFROMFILE)
                    {
                        fprintf(stderr,"Error: You can't change the gradient if the profile is read from a file\n");
                        exit(-1);
                    }        
                    ln=par;
                    for(i=0;i<N;i++)  density[i]= exp(-(rcor[i]*rcor[i])/(ln*ln));
                    for(i=0;i<N;i++) dx_density[i]= -2.*rcor[i]/(ln*ln);
                    break;
                case 6:
                    if(mode_file == MODE_READFROMFILE)
                    {
                        fprintf(stderr,"Error: You can't change plasma rotation if the profile is read from a file\n");
                        exit(-1);
                    }        
                    Vp=par;
                    for(i=0;i<N;i++)  velocity[i]= Vp;
                    break;     
            }
        }
        /* Ready setting up parameter scan */



        /******************* So, do we finally start to calculate something here ***************************/
        
        /* For all the m-modes we want to calculate */

        for(im=1;im<=LASTM;im++)
        {
            m = (double)im;

            fprintf(stderr,"m= %d, start searching....\t",im);
       
            /* Set up search area */
            
            if(omega[im-1][2] == FALSE) 
                {
                    /* Make a civilized guess for omega from the simple slab case */
                    omcr = omega[im-1][0]= .66;
                    omci = omega[im-1][1] = .5;
                    omega[im-1][2] = FALSE;
                    LOM = omega[im-1][3] = 4.*omcr;
                    fprintf(stderr,"fresh\n");
                }
            
            else
            {
                omcr = omega[im-1][0];
                omci =  omega[im-1][1];
                omega[im-1][2] = FALSE;
                LOM = omega[im-1][3]=1.25;
                fprintf(stderr,"from %g,%g width %g\n",omcr,omci,LOM);
            }
            
            emin =1.e100;
            h = 1.e-3;
            

            /* initialise dydt_in */
            GSL_ODEIV_JA_EVAL(&sys, rcor[0], y[0], dfdy, ys[0]);
            
            max_nu_in = 0.;
            for(i=0;i<N;i++)  max_nu_in = MAX(nin[i],max_nu_in);
            
            
            /* Determine the search area in real and imaginary part  */
            while( ((emin > MINERROR)  && (LOM > MINOM)) || LOM > MAXOM )
            {   
                omrs = MAX(omcr/(2.*(double)(NSAMPLE)),omcr-LOM*0.5);
                omre = omcr+LOM*0.5;
            
                omis = MAX(-max_nu_in,omci-LOM*0.5);
                omie = omci+LOM*0.5;

                FOUND = FALSE;
                omega[im-1][2] = FALSE;
                
                domr =  (omre-omrs)/(double)(NSAMPLE);
                domi =  (omie-omis)/(double)(NSAMPLE);   

                for(omr=omre;omr>=omrs;omr-=domr)
                {
                    for(omi=omie;omi>=omis;omi-=domi) /* search for mode with largest growthrate */
                    {
    
                        /* Initial values */

                        /* As modes with higher m start further out we say 
                           they are zero on  the first m gridpoints, otherwise the difficulty arises that they grow
                           above all limits 
                        */

                        for(i=0;i<3;i++)
                        {
                            t = rcor[i];
                            /* Values */
                            y[i][0] = pow(t,im);
                            y[i][1] = 0.0;
                        
                            y[i][2] = m*y[i][0];
                            y[i][3] = 0.0;
                        
                            /* Derivatives */
                            
                            rQ(t,&Qr,&Qi,&kappar,&kappai);

                            ys[i][0] = 1./t*y[i][2];
                            ys[i][1] = 1./t*y[i][3];

                            ys[i][2] = m*m/t* y[i][0] - Qr*y[i][0]+Qi*y[i][1] + kappar*y[i][2] - kappai*y[i][3] ;
                            ys[i][3] = m*m/t* y[i][1] - Qr*y[i][1]-Qi*y[i][0] + kappar*y[i][3] + kappai*y[i][2];
                            
                        }

                        /* initialise dydt_in */
#ifndef ADAPTIVE
                        GSL_ODEIV_JA_EVAL(&sys, rcor[0], y[0], dfdy, ys[0]);
#endif

                        in = 0;
                        ind = 0;
                        
                        for(i=0;i<N-1;i++) 
                        {
                            status= GSL_SUCCESS;
                            glob_index = i;
                            for (j=0;j<NDIM;j++)y[i+1][j]=y[i][j];                        
                            
#ifdef ADAPTIVE
                            t = rcor[i];
                            t1 = rcor[i+1];
                            h=1.e-4;

                            while (t < t1 && h > 1.e-7  && (status == GSL_SUCCESS ))
                            {
                                glob_delta = (t1-t)/(rcor[i+1]-rcor[i]);
                                status = gsl_odeiv_evolve_apply (e, c, s, 
                                                                 &sys,
                                                                 &t, t1,
                                                                 &h, y[i+1]);
                            }

                            /* Count number of zeros in realpart */
                            if(y[i+1][0]*y[i][0] < 0.0 ) in++;
                            if(in > seekn) {
                                status = !GSL_SUCCESS;
                                /* fprintf(stderr,"Mode at (%f,%f) has to many zeros. \n",omr,omi);*/
                            }
                            
                                
                            /* Count number of zeros in derivative of part */

                            if(y[i+1][2]*y[i][2] < 0.0 ) ind++;
                            if(ind > (seekn+1) ){
                                status = !GSL_SUCCESS;
                                /* fprintf(stderr,"Mode at (%f,%f) has to many zeros in derivative. \n",omr,omi);*/
                            }

                            /* Solution should be falling in amplitude at outermost 30 % of domain */
                               if(t>0.7*B)
                                if( (2.*dx_density[i]*density[i]*(y[i][0]*y[i][0]+y[i][1]*y[i][1])
                                                 +density[i]*density[i]*2./rcor[i]*(y[i][0]*y[i][2]+y[i][1]*y[i][3])) > 0. ) 
                                {
                                    status = !GSL_SUCCESS;
                                    /*fprintf(stderr,"Mode at (%f,%f) isn't decaying. \n",omr,omi);*/
                                }

                         /* Solution should be growing in amplitude at innermostmost 5 % of domain */
                               if(t<0.05*B)
                                if( (2.*dx_density[i]*density[i]*(y[i][0]*y[i][0]+y[i][1]*y[i][1])
                                                 +density[i]*density[i]*2./rcor[i]*(y[i][0]*y[i][2]+y[i][1]*y[i][3])) < 0. ) status = !GSL_SUCCESS;
                            


                            if(status  != GSL_SUCCESS) 
                            {
                                /*fprintf(stderr,"Integrated to r = %12.09g, Solution will be kept equal density for rest of domain.\n",t);*/
                                for(;i<N;i++) 
                                    for (j=0;j<NDIM;j++) 
                                    {
                                        glob_index = i; 
                                        y[i][j]= 1.e10*density[i];
                                    } 
                            }
#else              
                            status = gsl_odeiv_step_apply (s, rcor[i], dr, 
                                                           y[i+1], ys[i+2], ys[i], ys[i+1], &sys);
                            if (status != GSL_SUCCESS) break;
#endif     
                        }
                        
                        
                        /* Determine Error */
                        
                        /* Value shall be small at r+B */
                        
                        /* Contains absolute value of solution */
                        if (status != GSL_SUCCESS) 
                        {
                            error = 1000;
                        }
                        else
                        {
                            /* Count number of zeros */
                            

                            in = 0;
                            for(i=1;i<N-1;i++)if(y[i][0]*y[i+1][0] < 0.0 ) in++;

                            if(in == seekn)
                            {
                                for(i=0;i<N;i++)  ys[i][0] = density[i]*density[i]*(y[i][0]*y[i][0]+y[i][1]*y[i][1]);
                                norm = 0.;
                                for(i=0;i<N;i++) norm = MAX(norm ,ys[i][0]);
                                error = sqrt(ys[N-1][0])/norm;
                            }
                            else
                            {
                                error = 99.;
                                fprintf(stderr,"Mode zero count failed.\n");    
                            }
                                
                    
                        
                        /*
                          sprintf(lineout,"Error  %g found at (%.12g,%.12g). Mode (%d, %d) .",error,omcr,omci,(int)m,in);
                          printf("%- 120.90s\r",lineout);
                        */
                            
                        }
                        
                    if(error < emin)
                        {
                            emin = error;
                            omcr = omr;
                            omci = omi;
                            /* Primary output is fluctuation amplitude */
                            for(i=0;i<N;i++) result[im-1][i][0]= sqrt(ys[i][0]/norm);
                            /* Real and imaginary part of phi */
                            for(i=0;i<N;i++) result[im-1][i][1]= y[i][0];
                            for(i=0;i<N;i++) result[im-1][i][2]= y[i][1];



  
                          
#ifdef OUTPUT_PHASE
                            /* Phase relation n to phi */
                            for(i=0;i<N;i++) 
                                {
                                     kappa = -dx_density[i];
                                     r = rcor[i];
                                     wstar = m*kappa/r;  
                                     mvp =  m*velocity[i];
                                     tomr = omr - mvp;
                                     tomi = omi;
    
                                     
                                     /* calculate P */
                                     
                                     ZR = nei[i]+nen[i];
                                     ZI =  mvp;
                                     NORM = 1./(ZR*ZR+ZI*ZI);
                                     
                                     /* Ion Mass is at the moment fixed to Argon */
                                     
                                     P_r =   ( ARGON * PROTONMASS)*kp*kp*NORM*ZR;
                                     P_i =  -( ARGON * PROTONMASS)*kp*kp*NORM*ZI;    
                                     
                                     /* The n phi relationship */

                                     DR = wstar - P_i;
                                     DI =  P_r ;
                                     
                                     ZR = tomr - om1 - P_i;
                                     ZI =  tomi           + P_r;

                                     NORM = 1./(ZR*ZR+ZI*ZI);

                                     nfr =  NORM*(DR*ZR + DI*ZI);
                                     nfi =  NORM*(DI*ZR - DR*ZI);

                                     result[im-1][i][3]= nfr; 
                                     result[im-1][i][4]= nfi;
                                }
#else
                          /* Real and imaginary part of 1/r dr phi */
                            
                            for(i=0;i<N;i++) result[im-1][i][3]= y[i][2];
                            for(i=0;i<N;i++) result[im-1][i][4]= y[i][3];

#endif
                            /* Derivative of amplitude */

                            for(i=0;i<N;i++) result[im-1][i][3] = 2.*dx_density[i]*density[i]*(y[i][0]*y[i][0]+y[i][1]*y[i][1])
                                                 +density[i]*density[i]*2./rcor[i]*(y[i][0]*y[i][2]+y[i][1]*y[i][3]);


                            result[im-1][0][5] = result[im-1][2][5];
                            
                            /* And eigenvalue */

                            omega[im-1][0] = omr;
                            omega[im-1][1]= omi;
                            if(emin  < 1000)  omega[im-1][2]= TRUE;
                            omega[im-1][3]= 2.1*sqrt(domr*domr+domi*domi);
                            fprintf(stderr,"Minimal error  %10.07g found at (%10.07g,%10.07g). Mode (%.2d,%.2d).\r"
                                    ,emin,omcr,omci,(int)m,seekn);
                            if(emin  < 1000) FOUND = TRUE;
                            
                        }

                    }
                    
                }
                
                if(FOUND == TRUE) /* We found a minimal error and refine there */
                 {
                     LOM = omega[im-1][3];
                 }
                else  /* No minimal error was found, rescan with other grid  */
                {
                     omega[im-1][3] = LOM /= M_PI;
                     if (0. == omcr)
                     {
                         omcr=LOM/2.;
                         omci=LOM/2;
                     }
                 }
                
                domr = domi =  LOM/(double)(NSAMPLE-1);
                FOUND = FALSE;
                

                fprintf(stderr,"New search interval size: %10.07g at Omega= (%10.07g,%10.07g). Error %10.07g. Mode (%.2d,%.2d).\n",LOM,omcr,omci,emin,im,seekn);
            }
                
               fprintf(stderr,"\n\t Eigenvalue for m = %2.2d n = %2.2d: Omega= (%10.8g, %10.8g) eps %10.8g\n\n",(int)im,seekn,omcr,omci,emin);


            if(mode != MODE_PARAMETERSCAN)
            {
                sprintf(name,"eigenvalue_n=%d.dat",seekn);
                output = fopen(name,"a");
                fprintf(output,"%d %12.10g  %12.10g %12.10g\t ",(int)im,omcr,omci,emin);
                if(LASTM == im)  fprintf(output,"\n");       
                fclose(output);
            }
        
    }


/* Print OUTPUT */

        if(mode == MODE_EIGENFUNCTION || WriteEM)
        {
            sprintf(name,"eigenfunction_n=%d.dat",seekn); 
            
            if(FIRST)
            {
                
                output = fopen(name,"w");
                fprintf(output,"%%\tL_par\t= %g\n",M_PI/(2.*kp));
                fprintf(output,"%%\tV_par\t= %g\n",v_par);
                fprintf(output,"%%\tnu_ei\t= %g\n",nei[0]); 
                fprintf(output,"%%\tnu_en\t= %g\n",nen[0]); 
                fprintf(output,"%%\tnu_in\t= %g\n",nin[0]); 
                fprintf(output,"%%\tln\t= %g\n",ln); 
                fprintf(output,"%%\tradial order\t= %d\n\n",seekn); 
                if(mode == MODE_PARAMETERSCAN )
                {
                    switch(scan_kind)
                    {
                        case 0:
                            fprintf(output,"%% Lpar \t ");
                            break;
                        case 1:
                            fprintf(output,"%% Vpar \t ");
                            break;     
                        case 2:
                            fprintf(output,"%% nu_ei \t ");
                            break;     
                        case 3:
                            fprintf(output,"%% nu_en \t ");
                            break;     
                        case 4:
                            fprintf(output,"%% nu_in \t ");
                            break;     
                        case 5:
                            fprintf(output,"%% ln  \t ");
                            break;     
                        case 6:
                            fprintf(output,"%% Vp  \t ");
                            break;            
                    }
                    
                }
                else 
                    fprintf(output,"%% i \t ");
                
                fprintf(output,"r \t n \t nu_ei \t nu_en \t nu_in ");
                for(k=0;k<LASTM;k++) fprintf(output,"%d:amp\t y1R\t y1I\t PhaseR\t PhaseI\t ",k+1);
                fprintf(output,"\n");

                fprintf(output,"#%% 1 \t 2 \t 3 \t 4  \t 5 \t 6 \t 7 \t");
                for(k=0;k<LASTM;k++) fprintf(output," %d \t %d \t %d\t %d\t %d\t ", (7+k*5)+1,(7+k*5)+2,(7+k*5)+3,  (7+k*5)+4, (7+k*5)+5);
                fprintf(output,"\n");



            }
            else
                output = fopen(name,"a");

            if(mode == MODE_PARAMETERSCAN )
            {
                for(i=0;i<N-1;i++) 
                {
                    fprintf(output,"%10.08g \t%10.08g \t %10.08g \t%10.08g \t%10.08g \t%10.08g ",par,rcor[i],density[i],nei[i],nen[i],nin[i]);
                    for(k=0;k<LASTM;k++)  fprintf(output,"\t%10.08g \t%10.08g \t%10.08g \t%10.08g \t%10.08g ",
                                                  result[k][i][0],result[k][i][1],result[k][i][2],
                                                  result[k][i][3],result[k][i][4]);
                    fprintf(output,"\n");
                }
            }
            else
            {
                
                for(i=0;i<N-1;i++) 
                {
                    fprintf(output,"%5d \t%10.08g \t%10.08g \t%10.08g \t%10.08g \t%10.08g ",i,rcor[i],density[i],nei[i],nen[i],nin[i]);
                    for(k=0;k<LASTM;k++)  fprintf(output,"\t%10.08g \t%10.08g \t%10.08g \t%10.08g \t%10.08g ",
                                                  result[k][i][0],result[k][i][1],result[k][i][2],
                                                  result[k][i][3],result[k][i][4]);
                    fprintf(output,"\n");
                }
            }
            fprintf(output,"\n");
            fclose(output);
        }




        if(mode == MODE_PARAMETERSCAN )
        {
            switch(scan_kind)
            {
                case 0:
                    sprintf(name,"Lpar_scan_eigenvalues_n=%d.dat",seekn);
                    break;
                case 1:
                    sprintf(name,"Vpar_scan_eigenvalues_n=%d.dat",seekn);
                    break;     
                case 2:
                    sprintf(name,"nu_ei_scan_eigenvalues_n=%d.dat",seekn); 
                    break;     
                case 3:
                    sprintf(name,"nu_en_scan_eigenvalues_n=%d.dat",seekn);
                    break;     
                case 4:
                    sprintf(name,"nu_in_scan_eigenvalues_n=%d.dat",seekn);
                    break;     
                case 5:
                    sprintf(name,"ln_scan_eigenvalues_n=%d.dat",seekn);
                    break;     
                case 6:
                    sprintf(name,"Vp_scan_eigenvalues_n=%d.dat",seekn);
                    break;            
            }
            


            if(FIRST)
            {
                output = fopen(name,"w");
                fprintf(output,"%%\t L_par\t= %g\n",M_PI/(2.*kp));
                fprintf(output,"%%\t V_par\t= %g\n",v_par);
                fprintf(output,"%%\t nu_ei\t= %g\n",nei[0]); 
                fprintf(output,"%%\t nu_en\t= %g\n",nen[0]); 
                fprintf(output,"%%\t nu_in\t= %g\n",nin[0]); 
                fprintf(output,"%%\t ln\t= %g\n",ln); 
                fprintf(output,"%%\t VP\t= %g\n",Vp); 
                fprintf(output,"%%\t radial order\t= %d\n",seekn);

                switch(scan_kind)
                {
                    case 0:
                        fprintf(output,"%% Lpar\n ");
                        break;
                    case 1:
                        fprintf(output,"%% Vpar \n");
                        break;     
                    case 2:
                        fprintf(output,"%% nu_ei\t ");
                        break;     
                    case 3:
                        fprintf(output,"%% nu_en\t ");
                        break;     
                    case 4:
                        fprintf(output,"%% nu_in\t ");
                        break;     
                    case 5:
                        fprintf(output,"%% ln  \t");
                        break;     
                    case 6:
                        fprintf(output,"%% Vp  \t");
                        break;            
                }
                
                for(j=0;j<LASTM;j++) fprintf(output,"om(%d) gam(%d) \t",j+1,j+1); 
                fprintf(output,"\n");
                fclose(output);
            }

            output = fopen(name,"a");
            fprintf(output,"%g ",par);
            for(j=0;j<LASTM;j++) fprintf(output,"%g %g ",omega[j][0],omega[j][1]); 
            fprintf(output,"\n");
            fclose(output);
            printf("\nPrinted EW for Parameter %g\n",par);
            FIRST = FALSE;
        }
        FIRST = FALSE;
    }
    
    gsl_odeiv_evolve_free(e);
#ifdef ADAPTIVE
    gsl_odeiv_control_free(c);
#endif
    gsl_odeiv_step_free(s);
    return 0;
    
}
/******************/

void  rQ(double r,double *Qr,double *Qi,double *kappar,double *kappai)
{

    extern double m, nin[N+10],nei[N+10],nen[N+10],density[N+10],dx_density[N+10],rcor[N+10],
        velocity[N+10],S_p[N+10],kp,v_par, omr, omi;
    int i=glob_index;
    double wstar,kappa;
    double DR,DI,ZR,ZI,OR,OI,NORM;
    double nfr,nfi,rdi,rdr;
    
    double P_i,P_r,mvp,om1;    
    double tomr, tomi;
    


    kappa = -INTERP(dx_density,i);
    wstar = m*kappa/r;
    mvp =  m*(INTERP(velocity,i));
    
    tomr = omr - mvp;
    tomi = omi;
    

    /* calculate P */

    ZR = (INTERP(nei,i)+INTERP(nen,i));
    ZI =  mvp;
    NORM = 1./(ZR*ZR+ZI*ZI);

    /* Ion Mass is at the moment fixed to Argon */

    DR =  (ARGON*PROTONMASS)*kp*kp;

    P_r =   NORM*DR*ZR;
    P_i =  -NORM*DR*ZI; 

    /*Calculate W_1*/

    om1 = kp*v_par;
    
    /* Calculate terms in bracket */

    /* First the n phi relationship */

    DR = wstar - P_i;
    DI =  P_r ;
    
    ZR = tomr - om1 - P_i;
    ZI =  tomi           + P_r;

    NORM = 1./(ZR*ZR+ZI*ZI);

    nfr =  NORM*(DR*ZR + DI*ZI);
    nfi =  NORM*(DI*ZR - DR*ZI);

    /********************/


    OR = wstar - (tomr*nfr-tomi*nfi);
    OI =       -(tomi*nfr+tomr*nfi);

    /* Add rotation related part imSp */

    OR += m*INTERP(S_p,i)/r;

    /*Thats it, now multiply once with 1/(omega + i nu_in)  and r*/
  

    ZR = tomr ;
    ZI = tomi + INTERP(nin,i);

    NORM = 1./(ZR*ZR+ZI*ZI);

    *Qr =  r*NORM*(OR*ZR+OI*ZI) ;
    *Qi =  r*NORM*(OI*ZR-OR*ZI);


    /* Kappa Modified, reuse (wtilde + i nu_in)       */

    /* Calculate RD */

     DI =  r*INTERP(nin,i)*INTERP(velocity,i);
     
     rdr =  NORM*(DI*ZI) ;
     rdi =  NORM*(DI*ZR);

    *kappar = kappa - (rdr*nfr -rdi*nfi) ;
    *kappai = -(rdr*nfi+rdi*nfi);

}


/****************************************************/
/* GSL Part */


int func (double t, const double y[], double f[],
      void *params)
{
    extern double m;
    double kappar,kappai;
    double Qr,Qi,edt,m2edt;
    int i=glob_index;
    extern double dx_density[N+10];
    edt = 1./t;
    m2edt=m*m*edt;

    rQ(t,&Qr,&Qi,&kappar,&kappai);
    
            
    f[0]=  y[2]*edt;
    f[1] = y[3]*edt;
    
    f[2] = y[0]*m2edt  + (kappar*y[2] - kappai*y[3]-Qr*y[0] + Qi*y[1]);
    f[3] = y[1]*m2edt  + (kappar*y[3] + kappai*y[2]-Qr*y[1] - Qi*y[0]);

    return GSL_SUCCESS;
}

/*************************************************************************/
int jac (double t, const double y[], double *dfdy, 
     double dfdt[], void *params)
{
    extern double m,rcor[N+10],dx_density[N+10];
    double h = 0.001;
    double kappar,kappai;
    double Qr,Qi;
    double Qrp,Qip;
    double kip,krp;
    

    int i=glob_index;
   
    
    

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, NDIM, NDIM);
    gsl_matrix * mat = &dfdy_mat.matrix; 

    rQ(t,&Qr,&Qi,&kappar,&kappai);
     rQ(t+h,&Qrp,&Qip,&krp,&kip);

    gsl_matrix_set (mat, 0, 0, 
                    0.0);
    gsl_matrix_set (mat, 0, 1, 
                    0.0);
    gsl_matrix_set (mat, 0, 2, 
                    1.0/t);
    gsl_matrix_set (mat, 0, 3, 
                    0.0);



    gsl_matrix_set (mat, 1, 0, 
                    0.0);
    gsl_matrix_set (mat, 1, 1, 
                    0.0);
    gsl_matrix_set (mat, 1, 2, 
                    0.0);
    gsl_matrix_set (mat, 1, 3, 
                    1.0/t);


    gsl_matrix_set (mat, 2, 0, 
                    m*m/t-Qr);
    gsl_matrix_set (mat, 2, 1, 
                    Qi);
    gsl_matrix_set (mat, 2, 2, 
                    kappar);
    gsl_matrix_set (mat, 2, 3, 
                    -kappai);
    
    gsl_matrix_set (mat, 3, 0, 
                    -Qi);
    gsl_matrix_set (mat, 3, 1, 
                    m*m/t-Qr);
    gsl_matrix_set (mat, 3, 2, 
                    kappai);
    gsl_matrix_set (mat, 3, 3, 
                    kappar);
    

    dfdt[0] = -1./t/t *y[2];
    dfdt[1] = -1./t/t *y[3];
    dfdt[2] = -m*m/t/t *y[0]
        +(krp-kappar)/h*y[2]-(kip-kappai)/h*y[3]
        -(Qrp-Qr)/h*y[0]+(Qip-Qi)/h*y[1];
    dfdt[3] = -m*m/t/t *y[1]
        +(krp-kappar)/h*y[3]+(kip-kappai)/h*y[2]
        -(Qrp-Qr)/h*y[1]-(Qip-Qi)/h*y[0];
    
    return GSL_SUCCESS;
}




































































