

#include <mpi.h>
#include <utilities.h>
#include <par_helper.h>

#include "cyto.h"

#undef DEBUG
#undef COMM_STATUS

/*****************************************************/
void  Cyto_Restart(double ***n_0,double ***w_0, double ***U_0, double ***V_0,double ***t_0,
                   double *zcor,HDF_DS *d, PARA *p);

void  Cyto_Ini(double ***n_0,double ***w_0, double ***U_0, double ***V_0,double ***t_0,
               double *zcor,HDF_DS *d, PARA *p);



/*****************************************************/





int  Cyto_Start_Restart(double ***n_0,double ***w_0, double ***U_0, double ***V_0,double ***t_0,
                         double *zcor,HDF_DS *d, PARA *p)
{
    
    if(d->restart == RESTART)
        Cyto_Restart(n_0,w_0, U_0,V_0,t_0,zcor,d,p);
    else if (d->restart == START_FROM_INI || d->restart == DEFAULTSTART )
        Cyto_Ini(n_0,w_0, U_0,V_0,t_0,zcor,d,p);
    else if (d->restart == START_FROM_FILE ) 
    {
        fprintf(stderr, "Start from file  unsupported \n");
        MPI_Finalize();
        return 0;
    }
    else
    {
        fprintf(stderr, "Restart undefined \n");
        MPI_Finalize();
        return -1;
    }
    


    MPI_Bcast(p,1,para_type,ROOT,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    BUGREPORT;

    return 0;
    
    
}


/************************************************************************************************************/


void  Cyto_Restart(double ***n_0,double ***w_0, double ***U_0, double ***V_0,double ***t_0,
                   double *zcor,HDF_DS *d, PARA *p)
{
    int
        nx            = 0,
        ny            = 0,
        nz            = 0;
    
   int
        ir            = 0,
        ip            = 0,
        iz            = 0,
        counter       = 0;

   nx = d->lnx; 
   ny = d->lny;  
   nz = d->lnz;

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

     //Electron temperature
     if(PH_Read3DFieldbyName(t_0,"Te",d,p)> 0  ) 
     {
         if(ISROOT) fprintf(stderr,"Restart success, read field VElectron\n");
     }
     else
     {
         if(ISROOT) fprintf(stderr,"Could not read field Te. Initialising as constant\n");
         FORALL_BD   t_0[iz][ip][ir] =  1.;
     }
     COMM(PH_3D_Write( t_0 ,"Te","RESTART",4,d,p,TRUE);     );
     

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



/*****************************************************************************************************/

void  Cyto_Ini(double ***n_0,double ***w_0, double ***U_0, double ***V_0,double ***t_0,
                         double *zcor,HDF_DS *d, PARA *p)
{
    int
        nx            = 0,
        ny            = 0,
        nz            = 0;
    
   int
        ir            = 0,
        ip            = 0,
        iz            = 0;

   int 
       irg = 0,
       m = 0;
   

   double 
       Lz = 0.;
   

   double
       r  = 0.,
       phi_z = 0.,
       phi_l = 0.;
   
   double 
       x    = 0.,
       y    = 0.,
       z    = 0.;
   
   double
       r2 = 0.;

   double 
       phase = 0.;
   
   

   nx = d->lnx; 
   ny = d->lny;  
   nz = d->lnz;

 
       
   /******************START FROM INI FILE  ****************************************/
   
   Lz = (p->zmax - p->zmin);
   
   FORALL_BD  
   {       
       z = zcor[iz]/Lz;
       r = p->xmin+(d->lnx*d->grid_coords[2]+ir+0.5)*p->dx;
       
       // initial density helically around r0 = p->r0, phi_0 = 2 pi * z
       
       phi_z = 2.*M_PI*z;
       phi_l =  p->ymin+(d->lnx*d->grid_coords[1]+ip+0.5)*p->dy;
       
       x = p->r0 *cos(phi_z)  - r *cos(phi_l);
       y = p->r0*sin(phi_z) - r*sin(phi_l);
       r2 = x*x + y*y;
       
       n_0[iz][ip][ir] =  1. + 0.1*exp(-(r2*p->kappan*p->kappan + 1.*1.*z*z) );
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
   
   
   
   // Electron Temperature has a radial/axial profile, even if constant in time 
   
   FORALL_BD  
   {       
       z = zcor[iz]/Lz;
       r = p->xmin+(d->lnx*d->grid_coords[2]+ir+0.5)*p->dx;
       t_0[iz][ip][ir] =  log(exp(-(0.5*0.5*r*r*p->kappan*p->kappan + 0.1*0.1*z*z) ));
   }
   
   
   BUGREPORT;
   
   /* w   get initialized with zero values */
   
   
   FORALL_BD  
   {
       r = p->xmin+(d->lnx*d->grid_coords[2]+ir+0.5)*p->dx;
       w_0[iz][ip][ir] =   0.; // 0.1*exp(-(r*r*p->kappan*p->kappan) );
   }
   
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







