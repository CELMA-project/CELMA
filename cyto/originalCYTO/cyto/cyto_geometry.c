
#include <mpi.h>                                                                           
#include <utilities.h>                                                                     
#include <par_helper.h>                                                                    
                                                                                           
#include "cyto.h"                                                                                          
                                                                                           
#undef DEBUG                                                                               
#undef COMM_STATUS                                                                         
#define  CHECK_STRUCTURE  //Util_PrintStructures(d,p,__LINE__)                                                               
                                                                                           
                                                                                           
                                                                                           
int  Cyto_Geometry(int ispolar, double *rcor, double *edrcor, double *zhval, double *zgval, double *zcor,HDF_DS *d, PARA *p)
{
    int 	
        nx            = 0,
        ny            = 0,
        nz            = 0,
        i             = 0,
        ir            = 0,
        ip            = 0,
        iz            = 0,
        off           = 0;
    
    double 
        r  = 0.,
        zm = 0.,
        gm1= 0.,
        Lz = 0.,
        val = 0.,
        dz = 0.;

    double
        edz  = 0.,
        g    = 0.,
        rlen = 0.,
        xlen = 0.;  
    
 
    nx = d->lnx; 
    ny = d->lny;  
    nz = d->lnz;

    Lz = (p->zmax - p->zmin);

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

 
    if(ISROOT) {
        fprintf(stderr,"The geometry is ");
        if(ispolar)fprintf(stderr,"polar.\n");
        else fprintf(stderr,"rectangular.\n");
        
        COMM(fprintf(stderr,"%d coordsys: %s %d\t",__LINE__,d->coordsys, p->coordsys);
             fprintf(stderr,"with labels: %s %s %s \n",d->dim_label[2],d->dim_label[1],d->dim_label[0]););
        //fprintf(stderr,"%ld %ld %ld \n",d->dims[2],d->dims[1],d->dims[0]);
        //fprintf(stderr,"%ld %ld %ld \n",d->nx,d->ny,d->nz);
        }




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

   
   COMM(if(ISROOT) for(iz=0;iz<d->dims[0];iz++)  fprintf(stderr,"%d %.12f\n",iz,d->coordinate[0][iz] ););
   

   return 0;
   
}
