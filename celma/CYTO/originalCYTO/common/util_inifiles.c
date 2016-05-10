
#undef DEBUG
#undef COMM_STATUS
#include <utilities.h>
#include <string.h>



 
/**************************************************************************************/
/*! \brief  Gets a string and searches it for a description fitting keyword "search".
   Format is as ".......parameter= description[:;,]

   It returns the description of the keyword, if not found in the string it returns a pointer to default.
*/


// Changed to get result as a string allocated by calling  programm 
char* Util_GetDesc(const char *full,const char *search,const char *def,char *result)
{
    char strstart[DEFSTRLEN];
    char lsearch[DEFSTRLEN];
    char *lfull=NULL,*lfullstart=NULL;
    int l;
    
   
    COMM(fprintf(stderr,"%s %s %d: Begin to search for %s in\n\t >%s<\n\n",__FILE__,__func__,__LINE__,search,full););
    // Take a copy of full to modify it; remember to free it at end of function 
    lfullstart= lfull = strdup(full);
    
    
    // Set result to default and return if no search string is given
    snprintf(result,2*DEFSTRLEN,"%s",def);
    if(full ==  NULL)  goto done;

	// if empty full string return search string, e.g. definition is variable name
    if(strcmp(full,"empty") == 0) {
      snprintf(result,DEFSTRLEN,"%s",search);
      goto done;
    } 
    
    // Directly search for tokenised pattern 
    snprintf(lsearch,DEFSTRLEN,"%s:",search);
    
    while ( (lfull = strstr(lfull,lsearch)) != (char*)NULL)
    {
        BUGREPORT;
        /* String from first occurence of lsearch */
        snprintf(strstart,DEFSTRLEN,"%s",strstr(lfull,lsearch));
        
        COMM( fprintf(stderr,"%s %s %d: search >%s< at  >%s<\n",__FILE__,__func__,__LINE__,lsearch,lfull););
       
		// goto occurrence of ":" 
        if( (lfull = strchr(lfull,':')) == NULL) break;
         
        if( lfull++ == NULL) break;
        while (*lfull == ' ')  if( lfull++ == NULL) break;
        

        strncpy(result,lfull,l= MIN(strcspn(lfull,";,"),DEFSTRLEN)); // max DEFSTRLEN characters in result
        BUGREPORT;
        result[l]='\0'; /* set last charater to zero */
        BUGREPORT;
        COMM(fprintf(stderr,"%s %s %d: searched >%s< result is >%s<\n",__FILE__,__func__,__LINE__,lsearch,result););
        goto done;
    }
    
    COMM(fprintf(stderr,"%s %s %d: %s: is not identified as a token.\n",__FILE__,__func__,__LINE__,search););
    
    BUGREPORT;
    
    done:
    BUGREPORT;
    free(lfullstart);
    
    BUGREPORT;
    return result;
    
}
   /**************************************************************************************/

int Util_ReadIniFile(HDF_DS *data,PARA *para,char *ininame,int rw)
{
    /* 
       Reads in an *ini file and fills the structures data and para with the values
       data   : data structure
       para   : para structure to be filled
       ininame: name of .ini file
       rw     : FALSE: Read from .ini file 
       TRUE : Write DATA and PARA to ininame in .ini format 
    */ 
    FILE 
        *out   = stdout;
    int 
        xi     = 1,
        yi     = 0,
        zi     = 2,
        i      = 0,
        result = 0;

    char 
        var_name[DEFSTRLEN] = {0},
        desc_string[DEFSTRLEN] = {0};
    static char 
        corname[3][2] = {"x","y","z"};

    static char 
        Corname[3][2] = {"X","Y","Z"};

    int  
        hyperindex[] =  {0,1,2};

    char 
        *stringp = NULL;

    /* Check if file is there......*/
  
    if(!rw)
    {
        if(fopen(ininame,"r") == NULL)
        {
            fprintf(stderr,"Proc %d: Couldn't open %s for reading!\n",(int)data->this_process,ininame);
            return -1;
        }
    }

    BUGREPORT;

    if(rw) out = fopen(ininame,"w");
    fprintf(out,"#\n%s\n#\n",para->codedesc);
    fprintf(out,"#\n#Dimensions:\n#\n");
    if(rw) fclose(out);
    BUGREPORT;

  
    data->nx = (int)Util_ReadValueByName("nx",ininame,0,(double)data->nx,"Number of grid points in x",FALSE,rw);
    if(!rw) data->rank = 1;
  
    if((rw && (data->rank > 1)) || !rw)   
    {
        data->ny = (int)Util_ReadValueByName("ny",ininame,0,(double)data->ny,"Number of grid points in y",FALSE,rw);
        if(!rw) data->rank = 2;
    }
    if((rw && (data->rank > 2) ) || !rw)   
    {
        if( (data->nz = (int)Util_ReadValueByName("nz",ininame,0,(double)data->nz,"Number of grid points in z",FALSE,rw)) != 0)
        {
            data->rank = 3;
            xi = 0;yi =1;zi=2;
        }
    }


    COMM(fprintf(stderr,"Process %d: Util_ReadIniFile: datarank = %ld,(%d,%d,%d)\n",data->this_process,data->rank,
                 (int)data->nx,(int)data->ny,(int)data->nz););
  
 
    if(rw) out = fopen(ininame,"a");
    fprintf(out,"#\n#Units:\n#\n");
    if(rw) fclose(out);
    para->unit_time        = Util_ReadValueByName("time_unit",ininame,1.,para->unit_time,"Time unit in mu seconds",TRUE,rw);
    para->unit_length    = Util_ReadValueByName("lenght_unit",ininame,1.,para->unit_length,"Length unit in m",TRUE,rw);

    if(rw) out = fopen(ininame,"a");
    fprintf(out,"#\n#Domain:\n#\n");
    if(rw) fclose(out);

    hyperindex[0] = xi;
    hyperindex[1] = yi;
    hyperindex[2] = zi;

    BUGREPORT;

    para->xmin = Util_ReadValueByName("xmin",ininame,0.,para->xmin,"",TRUE,rw);
    para->xmax = Util_ReadValueByName("xmax",ininame,1.,para->xmax,"",TRUE,rw);
 

    if(data->rank > 1)
    {
        para->ymin = Util_ReadValueByName("ymin",ininame,0.,para->ymin,"times PI for poloidal coordinates and ymax = 1,2",TRUE,rw);
        para->ymax = Util_ReadValueByName("ymax",ininame,1.,para->ymax,"times PI for poloidal coordinates and ymax = 1,2",TRUE,rw);
    }


    if(data->rank > 2)
    {   
        para->zmin = Util_ReadValueByName("zmin",ininame,0.,para->zmin,"",TRUE,rw);
        para->zmax = Util_ReadValueByName("zmax",ininame,1.,para->zmax,"",TRUE,rw);
    }
  



    if(rw) out = fopen(ininame,"a");
    fprintf(out,"#\n#Processor Grid:\n#\n");
    if(rw) fclose(out);


    data->N[0] = (int)Util_ReadValueByName("NZ",ininame,0,(double)data->N[0],"Number of processes in Z (Zero if programme decides by itself)",FALSE,rw);
    data->N[1] = (int)Util_ReadValueByName("NY",ininame,1,(double)data->N[1],"Number of processes in Y (Always 1)",FALSE,rw);
    data->N[2] = (int)Util_ReadValueByName("NX",ininame,1,(double)data->N[2],"Number of processes in X (Zero if programme decides by itself)",FALSE,rw);

    if(rw) out = fopen(ininame,"a");
    fprintf(out,"#\n# Info on coordinate system:\n#\n");
    if(rw) fclose(out);
    BUGREPORT;

    para->coordsys  = 
        (int)Util_ReadValueByName("coordsys",ininame,0.,(double)para->coordsys,"Coordinate system 0:cartesian, 1:poloidal, 2:cylindrical, 3:spherical",FALSE,rw);

    para->r_spacing  = 
        (int)Util_ReadValueByName("r_spacing",ininame,0.,(double)para->r_spacing,"Equidistant [0] or cos [1]  points in x",FALSE,rw);
    
    para->r_offset  = 
        (int)Util_ReadValueByName("r_offset",ininame,0.,(double)para->r_offset,"r Offset when using cos distributed points",FALSE,rw);

    
    para->z_spacing  = 
        (int)Util_ReadValueByName("z_spacing",ininame,0.,(double)para->z_spacing,"Equidistant [0] or cos [1]  points in z",FALSE,rw);
    
    para->z_offset  = 
            (int)Util_ReadValueByName("z_offset",ininame,0.,(double)para->z_offset,"z Offset when using cos distributed points",FALSE,rw);
    
    

    if((para->coordsys == POLOIDAL || para->coordsys == CYLINDRICAL || para->coordsys == SPHERICAL) && ( para->ymax == 1. ||  para->ymax == 2.))
    {para->ymin *= M_PI; para->ymax *= M_PI;}

    if((para->coordsys == SPHERICAL)&& ( para->zmax == 1. ||  para->zmax == 2.))
    {para->zmin *= M_PI; para->zmax *= M_PI;}

    /*

 if(!rw) 
    {
        for(i=0;i<data->rank;i++)
        {
            snprintf(var_name,DEFSTRLEN,"N%s",Corname[i]);
            data->N[hyperindex[i]] = 
                MAX((int)Util_ReadValueByName(var_name,ininame,1,(double)data->N[hyperindex[i]],"",FALSE,rw),1);
        }
    }
    */

    if(rw) out = fopen(ininame,"a"); 
    fprintf(out,"#\n# Time Integration Parameters\n#\n");
    if(rw) fclose(out);
    BUGREPORT;
  

    para->dt       = Util_ReadValueByName("dt",ininame,0.001,para->dt,"Delta T",TRUE,rw);
    para->end_time = Util_ReadValueByName("end_time",ininame,1.,para->end_time,"Simulation stops at this time",TRUE,rw);
    para->out_time = Util_ReadValueByName("out_time",ininame,0.1,para->out_time,"Time between small outputs",TRUE,rw);
    para->otmult   = (int)Util_ReadValueByName("outmult",ininame,1.,(double)para->otmult,"Number of small outputs  before output of fields",FALSE,rw);


    /* Primary physics parameters */

   if(rw) out = fopen(ininame,"a");
   fprintf(out,"#\n# Primary Physics Parameters\n#\n");
   fprintf(out,"# if specified (Te > 0.) these will be used to calculate both\n#\n");
   fprintf(out,"#\t -secondary parameters and \n# \t -Parameters for equation unless initialised with value nonzero \n#\n");
   if(rw) fclose(out);
   
   
   
   stringp =  Util_GetDesc(para->Prim_Phys,"Ti","-",desc_string);  
   para->Ti  = Util_ReadValueByName("Ti",ininame,0.,para->Ti,stringp,TRUE,rw);


   stringp = Util_GetDesc(para->Prim_Phys,"Te","-",desc_string);
   para->Te  = Util_ReadValueByName("Te",ininame,0.,para->Te,stringp, TRUE,rw);
                                   
                                   
   stringp = Util_GetDesc(para->Prim_Phys,"B0","-",desc_string);
   para->B0  = Util_ReadValueByName("B0",ininame,0.,para->B0,stringp,TRUE,rw);
                                    
                                    
   stringp = Util_GetDesc(para->Prim_Phys,"n0","-",desc_string);
   para->n0  = Util_ReadValueByName("n0",ininame,0.,para->n0,stringp,TRUE,rw);
                                   
                                   
   stringp = Util_GetDesc(para->Prim_Phys,"Mi","-",desc_string);
   para->Mi  = Util_ReadValueByName("Mi",ininame,0.,para->Mi,stringp,TRUE,rw);


   stringp = Util_GetDesc(para->Prim_Phys,"Z","-",desc_string);
   para->Z   = Util_ReadValueByName("Z",ininame,1.,para->Z,stringp,TRUE,rw); // default 1
                                    
                                    
   stringp = Util_GetDesc(para->Prim_Phys,"p_n","-",desc_string);
   para->p_n = Util_ReadValueByName("p_n",ininame,0.,para->p_n,stringp,TRUE,rw);


   stringp = Util_GetDesc(para->Prim_Phys,"part_source","-",desc_string);
   para->part_source = Util_ReadValueByName("part_source",ininame,0.,para->part_source,stringp,TRUE,rw);

   stringp = Util_GetDesc(para->Prim_Phys,"temp_source","-",desc_string);
   para->temp_source = Util_ReadValueByName("temp_source",ininame,0.,para->temp_source,stringp,TRUE,rw);


   /* Secondary  parameters */
  
   if(rw) out = fopen(ininame,"a");
   fprintf(out,"#\n# Derived Parameters\n#\n");
   fprintf(out,"# only for information\n#\n");
   if(rw) fclose(out);

   stringp = Util_GetDesc(para->Sec_Phys,"rho_s","-",desc_string);
   para->rho_s    = Util_ReadValueByName("rho_s",ininame,0.,para->rho_s,stringp,TRUE,rw); 
                                        
   stringp = Util_GetDesc(para->Sec_Phys,"rho_i","-",desc_string);                                
   para->rho_i    = Util_ReadValueByName("rho_i",ininame,0.,para->rho_i,stringp,TRUE,rw);

   stringp = Util_GetDesc(para->Sec_Phys,"rho_e","-",desc_string);
   para->rho_e    = Util_ReadValueByName("rho_e",ininame,0.,para->rho_e,stringp,TRUE,rw); 

   stringp = Util_GetDesc(para->Sec_Phys,"omega_ci","-",desc_string);
   para->omega_ci    = Util_ReadValueByName("omega_ci",ininame,0.,para->omega_ci,stringp,TRUE,rw);
                                      
   stringp = Util_GetDesc(para->Sec_Phys,"omega_ce","-",desc_string);                                   
   para->omega_ce    = Util_ReadValueByName("omega_ce",ininame,0.,para->omega_ce,stringp,TRUE,rw); 
                                        
   stringp = Util_GetDesc(para->Sec_Phys,"omega_pi","-",desc_string);                            
   para->omega_pi    = Util_ReadValueByName("omega_pi",ininame,0.,para->omega_pi,stringp,TRUE,rw); 

   stringp = Util_GetDesc(para->Sec_Phys,"v_thi","-",desc_string);
   para->v_thi    = Util_ReadValueByName("v_thi",ininame,0.,para->v_thi,stringp,TRUE,rw); 
                                                                              
   stringp = Util_GetDesc(para->Sec_Phys,"v_the","-",desc_string);                                      
   para->v_the    = Util_ReadValueByName("v_the",ininame,0.,para->v_the,stringp,TRUE,rw); 
                                       
   stringp = Util_GetDesc(para->Sec_Phys,"c_s","-",desc_string);                                    
   para->c_s    = Util_ReadValueByName("c_s",ininame,0.,para->c_s,stringp,TRUE,rw); 
                                     
   stringp = Util_GetDesc(para->Sec_Phys,"v_alfven","-",desc_string);                                   
   para->v_alfven    = Util_ReadValueByName("v_alfven",ininame,0.,para->v_alfven,stringp,TRUE,rw); 

   stringp = Util_GetDesc(para->Sec_Phys,"n_n","-",desc_string);                                   
   para->n_n    = Util_ReadValueByName("n_n",ininame,0.,para->n_n,stringp,TRUE,rw); 



   /* Geometry  parameters */

   if(rw) out = fopen(ininame,"a");
   fprintf(out,"#\n# Primary Geometry Parameters\n#\n\n");
   if(rw) fclose(out);


   stringp = Util_GetDesc(para->Prim_Geom,"R0","-",desc_string);
   para->R0     = Util_ReadValueByName("R0",ininame,0.,para->R0,stringp,TRUE,rw);

   stringp = Util_GetDesc(para->Prim_Geom,"a","-",desc_string);
   para->a     = Util_ReadValueByName("a",ininame,0.,para->a, stringp,TRUE,rw);

   stringp = Util_GetDesc(para->Prim_Geom,"lpar","-",desc_string);
   para->lpar    = Util_ReadValueByName("lpar",ininame,0.,para->lpar,stringp,TRUE,rw);

   stringp = Util_GetDesc(para->Prim_Geom,"ln","-",desc_string);
   para->ln   = Util_ReadValueByName("ln",ininame,0.,para->ln,stringp,TRUE,rw);

   stringp = Util_GetDesc(para->Prim_Geom,"lTi","-",desc_string);
   para->lTi    = Util_ReadValueByName("lTi",ininame,0.,para->lTi,stringp,TRUE,rw);

   stringp = Util_GetDesc(para->Prim_Geom,"lTe","-",desc_string);
   para->lTe    = Util_ReadValueByName("lTe",ininame,0.,para->lTe,stringp,TRUE,rw);

   if(rw) out = fopen(ininame,"a");
   fprintf(out,"#\n# ----------------------------------------------------------------------------\n#\n\n");
   if(rw) fclose(out);


   if(rw) out = fopen(ininame,"a");
   fprintf(out,"#\n# Dimensionless geometry parameters that can be changed\n#\n\n");
   if(rw) fclose(out);

   stringp = Util_GetDesc(para->Prim_Geom,"q0","-",desc_string);
   para->q0	    = Util_ReadValueByName("q0",ininame,0.,para->q0,stringp,TRUE,rw);

   stringp = Util_GetDesc(para->Prim_Geom,"shat0","-",desc_string);
   para->shat0	    = Util_ReadValueByName("shat0",ininame,0.,para->shat0,stringp,TRUE,rw);

   /* Equation parameters */ 

    if(rw) out = fopen(ininame,"a");
    fprintf(out,"#\n#Parameters for equation, directly specified \n#\n");
    if(rw) fclose(out);
    BUGREPORT; 
    
    stringp = Util_GetDesc(para->desc,"k0","-",desc_string);
    para->k0     = Util_ReadValueByName("k0",ininame,0.,para->k0,stringp,TRUE,rw);
    BUGREPORT;

    stringp = Util_GetDesc(para->desc,"r0","-",desc_string);
    para->r0     = Util_ReadValueByName("r0",ininame,0.,para->r0,stringp,TRUE,rw); 

    BUGREPORT;
    stringp = Util_GetDesc(para->desc,"sigma","-",desc_string);
    para->sigma  = Util_ReadValueByName("sigma",ininame,0.,para->sigma,stringp,TRUE,rw);
    
    stringp = Util_GetDesc(para->desc,"alpha","-",desc_string);
    para->alpha  = Util_ReadValueByName("alpha",ininame,0.,para->alpha,stringp,TRUE,rw);
                                       
                                       
    stringp = Util_GetDesc(para->desc,"beta","-",desc_string);
    para->beta   = Util_ReadValueByName("beta",ininame,0.,para->beta,stringp,TRUE,rw);
                                        
    stringp = Util_GetDesc(para->desc,"gamma","-",desc_string);
    para->gamma  = Util_ReadValueByName("gamma",ininame,0.,para->gamma,stringp,TRUE,rw);
                                        
    stringp = Util_GetDesc(para->desc,"delta","-",desc_string);
    para->delta  = Util_ReadValueByName("delta",ininame,0.,para->delta,stringp,TRUE,rw);
                                        
    stringp = Util_GetDesc(para->desc,"muehat","-",desc_string);
    para->muehat = Util_ReadValueByName("muehat",ininame,0.,para->muehat,stringp,TRUE,rw);
                                      
    stringp = Util_GetDesc(para->desc,"betahat","-",desc_string);
    para->betahat= Util_ReadValueByName("betahat",ininame,0.,para->betahat,stringp,TRUE,rw);
                                      
    stringp = Util_GetDesc(para->desc,"adrhos","-",desc_string);
    para->adrhos = Util_ReadValueByName("adrhos",ininame,0.,para->adrhos,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"nu","-",desc_string);
    para->nu     = Util_ReadValueByName("nu",ininame,0.,para->nu,stringp,TRUE,rw);

    /* Profiles */
    stringp = Util_GetDesc(para->desc,"kappan","-",desc_string);
    para->kappan = Util_ReadValueByName("kappan",ininame,0.,para->kappan,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"kappat","-",desc_string);
    para->kappat = Util_ReadValueByName("kappat",ininame,0.,para->kappat,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"bprof","-",desc_string);
    para->bprof  = Util_ReadValueByName("bprof",ininame,0,para->bprof,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"nprof", "-",desc_string);
    para->nprof  = Util_ReadValueByName("nprof",ininame,0.,para->nprof,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"tprof","-",desc_string);
    para->tprof  = Util_ReadValueByName("tprof",ininame,0.,para->tprof,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"phiprof","-",desc_string);
    para->phiprof= Util_ReadValueByName("phiprof",ininame,0.,para->phiprof,stringp,TRUE,rw);
  

    stringp = Util_GetDesc(para->desc,"xbnd","-",desc_string);
    para->xbnd   = (int)Util_ReadValueByName("xbnd",ininame,0.,para->xbnd,stringp,FALSE,rw);

    stringp = Util_GetDesc(para->desc,"ybnd","-",desc_string);
    para->ybnd   = (int)Util_ReadValueByName("ybnd",ininame,0.,para->ybnd,stringp,FALSE,rw);

    stringp = Util_GetDesc(para->desc,"zbnd","-",desc_string);
    para->zbnd   = (int)Util_ReadValueByName("zbnd",ininame,0.,para->zbnd,stringp,FALSE,rw);

    stringp = Util_GetDesc(para->desc,"uvortex","-",desc_string);
    para->uvortex = Util_ReadValueByName("uvortex",ininame,0.,para->uvortex,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"radius","-",desc_string);
    para->radius  = Util_ReadValueByName("radius",ininame,0.,para->radius,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"epsilon","-",desc_string);
    para->epsilon = Util_ReadValueByName("epsilon",ininame,0.,para->epsilon,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"circulation","-",desc_string);
    para->circulation = Util_ReadValueByName("circulation",ininame,0.,para->circulation,stringp,TRUE,rw);


    if(rw) out = fopen(ininame,"a");
    fprintf(out,"#\n# ------------------------ numerical parameters ----------------------------------------------\n#\n\n");
    if(rw) fclose(out);

    stringp = Util_GetDesc(para->desc,"mue_w","-",desc_string);
    para->mue_w  = Util_ReadValueByName("mue_w",ininame,0.,para->mue_w,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"mue_n","-",desc_string);
    para->mue_n  = Util_ReadValueByName("mue_n",ininame,0.,para->mue_n,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"mue_t","-",desc_string);
    para->mue_t  = Util_ReadValueByName("mue_t",ininame,0.,para->mue_t,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"source","-",desc_string);
     para->source = Util_ReadValueByName("source",ininame,0.,para->source,stringp,TRUE,rw);

    stringp = Util_GetDesc(para->desc,"limiter","-",desc_string);
     para->limiter= Util_ReadValueByName("limiter",ininame,0.,para->limiter,stringp,TRUE,rw);


     if(rw) out = fopen(ininame,"a");
     fprintf(out,"#\n# ------------------------- numerical switches ------------------------------\n#\n\n");
     if(rw) fclose(out);

     stringp = Util_GetDesc(para->desc,"hm_nl","-",desc_string);
      para->hm_nl  = Util_ReadValueByName("hm_nl",ininame,0.,para->hm_nl,stringp,TRUE,rw);

      stringp = Util_GetDesc(para->desc,"exb_ll","-",desc_string);
      para->exb_ll = Util_ReadValueByName("n_nl",ininame,0.,para->exb_ll,stringp,TRUE,rw);


      stringp = Util_GetDesc(para->desc,"dt_pol","-",desc_string);
      para->dt_pol = Util_ReadValueByName("dt_polarization",ininame,0.,para->dt_pol,stringp,TRUE,rw);


     /* Erzeugen der Verteilung
     */

#ifdef OLD_STARTUP
    if(rw) out = fopen(ininame,"a");
    fprintf(out,"#\n# Initial condition\n#\n");
    fprintf(out,"#\t 0: Dipol: uvortex,radius\n"\
            "#\t 1: Wave: x,\t y,\t amp\n"\
            "#\t 2: GaussNoise:x,\t y,\t amp\n"\
            "#\t 3: Point:x,\t y,\t amp\n"\
            "#\t 4: Gauss: pos_x,\t pos_y,\t width_x,\t width_y,\t amp\n"\
            "#\t 5: White Noise: amp\n"\
            "#\t 6: Line: y,amp\n"\
            "#\t 7: Lamb Dipole: uvortex,radius\n"\
            "#\t 8: Cylinder: radius, amp\n"\
            "#\t10: read_asc\n"\
            "#\t11: drift_Util_Noise\n"\
            "#\t12: map_Util_Noise\n"\
            "#\t13: const_val\n"\
            "#\t14: read_bin\n"\
            "#\t15: read_simple_hdf\n");
    if(rw) fclose(out);
   
 
    para->verteilung = (int)Util_ReadValueByName("verteilung",ininame,0.,(double)para->verteilung,"select initial condition",FALSE,rw);

    if(rw && data->anzahl > 0) initialfields = data->anzahl;
    else initialfields = 3;
  
    for(i=0;i<initialfields;i++)
    {     
        if(rw) out = fopen(ininame,"a");
        fprintf(out,"# Data for field >%s< with identifier %d\n#\n",data->names[i],i);
        if(rw) fclose(out);

        snprintf(var_name,DEFSTRLEN,"randbedingung%d",i);
        para->boundary[i]     = (int)Util_ReadValueByName(var_name,ininame,0.,(double)para->boundary[i],"",FALSE,rw);

        snprintf(var_name,DEFSTRLEN,"bdvala%d",i);
        para->bdval[i][0]     = Util_ReadValueByName(var_name,ininame,0.,para->bdval[i][0],"",stringp,TRUE,rw);

        snprintf(var_name,DEFSTRLEN,"bdvalb%d",i);
        para->bdval[i][1]     = Util_ReadValueByName(var_name,ininame,0.,para->bdval[i][1],"",stringp,TRUE,rw);

        snprintf(var_name,DEFSTRLEN,"amp%d",i);
        para->amp[i]     = Util_ReadValueByName(var_name,ininame,0.,para->amp[i],"",stringp,TRUE,rw);

        snprintf(var_name,DEFSTRLEN,"width_x%d",i);
        para->width_x[i]     = Util_ReadValueByName(var_name,ininame,0.,para->width_x[i],"",stringp,TRUE,rw);

        snprintf(var_name,DEFSTRLEN,"width_y%d",i);
        para->width_y[i]     = Util_ReadValueByName(var_name,ininame,0.,para->width_y[i],"",stringp,TRUE,rw);

        snprintf(var_name,DEFSTRLEN,"pos_x%d",i);
        para->pos_x[i]     = Util_ReadValueByName(var_name,ininame,0.,para->pos_x[i],"",stringp,TRUE,rw);

        snprintf(var_name,DEFSTRLEN,"pos_y%d",i);
        para->pos_y[i]     = Util_ReadValueByName(var_name,ininame,0.,para->pos_y[i],"",stringp,TRUE,rw);

        snprintf(var_name,DEFSTRLEN,"amp_random%d",i);
        para->amp_random[i]     = Util_ReadValueByName(var_name,ininame,0.,para->amp_random[i],"",stringp,TRUE,rw);

        snprintf(var_name,DEFSTRLEN,"width_random_x%d",i);
        para->width_random_x[i]     = Util_ReadValueByName(var_name,ininame,0.,para->width_random_x[i],"",stringp,TRUE,rw);

        snprintf(var_name,DEFSTRLEN,"width_random_y%d",i);
        para->width_random_y[i]     = Util_ReadValueByName(var_name,ininame,0.,para->width_random_y[i],"",stringp,TRUE,rw);

        if(rw) out = fopen(ininame,"a");
        fprintf(out,"#\n");

        if(rw) fclose(out);
    }
#endif /* OLD_STARTUP */
    

    return result;

}



/**************************************************************************/

double Util_ReadValueByName(const char *name, char *file, double def,double hadval,const char *desc,int isfloat, int rw)
{
    /* Reads value from an ini file or writes out an default .ini file
       name    : name of variable
       file    : name of .ini file
       def     : default value
       hadval  : value already known
       desc    : description of variable
       isfloat : TRUE: Float, FALSE Integer value
       rw      : TRUE   write to file
       FALSE  read from file

       On exit function returns value of variable
    */ 

    FILE *in,*out;
    char buf[DEFSTRLEN];
    char ldesc[DEFSTRLEN];
    char *token;
    int search=TRUE;
    double val=0.;

    BUGREPORT;

    if(!rw)
    {
        val = def;

        in = fopen(file,"r");

        while( fgets(buf,DEFSTRLEN-1,in) != NULL && search)
        {
            COMM(fprintf(stderr,"%s %s %d (search >%s<: %s\n",__FILE__,__func__,__LINE__,name,buf));
            if( (*buf != '#') && strstr(buf,name)) 
                if((token =strtok(buf," =;,")) != NULL)
                {
                    COMM(fprintf(stderr,"%s %s %d :Search for >%s<. Actual token:  >%s<.\n",
                                 __FILE__,__func__,__LINE__,name,token););
                    if(strlen(name) == strlen(token) && !strcmp(name,token))
                    {
                        token = strtok(NULL," =;,");
                        COMM(fprintf(stderr,"%s %s %d : >%s< value found:  >%s<\n",
                                     __FILE__,__func__,__LINE__,name,token););
                        val =  Util_Str2Double(token);
                        search = FALSE;
                    }
                }
        }
        fclose(in);
    
        printf(" %-20s = %f; %s\n",name,val,desc);
        BUGREPORT;

    }
    else
    {
      snprintf(ldesc,DEFSTRLEN,"%s",desc);
        
        BUGREPORT;
        Util_TrimString(ldesc);
        if(strcmp(ldesc,"-"))
        {
            out = fopen(file,"a");
            val = hadval;
          
            if(isfloat) fprintf(out," %-20s = %f; %s\n",name,val,ldesc);
            else fprintf(out," %-20s = %d; %s\n",name,(int)val,ldesc);
            fclose(out); 
        }
        BUGREPORT;

    }

    return val;
}

/**************************************************************************/
double Util_Str2Double(char *str)
{
    double sign = 1.;
    double value = 0., fac = 1.;

    while(!isdigit(*str)) 
	{str++;if(*str == '\0') return 0.;}

/* First occurrence of a digit, check sign and if it's starting with a dot....*/

    if(*(str-1) == '.') str--;

    while (isspace(*(str--)));

    if(*str == '-') sign = -1.; 

    while(!isdigit(*str)) str++;

    if(*(str-1) == '.') str--;

    while(isdigit(*str)) {value *=10.; value+=*str-'0';str++;}


    if(*str == '.') {
        str++;	if(*str != '\0') 
                    while(isdigit(*str)) {fac *=.1;value+=fac*(*str-'0');str++;}	
	}

    value *=sign;
	
/*printf("%.5g\n",value);*/
	
    return(value);
}
/**************************************************************************/
int  Util_SetupGeom(double **rcor,double **edr, 
                    double **hval,double **vval,double **gval,
                    double **rs, double **edrs, 
                    double **hs, double **vs, double **gs,
                    double **zcor,double **zhval,double **zgval,
                    double **norm, HDF_DS *data,PARA *p,int rindexl)
{
    
    int 
        nr =  0,
        ir = 0,
        iz = 0,
        off= 1,
        oc,
        nz = 0,
        ispolar = FALSE;
    double 
        r0,val,g,gm1,rm,
        dz,edz,zm,
        rlen,xlen,
        dx,edx,stag;
    
  
    printf("x [%f;%f]   dx = %f  n*p->dx %f\n",p->xmin,p->xmax,p->dx,data->nx*p->dx);
    printf("y [%f;%f]   dy = %f  n*p->dy %f\n",p->ymin,p->ymax,p->dy,data->ny*p->dy);
 
    BUGREPORT;
    
    r0 = (p->dx*((double)data->start[rindexl]+0.5)) + p->xmin;

    nr = data->elements[rindexl];
    nz = data->elements[0];
   
    BUGREPORT;

    if(*rcor != NULL) free(*rcor-off);   BUGREPORT; *rcor  = Util_DVector(nr,off);
    if(*edr  != NULL) free(*edr-off) ;   BUGREPORT; *edr   = Util_DVector(nr,off);
    if(*hval != NULL) free(*hval-off);   BUGREPORT; *hval  = Util_DVector(nr,off);
    if(*vval != NULL) free(*vval-off);   BUGREPORT; *vval  = Util_DVector(nr,off);
    if(*gval != NULL) free(*gval-off);   BUGREPORT; *gval  = Util_DVector(nr,off);
    if(*norm != NULL) free(*norm-off);   BUGREPORT; *norm  = Util_DVector(nr,off);
    if(*hs   != NULL) free(*hs-off);     BUGREPORT; *hs    = Util_DVector(nr,off);
    if(*vs   != NULL) free(*vs-off);     BUGREPORT; *vs    = Util_DVector(nr,off);
    if(*gs   != NULL) free(*gs-off);     BUGREPORT; *gs    = Util_DVector(nr,off);
    if(*rs   != NULL) free(*rs-off);     BUGREPORT; *rs    = Util_DVector(nr,off);
    if(*edrs != NULL) free(*edrs-off);   BUGREPORT; *edrs  = Util_DVector(nr,off);

    if(*zcor  != NULL) free(*zcor-off);    BUGREPORT; *zcor   = Util_DVector(nz,off);
    if(*zhval != NULL) free(*zhval-off);   BUGREPORT; *zhval  = Util_DVector(nz,off);
    if(*zgval != NULL) free(*zgval-off);   BUGREPORT; *zgval  = Util_DVector(nz,off);


    BUGREPORT;

    // Check for z-Tschebischeff
    if (p->z_spacing == 1)
    {
        fprintf(stderr,"%s: Setup z-Tschebicheff\n",__FILE__);
        
        
        zm = 0.5*(p->zmax+p->zmin);
        off = p->z_offset;

        dz = (p->zmax-p->zmin)/(double)(nz+2*off);
        edz = 1./dz;
        

        BUGREPORT;
      
        rlen = p->zmax-p->zmin;
        xlen = cos((double)(off)*dz)-cos((double)(off+nz)*dz);
      
        g = xlen/rlen;
        gm1 = 1./g;

        fprintf(stderr,"#dz = %g, off = %d  gm1 = %f zm = %f\n",dz,off,gm1,zm);
        stag = 0.5;
      
        BUGREPORT;

        // calculate locations of z planes

        for(iz=-1;iz<=nz;iz++)
        {
            (*zcor)[iz]  = fabs(-cos(dz*((double)(iz+off)+stag))*gm1 + zm);
            COMM( fprintf(stderr,"%d %.12f\n",iz,(*zcor)[iz] ););
        }
     
        for(iz=-1;iz<=nz;iz++)
        {
            val = g*((*zcor)[iz]-rm);	
            (*zhval)[iz] = g/sqrt(1.-val*val)*edz;
            (*zgval)[iz] = g*g*val*edz/((1.-val*val)*sqrt(1.-val*val));
        }

    }
    


    // Check for r-Tschebischeff
    if (p->r_spacing == 1)
     {
        fprintf(stderr,"%s: Setup radial Tschebicheff\n",__FILE__);

        if(ispolar)
        {
            if(p->xmin == 0.) 
            {
                oc = nr ;
                rm = 0.;
            }
        }
        else
        {
            oc = 0;
            rm = 0.5*(p->xmax+p->xmin);
        }
        
        BUGREPORT;
        off = p->r_offset;
        
        
        dx = M_PI/(double)(nr+2*off+oc);
        edx = 1./dx;
        

        BUGREPORT;
      
        rlen = p->xmax-p->xmin;
        xlen = cos((double)(off+oc)*dx)-cos((double)(off+oc+nr)*dx);
      
        g = xlen/rlen;
        gm1 = 1./g;

        fprintf(stderr,"#dx = %g, off = %d oc=%d gm1 = %f rm = %f\n",dx,off,oc,gm1,rm);
        stag = 0.5;
      
        BUGREPORT;
        for(ir=-1;ir<=nr;ir++)
        {
            (*rcor)[ir]  = fabs(-cos(dx*((double)(ir+off+oc)+stag))*gm1 + rm);
            COMM( fprintf(stderr,"%d %.12f\n",ir,(*rcor)[ir] ););
        }
     
 
        for(ir=-1;ir<=nr;ir++)
            (*rs)[ir] = fabs(-cos(dx*((double)(ir+off+oc)+stag+0.5))*gm1 + rm);
   
        for(ir=-1;ir<=nr;ir++)
        {
            val = g*((*rcor)[ir]-rm);	
            (*hval)[ir] = g/sqrt(1.-val*val)*edx;
            (*gval)[ir] = g*g*val*edx/((1.-val*val)*sqrt(1.-val*val));
        }

        for(ir=-1;ir<=nr;ir++)
        {
            val = g*((*rs)[ir]-rm);
            (*hs)[ir] = g/sqrt(1.-val*val)*edx;
            (*gs)[ir] = g*g*val*edx/((1.-val*val)*sqrt(1.-val*val));
        }

        BUGREPORT;
        for(ir=-1;ir<=nr;ir++) (*edr)[ir]  = 1./((*rcor)[ir]);    
        for(ir=-1;ir<=nr;ir++) (*edrs)[ir] = 1./((*rs)[ir]);    
        BUGREPORT;

        for(ir=-1;ir<=nr;ir++) (*vval)[ir] = ((*edr)[ir]);
        for(ir=-1;ir<=nr;ir++) ( *vs )[ir] = ((*edrs)[ir]);
    }
    else if((strncmp(data->coordsys,"pol",3)== 0) ||
            (strncmp(data->coordsys,"cyl",3)== 0) && p->r_spacing == 0)
    {
        fprintf(stderr,"%s: Setup polar equidistant\n",__FILE__);
        ispolar = TRUE;
        for(ir=-1;ir<=nr;ir++) (*rcor)[ir] = p->dx*(double)(ir) + r0;
        for(ir=-1;ir<=nr;ir++) (*hval)[ir] = 1./p->dx;
        for(ir=-1;ir<=nr;ir++) (*gval)[ir] = 0.;
        for(ir=-1;ir<=nr;ir++) (*rs)[ir]   = ((*rcor)[ir])+p->dx*.5;
        for(ir=-1;ir<=nr;ir++) (*hs)[ir]   = 1./p->dx;
        for(ir=-1;ir<=nr;ir++) (*gs)[ir]   = 0.;
        for(ir=-1;ir<=nr;ir++) (*edr)[ir]  = 1./((*rcor)[ir]);    
        for(ir=-1;ir<=nr;ir++) (*edrs)[ir] = 1./((*rs)[ir]);
        for(ir=-1;ir<=nr;ir++) (*vval)[ir] = ((*edr)[ir]);
        for(ir=-1;ir<=nr;ir++) ( *vs )[ir] = ((*edrs)[ir]);
    }
    else 
    {
        fprintf(stderr,"%s: Setup Cartesian\n",__FILE__);
        for(ir=-1;ir<=nr;ir++) (*rcor)[ir] = 1.;     
        for(ir=-1;ir<=nr;ir++) (*hval)[ir] = 1./p->dx;
        for(ir=-1;ir<=nr;ir++) (*gval)[ir] = 0.;
        for(ir=-1;ir<=nr;ir++) (*rs)[ir]   = 1.;
        for(ir=-1;ir<=nr;ir++) (*hs)[ir]   = 1./p->dx;
        for(ir=-1;ir<=nr;ir++) (*gs)[ir]   = 0.;
        for(ir=-1;ir<=nr;ir++) (*edr)[ir]  = 1./((*rcor)[ir]);    
        for(ir=-1;ir<=nr;ir++) (*edrs)[ir] = 1./((*rs)[ir]);
        for(ir=-1;ir<=nr;ir++) (*vval)[ir] = 1.;
        for(ir=-1;ir<=nr;ir++) ( *vs )[ir] = 1.;
        ispolar = FALSE;
    }


    for(ir=-1;ir<=nr;ir++) (*vval)[ir] /= p->dy;
    for(ir=-1;ir<=nr;ir++) ( *vs )[ir] /= p->dy;

    for(ir=-1;ir<=nr;ir++) (*norm)[ir] = 1./fabs((*hval)[ir]*(*vval)[ir]);


    COMM(for(ir=-1;ir<=nr;ir++) fprintf(stderr,"Proc %d: norm[%d] = %f\n",
                                        data->this_process,ir,(*norm)[ir]););

    COMM(for(ir=-1;ir<=nr;ir++) fprintf(stderr,"Proc %d: hval[%d] = %f, vval[%d] = %f\n",
                                        data->this_process,ir,(*hval)[ir],ir,(*vval)[ir]););


    COMM(for(ir=-1;ir<=nr;ir++) fprintf(stderr,"%d %f %f %f %f \n",ir,(*rcor)[ir],(*hval)[ir],(*gval)[ir],(*vval)[ir]););
  


    return ispolar; 
}


// Remove leading whitespaces from string

void Util_TrimString(char *str)
{
    int i=0,j=0;
    // fprintf(stderr,"Trimming >%s<\n",str); 

    //while(str[i]!='\0')  
    //{
    //    fprintf(stderr,"%5d\t:>%c< \n",i,str[i]);
    //    i++;
    //    
    //}
    //i = 0;
    
    while((str[i]!='\0') && (str[i] == ' ' ||  str[i] == '\t')) i++;
    //fprintf(stderr,"Remove >%d< leading whitespaces\n",i); 

    while(str[i]!='\0') 
    {
        //    fprintf(stderr,"%5d -> %5d\t: >%c< -> >%c<\n",i,j,str[j],str[i]); 
        str[j] = str[i];
        j++;
        i++;
    }
    //fprintf(stderr,"Close string at  >%d< \n",j);   

str[j]='\0';
//fprintf(stderr,"Trimmed >%s<\n",str); 

} 
