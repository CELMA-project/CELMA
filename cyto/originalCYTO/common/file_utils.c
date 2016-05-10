/* 
  V. Naulin 
 
*/

 
#undef DEBUG 
#undef COMM_STATUS

#include <utilities.h>
#define DMATRIX Util_DMatrix
#define  JUST_WRITE_INI_FILE  1 
#define  READ_INI_FILE  0



/*****************************************************************************/
/*! \brief Parses commandline 
 *
 * This function reads in a HDF datafile including parameters
 * it fills the structures PARA and HDF_DS with the parameters
 * but will overwrite values with the ones given as optional arguments 
 * in the commandline. 
 * It adds the run to the databasefile and arranges a run numb
  er
 * 
 * ON output: It returns a flag indicating if the file read is an initial file, we start from an .ini file 
 * or if we restart a run
 *
 * Commandline parameters:
 *
 *  -F|-R|-B <name> <number>:              name of HDF file in form name.number, 
 *                                     options are  F :normal, R: Restart, B: Batch (not interactive)
 *  -I <inifilename>	    :              start directly from .ini file
 *  -P <inifilename>        :              calc secondary parameters from physical ones and write a new .ini file.
 *  -N <name <number>	    :              additional argument if new outputfilename to be specified
 *  -H 			            :              short help text and a sample .ini file
 *  -D 			            :              Run with buildin parameters
 *  -<name> <value>	     	:              change parameter <name> to <value>
 */ 


int FUtils_ReadArguments(int argc,char **argv,HDF_DS *data,PARA *para, void (*CalcParameters)(PARA *,HDF_DS *) )
{
    char	
        namelong[DEFSTRLEN] = "\0",
        parname[DEFSTRLEN]  = "\0",
        name2[DEFSTRLEN]    = "\0",
        ininame[DEFSTRLEN]  = "\0",
        *pname = NULL;

    int	
        numberin      = -1,
        i             = 0,
        restart       = START_FROM_FILE,
        PARA_FROM_INI = FALSE,
        verbose       = TRUE,
        NEW_FILE_NAME_REQUEST          = FALSE;

    double 
      dummy  = 0.;

  

    //	Read in a datafile and evaluate the commandline
    data->number   = -1;
 
    snprintf(data->name, DEFSTRLEN,"%s","EMPTY");
    snprintf(namelong,   DEFSTRLEN,"%s","EMPTY");
    snprintf(parname,    DEFSTRLEN,"%s","EMPTY");
    snprintf(name2,      DEFSTRLEN,"%s","EMPTY");
    snprintf(ininame,    DEFSTRLEN,"%s","EMPTY");
    snprintf(data->jobid,DEFSTRLEN,"%s",argv[0]);
    BUGREPORT;

 

    BUGREPORT;


/*
 * Evaluate commandline
 */

    i=0; 
    while(argv[++i] != NULL)
    {
        if(strcmp(argv[i],"-F") == 0) // Start from initial condition
        {
            if(argv[++i] == 0){USAGE;}
            strcpy(namelong,argv[i]);

            if(argv[++i] == 0) {USAGE;}
            sscanf(argv[i],"%i",&numberin);

            verbose = TRUE;
            restart = START_FROM_FILE;
        }
        else if(strcmp(argv[i],"-B") == 0) // Batch job. no questions asked 
        {
            if(argv[++i] == 0) {USAGE;}
            strcpy(namelong,argv[i]);

            if(argv[++i] == 0) {USAGE;}
            sscanf(argv[i],"%i",&numberin);

            verbose = FALSE;
            restart = START_FROM_FILE;

        }
        else if(strcmp(argv[i],"-R") == 0) // Restart from file 
        {
            if(argv[++i] == 0) {USAGE;}
            strcpy(namelong,argv[i]);

            if(argv[++i] == 0) {USAGE;}
            sscanf(argv[i],"%i",&numberin);

            fprintf(stderr,"Restarting from file %s.%03d.\n",namelong,numberin);

            verbose = FALSE;
            restart = RESTART;
        }
        else if(strcmp(argv[i],"-D") == 0) // Start from values in memory
        {    
            snprintf(namelong,DEFSTRLEN,"default");
            snprintf(ininame,DEFSTRLEN,"default");
            numberin = 0;
            fprintf(stderr,"Starting with default parameters.\n");
  
            PARA_FROM_INI = TRUE;  
            (*CalcParameters)(para,data);
            verbose = FALSE;            
            restart = DEFAULTSTART;

        }
        else if(strcmp(argv[i],"-N") == 0) // Give output files a distinct name 
        {
            if(argv[++i] == 0) {USAGE;}
            strcpy(name2,argv[i]);

            if(argv[++i] == 0) {USAGE;}
            sscanf(argv[i],"%d",&data->number_out);
  
            restart = START_FROM_FILE;
            verbose = FALSE;
            NEW_FILE_NAME_REQUEST = TRUE;
        }
        else if(strcmp(argv[i],"-I") == 0) // Start with parameters from .ini  file 
        { 
            if(argv[++i] == 0) {USAGE;}
            strcpy(parname,argv[i]);

            numberin = 0;
            Util_ReadIniFile(data,para,parname,READ_INI_FILE);
            PARA_FROM_INI = TRUE;
            restart = START_FROM_INI;
            verbose = FALSE;         
        }
        else if(strcmp(argv[i],"-H") == 0) // Write a sample .ini file and exit  
        {
            snprintf(ininame,DEFSTRLEN,"%s_sample.ini",para->codename); 
            Util_ReadIniFile(data,para,ininame,JUST_WRITE_INI_FILE);
            USAGE;
        }
        else if(strcmp(argv[i],"-P") == 0) // Read .ini file, calculate parameters from primary physics and write a secondary .ini file  
        {

            if(argv[++i] == 0) {USAGE;}
            strcpy(parname,argv[i]);
            numberin = 0;
            Util_ReadIniFile(data,para,parname,READ_INI_FILE);
            BUGREPORT;

            /* Needs call to  a function depending on caller to 
               calculate new parameters!!! */
            (*CalcParameters)(para,data);
            BUGREPORT;
            snprintf(ininame,DEFSTRLEN,"%s_second.ini",para->codename); 
            BUGREPORT;
            Util_ReadIniFile(data,para,ininame,JUST_WRITE_INI_FILE);
            USAGE; // ends here
        }
    }
       


 
    if(numberin == -1 ) USAGE;
    BUGREPORT;

    // Read attributes from file
    data->number = numberin;
    snprintf(data->name_in,DEFSTRLEN,"%s",namelong);
    snprintf(data->name,   DEFSTRLEN,"%s",namelong);
    
    // deal with output filenumber
    if(NEW_FILE_NAME_REQUEST)
    {
        snprintf(data->name_out,DEFSTRLEN,"%s",name2);
        data->number_out--;   
    }
    else
    {
        snprintf(data->name_out,DEFSTRLEN,"%s",namelong);
        data->number_out = numberin;
    }
      
    snprintf(data->filename,DEFSTRLEN,"%s",data->name_out);
    BUGREPORT;


    //Read parameters from an existing HDF file  
    if(!PARA_FROM_INI) 
    {
        fprintf(stderr,"Reading attributes from %s.%03d\n",data->name,(int)data->number);
        data->ReadAttributes = TRUE;
        FUtilsInt_ReadHDF4Attributes(data->name_in,data->number,data,para);
        // This fills para structure 
    }
       
         
    /*
     *
     *  Read parameters from commandline
     * here as well the .ini file as the sample file should have been read!
     * commandline parameters overwite the information retrieved
     *
     */
      
     
    if(verbose)  printf("\n");
  
#define CHANGEPAR(a)	{sscanf(argv[i],"%lf", &(a) );if(TRUE) printf("Change of parameter >>%s<< to a new value of >>%e<<.\n",pname,(a)) ;}
  
    i=0; 
    while(argv[++i] != NULL)
    {
        if(strcmp(argv[i],"-F") == 0)   i += 2;		// test for old stuff 
        else if(strcmp(argv[i],"-R") == 0) i += 2;
        else if(strcmp(argv[i],"-B") == 0) i += 2;
        else if(strcmp(argv[i],"-N") == 0) i += 2;
        else if(strcmp(argv[i],"-I") == 0) i += 1;
        else if(strcmp(argv[i],"-D") == 0) i += 0;
        else if(argv[i][0] != '-') {USAGE;}
        else
        {
            pname = argv[i]+1;
            if(argv[++i] == NULL) USAGE;
            if     (strcmp(pname,"time")    == 0) {CHANGEPAR(para->time); }
            else if(strcmp(pname,"dt")	    == 0) {CHANGEPAR(para->dt);   }
            else if(strcmp(pname,"xmax")    == 0) {CHANGEPAR(para->xmax); }
            else if(strcmp(pname,"xmin")    == 0) {CHANGEPAR(para->xmin); }
            else if(strcmp(pname,"ymax")    == 0) {CHANGEPAR(para->ymax); }
            else if(strcmp(pname,"ymin")    == 0) {CHANGEPAR(para->ymin); }
            else if(strcmp(pname,"zmax")    == 0) {CHANGEPAR(para->zmax); }
            else if(strcmp(pname,"zmin")    == 0) {CHANGEPAR(para->zmin); }
            else if(strcmp(pname,"mue_n")	== 0) {CHANGEPAR(para->mue_n);}
            else if(strcmp(pname,"mue_w")	== 0) {CHANGEPAR(para->mue_w);}
            else if(strcmp(pname,"mue_t")	== 0) {CHANGEPAR(para->mue_t);}
            else if(strcmp(pname,"nu")	    == 0) {CHANGEPAR(para->nu);   }
            else if(strcmp(pname,"sigma")   == 0) {CHANGEPAR(para->sigma);}
            else if(strcmp(pname,"kappan")  == 0) {CHANGEPAR(para->kappan);}
            else if(strcmp(pname,"kappat")  == 0) {CHANGEPAR(para->kappat);}
            else if(strcmp(pname,"source")  == 0) {CHANGEPAR(para->source);}
            else if(strcmp(pname,"limiter") == 0) {CHANGEPAR(para->limiter);}
            else if(strcmp(pname,"dt_pol")  == 0) {CHANGEPAR(para->dt_pol);}
            else if(strcmp(pname,"alpha")   == 0) {CHANGEPAR(para->alpha);}
            else if(strcmp(pname,"beta")    == 0) {CHANGEPAR(para->beta);}
            else if(strcmp(pname,"delta")   == 0) {CHANGEPAR(para->delta);}
            else if(strcmp(pname,"gamma")   == 0) {CHANGEPAR(para->gamma);}
            else if(strcmp(pname,"hm_nl")   == 0) {CHANGEPAR(para->hm_nl);}
            else if(strcmp(pname,"nprof")   == 0) {CHANGEPAR(para->nprof);}
            else if(strcmp(pname,"bprof")   == 0) {CHANGEPAR(para->bprof);}
            else if(strcmp(pname,"phiprof") == 0) {CHANGEPAR(para->phiprof);}
            else if(strcmp(pname,"tprof")   == 0) {CHANGEPAR(para->tprof);}
            else if(strcmp(pname,"r0")      == 0) {CHANGEPAR(para->r0);}
            else if(strcmp(pname,"k0")      == 0) {CHANGEPAR(para->k0);}
            else if(strcmp(pname,"exb_ll")  == 0) {CHANGEPAR(para->exb_ll);}
            else if(strcmp(pname,"nl_p")    == 0) {CHANGEPAR(para->exb_ll);}

            //Process Grid 

            else if(strcmp(pname,"NX")    == 0) {CHANGEPAR(dummy); 
                data->N[2] = (int)dummy; 
                if(data->N[0]== 0) data->N[0] =1;
            }
            else if(strcmp(pname,"NY")    == 0) {CHANGEPAR(dummy); data->N[1] = (int)dummy;}
            else if(strcmp(pname,"NZ")    == 0) {CHANGEPAR(dummy); data->N[0] = (int)dummy;  
                if(data->N[2]== 0) data->N[2] =1;}
   
            /*Primary Physics parameters */
           
            else if(strcmp(pname,"Ti")      == 0) {CHANGEPAR(para->Ti);}
            else if(strcmp(pname,"Te")      == 0) {CHANGEPAR(para->Te);}
            else if(strcmp(pname,"B0")      == 0) {CHANGEPAR(para->B0);}
            else if(strcmp(pname,"n0")      == 0) {CHANGEPAR(para->n0);}
            else if(strcmp(pname,"Mi")   	== 0) {CHANGEPAR(para->Mi);}
            else if(strcmp(pname,"Z")	    == 0) {CHANGEPAR(para->Z);}
            else if(strcmp(pname,"p_n")	    == 0) {CHANGEPAR(para->p_n);}

            /* Primary geometry parameters */ 

            else if(strcmp(pname,"R0")      == 0) {CHANGEPAR(para->R0);}
            else if(strcmp(pname,"a")       == 0) {CHANGEPAR(para->a);}
            else if(strcmp(pname,"shat")    == 0) {CHANGEPAR(para->shat0);}
            else if(strcmp(pname,"q0")      == 0) {CHANGEPAR(para->q0);}
            else if(strcmp(pname,"lpar")    == 0) {CHANGEPAR(para->lpar);}
            else if(strcmp(pname,"ln")      == 0) {CHANGEPAR(para->ln);}
            else if(strcmp(pname,"lTe")     == 0) {CHANGEPAR(para->lTe);}
            else if(strcmp(pname,"lTi")     == 0) {CHANGEPAR(para->lTi);}
 

            /* Time integration parameters */
            else if(strcmp(pname,"et")	    == 0) {CHANGEPAR(para->end_time);}
            else if(strcmp(pname,"ot")	    == 0) {CHANGEPAR(para->out_time);}
            else if(strcmp(pname,"otmult")  == 0) {CHANGEPAR(dummy); para->otmult = (int)dummy;}
            else if(strcmp(pname,"xbnd")    == 0) {CHANGEPAR(dummy); para->xbnd = (int)dummy;}
            else if(strcmp(pname,"ybnd")    == 0) {CHANGEPAR(dummy); para->ybnd = (int)dummy;}
            else if(strcmp(pname,"zbnd")    == 0) {CHANGEPAR(dummy); para->zbnd = (int)dummy;}
            else if(strcmp(pname,"wbdra")   == 0) {CHANGEPAR(para->bdval[0][0]);}
            else if(strcmp(pname,"wbdrb")   == 0) {CHANGEPAR(para->bdval[0][1]);}
            else if(strcmp(pname,"nbdra")   == 0) {CHANGEPAR(para->bdval[2][0]);}
            else if(strcmp(pname,"nbdrb")   == 0) {CHANGEPAR(para->bdval[2][1]);}
            else if(strcmp(pname,"fbdra")   == 0) {CHANGEPAR(para->bdval[1][0]);}
            else if(strcmp(pname,"fbdrb")   == 0) {CHANGEPAR(para->bdval[1][1]);}
            else if(strcmp(pname,"adrhos")  == 0) {CHANGEPAR(para->adrhos);}
            else if(strcmp(pname,"r_offset")  == 0) {CHANGEPAR(dummy); para->r_offset = (int)dummy;}
            else if(strcmp(pname,"z_offset")  == 0) {CHANGEPAR(dummy); para->z_offset = (int)dummy;}
            else if(strcmp(pname,"muehat")  == 0) {CHANGEPAR(para->muehat);}	
            else 	      
            { 
                {
                    fprintf(stderr,"Online change of parameter >%s< not supported. Ignoring.\n",pname);
                    fprintf(stderr,"Known parameters are f.x.:\n");
                    fprintf(stderr,"et,      ot,      dt,      otmult\n");
                    fprintf(stderr,"nu,      beta,    sigma,   delta,\n");
                    fprintf(stderr,"hm_nl,   nl_phi,  exb_ll,  nl_p\n");
                    fprintf(stderr,"mue_w,   mue_n,   mue_t,\n");
                    fprintf(stderr,"shat0,   q0,      source,  limiter,\n");
                    fprintf(stderr,"nprof,   bprof,   tprof,   phiprof, r0,\n");
                    fprintf(stderr,"xmax,    xmin,    ymax,    ymin,   zmax,     zmin\n");
                    fprintf(stderr,"shat0,   q0,      source,  limiter,\n");
                    fprintf(stderr,"dt_pol,  nprof,   bprof,   tprof,   phiprof, r0,\n");
                    fprintf(stderr,"kappan,  kappat,  xbnd,    ybnd,    zbndn\n\n" );
                    fprintf(stderr,"\n");
                    fprintf(stderr,"The usefulness of such changes is down to the user.\n");

                }  
            }
        }
    }


    // Here all parameters are now known, in case axes have changed we acknowledge this here
    FUtilsInt_SetupDataRange(data,para);

    // Calculates Dx, Dy, Allocates Coordinates If Not Done Yet 
    FUtils_SetupParaStructure(data,para);

    // Handles interaction with run database
    FUtilsInt_AddToRunDatabase(data,para,argv);

    // Make runname and inform user
    fprintf(stderr,"Starting run No. %04d. %s.\n",(int)data->run_no,data->jobid);
    snprintf(data->name_out,DEFSTRLEN,"%s.%d",para->codename,(int)data->run_no);
    

    return (restart);
}


/***************************************************************************/
/*! Packs 2d field in datarw structure for writing, converts to float at the same time */
void FUtils_Write2dSpace(double **a,const char *fname,const char *name,long num, HDF_DS *data,PARA *p,int create)
{  
    int iy,ix;
    float *element;    
    
    data->create = create;
    FUtils_AllocateReadMemory(data);
    

    element = (float*)data->datarw;
    BUGREPORT;
    
    for(iy=0;iy<data->elements[0];iy++) 
        for(ix=0;ix<data->elements[1];ix++,element++)
            *element = (float)a[iy][ix];



    BUGREPORT;
    FUtilsInt_WriteHDF4(name,num,data,p,fname);
    data->create = FALSE;
}

/***************************************************************************/

/*! Packs 3d field in datarw structure for writing, converts to float at the same time */
 
void FUtils_Write3dSpace(double ***a,const char *fname, const char *name,long num,HDF_DS *data,PARA *p,int create)
{
  int iz,iy,ix,l=0;
    float *fltp=NULL;    

    data->create = create;
    FUtils_AllocateReadMemory(data);

    fltp = (float*)data->datarw;

    for(iz=0;iz<data->elements[0];iz++) 
        for(iy=0;iy<data->elements[1];iy++) 
            for(ix=0;ix<data->elements[2];ix++,l++)
                fltp[l] = (float)a[iz][iy][ix];

    FUtilsInt_WriteHDF4(name,num,data,p,fname);    
    data->create=FALSE;
}

/********************************************************************/
/*! Internal routine, writes netcdf, data are floats packed int datarw field  */
 
void FUtilsInt_WriteHDF4(const char *name,int32 number,HDF_DS *data,PARA *para,const char *feldname)
{
    int 
      i;
    char 
        filename[DEFSTRLEN];
 
    int32 
        sd_id     = -1,
        sds_id    = -1,
        dim_id,
        d_type    = DFNT_FLOAT32,
        icreate;
  
 
 
    BUGREPORT;

    if(number < 0) snprintf(filename,DEFSTRLEN,"%s",name);
    else snprintf(filename,DEFSTRLEN,"%s.%03d",name,(int)number);

    // Set data->range from para
    FUtilsInt_SetupDataRange(data,para);

    COMM(fprintf(stderr,"Write >%s<. Data-create = >%d< (%d == create)\n",filename,data->create,TRUE););

    // Set access to HDF file
    if(data->create == NEW_HDF) icreate = DFACC_CREATE;
    else icreate = DFACC_RDWR;
    
    sd_id = SDstart(filename,icreate);
    BUGREPORT;

    if(sd_id == -1)
    {
        fprintf(stderr,"%s: Process %d: SDstart failed on >%s< with create = %d\n",
                __func__,(int)data->this_process,filename,(int)data->create);
        return;
    }
    

    //Now the SDS 
 
    COMM(fprintf(stderr,"Write rank %d,sd_id = %d\n",(int)data->rank,(int)sd_id););  
    BUGREPORT;
    
    // Check if sd_id already contains field fieldname
    if ( (sds_id = SDselect(sd_id,SDnametoindex(sd_id,feldname))) == -1)
    {
        sds_id = SDcreate(sd_id,feldname,(int32)d_type,(int32)data->rank,data->dims);
        COMM(fprintf(stderr,"created new field: %s id=%d\n",feldname,(int)sds_id););
    }
 
 
    if(sds_id == -1)

    {
        fprintf(stderr,"%s:%d Process %d: SDcreate/SDselect failed on >%s<\n",
                __func__,__LINE__,(int)data->this_process,filename);
        fprintf(stderr,"%s: d_type %d (%d), nx = %d ny = %d\n",
                __func__,(int)d_type,(int)DFNT_FLOAT32,(int)data->dims[0],(int)data->dims[1]);
        fprintf(stderr,"%s: sds_id = %d sd_id = %d, data_name = >%s<, rank = %d\n",
                __func__,(int)sds_id,(int)sd_id,feldname,(int)data->rank);
        return;
    }
   
    assert(SDwritedata(sds_id,data->start,NULL,data->elements,(VOIDP)data->datarw)+1);
    BUGREPORT;


    /* Attributes for SDS */
    COMM(fprintf(stderr,"cordsys = >%s< \n (IMPORTANT: name should not equal any fieldname!!!!!!!)\n",data->coordsys););
    assert(SDsetattr(sds_id,"cordsys",DFNT_CHAR8,(int32)strlen(data->coordsys),(VOIDP)data->coordsys)+1);
 
    /* Coordinates with attributes */
    BUGREPORT;	  	  
    for(i = 0;i<data->rank;i++)
    {
        dim_id = SDgetdimid(sds_id,i);
        BUGREPORT;
        COMM(fprintf(stderr,"dim_is %d coord is %d label is >%s<\n",(int)dim_id,i,data->dim_label[i]););
        assert(SDsetdimname(dim_id,data->dim_label[i])+1); 
        BUGREPORT;

        
        COMM(fprintf(stderr,"dim_id: >%d< coord is %d dims is %d\n",(int)dim_id,i,(int)data->dims[i]););
     
        if(data->dims[0] > 0)
        {
            if(data->coordinate[i] != NULL)
            {
                assert(SDsetdimscale(dim_id,(int32)data->dims[i],DFNT_FLOAT64,(VOIDP)data->coordinate[i])+1);
            }
            else
            {
                fprintf(stderr,"coordinate[%d] not defined !\n",i);
                exit(-1);
            }
            
        }
        
        BUGREPORT;       
        assert(SDsetattr(dim_id,"range",DFNT_FLOAT64,2,(VOIDP)&data->range[i][0])+1);
     
        COMM(printf("Set Range %f %f\n",data->range[i][0],data->range[i][1]););    
        BUGREPORT;
    }


    if(SDendaccess(sds_id) == FAIL)
    {    
        fprintf(stderr,"%s %s %d: end sds access failed \n",__FILE__,__func__,__LINE__);
    }
     
    if(SDend(sd_id)==FAIL)
    { 
      fprintf(stderr,"%s %s %d: end sd access failed \n",__FILE__,__func__,__LINE__);
    }

    COMM(fprintf(stderr,"%s %s %d: end sds and sd access\n",__FILE__,__func__,__LINE__));

    if(data->create)  
    {
        FUtilsInt_WriteHDF4Attributes(name,number,data,para);
        COMM(fprintf(stderr,"wrote attributes\n\n"));
    }
    
}

/********************************************************************/


/* Opens a HDF File and checks if it exists, uses different ways to determine a filename from a number/
   name combination 
   filename is an output */

int FUtilsInt_MakeHDFFilename(char *filename, const char *name,int32 number)
{
    FILE *pf;
 
    if(number >= 0) 
        snprintf(filename,DEFSTRLEN, "%s.%.03d", name,(int)number);
    else 
        snprintf(filename,DEFSTRLEN,"%s",name);

  
    if (!(pf = fopen(filename, "r"))) 
    {
        fprintf(stderr,"%s %s %d: Can't open data file: %s\n",__FILE__,__func__,__LINE__,filename);
        filename = NULL;
        return(-1);
    } 
    fclose(pf);
 
    if(DFishdf(filename)!=0)
    {
        fprintf(stderr,"%s %s %d: %s ist kein HDF-File !",__FILE__,__func__,__LINE__,filename);
        filename = NULL;
        return(-1);
    }

    return 0;
}

/********************************************************************/

int FUtilsInt_ReadFileInfo(const char *name,int32 number,HDF_DS *data,char *filename)
{
    register int 
        i;
  
    int 
        num_2d_sds=0,
        num_3d_sds=0,
        num_coord=0;

    int32 
        rank = 0,    
        sd_id = -1 ,   
        sds_id = -1,
        k,num_type,
        attributes;
  
    
    if (FUtilsInt_MakeHDFFilename(filename,name,number) != 0) return -1;

    if( (sd_id = SDstart(filename,DFACC_RDONLY)) == -1)
    {
        fprintf(stderr,"%s %s %d: Error on SDstart: >%s<\n",__FILE__,__func__,__LINE__,filename);
        return(-1); 
    }
  
    assert(SDfileinfo(sd_id,&data->data_sets,&data->file_attrs)+1);
  
    COMM(printf("Datasets = %d, FileAttr = %d \n",(int)data->data_sets,(int)data->file_attrs);)  ;
  
    k = 0;
  
    /* Check for number of 2D and 3D data sets */

    for(k=0; k <  data->data_sets;k++)
    {
      COMM(fprintf(stderr,"Probing %d out of %d datasets\n",(int)k,(int)data->data_sets);)  ;
        if((sds_id = SDselect(sd_id,k)) == -1)
        {
            fprintf(stderr,"%s: Error on SDselect: %d, %d\n",__func__,(int)sd_id,(int)k);
            SDend(sd_id);
            return(-1);
        } 
      
        if(!SDiscoordvar(sds_id)) 
        {
            SDgetinfo(sds_id,&data->name[0],&rank,data->dims,&num_type,&attributes);
            switch (rank)
            {
                case 3:
                    num_3d_sds++;
                    COMM(fprintf(stderr,"Dataset No. %d: Rank %d\n",(int)k,(int)rank););
                    break;
                case 2:
                    num_2d_sds++;
                    COMM(fprintf(stderr,"Dataset No. %d: Rank %d\n",(int)k,(int)rank););
                    break;
                default:
                    fprintf(stderr,"Dataset No. %d: Rank %d not supported\n",(int)k,(int)rank);
            }
            SDread_string(sds_id,"cordsys",data->coordsys,DEFSTRLEN-1,"cartesian"); 
        }
        else
        {
            COMM(fprintf(stderr,"Dataset No. %d is a COORDINATE\n",(int)k););
            num_coord++;
        }
        SDendaccess(sds_id);
    }
    BUGREPORT;
 
    SDend(sd_id);   
    BUGREPORT;      
  
    if(num_2d_sds == 0) 
    {
        data->rank = 3;
    }
    else
    {
        data->rank = 2;
    }
    
    COMM(fprintf(stderr,"Found %d 2D and %d 3D datasets, %d datasets are Coordinates \n", num_2d_sds,num_3d_sds,num_coord););
    COMM(for(i=0;i<data->rank;i++) 
	   fprintf(stderr,"Dim No. %d has  %d elements.\n",i,(int)data->dims[i]););


     for(i=0;i<data->rank;i++)
        if(data->elements[i]  == 0) 
        {
            data->start[i] = 0;
            data->end[i] = data->dims[i];
            data->elements[i] = data->end[i]-data->start[i];
        }
    BUGREPORT;
    return(0);
}


/****************************************************************************************************/

int FUtilsInt_ReadHDF4ByNumber( char *name,int32 number,int numsds, HDF_DS *data,PARA *para)
{
    register int 
        i;
    int 
        noe=0;

    int32
        rank=0,
        sd_id,
        sds_id = -1,
        k,
        num_type,
        attributes,
        found=0,
        ISCOORDINATE = TRUE;

    char	
        filename[DEFSTRLEN];
  
 
    BUGREPORT; 
    if (FUtilsInt_ReadFileInfo(name,number,data,filename) != 0) return -1;



    if(data->rank == 0 )
    {
        fprintf(stderr,"Error: Rank not set in main\n");
        return(-1); 
    }


    if( (sd_id = SDstart(filename,DFACC_RDONLY)) == -1)
    {
        fprintf(stderr,"Error on SDstart: >%s<\n",filename);
        return(-1); 
    }
  
    COMM(fprintf(stderr,"Reading data.field %d have %d\n", (int)numsds,(int)data->data_sets););
    if(data->read_data == FALSE) return(0);
    BUGREPORT;

    k=0;


    
  
    for(found = 0;found < numsds ;found++)
    {
        ISCOORDINATE = TRUE;
      
        // As long as we have co-ordinates or wrong dimension or wrong data-format: continue reading 
        while( ISCOORDINATE || rank != data->rank || 
               ( (num_type != DFNT_FLOAT32) && (num_type != DFNT_FLOAT64))) 
        {
            if((sds_id = SDselect(sd_id,k)) == -1) 
                printf("Error on SDselect: %d, %d,%d \n",(int)sd_id,(int)k,(int)numsds); 
	  
            /* Check if coordinate information */
	  
            ISCOORDINATE = SDiscoordvar(sds_id);
            if(!ISCOORDINATE) SDgetinfo(sds_id,&data->name[0],&rank,data->dims,&num_type,&attributes);

            k++;

            if(k > data->data_sets)
            {
                fprintf(stderr,"file_utils.c: Search for dataset %d ( %dD) failed\n",(int)found+1, (int)data->rank);
                SDendaccess(sds_id);
                SDend(sd_id);
                return(-1);
            } 
        }
    }



    noe=1;
    for(i=0;i<data->rank;i++) noe*=(data->elements[i]);

    if(data->datarw != NULL && data->datarw_size < 2*noe) {free(data->datarw);data->datarw_size = 0;}
    if( (data->datarw = (void*) calloc(2*noe,sizeof(double)))  == NULL) 
    {
        fprintf(stderr,"%s %s %d: could not allocate read field....\n",__FILE__,__func__,__LINE__);
        return(-1);
    }
    data->datarw_size = 2*noe; // size of doubles 

    // Workspace twice as big as needed
  
    COMM(
        for(i=0;i<data->rank;i++)
            fprintf(stderr,"Proc %d Field %d: Read field %s. Dim No%d in [%d,%d]; %d Elements (max %d)\n",
                    data->this_process,(int)found,data->name,(int)i,(int)data->start[i],(int)data->end[i],(int)data->elements[i],(int)data->dims[i]);
        );
  
    /* Read the Stuff */
    BUGREPORT;

    if(SDreaddata(sds_id,data->start,NULL,data->elements,(void*)data->datarw) == -1)
    {
      fprintf(stderr,"SDreaddata nicht erfolgreich in %s,process %d\n",filename,(int)data->this_process);
        SDendaccess(sds_id);
        SDend(sd_id);
        return(-1);
    }
    BUGREPORT;

    SDendaccess(sds_id);
    SDend(sd_id);


    BUGREPORT;
    if (data->ReadAttributes) 
    {
        BUGREPORT;
        FUtilsInt_ReadHDF4Attributes(name,number,data,para);
    }
    
    BUGREPORT;
    return(num_type);
} 

/****************************************************************************************************/

/* Setup of some para values after data structure is known */
/* sets xmin, xmax from data->range */

void FUtils_SetupParaStructure(HDF_DS *data, PARA *para)
{
    double 
        val             = 0.;
    
    int 
        i,
        SET_COORDINATES = FALSE;
    
    
    for(i=0;i<data->rank;i++)
        if(data->elements[i] < 1) data->elements[i] = data->dims[i];
    
    switch(data->rank)
    {
        case 2:
            para->xmax = data->range[1][1];
            para->xmin = data->range[1][0];
            para->ymax = data->range[0][1];
            para->ymin = data->range[0][0];
            para->zmax = 1.;
            para->zmin = 0.;

            data->nx = data->dims[1];
            data->ny = data->dims[0];
            data->nz = 1;


            data->lnx = data->elements[1];
            data->lny = data->elements[0];
            data->lnz = 0;
            break;
        case 3:
            para->xmax = data->range[2][1];
            para->xmin = data->range[2][0];
            para->ymax = data->range[1][1];
            para->ymin = data->range[1][0];
            para->zmax = data->range[0][1];
            para->zmin = data->range[0][0];
    
            data->nx = data->dims[2];
            data->ny = data->dims[1];
            data->nz = data->dims[0];
            
            data->lnx = data->elements[2];
            data->lny = data->elements[1];
            data->lnz = data->elements[0];
            break;
    }

    para->dx = (para->xmax-para->xmin)/(double)data->nx;
    para->dy = (para->ymax-para->ymin)/(double)data->ny;   
    para->dz = (para->zmax-para->zmin)/(double)data->nz;

    para->dkx = 2.0*M_PI/(para->xmax-para->xmin);
    para->dky = 2.0*M_PI/(para->ymax-para->ymin);
    para->dkz = 2.0*M_PI/(para->zmax-para->zmin);


      
    for(i=0;i<data->rank;i++) 
        if(data->coordinate[i] == NULL && data->dims[i] > 0 ) 
        {
            COMM(fprintf(stderr,"%s %s %d: Do allocate coordinate %d dim %d\n",__FILE__,__func__,__LINE__, i,(int)data->dims[i]););
            data->coordinate[i] = (double *)malloc(sizeof(double)*data->dims[i]);
            SET_COORDINATES = TRUE;//set equidistant coordinates 
            
        }

    if(SET_COORDINATES)
    {
        fprintf(stderr,"Setting up equidistant coordinates.\n");
        /* 
           Set coordinate points to equidistant default and name coordinates, 
           also in cases which do not make real sense 
        */

        if(data->rank == 3)
        {
            
            for(i=0,val = para->zmin+para->dz*.5;i<data->dims[0];i++,val+=para->dz)  data->coordinate[0][i] = val;
            for(i=0,val = para->ymin+para->dy*.5;i<data->dims[1];i++,val+=para->dy)  data->coordinate[1][i] = val;
            for(i=0,val = para->xmin+para->dx*.5;i<data->dims[2];i++,val+=para->dx)  data->coordinate[2][i] = val;

            switch(para->coordsys){
                case CARTESIAN:	
                    snprintf(data->coordsys,DEFSTRLEN,"%s","cartesian, equidistant");
                    snprintf(data->dim_label[0],DEFSTRLEN,"%s","z");
                    snprintf(data->dim_label[1],DEFSTRLEN,"%s","y");
                    snprintf(data->dim_label[2],DEFSTRLEN,"%s","x");
                break;
                case POLOIDAL:
                case CYLINDRICAL:	
                    snprintf(data->coordsys,DEFSTRLEN,"%s","cylindrical, equidistant");
                    snprintf(data->dim_label[0],DEFSTRLEN,"%s","z");
                    snprintf(data->dim_label[1],DEFSTRLEN,"%s","phi");
                    snprintf(data->dim_label[2],DEFSTRLEN,"%s","r");
                    break;
                case SPHERICAL:	
                    snprintf(data->coordsys,DEFSTRLEN,"%s","spherical, equidistant");
                    snprintf(data->dim_label[0],DEFSTRLEN,"%s","theta");
                    snprintf(data->dim_label[1],DEFSTRLEN,"%s","phi");
                    snprintf(data->dim_label[2],DEFSTRLEN,"%s","r");
                    break;
                case CART_TB_X:	
                    snprintf(data->coordsys,DEFSTRLEN,"%s","cartesian, x-Tschebischeff");
                    snprintf(data->dim_label[0],DEFSTRLEN,"%s","z");
                    snprintf(data->dim_label[1],DEFSTRLEN,"%s","y");
                    snprintf(data->dim_label[2],DEFSTRLEN,"%s","x");
                    fprintf(stderr,"ERROR:Tschebischeff is not implemented .\n");
                    break;
                case POL_TB_R:
                case CYL_TB_R:	
                    snprintf(data->coordsys,DEFSTRLEN,"%s","cylindrical, r-Tschebischeff");
                    snprintf(data->dim_label[0],DEFSTRLEN,"%s","z");
                    snprintf(data->dim_label[1],DEFSTRLEN,"%s","phi");
                    snprintf(data->dim_label[2],DEFSTRLEN,"%s","r");
                    fprintf(stderr,"ERROR:Tschebischeff is not implemented .\n");
                    break;
                case SPHER_TB_R:	
                    snprintf(data->coordsys,DEFSTRLEN,"%s","spherical, r-Tschebischeff");
                    snprintf(data->dim_label[0],DEFSTRLEN,"%s","theta");
                    snprintf(data->dim_label[1],DEFSTRLEN,"%s","phi");
                    snprintf(data->dim_label[2],DEFSTRLEN,"%s","r");
                    fprintf(stderr,"ERROR:Tschebischeff is not implemented .\n");
                    break;
                default: 
                    fprintf(stderr,"ERROR:No such coordinate in 3D.\n");
                    break;
            }
            
        }
        else if (data->rank == 2)
        {
            for(i=0,val = para->ymin+para->dy*.5;i<data->dims[0];i++,val+=para->dy)  data->coordinate[0][i] = val;
            for(i=0,val = para->xmin+para->dx*.5;i<data->dims[1];i++,val+=para->dx)  data->coordinate[1][i] = val;


            switch(para->coordsys){
                case CARTESIAN:	
                    snprintf(data->coordsys,DEFSTRLEN,"%s","cartesian, equidistant");
                    snprintf(data->dim_label[0],DEFSTRLEN,"%s","y");
                    snprintf(data->dim_label[1],DEFSTRLEN,"%s","x");
                    snprintf(data->dim_label[2],DEFSTRLEN,"%s","0");
                    break;
                case POLOIDAL:	
                    snprintf(data->coordsys,DEFSTRLEN,"%s","polar, equidistant");
                    snprintf(data->dim_label[0],DEFSTRLEN,"%s","phi");
                    snprintf(data->dim_label[1],DEFSTRLEN,"%s","r");
                    snprintf(data->dim_label[2],DEFSTRLEN,"%s","0");
                    break;
                case CART_TB_X:	
                    snprintf(data->coordsys,DEFSTRLEN,"%s","cartesian, x-Tschebischeff");
                    snprintf(data->dim_label[0],DEFSTRLEN,"%s","y");
                    snprintf(data->dim_label[1],DEFSTRLEN,"%s","x");
                    snprintf(data->dim_label[2],DEFSTRLEN,"%s","0");
                    fprintf(stderr,"ERROR:Tschebischeff is not implemented .\n");
                    break;
                case POL_TB_R:	
                    snprintf(data->coordsys,DEFSTRLEN,"%s","polar, r-Tschebischeff");
                    snprintf(data->dim_label[0],DEFSTRLEN,"%s","phi");
                    snprintf(data->dim_label[1],DEFSTRLEN,"%s","r");
                    snprintf(data->dim_label[2],DEFSTRLEN,"%s","0");
                    fprintf(stderr,"ERROR:Tschebischeff is not implemented .\n");
                    break;
                default: 
                    fprintf(stderr,"ERROR:No such coordinate in 2D.\n");
                    break;
            }
        }
        BUGREPORT;
    }
}




/*********************************************************************************************/

/*! Allocate memory to hold fields and coordinates, returns number of elements per process on cuccess, -1 on error */

int FUtils_AllocateReadMemory(HDF_DS *data)
{    
    int i=0,noe=1;
    
    
    BUGREPORT;
    COMM(printf("proc: %d\n",data->this_process););
    
    for(i=0;i<data->rank;i++) noe*=MAX(1,(data->elements[i]));
    if(data->datarw != NULL && data->datarw_size < 2*noe ) 
    {
        free(data->datarw);
        data->datarw = NULL;
	data->datarw_size = 0; 
    }
    

    if(data->datarw == NULL || data->datarw_size < 2*noe )
    {
      if( (data->datarw = (void*) calloc(2*noe,sizeof(double)))  == NULL) 
        {
            fprintf(stderr,"%s %s %d: could not allocate read field....\n",__FILE__,__func__,__LINE__);
            return(-1);
        }  
      data->datarw_size = 2*noe;
      
      BUGREPORT;
      COMM(fprintf(stderr,"Allocated datarw \n");); 
      
      /* Make sure that coordinates are allocated */
      
      for(i=0;i<data->rank;i++) 
      { 
          if(data->coordinate[i] == NULL && data->dims[i] > 0 ) 
          {
              COMM(fprintf(stderr,"Do allocate coordinate %d dim %d\n", i,(int)data->dims[i]););
              data->coordinate[i] = (double *)malloc(sizeof(double)*data->dims[i]);
          }
          if(data->coordinate[i] == NULL && data->dims[i] > 0 ) {
              fprintf(stderr,"Could not allocate read coordinate %d\n", i);
              return(-1);
          }
      }
      COMM(fprintf(stderr,"Allocate ready\n"););
    }
    BUGREPORT;
    return(noe);
}


/*********************************************************************************************/


int FUtilsInt_WriteHDF4Attributes(const char *name,int32 number,HDF_DS *data,PARA *para)
{
    time_t  
        tp = {0};
    char 
        *desc=NULL;
    int32 
        sd_id = -1;
    char	
        filename[DEFSTRLEN];
    char 
      desc_string[DEFSTRLEN];
    int32
        intval=0;
    


    BUGREPORT;
    if(number >= 0) 
        snprintf(filename,DEFSTRLEN, "%s.%.03d", name,(int)number);
    else 
        snprintf(filename,DEFSTRLEN,"%s",name);
  
    if( (sd_id = SDstart(filename,DFACC_RDWR)) == -1) // open SDS for RW 
    {
        if( (sd_id = SDstart(filename,DFACC_CREATE)) == -1) // open SDS for Write 
        {
            fprintf(stderr,"%s %s %d:Error on SDstart: >%s<\n",__FILE__,__func__,__LINE__,filename);
            return(-1); 
        }
        
    }

    BUGREPORT;
    /* Attributes for SD */
    COMM(fprintf(stderr,"%s %s %d:  >%s< %d\n",__FILE__,__func__,__LINE__,data->jobid,(int)strlen(data->jobid)););
    assert(SDsetattr(sd_id,"Commandline",DFNT_CHAR8,(int32)strlen(data->jobid),(VOIDP)data->jobid)+1);

    BUGREPORT;
    COMM(fprintf(stderr,"  >%s< %d\n",data->integrator,(int)strlen(data->integrator)););
    assert(SDsetattr(sd_id,"Programm",DFNT_CHAR8,(int32)strlen(data->integrator),(VOIDP)data->integrator)+1);


    COMM(fprintf(stderr,"  >%s< %d\n",data->revision,(int)strlen(data->revision)););
    assert(SDsetattr(sd_id,"Revision",DFNT_CHAR8,(int32)strlen(data->revision),(VOIDP)data->revision)+1);

 
    COMM(fprintf(stderr," >%s< %d\n",data->maschine,(int)strlen(data->maschine)););
    assert(SDsetattr(sd_id,"Machine",DFNT_CHAR8,(int32)strlen(data->maschine),(VOIDP)data->maschine)+1);
    BUGREPORT;

    COMM(fprintf(stderr," >%s< %d\n",data->compile_date,(int)strlen(data->compile_date)););
    assert(SDsetattr(sd_id,"Compile_date",DFNT_CHAR8,(int32)strlen(data->compile_date),(VOIDP)data->compile_date)+1);

    
    if(time(&tp) != -1) 
        strftime(data->write_date,21,"%d.%m.%Y %H.%M.%S",localtime(&tp));
    else
        snprintf(data->write_date,DEFSTRLEN,"%s","no_date_available");

    COMM(fprintf(stderr," >%s< %d\n",data->write_date,(int)strlen(data->write_date)););
    assert(SDsetattr(sd_id,"data_write_date",DFNT_CHAR8,(int32)strlen(data->write_date),(VOIDP)data->write_date)+1);
    BUGREPORT;

    COMM(fprintf(stderr,"  >%s<\n",data->desc););
    assert(SDsetattr(sd_id,"description",DFNT_CHAR8,(int32)strlen(data->desc),(VOIDP)data->desc)+1);


  
    BUGREPORT;

    COMM(fprintf(stderr,"  %d\n",(int)data->run_no););
    intval = (int32)data->run_no;
    assert(SDsetattr(sd_id,"run_no",DFNT_INT32,1,(VOIDP)&intval)+1);
    
    BUGREPORT;


    /****************************************************************************/  
    assert(SDsetattr(sd_id,"time",    DFNT_FLOAT64,1,(VOIDP)&para->time)+1);
    assert(SDsetattr(sd_id,"dt",      DFNT_FLOAT64,1,(VOIDP)&para->dt)+1);
    assert(SDsetattr(sd_id,"end_time",DFNT_FLOAT64,1,(VOIDP)&para->end_time)+1);
    assert(SDsetattr(sd_id,"out_time",DFNT_FLOAT64,1,(VOIDP)&para->out_time)+1);




    intval = (int32)para->otmult;
    assert(SDsetattr(sd_id,"ot_mult",DFNT_INT32,1,(VOIDP)&intval)+1);
    intval = (int32)para->r_offset;
    assert(SDsetattr(sd_id,"r_offset",DFNT_INT32,1,(VOIDP)&intval)+1);
    intval = (int32)para->z_offset;
    assert(SDsetattr(sd_id,"z_offset",DFNT_INT32,1,(VOIDP)&intval)+1);
    
    intval = (int32)para->r_spacing;
    assert(SDsetattr(sd_id,"r_spacing",DFNT_INT32,1,(VOIDP)&intval)+1);
    intval = (int32)para->z_spacing;
    assert(SDsetattr(sd_id,"z_spacing",DFNT_INT32,1,(VOIDP)&intval)+1);




  
    BUGREPORT;
 
    /****************************************************************************/
    if(strcmp(desc=Util_GetDesc(para->desc,"mue_n","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"mue_n",DFNT_FLOAT64,1,(VOIDP)&para->mue_n)+1);
        assert(SDsetattr(sd_id,"D_mue_n",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,"%s %s  %d:>%f< >%s<\n",__FILE__,__func__,__LINE__,para->mue_n,desc););
    }
     

    if(strcmp(desc=Util_GetDesc(para->desc,"mue_w","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"mue_w",DFNT_FLOAT64,1,(VOIDP)&para->mue_w)+1);
        assert(SDsetattr(sd_id,"D_mue_w",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->mue_w,desc););
    }
 
    if(strcmp(desc=Util_GetDesc(para->desc,"mue_t","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"mue_t",DFNT_FLOAT64,1,(VOIDP)&para->mue_t)+1);
        assert(SDsetattr(sd_id,"D_mue_t",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->mue_t,desc););
    }

    BUGREPORT;
 
    if(strcmp(desc=Util_GetDesc(para->desc,"alpha","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"alpha",DFNT_FLOAT64,1,(VOIDP)&para->alpha)+1);
        assert(SDsetattr(sd_id,"D_alpha",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->alpha,desc););
    }
  
    if(strcmp(desc=Util_GetDesc(para->desc,"beta","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"beta",DFNT_FLOAT64,1,(VOIDP)&para->beta)+1);
        assert(SDsetattr(sd_id,"D_beta",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->beta,desc););
    }
    if(strcmp(desc=Util_GetDesc(para->desc,"delta","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"delta",DFNT_FLOAT64,1,(VOIDP)&para->delta)+1);
        assert(SDsetattr(sd_id,"D_delta",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->delta,desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"gamma","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"gamma",DFNT_FLOAT64,1,(VOIDP)&para->gamma)+1);
        assert(SDsetattr(sd_id,"D_gamma",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->gamma,desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"sigma","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"sigma",DFNT_FLOAT64,1,(VOIDP)&para->sigma)+1);
        assert(SDsetattr(sd_id,"D_sigma",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->sigma,desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"nu","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"nu",DFNT_FLOAT64,1,(VOIDP)&para->nu)+1);
        assert(SDsetattr(sd_id,"D_nu",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->nu,desc););
    }
  
    if(strcmp(desc=Util_GetDesc(para->desc,"muehat","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"muehat",DFNT_FLOAT64,1,(VOIDP)&para->muehat)+1);
        assert(SDsetattr(sd_id,"D_muehat",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s< \n",para->muehat,desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"betahat","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"betahat",DFNT_FLOAT64,1,(VOIDP)&para->betahat)+1);
        assert(SDsetattr(sd_id,"D_betahat",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->betahat,desc););
    }
  
    if(strcmp(desc=Util_GetDesc(para->desc,"hm_nl","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"hm_nl",DFNT_FLOAT64,1,(VOIDP)&para->hm_nl)+1);
        assert(SDsetattr(sd_id,"D_hm_nl",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr," >%f< >%s<\n",para->hm_nl,desc););
    }
    if(strcmp(desc=Util_GetDesc(para->desc,"exb_ll","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"exb_ll",DFNT_FLOAT64,1,(VOIDP)&para->exb_ll)+1);
        assert(SDsetattr(sd_id,"D_exb_ll",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->exb_ll,desc););
    }


    BUGREPORT;


    if(strcmp(desc=Util_GetDesc(para->desc,"source","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"source",DFNT_FLOAT64,1,(VOIDP)&para->source)+1);
        assert(SDsetattr(sd_id,"D_source",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->source,desc););
    }


    if(strcmp(desc=Util_GetDesc(para->desc,"shat0","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"shat0",DFNT_FLOAT64,1,(VOIDP)&para->shat0)+1);
        assert(SDsetattr(sd_id,"D_shat0",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->shat0,desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"q0","-",desc_string),"-"))
     {
         assert(SDsetattr(sd_id,"q0",DFNT_FLOAT64,1,(VOIDP)&para->q0)+1);
         assert(SDsetattr(sd_id,"D_q0",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
         COMM(fprintf(stderr,">%f< >%s<\n",para->q0,desc););
     }

    if(strcmp(desc=Util_GetDesc(para->desc,"limiter","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"limiter",DFNT_FLOAT64,1,(VOIDP)&para->limiter)+1);
        assert(SDsetattr(sd_id,"D_limiter",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->limiter,desc););
    }
    if(strcmp(desc=Util_GetDesc(para->desc,"dt_pol","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"dt_pol",DFNT_FLOAT64,1,(VOIDP)&para->dt_pol)+1);
        assert(SDsetattr(sd_id,"D_dt_pol",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->dt_pol,desc););
    }
 
    if(strcmp(desc=Util_GetDesc(para->desc,"r0","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"r0",DFNT_FLOAT64,1,(VOIDP)&para->r0)+1);
        assert(SDsetattr(sd_id,"D_r0",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->r0,desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"k0","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"k0",DFNT_FLOAT64,1,(VOIDP)&para->k0)+1); 
        assert(SDsetattr(sd_id,"D_k0",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->k0,desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"kappan","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"kappan",DFNT_FLOAT64,1,(VOIDP)&para->kappan)+1);
        assert(SDsetattr(sd_id,"D_kappan",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->kappan,desc););
    }
 
    if(strcmp(desc=Util_GetDesc(para->desc,"kappat","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"kappat",DFNT_FLOAT64,1,(VOIDP)&para->kappat)+1);
        assert(SDsetattr(sd_id,"D_kappat",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->kappat,desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"adrhos","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"adrhos",DFNT_FLOAT64,1,(VOIDP)&para->adrhos)+1);
        assert(SDsetattr(sd_id,"D_adrhos",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->adrhos,desc););
    }


    if(strcmp(desc=Util_GetDesc(para->desc,"qprof","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"qprof",DFNT_FLOAT64,data->nx,(VOIDP)&para->qprof[0])+1);
        assert(SDsetattr(sd_id,"D_qprof",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);   
        COMM(fprintf(stderr,">%f< >%s<\n",para->qprof[0],desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"shat","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"shat",DFNT_FLOAT64,data->nx,(VOIDP)&para->shat[0])+1);
        assert(SDsetattr(sd_id,"D_shat",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->shat[0],desc););
    }
 
    if(strcmp(desc=Util_GetDesc(para->desc,"phiprof","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"phiprof",DFNT_FLOAT64,1,(VOIDP)&para->phiprof)+1);
        assert(SDsetattr(sd_id,"D_phiprof",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr," >%f< >%s<\n",para->phiprof,desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"bprof","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"bprof",DFNT_FLOAT64,1,(VOIDP)&para->bprof)+1);
        assert(SDsetattr(sd_id,"D_bprof",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr," >%f< >%s<\n",para->bprof,desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"nprof","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"nprof",DFNT_FLOAT64,1,(VOIDP)&para->nprof)+1);
        assert(SDsetattr(sd_id,"D_nprof",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr," >%f< >%s<\n",para->nprof,desc););
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"tprof","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"tprof",DFNT_FLOAT64,1,(VOIDP)&para->tprof)+1);
        assert(SDsetattr(sd_id,"D_tprof",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }

    BUGREPORT;
    assert(SDsetattr(sd_id,"energy",DFNT_FLOAT64,1,(VOIDP)&para->energy)+1);
    BUGREPORT;
    
    assert(SDsetattr(sd_id,"vorticity",DFNT_FLOAT64,1,(VOIDP)&para->vorticity)+1);
    BUGREPORT;
    assert(SDsetattr(sd_id,"reynolds",DFNT_FLOAT64,1,(VOIDP)&para->reynolds)+1);
    BUGREPORT;
    assert(SDsetattr(sd_id,"prandel",DFNT_FLOAT64,1,(VOIDP)&para->prandel)+1);
    assert(SDsetattr(sd_id,"Nusselt",DFNT_FLOAT64,1,(VOIDP)&para->nusselt)+1);
 
    BUGREPORT;
 
    intval = (int32)para->xbnd;
    assert(SDsetattr(sd_id,"xbnd",DFNT_INT32,1,(VOIDP)&intval)+1);
    intval = (int32)para->ybnd;
    assert(SDsetattr(sd_id,"ybnd",DFNT_INT32,1,(VOIDP)&intval)+1);
    intval = (int32)para->zbnd;
    assert(SDsetattr(sd_id,"zbnd",DFNT_INT32,1,(VOIDP)&intval)+1);
    assert(SDsetattr(sd_id,"bndcnd",DFNT_INT32,3,(VOIDP)&para->boundary[0])+1);
 
    BUGREPORT;
 
    if(strcmp(desc=Util_GetDesc(para->desc,"uvortex","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"uvortex",DFNT_FLOAT64,1,(VOIDP)&para->uvortex)+1);
        assert(SDsetattr(sd_id,"D_uvortex",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"radius","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"radius",DFNT_FLOAT64,1,(VOIDP)&para->radius)+1);
        assert(SDsetattr(sd_id,"D_radius",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }
 
    if(strcmp(desc=Util_GetDesc(para->desc,"epsilon","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"epsilon",DFNT_FLOAT64,1,(VOIDP)&para->epsilon)+1);
        assert(SDsetattr(sd_id,"D_epsilon",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }

    if(strcmp(desc=Util_GetDesc(para->desc,"v_pol","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"v_pol",DFNT_FLOAT64,1,(VOIDP)&para->v_pol)+1);
        assert(SDsetattr(sd_id,"D_v_pol",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }


/* Write primary and secondary physical parameters */

   /* Primary units in [eV],[T],[1/m**3],[u] */

 if(strcmp(desc=Util_GetDesc(para->Prim_Phys,"T_i","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"T_i",DFNT_FLOAT64,1,(VOIDP)&para->Ti)+1);
        assert(SDsetattr(sd_id,"D_T_i",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->Ti,desc););
    }

 if(strcmp(desc=Util_GetDesc(para->Prim_Phys,"T_e","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"T_e",DFNT_FLOAT64,1,(VOIDP)&para->Te)+1);
        assert(SDsetattr(sd_id,"D_T_e",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->Te,desc););
    }

 if(strcmp(desc=Util_GetDesc(para->Prim_Phys,"B0","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"B0",DFNT_FLOAT64,1,(VOIDP)&para->B0)+1);
        assert(SDsetattr(sd_id,"D_B0",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->B0,desc););
    }

 if(strcmp(desc=Util_GetDesc(para->Prim_Phys,"n0","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"n0",DFNT_FLOAT64,1,(VOIDP)&para->n0)+1);
        assert(SDsetattr(sd_id,"D_n0",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->n0,desc););

    }
 if(strcmp(desc=Util_GetDesc(para->Prim_Phys,"Mi","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"Mi",DFNT_FLOAT64,1,(VOIDP)&para->Mi)+1);
        assert(SDsetattr(sd_id,"D_Mi",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
        COMM(fprintf(stderr,">%f< >%s<\n",para->Mi,desc););
    }
 if(strcmp(desc=Util_GetDesc(para->Prim_Phys,"Z","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"Z",DFNT_FLOAT64,1,(VOIDP)&para->Z)+1);
        assert(SDsetattr(sd_id,"D_Z",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
       COMM(fprintf(stderr,">%f< >%s<\n",para->Z,desc););
    }


    /* Secondary, derived  units */

     if(strcmp(desc=Util_GetDesc(para->Sec_Phys,"rho_s","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"rho_s",DFNT_FLOAT64,1,(VOIDP)&para->rho_s)+1);
        assert(SDsetattr(sd_id,"D_rho_s",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }

   if(strcmp(desc=Util_GetDesc(para->Sec_Phys,"rho_i","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"rho_i",DFNT_FLOAT64,1,(VOIDP)&para->rho_i)+1);
        assert(SDsetattr(sd_id,"D_rho_i",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }

   if(strcmp(desc=Util_GetDesc(para->Sec_Phys,"rho_e","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"rho_e",DFNT_FLOAT64,1,(VOIDP)&para->rho_e)+1);
        assert(SDsetattr(sd_id,"D_rho_e",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }


  if(strcmp(desc=Util_GetDesc(para->Sec_Phys,"omega_ci","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"omega_ci",DFNT_FLOAT64,1,(VOIDP)&para->omega_ci)+1);
        assert(SDsetattr(sd_id,"D_omega_ci",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }

      if(strcmp(desc=Util_GetDesc(para->Sec_Phys,"omega_ce","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"omega_ce",DFNT_FLOAT64,1,(VOIDP)&para->omega_ce)+1);
        assert(SDsetattr(sd_id,"D_omega_ce",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }

  if(strcmp(desc=Util_GetDesc(para->Sec_Phys,"omega_pi","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"omega_pi",DFNT_FLOAT64,1,(VOIDP)&para->omega_pi)+1);
        assert(SDsetattr(sd_id,"D_omega_pi",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }

   BUGREPORT;
 if(strcmp(desc=Util_GetDesc(para->Sec_Phys,"v_thi","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"v_thi",DFNT_FLOAT64,1,(VOIDP)&para->v_thi)+1);
        assert(SDsetattr(sd_id,"D_v_thi",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }
 if(strcmp(desc=Util_GetDesc(para->Sec_Phys,"v_the","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"v_the",DFNT_FLOAT64,1,(VOIDP)&para->v_the)+1);
        assert(SDsetattr(sd_id,"D_v_the",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }
 if(strcmp(desc=Util_GetDesc(para->Sec_Phys,"c_s","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"c_s",DFNT_FLOAT64,1,(VOIDP)&para->c_s)+1);
        assert(SDsetattr(sd_id,"D_c_s",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }
 if(strcmp(desc=Util_GetDesc(para->Sec_Phys,"v_alfven","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"v_alfven",DFNT_FLOAT64,1,(VOIDP)&para->v_alfven)+1);
        assert(SDsetattr(sd_id,"D_v_alfven",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }
  if(strcmp(desc=Util_GetDesc(para->Sec_Phys,"n_n","-",desc_string),"-"))
    {
        assert(SDsetattr(sd_id,"n_n",DFNT_FLOAT64,1,(VOIDP)&para->n_n)+1);
        assert(SDsetattr(sd_id,"D_n_n",DFNT_CHAR8,(int32)strlen(desc),(VOIDP)desc)+1);
    }
  
 SDend(sd_id);      
 return(0);
 
 
}  
/*********************************************************************************************/
  
void FUtilsInt_ReadHDF4Attributes(const char *name,int32 number,HDF_DS *data, PARA *para)
{
    double tmp;
    double tarray[10];
    register int i;
    int32 sd_id = -1;
    char filename[DEFSTRLEN];


    if(!data->ReadAttributes){
        COMM(fprintf(stderr,"Do not read attribs\n"););
        return;
    }
    BUGREPORT;
    COMM(fprintf(stderr,"Do read attribs\n"););

    if (FUtilsInt_ReadHDF4Coordinates(name,number,data,para) == -1)
    {
        fprintf(stderr,"%s %s %d: Warning: Could not read coordinates\n",__FILE__,__func__,__LINE__);
    }
    else
    {
        FUtils_SetupParaStructure(data,para);
    }
    

    if (FUtilsInt_MakeHDFFilename(filename,name,number) != 0) return;
    COMM(fprintf(stderr,"JOBID: >%s<\n",data->jobid););

  
    if( (sd_id = SDstart(filename,DFACC_RDONLY)) == -1)
    {
        fprintf(stderr,"%s %s %d:Error on SDstart: >%s<\n",__FILE__,__func__,__LINE__,filename);
        return; 
    }

    
    BUGREPORT;
    
    SDread_string(sd_id,"file_contents",data->jobid,DEFSTRLEN-1,"Empty");
    BUGREPORT;
    SDread_string(sd_id,"Machine",data->maschine,DEFSTRLEN-1,"Empty");
    BUGREPORT;
 
    SDread_string(sd_id,"Programm",data->integrator,DEFSTRLEN-1,"Empty");
    SDread_string(sd_id,"Revision",data->revision,DEFSTRLEN-1,"Empty");
    SDread_string(sd_id,"Compile_date",data->compile_date,DEFSTRLEN-1,"Empty");
    BUGREPORT;
  
    SDread_string(sd_id,"data_write_date",data->write_date,DEFSTRLEN-1,"Empty");

    BUGREPORT;
    SDread_string(sd_id,"description",data->desc,DEFSTRLEN-1,"Empty");
  

    BUGREPORT;
    /*
      SDread_variable(sd_id,"run_no",&tmp,1,0); data->run_no = (int) tmp;  
    */

    BUGREPORT;
    /* Read File attributes */
  
    SDread_variable(sd_id,"time",&para->time,1,0.);
    BUGREPORT;
    SDread_variable(sd_id,"dt",&para->dt,1,0.);
    SDread_variable(sd_id,"end_time",&para->end_time,1,0.);
    SDread_variable(sd_id,"out_time",&para->out_time,1,0.);
    SDread_variable(sd_id,"ot_mult",&tmp,1,0.); para->otmult = (int) tmp;
 
    SDread_variable(sd_id,"r_offset",&tmp,1,0.); para->r_offset =  (int) tmp;
    SDread_variable(sd_id,"z_offset",&tmp,1,0.); para->z_offset =  (int) tmp;
    SDread_variable(sd_id,"r_spacing",&tmp,1,0.); para->r_spacing =  (int) tmp;
    SDread_variable(sd_id,"z_spacing",&tmp,1,0.); para->z_spacing =  (int) tmp;

    BUGREPORT;
    SDread_variable(sd_id,"mue_n",&para->mue_n,1,0.);
    SDread_variable(sd_id,"mue_w",&para->mue_w,1,0.);
    SDread_variable(sd_id,"mue_t",&para->mue_t,1,0.);
    BUGREPORT;  
    SDread_variable(sd_id,"alpha",&para->alpha,1,0.);
    SDread_variable(sd_id,"beta",&para->beta,1,0.);
    SDread_variable(sd_id,"delta",&para->delta,1,0.);
    SDread_variable(sd_id,"gamma",&para->gamma,1,0.);
    SDread_variable(sd_id,"sigma",&para->sigma,1,0.);
    SDread_variable(sd_id,"nu",&para->nu,1,0.);
    BUGREPORT;


    SDread_variable(sd_id,"muehat",&para->muehat,1,0.);
    SDread_variable(sd_id,"betahat",&para->betahat,1,0.);
    
    SDread_variable(sd_id,"hm_nl",&para->hm_nl,1,0.);
    SDread_variable(sd_id,"exb_ll",&para->exb_ll,1,0.);
  
    SDread_variable(sd_id,"source",&para->source,1,0.);
    SDread_variable(sd_id,"limiter",&para->limiter,1,0.);
    SDread_variable(sd_id,"shat0",&para->shat0,1,0.);
    SDread_variable(sd_id,"q0",&para->q0,1,0.);
    SDread_variable(sd_id,"dt_pol",&para->dt_pol,1,0.);
    SDread_variable(sd_id,"r0",&para->r0,1,0.);
    SDread_variable(sd_id,"k0",&para->k0,1,0.);
    BUGREPORT;

    
    SDread_variable(sd_id,"kappan",&para->kappan,1,0.);   
    SDread_variable(sd_id,"kappat",&para->kappat,1,0.); 
    SDread_variable(sd_id,"adrhos",&para->adrhos,1,0.);
    BUGREPORT;
  
    BUGREPORT;
    SDread_variable(sd_id,"qprof",&para->qprof[0],1024,0.);    
    SDread_variable(sd_id,"shat",&para->shat[0],1024,0.); 
    
    BUGREPORT;
  
    SDread_variable(sd_id,"phiprof",&para->phiprof,1,0.);   
    SDread_variable(sd_id,"bprof",&para->bprof,1,0.);   
    SDread_variable(sd_id,"nprof",&para->nprof,1,0.);   
    SDread_variable(sd_id,"tprof",&para->tprof,1,0.); 
  
    SDread_variable(sd_id,"energy",&para->energy,1,0.);   
    SDread_variable(sd_id,"vorticity",&para->vorticity,1,0.);   
    SDread_variable(sd_id,"reynolds",&para->reynolds,1,0.);   
    SDread_variable(sd_id,"prandel",&para->prandel,1,0.); 
    SDread_variable(sd_id,"Nusselt",&para->nusselt,1,0.); 
    BUGREPORT;

  
    SDread_variable(sd_id,"xbnd",&tmp,1,0.);    para->xbnd = (long) tmp;
    SDread_variable(sd_id,"ybnd",&tmp,1,0.);    para->ybnd = (long) tmp;
    SDread_variable(sd_id,"zbnd",&tmp,1,0.);    para->zbnd = (long) tmp;

    SDread_variable(sd_id,"bndcnd",&tarray[0],3,0.);    
    for(i = 0;i<3;i++) para->boundary[i] = (long) tarray[i];

    SDread_variable(sd_id,"bdval",&para->bdval[0][0],10,0.); 
    BUGREPORT;

    SDread_variable(sd_id,"amp",&para->amp[0],3,0.);  
    SDread_variable(sd_id,"pos_x",&para->pos_x[0],3,0.);  
    SDread_variable(sd_id,"pos_y",&para->pos_y[0],3,0.);  
    SDread_variable(sd_id,"width_x",&para->width_y[0],3,0.);  
    SDread_variable(sd_id,"width_random_x",&para->width_random_x[0],3,0.);  
    SDread_variable(sd_id,"width_random_y",&para->width_random_y[0],3,0.);  
    SDread_variable(sd_id,"amp_random",&para->amp_random[0],3,0.);  


    BUGREPORT;
    SDread_variable(sd_id,"uvortex",&para->uvortex,1,0.); 
    SDread_variable(sd_id,"radius",&para->radius,1,0.);
    SDread_variable(sd_id,"epsilon",&para->epsilon,1,0.);
    SDread_variable(sd_id,"v_pol",&para->v_pol,1,0.);

    /* Primary Physics parameters */

    /* Primary units in [eV],[T],[1/m**3],[u] */
    SDread_variable(sd_id,"Ti",&para->Ti,1,0.);
    SDread_variable(sd_id,"Te",&para->Te,1,0.);
    SDread_variable(sd_id,"B0",&para->B0,1,0.);
    SDread_variable(sd_id,"n0",&para->n0,1,0.);
    SDread_variable(sd_id,"Mi",&para->Mi,1,0.);
    SDread_variable(sd_id,"Z",&para->Z,1,1.); // default 1
    
    BUGREPORT;

    SDend(sd_id);   
}


/*********************************************************************************************/
void FUtilsInt_SetupDataRange(HDF_DS *data,PARA *para)
{
 
    if(data->rank == 3)
    {
        data->range[2][0] = para->xmin;
        data->range[2][1] = para->xmax;
        data->range[1][0] = para->ymin;
        data->range[1][1] = para->ymax;
        data->range[0][0] = para->zmin;
        data->range[0][1] = para->zmax;


        data->dims[0] = data->nz;
        data->dims[1] = data->ny;      
        data->dims[2] = data->nx;         
    }
    else
    {
        data->rank = 2;
        data->range[1][0] = para->xmin;
        data->range[1][1] = para->xmax;
        data->range[0][0] = para->ymin;
        data->range[0][1] = para->ymax;

        data->dims[0] = data->ny;
        data->dims[1] = data->nx;      
        data->dims[2] = 1;
        
        
        para->zmax = 1.;
        para->zmin = 0.;

    }
  
    COMM(fprintf(stderr,"rank %d\n", (int)data->rank););
}






/*********************************************************************************************/

int FUtilsInt_ReadHDF4ByName(char *name,int32 number,const char *sdsname, HDF_DS *data,PARA *para)
{
    register int 
        i;
    int 
        noe=0;
    float 
      *flt_pointer = NULL;
    double
      *dbl_pointer = NULL;
    VOIDP 
      data_pointer = NULL;

    int32
        rank=0,
        sd_id,
        sds_id = -1,
        k,
        num_type,
        attributes;

    char	
        filename[DEFSTRLEN];
  
 
    if(data->rank == 0 )
    {
        fprintf(stderr,"Error: Rank not set in main\n");
        return(-1); 
    }

    BUGREPORT;
    if (FUtilsInt_MakeHDFFilename(filename,name,number) != 0) return -1;
    BUGREPORT;
    if( (sd_id = SDstart(filename,DFACC_RDONLY)) == -1)
    {
        fprintf(stderr,"Error on SDstart: >%s<\n",filename);
        return(-1); 
    }
    
    BUGREPORT;
    if ((k = SDnametoindex(sd_id, sdsname)) == -1)
    {
      fprintf(stderr,"Dataset >%s< not found  in %s,process %d\n",sdsname,filename,(int)data->this_process);
        SDend(sd_id);
        return(-1);
    }
     
    BUGREPORT;
    sds_id = SDselect(sd_id,k);
    SDgetinfo(sds_id,&data->name[0],&rank,data->dims,&num_type,&attributes);
    if(rank != data->rank)
      fprintf(stderr,"%s %s %d:Warning Rank %d does not equal required rank %d\n",__FILE__,__func__,__LINE__, (int)rank, (int)data->rank);   
    BUGREPORT;
    
    noe=1;
    for(i=0;i<data->rank;i++) noe*=(data->elements[i]);
    
    
    BUGREPORT;
    if(data->datarw != NULL && data->datarw_size < 2*noe  ) { free(data->datarw), data->datarw_size = 0;}     
    if( (data->datarw = (void*) calloc(2*noe,sizeof(double)))  == NULL) 
      {
	fprintf(stderr,"%s %s %d: could not allocate read field....\n",__FILE__,__func__,__LINE__);
	return(-1);
      }
    data->datarw_size = 2*noe;
  
    
    BUGREPORT;

    if( num_type == DFNT_FLOAT32)
    {
        flt_pointer = (float*)data->datarw;
	flt_pointer += noe-1; 
	data_pointer = (VOIDP)flt_pointer;
        COMM(fprintf(stderr,"noe = %d\n",noe););
    }
    else
      data_pointer = (VOIDP)data->datarw;

  
    /* Read the Stuff */
    BUGREPORT;
  
    if(SDreaddata(sds_id,data->start,NULL,data->elements,data_pointer) == -1)
    {
      fprintf(stderr,"SDreaddata nicht erfolgreich in %s,process %d\n",filename,(int)data->this_process);
        SDendaccess(sds_id);
        SDend(sd_id);
        return(-1);
    }
    SDendaccess(sds_id); 
    SDend(sd_id);

    BUGREPORT;

    dbl_pointer = (double *)data->datarw;
    /* Float to double */
    if( num_type ==  DFNT_FLOAT32)        
      for(i=0;i<noe;i++,flt_pointer++) dbl_pointer[i] = (double) (*flt_pointer);
  


    BUGREPORT;

    if (data->ReadAttributes) 
    {
        BUGREPORT;
        FUtilsInt_ReadHDF4Attributes(name,number,data,para);
        COMM(fprintf(stderr,"JOBID: >%s<\n",data->jobid););
    }
    
  
    BUGREPORT;
    return(num_type);
} 
/*******************************************************************************************************/

void FUtils_Read3DFieldbyName(double ***f_0,const char *name,HDF_DS *data,PARA *para)
{
    int i,j,k,l;
    double *dblp=NULL;

    BUGREPORT;
    FUtilsInt_ReadHDF4ByName(data->name_in,data->number,name,data,para);
    BUGREPORT;
    if(data->elements[0]*data->elements[1]*data->elements[2] > data->datarw_size)
      {
	fprintf(stderr,"%s: Sizeof datarw array not compatible with dimensions\n",__func__);
	fprintf(stderr,"%s: [%d, %d, %d] and %d\n",__func__,
		(int)data->elements[0],(int)data->elements[1],(int)data->elements[2],(int)data->datarw_size);
	exit(-1);
      }

    l=0;
    dblp = (double*) data->datarw;
    for(i=0;i<data->elements[0];i++) 
        for(j=0;j<data->elements[1];j++)
            for(k=0;k<data->elements[2];k++,l++) 
                f_0[i][j][k] = dblp[l];
}
/*******************************************************************************************************/

void FUtils_Read2DFieldbyName(double **f_0,const char *name,HDF_DS *data,PARA *para)
{
    int i,j,l;
    double *dblp=NULL;
    BUGREPORT;
    FUtilsInt_ReadHDF4ByName(data->name_in,data->number,name,data,para);
    BUGREPORT;
    l=0;
    dblp = (double*) data->datarw;
    for(i=0;i<data->elements[0];i++) 
        for(j=0;j<data->elements[1];j++,l++)
            f_0[i][j] = dblp[l];
}


/*******************************************************************************************************/

void FUtils_Read3DFieldbyNumber(double ***f_0,int number,HDF_DS *data,PARA *para)
{
    int i,j,k,l;
    double *dblp=NULL;
    BUGREPORT;
    FUtilsInt_ReadHDF4ByNumber( data->name_in,data->number,number,data,para);
    BUGREPORT;

    if(data->elements[0]*data->elements[1]*data->elements[2] > data->datarw_size)
      {
	fprintf(stderr,"%s: Sizeof datarw array not compatible with dimensions\n",__func__);
	fprintf(stderr,"%s: [%d, %d, %d] and %d\n",__func__,
		(int)data->elements[0],(int)data->elements[1],(int)data->elements[2],(int)data->datarw_size);
	exit(-1);
      }

    l=0;
    dblp = (double*)data->datarw;
    for(i=0;i<data->elements[0];i++) 
        for(j=0;j<data->elements[1];j++)
            for(k=0;k<data->elements[2];k++,l++) 
                f_0[i][j][k] = dblp[l];
}
/*******************************************************************************************************/

void FUtils_Read2DFieldbyNumber(double **f_0,int number,HDF_DS *data,PARA *para)
{
    int i,j,l;
    double *dblp = NULL;
    BUGREPORT;
    FUtilsInt_ReadHDF4ByNumber( data->name_in,data->number,number,data,para);
    BUGREPORT;
    l=0;
    dblp = (double*)data->datarw;
    for(i=0;i<data->elements[0];i++) 
        for(j=0;j<data->elements[1];j++,l++)
            f_0[i][j] = dblp[l];
}



/*******************************************************************************************************/



int FUtilsInt_ReadHDF4Coordinates(const char *name,int32 number, HDF_DS *data,PARA *para)
{
   
    int 
        i;
    
    int32
        sd_id      = -1,
        sds_id     = -1,
        rank       = 0,
        attributes = 0,
        dim_id     = -1,
        count      = 0,
        num_type   = 0;
    char	
        filename[DEFSTRLEN]; 
    

    
    if (FUtilsInt_MakeHDFFilename(filename,name,number) != 0) return(-1);
    COMM(fprintf(stderr,"JOBID: >%s<\n",data->jobid););
    
    
    if( (sd_id = SDstart(filename,DFACC_RDONLY)) == -1)
    {
        fprintf(stderr,"%s %s %d: Error on SDstart: >%s<\n",__FILE__,__func__,__LINE__,filename);
        return(-1); 
    }
    
    if((sds_id = SDselect(sd_id,0)) == -1) 
    {
        fprintf(stderr,"%s %s %d:Error on SDselect\n",__FILE__,__func__,__LINE__); 
        return(-1); 
    }
    
    SDgetinfo(sds_id,&data->name[0],&rank,data->dims,&num_type,&attributes);
    if(rank != data->rank)
      fprintf(stderr,"%s %s %d:Warning Rank %d does not equal required rank %d\n",__FILE__,__func__,__LINE__, (int)rank, (int)data->rank); 

    SDread_string(sds_id,"cordsys",data->coordsys,DEFSTRLEN-1,"cartesian"); 

    /* Read cooordinate information */
    for(i = 0;i<rank;i++) 
    {
        dim_id = SDgetdimid(sds_id,i);
        assert(SDdiminfo(dim_id,data->dim_label[i],&count,&num_type,&attributes)+1);
        COMM(fprintf(stderr,"dimlabel[%d]= >%s<\n",i,data->dim_label[i]););
        
        if(data->coordinate[i]) free(data->coordinate[i]);
        data->coordinate[i] = (double *)malloc(sizeof(double)*data->dims[i]);
        
        SDgetdimscale(dim_id,(VOIDP)&data->coordinate[i][0]);
        SDread_variable(dim_id,"range",(double*)&data->range[i],2,0.);  
        COMM(fprintf(stderr,"%s %s %d: Range %d = [%f,%f]\n",__FILE__,__func__,__LINE__,i,data->range[i][0],data->range[i][1]););
        
        
    }
    SDendaccess(sds_id);
    SDend(sd_id);
    return(0);
}


/*********************************************************************/
void FUtils_IniStructure(HDF_DS *d,PARA *p)
{
    int i,j;
    
// holds string for primary physics parameters 
    const char Str_Prim_Phys[] ={
		 "Ti:  normalising ion temperature [eV];\n"\
		 "Te:  normalising electron temperature [eV];\n"\
         "B0:  magnetic field strength [T];\n"\
         "n0:  normalising density [1e19 m**-3];\n"\
         "Mi:  ion mass in [u];\n"\
         "Z:   charge of main ion species [e];\n"\
         "p_n: neutral pressure in [pa];\n"\
         "part_source: particle source in [10 19 m-3 s-1];\n"\
         "temp_source: energy source in [ev 10 19 m-3 s-1];\n"};     
         
// secondary parameter descriptions  
    const char Str_Sec_Phys[] = {
        "rho_s:   \t 1.02 *        sqrt(Mi)*sqrt(Te) / B    \t in [cm];\n" \
        "rho_i:   \t 1.02 *        sqrt(Mi)*sqrt(Ti) / B    \t in [cm];\n" \
        "omega_ci:\t 9.58 * 10^3 * Z/Mi  B                  \t in [10^6 rad/s];\n" \
        "omega_ce:\t 1.76 * 10^7 * B                        \t in [10^6 rad/s];\n" \
        "omega_pi:\t 1.32 *        Z/sqrt(Mi)*sqrt(n0)      \t in [10^6 rad/s];\n" \
        "v_thi:   \t 9.79 * 10^3 / sqrt(Mi)*sqrt(Ti)        \t in [km/s];\n" \
        "v_the:   \t 4.19 * 10^5 * sqrt(Te)                 \t in [km/s];\n" \
        "c_s:     \t 9.79 * 10^3 * sqrt(gamma  Z Te / Mi)   \t in [km/s];\n" \
        "v_alfven:\t 2.18 * 10^12/ sqrt(Mi) /sqrt(n_i)      \t in [km/s];\n" \
        "nu_ei:   \t                                        \t in [10^6 rad/s];\n" \
        "nu_en:   \t n_n* 5e-15 *sqrt(k_b Te/ME)*1e-4       \t in [10^6 rad/s];\n" \
        "nu_in:   \t n_n* 5e-15 *sqrt(k_b Te/Mi)*1e-4       \t in [10^6 rad/s];\n" \
        "n_n:     \t                                        \t in [10^19 m^3];\n"} ;
            
// holds string with primary geometry parameter descriptions 
    const char Str_Prim_Geom[] = {
        "R0: major radius in [cm];\n"\
        "a: radius of plasma column in [cm];\n"\
        "q0: q value in middle of domain (for tokamaks);\n"\
        "shat0: magnetic shear in middle of domain (for tokamaks);\n"\
        "lpar: parallel length in [cm] (for linear devices);\n"\
        "ln: density gradient length scale  in [cm];\n"\
        "lTe: electron Temperature gradient scale length in [cm]\n;" \
        "lTi: ion Temperature gradient scale length in [cm];\n"};
            
            
    // Null structures and only deal with elements non-zero
    
    memset(d, 0, sizeof(HDF_DS)); 
    memset(p, 0, sizeof(PARA));
    
    //
    // Data structure
    //

    // at least on ghostpoint
    d->offx = 1;
    d->offy = 1;
    d->offz = 1;

    // If user does not specify any setup for domain decomposition, assume free choice in z
    d->N[0] = 0;
    d->N[1] = 1;
    d->N[2] = 1;


    snprintf(d->name         ,DEFSTRLEN,"EMPTY");
    snprintf(d->coordsys     ,DEFSTRLEN,"EMPTY");
    snprintf(d->dim_label[0] ,DEFSTRLEN,"Coord 0");
    snprintf(d->dim_label[2] ,DEFSTRLEN,"Coord 2");
    snprintf(d->dim_label[1] ,DEFSTRLEN,"Coord 1");
    
    d->read_data        = TRUE;
    d->create           = TRUE; 
    d->ReadAttributes   = TRUE;
    BUGREPORT;  

    for(i=0;i<3;i++)
        for(j=0;j<2;j++)
            d->neighbour[i][j] = -1;

    BUGREPORT;
    
    snprintf(d->revision,     DEFSTRLEN,"%s","EMPTY");
    snprintf(d->compile_date, DEFSTRLEN,"%s %s %s",__DATE__,__TIME__,__func__);
    snprintf(d->filename   ,  DEFSTRLEN,"%s", "EMPTY");
    snprintf(d->desc,         DEFSTRLEN,"%s", "EMPTY");
    snprintf(d->maschine,     DEFSTRLEN,"%s", "EMPTY");
    snprintf(d->write_date,   DEFSTRLEN,"%s", "EMPTY"); 
    snprintf(d->name_out,     DEFSTRLEN,"%s", "EMPTY");
    snprintf(d->name_in,      DEFSTRLEN,"%s", "EMPTY");
    snprintf(d->erhname,      DEFSTRLEN,"%s", "EMPTY");
    //   snprintf(d->jobid,        DEFSTRLEN,"%s", "EMPTY");

    //
    //  Parameter structure 
    //

    p->xmin     = 0.;
    p->xmax     = 1.;
    p->ymin     = 0.;
    p->ymax     = 1.;
    p->zmin     = 0.;
    p->zmax     = 1.;
    p->dt       = .01;
    p->end_time = 1.;
    p->out_time = 0.1;
    p->otmult   = 1.;
	p->Z        = 1.;

    snprintf(p->Prim_Phys,2*DEFSTRLEN,"%s", Str_Prim_Phys);     
    snprintf(p->Sec_Phys ,2*DEFSTRLEN,"%s", Str_Sec_Phys); 
    snprintf(p->Prim_Geom,2*DEFSTRLEN,"%s", Str_Prim_Geom); 
    
}


/****************************************************************************************************/

void FUtilsInt_AddToRunDatabase(HDF_DS *d, PARA *p,char **argv)
{
    char   
        buffer[DEFSTRLEN]    = "\0",
        dbasename[DEFSTRLEN] = "\0";
    int 
        i                    = 0;
    FILE 
        *output              = NULL;
   time_t  
        tp;
     

    
    BUGREPORT;
 
    snprintf(dbasename,DEFSTRLEN,"%s/%s",DATABASEPATH,p->codename);
    COMM(fprintf(stderr,"%s\n",dbasename););
    
    // Read database entries to get run no  
    if( (output = fopen(dbasename,"r")) )
    {
        BUGREPORT;   

        while(fscanf(output,"%04d",&i) == 1)
        {
            if(fgets(buffer,DEFSTRLEN-1,output) == NULL) fprintf(stderr,"error in fgets\n ");        
        }
        fclose(output);
        BUGREPORT;   
        d->run_no= i+1;
    }



    BUGREPORT;
 
    // Directory
    snprintf(d->cwd,DEFSTRLEN,">%s<",getenv("PWD"));
    
    // on which host 
    if(getenv("HOSTNAME") != NULL)    
      snprintf(d->maschine,DEFSTRLEN,"%s",getenv("HOSTNAME"));
    else 
      snprintf(d->maschine,DEFSTRLEN,"%s","HOSTNAME not defined");

    COMM(fprintf(stderr,"cwd:%s machine:>%s<\n", d->cwd,d->maschine););



    // Time of day
    if(time(&tp) != -1) strftime(d->start_date,24,"%d.%m.%Y at %H.%M.%S",localtime(&tp));
    else                snprintf(d->start_date,DEFSTRLEN,"no_date_available");
    COMM(fprintf(stderr,">%s<\n", d->start_date););



    // Insert jobid into a single string
    i=0;
    while(argv[++i] != NULL) {
        COMM(fprintf(stderr,"jobid >%s< %d: >%s< %d\n",d->jobid,i,argv[i],(int)strlen(argv[i])););
        if(i>1) strncat(d->jobid," ",1);
        strncat(d->jobid,argv[i],strlen(argv[i]));
    }	
    strncat(d->jobid,"\0",1); // End string with NULL


    // Append to database
    if ((output = fopen(dbasename,"a")))
    {
        COMM(fprintf(stderr,"%04d on >%s< >%s< in >%s<: >%s<\n",(int) d->run_no,d->maschine,d->start_date,d->cwd,d->jobid););
        fprintf(output,"%04d on >%s< >%s< in >%s<: >%s<\n",(int)d->run_no,d->maschine,d->start_date,d->cwd,d->jobid);
        fclose(output);    
    }
    else 
    {
        fprintf(stderr,"%04d on >%s< >%s< in >%s<: >%s<\n",(int) d->run_no,d->maschine,d->start_date,d->cwd,d->jobid);
        fprintf(stderr,"Could not make use of run database! Check if DATABASE is defined at compiletime!\n");
    }
 
 
    BUGREPORT;
    }
