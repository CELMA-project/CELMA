#undef  DEBUG
#undef   COMM_STATUS
#define  COMM2(a)



#include "utilities.h"
#include "par_helper.h"


MPI_Datatype create_mpi_hdfds(HDF_DS *d, int root, MPI_Comm comm,int rank)
{
  int i;

  /* We have 5 markings, of the following types*/
  MPI_Datatype types[5] ={MPI_DOUBLE,MPI_LONG,MPI_INT,MPI_CHAR,MPI_CHAR};
  /* How many of the types do we have */
  int blockcounts[5] = {8,31,5,4*DEFSTRLEN,11*DEFSTRLEN};
  /* Where are the markers */
  MPI_Aint displs[5];
  //  MPI_Datatype hdfdata_type;

  /* initialize marker addresses */
  MPI_Address (&d->range[0][0],&displs[0]);
  MPI_Address (&d->offx,&displs[1]);
  MPI_Address (&d->zisperiodic,&displs[2]);
  MPI_Address (&d->jobid[0],&displs[3]);
  MPI_Address (&d->name[0],&displs[4]);



  /* Make relative addresses from absolute*/
  for(i=4;i>=0;i--) displs[i]-=displs[0];

  /* Define the structure */
  MPI_Type_struct(5,blockcounts,displs,types,&hdfdata_type);
  MPI_Type_commit( &hdfdata_type);

  return hdfdata_type;
}


/**
  ***************************************************************************/

MPI_Datatype create_mpi_para(PARA *p, int root,MPI_Comm comm, int rank )
{
  int i;

  int blockcounts[10] = {58,1024,10,24,13,3*1024,12*DEFSTRLEN,9,16,6};
  MPI_Datatype types[10] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_DOUBLE,MPI_CHAR,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
  MPI_Aint displs[10];
  // MPI_Datatype para_type;

  /* initialize types */
  MPI_Address (&p->time,&displs[0]);
  MPI_Address (&p->qprof[0],&displs[1]);
  MPI_Address (&p->bdval[0][0],&displs[2]);
  MPI_Address (&p->amp[0],&displs[3]);
  MPI_Address (&p->xbnd,&displs[4]);
  MPI_Address (&p->shat[0],&displs[5]);
  MPI_Address (&p->desc[0],&displs[6]);
  MPI_Address (&p->Ti,&displs[7]);
  MPI_Address (&p->rho_s,&displs[8]);
  MPI_Address (&p->R0,&displs[9]);


  for(i=9;i>=0;i--) displs[i]-=displs[0];

  MPI_Type_struct(10,blockcounts,displs,types,&para_type);
  MPI_Type_commit( &para_type);


  return para_type;

}





/*****************************************************************************/

void PH_3D_Write(double ***f_0,const char *fname,const char *name,int number, HDF_DS *data,PARA *p,int newfile)
{

int
    nx,ny,nz,
    i,j,k,
    dest,
    dest_coords[3],
    records,
    root=0,
    field_tag = 105;

 float
   *element=NULL;

 MPI_Status
   status;

 // Just to make sure
 FUtils_AllocateReadMemory(data);


 if(fname != NULL) snprintf(data->name,DEFSTRLEN,"%s",fname);

 nx = data->elements[2];
 ny = data->elements[1];
 nz = data->elements[0];

 COMM(fprintf(stderr,"Proc %d (%d,%d,%d)\n",data->this_process,nx,ny,nz););

 BUGREPORT;
 records = nx*ny*nz;
 if(2*records > data->datarw_size)
   {
     fprintf(stderr,"%s: Sizeof datarw array not compatible with dimensions\n",__func__);
     fprintf(stderr,"%s: [%d, %d, %d] and %d\n",__func__,
             (int)data->elements[0],(int)data->elements[1],(int)data->elements[2],(int)data->datarw_size);
     exit(-1);
   }

 element = (float*)data->datarw;

 for(i=0;i<data->elements[0];i++)
   for(j=0;j<data->elements[1];j++)
     for(k=0;k<data->elements[2];k++,element++)
       *element = (float)f_0[i][j][k];

 if(data->this_process == root)
   {
     data->create = newfile;
     BUGREPORT;
     for(dest=0; dest<data->num_procs ;dest++)
       {

         /* Determine position of process to receive from in coordinate grid */
         MPI_Cart_coords(data->cart_comm,dest,3,dest_coords);

         /* calculate start coordinates for the datafield to be received */
         for(i=0;i<3;i++) data->start[i] = data->elements[i]*dest_coords[i];

         BUGREPORT;

         COMM(fprintf(stderr,"Proc %d: --> dest[%d]: (%d,%d,%d)\n",data->this_process,dest,
                      data->start[0],data->start[1],data->start[2]););

         if(dest != root)
           MPI_Recv((void*)data->datarw,records,MPI_FLOAT,dest,field_tag,MPI_COMM_WORLD,&status);

         BUGREPORT;
         FUtilsInt_WriteHDF4(name,number,data,p,data->name);
         COMM(fprintf(stderr,"Proc %d \n",data->this_process););
         data->create = FALSE;
     }
     COMM(fprintf(stderr,"%d \n",data->this_process););
   }
 else
 {
   /* If not root, just dump your stuff to root */
   COMM(fprintf(stderr,"%d: send to root \n",data->this_process););
   MPI_Send((void*)data->datarw,records,MPI_FLOAT,root,field_tag,MPI_COMM_WORLD);
   COMM(fprintf(stderr,"%d: send to root done \n",data->this_process););
 }
}


/**************************************************************************/

int PHInt_ReadHDFbyName(const char *name, HDF_DS *data,PARA *para)
{
    int
        dest_coords[3],
        i,
        dest,root=0,
        records=1,
        success=-1,
        field_tag = 125;

    MPI_Status
        status;
    MPI_Request
        request;


    records = FUtils_AllocateReadMemory(data);


    if(records < 0)
    {
        fprintf(stderr,"failed to allocate\n");
        exit(-1);
    };


  BUGREPORT;



  if(data->this_process == root)
  {
      data->read_data = TRUE;

      for(dest=data->num_procs-1;dest>= 0 ;dest--)
      {

          /* Determine position of process to send to in coordinate grid */
          MPI_Cart_coords(data->cart_comm,dest,3,dest_coords);

          COMM(fprintf(stderr,"Proc %d:determined Proc %d position as [%d,%d,%d]\n",
                       data->this_process,dest,dest_coords[0],dest_coords[1],dest_coords[2]););

          for(i=0;i<data->rank;i++) data->start[i] = data->elements[i]*dest_coords[i];


          COMM(fprintf(stderr,"Proc %d:to Proc %d: Start from [%d,%d,%d]. %d Records\n",
                        data->this_process,dest, data->start[0],data->start[1],data->start[2],records););

          success = FUtilsInt_ReadHDF4ByName(data->name_in,data->number,name,data,para);

          COMM(fprintf(stderr,"Proc %d:to Proc %d: have read data\n",
                        data->this_process,dest););

          if(dest != root)
              {
                  COMM2(fprintf(stderr,"Proc %d: Send %d records to Proc %d, tag = %d\n",data->this_process,records,dest,field_tag););
                  MPI_Isend((void*)data->datarw,records, MPI_DOUBLE,dest,field_tag,MPI_COMM_WORLD,&request);
                  MPI_Wait (&request, &status);
              }
      }
  }
  else
  {
      COMM2(fprintf(stderr,"Proc %d: rank is %d\n", data->this_process,data->rank););
      COMM2(fprintf(stderr,"Proc %d: ready to receive %d records, tag = %d\n",data->this_process,records,field_tag););
      MPI_Irecv((void*)data->datarw,records,MPI_DOUBLE,root,field_tag,MPI_COMM_WORLD,&request);
      MPI_Wait (&request, &status);
      COMM2(fprintf(stderr,"Proc %d: received %d records, tag = %d\n",data->this_process,records,field_tag););
      success = 1;

  }


 COMM(fprintf(stderr,"Proc %d:leave routine %s with return value  %d\n",
              data->this_process,__func__,success););

  return success;
}



/**************************************************************************/

int PHInt_ReadHDFbyNumber(int number, HDF_DS *data,PARA *para)
{
    int
        dest_coords[3],
        i,
        dest,root=0,
        records=1,
        success=-1,
        field_tag = 125;

    MPI_Status
        status;


  BUGREPORT;
  FUtils_AllocateReadMemory(data);
  for(i=0;i<data->rank;i++) records *= data->elements[i];


  BUGREPORT;

  if(2*records > data->datarw_size)
  {
      fprintf(stderr,"%s: Sizeof datarw array not compatible with dimensions\n",__func__);
      fprintf(stderr,"%s: [%d, %d, %d] and %d\n",__func__,
              (int)data->elements[0],(int)data->elements[1],(int)data->elements[2],(int)data->datarw_size);
      exit(-1);
  }


  if(data->this_process == root)
  {
      data->read_data = TRUE;

      for(dest=data->num_procs-1;dest>= 0 ;dest--)
      {

          /* Determine position of process to send to in coordinate grid */
          MPI_Cart_coords(data->cart_comm,dest,3,dest_coords);

          COMM2(fprintf(stderr,"Proc %d:determined Proc %d position as [%d,%d,%d]\n",data->this_process,dest,dest_coords[0],dest_coords[1],dest_coords[2]););

          for(i=0;i<3;i++) data->start[i] = data->elements[i]*dest_coords[i];


          BUGREPORT;
          COMM2(fprintf(stderr," Proc %d: Read data for  Proc %d....\n",data->this_process,dest));

          success = FUtilsInt_ReadHDF4ByNumber(data->name_in,data->number,number,data,para);

          if(dest != root)
              {
                  COMM2(fprintf(stderr,"Proc %d: Send %d records to Proc %d, tag = %d\n",data->this_process,records,dest,field_tag););
                  MPI_Send((void*)data->datarw,records, MPI_DOUBLE,dest,field_tag,MPI_COMM_WORLD);
              }
      }
  }
  else
  {
      COMM2(fprintf(stderr,"Proc %d: rank is %d\n", data->this_process,data->rank););
      COMM2(fprintf(stderr,"Proc %d:  Before receive of %d records, tag = %d\n",data->this_process,records,field_tag););
      MPI_Recv((void*)data->datarw,records,MPI_DOUBLE,root,field_tag,MPI_COMM_WORLD,&status);
      success = 1;

  }
  return success;
}


/******************************************************************/

int  PH_Read3DFieldbyName(double ***f_0,const char *name,HDF_DS *data,PARA *para)
{
    int i,j,k,l,success;
    double *dblp = NULL;
    success = PHInt_ReadHDFbyName(name,data,para);
    l = 0;



    if(data->elements[0]*data->elements[1]*data->elements[2] > data->datarw_size)
      {
          fprintf(stderr,"%s: Sizeof datarw array not compatible with dimensions\n",__func__);
          fprintf(stderr,"%s: [%d, %d, %d] and %d\n",__func__,
                  (int)data->elements[0],(int)data->elements[1],(int)data->elements[2],(int)data->datarw_size);
          exit(-1);
      }

    dblp = (double*) data->datarw;
    for(i=0;i<data->elements[0];i++)
        for(j=0;j<data->elements[1];j++)
            for(k=0;k<data->elements[2];k++,l++)
                f_0[i][j][k] = dblp[l];
    return success;
}

/******************************************************************/

int  PH_Read3DFieldbyNumber(double ***f_0,int number,HDF_DS *data,PARA *para)
{
    int i,j,k,l,success;
    double *dblp = NULL;
    success = PHInt_ReadHDFbyNumber(number,data,para);
    l = 0;

    if(data->elements[0]*data->elements[1]*data->elements[2] > data->datarw_size)
      {
        fprintf(stderr,"%s: Sizeof datarw array not compatible with dimensions\n",__func__);
        fprintf(stderr,"%s: [%d, %d, %d] and %d\n",__func__,
                (int)data->elements[0],(int)data->elements[1],(int)data->elements[2],(int)data->datarw_size);
        exit(-1);
      }


    dblp =  (double*)data->datarw;
    for(i=0;i<data->elements[0];i++)
        for(j=0;j<data->elements[1];j++)
            for(k=0;k<data->elements[2];k++,l++)
                f_0[i][j][k] = dblp[l];
    return success;
}


/******************************************************************/

int  PH_Read2DFieldbyName(double **f_0,const char* name,HDF_DS *data,PARA *para)
{
    int i,j,l,success;
    double *dblp = NULL;

    success = PHInt_ReadHDFbyName(name,data,para);

    l = 0;
    dblp = (double*)data->datarw;
    for(i=0;i<data->elements[0];i++)
        for(j=0;j<data->elements[1];j++,l++)
                f_0[i][j] = dblp[l];
    return success;
}


/******************************************************************/

int  PH_Read2DFieldbyNumber(double **f_0,int number,HDF_DS *data,PARA *para)
{
    int i,j,l,success;
    double *dblp=NULL;

    success = PHInt_ReadHDFbyNumber(number,data,para);

    l = 0;
    dblp = (double*)data->datarw;
    for(i=0;i<data->elements[0];i++)
        for(j=0;j<data->elements[1];j++,l++)
                f_0[i][j] = dblp[l];
    return success;
}










/***********************************************************************************/

/*
  Determines Geometry of Processes in 3D,
  We split ONLY radial and parallel here ....although
  the routine supports different splitting......

*/
void  PH_3D_Geometry(HDF_DS *data)
{
int i,j;

int neighbour_coords[3];
int neighbour_rank;

// Capital letter variables refer to the meta-grid of processors

if(data->num_procs == 0)  {
    fprintf( stderr,"Incompatible Number %d of processes! Something is terribly wrong\n",(int)data->num_procs);
    MPI_Finalize();
    exit(-1);
}

BUGREPORT;

if((data->N[0]*data->N[1]*data->N[2]) == 0) // process grid not specified from outside
{
    if(data->num_procs <= data->dims[0])
    {
        /* If number of processors is same or less than number of points in Z direction
           then parallelisation is done only in Z direction                   */

        data->N[0] = 0; /*Z direction domain decomposited,
                          num of procs automatically filled in by MPI_Dims_create*/
        data->N[1] = 1; //Y direction is not
        data->N[2] = 1; //X direction is
    }
    else
    {
        /* If number of processors is greater than number of points in Z direction,
           parallellisation is done in both Z and X directions                    */
        if( data->num_procs%data->dims[0] != 0) {
            fprintf(stderr,"%s: Incompatible Number %d of processes!\n",__func__,(int)data->num_procs);
            MPI_Finalize();
            exit(-1);
        }

        if(data->dims[0] < data-> num_procs )
        {
            fprintf(stderr,"%s: Incompatible Number %d of gridpoints in [0]!\n",__func__,(int)data->dims[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        data->N[0] = data->dims[0]; /*Z direction is parallelised, 1 layer per proc*/
        data->N[1] = 1; /*Y direction is not parallelised*/
        data->N[2] = (int)(data->num_procs/data->dims[0]); /*X direction is parallelised*/
    } /*End if-else*/
}
else // Proc grid determined from outside
{

      data->N[1] = 1; //Never parallelize y direction

}




BUGREPORT;


if (data->zisperiodic == FALSE) data->periodic[0] = 0;
else data->periodic[0] = 1;  // Z direction is sometimes periodic

data->periodic[1] = 1; // Y direction is always periodic
data->periodic[2] = 0; // X direction isn't periodic

COMM(fprintf(stderr,"%d processes! (%d,%d,%d)\n",data->num_procs, (int)data->N[2],(int)data->N[1],(int)data->N[0]););


MPI_Dims_create(data->num_procs,3,data->N);
MPI_Cart_create (MPI_COMM_WORLD,3, data->N,data->periodic,1,&data->cart_comm);
MPI_Comm_rank (data->cart_comm, &data->grid_id);
MPI_Cart_coords (data->cart_comm, data->grid_id,3,data->grid_coords);

COMM(fprintf(stderr,"%s: %d processes! Process Grid: (%d,%d,%d)\n",__FILE__,data->num_procs,
             (int)data->N[2],(int)data->N[1],(int)data->N[0]););


/*Create communicator for each x-row*/
/*Ranks in this communicator are same as x-coordinate of processor in cart_comm*/
MPI_Comm_split(data->cart_comm,data->grid_coords[0],data->grid_coords[2],&data->xrow_comm);
MPI_Comm_rank (data->xrow_comm, &data->xrow_id);


/*Create communicator for each z-row*/
/*Ranks in this communicator are same as z-coordinate of processor in cart_comm*/
MPI_Comm_split(data->cart_comm,data->grid_coords[2],data->grid_coords[0],&data->zrow_comm);
MPI_Comm_rank (data->zrow_comm, &data->zrow_id);

BUGREPORT;

/*Determine neighbours for the each processors in Cart_comm communicator*/
for(j=0;j<3;j++) for(i=0;i<2;i++) data->neighbour[j][i]=MPI_PROC_NULL;
/*For Z axis*/

if(data->periodic[0] == 1)
{
    for(i=0;i<2;i++) {
        neighbour_coords[0] = data->grid_coords[0]-1+2*i;
        neighbour_coords[1] = data->grid_coords[1];
        neighbour_coords[2] = data->grid_coords[2];
        MPI_Cart_rank (data->cart_comm,neighbour_coords, &neighbour_rank);
        data->neighbour[0][i]=neighbour_rank;

        COMM(fprintf(stderr,"Z Axis: Mycoords (%d,%d,%d), neighbour (%d,%d,%d), id = %d\n",
                data->grid_coords[2],data->grid_coords[1], data->grid_coords[0],
                     neighbour_coords[2],neighbour_coords[1],neighbour_coords[0],neighbour_rank););


    } /* End for loop*/
}
else // data->periodic[0]
{
    if(data->N[0]>1)
    {
        if( data->grid_coords[0]== 0 ) {
            i = 0;
            /*Set rank of non-existing neighbour as MPI_NULL_PROC. One can send to and receive from here, but nothing happens then*/
            data->neighbour[0][i]=MPI_PROC_NULL;

            i=1;
            neighbour_coords[0] = data->grid_coords[0]-1+2*i;
            neighbour_coords[1] = data->grid_coords[1];
            neighbour_coords[2] = data->grid_coords[2];
            MPI_Cart_rank (data->cart_comm,neighbour_coords, &neighbour_rank);
            data->neighbour[0][i]=neighbour_rank;
        }
        else if ( data->grid_coords[0] == (data->N[0]-1))
        {
            i=0;
            neighbour_coords[0] = data->grid_coords[0]-1+2*i;
            neighbour_coords[1] = data->grid_coords[1];
            neighbour_coords[2] = data->grid_coords[2];
            MPI_Cart_rank (data->cart_comm,neighbour_coords, &neighbour_rank);
            data->neighbour[0][i]=neighbour_rank;

            i = 1;
            /*Set rank of non-existing neighbour as MPI_NULL_PROC. One can send to and receive from here, but nothing happens then*/
            data->neighbour[0][i]=MPI_PROC_NULL;
        }
        else {
            for(i=0;i<2;i++) {
                neighbour_coords[0] = data->grid_coords[0]-1+2*i;
                neighbour_coords[1] = data->grid_coords[1];
                neighbour_coords[2] = data->grid_coords[2];
                MPI_Cart_rank (data->cart_comm,neighbour_coords, &neighbour_rank);
                data->neighbour[0][i]=neighbour_rank;
            }
        } /*End if-elsi if-else*/
    }

}

/*For Y axis: Not necessary needed yet as no decomposition in this direction is implemented*/

for(i=0;i<2;i++) {
        neighbour_coords[0] = data->grid_coords[0];
        neighbour_coords[1] = data->grid_coords[1]-1+2*i;
        neighbour_coords[2] = data->grid_coords[2];
        MPI_Cart_rank (data->cart_comm,neighbour_coords, &neighbour_rank);
        data->neighbour[1][i]=neighbour_rank;

    COMM(fprintf(stderr,"Y Axis: Mycoords (%d,%d,%d), neighbour (%d,%d,%d), id = %d\n",
            data->grid_coords[2],data->grid_coords[1], data->grid_coords[0],
                 neighbour_coords[2],neighbour_coords[1],neighbour_coords[0],neighbour_rank););




} /* End for loop*/

/*For X axis.*/

if(data->N[2]>1)
{
/*Do only if there is more than 1 processor in X-dimension*/
        if( data->grid_coords[2]== 0 ) {
        i = 0;
                /*Set rank of non-existing neighbour as MPI_NULL_PROC. One can send to and receive from here, but nothing happens then*/
                data->neighbour[2][i]=MPI_PROC_NULL;

                i=1;
                neighbour_coords[0] = data->grid_coords[0];
                neighbour_coords[1] = data->grid_coords[1];
                neighbour_coords[2] = data->grid_coords[2]-1+2*i;
                MPI_Cart_rank (data->cart_comm,neighbour_coords, &neighbour_rank);
                data->neighbour[2][i]=neighbour_rank;


        }
        else if ( data->grid_coords[2] == data->N[2]-1){
                i=0;
                neighbour_coords[0] = data->grid_coords[0];
                neighbour_coords[1] = data->grid_coords[1];
                neighbour_coords[2] = data->grid_coords[2]-1+2*i;
                MPI_Cart_rank (data->cart_comm,neighbour_coords, &neighbour_rank);
                data->neighbour[2][i]=neighbour_rank;

        i = 1;
                /*Set rank of non-existing neighbour as MPI_NULL_PROC. One can send to and receive from here, but nothing happens then*/
                data->neighbour[2][i]=MPI_PROC_NULL;
        }
        else {
                for(i=0;i<2;i++) {
                        neighbour_coords[0] = data->grid_coords[0];
                        neighbour_coords[1] = data->grid_coords[1];
                        neighbour_coords[2] = data->grid_coords[2]-1+2*i;
                        MPI_Cart_rank (data->cart_comm,neighbour_coords, &neighbour_rank);
                        data->neighbour[2][i]=neighbour_rank;
                }
        } /*End if-elsi if-else*/

    COMM(fprintf(stderr,"X Axis: Mycoords (%d,%d,%d), neighbour (%d,%d,%d), id = %d\n",
            data->grid_coords[2],data->grid_coords[1], data->grid_coords[0],
                 neighbour_coords[2],neighbour_coords[1],neighbour_coords[0],neighbour_rank););


}/*End if */

/*Testing */
COMM(fprintf(stderr,"\n I %i NBR:Z- %i, Z+ %i, Y- %i, Y+ %i, X- %i, X+ %i\n",data->grid_id, data->neighbour[0][0],data->neighbour[0][1],data->neighbour[1][0],data->neighbour[1][1],data->neighbour[2][0],data->neighbour[2][1]);)

COMM(fprintf(stderr,"%s: PROC: %d dims: %d %d %d N: %d %d %d\n\n",
             __FILE__,data->this_process,(int)data->dims[0], (int)data->dims[1],
             (int)data->dims[2], (int)data->N[0], (int)data->N[1], (int)data->N[2]);)


/* Determine how many grid points there are for each processor */
for(i=0;i<data->rank;i++) data->elements[i] = data->dims[i]/data->N[i];

COMM(fprintf(stderr,"%s: PROC: %d dims: %d %d %d ELEMENTS: %d %d %d\n\n",
             __FILE__, data->this_process,data->dims[0], data->dims[1], data->dims[2], data->elements[0], data->elements[1], data->elements[2]););



/*Check that there isn't incompatible number of elements*/
for(i=0;i<data->rank;i++)
    if(data->elements[i]*data->N[i] != data->dims[i])
    {
        fprintf(stderr,"Incompatible Number %d of processes nl[%d] = %d\n",
                (int)data->num_procs,i,(int)data->elements[i] );
        MPI_Abort(MPI_COMM_WORLD, 1);
    }



switch (data->rank)
{
    case 1:
        data->lnx = data->elements[0];
        break;
    case 2:
        data->lny = data->elements[0];
        data->lnx = data->elements[1];
        break;
    case 3:
        data->lnz = data->elements[0];
        data->lny = data->elements[1];
        data->lnx = data->elements[2];
        break;
    default:
        fprintf(stderr, "%s, %s: Rank %d not supported \n",__FILE__,__func__, (int)data->rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
}

BUGREPORT;

COMM( if(data->this_process == 0)
   {
     fprintf(stderr,"%s %d:%d processes distributed as (%d,%d,%d).\n",__FILE__,__LINE__,data->num_procs,(int)data->N[0],(int)data->N[1],(int)data->N[2]);
     fprintf(stderr,"%s %d: Local dims (%d,%d,%d).\n", __FILE__,__LINE__,(int)data->elements[0],(int)data->elements[1],(int)data->elements[2]);
   });


/* Get info on local indices */

for(i=0;i<data->rank;i++)
 {
     data->start[i] = data->elements[i]*data->grid_coords[i];
     data->end[i]   = data->start[i]+data->elements[i];
     COMM2(fprintf(stderr,"Proc %d: Coord[%d] from %d to %d\n",data->this_process,i,data->start[i],data->end[i]););
 }



// Print a summary of information in the processes grid

COMM(
    fprintf(stderr,"%s %d: Position (%d, %d, %d).\n",__FILE__,(int)data->this_process,
            (int)data->grid_coords[2],(int)data->grid_coords[1],(int)data->grid_coords[0]);

    fprintf(stderr,"%s %d: Neighbours (%d %d, %d %d, %d %d) .\n",__FILE__,(int)data->this_process,
            (int)data->neighbour[2][0],(int)data->neighbour[2][1],
            (int)data->neighbour[1][0],(int)data->neighbour[1][1],
            (int)data->neighbour[0][0],(int)data->neighbour[0][1] );

    fprintf(stderr,"%s %d: local dims (%d, %d, %d).\n",__FILE__,data->this_process,
            data->elements[2],data->elements[1],data->elements[0]);
);


}

/*********************************************************/

int PH_CheckCFL3D(double ***f_0,double ***vr,double ***vp,HDF_DS *data,PARA *p,
                          double **hval,double **vval,double *cflr,double *cflp)
{
  long
      i      = 0,
      result = 0,
      val    = 0;

  double
      lcflr  = 0.,
      lcflp  = 0.;

  for(i=0;i<data->lnz;i++)
    {
      result = MIN(Util_2DCheckCFL(f_0[i],vr[i],vp[i],data,p,hval[i],vval[i],cflr,cflp),result);
      lcflr  = MAX(*cflr,lcflr);
      lcflp  = MAX(*cflp,lcflp);
    }

  MPI_Allreduce(&result,&val,1,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&lcflp,cflp,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&lcflr,cflr,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);

  return (int)val;

}

/*********************************************************/

int PH_CheckCFL3D_Metric(double ***f_0,double ***vr,double ***vp,HDF_DS *data,PARA *p,
                          double **hval,double **vval,double *kai,double **grr,double **gff, double **grf,
                          double **norm_gr,double **norm_gp,  double *cflr,double *cflp)
{
  long
      i      = 0,
      result = 0,
      val    = 0;

  double
      lcflr = 0.,
      lcflp = 0.;

  for(i=0;i<data->elements[0];i++)
    {
      result = MIN(Util_CheckCFLMetric(f_0[i],vr[i],vp[i],data,p,hval[i],vval[i],kai,grr[i],gff[i],grf[i],
                                    norm_gr[i],norm_gp[i],cflr,cflp),result);
      lcflr  = MAX((*cflr),lcflr);
      lcflp  = MAX((*cflp),lcflp);
      COMM2(fprintf(stderr," CR = %g CP = %g loc. CR = %g CP = %g\n", lcflr,lcflp, *cflr,*cflp  ););
    }

  MPI_Allreduce(&result,&val,1,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&lcflp,cflp,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&lcflr,cflr,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);

  return (int)val;
}


/**********************************************************************/
/*
     Calculates the Util_Integral of quantity f to order order

*/
double  PH_Integral(double **f, int order, double *norm,int ny,int nx)
{
  double
      result = 0.,
      val    = 0.;

  Util_Integral(f,order,norm,ny,nx,&result);
  MPI_Reduce((void *)&result,(void *)&val,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  return val;
}

/**********************************************************************/
double  PH_3DIntegral(double ***f, int order, double **norm,int nz,int ny,int nx)
{
    int i;
    double
        result = 0.0,
        val    = 0.0;

    BUGREPORT;

    for(i=0;i<nz;i++)
    {
        val = 0.;
        Util_Integral(f[i],order,norm[i],ny,nx,&val);
        result += val;
    }
    BUGREPORT;
    MPI_Allreduce((void *)&result,(void *)&val,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return val;
}

/*********************************************************************
 PH_3DIntegral without MPI_Allreduce. Reduction has to be done in the calling function!

**********************************************************************/
double  PH_3DIntegral_noreduce(double ***f, int order, double **norm,int nz,int ny,int nx)
{
    int i;
    double result=0.0,val=0.0;

    BUGREPORT;

    for(i=0;i<nz;i++)
    {
        val = 0.;
        Util_Integral(f[i],order,norm[i],ny,nx,&val);
        result += val;
    }
    BUGREPORT;

    return result;
}




/**********************************************************************/
/* Calculates flux surface averaged  components of field, without radial ghost points */
void  PH_3D_FS_Average(double ***f,double *result,HDF_DS *data)
{
    static  int nx,ny,nz,nz_global;
    static double fac,*help;
    static int FIRST = TRUE;

    register int i,j,k;


      BUGREPORT;

      if(FIRST)
      {

      nx = data->elements[2];
      ny = data->elements[1];
      nz = data->elements[0];

      nz_global= data->dims[0];
      fac =1./(double)(ny*nz_global);

      help =  Util_DVector(nx,data->offx);
      FIRST = FALSE;
      }

      BUGREPORT;
      for(k=0;k<nx;k++) help[k] =0.;
      BUGREPORT;
      for(i=0;i<nz;i++)  for(j=0;j<ny;j++)  for(k=0;k<nx;k++) help[k]+=f[i][j][k]*fac;

      BUGREPORT;
      MPI_Allreduce((void *)help,(void *)result,nx, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}
/**********************************************************************/
/*!
  Fluxsurface averaged  n_0 m_0 (yz-average) components of field, with radial ghost points
*/
void     PH_3D_n0m0_Component(double ***f,double *result_n0,HDF_DS *data)
{
    static  int nx,ny,nz,nz_global;
    static double fac_n0,*help_n0;
    static int FIRST = TRUE;

    register int i,j,k;


      BUGREPORT;

      if(FIRST)
      {

      nx = data->lnx;
      ny = data->lny;
      nz = data->lnz;

      nz_global= data->dims[0];
      fac_n0 = 1./(double)(nz_global*ny);
      help_n0 =  Util_DVector(nx,data->offx);
      FIRST = FALSE;
      }

      BUGREPORT;
      // Zero results vector
      for(k=-data->offx;k<nx+data->offx;k++) help_n0[k]=0.;

      BUGREPORT;
      for(i=0;i<nz;i++) for(j=0;j<ny;j++) for(k=-data->offx;k<nx+data->offx;k++)
          help_n0[k]+=f[i][j][k];

      for(k=-data->offx;k<nx+data->offx;k++) help_n0[k]*=fac_n0;

      BUGREPORT;

      MPI_Allreduce((void *)&help_n0[-data->offx],(void *)&result_n0[-data->offx],nx+data->offx+data->offx, MPI_DOUBLE, MPI_SUM,data->zrow_comm);



}


/**********************************************************************/
/* Calculates flux n_0 (z-average) components of field */
void  PH_3D_n0_Component(double ***f,double **result_n0,HDF_DS *data)
{
    static  int nx,ny,nz,nz_global;
    static double fac_n0,**help_n0;
    static int initialise= TRUE;

    register int i,j,k;


      BUGREPORT;

      if(initialise)
      {

      nx = data->lnx;
      ny = data->lny;
      nz = data->lnz;

      nz_global= data->dims[0];
      fac_n0 = 1./(double)nz_global;
      help_n0 =  Util_DMatrix(ny,data->offy,nx,data->offx);
      initialise = FALSE;
      }

      /* Zero Helpfield */
      for(j=0;j<ny;j++)  for(k=0;k<nx;k++) help_n0[j][k]=0.;
      BUGREPORT;
      for(i=0;i<nz;i++)  for(j=0;j<ny;j++)  for(k=0;k<nx;k++) help_n0[j][k]+=f[i][j][k]*fac_n0;

      BUGREPORT;

      /* Put together averages from z-slides, no averaging over r */
      MPI_Allreduce((void *)&help_n0[0][0],(void *)&result_n0[0][0],(nx+2*data->offx)*ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      BUGREPORT;
}

/******************************************************************/
/* Non-blocking version using IRecv and Isend                     */
/* Only exchange in z-direction                                   */
/* Exchanges data->offz layers                                    */


void PH_UpdateZ_FourierField(double ***f,HDF_DS *data,MPI_Comm comm)
{
    static int
        nlayer = 1,
        slayer = 1,
        tag    = 0,
        layer  = 0,
        j      = 0,
        ex_num = 0;

    int
        nz     = data->lnz,
        ny     = data->lny,
        nx     = data->lnx,
        hn     = 0;


  static MPI_Request
      handles[12];

  static MPI_Status
      statuses[12];

  /* There are 2 neighbours (u,o) ordered as (Previous, Next) */
  /* This is made to propagate a ffted version of the field */

  if(nz < data->offz)
  {
      nlayer = data->offz;
      slayer = 1;
  }
  else
  {
      nlayer = 1;
      slayer = data->offz;
  }

  /* Number of data to transfer to other processes */
  ex_num =(nx+2)*(ny+2)*slayer;
  BUGREPORT;

  /*
     exchange next neighbours only and progressivly  to make sure the stuff works
     also in the (stupid) case  when we need to exchange offz layers
     but have only a single layer per process
  */

  for(layer=0;layer<nlayer;layer++)
  {
      for(j = 0;j<2;j++) /* lower and upper boundary */
      {
          if(data->neighbour[0][j] != MPI_PROC_NULL) /* If we have a neighbour */
          {
              /* Receive and Send Calls */
              tag = j;
              /* Receive (ny+2*offset)*(nx+2*offset)*slayer =ex_num entries to locations
                 iz = -slayer:   j=0, layer=0
                 iz = nz:        j=1, layer =0
                 iz = -2:        j=0, layer=1
                 iz = nz+1:      j=1, layer =1
              */

            MPI_Irecv(&f[j*(nz+slayer+2*layer)-(layer+slayer)][0][0],
                      ex_num,MPI_DOUBLE,data->neighbour[0][j],tag,MPI_COMM_WORLD,&handles[hn]);hn++;

              /*Revert tag*/
              tag = 1-j;
              /* Send ex_num entries from location
                 iz = 0:      j=0, layer=0
                 iz = nz-slayer:   j=1, layer =0
                 iz = 1:      j=0, layer=1
                 iz = nz-2:   j=1, layer =1
              */

              MPI_Isend(&f[j*(nz-slayer-2*layer)+layer][0][0],
                        ex_num,MPI_DOUBLE,data->neighbour[0][j],tag,MPI_COMM_WORLD,&handles[hn]);hn++;
              COMM2( fprintf(stderr,"Proc %d: comm with  %d on axis %d, %d elements, layer %d\n",
                             data->this_process,data->neighbour[0][j],0,ex_num,layer););
          }
      }

      BUGREPORT;
      /* Now wait until everything is fine */
      if(hn != 0)
      {
          COMM2(fprintf(stderr,"\n Proc %d waiting, %d messages in transition\n",data->this_process,hn); );
          MPI_Waitall(hn,handles,statuses);
          COMM2(fprintf(stderr,"\n: Proc %d waiting over\n",data->this_process););
      }
  }
}

/******************************************************************/
/* Non-blocking version using IRecv and Isend                     */
/* Only exchange in z-direction                                   */
/* Does not exchange data->offz layers, but only one                                    */
#define HPC

#ifdef HPC
void PH_UpdateZ(double ***f,HDF_DS *data,MPI_Comm comm)
{
  int
      nx,
      ny,
      nz,
      tag,
      ex_num,
      hn = 0,
      i,
      offx,
      offy;

  MPI_Request
      handles[12];
  MPI_Status
      statuses[12];


  BUGREPORT;

  nz = data->lnz;
  ny = data->lny;
  nx = data->lnx;
  offx = data->offx;
  offy = data->offy;
  BUGREPORT;

  ex_num=(nx+2*offx)*(ny+2*offy);


  BUGREPORT;

  for(i=0;i<2;i++)
    {
      tag=i;
      //Send first to up and then down
      MPI_Isend(&f[nz-1-(nz-1)*i][-offy][-offx],ex_num,MPI_DOUBLE,data->neighbour[0][1-i],tag,data->cart_comm,&handles[hn]); hn++;
      //Recv first from down and then from up
      MPI_Irecv(&f[-1+(nz+1)*i][-offy][-offx], ex_num, MPI_DOUBLE, data->neighbour[0][i],tag, data->cart_comm,&handles[hn]);hn++;
    }
  BUGREPORT;
  MPI_Waitall(hn,handles,statuses);
  BUGREPORT;
}

#else

void PH_UpdateZ(double ***f,HDF_DS *data,MPI_Comm comm)
{
  static int
      nlayer=1,slayer = 1,
      tag,
      layer,
      nz,ny,nx,
      j,ex_num,
      hn=0;
  static MPI_Request
      handles[12];
  static MPI_Status
      statuses[12];


  /* There are 2 neighbours (u,o) ordered as (Previous, Next) */

  BUGREPORT;

  nz = data->lnz;
  ny = data->lny;
  nx = data->lnx;
  hn = 0;

        /* Careful: This is made to propagate a ffted version of the field */
  ex_num =(nx+2*data->offx)*(ny+2*data->offy);

  if(nz < data->offz)
  {
      nlayer = data->offz;
      slayer = 1;
  }
  else
  {
      nlayer = 1;
      slayer = data->offz;
  }

  /* Number of data to transfer to other processes */
  ex_num =(nx+2*data->offx)*(ny+2*data->offy)*slayer;
  BUGREPORT;

  /*
     exchange next neighbours only and progressivly  to make sure the stuff works
     also in the (stupid) case  when we need to exchange offz layers
     but have only a single layer per process
  */

  for(layer=0;layer<nlayer;layer++)
  {
      for(j = 0;j<2;j++) /* lower and upper boundary */
      {
          if(data->neighbour[0][j] !=  MPI_PROC_NULL) /* If we have a neighbour */
          {
              /* Receive and Send Calls */
              tag = j;
              /* Receive (ny+2*offset)*(nx+2*offset)*slayer =ex_num entries to locations
                 iz = -slayer:   j=0, layer=0
                 iz = nz:        j=1, layer =0
                 iz = -2:        j=0, layer=1
                 iz = nz+1:      j=1, layer =1
              */
            MPI_Irecv(&f[j*(nz+slayer+2*layer)-(layer+slayer)][-data->offy][-data->offx],
                      ex_num,MPI_DOUBLE,data->neighbour[0][j],tag,MPI_COMM_WORLD,&handles[hn]);hn++;

              /*Revert tag*/
              tag = 1-j;
              /* Send ex_num entries from location
                 iz = 0:      j=0, layer=0
                 iz = nz-slayer:   j=1, layer =0
                 iz = 1:      j=0, layer=1
                 iz = nz-2:   j=1, layer =1
              */
              MPI_Isend(&f[j*(nz-slayer-2*layer)+layer][-data->offy][-data->offx],
                        ex_num,MPI_DOUBLE,data->neighbour[0][j],tag,MPI_COMM_WORLD,&handles[hn]);hn++;
              COMM2( fprintf(stderr,"Proc %d: comm with  %d on axis %d, %d elements, layer %d\n",
                             data->this_process,data->neighbour[0][j],0,ex_num,layer););
          }
      }

      BUGREPORT;
      /* Now wait until everything is fine */

      if(hn != 0)
      {
          COMM2(fprintf(stderr,"\nProc %d waiting, %d messages in transition\n",data->this_process,hn); );
          MPI_Waitall(hn,handles,statuses);
          COMM2(fprintf(stderr,"\nProc %d waiting over\n",data->this_process););
      }
  }
}
#endif


/*
  New version of UpdateZ. While the communication takes place, the FFTs of
  the inner points (i.e. all points that are not in the ghost layer) are done.
*/

void PH_UpdateZ_fft(double ***f,HDF_DS *data,MPI_Comm comm)
{

  int
      nx,
      ny,
      nz,
      tag,
      ex_num,
      hn = 0,
      offx,
      offy,
      i;


  MPI_Request
      handles[12];

  MPI_Status
      statuses[12];

  nz = data->lnz;
  ny = data->lny;
  nx = data->lnx;
  offx = data->offx;
  offy = data->offy;

  ex_num=(nx+2*offx)*((ny+2)+2*offy);

  /* Calculate fft of exchanged data. At this point work only with nlayer=1*/
  Fft_1d_2d_f(f[0],nx,ny); // User has to make sure that field f has ny+2 elements, check in calling routine (Dparallel)

  if(nz > 1) Fft_1d_2d_f(f[nz-1],nx,ny); // Careful: If we have nz == 1 we need to avoid doing  the FFT here several times
  for(i=0;i<2;i++)
    {
      tag=i;
      /*Recv first from' down' secondly from 'up'*/
      MPI_Irecv(&f[-1+(nz+1)*i][-offy][-offx], ex_num, MPI_DOUBLE, data->neighbour[0][i],tag, data->cart_comm,&handles[hn]);hn++;

     /*Send first to the 'up' and then 'down'*/
      MPI_Isend(&f[nz-1-(nz-1)*i][-offy][-offx],ex_num,MPI_DOUBLE,data->neighbour[0][1-i],tag,data->cart_comm,&handles[hn]); hn++;

    }
  /* Calculate fft of data kept local (only if nz>2), while communication is done*/
  if (nz > 2) for(i=1;i<nz-1;i++)  Fft_1d_2d_f(f[i],nx,ny); /* Careful: If we have nz == 1 we need to avoid doing  the FFT here several times .....*/

  MPI_Waitall(hn,handles,statuses);
}


/***********************************************************************************************/

/*
  This procedure is not used in the code with the current 2D decomposition.
 */
void PH_UpdateY(double ***f,HDF_DS *data)
{
  static int tag,ex_num,nx,nz,ny,offz,offx,offy;
  static int i;
  static int hn;
  static int FIRST = TRUE;

  MPI_Request handles[12];
  MPI_Status statuses[12];

  if (FIRST){
    nz = data->lnz;
    nx = data->lnx;
    ny = data->lny;

    offx = data->offx;
    offy = data->offy;
    offz = data->offz;
    ex_num=(nx+2*offx)*(ny+2*offy);

    MPI_Type_vector(1,(nx+2*offx),ex_num,MPI_DOUBLE, &ysend);
    MPI_Type_commit(&ysend);

    FIRST = FALSE;
  }

  hn=0;

  for(i=0;i<2;i++)
    {
      tag=i;

      /*Send first to 'up' and then 'down'*/
      MPI_Isend(&f[-offz][ny-1-(ny-i)*i][-offx], (nz+2*offz), ysend, data->neighbour[1][1-i], tag, data->cart_comm, &handles[hn]); hn++;

      /*Recv first from 'down' and then from 'up'*/
      MPI_Irecv(&f[-offz][-1+(ny+1)*i][-offx], (nz+2*offz), ysend, data->neighbour[1][i],tag, data->cart_comm,&handles[hn]);hn++;
    }

  MPI_Waitall(hn,handles,statuses);

}

/********************************************************************************************/
void PH_UpdateX(double ***f,HDF_DS *data)
{

        static int i,tag,ex_num,nx,nz,ny,offz,offy,offx;
        static int FIRST = TRUE;
        static int hn;

        MPI_Request handles[12];
        MPI_Status statuses[12];

        if (FIRST){
                nz = data->lnz;
                nx = data->lnx;
                ny = data->lny;

                offz = data->offz;
                offy = data->offy;
                offx = data->offx;
                ex_num=(ny+2*offy)*(nz+2*offz);

                MPI_Type_vector(ex_num,1,nx+2*offx,MPI_DOUBLE, &xsend);
                MPI_Type_commit(&xsend);

                FIRST = FALSE;
        }

        hn=0;
    COMM2(fprintf(stderr,"Proc %d: in radial communication\n",data->this_process));
        for(i=0;i<2;i++)
        {
                tag=i;
                /*Send first to 'up' secondly to 'down'*/
                MPI_Isend(&f[-offz][-offy][nx-1-(nx-1)*i], 1, xsend, data->neighbour[2][1-i], tag, data->cart_comm, &handles[hn]); hn++;

                /*Recv first from 'down' secondly 'up' */
                MPI_Irecv(&f[-offz][-offy][-1+(nx+1)*i], 1, xsend, data->neighbour[2][i],tag, data->cart_comm,&handles[hn]);hn++;
        } /* End for loop*/
        MPI_Waitall(hn,handles,statuses);
}


/***********************************************************************************************/
void PH_Update2dBoundaries(double ***f,int mbdcnd, double **a, double **b,
                           double **hval,HDF_DS *data)
{
/* loeiten: Explaination to mbdcnd is given in laplace_solver*/
  int
      i;
  static int
      nz,nx,ny;
  static int FIRST = TRUE;

  if(FIRST)
  {
    nx = data->lnx;
    ny = data->lny;
    nz = data->lnz;

    FIRST = FALSE;

  }

  BUGREPORT;
  /* Do the periodic boundaries */

  for(i=0;i<nz;i++)  Util_BdPer(f[i],NULL,NULL,NULL,nx,ny,FALSE);
  /*PH_UpdateY(f, data);*/
  BUGREPORT;


  COMM2(fprintf(stderr,"Proc %d: Before Radial boundary\n",data->this_process));
  if(data->POS[2] == 0){

    COMM2(fprintf(stderr,"Proc %d: Before Radial boundary, GRID coordinate %d\n",data->this_process,data->POS[2]));
                switch(mbdcnd)
        {
            case 1:
            case 2:
                for(i=0;i<nz;i++)
                    Util_BdDirA(f[i],a[i],hval[i],nx,ny,FALSE);
                break;
            case 3:
            case 4:
                for(i=0;i<nz;i++)
                        Util_BdNeuA(f[i],a[i],hval[i],nx,ny,FALSE);
                break;
            case 5:
            case 6:
                /* R = 0 */
                /* loeiten: Sets the inner rho BC */
                for(i=0;i<nz;i++)
                        Util_BdCylA(f[i],a[i],hval[i],nx,ny,FALSE);
                break;
            default:
                fprintf(stderr,"Error in %s %s %d, boundary r at a  %d not supported!\n",__FILE__,__func__,__LINE__,mbdcnd);
                }
        }


    BUGREPORT;
    COMM2(fprintf(stderr,"Proc %d: Before r boundary 2 \n",data->this_process));
    if(data->POS[2] == (data->N[2] -1) ){
        COMM2(fprintf(stderr,"Proc %d: Before Radial boundary, GRID coordinate %d\n",data->this_process,data->POS[2]));
        switch(mbdcnd)
        {
            case 1:
            case 4:
            case 5:
                for(i=0;i<nz;i++)
                        /* loeiten: b is the value at the boundary point*/
                        Util_BdDirB(f[i],b[i],hval[i],nx,ny,FALSE);
                break;
            case 2:
            case 3:
            case 6:
                for(i=0;i<nz;i++)
                        Util_BdNeuB(f[i],b[i],hval[i],nx,ny,FALSE);
                break;
            default:
                fprintf(stderr,"Error in %s %s %d, r-boundary %d at b  not supported!\n",__FILE__,__func__,__LINE__,mbdcnd);
                }
    }

    COMM2(fprintf(stderr,"Proc %d: Boundaries ready\n",data->this_process));

    /* Boundary in radial direction */

    if (data->N[2]>1)   PH_UpdateX(f, data);


}


/***************************************************************************/
void PH_Untwist_Local(double **res,double **f,double *shift,double dtheta,int dist,
                       HDF_DS *data,int dir)
{
/*
   OUTPUT: res, back fft transformed and shifted 2d slice, without boundaries, NOTE = ny+2, nx+2 elements!
   INPUT f, fft forward transformed slice

*/



  static int
    initialized = FALSE;

  static double
    **fp  = NULL,
    **fm  = NULL,
    **fp2 = NULL,
    **fm2 = NULL,
    dtheta_old=0.;

  double
    **trafo = NULL;

  int
    i,j,
    nx, ny;

  double
    arg = 0.,
    ky  = 0.;


/*   printf("dist = %i dir = %i dtheta = %d \n",dist,dir,dtheta); */

  nx  = data->lnx  ; ny  = data->lny;

  COMM(fprintf(stderr," nx = %d, ny = %d,  xc = %d\n", nx,ny,data->grid_coords[2]););

  if(!initialized)
    {
      fp  = Util_DMatrix(ny+2,0,nx,0);
      fm  = Util_DMatrix(ny+2,0,nx,0);
      fp2 = Util_DMatrix(ny+2,0,nx,0);
      fm2 = Util_DMatrix(ny+2,0,nx,0);
      initialized = TRUE;
    }

  BUGREPORT;


  if(dtheta_old != dtheta)
  {
      for(i=0;i<ny+2;i+=2)
      {
          ky = (double)(i/2); // Ly = 2 PI  assumed
          for(j=0;j<nx;j++)
          {
                  // note that theta_trafo also contains dtheta!!!!

                  arg = ky*(shift[j]*dtheta);   // shift contains the q profile
              fp [i][j]   = cos(    arg);
              fp [i+1][j] = sin(    arg);

              fm [i][j]   = cos(   -arg);
              fm [i+1][j] = sin(   -arg);

              fp2[i][j]   = cos( 2.*arg);
              fp2[i+1][j] = sin( 2.*arg);

              fm2[i][j]   = cos(-2.*arg);
              fm2[i+1][j] = sin(-2.*arg);

          }
      }

      /*  for(j=0;j<nx;j++)
      {
          fm [1][j] = 0.;
          fp [1][j] = 0.;
          fm2[1][j] = 0.;
          fp2[1][j] = 0.;

          fm [ny+1][j] = 0.;
          fp [ny+1][j] = 0.;
          fm2[ny+1][j] = 0.;
          fp2[ny+1][j] = 0.;
      }
      */

      dtheta_old = dtheta;
  }


  /* Select the transform */

  if(dist == 1)
    trafo = (dir != -1) ? fp : fm;
  else if (dist == 2)
    trafo = (dir != -1) ? fp2 : fm2;
  else
    exit(-1);

  BUGREPORT;
  // Perform the multiplication  fres= f * trafo

  for(i=0;i<ny+2;i+=2)
      for(j=0;j<nx;j++)
      {
          res[i][j]      = f[i][j]*trafo[i][j]     - f[i+1][j]*trafo[i+1][j];/*Real part of f_k*/
          res[i+1][j]    = f[i][j]*trafo[i+1][j]   + f[i+1][j]*trafo[i][j];  /*Imaginary part of f_k*/
      }

BUGREPORT;

// FFT backwards in place
Fft_1d_2d_b(res,nx,ny);
BUGREPORT;
}
/* ******************************************************************** */

/*! Gather local profile data and write one file */

void  PH_write_radial_profile(double *bar,char *name,HDF_DS *d)
{
    FILE *file = NULL;
    int k = 0;
    static double
        *profile_buffer = NULL;
    static int
        first = TRUE;


    COMM(fprintf(stderr,"Doing profile ....\n"););


    if(first)
    {
        profile_buffer = (double *)malloc(d->dims[2]*sizeof(double));
    }


    // profiles
    COMM(fprintf(stderr,"sending %d to join into %d\n", d->lnx,d->dims[2]););



    MPI_Allgather(&bar[0],d->lnx, MPI_DOUBLE,&profile_buffer[0],d->lnx, MPI_DOUBLE,d->xrow_comm);

    if(ISROOT)
    {
        file = fopen(name,"a");
        for(k=0;k<d->dims[2];k++)  fprintf(file,"%g\n",profile_buffer[k]);
        fprintf(file,"\n");
        fclose(file);
    }
    COMM(fprintf(stderr,"ready\n"););

    first = FALSE;


}

/*! Gather local profile data and write one file */

void  PH_write_axial_profile(double *bar,char *name,HDF_DS *d)
{
    FILE *file = NULL;
    int k = 0;
    static double
        *profile_buffer = NULL;
    static int
        first = TRUE;


    COMM(fprintf(stderr,"Doing profile ....\n"););


    if(first)
    {
        profile_buffer = (double *)malloc(d->dims[0]*sizeof(double));
    }


    // profiles
    COMM(fprintf(stderr,"sending %d to join into %d\n", d->lnz,d->dims[0]););



    MPI_Allgather(&bar[0],d->lnz, MPI_DOUBLE,&profile_buffer[0],d->lnz, MPI_DOUBLE,d->zrow_comm);

    if(ISROOT)
    {
        file = fopen(name,"a");
        for(k=0;k<d->dims[0];k++)  fprintf(file,"%f %g\n",d->coordinate[0][k],profile_buffer[k]);
        fprintf(file,"\n");
        fclose(file);
    }
    COMM(fprintf(stderr,"ready\n"););

    first = FALSE;


}
