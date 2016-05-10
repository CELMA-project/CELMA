/* Include File for parallel Helpers */
/* Include after sles.h and da.h and utilities.h */





void par_write_vector(Vec source,char *name,int number, HDF_DS *data,PARA *p);

void par_matrix_Util_RectLaplace(Mat *A,HDF_DS *data,PARA *p,
			     int xbdcnd, double *xbdra, double *xbdrb,
			     int ybdcnd, double *ybdra, double *ybdrb,
			     double *lam);


void par_matrix_pol_laplace(Mat *A,HDF_DS *data,PARA *p,
			    int xbndcnd, double *cora,double *corb,
			    int ybndcnd, double lamd);

void par_matrix_ffd(Mat *A,HDF_DS *data,PARA *p,double *hval,double *gval,double *rcor,
			    int mbdcnd, double *cora,double *corb, double lamda,MPI_Comm comm);

void  par_prepare_boundaries(HDF_DS *data,PARA *para,Vec target,
			     double *bdra,double *bdrb,double *cora, double *corb);

void par_get_solver(SLES *sles,Mat A,MPI_Comm comm);

void par_field_to_global_vector(HDF_DS *data,Vec target, double **source,int nx,int ny);
void par_global_vector_to_field(HDF_DS *data,Vec source, double **target, int nx, int ny);

void par_field_to_da_vector(HDF_DS *data,Vec target, double **source);
void par_da_vector_to_field(HDF_DS *data,Vec target, double **source);

void par_write_da_vector(DA da,Vec lvec,Vec gvec,char *name,int number,HDF_DS *data,PARA *para,double **res);

void par_read_to_da_vector(DA da,Vec lvec,Vec gvec,HDF_DS *data,PARA *para,double **res);


int par_sh_ftfd(HDF_DS *data, PARA *p,double **f, double **w,double **help,double *bda, 
		double *bdb,double *cora, double *corb,double *rcor,
		SLES sles,Vec b,Vec approx,double lamda,double *dyya, double *dyyb);



int par_sh_ftfd_3d(HDF_DS *data,PARA *p,double ***f_0,double ***w_0,double **res, 
		   double **bdra,double **bdrb, double *cora, double *corb, double *rcor, 
		   SLES sles,Vec b, Vec approx, double lamda_mul,double **dyya, double **dyyb);












