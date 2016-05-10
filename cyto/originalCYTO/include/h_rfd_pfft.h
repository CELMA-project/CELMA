


void Laplace_Solve3D(HDF_DS *data, PARA *p,double ***f, double ***w,
	       double **bdra,double **bdrb,double *hval,double *gval,double *rcor,int mbdcnd,
	       double lamda,int mormalize);

void diskplplr_3d(HDF_DS *data, PARA *p,double ***f, double ***w,
	       double **bdra,double **bdrb,double *hval,double *gval,double *rcor,int mbdcnd,
	       double lamda,int mormalize);

void vnauplplr(HDF_DS *data, PARA *p,double **f, double **w,
	       double *bda,double *bdb,double *hval,double *gval,double *rcor,int mbdcnd,
	       double lamda, int normalize);

void HHSOLV_SM_GLOB_PSI(HDF_DS *data, PARA *p,double **f, double **w,
	       double *ibdra,double *ibdrb,double *hval,double *gffm1,
               int mbdcnd, double lamda,double *n,int nfd);


void HHSOLV_SM(HDF_DS *data, PARA *p,double **f, double **w,
	       double *ibdra,double *ibdrb,double *hval,double **gffm1,
               int mbdcnd, double lamda,int scale,int nfd,int iz);

void HHSOLV_SM_noparallel(HDF_DS *data, PARA *p,double **f, double **w,
               double *ibdra,double *ibdrb,double *hval,double *gffm1,
               int mbdcnd, double lamda,int scale,int nfd);
	       
void vnauplplr_metric(HDF_DS *data, PARA *p,double **f, double **w,
	       double *bdra,double *bdrb,
	       double **grr,double **gff,
	       double **grf,double **gffm1, double *shat, int iz,
	       double *edr, int mbdcnd,double *hval, double *vval,
	       double lamda, int normalize);

int  cgtsvn(int n, double *dl, double *d, double *du, double *b);

int  complex_cgtsvn(int n, double *dl, double *d, double *du,double *b);


void HRDF_HelmholtzVarCoeff(HDF_DS *data, PARA *p,double **f, double **w,
			    double *bdra,double *bdrb,double *hval,double *gval,double *rcor,int mbdcnd,
			    double *lamda, int normalize);

void diskplplr(HDF_DS *data, PARA *p,double **f, double **w,
	       double *bdra,double *bdrb,double *hval,double *gval,double *rcor,int mbdcnd,
	       double lamda, int normalize);


