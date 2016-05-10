#ifdef CPP
extern "C"
{
#endif
void laplace_2( double **f, double **g,long nx,long ny);

void AblY_2(double **f, double **g, long nx, long ny);

void AblX_2(double **f, double **g,long nx,long ny);

void laplace_4( double **f, double **g,long nx,long ny);

void AblY_4(double **f, double **g, long nx, long ny);

void AblX_4(double **f, double **g,long nx,long ny);
#ifdef CPP
}
#endif
