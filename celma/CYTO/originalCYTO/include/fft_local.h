#ifndef FFT_LOCAL_H_
#define FFT_LOCAL_H_


void  Fft_1d_transpose_b(double **f,double **res,int lnx,int lny);
void  Fft_1d_transpose_f(double **res,double **f,int lnx,int lny);

void  Fft_1d_2d_f(double **f,int lnx,int lny);

void  Fft_1d_2d_b(double **f,int lnx,int lny);
/* 1D backward (inverse) fft */ 
void  Fft_1d_b(double *f,int lny);
/* 1D forward fft  */ 
void  Fft_1d_f(double *f,int lny);
#endif /* !_FFT_LOCAL_H_ */
