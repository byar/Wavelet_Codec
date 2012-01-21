#ifndef _WAVELETLIB__H_
#define _WAVELETLIB__H_

void print1DArray(double *In,int dim);
void print2DArray(double *In,int dim);
void dwt_row(double *In, double* Out, int len, int SP, int AL);
void dwt_row_simd(double *In, double* Out, int len, int SP, int AL);
void idwt_row(double* In, double* Out, int len, int SP, int AL);
void dwt_colon(double* In, double* Out, int len, int SP, int AL);
void idwt_colon(double* In, double* Out, int len, int SP, int AL);
void wavedec(double* In, double* Out, int dim, int Level);
void waverec(double* In, double* Out, int dim, int Level);
void wavedecNOMEM(double* In, double* Out, int dim, int Level);
void waverecNOMEM(double* In, double* Out, int dim, int Level);
void wpfulldec2(double* In, double* Out, int dim, int Level);
void wpfullrec2(double* In, double* Out, int dim, int Level);

void fast_mul9(const double* In, double* Out, const double* IR);
void fast_mul8(const double* In, double* Out, const double* IR);

#endif //#ifndef __WAVELETLIB__H__