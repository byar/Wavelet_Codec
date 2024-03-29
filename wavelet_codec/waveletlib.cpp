#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "config.h"
#include <immintrin.h>
#include "opencv2\opencv.hpp"
#include "waveletlib.h"
#include "test.h"
#include "types.h"

static const double SIMDALIGN H[] =  {  
   3.782845550699535e-02,
  -2.384946501937986e-02,
  -1.106244044184226e-01,
   3.774028556126536e-01,
   8.526986790094022e-01,
   3.774028556126537e-01,
  -1.106244044184226e-01,
  -2.384946501937986e-02,
   3.782845550699535e-02
};
static const double SIMDALIGN H1[] = {
  -6.453888262869706e-02,
  -4.068941760916406e-02,
   4.180922732216172e-01,
   7.884856164055829e-01,
   4.180922732216172e-01,
  -4.068941760916406e-02,
  -6.453888262869706e-02,
  0
};
static const double SIMDALIGN G[] = {
  -6.453888262893856e-02, //-2
   4.068941760955867e-02,
   4.180922732222124e-01,
  -7.884856164056651e-01,
   4.180922732222124e-01,
   4.068941760955867e-02,
  -6.453888262893856e-02,
   0  // 4
};
static const double SIMDALIGN G1[] = {
  -3.782845550699535e-02, //-3
  -2.384946501937986e-02,
   1.106244044184226e-01,
   3.774028556126536e-01,
  -8.526986790094022e-01,
   3.774028556126537e-01,
   1.106244044184226e-01,
  -2.384946501937986e-02,
  -3.782845550699535e-02 // 5
};

void dwt_row(double* In, double* Out, int len,int SP, int AL)
{
  register int mid = len/2,i;
  register int L = mid-2;
  register int I1;
  In = In+0;
  Out = Out+0;
  Out[0] =   In[4]*H[0]+In[3]*H[1]+In[2]*H[2]+In[1]*H[3]+In[0]*H[4]+In[1]*H[5]+In[2]*H[6]+In[3]*H[7]+In[4]*H[8];
  Out[1] =   In[2]*H[0]+In[1]*H[1]+In[0]*H[2]+In[1]*H[3]+In[2]*H[4]+In[3]*H[5]+In[4]*H[6]+In[5]*H[7]+In[6]*H[8];
  Out[mid] = In[2]*G[0]+In[1]*G[1]+In[0]*G[2]+In[1]*G[3]+In[2]*G[4]+In[3]*G[5]+In[4]*G[6];
  Out[mid+1] = In[0]*G[0]+In[1]*G[1]+In[2]*G[2]+In[3]*G[3]+In[4]*G[4]+In[5]*G[5]+In[6]*G[6];
  //#pragma omp parallel for shared(Out, In, L,mid,H,G) private(i)
  for(i=2; i<L; ++i)
  {
    I1 = 2*i;
    Out[i] = In[I1-4]*H[0]+In[I1-3]*H[1]+In[I1-2]*H[2]+In[I1-1]*H[3]+In[I1]*H[4]+In[I1+1]*H[5]+In[I1+2]*H[6]+In[I1+3]*H[7]+In[I1+4]*H[8];
    Out[mid+i] = In[I1-2]*G[0]+In[I1-1]*G[1]+In[I1]*G[2]+In[I1+1]*G[3]+In[I1+2]*G[4]+In[I1+3]*G[5]+In[I1+4]*G[6];
  }
  Out[mid-2] = In[len-8]*H[0]+In[len-7]*H[1]+In[len-6]*H[2]+In[len-5]*H[3]+In[len-4]*H[4]+In[len-3]*H[5]+In[len-2]*H[6]+In[len-1]*H[7]+In[len-2]*H[8];
  Out[mid-1] = In[len-6]*H[0]+In[len-5]*H[1]+In[len-4]*H[2]+In[len-3]*H[3]+In[len-2]*H[4]+In[len-1]*H[5]+In[len-2]*H[6]+In[len-3]*H[7]+In[len-4]*H[8];
  Out[len-2] = In[len-6]*G[0]+In[len-5]*G[1]+In[len-4]*G[2]+In[len-3]*G[3]+In[len-2]*G[4]+In[len-1]*G[5]+In[len-2]*G[6];
  Out[len-1] = In[len-4]*G[0]+In[len-3]*G[1]+In[len-2]*G[2]+In[len-1]*G[3]+In[len-2]*G[4]+In[len-3]*G[5]+In[len-4]*G[6];
  /*
  1.886053  4.346256  7.071068  9.899495 12.727922        15.556349 18.384776 21.213203 23.965974 26.959734
 -0.176777 -0.000000 -0.000000 -0.000000  0.000000        -0.000000 -0.000000 -0.000000  0.129078 -0.611709
  */
}

void dwt_row_simd(double* In, double* Out, int len,int SP, int AL)
{
  register int mid = len/2,i;
  register int L = mid-2;
  register int I1;
  In = In+0;
  Out = Out+0;
  Out[0] =   In[4]*H[0]+In[3]*H[1]+In[2]*H[2]+In[1]*H[3]+In[0]*H[4]+In[1]*H[5]+In[2]*H[6]+In[3]*H[7]+In[4]*H[8];
  Out[1] =   In[2]*H[0]+In[1]*H[1]+In[0]*H[2]+In[1]*H[3]+In[2]*H[4]+In[3]*H[5]+In[4]*H[6]+In[5]*H[7]+In[6]*H[8];
  Out[mid] = In[2]*G[0]+In[1]*G[1]+In[0]*G[2]+In[1]*G[3]+In[2]*G[4]+In[3]*G[5]+In[4]*G[6];
  Out[mid+1] = In[0]*G[0]+In[1]*G[1]+In[2]*G[2]+In[3]*G[3]+In[4]*G[4]+In[5]*G[5]+In[6]*G[6];
  //#pragma omp parallel for shared(Out, In, L,mid,H,G) private(i)
  for(i=2; i<L; ++i)
  {
    I1 = 2*i;
    double *pIN = In+I1;
    fast_mul9(pIN-4,Out+i,H);
    fast_mul8(pIN-2,Out+mid+i,H);
  }
  Out[mid-2] = In[len-8]*H[0]+In[len-7]*H[1]+In[len-6]*H[2]+In[len-5]*H[3]+In[len-4]*H[4]+In[len-3]*H[5]+In[len-2]*H[6]+In[len-1]*H[7]+In[len-2]*H[8];
  Out[mid-1] = In[len-6]*H[0]+In[len-5]*H[1]+In[len-4]*H[2]+In[len-3]*H[3]+In[len-2]*H[4]+In[len-1]*H[5]+In[len-2]*H[6]+In[len-3]*H[7]+In[len-4]*H[8];
  Out[len-2] = In[len-6]*G[0]+In[len-5]*G[1]+In[len-4]*G[2]+In[len-3]*G[3]+In[len-2]*G[4]+In[len-1]*G[5]+In[len-2]*G[6];
  Out[len-1] = In[len-4]*G[0]+In[len-3]*G[1]+In[len-2]*G[2]+In[len-1]*G[3]+In[len-2]*G[4]+In[len-3]*G[5]+In[len-4]*G[6];
  /*
  1.886053  4.346256  7.071068  9.899495 12.727922        15.556349 18.384776 21.213203 23.965974 26.959734
 -0.176777 -0.000000 -0.000000 -0.000000  0.000000        -0.000000 -0.000000 -0.000000  0.129078 -0.611709
  */
}

void fast_mul9(const double* In,double *Out, const double* IR)
{
  //static double Out[4];
#ifdef SSE2PLUS
  __m128d A = _mm_load_pd(In);
  __m128d B = _mm_load_pd(IR);
  __m128d C = _mm_mul_pd(A,B);
  __m128d D = C;
  A = _mm_load_pd(In+2);
  B = _mm_load_pd(IR+2);
  C = _mm_mul_pd(A,B);
  D = _mm_add_pd(D,C);
  A = _mm_load_pd(In+4);
  B = _mm_load_pd(IR+4);
  C = _mm_mul_pd(A,B);
  D = _mm_add_pd(D,C);
  A = _mm_load_pd(In+6);
  B = _mm_load_pd(IR+6);
  C = _mm_mul_pd(A,B);
  D = _mm_add_pd(D,C);
  _mm_store_pd(Out,D);
  *Out = Out[0]+Out[1] + In[8]*IR[8];
#elif defined(AVX)
  __m256d InData1 =  _mm256_loadu_pd(In);
  __m256d IRdata1 =  _mm256_loadu_pd(IR);
  __m256d Out1 = _mm256_mul_pd(InData1,IRdata1);
  __m256d InData2 =  _mm256_loadu_pd(In+4);
  __m256d IRdata2 =  _mm256_loadu_pd(IR+4);
  __m256d Out2 = _mm256_mul_pd(InData2,IRdata2);
  Out1 = _mm256_add_pd(Out1,Out2);
  Out1 = _mm256_hadd_pd(Out1,Out1);
  _mm256_storeu_pd(Out,_mm256_add_pd(Out1,_mm256_permute2f128_pd(Out1,Out1, 0x01)));
  *Out +=  In[8]*IR[8];
#endif
}

void fast_mul8(const double* In, double* Out, const double* IR)
{
#ifdef SSE2PLUS
#error NOT IMPLEMENTED YET
#elif defined(AVX)
__m256d InData1 =  _mm256_loadu_pd(In);
  __m256d IRdata1 =  _mm256_loadu_pd(IR);
  __m256d Out1 = _mm256_mul_pd(InData1,IRdata1);
  __m256d InData2 =  _mm256_loadu_pd(In+4);
  __m256d IRdata2 =  _mm256_loadu_pd(IR+4);
  __m256d Out2 = _mm256_mul_pd(InData2,IRdata2);
  Out1 = _mm256_add_pd(Out1,Out2);
  Out1 = _mm256_hadd_pd(Out1,Out1);
  _mm256_storeu_pd(Out,_mm256_add_pd(Out1,_mm256_permute2f128_pd(Out1,Out1, 0x01)));
#endif
}

void idwt_row(double* In, double* Out, int len,int SP, int AL)
{
  register int mid = len>>1, i, j;
  register int L=len-2;
  register int I1;
  Out[0] = In[1]*H1[1]+In[0]*H1[3]+In[1]*H1[5] + In[mid+1]*G1[7]+In[mid]*G1[5]+In[mid  ]*G1[3]+In[mid+1]*G1[1];
  Out[1] = In[1]*H1[6]+In[0]*H1[4]+In[1]*H1[2]+In[2]*H1[0] + In[mid+1]*G1[8]+In[mid]*G1[6]+In[mid  ]*G1[4]+In[mid+1]*G1[2]+In[mid+2]*G1[0];
  Out[2] = In[0]*H1[5]+In[1]*H1[3]+In[2]*H1[1] + In[mid  ]*G1[7]+In[mid]*G1[5]+In[mid+1]*G1[3]+In[mid+2]*G1[1];
  Out[3] = In[0]*H1[6]+In[1]*H1[4]+In[2]*H1[2]+In[3]*H1[0] + In[mid  ]*G1[8]+In[mid]*G1[6]+In[mid+1]*G1[4]+In[mid+2]*G1[2]+In[mid+3]*G1[0];
  //#pragma omp parallel for shared(Out, In, L,mid,H1,G1) private(i,j)
  for(i=4;i<L;i=i+2)
  {
    j = i>>1;
    I1 = mid+j;
    Out[i]   = In[j-1]*H1[5]+In[j]*H1[3]+In[j+1]*H1[1]+In[I1-2]*G1[7]+In[I1-1]*G1[5]+In[I1]*G1[3]+In[I1+1]*G1[1];
    Out[i+1] = In[j-1]*H1[6]+In[j]*H1[4]+In[j+1]*H1[2]+In[j+2]*H1[0]+In[I1-2]*G1[8]+In[I1-1]*G1[6]+In[I1]*G1[4]+In[I1+1]*G1[2]+In[I1+2]*G1[0];
  }
  j = (len-4)>>1;
  Out[len-4] = In[j-1]*H1[5]+In[j]*H1[3]+In[j+1]*H1[1] + 
      In[mid+j-2]*G1[7]+In[mid+j-1]*G1[5]+In[mid+j]*G1[3]+In[mid+j+1]*G1[1];
  Out[len-3] = In[j-1]*H1[6]+In[j]*H1[4]+In[j+1]*H1[2]+In[j+1]*H1[0] + 
      In[mid+j-2]*G1[8]+In[mid+j-1]*G1[6]+In[mid+j]*G1[4]+In[mid+j+1]*G1[2]+In[mid+j]*G1[0];
  j = (len-2)>>1;
  Out[len-2] = In[j-1]*H1[5]+In[j]*H1[3]+In[j]*H1[1] + 
      In[mid+j-2]*G1[7]+In[mid+j-1]*G1[5]+In[mid+j]*G1[3]+In[mid+j-1]*G1[1];
  Out[len-1] = In[j-1]*H1[6]+In[j]*H1[4]+In[j]*H1[2]+In[j-1]*H1[0] + 
      In[mid+j-2]*G1[8]+In[mid+j-1]*G1[6]+In[mid+j]*G1[4]+In[mid+j-1]*G1[2]+In[mid+j-2]*G1[0];
  /* 
  8.356204  -6.849453 9.522092  -6.755343 10.960155 -6.717514 12.374369 -6.717514 13.788582 -6.717514 
  15.202796 -6.717514 16.617009 -6.717514 18.031223 -5.717514 19.445436 -6.577319 20.948038 -7.011925
  */
}

void dwt_colon(double* In, double* Out, int len, int SP, int AL)
{
  register int mid = len/2,i;
  register int L = mid-2;
  register int I1;
  Out[0] =   In[4*AL]*H[0]+In[3*AL]*H[1]+In[2*AL]*H[2]+In[1*AL]*H[3]+In[0]*H[4]+In[1*AL]*H[5]+In[2*AL]*H[6]+In[3*AL]*H[7]+In[4*AL]*H[8];
  Out[1*AL] =    In[2*AL]*H[0]+In[1*AL]*H[1]+In[0]*H[2]+In[1*AL]*H[3]+In[2*AL]*H[4]+In[3*AL]*H[5]+In[4*AL]*H[6]+In[5*AL]*H[7]+In[6*AL]*H[8];
  Out[mid*AL] = In[2*AL]*G[0]+In[1*AL]*G[1]+In[0]*G[2]+In[1*AL]*G[3]+In[2*AL]*G[4]+In[3*AL]*G[5]+In[4*AL]*G[6];
  Out[(mid+1)*AL] = In[0]*G[0]+In[1*AL]*G[1]+In[2*AL]*G[2]+In[3*AL]*G[3]+In[4*AL]*G[4]+In[5*AL]*G[5]+In[6*AL]*G[6];
//  #pragma omp parallel for shared(Out, In, L,mid,H,G) private(i)
  for(i=2; i<L; ++i)
  {
    I1 = 2*i*AL;
    //Out[i*AL] = In[(2*i-4)*AL]*H[0]+In[(2*i-3)*AL]*H[1]+In[(2*i-2)*AL]*H[2]+In[(2*i-1)*AL]*H[3]+In[(2*i)*AL]*H[4]+In[(2*i+1)*AL]*H[5]+In[(2*i+2)*AL]*H[6]+In[(2*i+3)*AL]*H[7]+In[(2*i+4)*AL]*H[8];
    //Out[(mid+i)*AL] = In[(2*i-2)*AL]*G[0]+In[(2*i-1)*AL]*G[1]+In[(2*i)*AL]*G[2]+In[(2*i+1)*AL]*G[3]+In[(2*i+2)*AL]*G[4]+In[(2*i+3)*AL]*G[5]+In[(2*i+4)*AL]*G[6];
    Out[i*AL] = In[I1-4*AL]*H[0]+In[I1-3*AL]*H[1]+In[I1-2*AL]*H[2]+In[I1-AL]*H[3]+In[I1]*H[4]+In[I1+AL]*H[5]+In[I1+2*AL]*H[6]+In[I1+3*AL]*H[7]+In[I1+4*AL]*H[8];
    Out[(mid+i)*AL] = In[I1-2*AL]*G[0]+In[I1-AL]*G[1]+In[I1]*G[2]+In[I1+AL]*G[3]+In[I1+2*AL]*G[4]+In[I1+3*AL]*G[5]+In[I1+4*AL]*G[6];
  }
  Out[(mid-2)*AL] = In[(len-8)*AL]*H[0]+In[(len-7)*AL]*H[1]+In[(len-6)*AL]*H[2]+In[(len-5)*AL]*H[3]+In[(len-4)*AL]*H[4]+In[(len-3)*AL]*H[5]+In[(len-2)*AL]*H[6]+In[(len-1)*AL]*H[7]+In[(len-2)*AL]*H[8];
  Out[(mid-1)*AL] = In[(len-6)*AL]*H[0]+In[(len-5)*AL]*H[1]+In[(len-4)*AL]*H[2]+In[(len-3)*AL]*H[3]+In[(len-2)*AL]*H[4]+In[(len-1)*AL]*H[5]+In[(len-2)*AL]*H[6]+In[(len-3)*AL]*H[7]+In[(len-4)*AL]*H[8];
  Out[(len-2)*AL] = In[(len-6)*AL]*G[0]+In[(len-5)*AL]*G[1]+In[(len-4)*AL]*G[2]+In[(len-3)*AL]*G[3]+In[(len-2)*AL]*G[4]+In[(len-1)*AL]*G[5]+In[(len-2)*AL]*G[6];
  Out[(len-1)*AL] = In[(len-4)*AL]*G[0]+In[(len-3)*AL]*G[1]+In[(len-2)*AL]*G[2]+In[(len-1)*AL]*G[3]+In[(len-2)*AL]*G[4]+In[(len-3)*AL]*G[5]+In[(len-4)*AL]*G[6];
  /*
  1.886053  4.346256  7.071068  9.899495 12.727922        15.556349 18.384776 21.213203 23.965974 26.959734
 -0.176777 -0.000000 -0.000000 -0.000000  0.000000        -0.000000 -0.000000 -0.000000  0.129078 -0.611709
  */
}

void idwt_colon(double* In, double* Out, int len, int SP, int AL)
{
  register int mid = len>>1, i, j;
  register int L=len-2;
  register int I1,I2;
  Out[0]   = In[1*AL]*H1[1]+In[0]*H1[3]+In[1*AL]*H1[5] + In[(mid+1)*AL]*G1[7]+In[(mid)*AL]*G1[5]+In[(mid)*AL  ]*G1[3]+In[(mid+1)*AL]*G1[1];
  Out[1*AL] = In[1*AL]*H1[6]+In[0]*H1[4]+In[1*AL]*H1[2]+In[(2)*AL]*H1[0] + In[(mid+1)*AL]*G1[8]+In[(mid)*AL]*G1[6]+In[(mid)*AL  ]*G1[4]+In[(mid+1)*AL]*G1[2]+In[(mid+2)*AL]*G1[0];
  Out[2*AL] = In[0]*H1[5]+In[1*AL]*H1[3]+In[2*AL]*H1[1] + In[mid*AL  ]*G1[7]+In[(mid)*AL]*G1[5]+In[(mid+1)*AL]*G1[3]+In[(mid+2)*AL]*G1[1];
  Out[3*AL] = In[0]*H1[6]+In[1*AL]*H1[4]+In[2*AL]*H1[2]+In[(3)*AL]*H1[0] + In[(mid)*AL  ]*G1[8]+In[(mid)*AL]*G1[6]+In[(mid+1)*AL]*G1[4]+In[(mid+2)*AL]*G1[2]+In[(mid+3)*AL]*G1[0];
//  #pragma omp parallel for shared(Out, In, L,mid,H1,G1) private(i,j)
  for(i=4;i<L;i=i+2)
  {
    j = i>>1;
    I1 = j*AL;
    I2 = (mid+j)*AL;
    //Out[i*AL]   = In[(j-1)*AL]*H1[5]+In[j*AL]*H1[3]+In[(j+1)*AL]*H1[1] + In[(mid+j-2)*AL]*G1[7]+In[(mid+j-1)*AL]*G1[5]+In[(mid+j)*AL]*G1[3]+In[(mid+j+1)*AL]*G1[1];
    //Out[(i+1)*AL] = In[(j-1)*AL]*H1[6]+In[(j)*AL]*H1[4]+In[(j+1)*AL]*H1[2]+In[(j+2)*AL]*H1[0] + In[(mid+j-2)*AL]*G1[8]+In[(mid+j-1)*AL]*G1[6]+In[(mid+j)*AL]*G1[4]+In[(mid+j+1)*AL]*G1[2]+In[(mid+j+2)*AL]*G1[0];
    Out[i*AL]   = In[I1-AL]*H1[5]+In[I1]*H1[3]+In[I1+AL]*H1[1]+In[I2-2*AL]*G1[7]+In[I2-AL]*G1[5]+In[I2]*G1[3]+In[I2+AL]*G1[1];
    Out[(i+1)*AL] = In[I1-AL]*H1[6]+In[I1]*H1[4]+In[I1+AL]*H1[2]+In[I1+2*AL]*H1[0] + In[I2-2*AL]*G1[8]+In[I2-AL]*G1[6]+In[I2]*G1[4]+In[I2+AL]*G1[2]+In[I2+2*AL]*G1[0];

  }
  j = (len-4)>>1;
  Out[(len-4)*AL] = In[(j-1)*AL]*H1[5]+In[j*AL]*H1[3]+In[(j+1)*AL]*H1[1] + 
      In[(mid+j-2)*AL]*G1[7]+In[(mid+j-1)*AL]*G1[5]+In[(mid+j)*AL]*G1[3]+In[(mid+j+1)*AL]*G1[1];
  Out[(len-3)*AL] = In[(j-1)*AL]*H1[6]+In[(j)*AL]*H1[4]+In[(j+1)*AL]*H1[2]+In[(j+1)*AL]*H1[0] + 
      In[(mid+j-2)*AL]*G1[8]+In[(mid+j-1)*AL]*G1[6]+In[(mid+j)*AL]*G1[4]+In[(mid+j+1)*AL]*G1[2]+In[(mid+j)*AL]*G1[0];
  j = (len-2)>>1;
  Out[(len-2)*AL] = In[(j-1)*AL]*H1[5]+In[j*AL]*H1[3]+In[j*AL]*H1[1] + 
      In[(mid+j-2)*AL]*G1[7]+In[(mid+j-1)*AL]*G1[5]+In[(mid+j)*AL]*G1[3]+In[(mid+j-1)*AL]*G1[1];
  Out[(len-1)*AL] = In[(j-1)*AL]*H1[6]+In[j*AL]*H1[4]+In[j*AL]*H1[2]+In[(j-1)*AL]*H1[0] + 
      In[(mid+j-2)*AL]*G1[8]+In[(mid+j-1)*AL]*G1[6]+In[(mid+j)*AL]*G1[4]+In[(mid+j-1)*AL]*G1[2]+In[(mid+j-2)*AL]*G1[0];
  /* 
  8.356204  -6.849453 9.522092  -6.755343 10.960155 -6.717514 12.374369 -6.717514 13.788582 -6.717514 
  15.202796 -6.717514 16.617009 -6.717514 18.031223 -5.717514 19.445436 -6.577319 20.948038 -7.011925
  */
}

void wavedec(double* In, double* Out, int dim, int Level)
{
  double *temp = (double*)malloc(dim*dim*sizeof(double));
  //double *temp = new double[dim*dim];
  int i,j,k,Length = dim;
  double *tempIn = In;
  assert(temp);
  assert(In);
  assert(Out);
  for(k=Level;k>0;--k)
  {
    for(i=0;i<dim;i++)
    {
      dwt_row(tempIn+Length*i,temp+Length*i,dim,0,Length);
    }
    
    //Transpose1
    for(i=0;i<dim;++i)
    {
      for(j=0;j<dim;j++)
      {
        Out[i*Length+j] = temp[j*Length+i];
      }
    }

    for(i=0;i<dim;i++)
    {
      dwt_row(Out+Length*i,temp+Length*i,dim,0,Length);
    }

    //Transpose2
    for(i=0;i<dim;++i)
    {
      for(j=0;j<dim;j++)
      {
        Out[i*Length+j] = temp[j*Length+i];
      }
    }
    dim >>= 1;
    tempIn = Out;
  }
  free(temp);
  //delete[] temp;
}

void waverec(double* In, double* Out, int dim, int Level)
{
  double *temp = (double*)malloc(dim*dim*sizeof(double));
  //double *temp = new double[dim*dim];
  int i,j,k, Length=dim;
  double *tempIn = In;
  assert(temp);
  assert(In);
  assert(Out);
  memcpy(Out,In,Length*Length*sizeof(double));
  dim >>= Level-1;
  for(k=Level;k>0;--k)
  {
    for(i=0;i<dim;++i)
    {
      for(j=0;j<dim;j++)
      {
        temp[i*Length+j] = tempIn[j*Length+i];
      }
    }

    for(i=0;i<dim;i++)
    {
      idwt_row(temp+Length*i,Out+Length*i,dim,0,Length);
    }
  
    for(i=0;i<dim;++i)
    {
      for(j=0;j<dim;j++)
      {
        temp[i*Length+j] = Out[j*Length+i];
      }
    }

    for(i=0;i<dim;i++)
    {
      idwt_row(temp+Length*i,Out+Length*i,dim,0,Length);
    }
    dim <<= 1;
    tempIn = Out;
  }
  free(temp);
  //delete[] temp;
}

void wavedecNOMEM(double* In, double* Out, int dim, int Level)
{
  int i,j,Length = dim;
  double *temp = (double*)malloc(dim*dim*sizeof(double));
  //double *temp = new double[dim*dim];
  double *tempIn = In;
  assert(temp);
  assert(In);
  assert(Out);
  for(j=Level;j>0;--j)
  {
    for(i=0;i<dim;i++)
    {
      dwt_row(tempIn+Length*i,temp+Length*i,dim,0,Length);
    }
    for(i=0;i<dim;i++)
    {
      dwt_colon(temp+i,Out+i,dim,0,Length);
    }
    dim >>= 1;
    tempIn = Out;
  }
  free(temp);
  //delete[] temp;
}

void waverecNOMEM(double* In, double* Out, int dim, int Level)
{
  int i,j,Length = dim;
  double *temp = (double*)malloc(dim*dim*sizeof(double));
  assert(temp);
  assert(In);
  assert(Out);
  //double *temp = new double[dim*dim];
  double *tempIn = In;
  memcpy(Out,In,Length*Length*sizeof(double));
  dim >>= Level-1;
  for(j=Level;j>0;--j)
  {
    for(i=0;i<dim;i++)
    {
      idwt_colon(tempIn+i,temp+i,dim,0,Length);
    }
    for(i=0;i<dim;i++)
    {
      idwt_row(temp+Length*i,Out+Length*i,dim,0,dim);
    }
    dim <<= 1;
    tempIn = Out;
  }
  free(temp);
  //delete[] temp;
}

void wpfulldec2(double* In, double* Out, int dim, int Level)
{
  double *temp = (double*)malloc(dim*dim*sizeof(double));
  //double *temp = new double[dim*dim];
  int i,j,k,m,n,p,Length = dim;
  double *tempIn = In;
  assert(temp);
  assert(In);
  assert(Out);
  cv::Mat frame(dim,dim,CV_64F);
  uchar* PushedData = frame.data;
  p=0;
  for(k=Level;k>0;--k)
  {
    int DecTimes = 1<<(Level-k);
    for(m=0; m < DecTimes; ++m)
    for(n=0; n < DecTimes; ++n)
    {
      int Ycurr = dim*m;
      int Xcurr = dim*n;

      printf("(%d,%d) - %d\n",Xcurr,Ycurr,dim);
      for(i=0;i<dim;i++)
      {
        dwt_row(tempIn+(Ycurr+i)*Length+Xcurr,temp+(Ycurr+i)*Length+Xcurr,dim,0,Length);
      }
      //frame.data =reinterpret_cast<uchar*>(temp);
      //cv::imshow("DEBUG1", 30*frame);
      //printf("%d\n",++p);
      //cv::waitKey();

      //Transpose1
      for(i=0;i<Length;++i)
      {
        for(j=0;j<Length;j++)
        {
          Out[i*Length+j] = temp[j*Length+i];
        }
      }    
      memcpy(temp,Out,Length*Length*sizeof(double));


      for(i=0;i<dim;i++)
      {
        dwt_row(Out+(Xcurr+i)*Length+Ycurr,temp+(Xcurr+i)*Length+Ycurr,dim,0,Length);
      }
      //frame.data =reinterpret_cast<uchar*>(temp);
      //cv::imshow("DEBUG2", 30*frame);
      //printf("%d\n",++p);
      //cv::waitKey();
      //Transpose2
      for(i=0;i<Length;++i)
      {
        for(j=0;j<Length;j++)
        {
          Out[i*Length+j] = temp[j*Length+i];
        }
      }
      memcpy(temp,Out,Length*Length*sizeof(double));

      //frame.data =reinterpret_cast<uchar*>(Out);
      //cv::imshow("FINAL", 30*frame);
      //printf("%d\n",++p);
      //cv::waitKey();
      tempIn = Out;
    }
    
      dim >>= 1;
  }
  frame.data = PushedData;
  free(temp);
  //delete[] temp;
}

void wpfullrec2(double* In, double* Out, int dim, int Level)
{
  double *temp = (double*)malloc(dim*dim*sizeof(double));
  //double *temp = new double[dim*dim];
  int i,j,k,m,n,p,Length = dim;
  double *tempIn = In;
  assert(temp);
  assert(In);
  assert(Out);
  cv::Mat frame(dim,dim,CV_64F);
  uchar* PushedData = frame.data;
  p=0;
  dim >>= Level-1;
  for(k=Level;k>0;--k)
  {
    int DecTimes = 1<<(k-1);
    for(m=0; m < DecTimes; ++m)
    for(n=0; n < DecTimes; ++n)
    {
      int Ycurr = dim*m;
      int Xcurr = dim*n;
      printf("(%d,%d) - %d\n",Xcurr,Ycurr,dim);
      //Transpose1

      //frame.data =reinterpret_cast<uchar*>(tempIn);
      //cv::imshow("DEBUG1", frame);
      //printf("%d\n",++p);
      //cv::waitKey();

      for(i=0;i<Length;++i)
      {
        for(j=0;j<Length;j++)
        {
          temp[i*Length+j] = tempIn[j*Length+i];
        }
      }    
      memcpy(Out,temp,Length*Length*sizeof(double));

      for(i=0;i<dim;i++)
      {
        idwt_row(temp+(Ycurr+i)*Length+Xcurr,Out+(Ycurr+i)*Length+Xcurr,dim,0,Length);
      }


      //frame.data =reinterpret_cast<uchar*>(Out);
      //cv::imshow("DEBUG2", frame);
      //printf("%d\n",++p);
      //cv::waitKey();
      //Transpose2

      for(i=0;i<Length;++i)
      {
        for(j=0;j<Length;j++)
        {
          temp[i*Length+j] = Out[j*Length+i];
        }
      }
      memcpy(Out,temp,Length*Length*sizeof(double));


      for(i=0;i<dim;i++)
      {
        idwt_row(temp+(Xcurr+i)*Length+Ycurr,Out+(Xcurr+i)*Length+Ycurr,dim,0,Length);
      }

      //frame.data =reinterpret_cast<uchar*>(Out);
      //cv::imshow("FINAL", frame);
      //printf("%d\n",++p);
      //cv::waitKey();
      tempIn = Out;
    }
    
      dim <<= 1;
  }
  frame.data = PushedData;
  free(temp);
  //delete[] temp;
}