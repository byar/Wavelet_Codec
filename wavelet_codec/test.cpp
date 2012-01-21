#include <stdio.h>
#include <math.h>
#include <Windows.h>
#include "config.h"
#include "types.h"
#include "waveletlib.h"
#pragma comment(lib, "user32.lib")

#ifdef OPENCV
#include "opencv2/opencv.hpp"
#endif


#ifdef DAVIS_CODEC
#include "WAVELET.HH"
extern FilterSet Antonini;
#endif

#define Nlen 1024
#define pi 3.141592653589793

double UpdateTime()
{
  double deltaTime;
  static int64 gTime,gLastTime;
  int64 freq;
  QueryPerformanceCounter((LARGE_INTEGER *)&gTime);  // Get current count
  QueryPerformanceFrequency((LARGE_INTEGER *)&freq); // Get processor freq
  deltaTime = (double)((double)(gTime - gLastTime)/(double)freq);
  gLastTime = gTime;
  return deltaTime;
}

void print2DArray(double *In,int dim)
{
  int i,j;
  for(i=0;i<dim;++i)
  {
    for(j=0;j<dim;j++)
    {
      printf("%5.2f ",In[i*dim+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void print1DArray(double *In,int dim)
{
  int j;
  for(j=0;j<dim;j++)
  {
    printf("%+5.2f ",In[j]);
  }
  printf("\n\n");
}

void ShowArray(double *Orig, double *Dec, double* Rec,int Dim)
{
  #ifdef OPENCV
  cv::Mat frameOrig(Dim,Dim,CV_64F,Orig);
  //cv::Mat frameOrig(Dim,Dim,CV_64F);
  //memcpy(frameOrig.data,Orig,Dim*Dim*sizeof(double));
  cv::Mat frameDec(Dim,Dim,CV_64F,Dec);
  cv::Mat frameRec(Dim,Dim,CV_64F,Rec);
  uint32 dwWidth = GetSystemMetrics(SM_CXSCREEN);
  uint32 dwHeight = GetSystemMetrics(SM_CYSCREEN);
  printf("start draw %d %d\n",dwWidth,dwHeight);
  cv::imshow("Orig", frameOrig);
  cv::imshow("Decomposition", (frameDec)*30);
  cv::imshow("Reconstruction", frameRec);
  cv::imshow("Difference", (1000000*(frameOrig-frameRec)));
  cvMoveWindow("Orig",0,0);
  cvMoveWindow("Decomposition",dwWidth/2,0);
  cvMoveWindow("Reconstruction",0,dwHeight/2);
  cvMoveWindow("Difference",dwWidth/2,dwHeight/2);
  cvWaitKey();
  cvDestroyWindow("Orig");
  cvDestroyWindow("Decomposition");
  cvDestroyWindow("Reconstruction");
  cvDestroyWindow("Difference");
  printf("end draw\n");
  #endif //OPENCV
}


double SIMDALIGN V1[Nlen] = {0.0}; 
double SIMDALIGN V2[Nlen] = {0.0};
double SIMDALIGN V3[Nlen] = {0.0};
double SIMDALIGN V4[Nlen*Nlen] = {0.0};
double SIMDALIGN V5[Nlen*Nlen] = {0.0};
double SIMDALIGN V6[Nlen*Nlen] = {0.0};

void InitTest()
{
  int i,j;
  for(i=0;i<Nlen;i++)
  {
    V1[i] = sin((double)(i)/100*pi);
  }
  for(i=0;i<Nlen;i++)
  for(j=0;j<Nlen;j++)
  {
    V4[i*Nlen+j] = sin((double)(i*13+j*17)/100*pi);
  }
#ifdef OPENCV
  cv::Mat Image = cv::imread("D:/lena512.bmp",CV_LOAD_IMAGE_GRAYSCALE);

  printf("WTF: %d %d %d\n",Image.size().width,Image.size().height,Image.step);
  //Image.convertTo(Image,CV_64F);
  for(i=0;i<Nlen;i++)
  for(j=0;j<Nlen;j++)
  {
    //V4[i*Nlen+j] = cos((double)(2*i+j)/100*pi)+sin((double)(i+2*j)/100*pi);
    V4[i*Nlen+j] = (double)Image.data[Nlen*i+j]/255;
  }
#endif
}

void WaveletTest()
{
  int i;

  //*****************************//
  //***1D DWT performance test***//
  //*****************************//
  UpdateTime();
  int Ntimes1D = 10000;
  for(i=0;i<Ntimes1D;i++)
  {
    V1[i]+=0.1;
    dwt_row(V1,V2,Nlen,0,Nlen);
    //idwt_row(V2,V3,Nlen,0,Nlen);
  }
  double tm_1d1 = UpdateTime()/Ntimes1D;
  printf("My row:         Elapsed: %f s; %+f %+f %+f\n",tm_1d1,V1[Nlen/3],V2[Nlen/3],V3[Nlen/3]);

  UpdateTime();
  for(i=0;i<Ntimes1D;i++)
  {
    dwt_colon(V4,V5,Nlen,0,Nlen);
    //idwt_colon(V6,V7,Nlen,0,Nlen);
  }
  double tm_1d2 = UpdateTime()/Ntimes1D;
  printf("My col:         Elapsed: %f s; %+f %+f %+f\n",tm_1d2,V1[Nlen/3],V2[Nlen/3],V3[Nlen/3]);
  //*****************************//

  UpdateTime();
  for(i=0;i<Ntimes1D;i++)
  {
    dwt_row_simd(V1,V2,Nlen,0,Nlen);
    //idwt_row(V2,V3,Nlen,0,Nlen);
  }
  double tm_1d3 = UpdateTime()/Ntimes1D;
  printf("My row simd:    Elapsed: %f s; %+f %+f %+f\n",tm_1d3,V1[Nlen/3],V2[Nlen/3],V3[Nlen/3]);
  
  
  //*****************************//
  //***2D DWT performance test***//
  //*****************************//
  int Ntimes2D = (int)sqrt((double)Ntimes1D)/5;
  UpdateTime();
  for(i=0;i<Ntimes2D;i++)
  {
    V4[i]+=0.1;
    wavedec(V4,V5,Nlen,4);
    waverec(V5,V6,Nlen,4); 
  }
  double tm_2d1 = UpdateTime()/Ntimes2D;
  printf("With transpose: Elapsed: %f s; %+f %+f %+f\n",tm_2d1,V4[Nlen/3],V5[Nlen/3],V6[Nlen/3]);

  UpdateTime();
  for(i=0;i<Ntimes2D;i++)
  {
    V4[i]+=0.1;
    wavedecNOMEM(V4,V5,Nlen,4);
    waverecNOMEM(V5,V6,Nlen,4);
  }
  double tm_2d2 = UpdateTime()/Ntimes2D;
  printf("No transpose:   Elapsed: %f s; %+f %+f %+f\n",tm_2d2,V4[Nlen/3],V5[Nlen/3],V6[Nlen/3]);

#ifdef DAVIS_CODEC
  Wavelet wav(&Antonini);

  UpdateTime();
  for(i=0;i<Ntimes1D;i++)
  {
    V1[i]+=0.1;
    wav.transform1d(V1,V2,Nlen,1,1);
    //wav.invert1d(V2,V3,Nlen,1,1);
  }
  double td_1d=UpdateTime()/Ntimes1D;
  printf("Davis 1d:       Elapsed: %f s; %+f %+f\n",td_1d,V1[Nlen/3],V2[Nlen/3],V3[Nlen/3]);

  UpdateTime();
  for(i=0;i<Ntimes2D;i++)
  {
    V4[i]+=0.1;
    wav.transform2d(V4,V5,Nlen,Nlen,4,1);
    wav.invert2d(V5,V6,Nlen,Nlen,4,1);
  }
  double td_2d = UpdateTime()/Ntimes2D;
  printf("Davis 2d:       Elapsed: %f s; %+f %+f %+f\n",td_2d,V4[Nlen/3],V5[Nlen/3],V6[Nlen/3]);
  
  printf("\nResults:\n");
  printf("1d_row         : my - %.2f%% %s\n",(tm_1d1>td_1d?tm_1d1/td_1d:td_1d/tm_1d1)*100,(tm_1d1>td_1d?"slower":"faster"));
  printf("1d_colon       : my - %.2f%% %s\n",(tm_1d2>td_1d?tm_1d2/td_1d:td_1d/tm_1d2)*100,(tm_1d2>td_1d?"slower":"faster"));
  printf("1d_row simd    : my - %.2f%% %s\n",(tm_1d3>td_1d?tm_1d3/td_1d:td_1d/tm_1d3)*100,(tm_1d3>td_1d?"slower":"faster"));
  printf("2d_tr          : my - %.2f%% %s\n",(tm_2d1>td_2d?tm_2d1/td_2d:td_2d/tm_2d1)*100,(tm_2d1>td_2d?"slower":"faster"));
  printf("2d_no tr       : my - %.2f%% %s\n",(tm_2d2>td_2d?tm_2d2/td_2d:td_2d/tm_2d2)*100,(tm_2d2>td_2d?"slower":"faster"));
#endif
  //*****************************//
}

void WaveletPacketTest()
{
  wpfulldec2(V4,V5,Nlen,4);
  wpfullrec2(V5,V6,Nlen,4);
  //wavedec(V4,V5,Nlen,4);
  //waverec(V5,V6,Nlen,4);
  ShowArray(V4,V5,V6,Nlen);
}

void intrinsicTest()
{
   double SIMDALIGN In[] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0};
   double SIMDALIGN IR[] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0};
   double SIMDALIGN Out[]= {0,0,0,0,0,0,0,0,0};
  fast_mul9(In,Out,IR);
  printf("fast_mul9 result = %f\n",Out[0]);
}
  