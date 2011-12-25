#include <stdio.h>
#include <math.h>
#include <Windows.h>
#include "types.h"
#include "waveletlib.h"
#pragma comment(lib, "user32.lib")
#define OPENCV
#ifdef OPENCV
#include "opencv2/opencv.hpp"
//#include "opencv/cv.h"
//#include "opencv/highgui.h"
#endif

#define DAVIS_CODEC
#ifdef DAVIS_CODEC
#include "WAVELET.HH"
extern FilterSet Antonini;
#endif

#define Nlen 512
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
}

double V1[Nlen] = {0.0}; 
double V2[Nlen] = {0.0};
double V3[Nlen] = {0.0};
double V4[Nlen*Nlen] = {0.0};
double V5[Nlen*Nlen] = {0.0};
double V6[Nlen*Nlen] = {0.0};

void InitTest()
{
  int i,j;
  for(i=0;i<Nlen;i++)
  {
    V1[i] = sin((double)(i)/100*pi);
  }
  cv::Mat Image = cv::imread("D:/lena512.bmp",CV_LOAD_IMAGE_GRAYSCALE);
  printf("WTF: %d %d %d\n",Image.size().width,Image.size().height,Image.step);
  //Image.convertTo(Image,CV_64F);
  for(i=0;i<Nlen;i++)
  for(j=0;j<Nlen;j++)
  {
    //V4[i*Nlen+j] = cos((double)(2*i+j)/100*pi)+sin((double)(i+2*j)/100*pi);
    V4[i*Nlen+j] = (double)Image.data[Nlen*i+j]/255;
  }
}

void WaveletTest()
{
  int i;

  //*****************************//
  //***1D DWT performance test***//
  //*****************************//
  UpdateTime();
  int Ntimes = 1;
  for(i=0;i<Ntimes;i++)
  {
    V1[i]+=0.1;
    dwt_row(V1,V2,Nlen,0,Nlen);
    idwt_row(V2,V3,Nlen,0,Nlen);
  }
  double tm_1d1 = UpdateTime()/Ntimes;
  printf("My row:         Elapsed: %f s; %+f %+f %+f\n",tm_1d1,V1[Nlen/3],V2[Nlen/3],V3[Nlen/3]);

  UpdateTime();
  for(i=0;i<Ntimes;i++)
  {
    dwt_colon(V4,V5,Nlen,0,Nlen);
    idwt_colon(V5,V6,Nlen,0,Nlen);
  }
  double tm_1d2 = UpdateTime()/Ntimes;
  printf("My col:         Elapsed: %f s; %+f %+f %+f\n",tm_1d2,V4[Nlen*(Nlen/4)],V5[Nlen*(Nlen/4)], V6[Nlen*(Nlen/4)]);
  //*****************************//
  

  //*****************************//
  //***2D DWT performance test***//
  //*****************************//
  UpdateTime();
  for(i=0;i<Ntimes;i++)
  {
    V4[i]+=0.1;
    wavedec(V4,V5,Nlen,4);
    waverec(V5,V6,Nlen,4); 
  }
  double tm_2d1 = UpdateTime()/Ntimes;
  printf("With transpose: Elapsed: %f s; %+f %+f %+f\n",tm_2d1,V4[Nlen*5+Nlen/4],V5[Nlen*5+Nlen/4],V6[Nlen*5+Nlen/4]);

  UpdateTime();
  for(i=0;i<Ntimes;i++)
  {
    V4[i]+=0.1;
    wavedecNOMEM(V4,V5,Nlen,4);
    waverecNOMEM(V5,V6,Nlen,4);
  }
  double tm_2d2 = UpdateTime()/Ntimes;
  printf("No transpose:   Elapsed: %f s; %+f %+f %+f\n",tm_2d2,V4[Nlen*5+Nlen/4],V5[Nlen*5+Nlen/4],V6[Nlen*5+Nlen/4]);

#ifdef DAVIS_CODEC
  Wavelet wav(&Antonini);

  UpdateTime();
  for(i=0;i<Ntimes;i++)
  {
    V1[i]+=0.1;
    wav.transform1d(V1,V2,Nlen,1,1);
    wav.invert1d(V2,V3,Nlen,1,1);
  }
  double td_1d=UpdateTime()/Ntimes;
  printf("Davis 1d:       Elapsed: %f s; %+f %+f\n",td_1d,V1[Nlen-1],V3[Nlen-1]);

  UpdateTime();
  for(i=0;i<Ntimes;i++)
  {
    V4[i]+=0.1;
    wav.transform2d(V4,V5,Nlen,Nlen,4,1);
    wav.invert2d(V5,V6,Nlen,Nlen,4,1);
  }
  double td_2d = UpdateTime()/Ntimes;
  printf("Davis 2d:       Elapsed: %f s; %+f %+f %+f\n",td_2d,V4[Nlen*5+Nlen/4],V5[Nlen*5+Nlen/4], V6[Nlen*5+Nlen/4]);
  
  printf("\nResults:\n");
  printf("1d_row  : my - %.2f%% %s\n",(tm_1d1>td_1d?tm_1d1/td_1d:td_1d/tm_1d1)*100,(tm_1d1>td_1d?"slower":"faster"));
  printf("1d_colon: my - %.2f%% %s\n",(tm_1d2>td_1d?tm_1d2/td_1d:td_1d/tm_1d2)*100,(tm_1d2>td_1d?"slower":"faster"));
  printf("2d_tr   : my - %.2f%% %s\n",(tm_2d1>td_2d?tm_2d1/td_2d:td_2d/tm_2d1)*100,(tm_2d1>td_2d?"slower":"faster"));
  printf("2d_no tr: my - %.2f%% %s\n",(tm_2d2>td_2d?tm_2d2/td_2d:td_2d/tm_2d2)*100,(tm_2d2>td_2d?"slower":"faster"));
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