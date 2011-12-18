#include <stdio.h>
#include <math.h>
#include <Windows.h>
//#include "WAVELET.HH"
#include "types.h"
#include "waveletlib.h"
#include "quantizer.h"
#include "stats.h"

#define Nlen 512
#define pi 3.141592653589793
//extern FilterSet Antonini;
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

double V1[Nlen] = {0.0}; 
double V2[Nlen] = {0.0};
double V3[Nlen] = {0.0};
double V4[Nlen*Nlen] = {0.0};
double V5[Nlen*Nlen] = {0.0};
double V6[Nlen*Nlen] = {0.0};
int main()
{
  int i,j;
  for(i=0;i<Nlen;i++)
  {
    V1[i] = sin((double)(i)/100*pi);
  }

  for(i=0;i<Nlen;i++)
  for(j=0;j<Nlen;j++)
  {
    V4[i*Nlen+j] = sin((double)(i+2*j)/100*pi);
  }

  //*****************************//
  //***1D DWT performance test***//
  //*****************************//
  UpdateTime();
  int Ntimes = 10;
  for(i=0;i<Ntimes;i++)
  {
    V1[i]+=0.1;
    dwt_row(V1,V2,Nlen,0,Nlen);
    idwt_row(V2,V3,Nlen,0,Nlen);
  }
  double tm_1d1 = UpdateTime()/Ntimes;
  printf("My row: Elapsed: %f %f %f\n",tm_1d1,V1[Nlen-1],V3[Nlen-1]);

  UpdateTime();
  for(i=0;i<Ntimes;i++)
  {
    dwt_colon(V4,V5,Nlen,0,Nlen);
    idwt_colon(V5,V6,Nlen,0,Nlen);
  }
  double tm_1d2 = UpdateTime()/Ntimes;
  printf("My col: Elapsed: %f %f %f\n",tm_1d2,V4[Nlen],V6[Nlen]);
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
  printf("With transpose: Elapsed: %f %f %f %f\n",tm_2d1,V4[Nlen*5+Nlen/4],V5[Nlen*5+Nlen/4],V6[Nlen*5+Nlen/4]);

  UpdateTime();
  for(i=0;i<Ntimes;i++)
  {
    V4[i]+=0.1;
    wavedecNOMEM(V4,V5,Nlen,4);
    waverecNOMEM(V5,V6,Nlen,4);
  }
  double tm_2d2 = UpdateTime()/Ntimes;
  printf("No transpose:   Elapsed: %f %f %f %f\n",tm_2d2,V4[Nlen*5+Nlen/4],V5[Nlen*5+Nlen/4],V6[Nlen*5+Nlen/4]);

  /*Wavelet wav(&Antonini);

  UpdateTime();
  for(i=0;i<Ntimes;i++)
  {
    V1[i]+=0.1;
    wav.transform1d(V1,V2,Nlen,1,1);
    wav.invert1d(V2,V3,Nlen,1,1);
  }
  double td_1d=UpdateTime()/Ntimes;
  printf("Davis 1d:       Elapsed: %f %f %f\n",td_1d,V1[Nlen-1],V3[Nlen-1]);

  UpdateTime();
  for(i=0;i<Ntimes;i++)
  {
    V4[i]+=0.1;
    wav.transform2d(V4,V5,Nlen,Nlen,4,1);
    wav.invert2d(V5,V6,Nlen,Nlen,4,1);
  }
  double td_2d = UpdateTime()/Ntimes;
  printf("Davis 2d:       Elapsed: %f %f %f %f\n",td_2d,V4[Nlen*5+Nlen/4],V5[Nlen*5+Nlen/4], V6[Nlen*5+Nlen/4]);
  
  printf("\nResults:\n");
  printf("1d_row  : my - %.2f%% %s\n",(tm_1d1>td_1d?tm_1d1/td_1d:td_1d/tm_1d1)*100,(tm_1d1>td_1d?"slower":"faster"));
  printf("1d_row  : my - %.2f%% %s\n",(tm_1d2>td_1d?tm_1d2/td_1d:td_1d/tm_1d2)*100,(tm_1d2>td_1d?"slower":"faster"));
  printf("2d_tr   : my - %.2f%% %s\n",(tm_2d1>td_2d?tm_2d1/td_2d:td_2d/tm_2d1)*100,(tm_2d1>td_2d?"slower":"faster"));
  printf("2d_no tr: my - %.2f%% %s\n",(tm_2d2>td_2d?tm_2d2/td_2d:td_2d/tm_2d2)*100,(tm_2d2>td_2d?"slower":"faster"));*/
  //*****************************//
  
  i=getchar();
  return 0;
}
