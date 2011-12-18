#include "stats.h"
double var0(double *D, int32 dim, int32 Length)
{
  int32 i,j;
  double Std = 0;
  for(i=0;i<dim;i++)
  for(j=0;j<dim;j++)
  {
    Std += D[i*Length+j]*D[i*Length+j];
  }
  return Std;
}