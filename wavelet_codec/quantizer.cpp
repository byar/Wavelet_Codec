#include "quantizer.h"
#ifndef abs
#define abs(x) (x>0?x:-x)
#endif
void Quantize(double *t, uint8 *levels, uint32 tresLen, double *In, uint8 *Out, uint32 dataLen)
{
  uint32 i,step,pos;
  double curr;
  for(i=0;i<dataLen;i++)
  {
    step = tresLen>>1;
    pos = step;
    curr = In[i];
    while(step!=1)
    {
      step >>= 1;    
      pos += (2*(curr > t[pos])-1)*step;
    }
    pos -= !(curr > t[pos]);
    Out[i] = pos;
  }
  int32 mid = tresLen/2-1;
  for(i=0;i<dataLen;i++)
  {
    if(Out[i]==mid)
      Out[i] = 0;
    else if(Out[i] < mid)
    {
      Out[i] = abs((int)Out[i]-mid+1)*2+2;
    }
    else if(Out[i]>mid)
    {
      Out[i] = (Out[i]-mid+1)*2+1;
    }
  }
}