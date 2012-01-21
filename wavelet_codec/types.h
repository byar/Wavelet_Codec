#ifndef _TYPES__H_
#define _TYPES__H_

typedef unsigned char uint8;
typedef signed char int8;
typedef unsigned short int uint16;
typedef signed short int int16;
typedef unsigned int uint32;
typedef signed int int32;
typedef unsigned long long uint64;
typedef signed long long int64;

#ifdef AVX 
#define SIMDALIGN __declspec(align(32))
#elif defined SSE2PLUS
#define SIMDALIGN __declspec(align(16))
#endif

#endif //#ifndef _TEST__H_