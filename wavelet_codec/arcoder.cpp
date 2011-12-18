// Arithmetic codec
// Almost original (miserable modifications) variant
#include <stdio.h>
#include <process.h>
#include <stdlib.h>
#include <math.h>
#include "arcoder.h"
#define SHIFT 14  // defines MAX_FREQUENCY and BITS_IN_REGISTER,
// to be empirically adjusted (<=15)
#define MAX_FREQUENCY ((unsigned)1<<SHIFT)
#define BITS_IN_REGISTER (31-SHIFT)
#define TOP_VALUE (((unsigned long) 1<<BITS_IN_REGISTER) -1)
#define FIRST_QTR ((TOP_VALUE>>2) +1)
#define HALF      (2*FIRST_QTR)
#define THIRD_QTR (3*FIRST_QTR)
int  NO_OF_CHARS;
int EOF_SYMBOL;// (NO_OF_CHARS+1)
int NO_OF_SYMBOLS;// (NO_OF_CHARS+1)
bool Initialized = false;
unsigned char         *index_to_char;//[NO_OF_SYMBOLS];
int                   *char_to_index;//[NO_OF_CHARS];
unsigned int          *cum_freq;//[NO_OF_SYMBOLS+1];
unsigned int          *freq;//[NO_OF_SYMBOLS+1];
unsigned long         low, high, value, bits_to_follow;
int                   buffer, bits_to_go, garbage_bits;
FILE *in, *out;
int Nlevels;
void update_model(int symbol);
void start_model(void)
{
	int i;
	for (i=0; i<NO_OF_CHARS; i++)
	{
		char_to_index[i]=i+1;
		index_to_char[i+1]=i;
	}
	for (i=0; i<=NO_OF_SYMBOLS; i++)
	{
		freq[i]=1;
		cum_freq[i]=NO_OF_SYMBOLS-i;
	}
	freq[0]=0;

	int Narr = MAX_FREQUENCY/4;
	if(Nlevels!=0)
	{
		for(int j=0;j<Narr;j++)
		{
			update_model(char_to_index[0]);
		}
		for(int i=1;i<Nlevels/3;i+=2)
		{
			int Ncurr = (int)(Narr*exp(-pow(60*(double)i/Nlevels,0.7)));
			for(int j=0;j<Ncurr;j++)
			{
				int num = char_to_index[i];
				update_model(num);
				num = char_to_index[i+1];
				update_model(num);
			}
		}
	}
}

void update_model(int symbol)
{
	int i, ch_i, ch_symbol, cum;
	if (cum_freq[0]==MAX_FREQUENCY)
	{
		cum=0;
		for (i=NO_OF_SYMBOLS; i>=0; i--)
		{
			freq[i]=(freq[i]+1)>>1;
			cum_freq[i]=cum;
			cum+=freq[i];
		}
	}
	for (i=symbol; freq[i]==freq[i-1]; i--);
	if (i<symbol)
	{
		ch_i                    =index_to_char[i];
		ch_symbol               =index_to_char[symbol];
		index_to_char[i]        =ch_symbol;
		index_to_char[symbol]   =ch_i;
		char_to_index[ch_i]     =symbol;
		char_to_index[ch_symbol]=i;
	}
	freq[i]++;
	while (i>0)
	{
		i--;
		cum_freq[i]++;
	}
}

void start_inputing_bits(void)
{
	bits_to_go=0;
	garbage_bits=0;
}

int input_bit(void)
{
	int t;
	if (bits_to_go==0)
	{
		buffer=getc(in);
		if (buffer==EOF)
		{
			garbage_bits++;
			if (garbage_bits>BITS_IN_REGISTER-2)
			{
				printf("ERROR IN COMPRESSED FILE !!! \n");
				exit(-1);
			}
		}
		bits_to_go=8;
	}
	t=buffer&1;
	buffer>>=1;
	bits_to_go--;
	return t;
}

void start_outputing_bits(void)
{
	buffer=0;
	bits_to_go=8;
}

void output_bit(int bit)
{
	buffer>>=1;
	if (bit) buffer|=0x80;
	bits_to_go--;
	if (bits_to_go==0)
	{
		putc(buffer,out);
		bits_to_go=8;
	}
}

void done_outputing_bits(void)
{
	putc(buffer>>bits_to_go,out);
}

void output_bit_plus_follow(int bit)
{
	output_bit(bit);
	while (bits_to_follow>0)
	{
		output_bit(!bit);
		bits_to_follow--;
	}
}

void start_encoding(void)
{
	low     =0l;
	high    =TOP_VALUE;
	bits_to_follow=0l;
}

void done_encoding(void)
{
	bits_to_follow++;
	if ( low < FIRST_QTR ) output_bit_plus_follow(0);
	else output_bit_plus_follow(1);
}

void start_decoding(void)
{
	int i;
	value=0l;
	for (i=1; i<=BITS_IN_REGISTER; i++)
		value=(value<<1) +input_bit();
	low=0l;
	high=TOP_VALUE;
}

void encode_symbol(int symbol)
{
	unsigned long range;
	range=high-low+1;
	if (symbol!=1) high=low+(range*cum_freq[symbol-1])/cum_freq[0]-1;
	low+=range*cum_freq[symbol]/cum_freq[0];
	for (;;)
	{
		if (high<HALF)
			output_bit_plus_follow(0);
		else if (low>=HALF)
		{
			output_bit_plus_follow(1);
			low-=HALF;
			high-=HALF;
		}
		else if (low>=FIRST_QTR && high < THIRD_QTR)
		{
			bits_to_follow++;
			low-=FIRST_QTR;
			high-=FIRST_QTR;
		}
		else break;
		low<<=1;
		(high<<=1)++;
	}
}

int decode_symbol(void)
{
	unsigned long range, cum;
	int symbol;
	range=high-low+1;
	cum=((value-low+1)*cum_freq[0]-1)/range;
	for (symbol=1; cum_freq[symbol]>cum; symbol++);
	if (symbol!=1) high=low +(range*cum_freq[symbol-1])/cum_freq[0] -1;
	low+=range*cum_freq[symbol]/cum_freq[0];
	for (;;)
	{
		if (high<HALF) {}
		else if (low>=HALF)
		{
			value-=HALF;
			low-=HALF;
			high-=HALF;
		}
		else if (low>=FIRST_QTR && high<THIRD_QTR)
		{
			value-=FIRST_QTR;
			low-=FIRST_QTR;
			high-=FIRST_QTR;
		}
		else break;
		low<<=1;
		(high<<=1)++;
		value=(value<<1)+input_bit();
	}
	return symbol;
}

void encode(char *infile, char *outfile)
{
	int ch,symbol;
	fopen_s(&in,infile,"rb+");
	fopen_s(&out,outfile,"wb+");
	if (in==NULL || out==NULL) return;
	start_model();
	start_outputing_bits();
	start_encoding();
	while ( (ch=getc(in))!=EOF )
	{
		symbol=char_to_index[ch];
		encode_symbol(symbol);
		update_model(symbol);
	}
	encode_symbol(EOF_SYMBOL);
	done_encoding();
	done_outputing_bits();
	fclose(in);
	fclose(out);
}

void decode(char *infile, char *outfile)
{
	int ch,symbol;
	fopen_s(&in,infile,"r+b");
	fopen_s(&out,outfile,"w+b");
	if (in==NULL || out==NULL) return;
	start_model();
	start_inputing_bits();
	start_decoding();
	while ((symbol=decode_symbol())!=EOF_SYMBOL)
	{
		ch=index_to_char[symbol];
		putc(ch,out);
		update_model(symbol);
	}
	fclose(in);
	fclose(out);
}

void initArCoder(int Nlevels)
{
  if(Initialized)
  {
    freeArCoder();
  }
	NO_OF_CHARS = Nlevels;
	EOF_SYMBOL=NO_OF_CHARS+1;
	NO_OF_SYMBOLS = (NO_OF_CHARS+1);
  index_to_char = new unsigned char[NO_OF_SYMBOLS];
	char_to_index =new int[NO_OF_CHARS];
	cum_freq = new unsigned int[NO_OF_SYMBOLS+1];
	freq = new unsigned int[NO_OF_SYMBOLS+1];
  Initialized = true;
}

void freeArCoder()
{
  delete [] index_to_char;
	delete [] char_to_index;
	delete [] cum_freq;
	delete [] freq;
}
/*
void main(int argc, char **argv)
{
	if (argc<5)
	{
		printf("\nUsage: ar e|d infile outfile \n");
		exit(1);
	}
	Nlevels = atoi(argv[4]);
	NO_OF_CHARS = Nlevels;
	EOF_SYMBOL=NO_OF_CHARS+1;
	NO_OF_SYMBOLS = (NO_OF_CHARS+1);
	index_to_char = new unsigned char[NO_OF_SYMBOLS];
	char_to_index =new int[NO_OF_CHARS];
	cum_freq = new unsigned int[NO_OF_SYMBOLS+1];
	freq = new unsigned int[NO_OF_SYMBOLS+1];
	if (argv[1][0]=='e')
		encode(argv[2],argv[3]);
	else if (argv[1][0]=='d')
		decode(argv[2],argv[3]);

	delete [] index_to_char;
	delete [] char_to_index;
	delete [] cum_freq;
	delete [] freq;

	exit(0);
}
*/