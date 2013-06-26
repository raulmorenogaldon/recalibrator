#include "aux_misc.h"

void printf_proc_features()
{
	int x64     = 0;
	int MMX     = 0;
	int SSE     = 0;
	int SSE2    = 0;
	int SSE3    = 0;
	int SSSE3   = 0;
	int SSE41   = 0;
	int SSE42   = 0;
	int SSE4a   = 0;
	int AVX     = 0;
	int XOP     = 0;
	int FMA3    = 0;
	int FMA4    = 0;
	int a,b,c,d;

	/*cpuid(0,a,b,c,d);
	int nIds = a;

	cpuid(0x80000000,a,b,c,d);
	int nExIds = a;

	if (nIds >= 1)
	{
		cpuid(0x00000001,a,b,c,d);
		MMX   = (d & ((int)1 << 23)) != 0;
		SSE   = (d & ((int)1 << 25)) != 0;
		SSE2  = (d & ((int)1 << 26)) != 0;
		SSE3  = (c & ((int)1 <<  0)) != 0;

		SSSE3 = (c & ((int)1 <<  9)) != 0;
		SSE41 = (c & ((int)1 << 19)) != 0;
		SSE42 = (c & ((int)1 << 20)) != 0;

		AVX   = (c & ((int)1 << 28)) != 0;
		FMA3  = (c & ((int)1 << 12)) != 0;
    }

    if (nExIds >= 0x80000001){
		cpuid(0x80000001,a,b,c,d);

		x64   = (d & ((int)1 << 29)) != 0;
		SSE4a = (c & ((int)1 <<  6)) != 0;
		FMA4  = (c & ((int)1 << 16)) != 0;
		XOP   = (c & ((int)1 << 11)) != 0;
	}

	printf("-----------------\nProcessor features:\n");*/

	if(x64)
		printf("64 bits\n");

	printf("Supported: ");

	if(MMX)
		printf("/MMX");

	if(SSE)
		printf("/SSE");

	if(SSE2)
		printf("/SSE2");

	if(SSE3)
		printf("/SSE3");

	printf("\n");
}

void print_binary(unsigned int num)
{
	unsigned int mask;

	/* mask = [0000 0000 1000 0000] */
	mask = 0x80;

	while(mask > 0)
	{
		if((num & mask) == 0 )
			printf("0");
		else
			printf("1");
		mask = mask >> 1 ;  // Right Shift
	}
}
