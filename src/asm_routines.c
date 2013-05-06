#include "asm_routines.h"


//SSE detection
int check_sse()
{
	__asm
	(
		"mov EAX,1"
		"cpuid"
		"test EDX, 4000000h" /*If CPUID.01H:EDX.SSE[bit 25] = 1, SSE extensions are present.*/
		"jnz yes"
	);
	
	return 0;
	
	yes: 
	
	return 1;
}
