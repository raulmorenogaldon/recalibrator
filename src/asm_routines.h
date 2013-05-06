#ifdef ASM_ROUTINES_H

#include <stdio.h>
#include <mmintrin.h> 	//- X86 MMX
#include <xmmintrin.h> 	// - X86 SSE1
#include <emmintrin.h> 	// - X86 SSE2 

//__MMX__ -- X86 MMX
//__SSE__ -- X86 SSE
//__SSE2__ -- X86 SSE2 

//SSE detection
int check_sse();

#endif
