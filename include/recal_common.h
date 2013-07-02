#ifndef COMMON_H_
#define COMMON_H_

#ifdef __MMX__
#include <mmintrin.h>
#endif

#ifdef __SSE__
#include <xmmintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifndef NULL
	#ifdef _cplusplus
		#define NULL 0
	#else
		#define NULL (void *)0
	#endif
#endif

#endif /* COMMON_H_ */
