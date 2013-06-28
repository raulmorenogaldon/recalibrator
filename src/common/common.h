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

#define PNULL (void *)0

#endif /* COMMON_H_ */
