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
	#ifdef __cplusplus
		#define NULL 0
	#else
		#define NULL (void *)0
	#endif
#endif

#ifdef __cplusplus
	#define EXTERNC extern "C"
#else
	#define EXTERNC extern
#endif

#define BOOL unsigned char
#define TRUE 1
#define FALSE 0

//#define ERROR_CODE unsigned char;
enum ERROR_C {
	NO_ERROR,

	INVALID_INPUT_PARAMS_NULL,
	INVALID_INPUT_PARAMS_0,
	INVALID_INPUT_SIZE_0,

	//Timestats specific
	INVALID_INPUT_SLOT,

	//Recalibrate specific
	INVALID_INPUT_QUAL,
	INVALID_DINUCLEOTIDE
};
typedef enum ERROR_C ERROR_CODE;


#endif /* COMMON_H_ */
