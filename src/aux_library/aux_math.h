#ifndef AUX_MATH_H_
#define AUX_MATH_H_

#include <math.h>

#include "common.h"

/**
 * Returns quality from probability.
 */
extern inline double Qvalue(double P);

/**
 * Return probability from quality.
 */
extern inline double Pvalue(double Q);

/**
 * Return Solexa quality from probability.
 */
extern inline double Qsolexa(double p);

/**
 * Return probability from Solexa quality.
 */
extern inline double Psolexa(double Q);

/**
 * Return Sanger quality from probability.
 */
extern inline double Qsanger(double p);

/**
 * Return probability from Solexa quality.
 */
extern inline double Psanger(double Q);


#endif /* AUX_MATH_H_ */
