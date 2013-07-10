#ifndef AUX_MATH_H_
#define AUX_MATH_H_

#include <math.h>
#include <recal_common.h>
#include <recal_config.h>

/**
 * Returns quality from probability.
 */
EXTERNC double Qvalue(double P);

/**
 * Return probability from quality.
 */
EXTERNC double Pvalue(double Q);

/**
 * Return Solexa quality from probability.
 */
EXTERNC double Qsolexa(double p);

/**
 * Return probability from Solexa quality.
 */
EXTERNC double Psolexa(double Q);

/**
 * Return Sanger quality from probability.
 */
EXTERNC double Qsanger(double p);

/**
 * Return probability from Solexa quality.
 */
EXTERNC double Psanger(double Q);

#endif /* AUX_MATH_H_ */
