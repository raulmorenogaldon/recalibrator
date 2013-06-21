#ifndef AUX_MATH_H_
#define AUX_MATH_H_

#include <math.h>

/**
 * Returns quality from probability.
 */
inline double Qvalue(double P);

/**
 * Return probability from quality.
 */
inline double Pvalue(double Q);

/**
 * Return Solexa quality from probability.
 */
inline double Qsolexa(double p);

/**
 * Return probability from Solexa quality.
 */
inline double Psolexa(double Q);

/**
 * Return Sanger quality from probability.
 */
inline double Qsanger(double p);

/**
 * Return probability from Solexa quality.
 */
inline double Psanger(double Q);


#endif /* AUX_MATH_H_ */
