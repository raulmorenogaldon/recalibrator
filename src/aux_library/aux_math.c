#include "aux_math.h"

/**
 * Returns quality from probability.
 */
inline double Qvalue(double P)
{
	#ifdef P_SOLEXA
		return Qsolexa(P);
	#else
		return Qsanger(P);
	#endif
}

/**
 * Return probability from quality.
 */
inline double Pvalue(double Q)
{
	#ifdef P_SOLEXA
		return Psolexa(Q);
	#else
		return Psanger(Q);
	#endif
}

/**
 * Return Solexa quality from probability.
 */
inline double Qsolexa(double p)
{
	return -10.0 * log10(p/(1-p));
}

/**
 * Return probability from Solexa quality.
 */
inline double Psolexa(double Q)
{
	return pow(10.0, -Q/10.0) / (pow(10.0, -Q/10.0) + 1);
}

/**
 * Return Sanger quality from probability.
 */
inline double Qsanger(double p)
{
	return -10.0 * log10(p);
}

/**
 * Return probability from Solexa quality.
 */
inline double Psanger(double Q)
{
	return pow(10.0, -Q/10.0);
}
