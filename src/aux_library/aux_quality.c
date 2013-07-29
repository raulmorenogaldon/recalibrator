#include "aux_quality.h"

/***************************
 * QUALITY OPERATIONS
 **************************/

ERROR_CODE
recal_get_estimated_Q(uint32_t *v_bases, uint32_t count, uint8_t start_quality, double *estimated_Q)
{
	double err0;
	double total_sum;
	uint32_t total_bases;
	double aux;
	uint8_t quality;
	int i;

	total_sum = 0.0;
	total_bases = 0.0;

	for(i = count-1; i >= 0; i--)
	{
		quality = (uint8_t)i + start_quality;
		err0 = (double)v_bases[i] * pow(10.0, (-((double)quality)*0.1));
		total_sum += err0;
		total_bases = total_bases + v_bases[i];
	}

	//Calc estimated Q
	aux = total_sum / (double)total_bases;

	aux = Qvalue(aux);
	*estimated_Q = aux;

	return NO_ERROR;
}

ERROR_CODE
recal_get_empirical_Q(double miss, uint32_t bases, double initial_quality, double *emp_qual)
{
	uint32_t mismatches;
	uint32_t observations;
	double log10[91];
	double norm[91];
	double sum_norm;
	uint16_t mle;
	double qempprior;
	double qemplike;
	double max;
	int i;

	mismatches = (uint32_t)(miss + 0.5) + SMOOTH_CONSTANT;
	observations = bases + (SMOOTH_CONSTANT*2);

	//Get logarithms
	for(i = 90; i >= 0; i--)
	{
		log10_Qemp_Reported((double)i, initial_quality, &qempprior);
		log10_Qemp_likelihood((double)i, observations, mismatches, &qemplike);
		log10[i] =  qempprior + qemplike;
	}

	//Normalize
	sum_norm = 0.0;
	max_value(log10, 91, &max);
	for(i = 90; i >= 0; i--)
	{
		norm[i] = pow(10, log10[i] - max);
		sum_norm += norm[i];
	}
	for(i = 90; i >= 0; i--)
	{
		norm[i] = norm[i] / sum_norm;
	}

	//Get maximum log index
	max_index(norm, 91, &mle);

	*emp_qual = (double)mle;

	return NO_ERROR;
}

ERROR_CODE
log10_Qemp_Reported(double Qemp, double Qreported, double *log)
{
	uint8_t difference;
	double gaussian;
	double local_log;

	difference = abs((int)(Qemp - Qreported));

	if(difference > 40)
		difference = 40;

	gaussian = gaussian_function(difference, 0.0, 0.9, 0.0, 0.5);

	local_log = log10(gaussian);

	if(isinf(local_log))
		local_log = -DBL_MAX;

	*log = local_log;

	return NO_ERROR;
}

ERROR_CODE
log10_Qemp_likelihood(double Qemp, uint32_t obs, uint32_t err, double *log)
{
	double log10p;
	double log10OneMinusP;

	if(obs == 0)
		return 0.0;

	log10p = -Qemp*0.1;


	log10OneMinusP = log10(1 - pow( 10, log10p) );

	if(isnan(log10OneMinusP) || isinf(log10OneMinusP))
		log10OneMinusP = -DBL_MAX;

	double factobs = log10_gamma(obs);
	double facterr = log10_gamma(err);
	double factdiff = log10_gamma(obs - err);

	*log = factobs - facterr - factdiff
			+ (log10p * err)
			+ (log10OneMinusP * (obs - err));

	return NO_ERROR;
}
