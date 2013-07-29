#include "check_recal.h"


START_TEST (check_empirical_quality)
{
	unsigned int test;
	size_t size = DATA_SIZE;
	double result;

	//Empirical
	printf("==============================\n");
	printf("CHECKING recal_get_empirical_Q\n");
	printf("==============================\n");
	for(test = 0; test < size; test++)
	{
		recal_get_empirical_Q(qual_errors[test], qual_obs[test], (double)quality[test], &result);

		//Is correct?
		printf("Test: %d, Result: %d, Expected: %d\n", test, (unsigned int)result, (unsigned int)qual_empirical[test]);
		ck_assert_msg(result == qual_empirical[test], "Empirical Quality test %d failed.", test);
	}
}
END_TEST

START_TEST (check_estimated_quality)
{
	double estimated_Q;

	//Empirical
	printf("==============================\n");
	printf("CHECKING recal_get_estimated_Q\n");
	printf("==============================\n");

	recal_get_estimated_Q(qual_obs, DATA_SIZE, MINIMUM_QUALITY, &estimated_Q);

	//Is correct?
	printf("Result: %.10f, Expected: %.10f, Diff: %.10f\n", estimated_Q, ESTIMATED_Q, estimated_Q - ESTIMATED_Q);
	ck_assert_msg(fabs(estimated_Q - ESTIMATED_Q) < 0.0001 , "Estimated Quality test failed.");
}
END_TEST

START_TEST (check_recalibration)
{
	double estimated_Q;
	recal_info_t *data;
	unsigned int test;
	unsigned int i;
	uint8_t quality;
	double test_global_emp;
	double data_global_emp;
	double test_estimated;
	double data_estimated;

	//Build data info
	recal_init_info(91, &data);

	//Global
	data->total_bases = GLOBAL_OBS;
	data->total_miss = GLOBAL_MISS;

	//Quality
	for(i = 0; i < DATA_SIZE; i++)
	{
		data->qual_bases[MINIMUM_QUALITY + i] = qual_obs[i];
		data->qual_miss[MINIMUM_QUALITY + i] = qual_errors[i];
	}

	//Calc deltas from data
	recal_calc_deltas(data);

	printf("======================================\n");
	printf("CHECKING DELTA PROCCESS FROM TEST DATA\n");
	printf("======================================\n");

	//Check estimated quality
	recal_get_estimated_Q(qual_obs, DATA_SIZE, MINIMUM_QUALITY, &test_estimated);
	recal_get_estimated_Q(&(data->qual_bases[6]), DATA_SIZE, MINIMUM_QUALITY, &data_estimated);
	printf("Global estimated quality check...\n");
	printf("Result: %.10f, Expected: %.10f, Diff: %.10f\n", data_estimated, test_estimated, data_estimated - test_estimated);
	ck_assert_msg(fabs(data_estimated - test_estimated) < 0.0001 , "Global estimated quality test failed.");

	//Check global empirical
	recal_get_empirical_Q(GLOBAL_MISS, GLOBAL_OBS, test_estimated, &test_global_emp);
	recal_get_empirical_Q(data->total_miss, data->total_bases, data_estimated, &data_global_emp);
	printf("Global empirical check...\n");
	printf("Result: %.10f, Expected: %.10f, Diff: %.10f\n", data_global_emp, test_global_emp, data_global_emp - test_global_emp);
	ck_assert_msg(fabs(data_global_emp - test_global_emp) < 0.0001 , "Global empirical test failed.");

	//Check global delta
	printf("Global Delta check...\n");
	printf("Result: %.10f, Expected: %.10f, Diff: %.10f\n", data->total_delta, GLOBAL_DELTA, data->total_delta - GLOBAL_DELTA);
	ck_assert_msg(fabs(data->total_delta - GLOBAL_DELTA) < 0.0001 , "Global Delta test failed.");

	//Check qualities deltas
	printf("Qualities Delta check...\n");
	for(test = 0; test < DATA_SIZE; test++)
	{
		quality = test + MINIMUM_QUALITY;
		printf("Quality: %d, \tResult: %.2f, \tExpected: %.2f, \t%.2f %d\n",
				quality, data->qual_delta[quality], qual_deltas[test], data->qual_miss[quality], data->qual_bases[quality]);
		ck_assert_msg(fabs(data->qual_delta[quality] - qual_deltas[test]) < 0.0001 , "Quality delta test failed.");
	}

	//recal_get_estimated_Q(obs, 33, 6, &estimated_Q);

	//Is correct?
	//printf("Result: %.10f, Expected: %.10f, Diff: %.10f\n", estimated_Q, expected, estimated_Q - expected);
	//ck_assert_msg(fabs(estimated_Q - expected) < 0.0001 , "Estimated Quality test failed.");

	//Destroy data
	recal_destroy_info(&data);
}
END_TEST


Suite *
recal_suite(void)
{
	Suite *s = suite_create("recal");

	/* Core test case */
	TCase *tc_core = tcase_create("core");
	tcase_add_test(tc_core, check_empirical_quality);
	tcase_add_test(tc_core, check_estimated_quality);
	tcase_add_test(tc_core, check_recalibration);
	//tcase_add_test(tc_core, check_timestats);
	suite_add_tcase(s, tc_core);

	return s;
}

int
main (void)
{
	int number_failed;

	Suite *s = recal_suite();
	SRunner *sr = srunner_create(s);
	srunner_set_log (sr, "test_empirical.log");
	srunner_run_all(sr, CK_NORMAL);
	number_failed = srunner_ntests_failed(sr);
	srunner_free(sr);
	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
