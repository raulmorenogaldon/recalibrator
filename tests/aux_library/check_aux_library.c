#include "check_aux_library.h"


START_TEST (check_aux_bam)
{
	// Not necessary
	// compare_bams_qual(const char* bamPath0, const char* bamPath1, int cycles);

	bam_header_t *header = NULL;
	int num_chrom;

	for(num_chrom = 1; num_chrom < 10; num_chrom++)
	{
		header = (bam_header_t *)malloc(sizeof(bam_header_t));
		init_empty_bam_header(num_chrom, header);

		/* Check output */
		ck_assert(header);//, "create_empty_bam_header returns NULL\n");
		ck_assert(header->n_targets == num_chrom);
		ck_assert(header->target_name);
		ck_assert(header->target_len);
		ck_assert(header->text);
		ck_assert(header->l_text == strlen(header->text));

		free(header);
	}

	num_chrom = 0;
	header = (bam_header_t *)malloc(sizeof(bam_header_t));
	if(!init_empty_bam_header(num_chrom, header))
	{
		ck_abort_msg("create_empty_bam_header must return ERROR");
	}
	free(header);
}
END_TEST

double
test_qvalue(double input)
{
	double ret;
	ret = Qvalue(input);
	#ifdef P_SOLEXA
		if(ret > P_SOLEXA_MAX || ret < P_SOLEXA_MIN)
			ck_abort_msg("Qvalue with P_SOLEXA must return in range [%d,%d]: input = %f output = %f",P_SOLEXA_MIN, P_SOLEXA_MAX,input, ret);
	#else
		if(ret > P_SANGER_MAX || ret < P_SANGER_MIN)
			ck_abort_msg("Qvalue must return in range [%d,%d]: input = %f output = %f",P_SANGER_MIN, P_SANGER_MAX, input, ret);
	#endif
	return ret;
}

double
test_pvalue(double input)
{
	double ret;
	ret = Pvalue(input);
	if(ret > 1.0 || ret < 0.0)
		ck_abort_msg("Pvalue must return in range [0,1]: input = %f output = %f", input, ret);
	return ret;
}

START_TEST (check_aux_math)
{
	double input;
	char ret;

	input = 0.000000000000000000001;
	test_qvalue(input);

	input = 0.0;
	test_qvalue(input);

	input = 1.0;
	test_qvalue(input);

	input = 1.1;
	test_qvalue(input);

	input = 0.5;
	ret = test_qvalue(input);
	#ifdef P_SOLEXA
		ck_assert(ret == 0);
	#else
		ck_assert(ret == 3);	// 3.01029...
	#endif

	input = 200;
	test_pvalue(input);

	input = -100;
	test_pvalue(input);

	input = 10.0;
	input = test_pvalue(input);
	#ifdef P_SOLEXA
		ck_assert(input > 0.09);
		ck_assert(input < 0.1);
	#else
		ck_assert_msg(input == 0.1,"%f",input);
	#endif

}
END_TEST

START_TEST (check_aux_vector)
{
	uint32_t *uret = NULL;
	double *dret = NULL;
	int i;

	new_vector_uint32(0, 0, &uret);
	if(uret)
		ck_abort_msg("new_vector must return NULL when size is vector is 0");

	/* new_vector*/
	new_vector_uint32(1, 0, &uret);
	ck_assert(uret && uret[0] == 0);
	free(uret);

	new_vector_uint32(100, 4, &uret);
	ck_assert(uret);
	for(i = 0; i < 100; i++)
		ck_assert(uret[i] == 4);

	/* initialize_vector */
	initialize_vector(uret, 100, -20);
	for(i = 0; i < 100; i++)
		ck_assert(uret[i] == -20);
	free(uret);

	/*new_vector_d */
	new_vector_double(0, 0, &dret);
	if(dret)
		ck_abort_msg("new_vector_d must return NULL when size is vector is 0");
	free(dret);

	new_vector_double(1, 0, &dret);
	ck_assert(dret && dret[0] == 0);
	free(dret);

	new_vector_double(100, 4, &dret);
	ck_assert(dret);
	for(i = 0; i < 100; i++)
		ck_assert(dret[i] == 4);
	free(dret);
}
END_TEST

START_TEST (check_timestats)
{
	double mean, mean2;
	p_timestats timestats = NULL;

	/* Check timestats creation */
	ck_assert(!TIME_GLOBAL_STATS);
	ck_assert(!timestats);
	time_new_stats(20, &timestats);
	ck_assert(timestats);
	ck_assert(TIME_GLOBAL_STATS);

	/* Check timing */

	/* Fails */
	if(!time_init_slot(-1, 20, TIME_GLOBAL_STATS))
	{
		ck_abort_msg("Must return error when slot is major than 20");
	}
	if(!time_init_slot(20, 20, TIME_GLOBAL_STATS))
	{
		ck_abort_msg("Must return error when slot is major than 20");
	}
	if(!time_set_slot(-1, 20, TIME_GLOBAL_STATS))
	{
		ck_abort_msg("Must return error when slot is major than 20");
	}
	if(!time_set_slot(20, 20, TIME_GLOBAL_STATS))
	{
		ck_abort_msg("Must return error when slot is major than 20");
	}

	ck_assert(time_init_slot(0, 0, TIME_GLOBAL_STATS) == NO_ERROR);
	ck_assert(time_set_slot(0, 20, TIME_GLOBAL_STATS) == NO_ERROR);

	ck_assert(time_get_mean_slot(0, TIME_GLOBAL_STATS, &mean) == NO_ERROR);
	ck_assert(time_get_mean_slot(0, timestats, &mean2) == NO_ERROR);
	ck_assert(mean == mean2);

	ck_assert(time_get_min_slot(0, TIME_GLOBAL_STATS, &mean) == NO_ERROR);
	ck_assert(time_get_min_slot(0, timestats, &mean2) == NO_ERROR);
	ck_assert(mean == mean2);

	ck_assert(time_get_max_slot(0, TIME_GLOBAL_STATS, &mean) == NO_ERROR);
	ck_assert(time_get_max_slot(0, timestats, &mean2) == NO_ERROR);
	ck_assert(mean == mean2);

	/* Check timestats destroy */
	ck_assert(time_destroy_stats(&TIME_GLOBAL_STATS) == NO_ERROR);

	/* Try to destroy null pointer */
	if(time_destroy_stats(NULL) == NO_ERROR)
	{
		ck_abort_msg("Destroy time stats must return error if NULL pointer is passed");
	}

}
END_TEST


Suite *
aux_library_suite(void)
{
	Suite *s = suite_create("aux_library");

	/* Core test case */
	TCase *tc_core = tcase_create("core");
	tcase_add_test(tc_core, check_aux_bam);
	tcase_add_test(tc_core, check_aux_math);
	tcase_add_test(tc_core, check_aux_vector);
	tcase_add_test(tc_core, check_timestats);
	suite_add_tcase(s, tc_core);

	return s;
}

int
main (void)
{
	int number_failed;

	Suite *s = aux_library_suite();
	SRunner *sr = srunner_create(s);
	srunner_run_all(sr, CK_NORMAL);
	number_failed = srunner_ntests_failed(sr);
	srunner_free(sr);
	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
