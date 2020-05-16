/*!
 *	\file tstest_uniqueness.hpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019-2020)
 */
#ifndef _TSTEST_UNIQUENESS_HPP_
#define _TSTEST_UNIQUENESS_HPP_





#include "util/common.hpp"
#include "util/bit_counters.hpp"





/**
 *  \brief
 *  This test validates component-wise uniqueness of all generated points.
 *
 *  Validation of component-wise uniqueness is performed using bit arrays:
 *  one bit array per each dimension. This allows reducing memory costs and,
 *  hence, extends the usage of this test on large-scaled data.
 *
 *  \warning
 *  This test is implemented for digital nets only! Its usage without defined
 *  TSTESTS_OPTIMISE_FOR_DIGITAL_NETS macro is impossible.
 *
 *  \param[in]  test_info   A valid pointer to \c TsTestsInfo.
 *
 *  \return
 *  \c TSTESTS_RETURNCODE_SUCCESS in case of successful pass of the test.
 *  \c TSTESTS_RETURNCODE_FAIL_GENERAL in case of existence of at least one
 *  repetetive coordinate for at least one component.
 *  \c TSTESTS_RETURNCODE_FAIL_INPUT in case of invalidity of \c test_info
 *  pointer or in case of an attempt to use this test without defined
 *  TSTESTS_OPTIMISE_FOR_DIGITAL_NETS macro.
 *  \c TSTESTS_RETURNCODE_FAIL_MEMORY in case of dynamic memory allocation
 *  fail.
 */
TsTestsReturnCode const tstest_uniqueness(TsTestsInfo *const test_info)
{
	TSTESTS_TEST_FUNCTION_BEGIN(TSTEST_UNIQUENESS, test_info->log_file)
	
#	ifdef TSTESTS_OPTIMISE_FOR_DIGITAL_NETS
	PUSHLOG_4("Test started.")
	
	TsTestsReturnCode   answer          = TSTESTS_RETURNCODE_SUCCESS;
	uint64_t            unique_points   = 0;
	uint8_t             s               = 0;
	uint64_t            amount          = 0;
	BitCounters        *counters        = NULL;
	
	if (test_info == NULL)
	{
		answer = TSTESTS_RETURNCODE_FAIL_INPUT;
		goto instant_death;
	}
	
	s        = test_info->s;
	amount   = 1ULL << test_info->m;
	counters = new BitCounters[s];
	
	/*
	 * Setup counters for each dimension
	 */
	PUSHLOG_4("Fetching resources...")
	for (uint8_t dim_i = 0; dim_i < s; ++dim_i)
	{
		answer = init_bit_counters(counters + dim_i, 1ULL << test_info->bitwidth, 1);
		if (answer != TSTESTS_RETURNCODE_SUCCESS)
			goto instant_death;
		
		if (dim_i + 1 == s >> 2)
			PUSHLOG_4("25% of rescources fetched.")
		if (dim_i + 1 == s >> 1)
			PUSHLOG_4("50% of rescources fetched.")
		if (dim_i + 1 == 3 * (s >> 2))
			PUSHLOG_4("75% of rescources fetched.")
	}
	PUSHLOG_4("100% of resources fetched.")
	
	/*
	 * Iterate over points
	 */
	PUSHLOG_4("Checking points...")
	for (uint64_t point_i = 0; point_i < amount; ++point_i)
	{
		// After the following line (point) is expected to be (s)-dimensional
		std::vector<TSTESTS_COORDINATE_TYPE>    point_tmp   = test_info->next_point_getter(point_i);
		std::vector<TSTESTS_DIGITAL_TYPE>       point(s, 0);
		std::transform(point_tmp.begin(), point_tmp.end(), point.begin(), [test_info](TSTESTS_COORDINATE_TYPE c){return (TSTESTS_DIGITAL_TYPE)(c * (1ULL << test_info->bitwidth));});
		
		// Check if it's unique; if it is, mark its components as already seen
		uint8_t is_unique = 1;
		for (uint8_t dim_i = 0; dim_i < s; ++dim_i)
		{
			if (verify_counter(counters + dim_i, point[dim_i], 0))
				increment_counter(counters + dim_i, point[dim_i]);
			else
			{
				answer = TSTESTS_RETURNCODE_FAIL_GENERAL;
				is_unique = 0;
			}
		}
		unique_points += is_unique;
		
		if (point_i + 1 == amount >> 2)
			PUSHLOG_4("25% of points checked.")
		if (point_i + 1 == amount >> 1)
			PUSHLOG_4("50% of points checked.")
		if (point_i + 1 == 3 * (amount >> 2))
			PUSHLOG_4("75% of points checked.")
	}
	
	PUSHLOG_4("100% of points checked.")
	
	instant_death:
	
	PUSHLOG_4("Test finished.")
	
	/*
	 * Print the result
	 */
	switch (answer)
	{
		case TSTESTS_RETURNCODE_SUCCESS:
			{
				PUSHLOG_1   ("+")
				PUSHLOGF_2  ("+ (%llu)", unique_points)
				PUSHLOG_3   ("Answer: POSITIVE.")
				APPENDLOGF_3("Component-wise unique points: %llu.", unique_points)
				break;
			}
		case TSTESTS_RETURNCODE_FAIL_GENERAL:
			{
				PUSHLOG_1   ("-")
				PUSHLOGF_2  ("- (%llu)", unique_points)
				PUSHLOG_3   ("Answer: NEGATIVE.")
				APPENDLOGF_3("Expected amount of unique points: %llu.", amount)
				APPENDLOGF_3("Actual   amount of unique points: %llu.", unique_points)
				APPENDLOGF_3("Loss: %llu.", amount - unique_points)
				break;
			}
		case TSTESTS_RETURNCODE_FAIL_INPUT:
			{
				PUSHLOG_1  ("- [!]")
				PUSHLOG_2  ("- [!]")
				PUSHLOG_3  ("Answer: NEGATIVE.")
				APPENDLOG_3("Invalid test info.")
				break;
			}
		case TSTESTS_RETURNCODE_FAIL_MEMORY:
			{
				PUSHLOG_1  ("- [!]")
				PUSHLOG_2  ("- [!]")
				PUSHLOG_3  ("Answer: NEGATIVE.")
				APPENDLOG_3("Operating system rejected memory allocation calls.")
				APPENDLOG_3("Try reducing values of m or s.")
				break;
			}
	}
	
	for (uint8_t i = 0; i < s; ++i)
		destroy_bit_counters(counters + i);
	delete counters;
#	else
	TsTestsReturnCode answer = TSTESTS_RETURNCODE_FAIL_INPUT;
	PUSHLOG_1("Unfortunately, this test is not yet implemented for non-digital nets.")
	PUSHLOG_2("Unfortunately, this test is not yet implemented for non-digital nets.")
	PUSHLOG_3("Unfortunately, this test is not yet implemented for non-digital nets.")
#	endif // TSTESTS_OPTIMISE_FOR_DIGITAL_NETS
	
	TSTESTS_TEST_FUNCTION_END
	
	return answer;
}





#endif // _TSTEST_UNIQUENESS_HPP_
