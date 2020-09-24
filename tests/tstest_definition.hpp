/*!
 *	\file tstest_definition.hpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019-2020)
 */
#ifndef _TSTEST_DEFINITION_HPP_
#define _TSTEST_DEFINITION_HPP_





#include "util/common.hpp"
#include "util/bit_counters.hpp"
#include <cstring>
#include <algorithm>
#include <math.h>





/*
 * Generate consecutive set of d's that are used for construction of elementary intervals
 */
bool const get_next_d_set(uint8_t *d_set, uint8_t const dim, uint8_t const d_sum)
{
	uint8_t    *borders     = new uint8_t[dim + 1];
	bool        success     = false;
	
	borders[0] = 0;
	for (uint8_t dim_i = 1; dim_i < dim + 1; ++dim_i)
		borders[dim_i] = borders[dim_i - 1] + d_set[dim_i - 1] + 1;
	for (uint8_t dim_i = dim - 1; dim_i > 0; --dim_i)
	{
		if (borders[dim_i] <= d_sum + dim_i - 1)
		{
			++borders[dim_i];
			for (uint8_t dim_j = dim_i + 1; dim_j < dim; ++dim_j)
				borders[dim_j] = borders[dim_j - 1] + 1;
			success = true;
			goto total_break;
		}
	}
	
	total_break:
	
	if (success)
	{
		for (uint8_t dim_i = 0; dim_i < dim; ++dim_i)
			d_set[dim_i] = borders[dim_i + 1] - borders[dim_i] - 1;
	}
	
	delete [] borders;
	
	return success;
}





/**
 *  \brief
 *  This test validates the definition of (t, m, s)-net for the generated
 *  set of points.
 *
 *  Validation of definition is performed by calculation of amount of points
 *  within each elementary interval. Counters of points occupy the least possible
 *  amount of memory due to bitwise packaging.
 *
 *  \param[in]  test_info   A valid pointer to \c TsTestsInfo.
 *
 *  \return
 *  \c TSTESTS_RETURNCODE_SUCCESS in case of successful definition assertion.
 *  \c TSTESTS_RETURNCODE_FAIL_GENERAL in case when generated points fail to
 *  meet the requirements of definition by not demonstrating the needed amount
 *  of points within at least one elementary interval.
 *  \c TSTESTS_RETURNCODE_FAIL_INPUT in case of invalidity of \c test_info
 *  pointer or in case when generated points fail to meet the requirements
 *  of definition by improperly set \c t and \c m parameters.
 *  \c TSTESTS_RETURNCODE_FAIL_MEMORY in case of dynamic memory allocation
 *  fail.
 */
TSTESTS_TEST_FUNCTION(tstest_definition)
{
	TSTESTS_TEST_FUNCTION_BEGIN(TSTEST_DEFINITION)
	
	PUSHLOG_4("Test started.")
	
	BitCounters *counters       = NULL;
	
	uint8_t      t              = 0;
	uint8_t      m              = 0;
	uint8_t      s              = 0;
	uint64_t     points_count   = 0;
	uint8_t     *d              = NULL;
	uint64_t    *a              = NULL;
	uint8_t      d_sum          = 0;
	uint64_t     numerator      = 1;
	uint64_t     denominator    = 1;
	
	if (test_info == NULL)
	{
		answer = TSTESTS_RETURNCODE_FAIL_INPUT;
		goto instant_death;
	}
	
	/*
	 * Specify net's parameters
	 */
	t            = test_info->t;
	m            = test_info->m;
	s            = test_info->s;
	points_count = 1ULL << m;
	d            = new uint8_t[s];
	a            = new uint64_t[s];
	d_sum        = m - t;
	
	/*
	 * Check relation between (t) and (m)
	 */
	if (t > m)
	{
		answer = TSTESTS_RETURNCODE_FAIL_INPUT;
		goto instant_death;
	}
	
	// All the following is pointless if we already have TSTESTS_RETURNCODE_FAIL_INPUT
	/*
	 * Setup counters
	 * Their amount is calculated as
	 *   / m - t + s - 1 \    m - t
	 *   |               | * 2
	 *   \     s - 1     /
	 */
	for (uint8_t i = 1; i <= s - 1; ++i)
	{
		numerator   *= m - t + s - i;
		denominator *= i;
	}
	counters = new BitCounters;
	answer = init_bit_counters(/*target          = */counters,
	                           /*num_of_counters = */(numerator / denominator) * (1ULL << d_sum),
	                           /*size_of_counter = */d_sum == 0 ? t + 1 : t + 2);
	if (answer != TSTESTS_RETURNCODE_SUCCESS)
		goto instant_death;
	
	/*
	 * Iterate over points
	 */
	PUSHLOG_4("Checking points...")
	for (uint64_t point_i = 0; point_i < points_count; ++point_i)
	{
		// After the following (point) is expected to be (s)-dimensional
#		ifdef TSTESTS_OPTIMISE_FOR_DIGITAL_NETS
		std::vector<TSTESTS_DIGITAL_TYPE> point = test_info->next_point_getter(point_i);
#		else
		std::vector<TSTESTS_COORDINATE_TYPE> point = test_info->next_point_getter(point_i);
#		endif // TSTESTS_OPTIMISE_FOR_DIGITAL_NETS
		
		memset(d, 0, s * sizeof(uint8_t));
		d[s - 1] = d_sum;
		memset(a, 0, s * sizeof(uint64_t));
		
		/*
		 * Find all elementary intervals in which (point) falls
		 */
		uint64_t    d_sets_counter  = 0;
		// When t == m
		if (d_sum == 0)
			increment_counter(counters, d_sets_counter++);
		// When t < m
		// (when t == m, this loop is skipped right after the first line)
		while (true)
		{
			for (uint8_t dim_i = 0; dim_i < s; ++dim_i)
#				ifdef TSTESTS_OPTIMISE_FOR_DIGITAL_NETS
				a[dim_i] = point[dim_i] >> (test_info->bitwidth - d[dim_i]);
#				else
				a[dim_i] = (uint64_t) std::floor(point[dim_i] * (1 << d[dim_i]));
#				endif // TSTESTS_OPTIMISE_FOR_DIGITAL_NETS
			
			// Find out index of corresponding bit counter
			uint64_t    counter_index   = (1ULL << d_sum) * d_sets_counter;
			uint8_t     exponent        = 0;
			for (int16_t i = s - 1; i >= 0; --i)
			{
				counter_index += a[i] * (1ULL << exponent);
				exponent += d[i];
			}
			
			increment_counter(counters, counter_index);
			
			++d_sets_counter;
			
			if (!get_next_d_set(d, s, d_sum)) break;
		}
		
		if (point_i + 1 == points_count >> 2)
			PUSHLOG_4("25% of points checked.")
		if (point_i + 1 == points_count >> 1)
			PUSHLOG_4("50% of points checked.")
		if (point_i + 1 == 3 * (points_count >> 2))
			PUSHLOG_4("75% of points checked.")
	}
	
	PUSHLOG_4("100% of points checked.")
	
	if (!verify_counters(counters, 1ULL << t))
		answer = TSTESTS_RETURNCODE_FAIL_GENERAL;
	
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
				PUSHLOGF_2  ("+ (%u, %u, %u)", +t, +m, +s)
				PUSHLOG_3   ("Answer: POSITIVE.")
				APPENDLOG_3 ("Verified parameters of (t, m, s)-net:")
				APPENDLOGF_3("t = %u\tm = %u\ts = %u", +t, +m, +s)
				break;
			}
		case TSTESTS_RETURNCODE_FAIL_GENERAL:
			{
				PUSHLOG_1   ("-")
				PUSHLOGF_2  ("- (%u, %u, %u)", +t, +m, +s)
				PUSHLOG_3   ("Answer: NEGATIVE.")
				APPENDLOG_3 ("Unverified parameters of (t, m, s)-net:")
				APPENDLOGF_3("t = %u\tm = %u\ts = %u", +t, +m, +s)
				APPENDLOG_3 ("Test failed at the elementary intervals with parameters:")
#				if TSTESTS_VERBOSITY_LEVEL >= 3
				uint64_t    d_sets_counter              = 0;
				uint64_t    failed_intervals_counter    = 0;
				for (uint8_t i = 0; i < s - 1; ++i)
					d[i] = 0;
				d[s - 1] = d_sum;
				// When t == m
				if (d_sum == 0)
				{
					fprintf(test_info->log_file, "\t1)\ta: ");
					for (uint8_t j = 0; j < s; ++j)
						fprintf(test_info->log_file, "0\t");
					fprintf(test_info->log_file, "\n\t\td: ");
					for (uint8_t j = 0; j < s; ++j)
						fprintf(test_info->log_file, "0\t");
					fprintf(test_info->log_file, "\n");
					APPENDLOGF_3("\tExpected amount of points inside: %llu.", 1ULL << t)
					APPENDLOGF_3("\tActual   amount of points inside: %llu.", get_counter(counters, 0))
				}
				// When t < m
				// (when t == m, this loop is skipped right after the first line)
				while (1)
				{
					for (uint64_t i = 0; i < (1ULL << d_sum); ++i)
					{
						uint64_t    i_copy      = i;
						uint8_t     exponent    = m - t;
						for (uint8_t j = 0; j < s; ++j)
						{
							exponent -= d[j];
							a[j] = i_copy / (1ULL << exponent);
							i_copy %= 1ULL << exponent;
						}
						
						if (!verify_counter(counters, (1ULL << d_sum) * d_sets_counter + i, 1ULL << t))
						{
							fprintf(test_info->log_file, "\t%llu)\ta: ", ++failed_intervals_counter);
							for (uint8_t j = 0; j < s; ++j)
								fprintf(test_info->log_file, "%llu\t", a[j]);
							fprintf(test_info->log_file, "\n\t\td: ");
							for (uint8_t j = 0; j < s; ++j)
								fprintf(test_info->log_file, "%u\t", d[j]);
							fprintf(test_info->log_file, "\n");
							APPENDLOGF_3("\tExpected amount of points inside: %llu.", 1ULL << t)
							APPENDLOGF_3("\tActual   amount of points inside: %llu.", get_counter(counters, (1ULL << d_sum) * d_sets_counter + i))
						}
					}
					
					++d_sets_counter;
					
					if (!get_next_d_set(d, s, d_sum)) break;
				}
#				endif // TSTESTS_VERBOSITY_LEVEL
				break;
			}
		case TSTESTS_RETURNCODE_FAIL_INPUT:
			{
				if (test_info != NULL)
				{
					PUSHLOG_1   ("-")
					PUSHLOGF_2  ("- (%u, %u, %u)", +t, +m, +s)
					PUSHLOG_3   ("Answer: NEGATIVE.")
					APPENDLOG_3 ("Unacceptable parameters of (t, m, s)-net:")
					APPENDLOGF_3("t = %u\tm = %u\ts = %u", +t, +m, +s)
					APPENDLOG_3 ("t must be less than or equal to m.")
				}
				else
				{
					PUSHLOG_1   ("- [!]")
					PUSHLOG_2   ("- [!]")
					PUSHLOG_3   ("Answer: NEGATIVE.")
					APPENDLOG_3 ("Invalid test info.")
				}
				break;
			}
		case TSTESTS_RETURNCODE_FAIL_MEMORY:
			{
				PUSHLOG_1   ("- [!]")
				PUSHLOGF_2  ("- [!]", +t, +m, +s)
				PUSHLOG_3   ("Answer: NEGATIVE.")
				APPENDLOG_3 ("Operating system rejected memory allocation calls.")
				APPENDLOG_3 ("Try reducing values of m or s.")
				break;
			}
	}
	
	delete [] d;
	delete [] a;
	destroy_bit_counters(counters);
	delete counters;
	
	TSTESTS_TEST_FUNCTION_END
}





#endif // _TSTEST_DEFINITION_HPP_
