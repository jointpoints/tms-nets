/*!
 *	\file tstest_principals.hpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2020)
 */
#ifndef _TSTEST_PRINCIPALS_HPP_
#define _TSTEST_PRINCIPALS_HPP_





#include "util/common.hpp"
#include "util/incremental_pca.hpp"





/**
 *  \brief
 *  This test performs the incremental principal component analysis for
 *  the set of generated points.
 *
 *  This test can be used to find the axes in multidimensional space
 *  along which the coordinates of points are variated the most. This is
 *  useful for recognising linear patterns in the mutual disposition of
 *  points.
 *
 *  \warning
 *  The IPCA method used in this function is a modification made for achieving
 *  the most accurate results and optimal performance in **case of (t,m,s)-nets**!
 *  This function will not work properly for any general dataset.
 *
 *  \param[in]  test_info   A valid pointer to \c TsTestsInfo.
 *
 *  \return
 *  \c TSTESTS_RETURNCODE_SUCCESS in case of successful completion of calculations.
 *  \c TSTESTS_RETURNCODE_FAIL_INPUT in case of invalidity of \c test_info
 *  pointer.
 *  \c TSTESTS_RETURNCODE_FAIL_MEMORY in case of dynamic memory allocation
 *  fail.
 */
TsTestsReturnCode const tstest_principals(TsTestsInfo *const test_info)
{
	TSTESTS_TEST_FUNCTION_BEGIN(TSTEST_PRINCIPALS, test_info->log_file)
	
	PUSHLOG_4("Test started.")
	
	TsTestsReturnCode   answer  = TSTESTS_RETURNCODE_SUCCESS;
	uint8_t             s       = 0;
	uint64_t            amount  = 0;
	
	PCAMatrix           points_matrix;
	IncrementalPCAInfo  pca_info;
	
	if (test_info == NULL)
	{
		answer = TSTESTS_RETURNCODE_FAIL_INPUT;
		goto instant_death;
	}
	
	s      = test_info->s;
	amount = 1ULL << test_info->m;
	
	/*
	 * Prepare incremental principal component analyser
	 */
	points_matrix = PCAMatrix(s, PCAVector(s, 0.0));
	init_incremental_pca(&pca_info, s);
	
	/*
	 * If (amount) is less than dimensionality (s), test is inapplicable
	 */
	if (amount < s)
	{
		answer = TSTESTS_RETURNCODE_FAIL_INPUT;
		goto instant_death;
	}
	
	PUSHLOG_4("Processing points...")
	for (uint64_t point_i = 0; point_i < amount; ++point_i)
	{
		// After the following line (point) is expected to be (s)-dimensional
		std::vector<TSTESTS_COORDINATE_TYPE> point = test_info->next_point_getter(point_i);
		for (uint8_t dim_i = 0; dim_i < s; ++dim_i)
		{
			points_matrix[point_i % s][dim_i] = point[dim_i];
		}
		
		if (point_i % s == (uint8_t)(s - 1))
			partial_fit(&pca_info, points_matrix);
		
		if (point_i + 1 == amount >> 2)
			PUSHLOG_4("25% of points processed.")
		if (point_i + 1 == amount >> 1)
			PUSHLOG_4("50% of points processed.")
		if (point_i + 1 == 3 * (amount >> 2))
			PUSHLOG_4("75% of points processed.")
	}
	
	PUSHLOG_4("100% of points checked.")
	
	instant_death:
	
	PUSHLOG_4("Test finished.")
	
	switch (answer)
	{
		case TSTESTS_RETURNCODE_SUCCESS:
			{
				PCAMatrix axes;
				PCAVector variance;
				if (s > 1)
				{
					axes      = get_principal_axes(&pca_info);
					variance  = get_components_normalised_variance(&pca_info);
				}
				else
				{
					axes      = {{1.0}};
					variance  = {1.0};
				}
				PUSHLOG_1  ("+")
				PUSHLOG_2  ("+")
				PUSHLOG_3  ("Answer: POSITIVE.")
				APPENDLOG_3("The following principal axes have been detected:")
#				if TSTESTS_VERBOSITY_LEVEL >= 3
				for (uint8_t axis_i = 0; axis_i < s; ++axis_i)
				{
					fprintf(test_info->log_file, "\t%u)\t", axis_i + 1);
					for (uint8_t dim_i = 0; dim_i < s; ++dim_i)
						fprintf(test_info->log_file, "%Lf\t", axes[axis_i][dim_i]);
					fprintf(test_info->log_file, "\n\t\tExplained variance ratio : %Lf.\n", variance[axis_i]);
				}
				
				uint8_t   main_axes = 0;
				for (TSTESTS_COORDINATE_TYPE accumulated_sum = 0.0; accumulated_sum <= 0.8; ++main_axes)
					accumulated_sum += variance[main_axes];
#				endif // TSTESTS_VERBOSITY_LEVEL
				if (variance.front() <= 0.8)
					APPENDLOGF_3("80%% of the variance is explained by the first %u axes.", main_axes)
				else
					APPENDLOG_3 ("80% of the variance is explained by the first axis alone.")
			}
		case TSTESTS_RETURNCODE_FAIL_GENERAL:
			{
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
	
	TSTESTS_TEST_FUNCTION_END
	
	return answer;
}





#endif // _TSTEST_PRINCIPALS_HPP_
