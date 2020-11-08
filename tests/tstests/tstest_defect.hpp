/*!
 *	\file tstest_defect.hpp
 *
 *	\author
 *		Arseny Zakharov (Russian Technological University, KMBO-01-17, Russia, 2020)
 */

#ifndef _TSTEST_DEFECT_HPP_
#define _TSTEST_DEFECT_HPP_

#include "util/raref.hpp"

/**
 *  \brief
 *  This test performs (t,m,s)-net's defect calculation.
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
TSTESTS_TEST_FUNCTION(tstest_defect)
{
	TSTESTS_TEST_FUNCTION_BEGIN(TSTEST_DEFECT)

	uint64_t defect = 0;

	PUSHLOG_4("Test started.")

	if (test_info == NULL)
	{
		answer = TSTESTS_RETURNCODE_FAIL_INPUT;
		goto instant_death;
	}

	answer = find_defect(defect, test_info->m, test_info->s, test_info->gamma_matrix_getter);

	instant_death:

	PUSHLOG_4("Test finished.")

	switch (answer)
	{
		case TSTESTS_RETURNCODE_SUCCESS:
			{
				PUSHLOG_1   ("+")
				PUSHLOGF_2  ("+ The following t was calculated: %u", test_info->m - defect)
				PUSHLOG_3   ("Answer: POSITIVE.")
				APPENDLOGF_3("The following t was calculated: %u", test_info->m - defect)
				break;
			}
		case TSTESTS_RETURNCODE_FAIL_GENERAL:
			{
				PUSHLOG_1  ("-")
				PUSHLOG_2  ("-")
				PUSHLOG_3  ("Answer: NEGATIVE.")
				APPENDLOG_3("Failed to calculate defect.")
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
}

#endif