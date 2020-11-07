/*!
 *	\file tstest_defect.hpp
 *
 *	\author
 *		Arseny Zakharov (Russian Technological University, KMBO-01-17, Russia, 2020)
 */

#ifdef TMS_EXPERIMENTAL
#ifndef _TSTEST_DEFECT_HPP_
#define _TSTEST_DEFECT_HPP_

#include "util/raref.hpp"

TSTESTS_TEST_FUNCTION(tstest_defect)
{
	TSTESTS_TEST_FUNCTION_BEGIN(TSTEST_DEFECT)

	uint defect = 0;

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
				PUSHLOG_1("+")
				PUSHLOG_2("+")
				PUSHLOG_3("Answer: POSITIVE.")
				APPENDLOGF_3("The following defect was calculated: %u", defect)
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
#endif