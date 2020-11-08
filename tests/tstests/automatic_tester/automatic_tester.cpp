/*!
 *	\file main.cpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2020)
 */
#include "../../../include/tms-nets/niederreiter2.hpp"



#define TSTESTS_VERBOSITY_LEVEL     1           // Try setting this value to 0, 1, 2, 3 or 4
#define TSTESTS_OPTIMISE_FOR_DIGITAL_NETS

#include "../tstest_uniqueness.hpp"
#include "../tstest_definition.hpp"
#include "../tstest_principals.hpp"
#include "../tstest_defect.hpp"



#define DIGITAL_POINT_GETTER(UIntType, gen)     [&gen](uint64_t const point_i){std::vector<UIntType> point = gen.generate_point_int(point_i); return std::vector<TSTESTS_DIGITAL_TYPE>(point.begin(), point.end());}
#define GAMMA_MATRIX_GETTER(gen)                [&gen](uint64_t const dim){return gen.get_gamma_matrix((tms::BasicInt)dim);}



typedef enum AutoTesterMode
{
	AUTOTESTERMODE_UNDEFINED,
	AUTOTESTERMODE_CRITICAL,
	AUTOTESTERMODE_REGULAR,
	AUTOTESTERMODE_EXHAUSTIVE
}
AutoTesterMode;





void perform_critical_tests(FILE *out, BitCounters *const tests_results)
{
	uint8_t const       num_of_tests = 5;
	init_bit_counters(tests_results, num_of_tests, 1);
	
	TsTestsInfo tests_info;
	
	// 1. One-dimensional nets
	{
		tms::Niederreiter<uint8_t> generator(0, 1);
		fprintf(out, "\n\n\n=== TESTS FOR SINGLE-POINT NETS BEGIN ===\n");
		fprintf(out, "\n\n\nTest case #1.\n\t(0, 0, 1)-net.\n\n");
		
		tests_info =
		{
			.t                 = 0,
			.m                 = 0,
			.s                 = 1,
			.bitwidth          = 0,
			.next_point_getter = DIGITAL_POINT_GETTER(uint8_t, generator),
			.gamma_matrix_getter = GAMMA_MATRIX_GETTER(generator),
			.log_file          = out
		};
		if (tstest_uniqueness(&tests_info) + tstest_definition(&tests_info) + tstest_principals(&tests_info) + tstest_defect(&tests_info))
			increment_counter(tests_results, 0);
	}
	// 2. Boundary values
	{
		tms::Niederreiter<uint8_t>  generator8 (8,  std::vector<tms::BasicInt>{5});
		tms::Niederreiter<uint16_t> generator16(16, std::vector<tms::BasicInt>{5});
		tms::Niederreiter<uint32_t> generator32(32, std::vector<tms::BasicInt>{18});
		tms::Niederreiter<uint64_t> generator64(64, std::vector<tms::BasicInt>{18});
		
		fprintf(out, "\n\n\n=== TESTS FOR BOUNDARY VALUES BEGIN ===\n");
		fprintf(out, "\n\n\nTest case #2.\n\t(4, 8, 1)-net with UIntType = uint8_t.\n\n");
		tests_info =
		{
			.t                 = 4,
			.m                 = 8,
			.s                 = 1,
			.bitwidth          = 8,
			.next_point_getter = DIGITAL_POINT_GETTER(uint8_t, generator8),
			.gamma_matrix_getter = GAMMA_MATRIX_GETTER(generator8),
			.log_file          = out
		};
		if (tstest_uniqueness(&tests_info) + tstest_definition(&tests_info) + tstest_principals(&tests_info) + tstest_defect(&tests_info))
			increment_counter(tests_results, 1);
		
		fprintf(out, "\n\n\nTest case #3.\n\t(4, 16, 1)-net with UIntType = uint16_t.\n\n");
		tests_info =
		{
			.t                 = 4,
			.m                 = 16,
			.s                 = 1,
			.bitwidth          = 16,
			.next_point_getter = DIGITAL_POINT_GETTER(uint16_t, generator16),
			.gamma_matrix_getter = GAMMA_MATRIX_GETTER(generator16),
			.log_file          = out
		};
		if (tstest_uniqueness(&tests_info) + tstest_definition(&tests_info) + tstest_principals(&tests_info) + tstest_defect(&tests_info))
			increment_counter(tests_results, 2);
		
		fprintf(out, "\n\n\nTest case #4.\n\tFirst 2^18 points of (17, 32, 1)-net with UIntType = uint32_t.\n\n");
		tests_info =
		{
			.t                 = 17,
			.m                 = 18,
			.s                 = 1,
			.bitwidth          = 32,
			.next_point_getter = DIGITAL_POINT_GETTER(uint32_t, generator32),
			.gamma_matrix_getter = GAMMA_MATRIX_GETTER(generator32),
			.log_file          = out
		};
		if (tstest_uniqueness(&tests_info) + tstest_definition(&tests_info) + tstest_principals(&tests_info) + tstest_defect(&tests_info))
			increment_counter(tests_results, 3);
		
		fprintf(out, "\n\n\nTest case #5.\n\tFirst 2^18 points of (17, 64, 1)-net with UIntType = uint64_t.\n\tUniqueness is skipped due to memory limitations.\n\n");
		tests_info =
		{
			.t                 = 17,
			.m                 = 18,
			.s                 = 1,
			.bitwidth          = 64,
			.next_point_getter = [&generator64](uint64_t const point_i){return generator64.generate_point_int(point_i);},
			.gamma_matrix_getter = GAMMA_MATRIX_GETTER(generator64),
			.log_file          = out
		};
		// Uniqueness would require 147 billion GB
		if (/*tstest_uniqueness(&tests_info) + */tstest_definition(&tests_info) + tstest_principals(&tests_info) + tstest_defect(&tests_info))
			increment_counter(tests_results, 4);
	}
	
	return;
}





void perform_regular_tests(FILE *out, BitCounters *const tests_results)
{
	uint8_t const       num_of_tests = 45;
	init_bit_counters(tests_results, num_of_tests, 1);
	
	tms::Niederreiter<uint64_t> generator(32, 1);
	
	TsTestsInfo tests_info;
	
	uint8_t const max_s = 9;
	std::vector<tms::BasicInt> degrees = {1, 1, 2, 3, 3, 4, 4, 4, 5, 5};
	
	uint64_t test_i = 0;
	for (uint8_t curr_s = 1; curr_s <= max_s; ++curr_s)
	{
		fprintf(out, "\n\n\n=== TESTS FOR %u-DIMENSIONAL NETS BEGIN ===\n", curr_s);
		for (uint8_t first_polynomial_i = 0; first_polynomial_i <= max_s - curr_s; ++first_polynomial_i)
		{
			std::vector<tms::BasicInt> curr_degrees(curr_s);
			std::copy(degrees.begin() + first_polynomial_i, degrees.begin() + first_polynomial_i + curr_s, curr_degrees.begin());
			uint8_t curr_t = std::accumulate(curr_degrees.begin(), curr_degrees.end(), 0) - curr_s;
			uint8_t curr_m = 5 * ((curr_t / 5) + 1);
			generator = tms::Niederreiter<uint64_t>(curr_m, curr_degrees);
			fprintf(out, "\n\n\nTest case #%llu.\n\t(%u, %u, %u)-net.\n\n", test_i + 1, curr_t, curr_m, curr_s);
			
			tests_info =
			{
				.t                 = curr_t,
				.m                 = curr_m,
				.s                 = curr_s,
				.bitwidth          = curr_m,
				.next_point_getter = [&generator](uint64_t const point_i){return generator.generate_point_int(point_i);},
				.gamma_matrix_getter = GAMMA_MATRIX_GETTER(generator),
				.log_file          = out
			};
			if (tstest_uniqueness(&tests_info) + tstest_definition(&tests_info) + tstest_principals(&tests_info) + tstest_defect(&tests_info))
				increment_counter(tests_results, test_i);
			
			++test_i;
		}
	}
	
	return;
}





void perform_exhaustive_tests(FILE *out, BitCounters *const tests_results)
{
	uint8_t const       num_of_tests = 1;
	init_bit_counters(tests_results, num_of_tests, 1);
	
	TsTestsInfo tests_info;
	
	// 1. Massive data test
	{
		tms::Niederreiter<uint32_t> generator(32, std::vector<tms::BasicInt>{1, 2, 3, 3, 4, 4, 4, 5, 5, 5});
		fprintf(out, "\n\n\n=== TESTS FOR MASSIVE NETS BEGIN ===\n");
		fprintf(out, "\n\n\nTest case #1.\n\t(26, 29, 10)-net.\n\n");
		
		tests_info =
		{
			.t                 = 26,
			.m                 = 29,
			.s                 = 10,
			.bitwidth          = 32,
			.next_point_getter = DIGITAL_POINT_GETTER(uint32_t, generator),
			.gamma_matrix_getter = GAMMA_MATRIX_GETTER(generator),
			.log_file          = out
		};
		if (tstest_uniqueness(&tests_info) + tstest_definition(&tests_info) + tstest_principals(&tests_info) + tstest_defect(&tests_info))
			increment_counter(tests_results, 0);
	}
	
	return;
}





int main(int argc, char const *const argv[])
{
	int                 log_manual                  = -1;
	AutoTesterMode      mode                        = AUTOTESTERMODE_UNDEFINED;
	void (*autotests[])(FILE*, BitCounters *const)  = {NULL, perform_critical_tests, perform_regular_tests, perform_exhaustive_tests};
	char const *const   names_of_modes[]            = {NULL, "\n\n\n                   CRITICAL TESTS\n", "\n\n\n                   REGULAR TESTS\n", "\n\n\n                 EXHAUSTIVE TESTS\n"};
	
	BitCounters    tests_results;
	
	FILE *out = TSTESTS_LOG_IN_CONSOLE;
	
	/*
	 * Interpret command line parameters
	 */
	{
		int arg_i = 1;
		while (arg_i < argc)
		{
			if (!strcmp(argv[arg_i], "-log") && (log_manual == -1))
			{
				if (arg_i + 1 < argc)
				{
					log_manual = arg_i + 1;
					arg_i += 2;
					continue;
				}
				else
				{
					printf("Invalid command line.\n");
					goto instant_death;
				}
			}
			if (!strcmp(argv[arg_i], "-c") && (mode == AUTOTESTERMODE_UNDEFINED))
			{
				mode = AUTOTESTERMODE_CRITICAL;
				++arg_i;
				continue;
			}
			if (!strcmp(argv[arg_i], "-r") && (mode == AUTOTESTERMODE_UNDEFINED))
			{
				mode = AUTOTESTERMODE_REGULAR;
				++arg_i;
				continue;
			}
			if (!strcmp(argv[arg_i], "-e") && (mode == AUTOTESTERMODE_UNDEFINED))
			{
				mode = AUTOTESTERMODE_EXHAUSTIVE;
				++arg_i;
				continue;
			}
			printf("Invalid command line.\n");
			goto instant_death;
		}
	}
	if (log_manual > 0)
	{
		out = fopen(argv[log_manual], "w");
		if (out == NULL)
		{
			printf("File cannot be opened.\n");
			goto instant_death;
		}
	}
	if (mode == AUTOTESTERMODE_UNDEFINED)
		mode = AUTOTESTERMODE_CRITICAL;
	
	/*
	 * Print the header
	 */
	{
		fprintf(out, "* * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
		fprintf(out, "* Automatic tester for digital (t, m, s)-nets       *\n");
		fprintf(out, "* Built on %s %s                     *\n", __DATE__, __TIME__);
		fprintf(out, "* * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
		time_t    start_time;
		struct tm start_ltime;
		time(&start_time);
		start_ltime = *localtime(&start_time);
		fprintf(out, names_of_modes[mode]);
		fprintf(out, "          started on %02d.%02d.%d %02d:%02d:%02d\n\n", start_ltime.tm_mday, start_ltime.tm_mon+1, start_ltime.tm_year+1900, start_ltime.tm_hour, start_ltime.tm_min, start_ltime.tm_sec);
		fprintf(out, "         scroll to the bottom of this log\n");
		fprintf(out, "             to see the brief summary\n");
	}
	
	/*
	 * Perform tests
	 */
	autotests[mode](out, &tests_results);
	
	/*
	 * Print the summary
	 */
	{
		uint64_t failed_tests = 0;
		for (uint64_t counter_i = 0; counter_i < tests_results.amount_of_counters; ++counter_i)
			failed_tests += get_counter(&tests_results, counter_i);
		fprintf(out, "\n\n\n                      SUMMARY\n");
		fprintf(out, "\n\n\nTotal      test cases : %llu\n", tests_results.amount_of_counters);
		fprintf(out, "Successful test cases : %llu\n", tests_results.amount_of_counters - failed_tests);
		fprintf(out, "Failed*    test cases : %llu\n", failed_tests);
		fprintf(out, "Full list of failed* test cases:\n");
		if (failed_tests)
		{
			for (uint64_t counter_i = 0; counter_i < tests_results.amount_of_counters; ++counter_i)
				if (verify_counter(&tests_results, counter_i, 1))
					fprintf(out, "\t%llu\n", counter_i + 1);
		}
		else
			fprintf(out, "\t-----\n");
		fprintf(out, "\n* a test case is considered to be failed if at least\none of the tests in it has failed.\n");
		time_t    end_time;
		struct tm end_ltime;
		time(&end_time);
		end_ltime = *localtime(&end_time);
		fprintf(out, "\n\n\n                    END OF LOG\n          finished on %02d.%02d.%d %02d:%02d:%02d\n\n\n", end_ltime.tm_mday, end_ltime.tm_mon+1, end_ltime.tm_year+1900, end_ltime.tm_hour, end_ltime.tm_min, end_ltime.tm_sec);
	}
	
	instant_death:
		
	fclose(out);
	destroy_bit_counters(&tests_results);
	
	return 0;
}
