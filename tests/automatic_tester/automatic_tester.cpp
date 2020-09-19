/*!
 *	\file main.cpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2020)
 */
#include "../../include/tms-nets/niederreiter2.hpp"



#define TSTESTS_VERBOSITY_LEVEL     4           // Try setting this value to 0, 1, 2, 3 or 4
#define TSTESTS_OPTIMISE_FOR_DIGITAL_NETS

#include "../tstest_uniqueness.hpp"
#include "../tstest_definition.hpp"
#include "../tstest_principals.hpp"





int main(int argc, char const *const argv[])
{
	int log_manual = -1;
	
	tms::Niederreiter<uint32_t> generator(32, 1);
	
	TsTestsInfo tests_info;
	FILE *out = TSTESTS_LOG_IN_CONSOLE;
	
	uint8_t max_s     = 10;
	std::vector<tms::BasicInt> degrees = {1, 1, 2, 3, 3, 4, 4, 4, 5, 5};
	
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
	
	/*
	 * Perform tests
	 */
	fprintf(out, "* * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
	fprintf(out, "* Automatic tester for Niederreiter sequences       *\n");
	fprintf(out, "* Built on %s %s                     *\n", __DATE__, __TIME__);
	fprintf(out, "* * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
	for (uint8_t curr_s = 1; curr_s <= max_s; ++curr_s)
	{
		fprintf(out, "\n\n\n=== TESTS FOR %u-DIMENSIONAL NETS BEGIN ===\n", curr_s);
		for (uint8_t first_polynomial_i = 0; first_polynomial_i <= max_s - curr_s; ++first_polynomial_i)
		{
			std::vector<tms::BasicInt> curr_degrees(curr_s);
			std::copy(degrees.begin() + first_polynomial_i, degrees.begin() + first_polynomial_i + curr_s, curr_degrees.begin());
			generator = tms::Niederreiter<uint32_t>(32, curr_degrees);
			uint8_t curr_m = 5 * ((generator.get_t() / 5) + 1);
			fprintf(out, "\n\n\nTest case #%u.%u. A priori parameters:\n", curr_s, first_polynomial_i + 1);
			fprintf(out, "\tm = %u\n", curr_m);
			fprintf(out, "\ts = %u\n", curr_s);
			fprintf(out, "Values to be verified:\n");
			fprintf(out, "\tt = %u\n\n", generator.get_t());
			
			tests_info =
			{
				.t                 = (uint8_t) generator.get_t(),
				.m                 = curr_m,
				.s                 = (uint8_t) generator.get_s(),
				.bitwidth          = 32,
				.next_point_getter = [&generator](uint64_t const point_i){return generator.generate_point(point_i);},
				.log_file          = out
			};
			
			tstest_uniqueness(&tests_info);
			tstest_definition(&tests_info);
			tstest_principals(&tests_info);
		}
	}
	
	fclose(out);
	
	instant_death:
	
	return 0;
}
