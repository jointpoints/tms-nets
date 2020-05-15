/*!
 *	\file main.cpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2020)
 */
#include "include/niederreiter2.hpp"



#define TSTESTS_VERBOSITY_LEVEL     4           // Try setting this value to 0, 1, 2, 3 or 4
#define TSTESTS_OPTIMISE_FOR_DIGITAL_NETS

#include "tests/tstest_uniqueness.hpp"
#include "tests/tstest_definition.hpp"
#include "tests/tstest_principals.hpp"



#include <time.h>





int main()
{
	uint32_t const nbits = 15;
	sequences::Niederreiter<uint64_t, nbits> generator(3);
	
	FILE *out = fopen("log.txt", "w");
	
	TsTestsInfo tests_info =
	{
		.t                 = (uint8_t) generator.get_t(),
		.m                 = nbits,
		.s                 = (uint8_t) generator.get_s(),
		.bitwidth          = nbits,
		.next_point_getter = [&generator](uint64_t const point_i){return generator.get_point_real(point_i);},
		.log_file          = out
	};
	
	tstest_uniqueness(&tests_info);
	tstest_definition(&tests_info);
	tstest_principals(&tests_info);
	
	fclose(out);
	
	return 0;
}
