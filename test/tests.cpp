/*
 * Author:  Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019)
 */
#include "tests.hpp"

#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <iostream>
#include <set>

#if ENABLE_LOG
#	define LOG(func, step, outOf, stage, message)           func << " [" << step << '/' << outOf << "] [" << stage << "] " << message << '\n'
#	define PUSHLOG(func, step, outOf, stage, message)       std::cout << LOG(func, step, outOf, stage, message);
#else
#	define PUSHLOG(func, step, outOf, stage, message)       (void)0;
#endif // ENABLE_LOG



/*template<typename UIntType = uint64_t, unsigned int NBITS = 63>
void nied2_check_definition(Nied2Generator<UIntType, NBITS> *generator)
{
	std::cout << "!\n";
	return;
}
*/
