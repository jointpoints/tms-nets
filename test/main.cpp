/*!
 *	\file main.cpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019)
 */



/*
 * Uncomment the following line to enable verbosity
 * Possible levels of verbosity: 1, 2, 3, 4
 */
//#define VERBOSITY_LEVEL 1



#include "tests.hpp"



int main(void)
{
	std::vector<uint32_t> degrees1 = {3, 6, 5, 5, 3};
	std::vector<uint32_t> degrees2 = {3, 6, 5, 5, 3, 4};
	
	using GeneratorType = sequences::Niederreiter<uint32_t, 32>;
	
	/*                       15
	 * 4-dimensional set of 2   points with the best defect possible
	 * This should be a (t, m, s)-net
	 */
	GeneratorType generator(4);
	niederreiter_check_uniqueness(&generator, 15);
	niederreiter_check_definition(&generator, 15);
	std::cout << '\n';
	
	/*                       20
	 * 5-dimensional set of 2   points with the manually specified degrees of polynoms
	 * This should be a (t, m, s)-net
	 */
	generator = GeneratorType(degrees1);
	niederreiter_check_uniqueness(&generator, 20);
	niederreiter_check_definition(&generator, 20);
	std::cout << '\n';
	
	/*                       15
	 * 6-dimensional set of 2   points with the manually specified degrees of polynoms
	 * This should NOT be a (t, m, s)-net (due to very high degrees)
	 */
	generator = GeneratorType(degrees2);
	niederreiter_check_uniqueness(&generator, 15);
	niederreiter_check_definition(&generator, 15);
	std::cout << '\n';
	
	return 0;
}
