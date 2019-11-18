/*!
 *	\file main.cpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019)
 */



/*
 * Uncomment the following line to enable logging
 */
//#define ENABLE_LOG



#include "tests.hpp"



int main(void)
{
	std::vector<uint32_t> degrees1 = {3, 6, 5, 5};
	std::vector<uint32_t> degrees2 = {3, 6, 5, 5, 3, 4, 3, 2};
	
	
	
	/*                       15
	 * 4-dimensional set of 2   points with the best defect possible
	 * This should be a (t, m, s)-net
	 */
	Nied2Generator<uint64_t, 63> generator(4);
	nied2_check_definition(&generator, 1ULL << 15);
	
	/*                       15
	 * 4-dimensional set of 2   points with the manually specified degrees of polynoms
	 * This should be a (t, m, s)-net
	 */
	generator = Nied2Generator<uint64_t, 63>(degrees1);
	nied2_check_definition(&generator, 1ULL << 15);
	
	
	
	/*                       15
	 * 8-dimensional set of 2   points with the manually specified degrees of polynoms
	 * This should NOT be a (t, m, s)-net (due to very high degrees)
	 */
	generator = Nied2Generator<uint64_t, 63>(degrees2);
	nied2_check_definition(&generator, 1ULL << 15);
	
	/*                       15
	 * 4-dimensional set of 2   + 1 points with the manually specified degrees of polynoms
	 * This should NOT be a (t, m, s)-net (due to incorrect amount of points)
	 */
	generator = Nied2Generator<uint64_t, 63>(degrees1);
	nied2_check_definition(&generator, (1ULL << 15) + 1);
	
	/*                       15
	 * 4-dimensional set of 2   + 1 points with the best defect possible
	 * This should NOT be a (t, m, s)-net (due to incorrect amount of points)
	 */
	generator = Nied2Generator<uint64_t, 63>(4);
	nied2_check_definition(&generator, (1ULL << 15) + 1);
	
	return 0;
}
