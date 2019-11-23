/*!
 *	\file tests.hpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019)
 */
#ifndef NETSTESTS_HPP
#define NETSTESTS_HPP



#include "../niederreiter2.hpp"
#include "tests_routines.hpp"
#include <cmath>



#ifdef ENABLE_LOG
#	define LOG(func, step, outOf, stage, message)           func << " [" << step << '/' << outOf << "] [" << stage << "] " << message << '\n'
#	define PUSHLOG(func, step, outOf, stage, message)       std::cout << LOG(func, step, outOf, stage, message);
#else
#	define PUSHLOG(func, step, outOf, stage, message)       (void)0;
#endif // ENABLE_LOG

using namespace sequences;

template<typename UIntType = uint64_t, unsigned int NBITS = 63>
void nied2_check_definition(Niederreiter<UIntType, NBITS> *generator, uint64_t amount)
{
#	define TEST_NAME            "CHECKDEFINITION"
#	define TEST_SUCCESS         0
#	define TEST_INTERVAL_ERROR  1
#	define TEST_BAD_PARAMS      2
	using GeneratorType = Niederreiter<UIntType, NBITS>;
	
	typename GeneratorType::IntPoint	point;
	uint32_t							point_i = 0;
	IntervalCounters					counters;
	uint8_t                             answer = TEST_SUCCESS;
	
	PUSHLOG(TEST_NAME, ' ', ' ', "MESSAGE", "Test started.")
	
	/*
	 * Specify net's parameters
	 */
	uint32_t    t = generator->get_t();
	uint32_t    m = (uint32_t) std::log2(amount);
	uint32_t    s = generator->get_s();
	uint32_t    *d          = new uint32_t[s];
	uint64_t    *a          = new uint64_t[s];
	uint32_t    d_sum       = m - t;
	
	/*
	 * Check relation between <t> and <m>
	 */
	if(t > m)
	{
		answer = TEST_BAD_PARAMS;
		PUSHLOG(TEST_NAME, ' ', ' ', "MESSAGE", "Test will be skipped.")
	}
	
	// Setup counters
	// Their amount is calculated as
	//   / m - t + s - 1 \    m - t
	//   |               | * 2
	//   \     s - 1     /
	// All the following is pointless if we already have TEST_BAD_PARAMS
	if(answer != TEST_BAD_PARAMS)
	{
		uint64_t    numerator   = 1;
		uint64_t    denominator = 1;
		for(uint64_t i = 1; i <= s - 1; ++i)
		{
			numerator   *= m - t + s - i;
			denominator *= i;
		}
		init_interval_counters
		(
			&counters,
			(numerator / denominator)*(1U << d_sum),
			(d_sum == 0) ? (t + 1) : (t + 2)
		);
		
		/*
		 * Iterate over points
		 */
		for(point_i = 0; point_i < amount; ++point_i)
		{
			point = generator->get_point_int(point_i);
			
			/*
			 * Check current point
			 */
			bool        all_d_sets_generated    = false;
			uint64_t    d_sets_counter          = 0;
			for(uint32_t i = 0; i < s; ++i)
				d[i] = 0;
			// When t == m
			if(d_sum == 0)
			{
				for(uint64_t i = 0; i < s; ++i)
					a[i] = 0;
				
				increment_counter(&counters, s, d_sets_counter, 1, d, a);
				
				++d_sets_counter;
			}
			// When t < m, this loop is skipped right after the first line
			while(true)
			{
				get_next_param_set(d, s, d_sum, &all_d_sets_generated);
				if(all_d_sets_generated) break;
				
				for(uint64_t i = 0; i < s; ++i)
					a[i] = (point[i] >> (NBITS - d[i] - 1)) >> 1;

				increment_counter(&counters, s, d_sets_counter, 1ULL << d_sum, d, a);
				
				++d_sets_counter;
			}
			
			if(point_i + 1 == amount >> 2)
				PUSHLOG(TEST_NAME, point_i + 1, amount, "MESSAGE", "25% of points checked.")
			if(point_i + 1 == amount >> 1)
				PUSHLOG(TEST_NAME, point_i + 1, amount, "MESSAGE", "50% of points checked.")
			if(point_i + 1 == 3 * (amount >> 2))
				PUSHLOG(TEST_NAME, point_i + 1, amount, "MESSAGE", "75% of points checked.")
		}
		
		PUSHLOG(TEST_NAME, amount, amount, "MESSAGE", "100% of points checked.")
		
		if(!verify_counters(&counters, 1ULL << t))
			answer = TEST_INTERVAL_ERROR;
	}
	
	PUSHLOG(TEST_NAME, ' ', ' ', "MESSAGE", "Test finished.")
	
	/*
	 * Print the result
	 */
	switch(answer)
	{
		case TEST_SUCCESS:
			{
				std::cout << "\nCHECKDEFINITION\n\tAnswer: POSITIVE.\n\tParameters of (t, m, s)-net:\n";
				std::cout << "\tt = " << t << "\tm = " << m << "\ts = " << s << "\n\n";
				break;
			}			
		case TEST_INTERVAL_ERROR:
			{
				std::cout << "\nCHECKDEFINITION\n\tAnswer: NEGATIVE.\n\tCalculated parameters of (t, m, s)-net:\n";
				std::cout << "\tt = " << t << "\tm = " << m << "\ts = " << s << "\n";
				std::cout << "\tTest failed at the elementary intervals with parameters:\n";
				uint64_t failed_intervals_counter(0);
				bool all_d_sets_generated(false);
				uint64_t d_sets_counter(0);
				for(uint64_t i = 0; i < s; ++i)
					d[i] = 0;
				// When t == m
				if(d_sum == 0)
				{
					std::cout << '\t' << ++failed_intervals_counter << ")\ta: ";
					for(uint64_t j = 0; j < s; ++j)
						std::cout << 0 << '\t';
					std::cout << "\n\t\td: ";
					for(uint64_t j = 0; j < s; ++j)
						std::cout << 0 << '\t';
					std::cout << "\n\t\tExpected amount of points inside: " << (1 << t) << ".\n";
					std::cout <<   "\t\tActual   amount of points inside: " << get_counter(&counters, 0) << ".\n\n";
				}
				// When t < m, this loop is skipped right after the first line
				while(true)
				{
					get_next_param_set(d, s, d_sum, &all_d_sets_generated);
					if(all_d_sets_generated) break;
					
					for(uint64_t i = 0; i < (1ULL << d_sum); ++i)
					{
						uint64_t i_copy(i);
						uint64_t exponent(m - t);
						for(uint64_t j = 0; j < s; ++j)
						{
							exponent -= d[j];
							a[j] = i_copy / (1ULL << exponent);
							i_copy %= 1ULL << exponent;
						}
						
						if(!verify_counter(&counters, (1ULL << d_sum) * d_sets_counter + i, 1ULL << t))
						{
							std::cout << '\t' << ++failed_intervals_counter << ")\ta: ";
							for(uint64_t j = 0; j < s; ++j)
								std::cout << a[j] << '\t';
							std::cout << "\n\t\td: ";
							for(uint64_t j = 0; j < s; ++j)
								std::cout << d[j] << '\t';
							std::cout << "\n\t\tExpected amount of points inside: " << (1 << t) << ".\n";
							std::cout <<   "\t\tActual   amount of points inside: " << get_counter(&counters, (1ULL << d_sum) * d_sets_counter + i) << ".\n\n";
						}
					}
					
					++d_sets_counter;
				}
				break;
			}
		case TEST_BAD_PARAMS:
			{
				std::cout << "\nCHECKDEFINITION\n\tAnswer: NEGATIVE.\n\tUnacceptable calculated parameters of (t, m, s)-net:\n";
				std::cout << "\tt = " << t << "\tm = " << m << "\ts = " << s << "\n";
				std::cout << "\tt must be less than or equal to m.\n\n";
				break;
			}
	}
	
	delete [] d;
	delete [] a;
	if(answer != TEST_BAD_PARAMS)
		delete [] counters.counters;
	
	return;
#	undef TEST_NAME
#	undef TEST_SUCCESS
#	undef TEST_INTERVAL_ERROR
#	undef TEST_BAD_PARAMS
}



template<typename UIntType = uint64_t, unsigned int NBITS = 63>
void nied2_check_uniqueness(Niederreiter<UIntType, NBITS> *generator, uint64_t amount)
{
#	define TEST_NAME            "CHECKUNIQUENESS"
#	define TEST_SUCCESS         0
#	define TEST_FAIL            1
	using GeneratorType = Niederreiter<UIntType, NBITS>;
	
	typename GeneratorType::IntPoint	point;
	uint32_t							point_i = 0;
	uint8_t                             answer = TEST_SUCCESS;
	
	std::cout << "\nCHECKUNIQUENESS\n\tUniqueness test is under current development.\n\n";
	
	return;
#	undef TEST_NAME
#	undef TEST_SUCCESS
#	undef TEST_FAIL
}



#endif // NETSTESTS_HPP
