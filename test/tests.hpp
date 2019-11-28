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



#ifndef VERBOSITY_LEVEL
#	define VERBOSITY_LEVEL              0
#	define VERBOSITY_AUTO_DEFINED
#endif // VERBOSITY_LEVEL

#if VERBOSITY_LEVEL == 1
#	define LOG(func, message)           func << " : " << message << '\n'
#	define PUSHLOG_1(func, message)     std::cout << LOG(func, message);
#else
#	define PUSHLOG_1(func, message)     (void) 0;
#endif // VERBOSITY_LEVEL 1

#if VERBOSITY_LEVEL == 2
#	define LOG(func, message)           func << " : " << message << '\n'
#	define PUSHLOG_2(func, message)     std::cout << LOG(func, message);
#else
#	define PUSHLOG_2(func, message)     (void) 0;
#endif // VERBOSITY_LEVEL 2

#if VERBOSITY_LEVEL >= 3
#	define LOG(func, delim, message)    func << delim << message << '\n'
#	define LOGC(message)                '\t' << message << '\n'
#	define PUSHLOG_3(func, message)     std::cout << LOG(func, "\n\t", message);
#	define PUSHLOGC_3(message)          std::cout << LOGC(message);
#else
#	define PUSHLOG_3(func, message)     (void) 0;
#	define PUSHLOGC_3(message)          (void) 0;
#endif // VERBOSITY_LEVEL 3

#if VERBOSITY_LEVEL >= 4
#	define PUSHLOG_4(func, message)     std::cout << LOG(func, " : ", message);
#else
#	define PUSHLOG_4(func, message)     (void) 0;
#endif // VERBOSITY_LEVEL 4



template<typename UIntType = uint64_t, unsigned int NBITS = 63>
const bool niederreiter_check_definition(sequences::Niederreiter<UIntType, NBITS> *generator, uint8_t m)
{
#	define TEST_NAME            "CHECKDEFINITION"
#	define TEST_SUCCESS         0
#	define TEST_INTERVAL_ERROR  1
#	define TEST_BAD_PARAMS      2
	using GeneratorType = sequences::Niederreiter<UIntType, NBITS>;
	
	uint8_t                             answer = TEST_SUCCESS;
	typename GeneratorType::IntPoint	point;
	uint64_t							point_i = 0;
	BitCounters                         counters;
	
	PUSHLOG_4(TEST_NAME, "Test started.")
	
	/*
	 * Specify net's parameters
	 */
	uint32_t     t          = generator->get_t();
	uint32_t     s          = generator->get_s();
	uint64_t     amount     = 1ULL << m;
	uint32_t    *d          = new uint32_t[s];
	uint64_t    *a          = new uint64_t[s];
	uint32_t     d_sum      = m - t;
	
	/*
	 * Check relation between t and m
	 */
	if (t > m)
	{
		answer = TEST_BAD_PARAMS;
		PUSHLOG_4(TEST_NAME, "Test will be skipped.")
	}
	
	// All the following is pointless if we already have TEST_BAD_PARAMS
	if (answer != TEST_BAD_PARAMS)
	{
		/*
		 * Setup counters
		 * Their amount is calculated as
		 *   / m - t + s - 1 \    m - t
		 *   |               | * 2
		 *   \     s - 1     /
		 */
		uint64_t    numerator   = 1;
		uint64_t    denominator = 1;
		for (uint64_t i = 1; i <= s - 1; ++i)
		{
			numerator   *= m - t + s - i;
			denominator *= i;
		}
		init_bit_counters
		(
			&counters,
			(numerator / denominator) * (1UL << d_sum),
			(d_sum == 0) ? (t + 1) : (t + 2)
		);
		
		/*
		 * Iterate over points
		 */
		for (point_i = 0; point_i < amount; ++point_i)
		{
			point = generator->get_point_int(point_i);
			
			/*
			 * Check current point
			 */
			bool        all_d_sets_generated    = false;
			uint64_t    d_sets_counter          = 0;
			for (uint32_t i = 0; i < s; ++i)
				d[i] = 0;
			// When t == m
			if (d_sum == 0)
			{
				for (uint64_t i = 0; i < s; ++i)
					a[i] = 0;
				
				increment_counter(&counters, d_sets_counter);
				
				++d_sets_counter;
			}
			// When t < m
			// (when t == m, this loop is skipped right after the first line)
			while (true)
			{
				get_next_d_set(d, s, d_sum, &all_d_sets_generated);
				if (all_d_sets_generated) break;
				
				for (uint64_t i = 0; i < s; ++i)
					a[i] = ((uint64_t) point[i]) >> (NBITS - d[i]);
				
				// Find out counter index
				uint64_t counter_index = (1ULL << d_sum) * d_sets_counter;
				uint32_t exponent = 0;
				for (int32_t i = s - 1; i >= 0; --i)
				{
					counter_index += a[i] * (1ULL << exponent);
					exponent += d[i];
				}
				
				increment_counter(&counters, counter_index);
				
				++d_sets_counter;
			}
			
			if (point_i + 1 == amount >> 2)
				PUSHLOG_4(TEST_NAME, "25% of points checked.")
			if (point_i + 1 == amount >> 1)
				PUSHLOG_4(TEST_NAME, "50% of points checked.")
			if (point_i + 1 == 3 * (amount >> 2))
				PUSHLOG_4(TEST_NAME, "75% of points checked.")
		}
		
		PUSHLOG_4(TEST_NAME, "100% of points checked.")
		
		if (!verify_counters(&counters, 1ULL << t))
			answer = TEST_INTERVAL_ERROR;
	}
	
	PUSHLOG_4(TEST_NAME, "Test finished.")
	
	/*
	 * Print the result
	 */
	switch (answer)
	{
		case TEST_SUCCESS:
			{
				PUSHLOG_1(TEST_NAME, "+")
				PUSHLOG_2(TEST_NAME, "+ (" << t << ", " << +m << ", " << s << ")")
				PUSHLOG_3(TEST_NAME, "Answer: POSITIVE.")
				PUSHLOGC_3("Parameters of (t, m, s)-net:")
				PUSHLOGC_3("t = " << t << "\tm = " << +m << "\ts = " << s << '\n')
				break;
			}			
		case TEST_INTERVAL_ERROR:
			{
				PUSHLOG_1(TEST_NAME, "-")
				PUSHLOG_2(TEST_NAME, "- (" << t << ", " << +m << ", " << s << ")")
				PUSHLOG_3(TEST_NAME, "Answer: NEGATIVE.")
				PUSHLOGC_3("Calculated parameters of (t, m, s)-net:")
				PUSHLOGC_3("t = " << t << "\tm = " << +m << "\ts = " << s)
				PUSHLOGC_3("Test failed at the elementary intervals with parameters:")
				bool        all_d_sets_generated        = false;
				uint64_t    d_sets_counter              = 0;
#				if VERBOSITY_LEVEL >= 3
				uint64_t    failed_intervals_counter    = 0;
#				endif // VERBOSITY_LEVEL
				for (uint64_t i = 0; i < s; ++i)
					d[i] = 0;
				// When t == m
				if (d_sum == 0)
				{
#					if VERBOSITY_LEVEL >= 3
					std::cout << "\t1)\ta: ";
					for (uint64_t j = 0; j < s; ++j)
						std::cout << 0 << '\t';
					std::cout << "\n\t\td: ";
					for (uint64_t j = 0; j < s; ++j)
						std::cout << 0 << '\t';
					std::cout << '\n';
#					endif // VERBOSITY_LEVEL
					PUSHLOGC_3("\tExpected amount of points inside: " << (1ULL << t) << ".")
					PUSHLOGC_3("\tActual   amount of points inside: " << get_counter(&counters, 0) << ".\n")
				}
				// When t < m
				// (when t == m, this loop is skipped right after the first line)
				while (true)
				{
					get_next_d_set(d, s, d_sum, &all_d_sets_generated);
					if (all_d_sets_generated) break;
					
					for (uint64_t i = 0; i < (1ULL << d_sum); ++i)
					{
						uint64_t i_copy = i;
						uint64_t exponent = m - t;
						for (uint64_t j = 0; j < s; ++j)
						{
							exponent -= d[j];
							a[j] = i_copy / (1ULL << exponent);
							i_copy %= 1ULL << exponent;
						}
						
						if (!verify_counter(&counters, (1ULL << d_sum) * d_sets_counter + i, 1ULL << t))
						{
#							if VERBOSITY_LEVEL >= 3
							std::cout << '\t' << ++failed_intervals_counter << ")\ta: ";
							for (uint64_t j = 0; j < s; ++j)
								std::cout << a[j] << '\t';
							std::cout << "\n\t\td: ";
							for (uint64_t j = 0; j < s; ++j)
								std::cout << d[j] << '\t';
							std::cout << '\n';
#							endif // VERBOSITY_LEVEL
							PUSHLOGC_3("\tExpected amount of points inside: " << (1ULL << t) << ".")
							PUSHLOGC_3("\tActual   amount of points inside: " << get_counter(&counters, (1ULL << d_sum) * d_sets_counter + i) << ".\n")
						}
					}
					
					++d_sets_counter;
				}
				break;
			}
		case TEST_BAD_PARAMS:
			{
				PUSHLOG_1(TEST_NAME, "-")
				PUSHLOG_2(TEST_NAME, "- (" << t << ", " << +m << ", " << s << ")")
				PUSHLOG_3(TEST_NAME, "Answer: NEGATIVE.")
				PUSHLOGC_3("Unacceptable calculated parameters of (t, m, s)-net:")
				PUSHLOGC_3("t = " << t << "\tm = " << +m << "\ts = " << s)
				PUSHLOGC_3("t must be less than or equal to m.\n")
				break;
			}
	}
	
	delete [] d;
	delete [] a;
	if(answer != TEST_BAD_PARAMS)
		delete [] counters.counters;
	
	return answer == TEST_SUCCESS;
#	undef TEST_NAME
#	undef TEST_SUCCESS
#	undef TEST_INTERVAL_ERROR
#	undef TEST_BAD_PARAMS
}



template<typename UIntType = uint64_t, unsigned int NBITS = 32>
const bool niederreiter_check_uniqueness(sequences::Niederreiter<UIntType, NBITS> *generator, uint8_t m)
{
#	define TEST_NAME            "CHECKUNIQUENESS"
#	define TEST_SUCCESS         0
#	define TEST_FAIL            1
	using GeneratorType = sequences::Niederreiter<UIntType, NBITS>;
	
	uint8_t                              answer         = TEST_SUCCESS;
	typename GeneratorType::IntPoint     point;
	uint64_t                             point_i        = 0;
	uint64_t                             unique_points  = 0;
	uint32_t                             s              = generator->get_s();
	uint64_t                             amount         = 1ULL << m;
	BitCounters                         *counters       = new BitCounters[s];
	
	PUSHLOG_4(TEST_NAME, "Test started.")
	
	/*
	 * Setup counters for each dimension
	 */
	for (uint32_t i = 0; i < s; ++i)
	{
		init_bit_counters(counters + i, 1ULL << generator->get_nbits(), 1);
	}
	
	/*
	 * Iterate over points
	 */
	for (point_i = 0; point_i < amount; ++point_i)
	{
		point = generator->get_point_int(point_i);
		
		// Check if it's unique; if it is, mark its components as already seen
		bool is_unique = true;
		for (uint32_t dim_i = 0; dim_i < s; ++dim_i)
		{
			if (verify_counter(counters + dim_i, point[dim_i], 0))
				increment_counter(counters + dim_i, point[dim_i]);
			else
			{
				answer = TEST_FAIL;
				is_unique = false;
			}
		}
		unique_points = is_unique ? unique_points + 1 : unique_points;
		
		if (point_i + 1 == amount >> 2)
			PUSHLOG_4(TEST_NAME, "25% of points checked.")
		if (point_i + 1 == amount >> 1)
			PUSHLOG_4(TEST_NAME, "50% of points checked.")
		if (point_i + 1 == 3 * (amount >> 2))
			PUSHLOG_4(TEST_NAME, "75% of points checked.")
	}
		
	PUSHLOG_4(TEST_NAME, "100% of points checked.")
	PUSHLOG_4(TEST_NAME, "Test finished.")
	
	/*
	 * Print the result
	 */
	switch (answer)
	{
		case TEST_SUCCESS:
			{
				PUSHLOG_1(TEST_NAME, "+")
				PUSHLOG_2(TEST_NAME, "+ (" << unique_points << ")")
				PUSHLOG_3(TEST_NAME, "Answer: POSITIVE.")
				PUSHLOGC_3("Component-wise unique points: " << unique_points << '\n')
				break;
			}
		case TEST_FAIL:
			{
				PUSHLOG_1(TEST_NAME, "-")
				PUSHLOG_2(TEST_NAME, "+ (" << unique_points << ")")
				PUSHLOG_3(TEST_NAME, "Answer: NEGATIVE.")
				PUSHLOGC_3("Expected amount of unique points: " << amount)
				PUSHLOGC_3("Actual   amount of unique points: " << unique_points)
				PUSHLOGC_3("Loss: " << amount - unique_points << '\n')
				break;
			}
	}
	
	for (uint32_t i = 0; i < s; ++i)
		delete [] counters[i].counters;
	delete [] counters;
	
	return answer == TEST_SUCCESS;
#	undef TEST_NAME
#	undef TEST_SUCCESS
#	undef TEST_FAIL
}



#undef LOG
#undef LOGC
#undef PUSHLOG_1
#undef PUSHLOG_2
#undef PUSHLOG_3
#undef PUSHLOGC_3
#undef PUSHLOG_4
#ifdef VERBOSITY_AUTO_DEFINED
#	undef VERBOSITY_LEVEL
#	undef VERBOSITY_AUTO_DEFINED
#endif // VERBOSITY_AUTO_DEFINED



#endif // NETSTESTS_HPP
