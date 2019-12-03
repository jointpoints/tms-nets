/*!
 *	\file tests_routines.hpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019)
 */
#ifndef NETSTESTSROUTINES_HPP
#define NETSTESTSROUTINES_HPP

#include <stdint.h>



/*
 * BITCOUNTERS ====================================================================================
 */
typedef struct
{
	uint8_t         *counters;
	uint64_t         amount_of_counters;
	uint64_t         counter_size;
} BitCounters;
void                init_bit_counters   (      BitCounters *target, const uint64_t num_of_counters, const uint64_t size_of_counter);
void                set_bit             (      BitCounters *target, const uint64_t bit_i, const uint8_t value);
const uint8_t       get_bit             (const BitCounters *target, const uint64_t bit_i);
void                increment_counter   (      BitCounters *target, const uint64_t i);
const uint64_t      get_counter         (const BitCounters *target, const uint64_t i);
const bool          verify_counter      (const BitCounters *target, const uint64_t i, const uint64_t desired_value);
const bool          verify_counters     (const BitCounters *target,                   const uint64_t desired_value);
/*
 * ================================================================================================
 */

void get_next_d_set(uint32_t *sequence, const uint64_t dim, const uint64_t d_sum, bool *stop_flag);



#endif // NETSTESTSROUTINES_HPP
