/*!
 *	\file tests_routines.hpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019)
 */
#ifndef NETSTESTSROUTINES_HPP
#define NETSTESTSROUTINES_HPP

#include <stdint.h>
#include <set>



typedef struct
{
	unsigned char   *counters;
	uint64_t        amount_of_counters;
	uint64_t        counter_size;
} IntervalCounters;
void                init_interval_counters(IntervalCounters *target, const uint64_t num_of_counters, const uint64_t size_of_counter);
void                set_bit             (      IntervalCounters *target, const uint64_t i, const unsigned char value);
const uint8_t       get_bit             (const IntervalCounters *target, const uint64_t i);
void                increment_counter   (      IntervalCounters *target, const uint64_t dim, const uint64_t d_index, const uint64_t d_step, const uint32_t *d, const uint64_t *a);
const uint64_t      get_counter         (const IntervalCounters *target, const uint64_t i);
const bool          verify_counter      (const IntervalCounters *target, const uint64_t i, const uint64_t desired_value);
const bool          verify_counters     (const IntervalCounters *target,                   const uint64_t desired_value);

void get_next_param_set(uint32_t *sequence, const uint64_t dim, const uint64_t d_sum, bool *stop_flag);

void clean_unique_net(std::set<uint64_t *> *unique_net_points);



#endif // NETSTESTSROUTINES_HPP
