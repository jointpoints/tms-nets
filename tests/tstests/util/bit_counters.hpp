/*!
 *	\file bit_counters.hpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019-2020)
 */
#ifndef BIT_COUNTERS_H
#define BIT_COUNTERS_H



#include "common.hpp"
#include <cstdint>



typedef struct BitCounters
{
	uint8_t         *counters;
	uint64_t         amount_of_counters;
	uint64_t         counter_size;
}
BitCounters;



TsTestsReturnCode   init_bit_counters   (BitCounters       *target, uint64_t const num_of_counters, uint64_t const size_of_counter);
void                set_bit             (BitCounters       *target, uint64_t const bit_i, uint8_t const value);
uint8_t  const      get_bit             (BitCounters const *target, uint64_t const bit_i);
void                increment_counter   (BitCounters       *target, uint64_t const i);
uint64_t const      get_counter         (BitCounters const *target, uint64_t const i);
uint8_t  const      verify_counter      (BitCounters const *target, uint64_t const i, uint64_t const desired_value);
uint8_t  const      verify_counters     (BitCounters const *target,                   uint64_t const desired_value);
void                destroy_bit_counters(BitCounters       *target);



#endif // BIT_COUNTERS_H
