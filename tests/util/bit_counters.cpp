/*!
 *	\file bit_counters.cpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019-2020)
 */
#include "bit_counters.hpp"
#include <cstdlib>





/*
 * Allocate memory for counters
 */
TsTestsReturnCode init_bit_counters(BitCounters target, uint64_t const num_of_counters, uint64_t const size_of_counter)
{
	/*
	 * We define the size of the array of bit counters by the following scheme.
	 * <1>. (num_of_counters * size_of_counter) accounts for the total amount of bits that's needed.
	 * <2>. To identify the number of needed bytes we divide this number by 8: <1> >> 3.
	 * <3>. If 8 doesn't divide <1>, then we need to add an extra byte since not all of the bits managed
	 *      to fit into <2> bytes.
	 */
	uint64_t container_size = ((num_of_counters * size_of_counter) >> 3) + ((num_of_counters * size_of_counter) % 8 != 0);
	
	target->amount_of_counters = num_of_counters;
	target->counter_size = size_of_counter;
	target->counters = (uint8_t *) malloc(container_size * sizeof(uint8_t));
	if (target->counters == NULL)
		return TSTESTS_RETURNCODE_FAIL_MEMORY;
	
	for (uint64_t i = 0; i < container_size; ++i)
		target->counters[i] = 0;
	
	return TSTESTS_RETURNCODE_SUCCESS;
}

/*
 * Set bit_i-th bit in target to value
 */
void set_bit(BitCounters *target, uint64_t const bit_i, uint8_t const value)
{
	if (value)
		target->counters[bit_i / 8] |=   1U << (bit_i % 8);
	else
		target->counters[bit_i / 8] &= ~(1U << (bit_i % 8));
	
	return;
}

/*
 * Get bit_i-th bit in target
 */
uint8_t const get_bit(BitCounters const *target, uint64_t const bit_i)
{
	return (target->counters[bit_i / 8] & (1U << (bit_i % 8))) >> (bit_i % 8);
}

/*
 * Increment i-th counter
 */
void increment_counter(BitCounters *target, uint64_t const i)
{
	uint64_t bit_index = i * target->counter_size;
	
	while (1)
	{
		if (bit_index >= (i + 1) * target->counter_size)
			break;
		if (get_bit(target, bit_index) == 0)
		{
			set_bit(target, bit_index, 1);
			for (uint64_t bit_i = i * target->counter_size; bit_i < bit_index; ++bit_i)
				set_bit(target, bit_i, 0);
			break;
		}
		else
			++bit_index;
	}
	
	return;
}

/*
 * Get numeric value of an i-th counter
 */
uint64_t const get_counter(BitCounters const *target, uint64_t const i)
{
	uint64_t result = 0;
	for (uint64_t j = 0; j < target->counter_size; ++j)
	{
		result |= ((uint64_t) get_bit(target, i * target->counter_size + j)) << j;
	}
	
	return result;
}
 
/*
 * Check if the i-th counter equals to desired_value
 */
uint8_t const verify_counter(BitCounters const *target, uint64_t const i, uint64_t const desired_value)
{
	return get_counter(target, i) == desired_value;
}

/*
 * Check if all counters equal to desired_value
 */
uint8_t const verify_counters(BitCounters const *target, uint64_t const desired_value)
{
	for (uint64_t i = 0; i < target->amount_of_counters; ++i)
		if (!verify_counter(target, i, desired_value))
			return 0;
	
	return 1;
}

/*
 * Destroy unneeded bit counters
 */
void destroy_bit_counters(BitCounters *target)
{
	if (target != NULL)
		free(target->counters);
	
	return;
}
