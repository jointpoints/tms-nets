/*!
 *	\file tests_routines.cpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019)
 */
#include "tests_routines.hpp"



/*
 * BITCOUNTERS ====================================================================================
 */

/*
 * Allocate memory for counters
 */
void init_bit_counters(BitCounters *target, const uint64_t num_of_counters, const uint64_t size_of_counter)
{
	uint64_t container_size = (num_of_counters * size_of_counter) / 8 + ( (num_of_counters * size_of_counter) % 8 != 0 );
	
	target->amount_of_counters = num_of_counters;
	target->counter_size = size_of_counter;
	target->counters = new uint8_t[container_size];
	for (uint64_t i = 0; i < container_size; ++i)
		target->counters[i] = 0;
	
	return;
}

/*
 * Set bit_i-th bit in target to value
 */
void set_bit(BitCounters *target, const uint64_t bit_i, const uint8_t value)
{
	if (value)
		target->counters[bit_i / 8] |= 1U << (bit_i % 8);
	else
		target->counters[bit_i / 8] &= ~(1U << (bit_i % 8));
	
	return;
}

/*
 * Get bit_i-th bit in target
 */
const uint8_t get_bit(const BitCounters *target, const uint64_t bit_i)
{
	return (target->counters[bit_i / 8] & (1U << (bit_i % 8))) >> (bit_i % 8);
}

/*
 * Increment i-th counter
 */
void increment_counter(BitCounters *target, const uint64_t i)
{
	uint64_t bit_index = i * target->counter_size;
	
	while (true)
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
const uint64_t get_counter(const BitCounters *target, const uint64_t i)
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
const bool verify_counter(const BitCounters *target, const uint64_t i, const uint64_t desired_value)
{
	return get_counter(target, i) == desired_value;
}

/*
 * Check if all counters equal to desired_value
 */
const bool verify_counters(const BitCounters *target, const uint64_t desired_value)
{
	for (uint64_t i = 0; i < target->amount_of_counters; ++i)
		if (!verify_counter(target, i, desired_value))
			return false;
	
	return true;
}

/*
 * ================================================================================================
 */



/*
 * Generate consecutive set of d's used for construction of elementary intervals
 */
void get_next_d_set(uint32_t *sequence, const uint64_t dim, const uint64_t d_sum, bool *stop_flag)
{
	uint64_t    curr_dim            = 0;
	uint64_t    curr_sum            = 0;
	bool        full_recalculate    = true;
	bool        success             = false;
	
	while (!success)
	{
		if (++sequence[curr_dim] > d_sum)
		{
			if (curr_dim + 1 == dim)
			{
				break;
			}
			sequence[curr_dim] = 0;
			++curr_dim;
			full_recalculate = true;
		}
		else
		{
			curr_dim = 0;
			if (full_recalculate)
			{
				curr_sum = 0;
				for(uint64_t i = 0; i < dim; ++i)
					curr_sum += sequence[i];
				full_recalculate = false;
			}
			else
				++curr_sum;
			if (curr_sum == d_sum)
				success = true;
		}
	}
	
	if (!success)
		*stop_flag = true;
	
	return;
}
