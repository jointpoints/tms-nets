/*!
 *	\file tests_routines.cpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2019)
 */
#include "tests_routines.hpp"



/*
 * FUNCTIONS FOR INTERVAL_COUNTERS ================================================================
 */

/*
 * Allocate memory for counters
 */
void init_interval_counters(IntervalCounters *target, const uint64_t num_of_counters, const uint64_t size_of_counter)
{
	uint64_t container_size((num_of_counters * size_of_counter) / 8 + ( (num_of_counters * size_of_counter) % 8 != 0 ));
	target->amount_of_counters = num_of_counters;
	target->counter_size = size_of_counter;
	target->counters = new unsigned char[container_size];
	for(uint64_t i = 0; i < container_size; ++i)
		target->counters[i] = 0;
	
	return;
}

/*
 * Set <i>-th bit in <target> to <value>
 */
void set_bit(IntervalCounters *target, const uint64_t i, const unsigned char value)
{
	if(value)
		target->counters[i / 8] |= 1U << (i % 8);
	else
		target->counters[i / 8] &= ~(1U << (i % 8));
	
	return;
}

/*
 * Get <i>-th bit in <target>
 */
const uint8_t get_bit(const IntervalCounters *target, const uint64_t i)
{
	return (target->counters[i / 8] & (1U << (i % 8))) >> (i % 8);
}

/*
 * Increment counter for specified by <a>'s and <d>'s elementary interval
 */
void increment_counter(IntervalCounters *target, const uint64_t dim, const uint64_t d_index, const uint64_t d_step, const uint32_t *d, const uint64_t *a)
{
	uint64_t counter_index(d_step * d_index);
	uint64_t exponent(0);
	for(int64_t i = dim - 1; i >= 0; --i)
	{
		counter_index += a[i] * (1 << exponent);
		exponent += d[i];
	}
	
	uint64_t bit_index(counter_index * target->counter_size);
	while(true)
	{
		if(bit_index >= counter_index * target->counter_size + target->counter_size)
			break;
		if(get_bit(target, bit_index) == 0)
		{
			set_bit(target, bit_index, 1);
			for(uint64_t i = counter_index * target->counter_size; i < bit_index; ++i)
				set_bit(target, i, 0);
			break;
		}
		else
			++bit_index;
	}
	
	return;
}

/*
 * Get numeric value of an <i>-th counter
 */
const uint64_t get_counter(const IntervalCounters *target, const uint64_t i)
{
	uint64_t result(0);
	for(uint64_t j = 0; j < target->counter_size; ++j)
	{
		result |= (uint64_t)get_bit(target, i * target->counter_size + j) << j;
	}
	
	return result;
}
 
/*
 * Check if the needed amount of points got inside the <i>-th interval
 */
const bool verify_counter(const IntervalCounters *target, const uint64_t i, const uint64_t desired_value)
{
	return get_counter(target, i) == desired_value;
}

/*
 * Check if the needed amount of points got inside all the intervals
 */
const bool verify_counters(const IntervalCounters *target, const uint64_t desired_value)
{
	for(uint64_t i = 0; i < target->amount_of_counters; ++i)
		if(!verify_counter(target, i, desired_value))
			return false;
	
	return true;
}

/*
 * ================================================================================================
 */



/*
 * Generate consecutive set of d's used for construction of elementary intervals
 */
void get_next_param_set(uint32_t *sequence, const uint64_t dim, const uint64_t d_sum, bool *stop_flag)
{
	uint64_t    curr_dim            (0);
	uint64_t    curr_sum            (0);
	bool        full_recalculate    (true);
	bool        success             (false);
	
	while(!success)
	{
		if(++sequence[curr_dim] > d_sum)
		{
			if(curr_dim + 1 == dim)
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
			if(full_recalculate)
			{
				curr_sum = 0;
				for(uint64_t i = 0; i < dim; ++i)
					curr_sum += sequence[i];
				full_recalculate = false;
			}
			else
				++curr_sum;
			if(curr_sum == d_sum)
				success = true;
		}
	}
	
	if(!success)
		*stop_flag = true;
	
	return;
}
