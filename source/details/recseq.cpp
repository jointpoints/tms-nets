#include "../../include/tms-nets/details/recseq.hpp"
#include "../../include/tms-nets/details/gf2poly.hpp"





uintmax_t tms::recseq::pack_initial_values(std::vector<BasicInt> const &init_values)
{
	uintmax_t packed_init_values = 0;
	
	for (BasicInt i = 0; i < init_values.size() && i < max_nbits; ++i)
	{
		packed_init_values |= init_values[i] << i;
	}
	
	return packed_init_values;
}

std::vector<tms::BasicInt> tms::recseq::unpack_initial_values(uintmax_t packed_init_values, BasicInt count)
{
	std::vector<BasicInt> unpacked_init_values(count);
	
	count != 0 ?
		static_cast<void>(unpacked_init_values[0] = packed_init_values & 1), --count :
		0;
	
	while ( count != 0 )
	{
		*(unpacked_init_values.end() - count) = (packed_init_values >>= 1) & 1;
		--count;
	}
	
	return unpacked_init_values;
}

void tms::recseq::fill_vector_recursively(std::vector<BasicInt> &container, uintmax_t init_values, Polynomial const &char_poly)
{
	if ( char_poly.size() != 1 )
	{
		CountInt seq_i = 0;
		
		CountInt const deg = static_cast<CountInt>(char_poly.degree());
		
		while ( seq_i < deg && seq_i < container.size() )
		{
			container[seq_i] = init_values & 1;
			init_values >>= 1;
			++seq_i;
		}
		//i == deg
		
		while ( seq_i < container.size() )
		{
			container[seq_i] = 0;
			for (BasicInt poly_i = 0; poly_i < deg; ++poly_i)
			{
				container[seq_i] ^= char_poly[poly_i] & container[seq_i - deg + poly_i];
			}
			++seq_i;
		}
	}
	else
	{
		throw std::logic_error("Constant polynomial can't be a characteristic polynomial\n");
	}
}

uintmax_t tms::recseq::initial_poly_to_initial_values(Polynomial const &init_poly, Polynomial const &char_poly)
{
	CountInt const deg = static_cast<CountInt>(char_poly.degree());
	
	uintmax_t init_values = 0;
	
	std::vector<BasicInt> section((deg << 1) - 1);
	fill_vector_recursively(section, 1ULL << (deg - 1), char_poly);
	
	for (BasicInt i = 0; i < init_poly.size(); ++i)
	{
		for (BasicInt j = i; j < init_poly.size(); ++j)
		{
			init_values ^= (section[deg - 1 - i + j] & init_poly[j]) << (deg - 1 - i);
		}
	}
	
	return init_values;
}

tms::Polynomial tms::recseq::initial_values_to_initial_poly(uintmax_t init_values, Polynomial const &char_poly)
{
	BasicInt const deg = static_cast<BasicInt>(char_poly.degree());
	
	std::vector<uintmax_t> init_poly_coeffs(deg, 0);
	
	for (BasicInt i = 0; i < deg; ++i)
	{
		for (BasicInt j = i; j < deg; ++j)
		{
			init_poly_coeffs[i] ^= ((init_values >> (deg - 1 - j)) & 1) & char_poly[deg + i - j];
		}
	}
	
	return gf2poly::make_gf2poly(init_poly_coeffs);
}
