#ifndef TMS_NETS_RECSEQ_HPP
#define TMS_NETS_RECSEQ_HPP

#include "common.hpp"


namespace tms::recseq
{
	
	uintmax_t             pack_initial_values(std::vector<BasicInt> const &init_values);
	
	std::vector<BasicInt> unpack_initial_values(uintmax_t packed_init_values, BasicInt count);
	
	void fill_vector_recursively(std::vector<BasicInt> &container, uintmax_t init_values, Polynomial const &char_poly);
	
	uintmax_t  initial_poly_to_initial_values(Polynomial const &init_poly, Polynomial const &char_poly);
	
	Polynomial initial_values_to_initial_poly(uintmax_t init_values, Polynomial const &char_poly);
	
}

#endif
