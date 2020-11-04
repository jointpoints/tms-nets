/*!
 *	\file raref.hpp
 *
 *	\author
 *		Arseny Zakharov (Russian Technological University, KMBO-01-17, Russia, 2020)
 */
#ifndef _RAREF_HPP_
#define _RAREF_HPP_

#include "common.hpp"

typedef std::vector<bool>           RAREFVector;
typedef std::vector<RAREFVector>    RAREFMatrix;

typedef struct RAREF
{
	RAREFMatrix L;
	RAREFMatrix T;
	std::vector<size_t> p;
} RAREF;

RAREF               compute_RAREF(RAREFMatrix const &C);
RAREF               compute_RAREF(RAREFMatrix const &C, RAREFMatrix const &L);
RAREF               update_RAREF (RAREFMatrix const &C, RAREFMatrix const &C2, RAREF const &src);
std::vector<size_t> find_deffect (size_t t, size_t k, size_t s, size_t dmax,
                                  std::vector<RAREFMatrix> const &gen_mat);

#endif