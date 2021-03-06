/**
 * @file    gf2poly.hpp
 *
 * @author Vadim Piven
 * @author Alexey Burimov
 */

#ifndef TMS_NETS_GF2POLY_HPP
#define TMS_NETS_GF2POLY_HPP

#include "common.hpp"

#include <map>


// coefficient number of polynomial over GF(2) is a such integer number that it's n-th bit is equal to n-th coefficient of polynomial.

/** @namespace tms::gf2poly
 *  @brief Contains specific polynomial-related functions that are necessary for (t,m,s)-net generation */
namespace tms::gf2poly
{
	/** Returns polynomial over GF(2) with the specified coefficients (represents a wrapper above irrpoly::gfpoly constructor).
	 *  @param [in] coeffs - desired polynomial coefficients */
	Polynomial              make_gf2poly(std::vector<uintmax_t> const &coeffs);
	/** Generates vector of first least-degree irreducible polynomials over GF(2).
	 *  @param [in] amount - amount of irreducible polynomials to generate
	 *  @param [in] max_defect - upper limit for sum of degrees of polynomials */
	std::vector<Polynomial> generate_irrpolys(unsigned int const amount, unsigned int const max_defect = ~(unsigned int)(0));
	/** Generates vector of first least-degree irreducible polynomials over GF(2) using multithreading.
	 *  @param [in] amount - amount of irreducible polynomials to generate
	 *  @param [in] max_defect - upper limit for sum of degrees of polynomials */
	std::vector<Polynomial> generate_irrpolys_in_parallel(unsigned int const amount, unsigned int const max_defect = ~(unsigned int)(0));
	/** Generates vector of first irreducible polynomials over GF(2) with specified degrees.
	 *  @param [in] degrees - vector of desired degrees of irreducible polynomials
	 *  @param [in] max_defect - upper limit for sum of degrees of polynomials */
	std::vector<Polynomial> generate_irrpolys_with_degrees(std::vector<unsigned int> const &degrees, unsigned int const max_defect = ~(unsigned int)(0));
	
	/** Generates vector of all first irreducible polynomials over GF(2) with degrees <= degree in "ascending" order
	 *  @param [in] degree - greatest degree of generated polynomials */
	std::vector<Polynomial> generate_irrpolys_until_degree(unsigned int const degree);
	
};


#endif /* gf2poly_hpp */
