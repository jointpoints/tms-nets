#ifndef TMS_NETS_MODIFIED_NIEDERREITER_HPP
#define TMS_NETS_MODIFIED_NIEDERREITER_HPP

#include "niederreiter.hpp"

namespace tms
{
	class ModifiedNiederreiter : public Niederreiter
	{
	public:
		
		ModifiedNiederreiter(BasicInt nbits,
							 BasicInt dim,
							 bool     in_parallel = false);
		/** Constructs the generator of (t,m,s)-net with specified m, s degrees of initial irreducible polynomials,
		 *  and (optional, for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] degrees_of_irrpolys - vector of the s degrees of initial irreducible polynomials*/
		ModifiedNiederreiter(BasicInt                     nbits,
							 std::vector<BasicInt> const &degrees_of_irrpolys);
		/** Constructs the generator of (t,m,s)-net with specified m, s degrees of initial irreducible polynomials,
		 *  and (optional, for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] degrees_of_irrpolys - vector of the s degrees of initial irreducible polynomials*/
		ModifiedNiederreiter(BasicInt                               nbits,
							 std::initializer_list<BasicInt> const &degrees_of_irrpolys);
		/** Constructs the generator of (t,m,s)-net with specified m, s initial irreducible polynomials, and (optional,
		 *  for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] irrpolys_coeffs - initializer list of the s initial irreducible polynomials coefficients*/
		ModifiedNiederreiter(BasicInt                                     nbits,
							 std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs);
		/** Constructs the generator of (t,m,s)-net with specified m, s initial irreducible polynomials, and (optional,
		 *  for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] irrpolys_coeffs - initializer list of the s initial irreducible polynomials coefficients*/
		ModifiedNiederreiter(BasicInt                                               nbits,
							 std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs);
		
		DirNum get_inversed_direction_numbers(BasicInt dim) const;
		
		GenMat get_inversed_generating_matrix(BasicInt dim) const;
		
		
	protected:
		
		void initialize_direction_numbers(void) override;
	};
	
	
	
	
	
	
	inline GenMat
	ModifiedNiederreiter::get_inversed_generating_matrix(BasicInt dim) const
	{ return GenMat(get_inversed_direction_numbers(dim)); }
	
}

#endif
