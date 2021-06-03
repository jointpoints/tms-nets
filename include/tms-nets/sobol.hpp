#ifndef TMS_NETS_SOBOL_HPP
#define TMS_NETS_SOBOL_HPP

#include "niederreiter.hpp"


namespace tms
{
	/** Represents digital Sobol (t, m, s)-net in base 2 */
	class Sobol : public Niederreiter
	{
	public:
		
		Sobol(void);
		
		Sobol(BasicInt nbits,
			  BasicInt dim,
			  bool     in_parallel = false);
		/** Constructs the generator of (t,m,s)-net with specified m, s degrees of initial irreducible polynomials,
		 *  and (optional, for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] degrees_of_irrpolys - vector of the s degrees of initial irreducible polynomials*/
		Sobol(BasicInt                     nbits,
			  std::vector<BasicInt> const &degrees_of_irrpolys);
		/** Constructs the generator of (t,m,s)-net with specified m, s degrees of initial irreducible polynomials,
		 *  and (optional, for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] degrees_of_irrpolys - vector of the s degrees of initial irreducible polynomials*/
		Sobol(BasicInt                               nbits,
			  std::initializer_list<BasicInt> const &degrees_of_irrpolys);
		/** Constructs the generator of (t,m,s)-net with specified m, s initial irreducible polynomials, and (optional,
		 *  for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] irrpolys_coeffs - initializer list of the s initial irreducible polynomials coefficients*/
		Sobol(BasicInt                                     nbits,
			  std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs);
		/** Constructs the generator of (t,m,s)-net with specified m, s initial irreducible polynomials, and (optional,
		 *  for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] irrpolys_coeffs - initializer list of the s initial irreducible polynomials coefficients*/
		Sobol(BasicInt                                               nbits,
			  std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs);
		
		GenNum inversed_generating_numbers(BasicInt dim) const;
		
		GenMat inversed_generating_matrix(BasicInt dim) const;
		
		
	protected:
		
		void initialize_generating_numbers(void) override;
		
	};
	
	
	
	
	
	
	inline GenMat
	Sobol::inversed_generating_matrix(BasicInt dim) const
	{ return GenMat(inversed_generating_numbers(dim)); }
	
}

#endif
