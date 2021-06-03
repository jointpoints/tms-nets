/**
 *	@file niederreiter.hpp
 *
 *	@brief Includes the generator of a base 2 digital Niederreiter's (t,m,s)-nets.
 *
 *	@author Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter
 *	@author Previous C++ version by John Burkardt 2003-2009 [https://people.sc.fsu.edu/~jburkardt/cpp_src/niederreiter2/niederreiter2.html]
 * 	@author	Alexey Burimov, 2019-2021
 *
 */

#ifndef TMS_NETS_NIEDERREITER_HPP
#define TMS_NETS_NIEDERREITER_HPP

#include "digital_net.hpp"
#include "details/recseq.hpp"


namespace tms
{
	/** Represents digital classical Niederreiter (t, m, s)-net in base 2 */
	class Niederreiter : public DigitalNet
	{
	public:
		
		Niederreiter(void);
		
		/** Constructs the generator of (t,m,s)-net with specified m, s, and with induced least possible t.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] dim - s parameter of the net
		 *  @param [in] in_parallel - flag, defining the method that will be used to generate irreducible polynomials */
		Niederreiter(BasicInt nbits,
		             BasicInt dim,
		             bool     in_parallel = false);
		
		/** Constructs the generator of (t,m,s)-net with specified m, s degrees of initial irreducible polynomials,
		 *  and (optional, for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] degrees_of_irrpolys - vector of the s degrees of initial irreducible polynomials*/
		Niederreiter(BasicInt                     nbits,
		             std::vector<BasicInt> const &degrees_of_irrpolys);
		
		/** Constructs the generator of (t,m,s)-net with specified m, s degrees of initial irreducible polynomials,
		 *  and (optional, for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] degrees_of_irrpolys - vector of the s degrees of initial irreducible polynomials*/
		Niederreiter(BasicInt                               nbits,
		             std::initializer_list<BasicInt> const &degrees_of_irrpolys);
		
		/** Constructs the generator of (t,m,s)-net with specified m, s initial irreducible polynomials, and (optional,
		 *  for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] irrpolys_coeffs - initializer list of the s initial irreducible polynomials coefficients*/
		Niederreiter(BasicInt                                     nbits,
		             std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs);
		
		/** Constructs the generator of (t,m,s)-net with specified m, s initial irreducible polynomials, and (optional,
		 *  for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] irrpolys_coeffs - initializer list of the s initial irreducible polynomials coefficients*/
		Niederreiter(BasicInt                                               nbits,
		             std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs);
		
		~Niederreiter(void);
		
		/// Retutns t-value of a corresponding digital (t, s)-sequence
		BasicInt t_estimate(void) const;
	

	protected:
		
		/// Vector of s initial irreducible polynomials of a (t,m,s)-net.
		std::vector<Polynomial>	m_irrpolys;
		
		/** */
		Niederreiter(BasicInt                       nbits,
					 BasicInt                       dim,
					 std::vector<GenNum>     const &generating_numbers,
					 Real                           recip,
					 std::vector<Polynomial> const &irrpolys,
					 void           (Niederreiter::*ptr_check)(void const *),
					 void                    const *ptr_arg);
		
		void check_init1(void const *ptr_arg);
		void check_init2(void const *ptr_arg);
		void check_init3(void const *ptr_arg);
		
		/** Initializes (t,m,s)-net direction numbers.*/
		virtual void initialize_generating_numbers(void);
	};
	
};// namespace tms


#endif // #ifndef NIEDERREITER2_HPP
