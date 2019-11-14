/*!
 *	\file niederreiter2.hpp
 *
 *	\author
 *		Original C++ version by John Burkardt 2003-2009, __ 2019
 */

#ifndef NIEDERREITER2_HPP
#define NIEDERREITER2_HPP

#include "irrpoly/polynomialgf.hpp"


//! Provides an interface to generate an elements of base 2 \f$(t,m,s)\f$-nets in a fixed-dimensional real space.
class Nied2Generator
{
	
public:

	//! Type used for integer values lesser than #c_nbits.
	using 	BasicInt 	= unsigned int;
	//! Type used for quantity-values and indexation-values, i.e. the total number of points to generate.
	using	CountInt	= std::size_t;
	//! Type used to store integer numerators \f$\{Q_n\}\f$ of \f$(t,s)\f$-sequence elements, which bitwitdth defines the accuracy of generation.
	using 	NextInt 	= uint64_t;
	using 	Real		= long double;
	
	using 	IntPoint	= std::vector<NextInt>;
	using 	RealPoint	= std::vector<Real>;
	using 	Polynom 	= irrpoly::polynomialgf<2>;
	
	
	Nied2Generator(void) = delete;
	Nied2Generator(Nied2Generator const &) = delete;
	Nied2Generator& operator =(Nied2Generator &) = delete;
	
	//! Constructor, preparing the generation of uncustomisable dim-dimensional points.
	Nied2Generator(BasicInt const dim, bool const in_parallel = false);
	//! Constructor, preparing the generation of \f$(t, m, s)\f$-net with parameters, specified by vector of degrees of irreducible polynomials.
	Nied2Generator(std::vector<BasicInt> const &degrees_of_irred);
	//! Move constructor.
	Nied2Generator(Nied2Generator &&tmp_generator);
	//! Default destructor.
	~Nied2Generator(void) = default;
	
	//! Move assignment operator.
	Nied2Generator& operator =(Nied2Generator &&tmp_generator) = default;
	
	//! Returns s parameter of possible (t,m,s)-nets.
	BasicInt get_s(void) const;
	//! Returns t parameter of possible (t,m,s)-nets.
	BasicInt get_t(void) const;
	
	//! Loads \f$x_\textbf{pos}\f$ into \b point.
	void load_point_real(RealPoint &point, CountInt const pos);
	//! Loads \f$Q_\textbf{pos}\f$ into \b point.
	void load_point_int (IntPoint  &point, CountInt const pos);
	
	//! Generates \f$x_\textbf{pos}\f$ -- \b pos'th point of (t,s)-sequence.
	RealPoint get_point_real(CountInt const pos);
	//! Generates integer point \f$Q_n = x_\textbf{pos}\cdot 2^{\textbf{c_nbits}}\f$
	 IntPoint get_point_int (CountInt const pos);
	
	//! Loads the \b points vector with sequential integer points \f$\{Q_k\}_{k=\textbf{pos}}^{\textbf{pos}-1+\textbf{points.size()}}\f$.
	void load_points_int (std::vector< IntPoint> &points, CountInt const pos);
	//! Loads the \b points vector with sequential points \f$\{x_k\}_{k=\textbf{pos}}^{\textbf{pos}-1+\textbf{points.size()}}\f$.
	void load_points_real(std::vector<RealPoint> &points, CountInt const pos);
	
	//! Generates \b amount sequential integer points \f$\{Q_k\}_{k=\textbf{pos}}^{n-1+\textbf{amount}}\f$.
	std::vector<RealPoint> get_points_real(CountInt const pos, CountInt const amount);
	//! Generates \b amount sequential points \f$\{x_k\}_{k=\textbf{pos}}^{n-1+\textbf{amount}}\f$.
	std::vector< IntPoint> get_points_int (CountInt const pos, CountInt const amount);

	//! Loads \f$Q_\textbf{pos}\f$ into \b point using the previous point \f$Q_{\textbf{pos}-1}\f$.
	void load_next_int(IntPoint &point, CountInt pos, IntPoint const &prev_point); //MEANT TO BE PRIVATE!!
	
private:
	
	//! Number of available digits to store the integer \f$Q_n = x_n\cdot 2^\textbf{c_nbits+1},\ x_n^{(i)} < 1\ \ \forall i\f$.
	/*! Hence the maximum value of \f$m\f$ parameter of base 2 \f$(t,m,s)\f$-net equals to \b c_nbits as there can be
	 	generated only \f$2^\textbf{c_nbits}\f$ unique numbers. */
	static BasicInt const c_nbits;
	
	static Real		const c_recip; //!< Computational normalization constant equal to \f$2^\textbf{c_nbits}\f$.
		
	BasicInt 			 	c_dim;			//!< Spatial dimensionality, that's equal to s parameter of (t,m,s)-net.
	std::vector<Polynom>	c_irred_polys;	//!< Lookup table of \b c_dim least-degree irreducible polynomials over GF(2).
	BasicInt				c_defect;		//!< t paramteter of (t,m,s)-net.
	std::vector<IntPoint> 	c_cj;			//!< Lookup tables of constants \f$c^{(i)}_{jr}\f$.

//	//! Specifically compares polynomials over GF(2).
//	static bool          lt_poly2(Polynom const &lhs, Polynom const &rhs);
	
	//! Multiplies two polynomials over GF(2).
	static Polynom multiply_poly2(Polynom const &poly_pa, Polynom const &poly_pb);
	
	//!	Generates vector of first \b s_parameter least-degree irreducible polynomials over GF(2).
	std::vector<Polynom> generate_irreducible(std::size_t const s_parameter);
	//!	Generates vector of first \b s_parameter least-degree irreducible polynomials over GF(2) using multythreading.
	std::vector<Polynom> generate_irreducible_in_parallel(std::size_t const amount);
	

	
	//!	Generates vector of first \b s_parameter least-degree irreducible polynomials over GF(2).
	std::vector<Polynom> generate_irreducible_with_degrees(std::vector<BasicInt> const &degrees);
	
	//! Initializes lookup tables of \f$c^{(i)}_{jr}\f$ coefficients.
	void initialize_c(void);
	//! Calculates the value of the constants \f$\{v_n\}\f$.
	void calculate_v(Polynom const &poly_pi, Polynom &poly_b, BasicInt v[]);
	
};


#endif // #ifndef NIEDERREITER2_HPP
