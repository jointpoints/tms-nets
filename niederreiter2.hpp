/*!
 *	\file niederreiter2.hpp
 *
 *	\author
 *		Original C++ version by John Burkardt 2003-2009, __ 2019
 */

#ifndef NIEDERREITER2_HPP
#define NIEDERREITER2_HPP

#include "irrpoly/polynomialgf.hpp"
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <map>
#include <type_traits>

//! Provides an interface to generate an elements of base 2 \f$(t,m,s)\f$-nets in a fixed-dimensional real space.
template <
	//! Type used to store integer numerators \f$\{Q_n\}\f$ of \f$(t,s)\f$-sequence elements, which bitwitdth defines the accuracy of generation.
	typename UIntType = uint64_t,
	//! Number of available digits to store the integer \f$Q_n = x_n\cdot 2^\textbf{NBITS+1},\ x_n^{(i)} < 1\ \ \forall i\f$.
	/*! Hence the maximum value of \f$m\f$ parameter of base 2 \f$(t,m,s)\f$-net equals to \b NBITS as there can be
	 generated only \f$2^\textbf{NBITS}\f$ unique numbers. */
	unsigned int NBITS = 63>
class Nied2Generator
{
static_assert(std::is_unsigned<UIntType>::value, "UIntType template parameter must be an unsigned integer type.");
static_assert(NBITS && NBITS <= sizeof(UIntType)*8, "NBITS template parameter must be nonzero <= bitwidth of UIntType.");
	
public:

	//! Type used for integer values lesser than #NBITS.
	using 	BasicInt 	= unsigned int;
	//! Type used for quantity-values and indexation-values, i.e. the total number of points to generate.
	using	CountInt	= std::size_t;
	//! Type used to store integer numerators \f$\{Q_n\}\f$ of \f$(t,s)\f$-sequence elements, which bitwitdth defines the accuracy of generation.
	using 	NextInt 	= UIntType;
	using 	Real		= long double;
	
	using 	IntPoint	= std::vector<NextInt>;
	using 	RealPoint	= std::vector<Real>;
	using 	Polynom 	= irrpoly::polynomialgf<2>;
	
	
	Nied2Generator(void) = delete;
	Nied2Generator(Nied2Generator const &) = delete;
	Nied2Generator& operator =(Nied2Generator const &) = delete;
	
	//! Constructor, preparing the generation of uncustomisable dim-dimensional points.
	Nied2Generator(BasicInt const dim, bool const in_parallel = false);
	//! Constructor, preparing the generation of \f$(t, m, s)\f$-net with parameters, specified by vector of degrees of irreducible polynomials.
	Nied2Generator(std::vector<BasicInt> const &degrees_of_irred);
	//! Default destructor.
	~Nied2Generator(void) = default;
	
	//! Move assignment operator.
	Nied2Generator& operator =(Nied2Generator &&tmp_generator) = default;
	
	//! Returns s parameter of possible (t,m,s)-nets.
	BasicInt get_s(void) const;
	//! Returns t parameter of possible (t,m,s)-nets.
	BasicInt get_t(void) const;
	//! Returns #NBITS -- bitwidth of Q_n integers.
	BasicInt get_nbits(void) const;
	
	//! Loads \f$x_\textbf{pos}\f$ into \b point.
	void load_point_real(RealPoint &point, CountInt const pos);
	//! Loads \f$Q_\textbf{pos}\f$ into \b point.
	void load_point_int (IntPoint  &point, CountInt const pos);
	
	//! Generates \f$x_\textbf{pos}\f$ -- \b pos'th point of (t,s)-sequence.
	RealPoint get_point_real(CountInt const pos);
	//! Generates integer point \f$Q_n = x_\textbf{pos}\cdot 2^{\textbf{NBITS}}\f$
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
	
	static Real		const c_recip; //!< Computational normalization constant equal to \f$2^\textbf{NBITS}\f$.
		
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


//=========================================================================================================
// Implementation:
//=========================================================================================================

template <typename T, unsigned int NBITS>
typename Nied2Generator<T,NBITS>::Real     const Nied2Generator<T,NBITS>::c_recip   = static_cast<Real>(1.0)/pow(2, NBITS);

/*!
 *
 *	\author Vadim Piven
 *
 */
template <typename T, unsigned int NBITS>
std::vector<typename Nied2Generator<T,NBITS>::Polynom> Nied2Generator<T,NBITS>::generate_irreducible_in_parallel(std::size_t const amount)
{
	// возвращаемое значение
	std::vector<Polynom> irred_polys_table;
	if (amount == 0) { throw std::logic_error("Spatial dimensionality can't be equal to 0."); }
	irred_polys_table.reserve(amount);
	irred_polys_table.emplace_back(Polynom({0, 1}));
	if (amount == 1) { return irred_polys_table; }
	
	// создаём всё необходимое для многопоточности
	irrpoly::multithread::polychecker<2> ch;
	
	// функция, генерирующая многочлены для проверки
	auto input = []() -> Polynom {
		static CountInt index = 1;
		Polynom poly;
		BasicInt degree = 0;
		for (BasicInt i = 1; index >> i; ++i, ++degree);
		poly.data().reserve(degree + 2);
		poly.data().emplace_back(1);
		for (BasicInt i = 0; i <= degree; ++i) {
			poly.data().emplace_back((index & (1ull << i)) ? 1 : 0);
		}
		++index;
		return poly;
	};
	
	auto check_func = irrpoly::multithread::make_check_func<2>(
															   irrpoly::multithread::irreducible_method::berlekamp,
															   irrpoly::multithread::primitive_method::nil);
	
	c_defect = 0;
	
	// функция, вызываемая по окончании проверки, если результат нам подходит - сохраняем и возвращаем true, иначе false
	auto callback = [&](const Polynom &poly, const typename irrpoly::multithread::result_type &result) -> bool {
		if (!result.irreducible) { return false; }
		c_defect += poly.degree() - 1;
		irred_polys_table.emplace_back(poly);
		return true;
	};
	
	// запускаем генерацию num - 1 многочленов (т.к. один у нас уже есть)
	// последний false говорит, что нам нужны все результаты проверки, включая лишние, поскольку мы загружаем
	// многочлены на проверку последовательно, но многопоточность не гарантирует строгого порядка, и какие-то
	// многочлены из начала последовательности могли провериться только после того, как мы уже набрали необходимое
	// количество подходящих, поэтому нужно обработать и их, т.е. последовательность будет избыточна
	ch.check(amount - 1, input, check_func, callback, false);
	
	// сортируем многочлены в лексико-графическом порядке для получения правильной последовательности
	std::sort(irred_polys_table.begin(), irred_polys_table.end(), [](const Polynom &a, const Polynom &b) {
		if (a.degree() == b.degree())
		{
			for (BasicInt i = static_cast<BasicInt>(a.degree()); i > 0; --i)
			{
				if (a[i] != b[i]) { return a[i] < b[i]; }
			}
			return a[0] < b[0];
		} else { return a.degree() < b.degree(); }
	});
	// выкидываем лишние с конца, чтобы осталось только требуемое число
	while (irred_polys_table.size() > amount)
	{
		c_defect -= irred_polys_table.back().degree() - 1;
		irred_polys_table.pop_back();
	}
	
	if ( c_defect >= NBITS )
	{
		throw std::logic_error("Spatial dimensionality is too big.");
	}
	
	return irred_polys_table;
}

/*!
 *
 *	Let \f$\{p_n(x)\}_{n=0}^{\infty}\f$ be a sequence of irreducible polynomials over GF(2) sorted by increasing of degrees.
 *	To generate a base 2 \f$(t,m,s)\f$-net using \f$\{p_n(x)\}\f$, we need the sum \f$\sum_{k=0}^{s-1}(deg(p_k) - 1) = t\f$
 *	to be smaller than \f$\max\{m\} =\f$ #NBITS.\n
 *	Let's say, that \f$p(x) = \sum_{i=0}^{deg(p)}p_i\cdot x^i,\ \ p_i \in\f$ GF(2) is equal to the integer
 *	\f$\overline{p} = \sum_{i=0}^{\log_2 \overline{p}}p_i\cdot 2^i\f$.
 *	To generate \f$\{p_n(x)\}_{n=0}^{\textbf{s_parameter} - 1}\f$ we'll test every polynomial, equal to integers
 *	\f$\overline{p}\f$ running increasingly from 3 with step 2, with Berlekamp's irreducibility test, supposing that
 *	\f$p_0(x) = x\f$.
 *
 *	\return
 *		Vector of \b s_parameter least-degree irreducible polynomials.
 *
 *	\author
 *		Vadim Piven\n
 *		Alexey Burimov
 *
 */
template <typename T, unsigned int NBITS>
std::vector< typename Nied2Generator<T,NBITS>::Polynom > Nied2Generator<T,NBITS>::generate_irreducible(
																			//! [in] amount of least-degree irreducible to generate
																			std::size_t const s_parameter)
{
	CountInt coeffs_number = 3;
	c_defect = 0;
	
	if ( s_parameter == 0 )
	{ throw std::logic_error("Spatial dimensionality can't be equal to 0."); }
	
	std::vector<Polynom> irred_polys_table;
	irred_polys_table.reserve(s_parameter);
	irred_polys_table.emplace_back(Polynom({0, 1}));
	
	auto number_to_poly = [](CountInt coeffs_number) {
		Polynom poly;
		BasicInt degree = 0;
		while ( coeffs_number >> degree != 0 ) { ++degree; }
		poly.data().reserve(degree--);
		for (BasicInt i = 0; i <= degree; ++i)
		{
			poly.data().emplace_back((coeffs_number >> i) & 1);
		}
		return poly;
	};
	
	while ( irred_polys_table.size() < s_parameter && c_defect < NBITS )
	{
		irred_polys_table.emplace_back(number_to_poly(coeffs_number));
		
		while ( !is_irreducible_berlekamp(irred_polys_table.back()) )
		{
			coeffs_number += 2;
			irred_polys_table.back() = number_to_poly(coeffs_number);
		}
		coeffs_number += 2;
		c_defect += irred_polys_table.back().degree() - 1;
	}
	
	if ( c_defect >= NBITS )
	{
		//Then we can't generate a (t_parameter,m,s_parameter)-net, because t_parameter will be >= max{m} = NBITS.
		throw std::logic_error("Spatial dimensionality is too big.");
	}
	
	return irred_polys_table;
}

/*!
 *	\todo Implement error handling
 *
 *	Let \f$[p_n(x)]_{n=0}^{s-1}\f$ be a vector of irreducible polynomials over GF(2).
 *	To generate a base 2 \f$(t,m,s)\f$-net using \f$\{p_n(x)\}\f$, we need the sum \f$\sum_{n=0}^{s-1}(deg(p_n) - 1) = t\f$
 *	to be smaller than \f$\max\{m\} =\f$ #NBITS.\n
 *	Let \f$[d_n]\f$ be equal to \b degrees and let's say, that \f$p(x) = \sum_{i=0}^{deg(p)}p_i\cdot x^i,\ \ p_i \in\f$ GF(2) is
 *	equal to the integer \f$\overline{p} = \sum_{i=0}^{\log_2 \overline{p}}p_i\cdot 2^i\f$.
 *	To generate \f$[p_n(x)]\f$ for each \f$n\f$ we will test polynomials, equal to every integer
 *	\f$\overline{p_n}\f$ such that \lfloor \log_2 \overline{p_n}\rfloor = d_n\f$ with Berlekamp's irreducibility test.
 *
 *	\return
 *		Vector of \b s_parameter least-degree irreducible polynomials.
 *
 *	\author
 *		Vadim Piven\n
 *		Alexey Burimov
 *
 */
template <typename T, unsigned int NBITS>
std::vector< typename Nied2Generator<T,NBITS>::Polynom > Nied2Generator<T,NBITS>::generate_irreducible_with_degrees(
																						 //! [in] amount of least-degree irreducible to generate
																						 std::vector<Nied2Generator<T,NBITS>::BasicInt> const &degrees)
{
	// counter of possible t values for (t,m,s)-nets with sush irred. polynomials degrees.
	c_defect = 0;
	// coefficient number of polynomial over GF(2) such that n-th bit is equal to n-th coefficient of polynomial.
	
	auto number_to_poly = [](CountInt coeffs_number) {
		Polynom poly;
		BasicInt degree = 0;
		while ( coeffs_number >> degree != 0 ) { ++degree; }
		poly.data().reserve(degree--);
		for (BasicInt i = 0; i <= degree; ++i)
		{
			poly.data().emplace_back((coeffs_number >> i) & 1);
		}
		return poly;
	};
	
	// map of elements { key, value }:
	//		key == degree;
	//		value == {degree's appearances count in degrees, current coefficient number } },
	//	where degree is taken from degrees vector.
	std::map< BasicInt, std::pair<BasicInt, BasicInt> > degrees_counter;
	
	BasicInt i = 0;
	for (i = 0; i < c_dim && c_defect < NBITS && degrees[i] != 0; ++i)
	{
		// if there's no key == degree, it will be added with value 0 and then will be incremented to 1;
		//		else the value associated with key == degree will be just incremented.
		degrees_counter[degrees[i]].first += 1;
		degrees_counter[degrees[i]].second = (1 << degrees[i]) + (degrees[i] != 1);
		c_defect += degrees[i] - 1;
	}
	
	if ( c_defect >= NBITS || (i < c_dim && degrees[i] == 0) )
	{
		throw std::logic_error("There can't be generated any (t, m, degrees.size())-nets with such polynomial's degrees for any m <= " + std::to_string(NBITS));
	}
	
	std::vector<Polynom> irred_polys_table;
	irred_polys_table.reserve(c_dim);
	
	for (i = 0; i < c_dim && degrees_counter[degrees[i]].second < (2 << degrees[i]); ++i)
	{
		irred_polys_table.emplace_back( number_to_poly(degrees_counter[degrees[i]].second) );
		
		while ( degrees_counter[degrees[i]].second < (2 << degrees[i]) && !is_irreducible_berlekamp(irred_polys_table.back()) )
		{
			degrees_counter[degrees[i]].second += 2;
			irred_polys_table.back() = number_to_poly(degrees_counter[degrees[i]].second);
		}
		degrees_counter[degrees[i]].second += (degrees[i] != 1) + 1;
	}
	
	if ( degrees_counter[degrees[i]].second >= (2 << degrees[i]) )
	{
		throw std::logic_error("There doesn't exist " + std::to_string(degrees_counter[degrees[i]].first) + " irreducible polynomials of degree " + std::to_string(degrees[i]) + " over GF(2)");
	}
	
	return irred_polys_table;
}

/*!
 *	Function performs multiplication of two polynomials over GF(2).\n
 *	Polynomials are stored as arrays of coefficients and have
 *	the coefficient of degree N as the N-th element of an array.\n
 *
 *	\par Modified
 *		29 March 2003
 *
 *	\author
 *		Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.\n
 *		C++ version by John Burkardt
 *
 */
template <typename T, unsigned int NBITS>
typename Nied2Generator<T,NBITS>::Polynom Nied2Generator<T,NBITS>::multiply_poly2(
													   //! [in] the first polynomial
													   Nied2Generator<T,NBITS>::Polynom const &poly_pa,
													   //! [in] the second polynomial
													   Nied2Generator<T,NBITS>::Polynom const &poly_pb)
{
	std::size_t jlo, jhi;
	BasicInt term;
	Polynom poly_pc(std::vector< irrpoly::gf<2> >(poly_pa.degree() + poly_pb.size(), 1));
	
	for (std::size_t i = 0; i <= poly_pa.degree() + poly_pb.degree(); ++i)
	{
		jlo = ( i < poly_pa.degree() ) ? 0 : static_cast<BasicInt>(i - poly_pa.degree());
		jhi = ( i < poly_pb.degree() ) ? i : static_cast<BasicInt>(poly_pb.degree());
		//
		// Find pc_i as a mod2 sum over all pa_j*pb_k such that j + k = i
		//
		term = 0;
		for (std::size_t j = jlo; j <= jhi; ++j)
		{
			term ^= poly_pa[i - j].data() & poly_pb[j].data();
		}
		poly_pc[i] = term;
	}
	
	return poly_pc;
}


template <typename T, unsigned int NBITS>
Nied2Generator<T,NBITS>::Nied2Generator(
							   //! [in] spatial dimensionality
							   Nied2Generator<T,NBITS>::BasicInt const dim,
							   //! [in] flag, defining usage of a single-thread or a multy-thead #c_irred_polys generation
							   bool const in_parallel) :
	c_dim(dim),
	c_irred_polys( (this->*(in_parallel ? &Nied2Generator<T,NBITS>::generate_irreducible_in_parallel : &Nied2Generator<T,NBITS>::generate_irreducible))(dim) ),
	c_cj(std::vector<IntPoint>(c_dim, IntPoint(NBITS, 0)))
{
	initialize_c();
}

/*! Let \f$\{d_i\}_{i=0}^{s-1}\f$ is equal to \b degrees_of_irred, where \f$s = \f$\b degrees_of_irred.size().\n
 *	Then \f$t = \sum_{i=0}^{s-1} (d_i - 1)\f$.
 *
 *	\throws std::logic_error if \f$t\f$ becomes larger than (#NBITS-1) or if ???.
 */
template <typename T, unsigned int NBITS>
Nied2Generator<T,NBITS>::Nied2Generator(
							   //! [in] vector of degrees of irreducible polynomials, which size must be equal to spatial dimensionality.
							   std::vector<BasicInt> const &degrees_of_irred) :
	c_dim(static_cast<BasicInt>(degrees_of_irred.size())),
	c_irred_polys(generate_irreducible_with_degrees(degrees_of_irred)),
	c_cj(std::vector<IntPoint>(c_dim, IntPoint(NBITS, 0)))
{
	initialize_c();
}

template <typename T, unsigned int NBITS>
typename Nied2Generator<T, NBITS>::BasicInt Nied2Generator<T, NBITS>::get_s(void) const
{
	return c_dim;
}

template <typename T, unsigned int NBITS>
typename Nied2Generator<T, NBITS>::BasicInt Nied2Generator<T, NBITS>::get_t(void) const
{
	return c_defect;
}

template <typename T, unsigned int NBITS>
typename Nied2Generator<T, NBITS>::BasicInt Nied2Generator<T, NBITS>::get_nbits(void) const
{
	return NBITS;
}

/*!
 *	Calculation of the constants \f$\{v_n\}\f$ implemented as
 *	described in the reference (BFN) section 3.3.\n
 *
 *	\par Modified
 *		29 March 2003
 *
 *	\par Author
 *		Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.\n
 *    	C++ version by John Burkardt.
 *
 *	\par Reference
 *		Paul Bratley, Bennett Fox, Harald Niederreiter,
 *		Algorithm 738:
 *		Programs to Generate Niederreiter's Low-Discrepancy Sequences,
 *		ACM Transactions on Mathematical Software,
 *		Volume 20, Number 4, pages 494-495, 1994.
 *
 */
template <typename T, unsigned int NBITS>
void Nied2Generator<T,NBITS>::calculate_v(
								 //! [in] the appropriate irreducible polynomial \f$p_i(x)\f$ for the dimension currently being considered
								 Polynom const &poly_pi,
								 /*! [in, out] on input, \f$b(x)\f$ is the polynomial
								  defined in section 2.3 of BFN. The degree \f$deg(b(x))\f$ implicitly defines
								  the parameter j of section 3.3, by \f$deg(b(x)) = e \cdot (j-1)\f$.  On output,
								  \f$b(x)\f$ has been multiplied by \f$p_i(x)\f$, so its degree is now \f$e \cdot j\f$*/
								 Polynom &poly_b,
								 //! [out] the computed \f$\{v_n\}\f$ array
								 BasicInt v[])
{
	Polynom poly_h;
	BasicInt term;
	BasicInt const v_size = static_cast<BasicInt>(NBITS + c_irred_polys.back().degree() + 1);
	//
	//  The polynomial h(x) = p_i(x)**(j-1) = b(x) on arrival.
	//
	//  In section 3.3, the values of h_i are defined with a minus sign,
	//  but in GF(2): h_i = -h_i for any value.
	//
	poly_h = poly_b;
	//
	//  Now choose a value of K_j as defined in section 3.3.
	//  We must have 0 <= K_j < e*j = m.
	//  The limit condition on K_j does not seem very relevant
	//  in this program.
	//	Let K_j = e*(j - 1) = deg(h).
	//
	//  Multiply b(x) by p_i(x) so b(x) becomes p(x)**j.
	//  In section 2.3, the values of b_i are defined with a minus sign,
	//  but in GF(2): b_i = -b_i for any value.
	//
	poly_b = multiply_poly2(poly_pi, poly_b);
	//
	//  Choose values of {v_n} in accordance with the conditions in section 3.3.
	//
	for (BasicInt r = 0; r < poly_h.degree(); ++r)
	{
		v[r] = 0;
	}
	v[poly_h.degree()] = 1;
	
	for (BasicInt r = static_cast<BasicInt>(poly_h.degree() + 1); r <= poly_b.degree() - 1; ++r)
	{
		v[r] = 1;
	}
	//
	//  Calculate the remaining v_n's using the recursion of section 2.3,
	//  remembering that the b_i's have the opposite sign.
	//
	for (BasicInt r = 0; r < v_size - poly_b.degree(); ++r)
	{
		term = 0;
		for (BasicInt i = 0; i <= poly_b.degree() - 1; ++i)
		{
			term ^= poly_b[i].data() & v[r + i];
		}
		v[r + poly_b.degree()] = term;
	}
	
	return;
}

/*!
 *	This program calculates the values of the constants \f$c^{(i)}_{jr}\f$ denoted as c(i, j, r).
 *	As far as possible, Niederreiter's notation is used.\n
 *	For each value of i, we first calculate all the corresponding
 *	values of c (all these values are either 0 or 1) and pack them into the 2D-array \b cj thats
 *	denoting \f$С^{(i)}_r\f$ \n
 *	in such a way that \b cj[i][r] holds the values of c(i, j, r)
 *	for the indicated i and r and for every j from 0 to (\b NBITS - 1).\n
 *	The most significant bit of \b cj[i][r] (not counting the sign bit) is c(i, 0, r) and the least
 *	significant bit is c(i, \b NBITS - 1, r) that is equivalent to\n
 *	\f$С^{(i)}_r = \sum\limits_{j=0}^{\textbf{NBITS}-1} c^{(i)}_{jr}\cdot 2^{\textbf{NBITS}-1-j}\f$.
 *
 *	\par Modified
 *		29 March 2003
 *
 *	\author
 *		Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.\n
 *		C++ version by John Burkardt.
 *
 *	\par Reference
 *		R Lidl, Harald Niederreiter,
 *		Finite Fields,
 *		Cambridge University Press, 1984, page 553.\n
 *		Harald Niederreiter,
 *		Low-discrepancy and low-dispersion sequences,
 *		Journal of Number Theory,
 *		Volume 30, 1988, pages 51-70.
 *
 */
template <typename T, unsigned int NBITS>
void Nied2Generator<T,NBITS>::initialize_c(void)
{
	Polynom const poly_ident{1};
	Polynom poly_b;
	Polynom poly_pi;
	BasicInt v[NBITS + c_irred_polys.back().degree() + 1];
	
	for (BasicInt i = 0, u = 0; i < c_dim; ++i)
	{
		//
		//  For each dimension, we need to calculate powers of an
		//  appropriate irreducible polynomial:  see Niederreiter
		//  page 65, just below equation (19).
		//
		//  Copy the appropriate irreducible polynomial p_i(x) into pi,
		//  and its degree into e.  Set polynomial b(x) = 1 and copy it into b.
		//  m is the degree of b(x).  Subsequently b(x) will hold higher
		//  powers of p_i(x).
		//
		poly_pi = c_irred_polys[i];
		poly_b = poly_ident;
		//
		//  Niederreiter (page 56, after equation (7), defines two
		//  variables q(i, j) and u(i, j).  We do not need q explicitly, but we do need u.
		//
		u = 0;
		
		for (BasicInt j = 0; j < NBITS; ++j)
		{
			//
			//  If u = 0, we need to set B to the next power of p_i(x)
			//  and recalculate {v_n}. This is done by subroutine CALCV2.
			//
			if ( u == 0 )
			{
				calculate_v(poly_pi, poly_b, v);
			}
			//
			//  Now c is obtained from v.
			//	!!! We can think about c(i0,j0,r) as a (NBITS - j0 - 1)-th bit in binary
			//		representation of a number cj(i0, r).
			//
			for (BasicInt r = 0; r < NBITS; ++r)
			{
				c_cj[i][r] |= ((NextInt)v[r + u]) << (NBITS - 1 - j); //???
			}
			//
			//  Increment u. If u = e, then u = 0 and in Niederreiter's
			//  paper q = q + 1.  Here, however, q is not used explicitly.
			//
			u = u + 1;
			if ( u == c_irred_polys[i].degree() )
			{
				u = 0;
			}
		}
	}
	
	return;
}

/*!
 *
 *
 */
template <typename T, unsigned int NBITS>
void Nied2Generator<T,NBITS>::load_point_int(
									//! [in, out] #c_dim-sized vector, which will contain the generated \f$Q_\textbf{pos}\f$
									Nied2Generator<T,NBITS>::IntPoint &point,
									//! [in] sequence number of generated point
									Nied2Generator<T,NBITS>::CountInt const pos)
{
	CountInt pos_gray_code;
	BasicInt r;
	
	pos_gray_code = pos ^ (pos >> 1);
	
	for (BasicInt i = 0; i < c_dim; ++i)
	{
		point[i] = 0;
	}
	
	r = 0;
	while ( pos_gray_code != 0 && r < NBITS )
	{
		if ( (pos_gray_code & 1) != 0 )
		{
			for (BasicInt i = 0; i < c_dim; ++i)
			{
				point[i] ^= c_cj[i][r];
			}
		}
		pos_gray_code >>= 1;
		++r;
	}
	
}

/*!
 *	Generation based on Gray code representation of \b pos.
 *
 *  \par Reference
 *		Harald Niederreiter,
 *		Low-discrepancy and low-dispersion sequences,
 *		Journal of Number Theory,
 *		Volume 30, 1988, pages 51-70.
 *
 */
template <typename T, unsigned int NBITS>
typename Nied2Generator<T,NBITS>::IntPoint Nied2Generator<T,NBITS>::get_point_int(
													   //! [in] sequence number of generated point
													   Nied2Generator<T,NBITS>::CountInt const pos)
{
	IntPoint point_Qn(c_dim);
	load_point_int(point_Qn, pos);
	
	return point_Qn;
}

/*!
 *
 *
 */
template <typename T, unsigned int NBITS>
void Nied2Generator<T,NBITS>::load_point_real(
									 //! [in, out] #c_dim-sized #RealPoint, which will contain the generated \f$x_\textbf{pos}\f$
									 Nied2Generator<T,NBITS>::RealPoint &point,
									 //! [in] sequence number of generated point
									 Nied2Generator<T,NBITS>::CountInt const pos)
{
	IntPoint  point_Qn(c_dim);
	load_point_int(point_Qn, pos);
	
	for (BasicInt i = 0; i < c_dim; ++i)
	{
		point[i] = static_cast<Real>(point_Qn[i])*c_recip;
	}
}

/*!
 *
 *
 */
template <typename T, unsigned int NBITS>
typename Nied2Generator<T,NBITS>::RealPoint Nied2Generator<T,NBITS>::get_point_real(
														 //! [in] sequence number of generated point
														 Nied2Generator<T,NBITS>::CountInt const pos)
{
	RealPoint point_xn(c_dim);
	load_point_real(point_xn, pos);
	
	return point_xn;
}

/*!
 *	Calling this function with \b pos = 0 causes bad results.
 */
template <typename T, unsigned int NBITS>
void Nied2Generator<T,NBITS>::load_next_int(
								   //! [in, out] #c_dim-sized #IntPoint, which will contain the generated \f$Q_\textbf{pos}\f$
								   Nied2Generator<T,NBITS>::IntPoint &point,
								   //! [in] sequence number of generated point
								   Nied2Generator<T,NBITS>::CountInt pos,
								   //! [in] previous point \f$Q_{\textbf{pos} - 1}\f$
								   Nied2Generator<T,NBITS>::IntPoint const &prev_point)
{
	BasicInt rightmost_zero_bit_pos = 0;
	pos -= (pos != 0);
	while ( (pos & 1) != 0 && rightmost_zero_bit_pos < NBITS - 1 )
	{
		++rightmost_zero_bit_pos;
		pos >>= 1;
	}
	//
	//  Compute the new numerators in vector Q.
	//
	for (BasicInt i = 0; i < c_dim; ++i)
	{
		point[i] = prev_point[i] ^ c_cj[i][rightmost_zero_bit_pos];
	}
}

/*!
 *
 *
 */
template <typename T, unsigned int NBITS>
void Nied2Generator<T,NBITS>::load_points_int(
									 //! [in, out] non-empty vector of #c_dim-sized #IntPoint, which will contain generated points
									 std::vector<Nied2Generator<T,NBITS>::IntPoint> &points,
									 //! [in] beginning sequence number of generated points
									 Nied2Generator<T,NBITS>::CountInt const pos)
{
	load_point_int(points[0], pos);
	
	for (CountInt seq_num = 1; seq_num < points.size(); ++seq_num)
	{
		load_next_int(points[seq_num], pos + seq_num, points[seq_num - 1]);
	}
}

/*!
 *
 *
 */
template <typename T, unsigned int NBITS>
void Nied2Generator<T,NBITS>::load_points_real(
									  //! [in, out] non-empty vector of #c_dim-sized #RealPoint, which will contain generated points
									  std::vector<Nied2Generator<T,NBITS>::RealPoint> &points,
									  //! [in] beginning sequence number of generated points
									  Nied2Generator<T,NBITS>::CountInt const pos)
{
	IntPoint point_Qn(c_dim);
	load_point_int(point_Qn, pos);
	
	for (BasicInt i = 0; i < c_dim; ++i)
	{
		points[0][i] = static_cast<Real>(point_Qn[i])*c_recip;
	}
	
	for (CountInt seq_num = 1; seq_num < points.size(); ++seq_num)
	{
		load_next_int(point_Qn, pos + seq_num, point_Qn);
		for (BasicInt i = 0; i < c_dim; ++i)
		{
			points[seq_num][i] = static_cast<Real>(point_Qn[i])*c_recip;
		}
	}
}

/*!
 *
 *
 */
template <typename T, unsigned int NBITS>
std::vector<typename Nied2Generator<T,NBITS>::IntPoint>  Nied2Generator<T,NBITS>::get_points_int(
																	  //! [in] beginning sequence number of generated points
																	  Nied2Generator<T,NBITS>::CountInt const pos,
																	  //! [in] amount of points being generated
																	  Nied2Generator<T,NBITS>::CountInt const amount)
{
	std::vector<IntPoint> points_Qn(amount, IntPoint(c_dim));
	load_points_int(points_Qn, pos);
	
	return points_Qn;
}

/*!
 *
 *
 */
template <typename T, unsigned int NBITS>
std::vector<typename Nied2Generator<T,NBITS>::RealPoint> Nied2Generator<T,NBITS>::get_points_real(
																	   //! [in] beginning sequence number of generated points
																	   Nied2Generator<T,NBITS>::CountInt const pos,
																	   //! [in] amount of points being generated
																	   Nied2Generator<T,NBITS>::CountInt const amount)
{
	std::vector<RealPoint> points_xn(amount, RealPoint(c_dim));
	load_points_real(points_xn, pos);
	
	return points_xn;
}



#endif // #ifndef NIEDERREITER2_HPP
