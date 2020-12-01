/**
 *	@file niederreiter2.hpp
 *
 *	@brief Includes the generator of a base 2 digital Niederreiter's (t,m,s)-nets.
 *
 *	@author Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter
 *	@author Previous C++ version by John Burkardt 2003-2009 [https://people.sc.fsu.edu/~jburkardt/cpp_src/niederreiter2/niederreiter2.html]
 * 	@author	__ 2019
 *
 */

#ifndef NIEDERREITER2_HPP
#define NIEDERREITER2_HPP

#include "./details/gf2poly.hpp"

#include <vector>
#include <cmath>		//for pow function
#include <type_traits>	//for type conditions in static_assert
#include <algorithm>	//for std::max_element function (is already included in "irrpoly/gfpoly.hpp")
#include <functional>	//for unified for_each_point* methods (is already included in "irrpoly/gfpoly.hpp")
#include <map>          //for usage in one constructor (is already included in "irrpoly/gf2poly.hpp")


/** @namespace tms
 *  @brief All the entities in the library are defined in this namespace. */
namespace tms
{
	/// Type used for storing integer values that are less than word size.
	using 	BasicInt 	= unsigned int;
	/// Type used for storing quantity-values and indexation-values, i.e. the total amount of points to generate.
	using	CountInt	= uintmax_t;
	/// Type used for storing coordinates of points of (t,m,s)-net.
	using 	Real		= long double;
	/// Type used for storing points of (t,m,s)-net.
	using 	Point		= std::vector<Real>;
	/// Type used for storing polynomials.
	using 	Polynomial 	= irrpoly::gfpoly;

	
	/** Provides an interface to create various generators of a digital base 2 Niederreiter's (t,m,s)-nets.	*/
	template <typename UIntType>
	class Niederreiter
	{
	static_assert(std::is_unsigned<UIntType>::value, "UIntType template parameter must be an unsigned integer type.");
		
		
	public:
		
		using IntPoint = std::vector<UIntType>;
		
		// Prevent generation of implicit-defined default constructors by compiler:
		Niederreiter(void)                 = delete;
		Niederreiter(Niederreiter const &) = delete;
		Niederreiter(Niederreiter &&)      = delete;
		// Prevent generation of implicit-defined copy assignment by compiler:
		Niederreiter& operator =(Niederreiter const &) = delete;
		
		/** Constructs the generator of (t,m,s)-net with specified m, s, and with induced least possible t.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] dim - s parameter of the net
		 *  @param [in] in_parallel - flag, defining the method that will be used to generate irreducible polynomials */
		Niederreiter(BasicInt const nbits,
					 BasicInt const dim,
					 bool     const in_parallel = false);
		/** Constructs the generator of (t,m,s)-net with specified m, s degrees of initial irreducible polynomials,
		 *  and (optional, for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] degrees_of_irrpolys - vector of the s degrees of initial irreducible polynomials
		 *  @param [in] matrix_of_initial_values - matrix of initial values for all recursive sequences */
		Niederreiter(BasicInt                              const  nbits,
					 std::vector<BasicInt>                 const &degrees_of_irrpolys,
					 std::vector< std::vector<uintmax_t> > const &matrix_of_initial_values = {});
		/** Constructs the generator of (t,m,s)-net with specified m, s degrees of initial irreducible polynomials,
		 *  and (optional, for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] degrees_of_irrpolys - vector of the s degrees of initial irreducible polynomials
		 *  @param [in] matrix_of_initial_values - matrix of initial values for all recursive sequences */
		Niederreiter(BasicInt                                        const  nbits,
					 std::initializer_list<BasicInt>                 const &degrees_of_irrpolys,
					 std::initializer_list< std::vector<uintmax_t> > const &matrix_of_initial_values = {});
		/** Constructs the generator of (t,m,s)-net with specified m, s initial irreducible polynomials, and (optional,
		 *  for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] irrpolys_coeffs - initializer list of the s initial irreducible polynomials coefficients
		 *  @param [in] matrix_of_initial_values - initializer list of vectors of initial values for all recursive sequences */
		Niederreiter(BasicInt                              const  nbits,
					 std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs,
					 std::vector< std::vector<uintmax_t> > const &matrix_of_initial_values = {});
		/** Constructs the generator of (t,m,s)-net with specified m, s initial irreducible polynomials, and (optional,
		 *  for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] irrpolys_coeffs - initializer list of the s initial irreducible polynomials coefficients
		 *  @param [in] matrix_of_initial_values - initializer list of vectors of initial values for all recursive sequences */
		Niederreiter(BasicInt                                        const  nbits,
					 std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs,
					 std::initializer_list< std::vector<uintmax_t> > const &matrix_of_initial_values = {});
		
		~Niederreiter(void) = default;
		
		Niederreiter& operator =(Niederreiter &&tmp_generator) = default;
		
		/** Returns t parameter of (t,m,s)-net. */
		BasicInt    get_t(void) const;
		/** Returns m parameter of (t,m,s)-net. */
		BasicInt    get_m(void) const;
		/** Returns s parameter of (t,m,s)-net. */
		BasicInt    get_s(void) const;
		
		
		/** Returns vector of direction numbers corresponging to certain dimension.
		 *  @param [in] dim - dimension */
		std::vector<UIntType> get_direction_numbers(BasicInt const dim) const;
		/** Returns generating matrix corresponging to certain dimension.
		 *  @param [in] dim - dimension */
		std::vector< std::vector<BasicInt> > get_gamma_matrix(BasicInt const dim) const;
		
		/** Generates scaled point of (t,m,s)-net with certain number, enumerated according to Gray's code.
		 *  @param [in] pos - sequence number of scaled generated point */
		IntPoint    generate_point_int(CountInt const pos) const;
		/** Generates point of (t,m,s)-net with certain number, enumerated according to Gray's code.
		 *  @param [in] pos - sequence number of generated point */
		Point       generate_point    (CountInt const pos) const;
		
		/** Generates point of (t,m,s)-net with certain number.
		 *  @param [in] pos - sequence number of generated point */
		Point       generate_point_classical(CountInt const pos) const;
		
		/** Sequentially generates a section of scaled (t,m,s)-net points and applies the handler function to each pair:
		 *  (point, point's number).
		 *  @param [in] handler - handler function to apply
		 *  @param [in] amount - amount of points in the section of the net
		 *  @param [in] pos - number of the first point in the section of the net */
		void        for_each_point_int(std::function<void (IntPoint const &, CountInt)> handler,
									   CountInt                                         amount,
									   CountInt                                         pos = 0) const;
		/** Sequentially generates a section of (t,m,s)-net points and applies the handler function to each pair:
		 *  (point, point's number).
		 *  @param [in] handler - handler function to apply
		 *  @param [in] amount - amount of points in the section of the net
		 *  @param [in] pos - number of the first point in the section of the net */
		void        for_each_point    (std::function<void (Point const &, CountInt)> handler,
									   CountInt                                      amount,
									   CountInt                                      pos = 0) const;
		

	private:
		
		/// t parameter of net (aka quality parameter).
		BasicInt                             m_quality_param;
		/// m parameter of net that defines net cardinality and bitwidth of computations.
		BasicInt                             m_nbits;
		/// s parameter of net that defines spatial dimensionality.
		BasicInt                             m_dim;
		/// Coefficient, equal to \f$ 2^{-m} \f$.
		Real                                 m_recip;
		/// Vector of s initial irreducible polynomials of a (t,m,s)-net.
		std::vector<Polynomial>	             m_irrpolys;
		/// Matrix of a (t,m,s)-net direction numbers g[i](k).
		std::vector< std::vector<UIntType> > m_direction_numbers;
		
		/** Verifies the initial values of recursive sequences supplied by user during generator construction time right
		 *  after irredusible polynomials were successfully verified.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] valid_irrpolys - vector of successfully verified initial irreducible polynomials
		 *  @param [in] matrix_of_initial_values - matrix of initial values to be verified */
		static bool is_matrix_of_initial_values_valid(BasicInt                              const  nbits,
													  std::vector<Polynomial>               const &valid_irrpolys,
													  std::vector< std::vector<uintmax_t> > const &matrix_of_initial_values);
		
		/** Fills vector sequentially with elements of the recursive sequence over GF(2), defined by certain characteristic
		 *  polynomial and initial values.
		 *  @param [out] container - vector to be filled
		 *  @param [in] initial_values - initial values of recursive sequence packed into an integer
		 *  @param [in] char_poly - characteristic polynomial of the sequence */
		static void fill_container_recursively(std::vector<BasicInt> &container, uintmax_t const initial_values, Polynomial const &char_poly);
		/** Initializes (t,m,s)-net direction numbers g[i](k).
		 *  @param [in] matrix_of_initial_values - initial values provided by user (empty matrix otherwise) */
		       void initialize_direction_numbers(std::vector< std::vector<uintmax_t> > const &matrix_of_initial_values);
		
		/** Returns transformed given integer point into a real point from a unit hypercube.
		 *  @param [in] point_int - point to be transformed (multiplied by \f$ 2^{-m} \f$) */
		Point cast_point_int_to_real(IntPoint const &point_int) const;
		
		/** Stores into the integer vector scaled (t,m,s)-net point with certain number, enumerated according to Gray's code.
		 *  @param [out] point - storage vector
		 *  @param [in] pos - generated scaled (t,m,s)-net point number */
		void  store_point_int     (IntPoint &point, CountInt pos) const;
		/** Stores into the integer vector scaled (t,m,s)-net point with certain number, enumerated according to Gray's code,
		 *  computed using the previous point.
		 *  @param [out] point - storage vector
		 *  @param [in] pos - generated scaled (t,m,s)-net point number (should be greater then 0)
		 *  @param [in] prev_point - scaled (t,m,s)-net point with the previous point number */
		void  store_next_point_int(IntPoint &point, CountInt pos, IntPoint const &prev_point) const;
		
#ifdef TMS_EXPERIMENTAL
		/* Experimental features{ */
		
	private:
		
		std::vector< std::vector<uintmax_t> > m_matrix_of_initial_values;
		/** Resets specified section in a certain (t,m,s)-net generating matrix.
		 *  */
		void reset_gamma_matrix_section(BasicInt const dim, BasicInt const q, Polynomial const &char_poly, std::vector<BasicInt> &alpha);
		
		
	public:
		
		uintmax_t get_irrpoly_degree(BasicInt const dim) const;
		
		/// Generates scaled by 2**m point with sequence number 'pos' of naturally enumerated (t,m,s)-net.
		IntPoint	     generate_point_int_classical(CountInt const pos) const;
		
		/// Compute rank of  Gamma[dim] matrix of generator.
		BasicInt		 get_rank_of_gamma_matrix(BasicInt const dim) const;
		
		void set_gamma_initial_values_in_section(BasicInt const dim, BasicInt const sec_pos, uintmax_t const initial_values);
		
		void set_gamma_initial_values(BasicInt const dim, std::vector<uintmax_t> const &vector_of_initial_values);
		
		/// Reforms the generator of (t,m,s)-net to the generator of (t,m',s)-net, where m' = 'nbits'
		void decrease_nbits(BasicInt const nbits);
		
		/* }Experimental features */
#endif
	};
	
	
	template <typename UIntType>
	Niederreiter<UIntType>::Niederreiter(BasicInt const nbits,
										 BasicInt const dim,
										 bool     const in_parallel) :
	    m_quality_param(0),
		m_nbits(nbits),
	    m_dim(dim),
		m_recip( pow(2, -static_cast<Real>(nbits)) ),
		m_irrpolys( (in_parallel ? tms::gf2poly::generate_irrpolys_in_parallel : tms::gf2poly::generate_irrpolys)(dim, nbits) ),
		m_direction_numbers( std::vector< std::vector<UIntType> >(m_dim, std::vector<UIntType>(nbits, 0)) )
	{
		if ( nbits > sizeof(UIntType)*8 )
		{
			throw std::logic_error("\nnbits is too high");
		}
		
		// if m_dim == 0 or
		//  if m (m_nbits) is too small for desired s (m_dim)
		if ( m_dim == 0 || m_irrpolys.size() != m_dim )
		{
			throw std::logic_error("\nWrong net's parameters");
		}
		// computing t parameter of the net:
		for (Polynomial const &poly : m_irrpolys)
		{
			m_quality_param += poly.size();
		}
		m_quality_param -= 2*m_dim;
		
#ifdef TMS_EXPERIMENTAL
		m_matrix_of_initial_values = std::vector< std::vector<uintmax_t> >(m_dim);
#endif
		
		initialize_direction_numbers(std::vector< std::vector<uintmax_t> >());
	}
	
	
	template <typename UIntType>
	Niederreiter<UIntType>::Niederreiter(BasicInt                              const  nbits,
										 std::vector<BasicInt>                 const &degrees_of_irrpolys,
										 std::vector< std::vector<uintmax_t> > const &matrix_of_initial_values) :
		m_quality_param(0),
	    m_nbits(nbits),
		m_dim(static_cast<BasicInt>(degrees_of_irrpolys.size())),
	    m_recip( pow(2, -static_cast<Real>(nbits)) ),
		m_irrpolys(gf2poly::generate_irrpolys_with_degrees(degrees_of_irrpolys, nbits)),
		m_direction_numbers(std::vector< std::vector<UIntType> >(m_dim, std::vector<UIntType>(nbits, 0)))
	{
		if ( nbits > sizeof(UIntType)*8 )
		{
			throw std::logic_error("\nnbits is too high");
		}
		
		for (BasicInt const degree : degrees_of_irrpolys)
		{
			m_quality_param += degree;
		}
		m_quality_param -= degrees_of_irrpolys.size();
		// if degrees_of_irrpolys is empty or
		//  if there is no polynomials with such degrees that induced t (m_quality_param) <= m (nbits)
		if ( degrees_of_irrpolys.size() == 0 || m_irrpolys.size() != degrees_of_irrpolys.size() )
		{
			throw std::logic_error("\nWrong polynomial degrees or the nbits parameter");
		}
		// Verification of user-supplied initial values:
		if ( !is_matrix_of_initial_values_valid(m_nbits, m_irrpolys, matrix_of_initial_values) )
		{
			throw std::logic_error("\nNot all initial values provided");
		}
		
#ifdef TMS_EXPERIMENTAL
		m_matrix_of_initial_values = std::vector< std::vector<uintmax_t> >(m_dim);
#endif
		
		initialize_direction_numbers(matrix_of_initial_values);
	}


	template <typename UIntType>
	Niederreiter<UIntType>::Niederreiter(BasicInt                                        const  nbits,
										 std::initializer_list<BasicInt>                 const &degrees_of_irrpolys,
										 std::initializer_list< std::vector<uintmax_t> > const &matrix_of_initial_values) :
		Niederreiter(nbits,
					 std::vector<BasicInt>{degrees_of_irrpolys},
					 std::vector< std::vector<uintmax_t> >{matrix_of_initial_values})
	{}
	
	
	template <typename UIntType>
	Niederreiter<UIntType>::Niederreiter(BasicInt                              const  nbits,
										 std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs,
										 std::vector< std::vector<uintmax_t> > const &matrix_of_initial_values) :
	    m_quality_param(0),
	    m_nbits(nbits),
		m_dim( static_cast<BasicInt>(irrpolys_coeffs.size()) ),
		m_recip( pow(2, -static_cast<Real>(nbits)) ),
		m_direction_numbers(std::vector< std::vector<UIntType> >(m_dim, std::vector<UIntType>(nbits, 0)))
	{
		if ( nbits > sizeof(UIntType)*8 )
		{
			throw std::logic_error("\nnbits is too high");
		}
		
		//checking the polynomials:
		{
			m_irrpolys.reserve(m_dim);
			for (BasicInt i = 0; i < m_dim; ++i)
			{
				m_irrpolys.push_back( gf2poly::make_gf2poly(*(irrpolys_coeffs.begin() + i)) );
			}
			
			BasicInt i = 0;
			// In order to not include another header, uniqueness of polynomials checked with std::map container.
			// std::map was chosen for two reasons:
			// 1. it has been already included in gf2poly.hpp;
			// 2. we need to keep original order of polynomials.
			std::map<Polynomial, BasicInt, bool(*)(Polynomial const &, Polynomial const &)> just_set(irrpoly::operator!=);
			while ( i < m_irrpolys.size()  &&  just_set.size() == i && \
				    m_quality_param <= m_nbits    &&  irrpoly::is_irreducible_berlekamp(m_irrpolys[i]) )
			{
				just_set.insert(std::make_pair(m_irrpolys[i], i));
				m_quality_param += m_irrpolys[i].size() - 2;
				++i;
			}
			// if irrpolys.empty() or
			//  if irrpolys contains not unique polynomials
			//  if induced t (m_quality_param) is greater than m (m_nbits)
			if ( m_irrpolys.empty() || just_set.size() != m_irrpolys.size() || m_quality_param > m_nbits )
			{
				throw std::logic_error("\nWrong polynomials");
			}
		}
		
		//checking user-supplied initial values:
		if ( !is_matrix_of_initial_values_valid(m_nbits, m_irrpolys, matrix_of_initial_values) )
		{
			throw std::logic_error("\nNot all initial values provided");
		}
		
#ifdef TMS_EXPERIMENTAL
		m_matrix_of_initial_values = std::vector< std::vector<uintmax_t> >(m_dim);
#endif
		
		initialize_direction_numbers(matrix_of_initial_values);
	}
	
	
	template <typename UIntType>
	Niederreiter<UIntType>::Niederreiter(BasicInt                                        const  nbits,
										 std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs,
										 std::initializer_list< std::vector<uintmax_t> > const &matrix_of_initial_values) :
		Niederreiter(nbits,
					 std::vector< std::vector<uintmax_t> >{irrpolys_coeffs},
					 std::vector< std::vector<uintmax_t> >{matrix_of_initial_values})
	{}
	
	
	template <typename UIntType>
	BasicInt
	Niederreiter<UIntType>::get_s(void) const
	{
		return m_dim;
	}
	
	
	template <typename UIntType>
	BasicInt
	Niederreiter<UIntType>::get_t(void) const
	{
		return m_quality_param;
	}
	
	
	template <typename UIntType>
	BasicInt
	Niederreiter<UIntType>::get_m(void) const
	{
		return m_nbits;
	}
	
	
	template <typename UIntType>
	std::vector<UIntType>
	Niederreiter<UIntType>::get_direction_numbers(BasicInt const dim) const
	{
		return m_direction_numbers[dim];
	}
	
	
	template <typename UIntType>
	std::vector< std::vector<BasicInt> >
	Niederreiter<UIntType>::get_gamma_matrix(BasicInt const dim) const
	{
		std::vector< std::vector<BasicInt> > gamma_matrix(m_nbits);
		for (int j = 0; j < m_nbits; ++j)
		{
			gamma_matrix[j] = std::vector<BasicInt>(m_nbits);
			for (int k = 0; k < m_nbits; ++k)
			{
				gamma_matrix[j][k] = (m_direction_numbers[dim][k] >> (m_nbits - 1 - j)) & 1;
				//((UIntType)v[k + u]) << (m_nbits - 1 - j) << ' ';
			}
		}
		return gamma_matrix;
	}
	
	
	template <typename UIntType>
	typename Niederreiter<UIntType>::IntPoint
	Niederreiter<UIntType>::generate_point_int(CountInt const pos) const
	{
		IntPoint point_int(m_dim);
		store_point_int(point_int, pos);
		
		return point_int;
	}
	
	
	template <typename UIntType>
	Point
	Niederreiter<UIntType>::generate_point(CountInt const pos) const
	{
		IntPoint point_int(m_dim);
		store_point_int(point_int, pos);
		
		return cast_point_int_to_real(point_int);
	}
	
	
	template <typename UIntType>
	Point
	Niederreiter<UIntType>::generate_point_classical(CountInt const pos) const
	{
		Point point(m_dim, 0);
		for (int i = 0; i < m_dim; ++i)
		{
			UIntType acc = 0;
			for (int k = 0; k < m_nbits; ++k)
			{
				acc ^= m_direction_numbers[i][k] * ((pos >> k) & 1);
			}
			point[i] = static_cast<Real>(acc) * pow(2.0, -static_cast<Real>(m_nbits));
		}
		return point;
	}
	
	
	template <typename UIntType>
	void
	Niederreiter<UIntType>::for_each_point_int(std::function<void (IntPoint const &, CountInt)> handler,
											   CountInt amount,
											   CountInt pos) const
	{
		if ( amount != 0 )
		{
			IntPoint curr_int(m_dim);
			store_point_int(curr_int, pos);
			handler(curr_int, pos);
			while ( --amount )
			{
				++pos;
				store_next_point_int(curr_int, pos, curr_int);
				handler(curr_int, pos);
			}
		}
	}
	
	
	template <typename UIntType>
	void
	Niederreiter<UIntType>::for_each_point(std::function<void (Point const &, CountInt)> handler,
										   CountInt amount,
										   CountInt pos) const
	{
		if ( amount != 0 )
		{
			IntPoint curr_int(m_dim);
			store_point_int(curr_int, pos);
			handler(cast_point_int_to_real(curr_int), pos);
			while ( --amount )
			{
				++pos;
				store_next_point_int(curr_int, pos, curr_int);
				handler(cast_point_int_to_real(curr_int), pos);
			}
		}
	}
	
	
	
	//static
	template <typename UIntType>
	bool
	Niederreiter<UIntType>::is_matrix_of_initial_values_valid(BasicInt                const  nbits,
															  std::vector<Polynomial> const &valid_irrpolys,
															  std::vector< std::vector<uintmax_t> >       const &matrix_of_initial_values)
	{
		if ( matrix_of_initial_values.size() != 0 )
		{
			BasicInt const dim = static_cast<BasicInt>(valid_irrpolys.size());
			
			if ( matrix_of_initial_values.size() < dim )
			{
				return false;
			}
			BasicInt i = 0;
			
			while ( i < dim && \
				   matrix_of_initial_values[i].size() >= (nbits - 1)/valid_irrpolys[i].degree() + 1 )
			{
				++i;
			}
			// if number of gived initial values for any dimension is less than numnber of recursive sequences
			if ( i != dim )
			{
				return false;
			}
		}
		return true;
	}

	
	//static
	template <typename UIntType>
	void
	Niederreiter<UIntType>::fill_container_recursively(std::vector<BasicInt>       &container,
													   uintmax_t             const  initial_values,
													   Polynomial            const &char_poly)
	{
		CountInt const deg = static_cast<CountInt>(char_poly.degree());
		
		// Setting the initial values of the recursive sequence
		for (CountInt poly_i = 0; poly_i < deg && poly_i < container.size(); ++poly_i)
		{
			container[poly_i] = (initial_values >> (deg - 1 - poly_i)) & 1;
		}
		
		// Computing recurrently remaining elements, if needed
		for (CountInt seq_i = deg; seq_i < container.size(); ++seq_i)
		{
			BasicInt sum = 0;
			for (CountInt poly_i = 0; poly_i < deg; ++poly_i)
			{
				sum ^= char_poly[poly_i] & container[seq_i - deg + poly_i];
			}
			container[seq_i] = sum;
		}
	}
	
	
	template <typename UIntType>
	void
	Niederreiter<UIntType>::initialize_direction_numbers(std::vector< std::vector<uintmax_t> > const &matrix_of_initial_values)
	{
		std::vector<BasicInt> alpha(m_nbits + \
									std::max_element(
													 m_irrpolys.begin(), \
													 m_irrpolys.end(), \
													 [](Polynomial const &lpoly, Polynomial const &rpoly) { return lpoly.degree() < rpoly.degree(); } \
													 )->degree() - 1);
		
		for (BasicInt i = 0; i < m_dim; ++i)
		{
			BasicInt const e = static_cast<BasicInt>(m_irrpolys[i].degree());
			
			Polynomial poly_mu(tms::gf2poly::make_gf2poly({1}));
			
#ifdef TMS_EXPERIMENTAL
			m_matrix_of_initial_values[i].resize((m_nbits - 1)/e + 1);
#endif
			
			for (BasicInt j = 0; j < m_nbits; )
			{
				poly_mu = poly_mu * m_irrpolys[i];
				
				BasicInt const r_nbits   = m_nbits % e;
				BasicInt       rows_remaining_in_section = ( (j/e + 1)*e > m_nbits ) ? r_nbits : e;
				
				// If there wasn't any initial values provided by user, then default initial value setted.
				// Default initial values guarantee zero defect of (0,m,1)-subnet over i-th dimension.
				uintmax_t const initial_values = \
					( !matrix_of_initial_values.empty() ) ? \
						matrix_of_initial_values[i][j/e] \
					: \
						( r_nbits != 0 && j/e == (m_nbits - 1)/e ) ? \
							1 << (e - r_nbits) \
						: 	1;
#ifdef TMS_EXPERIMENTAL
				m_matrix_of_initial_values[i][j/e] = initial_values;
#endif
				alpha.resize(m_nbits - 1 + rows_remaining_in_section);
				fill_container_recursively(alpha, initial_values, poly_mu);
				
				// Here we interpret the j-th digit of direction number g[i](k) as gamma[i](j,k) - an
				// element of i-th generating matrix Gamma[i]
				while ( rows_remaining_in_section != 0 )
				{
					for (BasicInt k = 0, r = j % e; k < m_nbits; ++k)
					{
						m_direction_numbers[i][k] |= static_cast<UIntType>(alpha[k + r]) << (m_nbits - 1 - j);
					}
					++j;
					--rows_remaining_in_section;
				}
			}
		}
	}
	

	
	
	template <typename UIntType>
	Point
	Niederreiter<UIntType>::cast_point_int_to_real(Niederreiter<UIntType>::IntPoint const &point_int) const
	{
		Point point_real(m_dim);
		for (BasicInt i = 0; i < m_dim; ++i)
		{
			point_real[i] = static_cast<Real>(point_int[i])*m_recip;
		}
		return point_real;
	}
	

	template <typename UIntType>
	void
	Niederreiter<UIntType>::store_point_int(IntPoint &point,
											CountInt  pos) const
	{
		for (BasicInt i = 0; i < m_dim; ++i)
		{
			point[i] = 0;
		}
		
		CountInt pos_gray_code = (pos ^ (pos >> 1));
		for (BasicInt k = 0; pos_gray_code != 0 && k < m_nbits; ++k)
		{
			if ( pos_gray_code & 1 )
			{
				for (BasicInt i = 0; i < m_dim; ++i)
				{
					point[i] ^= m_direction_numbers[i][k];
				}
			}
			pos_gray_code >>= 1;
		}
	}
	

	template <typename UIntType>
	void
	Niederreiter<UIntType>::store_next_point_int(IntPoint       &point,
												 CountInt        pos,
												 IntPoint const &prev_point) const
	{
		// here we get pos = pow(2, rightmost zero bit position in pos), pos should be greater than 0
		pos = ~(pos - 1) & pos;
		
		// count trailing zeros of pos ( same that floor(log2(pos)) )
		BasicInt rightmost_zero_bit_pos = 0;
		while ( pos >>= 1 && rightmost_zero_bit_pos < (sizeof(CountInt)*8 - 1) )
		{
			++rightmost_zero_bit_pos;
		}

		for (BasicInt i = 0; i < m_dim; ++i)
		{
			point[i] = prev_point[i] ^ m_direction_numbers[i][rightmost_zero_bit_pos];
		}
	}

	
#ifdef TMS_EXPERIMENTAL
	//=========================================================================================================
	// 	Experimental features
	//=========================================================================================================
	
	template <typename UIntType>
	void
	Niederreiter<UIntType>::reset_gamma_matrix_section(BasicInt              const  dim,
													   BasicInt              const  q,
													   Polynomial            const &char_poly,
													   std::vector<BasicInt>       &alpha)
	{
		BasicInt const e = static_cast<BasicInt>(m_irrpolys[dim].degree());
		BasicInt       rows_remaining_in_section = ( (q + 1)*e > m_nbits ) ? m_nbits % e : e;
		
		alpha.resize(m_nbits - 1 + rows_remaining_in_section);
		fill_container_recursively(alpha, m_matrix_of_initial_values[dim][q], char_poly);
		
		UIntType const mask = ((((UIntType)1 << (m_nbits - e)) - 1) | ~(((UIntType)1 << m_nbits) - 1)) >> q*e;
		for (BasicInt k = 0; k < m_nbits; ++k)
		{
			m_direction_numbers[dim][k] &= mask;
			for (BasicInt r = 0; r < rows_remaining_in_section; ++r)
			{
				m_direction_numbers[dim][k] |= static_cast<UIntType>(alpha[k + r]) << (m_nbits - 1 - q*e - r);
			}
		}
	}
	
	
	template <typename UIntType>
	uintmax_t
	Niederreiter<UIntType>::get_irrpoly_degree(BasicInt const dim) const
	{
		return m_irrpolys[dim].degree();
	}
	
	
	template <typename UIntType>
	typename Niederreiter<UIntType>::IntPoint
	Niederreiter<UIntType>::generate_point_int_classical(CountInt const pos) const
	{
		IntPoint point(m_dim);
		for (int i = 0; i < m_dim; ++i)
		{
			point[i] = 0;
			for (int k = 0; k < m_nbits; ++k)
			{
				point[i] ^= m_direction_numbers[i][k] * ((pos >> k) & 1);
			}
		}
		return point;
	}
	
	
	template <typename UIntType>
	BasicInt
	Niederreiter<UIntType>::get_rank_of_gamma_matrix(BasicInt const dim) const
	{
		BasicInt rank = 0;
		int size = m_nbits;
		std::vector<UIntType> vecs = m_direction_numbers[dim];
		int digit = m_nbits - 1;
		int curr_row;
		while ( digit >= 0 && rank < size )
		{
			int i = curr_row = size - 1 - rank;
			while ( i >= 0 && (vecs[i] >> digit & 1) == 0 )
			{
				--i;
			}
			
			if ( i >= 0 )
			{
				std::swap(vecs[curr_row], vecs[i]);
				i = curr_row - 1;
				while ( i >= 0 )
				{
					if ( (vecs[i] >> digit & 1) == 1 )
					{
						vecs[i] ^= vecs[curr_row];
					}
					--i;
				}
				++rank;
			}
			--digit;
		}
		
		return rank;
	}
	
	
	template <typename UIntType>
	void
	Niederreiter<UIntType>::set_gamma_initial_values_in_section(BasicInt const dim, BasicInt const sec_pos, uintmax_t const initial_values)
	{
		Polynomial const &cur_irrpoly = m_irrpolys[dim];
		Polynomial char_poly = gf2poly::make_gf2poly({1});
		std::vector<BasicInt> alpha(m_nbits + cur_irrpoly.degree() - 1);
		
		m_matrix_of_initial_values[dim][sec_pos] = initial_values;
		for (BasicInt q = 0; q < sec_pos + 1; ++q)
		{
			char_poly = char_poly * cur_irrpoly;
		}
		reset_gamma_matrix_section(dim, sec_pos, char_poly, alpha);
	}
	
	
	template <typename UIntType>
	void
	Niederreiter<UIntType>::set_gamma_initial_values(BasicInt const dim, std::vector<uintmax_t> const &vector_of_initial_values)
	{
		if ( m_matrix_of_initial_values[dim].size() == vector_of_initial_values.size() )
		{
			Polynomial const &cur_irrpoly = m_irrpolys[dim];
			Polynomial char_poly = gf2poly::make_gf2poly({1});
			std::vector<BasicInt> alpha(m_nbits + cur_irrpoly.degree() - 1);
			
			for (BasicInt q = 0; q < vector_of_initial_values.size(); ++q)
			{
				char_poly = char_poly * cur_irrpoly;
				if ( m_matrix_of_initial_values[dim][q] != vector_of_initial_values[q] )
				{
					m_matrix_of_initial_values[dim][q] = vector_of_initial_values[q];
					reset_gamma_matrix_section(dim, q, char_poly, alpha);
				}
			}
		}
	}
	
	
	template <typename UIntType>
	void
	Niederreiter<UIntType>::decrease_nbits(const BasicInt nbits)
	{
		if ( nbits <= m_nbits && m_quality_param <= nbits )
		{
			for (BasicInt i = 0; i < m_dim; ++i)
			{
				m_direction_numbers[i].resize(nbits);
				for (BasicInt k = 0; k < nbits; ++k)
				{
					m_direction_numbers[i][k] >>= (m_nbits - nbits);
				}
				m_nbits = nbits;
			}
		}
	}
	

#endif // #ifdef TMS_EXPERIMENTAL
	
};// namespace tms


#endif // #ifndef NIEDERREITER2_HPP
