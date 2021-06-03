#ifndef TMS_NETS_DIGITAL_NET_HPP
#define TMS_NETS_DIGITAL_NET_HPP

#include "details/gf2poly.hpp"

#include <vector>
#include <cmath>		//for pow function
#include <algorithm>	//for std::max_element function (is already included in "irrpoly/gfpoly.hpp")
#include <functional>	//for unified for_each_point* methods (is already included in "irrpoly/gfpoly.hpp")


namespace tms
{
	/** Represents digital \f$(t, m, s)\f$-net over \f$\mathbb{F}_2\f$ */
	class DigitalNet
	{
	public:
		
		DigitalNet(DigitalNet const &) = default;
		DigitalNet(DigitalNet &&)      = default;
		DigitalNet& operator =(DigitalNet const &) = default;
		DigitalNet& operator =(DigitalNet &&)      = default;
		
		/// Creates empty object
		DigitalNet(void);
		
		/// Creates digital net with given generating numbers
		DigitalNet(std::vector<GenNum> const &generating_numbers);
		
		/// Creates digital net with given generating matrices
		DigitalNet(std::vector<GenMat> const &generating_matrices);
		
		virtual ~DigitalNet(void);
		
		/// Returns \f$m\f$ parameter of the net
		BasicInt m(void) const;
		
		/// Returns \f$s\f$ parameter of the net
		BasicInt s(void) const;
		
		/** Returns generating numbers corresponging to certain dimension
		 *  @param dim – dimension */
		GenNum   generating_numbers(BasicInt dim) const;
		
		/** Returns generating matrix corresponding to certain dimetnsion
		 *  @param dim – dimension */
		GenMat   generating_matrix(BasicInt dim) const;
		
		/** Generates point of a digital net with the certain number
		 *  @param pos - number of generated point */
		Point    generate_point_classical(CountInt pos) const;
		
		/** Generates point of a digital net with certain Gray's code number
		 *  @param pos - sequence number of generated point */
		Point    generate_point(CountInt pos) const;
		
		/** Generates scaled point of a digital net with certain Gray's code number
		 *  @param pos - sequence number of scaled generated point */
		IntPoint generate_int_point(CountInt pos) const;
		
		/** Sequentially generates a section of reordered net points and applies the handler function to each pair:
		 *  (point, point's number)
		 *  @param handler - handler function to apply
		 *  @param amount - amount of points in the section of the net
		 *  @param pos - number of the first point in the section of the net */
		void        for_each_point    (std::function<void (Point const &, CountInt)> handler,
									   CountInt                                      amount,
									   CountInt                                      pos = 0) const;
		
		/** Sequentially generates a section of reordered scaled net points and applies the handler function to each pair:
		 *  (point, point's number)
		 *  @param handler - handler function to apply
		 *  @param amount - amount of points in the section of the net
		 *  @param pos - number of the first point in the section of the net */
		void        for_each_int_point(std::function<void (IntPoint const &, CountInt)> handler,
									   CountInt                                         amount,
									   CountInt                                         pos = 0) const;
		
		/** Casts scaled integer point to a point by multiplying it by \f$2^{-m}\f$
		 *  @param int_point - point to cast */
		Point cast_int_point_to_real(IntPoint const &int_point) const;
		
		
	protected:
		
		/// \f$m\f$ parameter of the digital net
		BasicInt m_nbits;
		/// \f$s\f$ parameter of the digital net
		BasicInt m_dim;
		/// Coefficient equal to \f$2^{-m}\f$
		Real     m_recip;
		/// Vector of a generating numbers of the digital net
		std::vector<GenNum> m_generating_numbers;
		
		/**
		 */
		DigitalNet(BasicInt                   nbits,
				   BasicInt                   dim,
				   std::vector<GenNum> const &generating_numbers);
		
		/** Stores into the integer vector scaled (t,m,s)-net point with certain number, enumerated according to Gray's code.
		 *  @param [out] point - storage vector
		 *  @param [in] pos - generated scaled (t,m,s)-net point number */
		void  store_int_point     (IntPoint &point,
								   CountInt  pos)   const;
		
		/** Stores into the integer vector scaled (t,m,s)-net point with certain number, enumerated according to Gray's code,
		 *  computed using the previous point.
		 *  @param [out] point - storage vector
		 *  @param [in] pos - generated scaled net point number (should be greater then 0)
		 *  @param [in] prev_point - scaled net point with the previous point number */
		void  store_next_int_point(IntPoint       &point,
								   CountInt        pos,
								   IntPoint const &prev_point) const;
		
	};






	inline BasicInt
	DigitalNet::m(void) const
	{ return m_nbits; }
	
	inline BasicInt
	DigitalNet::s(void) const
	{ return m_dim; }
	
	inline GenNum
	DigitalNet::generating_numbers(BasicInt dim) const
	{ return m_generating_numbers[dim]; }
	
	inline GenMat
	DigitalNet::generating_matrix(BasicInt dim) const
	{ return GenMat(m_generating_numbers[dim]); }
	
}


#endif

