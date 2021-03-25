#ifndef TMS_NETS_DIGITAL_NET_HPP
#define TMS_NETS_DIGITAL_NET_HPP

#include "details/gf2poly.hpp"

#include <vector>
#include <cmath>		//for pow function
#include <algorithm>	//for std::max_element function (is already included in "irrpoly/gfpoly.hpp")
#include <functional>	//for unified for_each_point* methods (is already included in "irrpoly/gfpoly.hpp")


namespace tms
{
	class DigitalNet
	{
	public:
		
		// Prevent implementation by compiler of implicit-defined copy/move constructors, assignment operators
		DigitalNet(DigitalNet const &) = delete;
		DigitalNet(DigitalNet &&)      = delete;
		DigitalNet& operator =(DigitalNet const &) = delete;
		DigitalNet& operator =(DigitalNet&&)       = delete;
		
		
		DigitalNet(void);
		
		DigitalNet(std::vector<DirNum> const &direction_numbers);
		
		DigitalNet(std::vector<GenMat> const &generating_matrices);
		
		virtual ~DigitalNet(void);
		
		/** Returns m parameter of generatred (t,m,s)-net.*/
		BasicInt get_m(void) const;
		/** Returns s parameter of generatred (t,m,s)-net.*/
		BasicInt get_s(void) const;
		/** Returns vector of direction numbers corresponging to certain dimension.
		 *  @param [in] dim - dimension */
		DirNum   get_direction_numbers(BasicInt dim) const;
		
		GenMat   get_generating_matrix(BasicInt dim) const;
		
		/** Generates point of (t,m,s)-net with certain number.
		 *  @param [in] n - sequence number of generated point */
		Point    generate_point_classical(CountInt pos) const;
		/** Generates point of (t,m,s)-net with certain number, enumerated according to Gray's code.
		 *  @param [in] pos - sequence number of generated point */
		Point    generate_point(CountInt pos) const;
		/** Generates scaled point of (t,m,s)-net with certain number, enumerated according to Gray's code.
		 *  @param [in] pos - sequence number of scaled generated point */
		IntPoint generate_point_int(CountInt pos) const;
		/** Sequentially generates a section of reordered (t,m,s)-net points and applies the handler function to each pair:
		 *  (point, point's number).
		 *  @param [in] handler - handler function to apply
		 *  @param [in] amount - amount of points in the section of the net
		 *  @param [in] pos - number of the first point in the section of the net */
		void        for_each_point    (std::function<void (Point const &, CountInt)> handler,
									   CountInt                                      amount,
									   CountInt                                      pos = 0) const;
		/** Sequentially generates a section of reordered scaled (t,m,s)-net points and applies the handler function to each pair:
		 *  (point, point's number).
		 *  @param [in] handler - handler function to apply
		 *  @param [in] amount - amount of points in the section of the net
		 *  @param [in] pos - number of the first point in the section of the net */
		void        for_each_point_int(std::function<void (IntPoint const &, CountInt)> handler,
									   CountInt                                         amount,
									   CountInt                                         pos = 0) const;
		
		/** Returns transformed given integer point into a real point from a unit hypercube.
		 *  @param [in] point_int - point to be transformed (multiplied by \f$ 2^{-m} \f$) */
		Point cast_point_int_to_real(IntPoint const &point_int) const;
		
		
	protected:
		
		/// m parameter of net that defines net cardinality and bitwidth of computations.
		BasicInt m_nbits;
		/// s parameter of net that defines spatial dimensionality.
		BasicInt m_dim;
		/// Coefficient, equal to \f$ 2^{-m} \f$.
		Real     m_recip;
		/// Vector of a (t,m,s)-net's direction numbers.
		std::vector<DirNum> m_direction_numbers;
		
		DigitalNet(BasicInt                   nbits,
				   BasicInt                   dim,
				   std::vector<DirNum> const &direction_numbers);
		/** Stores into the integer vector scaled (t,m,s)-net point with certain number, enumerated according to Gray's code.
		 *  @param [out] point - storage vector
		 *  @param [in] pos - generated scaled (t,m,s)-net point number */
		void  store_point_int     (IntPoint &point,
								   CountInt  pos)   const;
		/** Stores into the integer vector scaled (t,m,s)-net point with certain number, enumerated according to Gray's code,
		 *  computed using the previous point.
		 *  @param [out] point - storage vector
		 *  @param [in] pos - generated scaled (t,m,s)-net point number (should be greater then 0)
		 *  @param [in] prev_point - scaled (t,m,s)-net point with the previous point number */
		void  store_next_point_int(IntPoint       &point,
								   CountInt        pos,
								   IntPoint const &prev_point) const;
	};






	inline BasicInt
	DigitalNet::get_m(void) const
	{ return m_nbits; }
	
	inline BasicInt
	DigitalNet::get_s(void) const
	{ return m_dim; }
	
	inline DirNum
	DigitalNet::get_direction_numbers(BasicInt dim) const
	{ return m_direction_numbers[dim]; }
	
	inline GenMat
	DigitalNet::get_generating_matrix(BasicInt dim) const
	{ return GenMat(m_direction_numbers[dim]); }
	
}


#endif

