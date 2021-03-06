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






	DigitalNet::DigitalNet(void) :
	    m_nbits(0),
	    m_dim(0),
	    m_recip(1),
	    m_direction_numbers()
	{}
	
	DigitalNet::DigitalNet(std::vector<DirNum> const &direction_numbers) :
	    m_nbits(direction_numbers.empty() ? 0 : direction_numbers[0].size()),
	    m_dim(static_cast<BasicInt>(direction_numbers.size())),
	    m_recip( pow(2, -static_cast<Real>(m_nbits)) ),
	    m_direction_numbers(direction_numbers)
	{}
	
	DigitalNet::DigitalNet(std::vector<GenMat> const &generating_matrices) :
	    m_nbits(generating_matrices.empty() ? 0 : generating_matrices[0].size()),
	    m_dim(static_cast<BasicInt>(generating_matrices.size())),
	    m_recip( pow(2, -static_cast<Real>(m_nbits)) ),
		m_direction_numbers(m_dim)
	{
		for (BasicInt i = 0; i < m_dim; ++i)
		{
			m_direction_numbers[i] = (DirNum)generating_matrices[i];
		}
	}
	
	DigitalNet::~DigitalNet(void)
	{}
	
	
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
	
	
	Point
	DigitalNet::generate_point_classical(CountInt pos) const
	{
		Point point(m_dim, 0);
		for (int i = 0; i < m_dim; ++i)
		{
			DirNumInt acc = 0;
			for (int k = 0; k < m_nbits; ++k)
			{
				acc ^= m_direction_numbers[i][k] * ((pos >> k) & 1);
			}
			point[i] = static_cast<Real>(acc) * m_recip;
		}
		return point;
	}
	
	Point
	DigitalNet::generate_point(CountInt pos) const
	{
		//std::cout << "DERgen called\n";
		IntPoint point_int(m_dim);
		store_point_int(point_int, pos);
		
		return cast_point_int_to_real(point_int);
	}
	
	IntPoint
	DigitalNet::generate_point_int(CountInt pos) const
	{
		IntPoint point_int(m_dim);
		store_point_int(point_int, pos);
		
		return point_int;
	}
	
	void
	DigitalNet::for_each_point(std::function<void (Point const &, CountInt)> handler,
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
	
	void
	DigitalNet::for_each_point_int(std::function<void (IntPoint const &, CountInt)> handler,
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
	
	Point
	DigitalNet::cast_point_int_to_real(IntPoint const &point_int) const
	{
		Point point_real(m_dim);
		for (BasicInt i = 0; i < m_dim; ++i)
		{
			point_real[i] = static_cast<Real>(point_int[i])*m_recip;
		}
		return point_real;
	}
	
	
	
	DigitalNet::DigitalNet(BasicInt nbits,
						   BasicInt dim,
						   std::vector<DirNum> const &direction_numbers) :
	m_nbits(nbits),
	m_dim(dim),
	m_recip( pow(2, -static_cast<Real>(m_nbits)) ),
	m_direction_numbers(direction_numbers)
	{}
	
	void
	DigitalNet::store_point_int(IntPoint &point,
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
	
	void
	DigitalNet::store_next_point_int(IntPoint       &point,
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
	
	
}


#endif

