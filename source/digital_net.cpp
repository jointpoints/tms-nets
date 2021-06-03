#include "../include/tms-nets/digital_net.hpp"


namespace tms
{
	
	DigitalNet::DigitalNet(void) :
	    m_nbits(0),
	    m_dim(0),
	    m_recip(1),
	    m_generating_numbers()
	{}
	
	DigitalNet::DigitalNet(std::vector<GenNum> const &generating_numbers) :
	    m_nbits(generating_numbers.empty() ? 0 : generating_numbers[0].size()),
	    m_dim(static_cast<BasicInt>(generating_numbers.size())),
	    m_recip( pow(2, -static_cast<Real>(m_nbits)) ),
	    m_generating_numbers(generating_numbers)
	{
		if ( !generating_numbers.empty() && \
			 !std::all_of(generating_numbers.begin(),
						  generating_numbers.end(),
						  [&](GenNum const &cgenmat) { return cgenmat.size() == generating_numbers[0].size(); } ) )
		{
			throw std::logic_error("\nDirection numbers have different sizes\n");
		}
	}
	
	DigitalNet::DigitalNet(std::vector<GenMat> const &generating_matrices) :
		m_nbits(generating_matrices.empty() ? 0 : generating_matrices[0].size()),
	    m_dim(static_cast<BasicInt>(generating_matrices.size())),
	    m_recip( pow(2, -static_cast<Real>(m_nbits)) ),
	    m_generating_numbers(m_dim)
	{
		if ( !generating_matrices.empty() && \
			 std::all_of(generating_matrices.begin(),
						 generating_matrices.end(),
						 [&](GenMat const &genmat) { return genmat.size() == generating_matrices[0].size(); } ) )
		{
			for (BasicInt i = 0; i < m_dim; ++i)
			{
				m_generating_numbers[i] = (GenNum)generating_matrices[i];
			}
		}
		else if ( !generating_matrices.empty() )
		{
			throw std::logic_error("\nGenerating matrices have different sizes\n");
		}
	}
	
	DigitalNet::~DigitalNet(void)
	{}
	
	
	Point
	DigitalNet::generate_point_classical(CountInt pos) const
	{
		Point point(m_dim, 0);
		for (int i = 0; i < m_dim; ++i)
		{
			uintmax_t acc = 0;
			for (int k = 0; k < m_nbits; ++k)
			{
				acc ^= m_generating_numbers[i][k] * ((pos >> k) & 1);
			}
			point[i] = static_cast<Real>(acc) * m_recip;
		}
		return point;
	}
	
	Point
	DigitalNet::generate_point(CountInt pos) const
	{
		//std::cout << "DERgen called\n";
		IntPoint int_point(m_dim);
		store_int_point(int_point, pos);
		
		return cast_int_point_to_real(int_point);
	}
	
	IntPoint
	DigitalNet::generate_int_point(CountInt pos) const
	{
		IntPoint int_point(m_dim);
		store_int_point(int_point, pos);
		
		return int_point;
	}
	
	void
	DigitalNet::for_each_point(std::function<void (Point const &, CountInt)> handler,
							   CountInt amount,
							   CountInt pos) const
	{
		if ( amount != 0 )
		{
			IntPoint curr_int(m_dim);
			store_int_point(curr_int, pos);
			handler(cast_int_point_to_real(curr_int), pos);
			while ( --amount )
			{
				++pos;
				store_next_int_point(curr_int, pos, curr_int);
				handler(cast_int_point_to_real(curr_int), pos);
			}
		}
	}
	
	void
	DigitalNet::for_each_int_point(std::function<void (IntPoint const &, CountInt)> handler,
								   CountInt amount,
								   CountInt pos) const
	{
		if ( amount != 0 )
		{
			IntPoint curr_int(m_dim);
			store_int_point(curr_int, pos);
			handler(curr_int, pos);
			while ( --amount )
			{
				++pos;
				store_next_int_point(curr_int, pos, curr_int);
				handler(curr_int, pos);
			}
		}
	}
	
	Point
	DigitalNet::cast_int_point_to_real(IntPoint const &int_point) const
	{
		Point point_real(m_dim);
		for (BasicInt i = 0; i < m_dim; ++i)
		{
			point_real[i] = static_cast<Real>(int_point[i])*m_recip;
		}
		return point_real;
	}
	
	
	
	DigitalNet::DigitalNet(BasicInt nbits,
						   BasicInt dim,
						   std::vector<GenNum> const &generating_numbers) :
	    m_nbits(nbits),
	    m_dim(dim),
	    m_recip( pow(2, -static_cast<Real>(m_nbits)) ),
	    m_generating_numbers(generating_numbers)
	{}
	
	void
	DigitalNet::store_int_point(IntPoint &point,
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
					point[i] ^= m_generating_numbers[i][k];
				}
			}
			pos_gray_code >>= 1;
		}
	}
	
	void
	DigitalNet::store_next_int_point(IntPoint       &point,
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
			point[i] = prev_point[i] ^ m_generating_numbers[i][rightmost_zero_bit_pos];
		}
	}

};
