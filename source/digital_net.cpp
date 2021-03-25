#include "../include/tms-nets/digital_net.hpp"





tms::DigitalNet::DigitalNet(void) :
	m_nbits(0),
	m_dim(0),
	m_recip(1),
	m_direction_numbers()
{}

tms::DigitalNet::DigitalNet(std::vector<DirNum> const &direction_numbers) :
	m_nbits(direction_numbers.empty() ? 0 : direction_numbers[0].size()),
	m_dim(static_cast<BasicInt>(direction_numbers.size())),
	m_recip( pow(2, -static_cast<Real>(m_nbits)) ),
	m_direction_numbers(direction_numbers)
{}

tms::DigitalNet::DigitalNet(std::vector<GenMat> const &generating_matrices) :
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

tms::DigitalNet::~DigitalNet(void)
{}


tms::Point
tms::DigitalNet::generate_point_classical(CountInt pos) const
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

tms::Point
tms::DigitalNet::generate_point(CountInt pos) const
{
	//std::cout << "DERgen called\n";
	IntPoint point_int(m_dim);
	store_point_int(point_int, pos);
	
	return cast_point_int_to_real(point_int);
}

tms::IntPoint
tms::DigitalNet::generate_point_int(CountInt pos) const
{
	IntPoint point_int(m_dim);
	store_point_int(point_int, pos);
	
	return point_int;
}

void
tms::DigitalNet::for_each_point(std::function<void (Point const &, CountInt)> handler,
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
tms::DigitalNet::for_each_point_int(std::function<void (IntPoint const &, CountInt)> handler,
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

tms::Point
tms::DigitalNet::cast_point_int_to_real(IntPoint const &point_int) const
{
	Point point_real(m_dim);
	for (BasicInt i = 0; i < m_dim; ++i)
	{
		point_real[i] = static_cast<Real>(point_int[i])*m_recip;
	}
	return point_real;
}



tms::DigitalNet::DigitalNet(BasicInt nbits,
						BasicInt dim,
						std::vector<DirNum> const &direction_numbers) :
m_nbits(nbits),
m_dim(dim),
m_recip( pow(2, -static_cast<Real>(m_nbits)) ),
m_direction_numbers(direction_numbers)
{}

void
tms::DigitalNet::store_point_int(IntPoint &point,
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
tms::DigitalNet::store_next_point_int(IntPoint       &point,
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
