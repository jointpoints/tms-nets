#include "../../include/tms-nets/details/common.hpp"


namespace tms
{
	
	// class GenNum
	
	GenNum::GenNum(void) = default;
	GenNum::GenNum(GenNum const&) = default;
	GenNum::GenNum(GenNum &&)  = default;
	
	GenNum::~GenNum(void) = default;
	
	GenNum& GenNum::operator =(GenNum const &) = default;
	GenNum& GenNum::operator =(GenNum &&)      = default;
	
	
	GenNum::GenNum(std::vector<uintmax_t> const &values) :
	    m_nbits(static_cast<BasicInt>(values.size())),
	    m_numbers(values)
	{
		if ( m_nbits > max_nbits )
		{
			throw std::length_error("\nGenNum can't hold more than " + std::to_string(max_nbits) + " elements\n");
		}
	}
	
	GenNum::GenNum(BasicInt size) :
	    m_nbits(size > max_nbits ? 0 : size),
	    m_numbers(size)
	{
		if ( m_nbits > max_nbits )
		{
			throw std::length_error("\nGenNum can't hold more than " + std::to_string(max_nbits) + " elements\n");
		}
	}
	
	GenNum::operator GenMat(void) const
	{
		GenMat gamma_matrix(m_nbits);
		for (BasicInt j = 0; j < m_nbits; ++j)
		{
			for (BasicInt k = 0; k < m_nbits; ++k)
			{
				gamma_matrix[j][k] = this->get_bit(j, k);
			}
		}
		
		return gamma_matrix;
	}
	
	bool
	GenNum::is_toeplitz(void) const
	{	
		BasicInt i = m_numbers.empty() ? 0 : m_nbits - 1;
		while ( i > 0 && ((m_numbers[i] << 1) & ((1 << m_nbits) - 1)) == m_numbers[i - 1] )
		{
			--i;
		}
		
		return i == 0 && !m_numbers.empty(); 
	}
	
	GenNum&
	GenNum::operator *=(GenNum const &l)
	{
		GenNum c = *this;
		*this = GenNum(m_nbits);
		for (BasicInt j = 0; j < m_nbits; ++j)
		{
			for (BasicInt i = 0; i < m_nbits; ++i)
			{
				m_numbers[j] ^= c.m_numbers[i] * l.get_bit(i, j);
			}
		}
		
		return *this;
	}
	
	bool
	operator ==(GenNum const &l, GenNum const &r)
	{
		return l.m_numbers == r.m_numbers;
	}
	
	
	GenNum
	operator * (GenNum r, GenNum const &l)
	{
		return r *= l;
	}
	
	bool
	operator !=(GenNum const &l, GenNum const &r)
	{
		return !(l == r);
	}
	
	
	
	
	
	
	// class GenMatRow
	
	GenMatRow::GenMatRow(void) = default;
	GenMatRow::GenMatRow(GenMatRow const&) = default;
	GenMatRow::GenMatRow(GenMatRow &&)     = default;
	
	GenMatRow::~GenMatRow(void) = default;
	
	GenMatRow& GenMatRow::operator =(GenMatRow const &) = default;
	GenMatRow& GenMatRow::operator =(GenMatRow &&)      = default;
	
	
	GenMatRow::GenMatRow(BasicInt size) :
	    m_values(std::vector<uint8_t>(size > max_nbits ? 0 : size))
	{
		if ( size > max_nbits )
		{
			throw std::length_error("\nGenMatRow can't hold more than " + std::to_string(max_nbits) + " elements\n");
		}
	}
	
	GenMatRow::GenMatRow(std::vector<uint8_t> const &values_vector) :
	    m_values(values_vector)
	{
		if ( values_vector.size() > max_nbits )
		{
			throw std::length_error("\nGenMatRow can't hold more than " + std::to_string(max_nbits) + " elements\n");
		}
	}
	
	GenMatRow::GenMatRow(std::initializer_list<uint8_t> const &values_list) :
	    GenMatRow(std::vector<uint8_t>(values_list))
	{}
	
	
	GenMatRow&
	GenMatRow::operator ^=(GenMatRow const &r)
	{
		if ( m_values.size() == r.size() )
		{
			for (BasicInt i = 0; i < m_values.size(); ++i)
			{
				m_values[i] ^= r[i];
				m_values[i] &= 1;
			}
			
			return *this;
		}
		else
		{
			throw std::length_error("\nArguments sizes mismatching\n");
		}
	}
	
	GenMatRow&
	GenMatRow::operator >>=(BasicInt s)
	{
		BasicInt const size = static_cast<BasicInt>(m_values.size());
		s = std::min(s, size);
		std::move(m_values.data(), m_values.data() + size - s, m_values.data() + s);
		std::fill(m_values.data(), m_values.data() + s, 0);
		
		return *this;
	}
	
	GenMatRow&
	GenMatRow::operator <<=(BasicInt s)
	{
		BasicInt const size = static_cast<BasicInt>(m_values.size());
		s = std::min(s, size);
		std::move(m_values.data() + s, m_values.data() + size, m_values.data());
		std::fill(m_values.data() + size - s, m_values.data() + size, 0);
		
		return *this;
	}
	
	GenMatRow&
	GenMatRow::operator *=(bool m)
	{
		uint8_t n = static_cast<uint8_t>(m);
		for (uint8_t &value : m_values)
		{
			value &= n;
		}
		
		return *this;
	}
	
	
	bool
	operator ==(GenMatRow const &l, GenMatRow const &r)
	{
		if ( l.size() == r.size() )
		{
			return l.m_values == r.m_values;
		}
		else
		{
			throw std::length_error("\nArguments sizes mismatching\n");
		}
	}
	
	bool
	operator !=(GenMatRow const &l, GenMatRow const &r)
	{
		return !(l == r);
	}
	
	GenMatRow
	operator ^(GenMatRow l, GenMatRow const &r)
	{
		return l ^= r;
	}
	
	GenMatRow
	operator <<(GenMatRow l, BasicInt s)
	{
		return l <<= s;
	}
	
	GenMatRow
	operator >>(GenMatRow l, BasicInt s)
	{
		return l >>= s;
	}
	
	GenMatRow
	operator * (GenMatRow l, bool m)
	{
		return l *= m;
	}
	
	
	
	
	
	// class GenMat
	
	GenMat::GenMat(void) = default;
	GenMat::GenMat(GenMat const&) = default;
	GenMat::GenMat(GenMat &&)     = default;
	
	GenMat::~GenMat(void) = default;
	
	GenMat& GenMat::operator =(GenMat const &) = default;
	GenMat& GenMat::operator =(GenMat &&)      = default;
	
	
	GenMat::GenMat(BasicInt size) :
	    m_nbits(size > max_nbits ? 0 : size),
	    m_rows(m_nbits, GenMatRow(m_nbits))
	{
		if ( size > max_nbits )
		{
			throw std::length_error("\nGenMat can't hold more than " + std::to_string(max_nbits) + " rows\n");
		}
	}
	
	GenMat::GenMat(std::vector<GenMatRow> const &rows_vector) :
	    m_nbits(static_cast<BasicInt>(rows_vector.size())),
	    m_rows(rows_vector)
	{
		if ( m_nbits > max_nbits )
		{
			throw std::length_error("\nGenMat can't hold more than " + std::to_string(max_nbits) + " rows\n");
		}
	}
	
	GenMat::GenMat(std::initializer_list<GenMatRow> const &rows_list) :
	    GenMat(std::vector<GenMatRow>(rows_list))
	{}
	
	GenMat::operator GenNum(void) const
	{
		GenNum gen_num(m_nbits);
		for (BasicInt j = 0; j < m_nbits; ++j)
		{
			for (BasicInt i = 0; i < m_nbits; ++i)
			{
				gen_num.set_bit(i, j, m_rows[i][j] & 1);
			}
		}
		
		return gen_num;
	}
	
	
	bool
	GenMat::is_toeplitz(void) const
	{
		BasicInt i = 1;
		while ( i < m_nbits && \
			   std::equal(m_rows[0].m_values.data(),
						  m_rows[0].m_values.data() + m_nbits - i,
						  m_rows[i].m_values.data() + i) )
		{
			++i;
		}
		
		return i == m_nbits;
	}
	
	GenMat
	GenMat::inverse(void) const
	{
		GenMat matrix = *this;
		
		GenMat ident = GenMat(m_nbits);
		for (BasicInt i = 0; i < m_nbits; ++i)
		{
			ident[i][i] = 1;
		}
		
		//gauss down:
		for (BasicInt base_i = 0; base_i < m_nbits; ++base_i)
		{
			BasicInt row_i = base_i;
			
			while ( matrix[row_i][base_i] == 0 )
			{
				++row_i;
			}
			std::swap(matrix[row_i], matrix[base_i]);
			std::swap( ident[row_i],  ident[base_i]);
			
			for (row_i = base_i + 1; row_i < m_nbits; ++row_i)
			{
				if ( matrix[row_i][base_i] != 0 )
				{
					matrix[row_i] ^= matrix[base_i];
					ident[row_i] ^=  ident[base_i];
				}
			}
		}
		
		//gauss up:
		for (BasicInt base_i = m_nbits - 1; base_i > 0; --base_i)
		{
			for (BasicInt row_i = 0; row_i < base_i; ++row_i)
			{
				if ( matrix[row_i][base_i] != 0 )
				{
					ident[row_i] ^= ident[base_i];
				}
			}
		}
		
		return ident;
	}
	
	
	GenMat&
	GenMat::operator *=(GenMat const &r)
	{
		GenMat c = *this;
		*this = GenMat(m_nbits);
		for (BasicInt i = 0; i < m_nbits; ++i)
		{
			for (BasicInt j = 0; j < m_nbits; ++j)
			{
				for (BasicInt k = 0; k < m_nbits; ++k)
				{
					(*this)[i][j] ^= c[i][k] & r[k][j];
				}
			}
		}
		
		return *this;
	}
	
	GenMat
	operator *(GenMat l, GenMat const &r)
	{
		return l *= r;
	}
	
	bool
	operator ==(GenMat const &a, GenMat const &b)
	{
		return a.m_rows == b.m_rows;
	}
	
	
	
	
	
	// Verbosity for GenMatRow and GenMat
	
	std::ostream& operator <<(std::ostream& out, GenMatRow const &row)
	{
		out << "{";
		if ( !row.empty() )
		{
			out << +row[0];
		}
		for (BasicInt i = 1; i < row.size(); ++i)
		{
			out << ", " << +row[i];
		}
		out << "}";
		
		return out;
	}
	
	std::ostream& operator <<(std::ostream& out, GenMat const &gen_mat)
	{
		out << "{";
		for (BasicInt i = 0; i < gen_mat.size() - 1; ++i)
		{
			out << gen_mat[i] << ",\n";
		}
		out << gen_mat[gen_mat.size() - 1] << "}";
		
		return out;
	}

	
};
