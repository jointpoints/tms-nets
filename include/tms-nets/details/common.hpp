#ifndef TMS_NETS_COMMON_HPP
#define TMS_NETS_COMMON_HPP

#include "../thirdparty/irrpoly/gfcheck.hpp"

#include <vector>
#include <stdexcept>


namespace tms
{
	using DirNumInt  = uintmax_t;
	
	using IntPoint   = std::vector<DirNumInt>;
	/// Type used for storing integer values that are less than word size.
	using BasicInt   = unsigned int;
	/// Type used for storing quantity-values and indexation-values, i.e. the total amount of points to generate.
	using CountInt   = uintmax_t;
	/// Type used for storing coordinates of points of (t,m,s)-net.
	using Real	     = long double;
	/// Type used for storing points of (t,m,s)-net.
	using Point	     = std::vector<Real>;
	/// Type used for storing polynomials.
	using Polynomial = irrpoly::gfpoly;
	
	class DirNum;
	
	class GenMatRow;
	
	class GenMat;
	
	BasicInt const max_nbits = sizeof(DirNumInt)*8;
	
	
	class DirNum
	{
	public:
		
		DirNum(void) = default;
		DirNum(DirNum const&) = default;
		DirNum(DirNum &&)  = default;
		
		~DirNum(void) = default;
		
		DirNum& operator =(DirNum const &) = default;
		DirNum& operator =(DirNum &&)      = default;
		
		DirNum(std::vector<DirNumInt> const &values);
		DirNum(BasicInt size);
		//DirNum(GenMat   const &gen_mat);
	
		explicit operator GenMat(void) const;
		
		BasicInt   size(void) const;
		bool       get_bit(BasicInt i, BasicInt j) const;
		DirNumInt  operator[](BasicInt n) const;
		
		void       set_bit(BasicInt i, BasicInt j, bool value);
		DirNumInt& operator[](BasicInt n);
		
		
	private:
		
		BasicInt               m_nbits;
		std::vector<DirNumInt> m_numbers;
	};
	
	
	class GenMatRow
	{
	public:
		
		GenMatRow(void) = default;
		GenMatRow(GenMatRow const&) = default;
		GenMatRow(GenMatRow &&)     = default;
		
		~GenMatRow(void) = default;
		
		GenMatRow& operator =(GenMatRow const &) = default;
		GenMatRow& operator =(GenMatRow &&)      = default;
		
		GenMatRow(BasicInt size);
		GenMatRow(std::vector<uint8_t> const &values_vector); //?
		GenMatRow(std::initializer_list<uint8_t> const &values_list);
		
		BasicInt size(void) const;
		uint8_t  operator [](BasicInt n) const;
		uint8_t& operator [](BasicInt n);
		
		GenMatRow& operator ^= (GenMatRow const &r);
		GenMatRow& operator >>=(BasicInt s);
		GenMatRow& operator <<=(BasicInt s);
		GenMatRow& operator *= (bool m);

		friend bool operator ==(GenMatRow const &l, GenMatRow const &r);
		
		friend GenMat;

		
	private:
		
		std::vector<uint8_t> m_values;
	};
	
	bool      operator !=(GenMatRow const &l, GenMatRow const &r);
	GenMatRow operator ^ (GenMatRow l, GenMatRow const &r);
	GenMatRow operator <<(GenMatRow l, BasicInt s);
	GenMatRow operator >>(GenMatRow l, BasicInt s);
	GenMatRow operator * (GenMatRow l, bool m);
	
	
	class GenMat
	{
	public:
		
		GenMat(void) = default;
		GenMat(GenMat const&) = default;
		GenMat(GenMat &&)     = default;
		
		~GenMat(void) = default;
		
		GenMat& operator =(GenMat const &) = default;
		GenMat& operator =(GenMat &&)      = default;
		
		GenMat(BasicInt size);
		GenMat(std::vector<GenMatRow> const &matrix);
		GenMat(std::initializer_list<GenMatRow> const &lines_list);
		
		explicit operator DirNum(void) const;
		
		BasicInt size(void) const;
		GenMatRow  operator [](BasicInt n) const;
		GenMatRow& operator [](BasicInt n);
		
		bool   is_shifted(void) const;
		GenMat inverse(void) const;
		
		GenMat& operator *=(GenMat const &r);
		
		friend bool operator ==(GenMat const &l, GenMat const &r);
		
		
	private:
		
		BasicInt               m_nbits;
		std::vector<GenMatRow> m_rows;
		
	};
	
	GenMat operator * (GenMat l, GenMat const &r);
	
	
	std::ostream& operator <<(std::ostream& out, GenMatRow const &row);
	std::ostream& operator <<(std::ostream& out, GenMat const &gen_mat);
	
	
	
	
	
	
	DirNum::DirNum(std::vector<DirNumInt> const &values) :
	m_nbits(values.size() > max_nbits ? 0 : static_cast<BasicInt>(values.size())),
	m_numbers(values)
	{
		if ( m_nbits == 0 )
		{
			throw std::length_error("\n");
		}
	}

	DirNum::DirNum(BasicInt size) :
	m_numbers(size > max_nbits ? 0 : size),
	m_nbits(size)
	{
		if ( m_nbits == 0 )
		{
			throw std::length_error("\n");
		}
	}

	DirNum::operator GenMat(void) const
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


	inline BasicInt
	DirNum::size(void) const
	{ return m_nbits; }

	inline bool
	DirNum::get_bit(BasicInt i, BasicInt j) const
	{ return (m_numbers[j] >> (m_nbits - 1 - i)) & 1; }

	inline DirNumInt
	DirNum::operator[](BasicInt n) const
	{ return m_numbers[n]; }


	inline void
	DirNum::set_bit(BasicInt i, BasicInt j, bool value)
	{ m_numbers[j] |= static_cast<DirNumInt>(value) << (m_nbits - 1 - i); }

	inline DirNumInt&
	DirNum::operator[](BasicInt n)
	{ return m_numbers[n]; }




	GenMatRow::GenMatRow(BasicInt size) :
	m_values(std::vector<uint8_t>(size > max_nbits ? 0 : size))
	{
		if ( size > max_nbits )
		{
			throw std::length_error("\n");
		}
	}

	GenMatRow::GenMatRow(std::vector<uint8_t> const &values_vector) :
	m_values(values_vector)
	{
		if ( values_vector.size() > max_nbits )
		{
			throw std::length_error("\n");
		}
	}

	GenMatRow::GenMatRow(std::initializer_list<uint8_t> const &values_list) :
	GenMatRow(std::vector<uint8_t>(values_list))
	{}


	inline BasicInt
	GenMatRow::size(void) const
	{ return static_cast<BasicInt>(m_values.size()); }

	inline uint8_t
	GenMatRow::operator [](BasicInt n) const
	{ return m_values[n]; }

	inline uint8_t&
	GenMatRow::operator [](BasicInt n)
	{ return m_values[n]; }


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
			throw std::length_error("\n");
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


	bool operator ==(GenMatRow const &l, GenMatRow const &r)
	{
		if ( l.size() == r.size() )
		{
			return l.m_values == r.m_values;
		}
		else
		{
			throw std::length_error("\n");
		}
	}

	bool operator !=(GenMatRow const &l, GenMatRow const &r)
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




	GenMat::GenMat(BasicInt size) :
	m_nbits(size > max_nbits ? 0 : size),
	m_rows(m_nbits, GenMatRow(m_nbits))
	{
		if ( m_nbits == 0 )
		{
			throw std::length_error("\n");
		}
	}

	GenMat::GenMat(std::vector<GenMatRow> const &rows_vector) :
	m_nbits(rows_vector.size() > max_nbits ? 0 : static_cast<BasicInt>(rows_vector.size())),
	m_rows(rows_vector)
	{
		if ( m_nbits == 0 )
		{
			throw std::length_error("\n");
		}
	}

	GenMat::GenMat(std::initializer_list<GenMatRow> const &rows_list) :
	GenMat(std::vector<GenMatRow>(rows_list))
	{}

	GenMat::operator DirNum(void) const
	{
		DirNum dir_num(m_nbits);
		for (BasicInt j = 0; j < m_nbits; ++j)
		{
			for (BasicInt i = 0; i < m_nbits; ++i)
			{
				dir_num.set_bit(i, j, m_rows[i][j] & 1);
			}
		}

		return dir_num;
	}


	inline BasicInt
	GenMat::size(void) const
	{ return m_nbits; }

	inline GenMatRow
	GenMat::operator [](BasicInt n) const
	{ return m_rows[n]; }

	inline GenMatRow&
	GenMat::operator [](BasicInt n)
	{ return m_rows[n]; }

	bool
	GenMat::is_shifted(void) const
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




	std::ostream& operator <<(std::ostream& out, GenMatRow const &row)
	{
		out << "{" << +row[0];
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
	
}


#endif
