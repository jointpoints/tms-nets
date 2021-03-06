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
		
		DirNum(void);
		DirNum(DirNum const&);
		DirNum(DirNum &&);
		
		~DirNum(void);
		
		DirNum& operator =(DirNum const &);
		DirNum& operator =(DirNum &&);
		
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
		
		GenMatRow(void);
		GenMatRow(GenMatRow const&);
		GenMatRow(GenMatRow &&);
		
		~GenMatRow(void);
		
		GenMatRow& operator =(GenMatRow const &);
		GenMatRow& operator =(GenMatRow &&);
		
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
	
	bool      operator ==(GenMatRow const &l, GenMatRow const &r);
	bool      operator !=(GenMatRow const &l, GenMatRow const &r);
	GenMatRow operator ^ (GenMatRow l, GenMatRow const &r);
	GenMatRow operator <<(GenMatRow l, BasicInt s);
	GenMatRow operator >>(GenMatRow l, BasicInt s);
	GenMatRow operator * (GenMatRow l, bool m);
	
	
	class GenMat
	{
	public:
		
		GenMat(void);
		GenMat(GenMat const&);
		GenMat(GenMat &&);
		
		~GenMat(void);
		
		GenMat& operator =(GenMat const &);
		GenMat& operator =(GenMat &&);
		
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
	bool   operator ==(GenMat const &l, GenMat const &r);
	
	
	std::ostream& operator <<(std::ostream& out, GenMatRow const &row);
	std::ostream& operator <<(std::ostream& out, GenMat const &gen_mat);
	
	
	
	
	
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




	inline BasicInt
	GenMatRow::size(void) const
	{ return static_cast<BasicInt>(m_values.size()); }

	inline uint8_t
	GenMatRow::operator [](BasicInt n) const
	{ return m_values[n]; }

	inline uint8_t&
	GenMatRow::operator [](BasicInt n)
	{ return m_values[n]; }




	inline BasicInt
	GenMat::size(void) const
	{ return m_nbits; }

	inline GenMatRow
	GenMat::operator [](BasicInt n) const
	{ return m_rows[n]; }

	inline GenMatRow&
	GenMat::operator [](BasicInt n)
	{ return m_rows[n]; }
	
}


#endif
