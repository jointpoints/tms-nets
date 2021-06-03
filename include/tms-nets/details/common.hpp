#ifndef TMS_NETS_COMMON_HPP
#define TMS_NETS_COMMON_HPP

#include "../thirdparty/irrpoly/gfcheck.hpp"

#include <vector>
#include <stdexcept>


/** @namespace tms
 *  @brief All the entities in the library are defined in this namespace. */
namespace tms
{
	/// A type for integer values that are less than word size (e.g. m, s parameters of the net)
	using BasicInt   = unsigned int;
	/// A type for integer counting values, (e.g. number of net's point)
	using CountInt   = uintmax_t;
	/// A type of coordinates of points of (t, m, s)-net
	using Real	     = long double;
	/// Represents a generating number of a digital (t, m, s)-net
	using GenNumInt = uintmax_t;
	/// Represents a point of a scaled (t, m, s)-net
	using IntPoint   = std::vector<uintmax_t>;
	/// Represents a point of a (t, m, s)-net
	using Point	     = std::vector<Real>;
	/// Represents a polynomial over GF[2]
	using Polynomial = irrpoly::gfpoly;
	
	
	/// Highest allowed m parameter value of created nets, i.e. highest bit depth value
	BasicInt const max_nbits = sizeof(uintmax_t)*8;
	
	
	class GenNum;
	
	class GenMatRow;
	
	class GenMat;
	
	
	/** @class GenNum
	 *  @brief Represents container of generating numbers of a digital net.
	 *         Can be used as a shortened version of generating matrix  */
	class GenNum
	{
	public:
		
		/// Creates empty container with no generating numbers
		GenNum(void);
		GenNum(GenNum const&);
		GenNum(GenNum &&);
		
		~GenNum(void);
		
		GenNum& operator =(GenNum const &);
		GenNum& operator =(GenNum &&);
		
		/** Creates generating numbers with given values
		 *  @param values – vector of generating numbers values */
		GenNum(std::vector<GenNumInt> const &values);
		
		/** Creates certain amount of zero defined generating numbers
		 *  @param amount – amount of generating numbers */
		GenNum(BasicInt amount);
	
		/// Casts generating numbers to the corresponding generating matrix
		explicit operator GenMat(void) const;
		
		/// Checks whether generating numbers container is empty
		bool       empty(void) const;
		
		/// Checks whether the corresponding generating matrix is toeplitz
		bool       is_toeplitz(void) const;
		
		/// Returns amount of generating numbers in the containter
		BasicInt   size(void) const;
		
		/** Returns certain bit of a generating number addressed as an element of generating matrix
		 *  @param i – row number
		 *  @param j – column number */
		bool       get_bit(BasicInt i, BasicInt j) const;
		
		/** Returns certain generating number
		 *  @param n – index of a number */
		uintmax_t  operator[](BasicInt n) const;
		
		/** Sets new bit value of a certain generating number addressed as an element of generating matrix
		 *  @param i – row number
		 *  @param j – column number
		 *  @param value – new bit value */
		void       set_bit(BasicInt i, BasicInt j, bool value);
		
		/** Returns reference to a certain generating number
		 *  @param n – index of a number */
		uintmax_t& operator[](BasicInt n);
		
		GenNum& operator *=(GenNum const& l);
		
		friend bool operator ==(GenNum const &r, GenNum const &l);
		
		
	private:
		
		BasicInt               m_nbits;
		std::vector<GenNumInt> m_numbers;
	};
	
	GenNum operator * (GenNum r, GenNum const &l);
	bool   operator !=(GenNum const &l, GenNum const &r);
	
	
	/** Represents a row of a generating matrix */
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
		
		bool     empty(void) const;
		
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
		
		explicit operator GenNum(void) const;
		
		bool       empty(void) const;
		
		BasicInt   size(void) const;
		GenMatRow  operator [](BasicInt n) const;
		GenMatRow& operator [](BasicInt n);
		
		bool   is_toeplitz(void) const;
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
	
	
	
	
	
	inline bool
	GenNum::empty() const
	{ return m_numbers.empty(); }
	
	inline BasicInt
	GenNum::size(void) const
	{ return m_nbits; }

	inline bool
	GenNum::get_bit(BasicInt i, BasicInt j) const
	{ return (m_numbers[j] >> (m_nbits - 1 - i)) & 1; }

	inline uintmax_t
	GenNum::operator[](BasicInt n) const
	{ return m_numbers[n]; }


	inline void
	GenNum::set_bit(BasicInt i, BasicInt j, bool value)
	{ m_numbers[j] |= static_cast<uintmax_t>(value) << (m_nbits - 1 - i); }

	inline uintmax_t&
	GenNum::operator[](BasicInt n)
	{ return m_numbers[n]; }




	inline bool
	GenMatRow::empty(void) const
	{ return m_values.empty(); }
	
	inline BasicInt
	GenMatRow::size(void) const
	{ return static_cast<BasicInt>(m_values.size()); }

	inline uint8_t
	GenMatRow::operator [](BasicInt n) const
	{ return m_values[n]; }

	inline uint8_t&
	GenMatRow::operator [](BasicInt n)
	{ return m_values[n]; }




	inline bool
	GenMat::empty() const
	{ return m_rows.empty(); }
	
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
