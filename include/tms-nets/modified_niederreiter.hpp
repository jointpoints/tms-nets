#ifndef TMS_NETS_MODIFIED_NIEDERREITER_HPP
#define TMS_NETS_MODIFIED_NIEDERREITER_HPP

#include "niederreiter.hpp"

namespace tms
{
	class ModifiedNiederreiter : public Niederreiter
	{
	public:
		
		ModifiedNiederreiter(BasicInt nbits,
							 BasicInt dim,
							 bool     in_parallel = false);
		/** Constructs the generator of (t,m,s)-net with specified m, s degrees of initial irreducible polynomials,
		 *  and (optional, for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] degrees_of_irrpolys - vector of the s degrees of initial irreducible polynomials*/
		ModifiedNiederreiter(BasicInt                     nbits,
							 std::vector<BasicInt> const &degrees_of_irrpolys);
		/** Constructs the generator of (t,m,s)-net with specified m, s degrees of initial irreducible polynomials,
		 *  and (optional, for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] degrees_of_irrpolys - vector of the s degrees of initial irreducible polynomials*/
		ModifiedNiederreiter(BasicInt                               nbits,
							 std::initializer_list<BasicInt> const &degrees_of_irrpolys);
		/** Constructs the generator of (t,m,s)-net with specified m, s initial irreducible polynomials, and (optional,
		 *  for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] irrpolys_coeffs - initializer list of the s initial irreducible polynomials coefficients*/
		ModifiedNiederreiter(BasicInt                                     nbits,
							 std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs);
		/** Constructs the generator of (t,m,s)-net with specified m, s initial irreducible polynomials, and (optional,
		 *  for advanced users) initial values of ALL recursive sequences, defining the generation matrices.
		 *  @param [in] nbits - m parameter of the net
		 *  @param [in] irrpolys_coeffs - initializer list of the s initial irreducible polynomials coefficients*/
		ModifiedNiederreiter(BasicInt                                               nbits,
							 std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs);
		
		DirNum get_inversed_direction_numbers(BasicInt dim) const;
		
		GenMat get_inversed_generating_matrix(BasicInt dim) const;
		
		
	protected:
		
		void initialize_direction_numbers(void) override;
	};
	
	
	
	
	
	
	ModifiedNiederreiter::ModifiedNiederreiter(BasicInt nbits,
											   BasicInt dim,
											   bool     in_parallel) :
	Niederreiter(nbits,
				 dim,
				 std::vector<DirNum>(dim, DirNum(nbits)),
				 0,
				 (in_parallel ? tms::gf2poly::generate_irrpolys_in_parallel : tms::gf2poly::generate_irrpolys)(dim, nbits),
				 &ModifiedNiederreiter::check_init1,
				 (void *)0)
	{
		initialize_direction_numbers();
	}
	
	ModifiedNiederreiter::ModifiedNiederreiter(BasicInt              const  nbits,
											   std::vector<BasicInt> const &degrees_of_irrpolys) :
	Niederreiter(nbits,
				 static_cast<BasicInt>(degrees_of_irrpolys.size()),
				 std::vector<DirNum>(degrees_of_irrpolys.size(), DirNum(nbits)),
				 0,
				 gf2poly::generate_irrpolys_with_degrees(degrees_of_irrpolys, nbits),
				 &ModifiedNiederreiter::check_init2,
				 (void *)&degrees_of_irrpolys)
	{
		initialize_direction_numbers();
	}
	
	ModifiedNiederreiter::ModifiedNiederreiter(BasicInt                        const  nbits,
											   std::initializer_list<BasicInt> const &degrees_of_irrpolys) :
	ModifiedNiederreiter(nbits, std::vector<BasicInt>{degrees_of_irrpolys})
	{}
	
	ModifiedNiederreiter::ModifiedNiederreiter(BasicInt                              const  nbits,
											   std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs) :
	Niederreiter(nbits,
				 static_cast<BasicInt>(irrpolys_coeffs.size()),
				 std::vector<DirNum>(irrpolys_coeffs.size(), DirNum(nbits)),
				 0,
				 std::vector<Polynomial>(),
				 &ModifiedNiederreiter::check_init3,
				 (void *)&irrpolys_coeffs)
	{
		initialize_direction_numbers();
	}
	
	ModifiedNiederreiter::ModifiedNiederreiter(BasicInt                                        const  nbits,
											   std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs) :
	ModifiedNiederreiter(nbits, std::vector< std::vector<uintmax_t> >{irrpolys_coeffs})
	{}
	
	
	
	void
	ModifiedNiederreiter::initialize_direction_numbers(void)
	{
		std::vector<BasicInt> alpha(m_nbits + \
									std::max_element(
													 m_irrpolys.begin(), \
													 m_irrpolys.end(), \
													 [](Polynomial const &lpoly, Polynomial const &rpoly) { return lpoly.degree() < rpoly.degree(); } \
													 )->degree() - 1);
		
		for (BasicInt i = 0; i < m_dim; ++i)
		{
			BasicInt const e       = static_cast<BasicInt>(m_irrpolys[i].degree());
			BasicInt const r_nbits = m_nbits % e;
			
			Polynomial poly_mu(tms::gf2poly::make_gf2poly({1}));
			
			alpha.resize(m_nbits - 1 + e);
			
			for (BasicInt j = 0; j < m_nbits; )
			{
				BasicInt rows_remaining_in_section = ( (j/e + 1)*e > m_nbits ) ? r_nbits : e;
				
				poly_mu = poly_mu * m_irrpolys[i];
				
				
				recseq::fill_vector_recursively(alpha, 1ULL << ((j/e + 1)*e - 1), poly_mu);
				
				// Here we interpret the j-th digit of direction number g[i](k) as gamma[i](j,k) - an
				// element of i-th generating matrix Gamma[i]
				while ( rows_remaining_in_section != 0 )
				{
					for (BasicInt k = 0, r = e - 1 - (j % e);
						 k < m_nbits;
						 ++k)
					{
						m_direction_numbers[i].set_bit(j, k, alpha[k + r]);
						//m_direction_numbers[i][k] |= static_cast<DirNumInt>(alpha[k + r]) << (m_nbits - 1 - j);
					}
					++j;
					--rows_remaining_in_section;
				}
			}
		}
	}
	
	DirNum
	ModifiedNiederreiter::get_inversed_direction_numbers(BasicInt dim) const
	{
		DirNum dir_num(m_nbits);
		Polynomial omega = gf2poly::make_gf2poly({0, 1});
		
		Polynomial const &poly = m_irrpolys[dim];
		
		auto exp_poly = [&](irrpoly::gfpoly const &poly, BasicInt degree) -> Polynomial {
			Polynomial cur_poly = gf2poly::make_gf2poly({1});
			while ( degree > 0 )
			{
				cur_poly = cur_poly * poly;
				--degree;
			}
			return cur_poly;
		};
		
		Polynomial cur_poly = gf2poly::make_gf2poly({0, 1});
		for (BasicInt i = 0; i < m_nbits; ++i)
		{
			cur_poly = exp_poly(omega, i % poly.degree())*exp_poly(poly, i/poly.degree());
			for (BasicInt j = 0; j < cur_poly.value().size(); ++j)
			{
				dir_num[i] |= (cur_poly.value()[j]) << (m_nbits - 1 - j);
			}
		}
		
		BasicInt const e = static_cast<BasicInt>(poly.degree());
		
		//		BasicInt q = 0;
		//		while ( q < (m_nbits - 1)/e )
		//		{
		//			for (BasicInt r = e - 1; r > 0; --r)
		//			{
		//				for (BasicInt ri = 0; ri < r; ++ri)
		//				{
		//					dir_num[q*e + r] ^= poly[e - r + ri]*dir_num[q*e + ri];
		//				}
		//			}
		//			++q;
		//		}
		//
		//		for (BasicInt r = (m_nbits - 1) % e; r > 0; --r)
		//		{
		//			for (BasicInt ri = 0; ri < r; ++ri)
		//			{
		//				dir_num[q*e + r] ^= poly[e - r + ri]*dir_num[q*e + ri];
		//			}
		//		}
		
		int k = m_nbits - 1;
		
		while ( k >= 0 )
		{
			int const r_max = (k % e) + 1;
			
			for (int r = r_max - 1; r > 0; --r)
			{
				for (int ri = 0; ri < r; ++ri)
				{
					dir_num[k - r_max + 1 + r] ^= poly[e - r + ri]*dir_num[k - r_max + 1 + ri];
				}
			}
			k -= r_max;
		}
		
		return dir_num;
	}
	
	inline GenMat
	ModifiedNiederreiter::get_inversed_generating_matrix(BasicInt dim) const
	{ return GenMat(get_inversed_direction_numbers(dim)); }
	
}

#endif
