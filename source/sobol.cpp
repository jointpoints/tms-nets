#include "../include/tms-nets/sobol.hpp"


namespace tms
{
	
	Sobol::Sobol(void) :
	    Niederreiter()
	{}
	
	Sobol::Sobol(BasicInt nbits,
				 BasicInt dim,
				 bool     in_parallel) :
	Niederreiter(nbits,
	             dim,
			     std::vector<GenNum>(dim, GenNum(nbits)),
				 0,
				 (in_parallel ? tms::gf2poly::generate_irrpolys_in_parallel : tms::gf2poly::generate_irrpolys)(dim, nbits),
				 &Sobol::check_init1,
				 (void *)0)
	{
		initialize_generating_numbers();
	}
	
	Sobol::Sobol(BasicInt              const  nbits,
				 std::vector<BasicInt> const &degrees_of_irrpolys) :
	Niederreiter(nbits,
	             static_cast<BasicInt>(degrees_of_irrpolys.size()),
				 std::vector<GenNum>(degrees_of_irrpolys.size(), GenNum(nbits)),
				 0,
				 gf2poly::generate_irrpolys_with_degrees(degrees_of_irrpolys, nbits),
				 &Sobol::check_init2,
				 (void *)&degrees_of_irrpolys)
	{
		initialize_generating_numbers();
	}
	
	Sobol::Sobol(BasicInt                        const  nbits,
				 std::initializer_list<BasicInt> const &degrees_of_irrpolys) :
	    Sobol(nbits, std::vector<BasicInt>{degrees_of_irrpolys})
	{}
	
	Sobol::Sobol(BasicInt                              const  nbits,
				 std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs) :
	Niederreiter(nbits,
	             static_cast<BasicInt>(irrpolys_coeffs.size()),
				 std::vector<GenNum>(irrpolys_coeffs.size(), GenNum(nbits)),
				 0,
				 std::vector<Polynomial>(),
				 &Sobol::check_init3,
				 (void *)&irrpolys_coeffs)
	{
		initialize_generating_numbers();
	}
	
	Sobol::Sobol(BasicInt                                        const  nbits,
				 std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs) :
		Sobol(nbits, std::vector< std::vector<uintmax_t> >{irrpolys_coeffs})
	{}
	
	
	
	void
	Sobol::initialize_generating_numbers(void)
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
						m_generating_numbers[i].set_bit(j, k, alpha[k + r]);
						//m_generating_numbers[i][k] |= static_cast<GenNumInt>(alpha[k + r]) << (m_nbits - 1 - j);
					}
					++j;
					--rows_remaining_in_section;
				}
			}
		}
	}
	
	GenNum
	Sobol::inversed_generating_numbers(BasicInt dim) const
	{
		GenNum dir_num(m_nbits);
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

};
