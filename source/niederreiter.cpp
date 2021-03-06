#include "../include/tms-nets/niederreiter.hpp"





tms::Niederreiter::~Niederreiter(void) = default;


tms::Niederreiter::Niederreiter(BasicInt const nbits,
                                BasicInt const dim,
                                bool     const in_parallel) :
	DigitalNet(nbits, dim, std::vector<DirNum>(dim, DirNum(nbits))),
	m_irrpolys( (in_parallel ? tms::gf2poly::generate_irrpolys_in_parallel : tms::gf2poly::generate_irrpolys)(dim, nbits) )
{
	check_init1(static_cast<void *>(0));
	initialize_direction_numbers();
}

tms::Niederreiter::Niederreiter(BasicInt              const  nbits,
                                std::vector<BasicInt> const &degrees_of_irrpolys) :
	DigitalNet(nbits,
	           static_cast<BasicInt>(degrees_of_irrpolys.size()),
	           std::vector<DirNum>(degrees_of_irrpolys.size(), DirNum(nbits))),
	m_irrpolys{gf2poly::generate_irrpolys_with_degrees(degrees_of_irrpolys, nbits)}
{
	check_init2((void *)(&degrees_of_irrpolys));
	initialize_direction_numbers();
}

tms::Niederreiter::Niederreiter(BasicInt                        const  nbits,
                                std::initializer_list<BasicInt> const &degrees_of_irrpolys) :
	Niederreiter(nbits,
	             std::vector<BasicInt>{degrees_of_irrpolys})
{}

tms::Niederreiter::Niederreiter(BasicInt                              const  nbits,
                                std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs) :
	DigitalNet(nbits,
	           static_cast<BasicInt>(irrpolys_coeffs.size()),
	           std::vector<DirNum>(irrpolys_coeffs.size(), DirNum(nbits)))
{
	check_init3((void *)(&irrpolys_coeffs));
	initialize_direction_numbers();
}

tms::Niederreiter::Niederreiter(BasicInt                                        const  nbits,
                                std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs) :
	Niederreiter(nbits,
	             std::vector< std::vector<uintmax_t> >{irrpolys_coeffs})
{}


tms::BasicInt
tms::Niederreiter::get_t_estimate(void) const
{
	BasicInt t = 0;
	std::for_each(m_irrpolys.begin(), m_irrpolys.end(), [&](Polynomial const &poly) { t += poly.degree() - 1; });
	return t;
}



tms::Niederreiter::Niederreiter(BasicInt                       nbits,
                                BasicInt                       dim,
                                std::vector<DirNum>     const &direction_numbers,
                                Real                           recip,
                                std::vector<Polynomial> const &irrpolys,
                                void           (Niederreiter::*ptr_check)(void*),
                                void                          *ptr_arg) :
DigitalNet(nbits, dim, direction_numbers),
m_irrpolys(irrpolys)
{
	(this->*ptr_check)(ptr_arg);
}

void tms::Niederreiter::check_init1(void *ptr_arg)
{
	if ( m_nbits > sizeof(DirNumInt)*8 )
	{
		throw std::logic_error("\nnbits is too high");
	}
	
	// if m_dim == 0 or
	//  if m (m_nbits) is too small for desired s (m_dim)
	if ( m_dim == 0 || m_irrpolys.size() != m_dim )
	{
		throw std::logic_error("\nWrong net's parameters");
	}
}

void tms::Niederreiter::check_init2(void *ptr_arg)
{
	std::vector<BasicInt> *ptr_degrees_of_irrpolys = static_cast<std::vector<BasicInt>*>(ptr_arg);
	
	if ( m_nbits > sizeof(DirNumInt)*8 )
	{
		throw std::logic_error("\nnbits is too high");
	}
	
	// if degrees_of_irrpolys is empty or
	//  if there is no polynomials with such degrees that induced t (m_quality_param) <= m (nbits)
	if ( ptr_degrees_of_irrpolys->size() == 0 || m_irrpolys.size() != ptr_degrees_of_irrpolys->size() )
	{
		throw std::logic_error("\nWrong polynomial degrees or the nbits parameter");
	}
}

void tms::Niederreiter::check_init3(void *ptr_arg)
{
	std::vector< std::vector<uintmax_t> > *ptr_irrpolys_coeffs = static_cast<std::vector< std::vector<uintmax_t> >*>(ptr_arg);
	
	BasicInt quality_param = 0;
	
	if ( m_nbits > sizeof(DirNumInt)*8 )
	{
		throw std::logic_error("\nnbits is too high");
	}
	
	//checking the polynomials:
	m_irrpolys.reserve(m_dim);
	for (BasicInt i = 0; i < m_dim; ++i)
	{
		m_irrpolys.push_back( gf2poly::make_gf2poly(*(ptr_irrpolys_coeffs->begin() + i)) );
	}
	
	BasicInt i = 0;
	// In order to not include another header, uniqueness of polynomials checked with std::map container.
	// std::map was chosen for two reasons:
	// 1. it has been already included in gf2poly.hpp;
	// 2. we need to keep original order of polynomials.
	std::map<Polynomial, BasicInt, bool(*)(Polynomial const &, Polynomial const &)> just_set(irrpoly::operator!=);
	while ( i < m_irrpolys.size()  &&  just_set.size() == i && \
			quality_param <= m_nbits    &&  irrpoly::is_irreducible_berlekamp(m_irrpolys[i]) )
	{
		just_set.insert(std::make_pair(m_irrpolys[i], i));
		quality_param += m_irrpolys[i].size() - 2;
		++i;
	}
	// if irrpolys.empty() or
	//  if irrpolys contains not unique polynomials
	//  if induced t (m_quality_param) is greater than m (m_nbits)
	if ( m_irrpolys.empty() || just_set.size() != m_irrpolys.size() || quality_param > m_nbits )
	{
		throw std::logic_error("\nWrong polynomials");
	}
}

void
tms::Niederreiter::initialize_direction_numbers(void)
{
	//std::cout << "Classical called\n";
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
			

			recseq::fill_vector_recursively(alpha,
											( r_nbits != 0 && j/e == (m_nbits - 1)/e ) ? 1ULL << (m_nbits - 1) : 1ULL << ((j/e + 1)*e - 1),
											poly_mu);
			
			// Here we interpret the j-th digit of direction number g[i](k) as gamma[i](j,k) - an
			// element of i-th generating matrix Gamma[i]
			while ( rows_remaining_in_section != 0 )
			{
				for (BasicInt k = 0, r = j % e; k < m_nbits; ++k)
				{
					m_direction_numbers[i][k] |= static_cast<DirNumInt>(alpha[k + r]) << (m_nbits - 1 - j);
				}
				++j;
				--rows_remaining_in_section;
			}
		}
	}
}
