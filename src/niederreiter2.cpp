# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "niederreiter2.h"

//! \todo Extend irred list (#1), irred_deg list (#2), set MAXE to max degree of added polynomials (#3)
void calculate_c(//! [in] the dimension of the sequence to be generated
			 int dim,
			 //! [out] the packed values of Niederreiter's \f$c^{(i)}_{jr}\f$ into numbers \f$С^{(i)}_r\f$
			 Next_int cj[DIM_MAX][NBITS])

/*! \brief Computes values of the constants \f$c^{(i)}_{jr}\f$.
 //
 //	\callgraph
 //
 //	\par Discussion
 //		This program calculates the values of the constants \f$c^{(i)}_{jr}\f$ denoted as c(i, j, r).
 //		As far as possible, Niederreiter's notation is used.\n
 //		For each value of i, we first calculate all the corresponding
 //		values of c (all these values are either 0 or 1) and pack them into the 2D-array \b cj thats
 //		denoting \f$С^{(i)}_r\f$ \n
 //		in such a way that \b cj[i][r] holds the values of c(i, j, r)
 //    	for the indicated i and r and for every j from 0 to (NBITS - 1).\n
 //		The most significant bit of \b cj[i][r] (not counting the sign bit) is c(i, 0, r) and the least
 //    	significant bit is c(i, NBITS - 1, r) that is equivalent to\n
 //		\f$С^{(i)}_r = \sum\limits_{j=0}^{R-1} c^{(i)}_{jr}\cdot 2^{R-1-j}\f$ where \f$R\f$ is equal to NBITS.
 //
 // \par Local Parameters
 //    	MACRO \p int \b MAXE, the highest degree among DIM_MAX irreducible polynomials over GF(2)
 //		listed in \b irred.
 //
 //	\copyright This code is distributed under the GNU LGPL license.
 //
 //	\par Modified
 //		29 October 2019
 //
 // \author Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
 //    		C++ version by:\n
 //				John Burkardt (Florida State University, USA, 2003),\n
 //				Yekaterina Listyukhina (Russian Technological University, KMBO-03-16, Russia, 2019),\n
 //				Alexander Smekhov (Russian Technological University, KMBO-03-16, Russia, 2019),\n
 //				Alexey Burimov (Russian Technological University, KMBO-03-16, Russia, 2019).
 //
 // \par Reference
 //    R Lidl, Harald Niederreiter,
 //    Finite Fields,
 //    Cambridge University Press, 1984, page 553.\n
 //    Harald Niederreiter,
 //    Low-discrepancy and low-dispersion sequences,
 //    Journal of Number Theory,
 //    Volume 30, 1988, pages 51-70.
 //
 */
{
# define MAXE 6
	
	int poly_b[MAXDEG + 1];
	int poly_b_deg;
	int e;
	static int irred_polys[DIM_MAX][MAXE + 1] =
	{
		{ 0,1,0,0,0,0,0 },
		{ 1,1,0,0,0,0,0 },
		{ 1,1,1,0,0,0,0 },
		{ 1,1,0,1,0,0,0 },
		{ 1,0,1,1,0,0,0 },
		{ 1,1,0,0,1,0,0 },
		{ 1,0,0,1,1,0,0 },
		{ 1,1,1,1,1,0,0 },
		{ 1,0,1,0,0,1,0 },
		{ 1,0,0,1,0,1,0 },
		{ 1,1,1,1,0,1,0 },
		{ 1,1,1,0,1,1,0 },
		{ 1,1,0,1,1,1,0 },
		{ 1,0,1,1,1,1,0 },
		{ 1,1,0,0,0,0,1 },
		{ 1,0,0,1,0,0,1 },
		{ 1,1,1,0,1,0,1 },
		{ 1,1,0,1,1,0,1 },
		{ 1,0,0,0,0,1,1 },
		{ 1,1,1,0,0,1,1 }
	};
	int irred_polys_deg[DIM_MAX] = { 1, 1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6 };
	int poly_pi[MAXDEG + 1];
	int poly_pi_deg;
	int i, j, r, u;
	int v[NBITS + MAXE + 1];
	
	for (i = 0; i < dim; ++i)
	{
		//
		//  For each dimension, we need to calculate powers of an
		//  appropriate irreducible polynomial:  see Niederreiter
		//  page 65, just below equation (19).
		//
		//  Copy the appropriate irreducible polynomial p_i(x) into pi,
		//  and its degree into e.  Set polynomial b(x) = 1 and copy it into b.
		//  m is the degree of b(x).  Subsequently b(x) will hold higher
		//  powers of p_i(x).
		//
		e = irred_polys_deg[i];
		poly_pi_deg = irred_polys_deg[i];
		
		for (j = 0; j <= poly_pi_deg; ++j)
		{
			poly_pi[j] = irred_polys[i][j];
		}
		
		poly_b_deg = 0;
		poly_b[0] = 1;
		//
		//  Niederreiter (page 56, after equation (7), defines two
		//  variables q(i, j) and u(i, j).  We do not need q explicitly, but we do need u.
		//
		u = 0;
		
		for (j = 0; j < NBITS; ++j)
		{
			//
			//  If u = 0, we need to set B to the next power of p_i(x)
			//  and recalculate {v_n}. This is done by subroutine CALCV2.
			//
			if ( u == 0 )
			{
				calculate_v(NBITS + MAXE, poly_pi_deg, poly_pi, &poly_b_deg, poly_b, v);
			}
			//
			//  Now c is obtained from v.
			//	!!! We can think about c(i0,j0,r) as a (NBITS - j0 - 1)-th bit in binary
			//		representation of a number cj(i0, r).
			//
			for (r = 0; r < NBITS; ++r)
			{
				cj[i][r] |= ((Next_int)v[r + u]) << (NBITS - 1 - j); //???
			}
			//
			//  Increment u. If u = e, then u = 0 and in Niederreiter's
			//  paper q = q + 1.  Here, however, q is not used explicitly.
			//
			u = u + 1;
			if ( u == e )
			{
				u = 0;
			}
		}
	}
	
	return;
# undef MAXE
}

//****************************************************************************80

void calculate_v(//! [in] the size of vector \f$\{v_n\}\f$
			 int size_of_v,
			 //! [in] the degree of \f$p_i(x)\f$
			 int poly_pi_deg,
			 //! [in] the appropriate irreducible polynomial \f$p_i(x)\f$ for the dimension currently being considered
			 int poly_pi[MAXDEG + 1],
			 //! [in, out] pointer to a degree of the polynomial \f$b(x)\f$
			 int *pt_poly_b_deg,
			 /*! [in, out] on input, \f$b(x)\f$ is the polynomial
			  defined in section 2.3 of BFN. The degree \f$deg(b(x))\f$ implicitly defines
			  the parameter j of section 3.3, by \f$deg(b(x)) = e \cdot (j-1)\f$.  On output,
			  \f$b(x)\f$ has been multiplied by \f$p_i(x)\f$, so its degree is now \f$e \cdot j\f$*/
			 int poly_b[MAXDEG + 1],
			 //! [out] the computed \f$\{v_n\}\f$ array
			 int v[])

/*! \brief Calculates the value of the constants V(J,R).
 //
 //  \par Discussion
 //    This program calculates the values of the constants \f$\{v_n\}\f$ as
 //    described in the reference (BFN) section 3.3.\n
 //    Polynomials stored as arrays have the coefficient of degree N
 //    in POLY(N).\n
 //    A polynomial which is identically 0 is given degree -1.
 //
 //  \par Licensing
 //		This code is distributed under the GNU LGPL license.
 //
 //	 \par Modified
 //		29 October 2019
 //
 //  \par Author
 //		Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.\n
 //    	C++ version by:\n
 //			John Burkardt (Florida State University, USA, 2003),\n
 //			Yekaterina Listyukhina (Russian Technological University, KMBO-03-16, Russia, 2019),\n
 //			Alexander Smekhov (Russian Technological University, KMBO-03-16, Russia, 2019),\n
 //			Alexey Burimov (Russian Technological University, KMBO-03-16, Russia, 2019).
 //
 //  \par Reference
 //    Paul Bratley, Bennett Fox, Harald Niederreiter,
 //    Algorithm 738:
 //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
 //    ACM Transactions on Mathematical Software,
 //    Volume 20, Number 4, pages 494-495, 1994.
 //
 \callgraph
 */
{
	int poly_h[MAXDEG + 1];
	int poly_h_deg, poly_b_deg;
	int i, r;
	int K_j, e, m;
	int term;
	//
	e = poly_pi_deg;
	//
	//  The polynomial h(x) = p_i(x)**(j-1) = b(x) on arrival.
	//
	//  In section 3.3, the values of h_i are defined with a minus sign,
	//  but in GF(2): h_i = -h_i for any value.
	//
	poly_h_deg = *pt_poly_b_deg;
	
	for (i = 0; i <= poly_h_deg; ++i)
	{
		poly_h[i] = poly_b[i];
	}
	//
	//  Now choose a value of K_j as defined in section 3.3.
	//  We must have 0 <= K_j < e*j = m.
	//  The limit condition on K_j does not seem very relevant
	//  in this program.
	//	Let K_j = e*(j - 1).
	//
	K_j = poly_h_deg;
	//
	//  Multiply b(x) by p_i(x) so b(x) becomes p(x)**j.
	//  In section 2.3, the values of b_i are defined with a minus sign,
	//  but in GF(2): b_i = -b_i for any value.
	//
	poly_b_deg = *pt_poly_b_deg;
	multiply_poly2(poly_pi_deg, poly_pi, poly_b_deg, poly_b, &poly_b_deg, poly_b);
	*pt_poly_b_deg = poly_b_deg;
	m = *pt_poly_b_deg;
	//
	//  Choose values of {v_n} in accordance with the conditions in section 3.3.
	//
	for (r = 0; r < K_j; ++r)
	{
		v[r] = 0;
	}
	v[K_j] = 1;
	
	for (r = K_j + 1; r <= m - 1; ++r)
	{
		v[r] = 1;
	}
	//
	//  Calculate the remaining v_n's using the recursion of section 2.3,
	//  remembering that the b_i's have the opposite sign.
	//
	for (r = 0; r <= size_of_v - m; ++r)
	{
		term = 0;
		for (i = 0; i <= m - 1; ++i)
		{
			term ^= poly_b[i] & v[r + i];
		}
		v[r + m] = term;
	}
	
	return;
}
//****************************************************************************80

void generate_next_nied2_real(//! [in] the dimension of the sequence to be generated
					int dim,
					//! [in, out] the index of the element entry to compute.  On output, \b seed is typically reset by this routine to \b seed+1
					uint64_t *seed,
					//! [out] the next quasirandom vector
					Real next_elem[])

/*! \brief Writes into \b next_elem an element of the Niederreiter sequence base 2.
 //
 //  \par Licensing
 //		This code is distributed under the GNU LGPL license.
 //
 //	 \par Modified
 //		29 October 2019
 //
 //  \par Author
 //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.\n
 //    	C++ version by:\n
 //			John Burkardt (Florida State University, USA, 2003),\n
 //			Yekaterina Listyukhina (Russian Technological University, KMBO-03-16, Russia, 2019),\n
 //			Alexander Smekhov (Russian Technological University, KMBO-03-16, Russia, 2019),\n
 //			Alexey Burimov (Russian Technological University, KMBO-03-16, Russia, 2019).
 //
 //  \par Reference
 //    Harald Niederreiter,
 //    Low-discrepancy and low-dispersion sequences,
 //    Journal of Number Theory,
 //    Volume 30, 1988, pages 51-70.
 //
 //  \par Local Parameters
 //    \p int \b cj[DIM_MAX][NBITS], the packed values \f$С^{(i)}_r\f$ of
 //    Niederreiter's \f$c^{(i)}_{jr}\f$.\n
 //		\p int \b dim_save, the spatial dimension of the sequence
 //    as specified on an initialization call.\n
 //		\p int \b Q[DIM_MAX], the numerators of the next item in the
 //    series.  These are like Niederreiter's \f$x^{(i)}_n\f$ (page 54) except that
 //    n is implicit, and the Q are integers.  To obtain
 //    the values of \f$x^{(i)}_n\f$, multiply by RECIP.
 //
 //  \callgraph
 */
{
	static Next_int cj[DIM_MAX][NBITS];
	static int dim_save = 0;
	static Next_int Q[DIM_MAX];
	static uint64_t seed_save = 0;
	Next_int seed_gray_code;
	Next_int i;
	int rightmost_zero_bit_pos;
	//
	//  Initialization.
	//
	if ( dim_save < 1 || dim != dim_save || *seed <= 0 )
	{
		if ( dim <= 0 || DIM_MAX < dim )
		{
			cout << "\n";
			cout << "NIEDERREITER2 - Fatal error!\n";
			cout << "  Bad spatial dimension.\n";
			exit ( 1 );
		}
		
		dim_save = dim;
		
		if ( *seed < 0 )
		{
			*seed = 0;
		}
		
		seed_save = *seed;
		//
		//  Calculate the C array.
		//
		calculate_c(dim_save, cj);
	}
	//
	//  Set up Q appropriately, depending on the Gray code of SEED.
	//
	//  You can do this every time, starting Q back at 0,
	//  or you can do it once, and then carry the value of Q
	//  around from the previous computation.
	//
	if ( *seed != seed_save + 1 )
	{
		seed_gray_code = *seed ^ (*seed >> 1);
		
		for (i = 0; i < dim_save; ++i)
		{
			Q[i] = 0;
		}
		
		rightmost_zero_bit_pos = 0;
		
		while ( seed_gray_code != 0 )
		{
			if ( (seed_gray_code & 1) != 0 )
			{
				for (i = 0; i < dim_save; ++i)
				{
					Q[i] ^= cj[i][rightmost_zero_bit_pos];
				}
			}
			seed_gray_code >>= 1;
			++rightmost_zero_bit_pos;
		}
	}
	//
	//  Multiply the numerators in Q by RECIP to get the next
	//  quasi-random vector.
	//
	for (i = 0; i < dim_save; ++i)
	{
		next_elem[i] = ((Real)Q[i]) * RECIP;
	}
	//
	//  Find the position of the right-hand zero in seed.  This
	//  is the bit that changes in the Gray-code representation as
	//  we go from seed to seed+1.
	//
	rightmost_zero_bit_pos = 0;
	i = *seed;
	
	while ( (i & 1) != 0 )
	{
		++rightmost_zero_bit_pos;
		i >>= 1;
	}
	//
	//  Check that we have not passed 2**NBITS calls.
	//
	if ( NBITS <= rightmost_zero_bit_pos )
	{
		cout << "\n";
		cout << "NIEDERREITER2 - Fatal error!\n";
		cout << "  Too many calls!\n";
		exit ( 1 );
	}
	//
	//  Compute the new numerators in vector Q.
	//
	for (i = 0; i < dim_save; ++i)
	{
		Q[i] ^= cj[i][rightmost_zero_bit_pos];
	}
	
	seed_save = *seed;
	++(*seed);
	
	return;
}
//****************************************************************************80

void generate_next_nied2_int(//! [in] the dimension of the sequence to be generated
					int dim,
					//! [in, out] the index of the element entry to compute. On output, \b seed is typically reset by this routine to \b seed+1
					uint64_t *seed,
					//! [out] the next quasirandom vector multiplied by (1 << NBITS)
					Next_int next_elem[])

/*! \brief Writes into \b next_elem nominator of an element of the Niederreiter sequence base 2 multiplied by (1 << NBITS).
 //
 //  \par Licensing
 //		This code is distributed under the GNU LGPL license.
 //
 //	 \par Modified
 //		29 October 2019
 //
 //  \par Author
 //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.\n
 //    	C++ version by:\n
 //			John Burkardt (Florida State University, USA, 2003),\n
 //			Yekaterina Listyukhina (Russian Technological University, KMBO-03-16, Russia, 2019),\n
 //			Alexander Smekhov (Russian Technological University, KMBO-03-16, Russia, 2019),\n
 //			Alexey Burimov (Russian Technological University, KMBO-03-16, Russia, 2019).
 //
 //  \par Reference
 //    Harald Niederreiter,
 //    Low-discrepancy and low-dispersion sequences,
 //    Journal of Number Theory,
 //    Volume 30, 1988, pages 51-70.
 //
 //  \par Local Parameters
 //    \p int \b cj[DIM_MAX][NBITS], the packed values \f$С^{(i)}_r\f$ of
 //    Niederreiter's \f$c^{(i)}_{jr}\f$.\n
 //		\p int \b dim_save, the spatial dimension of the sequence
 //    as specified on an initialization call.\n
 //		\p int \b Q[DIM_MAX], the numerators of the next item in the
 //    series.  These are like Niederreiter's \f$x^{(i)}_n\f$ (page 54) except that
 //    n is implicit, and the Q are integers.  To obtain
 //    the values of \f$x^{(i)}_n\f$, multiply by RECIP.
 //
 //  \callgraph
 */
{
	static Next_int cj[DIM_MAX][NBITS];
	static int dim_save = 0;
	static Next_int Q[DIM_MAX];
	static uint64_t seed_save = 0;
	Next_int seed_gray_code;
	Next_int i;
	int rightmost_zero_bit_no;
	//
	//  Initialization.
	//
	if ( dim_save < 1 || dim != dim_save || *seed <= 0 )
	{
		if ( dim <= 0 || DIM_MAX < dim )
		{
			cout << "\n";
			cout << "NIEDERREITER2 - Fatal error!\n";
			cout << "  Bad spatial dimension.\n";
			exit ( 1 );
		}
		
		dim_save = dim;
		
		if ( *seed < 0 )
		{
			*seed = 0;
		}
		
		seed_save = *seed;
		//
		//  Calculate the C array.
		//
		calculate_c(dim_save, cj);
	}
	//
	//  Set up Q appropriately, depending on the Gray code of SEED.
	//
	//  You can do this every time, starting Q back at 0,
	//  or you can do it once, and then carry the value of Q
	//  around from the previous computation.
	//
	if ( *seed != seed_save + 1 )
	{
		seed_gray_code = *seed ^ (*seed >> 1);
		
		for (i = 0; i < dim_save; ++i)
		{
			Q[i] = 0;
		}
		
		rightmost_zero_bit_no = 0;
		
		while ( seed_gray_code != 0 )
		{
			if ( (seed_gray_code & 1) != 0 )
			{
				for (i = 0; i < dim_save; ++i)
				{
					Q[i] ^= cj[i][rightmost_zero_bit_no];
				}
			}
			seed_gray_code >>= 1;
			++rightmost_zero_bit_no;
		}
	}
	//
	//  Multiply the numerators in Q by RECIP to get the next
	//  quasi-random vector.
	//
	for (i = 0; i < dim_save; ++i)
	{
		next_elem[i] = Q[i];
	}
	//
	//  Find the position of the right-hand zero in seed.  This
	//  is the bit that changes in the Gray-code representation as
	//  we go from seed to seed+1.
	//
	rightmost_zero_bit_no = 0;
	i = *seed;
	
	while ( (i & 1) != 0 )
	{
		++rightmost_zero_bit_no;
		i >>= 1;
	}
	//
	//  Check that we have not passed 2**NBITS calls.
	//
	if ( NBITS <= rightmost_zero_bit_no )
	{
		cout << "\n";
		cout << "NIEDERREITER2 - Fatal error!\n";
		cout << "  Too many calls!\n";
		exit ( 1 );
	}
	//
	//  Compute the new numerators in vector Q.
	//
	for (i = 0; i < dim_save; ++i)
	{
		Q[i] ^= cj[i][rightmost_zero_bit_no];
	}
	
	seed_save = *seed;
	++(*seed);
	
	return;
}
//****************************************************************************80


//Real *niederreiter2_generate (//! [in] the spatial dimension
//								int dim,
//								//! [in] the number of points desired
//								uint64_t amount,
//								//! [in. out] a seed for the random number generator
//								uint64_t *seed )
//
///*! \brief Generates a set of Niederreiter values.
// //	 \callgraph
// //
// //  \par Licensing:
// //		This code is distributed under the GNU LGPL license.
// //
// //	\par Modified
// //		29 October 2019
// //
// //  \par Author
// //		John Burkardt (Florida State University, USA, 2003),\n
// //		Alexander Smekhov (Russian Technological University, KMBO-03-16, Russia, 2019),\n
// //
// //	\return
// //   \p double \b sequence[dim_num*n], the points.
// */
//{
//	uint64_t j;
//	Real *sequence;
//	
//	sequence = new Real[dim*amount];
//	
//	for ( j = 0; j < amount; j++ )
//	{
//		generate_next_nied2_real(dim, seed, sequence+j*dim);
//	}
//	
//	return sequence;
//}
////****************************************************************************80

void multiply_poly2(
				//! [in] the degree of \f$p_a(x)\f$
				int poly_pa_deg,
				//! [in] the first polynomial \f$p_a(x) = \sum_{i=0}^\textbf{pa_deg} a_i\cdot x^i\f$
				int poly_pa[MAXDEG + 1],
				//! [in] the degree of \f$p_b(x)\f$
				int poly_pb_deg,
				//! [in] the second polynomial \f$p_b(x) = \sum_{j=0}^\textbf{pb_deg} b_j\cdot x^j\f$
				int poly_pb[MAXDEG + 1],
				//! [out] the pointer to a degree of \f$p_c(x)\f$
				int *pt_poly_pc_deg,
				//! [out] the factor polynomial \f$p_c(x) = \sum_{k=i+j}^{\textbf{pa_deg}+\textbf{pb_deg}} a_i b_j\cdot x^k\f$
				int poly_pc[MAXDEG + 1]
				)

/*! \brief
 //		Multiplies two polynomials in GF(2).
 //
 //	\par Discussion
 //		Function performs multiplication of polynomials \f$p_c(x) = p_a(x) \cdot p_b(x)\f$.\n
 //		Polynomials are stored as arrays of coefficients and have
 //		the coefficient of degree N as the N-th element of an array.\n
 //		A polynomial which is identically 0 is given degree -1.
 //
 //	\copyright
 //		This code is distributed under the GNU LGPL license. 
 //
 //	\par Modified
 //		29 October 2019
 //
 //	\author
 //		Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
 //    	C++ version by:\n
 //			John Burkardt (Florida State University, USA, 2003),\n
 //			Yekaterina Listyukhina (Russian Technological University, KMBO-03-16, Russia, 2019),\n
 //			Alexander Smekhov (Russian Technological University, KMBO-03-16, Russia, 2019),\n
 //			Alexey Burimov (Russian Technological University, KMBO-03-16, Russia, 2019).
 //
 */
{
	int i, j;
	int jlo, jhi;
	int poly_pc_clone[MAXDEG + 1];
	int term;

	*pt_poly_pc_deg = ( poly_pa_deg != -1 && poly_pb_deg != -1 ) ?  poly_pa_deg + poly_pb_deg  :  -1;

	if ( MAXDEG < *pt_poly_pc_deg )
	{
		cout << "\n";
		cout << "PLYMUL2 - Fatal error!\n";
		cout << "	Degree of the product exceeds MAXDEG.\n";
		exit ( 1 );
	}

	for (i = 0; i <= *pt_poly_pc_deg; ++i)
	{
		jlo = ( i - poly_pa_deg < 0 ) ?  0  :  i - poly_pa_deg;
		jhi = 	  ( i < poly_pb_deg ) ?  i  :  poly_pb_deg;
		//
		// Find pc_i as a mod2 sum over all pa_j*pb_k such that j + k = i
		//
		term = 0;
		for (j = jlo; j <= jhi; ++j)
		{
			term ^= poly_pa[i - j] & poly_pb[j];
		}
		poly_pc_clone[i] = term;
	}

	for (i = 0; i <= *pt_poly_pc_deg; ++i)
	{
		poly_pc[i] = poly_pc_clone[i];
	}
	
	for (i = *pt_poly_pc_deg + 1; i <= MAXDEG; ++i)
	{
		poly_pc[i] = 0;
	}

	return;
}


void timestamp (void)

/*!
 // \brief Prints the current YMDHMS date as a time stamp.
 //
 //	\par Example
 //		May 31 2001 09:45:54 AM
 //
 //	\copyright
 //		This code is distributed under the GNU LGPL license. 
 //
 //	\par Modified
 //		03 October 2003
 //
 //	\author
 //		John Burkardt
 //
 */
{
#	define TIME_SIZE 40
	
	static char time_buffer[TIME_SIZE];
	const struct tm *tm;
	size_t len;
	time_t now;
	
	now = time ( NULL );
	tm = localtime ( &now );
	
	len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );
	
	cout << time_buffer << "\n";
	
	return;
#	undef TIME_SIZE
}
