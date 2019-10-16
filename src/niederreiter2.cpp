# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "niederreiter2.hpp"


void calcc2 (//! [in] the dimension of the sequence to be generated.
			 int dim_num,
			 //! [out] the packed values of Niederreiter's C(I,J,R)
			 int cj[DIM_MAX][NBITS])
/*! \brief CALCC2 computes values of the constants \f$c^{(i)}_{jr}\f$.
 //
 //	\callgraph
 //
 //	\par Discussion
 //		This program calculates the values of the constants C(I,J,R).
 //		As far as possible, Niederreiter's notation is used.\n
 //		For each value of I, we first calculate all the corresponding
 //		values of C.  These are held in the array CI.  All these
 //		values are either 0 or 1.\n
 //		Next we pack the values into the
 //    	array CJ, in such a way that CJ(I,R) holds the values of C
 //    	for the indicated values of I and R and for every value of
 //    	J from 1 to NBITS.  The most significant bit of CJ(I,R)
 //    	(not counting the sign bit) is C(I,1,R) and the least
 //    	significant bit is C(I,NBITS,R).
 //
 // \par Local Parameters
 //    	\p int \b MAXE, the highest degree among \b DIM_MAX irreducible polynomials over GF(2).\n
 //	 	\p int \b MAXV, the maximum possible index used in V.
 //
 //	\copyright This code is distributed under the GNU LGPL license.
 //
 // \par Modified
 //		29 March 2003
 //
 // \author Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
 //    		 C++ version by John Burkardt.
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
	
	int add[2][2];
	int b[MAXDEG+1];
	int b_deg;
	int ci[NBITS][NBITS];
	int count;
	int e;
	int i;
	static int irred[DIM_MAX][MAXE+1] =
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
	int irred_deg[DIM_MAX] =
	{ 1, 1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6 };
	int j;
	int maxv = NBITS + MAXE;
	int mul[2][2];
	int nextq[DIM_MAX];
	int p;
	int px[MAXDEG+1];
	int px_deg;
	int q;
	int r;
	int sub[2][2];
	int term;
	int u;
	int v[NBITS+MAXE+1];
	//
	//  Prepare to work in Z2.
	//
	setfld2 ( add, mul, sub );
	
	for ( i = 0; i < dim_num; i++ )
	{
		//
		//  For each dimension, we need to calculate powers of an
		//  appropriate irreducible polynomial:  see Niederreiter
		//  page 65, just below equation (19).
		//
		//  Copy the appropriate irreducible polynomial into PX,
		//  and its degree into E.  Set polynomial B = PX ** 0 = 1.
		//  M is the degree of B.  Subsequently B will hold higher
		//  powers of PX.
		//
		e = irred_deg[i];
		
		px_deg = irred_deg[i];
		
		for ( j = 0; j <= px_deg; j++ )
		{
			px[j] = irred[i][j];
		}
		
		b_deg = 0;
		b[0] = 1;
		//
		//  Niederreiter (page 56, after equation (7), defines two
		//  variables Q and U.  We do not need Q explicitly, but we do need U.
		//
		u = 0;
		
		for ( j = 0; j < NBITS; j++ )
		{
			//
			//  If U = 0, we need to set B to the next power of PX
			//  and recalculate V.  This is done by subroutine CALCV.
			//
			if ( u == 0 )
			{
				calcv2 ( maxv, px_deg, px, add, mul, sub, &b_deg, b, v );
			}
			//
			//  Now C is obtained from V.  Niederreiter obtains A from V (page 65,
			//  near the bottom), and then gets C from A (page 56, equation (7)).
			//  However this can be done in one step.  Here CI(J,R) corresponds to
			//  Niederreiter's C(I,J,R).
			//
			for ( r = 0; r < NBITS; r++ )
			{
				ci[j][r] = v[r+u];
			}
			//
			//  Increment U.
			//
			//  If U = E, then U = 0 and in Niederreiter's
			//  paper Q = Q + 1.  Here, however, Q is not used explicitly.
			//
			u = u + 1;
			if ( u == e )
			{
				u = 0;
			}
			
		}
		//
		//  The array CI now holds the values of C(I,J,R) for this value
		//  of I.  We pack them into array CJ so that CJ(I,R) holds all
		//  the values of C(I,J,R) for J from 1 to NBITS.
		//
		for ( r = 0; r < NBITS; r++ )
		{
			term = 0;
			for ( j = 0; j < NBITS; j ++ )
			{
				term = 2 * term + ci[j][r];
			}
			cj[i][r] = term;
		}
		
	}
	
	return;
# undef MAXE
}
//****************************************************************************80

void calcv2 (//! [in] the dimension of the array V
			 int maxv,
			 //! [in] the degree of PX.
			 int px_deg,
			 //! [in] the appropriate irreducible polynomial for the dimension currently being considered.
			 int px[MAXDEG+1],
			 //! [in] add, mul, sub: the addition, multiplication and subtraction tables, mod 2.
			 int add[2][2],
			 int mul[2][2],
			 int sub[2][2],
			 //! [in, out] the degree of the polynomial B
			 int *b_deg,
			 /*! [in, out] on input, B is the polynomial
			  defined in section 2.3 of BFN. The degree of B implicitly defines
			  the parameter J of section 3.3, by degree(B) = E*(J-1).  On output,
			  B has been multiplied by PX, so its degree is now E * J*/
			 int b[MAXDEG+1],
			 //! [out] the computed V array
			 int v[] )

/*! \brief CALCV2 calculates the value of the constants V(J,R).
 //
 // 	\par Discussion
 //    This program calculates the values of the constants V(J,R) as
 //    described in the reference (BFN) section 3.3.  It is called from CALCC2.\n
 //    Polynomials stored as arrays have the coefficient of degree N
 //    in POLY(N).\n
 //    A polynomial which is identically 0 is given degree -1.
 //
 //  \par Licensing
 //		This code is distributed under the GNU LGPL license.
 //
 //  \par Modified
 //		29 March 2003
 //
 //  \par Author
 //		Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.\n
 //    	C++ version by John Burkardt.
 //
 //  \par Reference
 //    Paul Bratley, Bennett Fox, Harald Niederreiter,
 //    Algorithm 738:
 //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
 //    ACM Transactions on Mathematical Software,
 //    Volume 20, Number 4, pages 494-495, 1994.
 //
 //  \par Local Parameters
 //		\p int \b arbit, indicates where the user can place
 //    an arbitrary element of the field of order 2.  This means
 //    0 <= \b arbit < 2.\n
 //    \p int \b bigm, is the M used in section 3.3.
 //    It differs from the [little] m used in section 2.3,
 //    denoted here by M.\n
 //    \p int \b nonzer, shows where the user must put an arbitrary
 //    non-zero element of the field.  For the code, this means
 //    0 < \b nonzer < 2.
 \callgraph
 */
{
	static int arbit = 1;
	int bigm;
	int e;
	int h[MAXDEG+1];
	int h_deg;
	int i;
	int j;
	int kj;
	int m;
	static int nonzer = 1;
	static int p = 2;
	int pb_deg;
	static int q = 2;
	int r;
	int term;
	//
	e = px_deg;
	//
	//  The polynomial H is PX**(J-1), which is the value of B on arrival.
	//
	//  In section 3.3, the values of Hi are defined with a minus sign:
	//  don't forget this if you use them later!
	//
	h_deg = *b_deg;
	
	for ( i = 0; i <= h_deg; i++ )
	{
		h[i] = b[i];
	}
	
	bigm = h_deg;
	//
	//  Multiply B by PX so B becomes PX**J.
	//  In section 2.3, the values of Bi are defined with a minus sign:
	//  don't forget this if you use them later!
	//
	pb_deg = *b_deg;
	
	plymul2 ( add, mul, px_deg, px, pb_deg, b, &pb_deg, b );
	
	*b_deg = pb_deg;
	m = *b_deg;
	//
	//  We don't use J explicitly anywhere, but here it is just in case.
	//
	j = m / e;
	//
	//  Now choose a value of Kj as defined in section 3.3.
	//  We must have 0 <= Kj < E*J = M.
	//  The limit condition on Kj does not seem very relevant
	//  in this program.
	//
	kj = bigm;
	//
	//  Choose values of V in accordance with the conditions in section 3.3.
	//
	for ( r = 0; r < kj; r++ )
	{
		v[r] = 0;
	}
	v[kj] = 1;
	
	if ( kj < bigm )
	{
		term = sub [ 0 ] [ h[kj] ];
		
		for ( r = kj+1; r <= bigm-1; r++ )
		{
			v[r] = arbit;
			//
			//  Check the condition of section 3.3,
			//  remembering that the H's have the opposite sign.
			//
			term = sub [ term ] [ mul [ h[r] ] [ v[r] ] ];
			
		}
		//
		//  Now V(BIGM) is anything but TERM.
		//
		v[bigm] = add [ nonzer] [ term ];
		
		for ( r = bigm+1; r <= m-1; r++ )
		{
			v[r] = arbit;
		}
	}
	else
	{
		for ( r = kj+1; r <= m-1; r++ )
		{
			v[r] = arbit;
		}
		
	}
	//
	//  Calculate the remaining V's using the recursion of section 2.3,
	//  remembering that the B's have the opposite sign.
	//
	for ( r = 0; r <= maxv - m; r++ )
	{
		term = 0;
		for ( i = 0; i <= m-1; i++ )
		{
			term = sub [ term] [ mul [ b[i] ] [ v[r+i] ] ];
		}
		v[r+m] = term;
	}
	
	return;
}
//****************************************************************************80

void niederreiter2 (//! [in] the dimension of the sequence to be generated.
					int dim_num,
					//! [in, out] the index of the element entry to compute.  On output, SEED is typically reset by this routine to SEED+1.
					int *seed,
					//! [out] the next quasirandom vector.
					double quasi[] )

/*! \brief NIEDERREITER2 returns an element of the Niederreiter sequence base 2.
 //
 //  \par Licensing
 //		This code is distributed under the GNU LGPL license.
 //
 //  \par Modified
 //    29 March 2003
 //
 //  \par Author
 //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.\n
 //    C++ version by John Burkardt.
 //
 //  \par Reference
 //    Harald Niederreiter,
 //    Low-discrepancy and low-dispersion sequences,
 //    Journal of Number Theory,
 //    Volume 30, 1988, pages 51-70.
 //
 //  \par Local Parameters
 //    \p int CJ(DIM_MAX,0:NBITS-1), the packed values of
 //    Niederreiter's C(I,J,R).\n
 //		\p int \b DIM_SAVE, the spatial dimension of the sequence
 //    as specified on an initialization call.\n
 //		\p int \b COUNT, the index of the current item in the sequence,
 //    expressed as an array of bits.  COUNT(R) is the same as Niederreiter's
 //    AR(N) (page 54) except that N is implicit.\n
 //		\p int \b NEXTQ[DIM_MAX], the numerators of the next item in the
 //    series.  These are like Niederreiter's XI(N) (page 54) except that
 //    N is implicit, and the NEXTQ are integers.  To obtain
 //    the values of XI(N), multiply by RECIP.
 //
 //  \callgraph
 */
{
	static int cj[DIM_MAX][NBITS];
	static int dim_save = 0;
	int gray;
	int i;
	static int nextq[DIM_MAX];
	int r;
	int skip;
	static int seed_save = 0;
	//
	//  Initialization.
	//
	if ( dim_save < 1 || dim_num != dim_save || *seed <= 0 )
	{
		if ( dim_num <= 0 || DIM_MAX < dim_num )
		{
			cout << "\n";
			cout << "NIEDERREITER2 - Fatal error!\n";
			cout << "  Bad spatial dimension.\n";
			exit ( 1 );
		}
		
		dim_save = dim_num;
		
		if ( *seed < 0 )
		{
			*seed = 0;
		}
		
		seed_save = *seed;
		//
		//  Calculate the C array.
		//
		calcc2 ( dim_save, cj );
	}
	//
	//  Set up NEXTQ appropriately, depending on the Gray code of SEED.
	//
	//  You can do this every time, starting NEXTQ back at 0,
	//  or you can do it once, and then carry the value of NEXTQ
	//  around from the previous computation.
	//
	if ( *seed != seed_save + 1 )
	{
		gray = ( *seed ) ^ ( *seed / 2 );
		
		for ( i = 0; i < dim_save; i++ )
		{
			nextq[i] = 0;
		}
		
		r = 0;
		
		while ( gray != 0 )
		{
			if ( ( gray % 2 ) != 0 )
			{
				for ( i = 0; i < dim_save; i++ )
				{
					nextq[i] = ( nextq[i] ) ^ ( cj[i][r] );
				}
			}
			gray = gray / 2;
			r = r + 1;
		}
	}
	//
	//  Multiply the numerators in NEXTQ by RECIP to get the next
	//  quasi-random vector.
	//
	for ( i = 0; i < dim_save; i++ )
	{
		quasi[i] = ( ( double ) nextq[i] ) * RECIP;
	}
	//
	//  Find the position of the right-hand zero in SEED.  This
	//  is the bit that changes in the Gray-code representation as
	//  we go from SEED to SEED+1.
	//
	r = 0;
	i = *seed;
	
	while ( ( i % 2 ) != 0 )
	{
		r = r + 1;
		i = i / 2;
	}
	//
	//  Check that we have not passed 2**NBITS calls.
	//
	if ( NBITS <= r )
	{
		cout << "\n";
		cout << "NIEDERREITER2 - Fatal error!\n";
		cout << "  Too many calls!\n";
		exit ( 1 );
	}
	//
	//  Compute the new numerators in vector NEXTQ.
	//
	for ( i = 0; i < dim_save; i++ )
	{
		nextq[i] = ( nextq[i] ) ^ ( cj[i][r] );
	}
	
	seed_save = *seed;
	*seed = *seed + 1;
	
	return;
}
//****************************************************************************80

double *niederreiter2_generate (//! [in] the spatial dimension.
								int dim_num,
								//! [in] the number of points desired.
								int n,
								//! [in. out] a seed for the random number generator.
								int *seed )

/*! \brief NIEDERREITER2_GENERATE generates a set of Niederreiter values.
 \callgraph
 //
 //  \par Licensing:
 //		This code is distributed under the GNU LGPL license.
 //
 //  \par Modified
 //		11 December 2009
 //
 //  \par Author
 //		John Burkardt
 //
 //	\return
 //   \p double \b r[dim_num*n], the points.
 */
{
	int j;
	double *r;
	
	r = new double[dim_num*n];
	
	for ( j = 0; j < n; j++ )
	{
		niederreiter2 ( dim_num, seed, r+j*dim_num );
	}
	
	return r;
}
//****************************************************************************80

void plymul2	(//! [in] the addition table for GF(2)
				 int add[2][2],
				 //! [in] the multiplication table for GF(2)
				 int mul[2][2],
				 //! [in] the degree of \f$p_a(x)\f$
				 int pa_deg,
				 //! [in] the first polynomial \f$p_a(x)\f$
				 int pa[MAXDEG+1],
				 //! [in] the degree of \f$p_b(x)\f$
				 int pb_deg,
				 //! [in] the second polynomial \f$p_b(x)\f$
				 int pb[MAXDEG+1],
				 //! [out] the degree of \f$p_c(x)\f$
				 int *pc_deg,
				 //! [out] the factor polynomial \f$p_c(x)\f$
				 int pc[MAXDEG+1]
				)

/*!
 //	\brief
 //		PLYMUL2 multiplies two polynomials in GF(2)
 //
 //	\callgraph
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
 //		29 March 2003
 //
 //	\author
 //		Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
 //		C++ version by John Burkardt.
 //
 */
{
	int i;
	int j;
	int jhi;
	int jlo;
	int pt[MAXDEG+1];
	int term;

	if ( pa_deg == -1 || pb_deg == -1 )
	{
		*pc_deg = -1;
	}
	else
	{
		*pc_deg = pa_deg + pb_deg;
	}

	if ( MAXDEG < *pc_deg )
	{
		cout << "\n";
		cout << "PLYMUL2 - Fatal error!\n";
		cout << "	Degree of the product exceeds MAXDEG.\n";
		exit ( 1 );
	}

	for ( i = 0; i <= *pc_deg; i++ )
	{

		jlo = i - pa_deg;
		if ( jlo < 0 ) 
		{
			jlo = 0;
		}

		jhi = pb_deg;
		if ( i < jhi ) 
		{
			jhi = i;
		}

		term = 0;

		for ( j = jlo; j <= jhi; j++ ) 
		{
			term = add [ term ] [ mul [ pa[i-j] ] [ pb[j] ] ];
		}
		pt[i] = term;
	}

	for ( i = 0; i <= *pc_deg; i++ )
	{
		pc[i] = pt[i];
	}

	for ( i = *pc_deg + 1; i <= MAXDEG; i++ )
	{
		pc[i] = 0;
	}

	return;
}


void r8mat_write	(//! [in] the output filename
					 string output_filename,
					 //! [in] the spatial dimension (amount of variables in target function)
					 int m,
					 //! [in] the number of points
					 int n,
					 //! [in] the m*n table data
					 double table[]
					)

/*!
 //	\brief
 //		R8MAT_WRITE writes an R8MAT file.
 //
 //	\callgraph
 //
 //	\par Discussion
 //		An R8MAT is an array of R8's.
 //
 //	\copyright
 //		This code is distributed under the GNU LGPL license. 
 //
 //	\par Modified
 //		29 June 2009
 //
 //	\author
 //		John Burkardt
 //
 */
{
	int i;
	int j;
	ofstream output;
//
//	Open the file.
//
	output.open ( output_filename.c_str ( ) );

	if ( !output )
	{
		cerr << "\n";
		cerr << "R8MAT_WRITE - Fatal error!\n";
		cerr << "	Could not open the output file.\n";
		return;
	}
//
//	Write the data.
//
	for ( j = 0; j < n; j++ )
	{
		for ( i = 0; i < m; i++ )
		{
			output << "	" << setw(24) << setprecision(16) << table[i+j*m];
		}
		output << "\n";
	}
//
//	Close the file.
//
	output.close ( );

	return;
}


void setfld2	(//! [in] the addition table for GF(2)
				 int add[2][2],
				 //! [in] the multiplication table for GF(2)
				 int mul[2][2],
				 //! [in] the subtraction table for GF(2)
				 int sub[2][2]
				)

/*! 
 // \brief
 //		SETFLD2 sets up arithmetic tables for GF(2).
 //
 //	\callgraph
 //
 //	\par Discussion
 //		SETFLD2 sets up addition, multiplication, and subtraction tables 
 //		for GF(2).
 //
 //	\copyright
 //		This code is distributed under the GNU LGPL license. 
 //
 //	\par Modified
 //		29 March 2003
 //
 //	\author
 //		Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
 //		C++ version by John Burkardt.
 //
 */
{
	int i;
	int j;
	int p = 2;
	int q = 2;
	
	for ( i = 0; i < q; i++ )
	{
		for ( j = 0; j < q; j++ )
		{
			add[i][j] = ( i + j ) % p;
			mul[i][j] = ( i * j ) % p;
		}
	}
//
//	Use the addition table to set the subtraction table.
//
	for ( i = 0; i < q; i++ )
	{
		for ( j = 0; j < q; j++ )
		{
			sub[ add[i][j] ] [i] = j;
		}
	}
	
	return;
}


void timestamp (void)

/*!
 //	\brief
 //		TIMESTAMP prints the current YMDHMS date as a time stamp.
 //
 //	\callgraph
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
