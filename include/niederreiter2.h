//  NIEDERREITER2.H
//
//  Local Parameters:
//
//    int NBITS, the number of bits (not counting the sign) in a
//    fixed-point integer.
//
//    int MAXDEG, the highest degree of polynomial
//    to be handled. Must be greater or equal than NBITS.
//
//    int DIM_MAX, the maximum dimension that will be used.
//
//    double RECIP, the multiplier which changes the
//    integers in NEXTQ into the required real values in QUASI.
//
# define NBITS 63U
# define MAXDEG 67U
# define DIM_MAX 20U

//
//  Thanks to Bradley Keister for pointing out that the compiler might
//  have trouble using POW since it's not declared yet...
//
//static double RECIP = pow ( 2.0, -NBITS );

typedef int64_t Next_int;

typedef long double Real;

//! \todo Remove all global constants
static Real const RECIP = 1.0 / (Real)(1ULL << NBITS);

void calculate_c(int dim, Next_int cj[DIM_MAX][NBITS]);

void calculate_v(int size_of_v, int poly_p_deg, int poly_p[MAXDEG+1], int *pt_poly_b_deg, int poly_b[MAXDEG+1], int v[]);

void generate_next_nied2_real(int dim, uint64_t *seed, Real next_elem[]);

void generate_next_nied2_int(int dim, uint64_t *seed, Next_int next_elem[]);

Real *niederreiter2_generate(int dim, uint64_t *seed, uint64_t amount);

void multiply_poly2(int poly_pa_deg, int poly_pa[MAXDEG+1], int poly_pb_deg, int poly_pb[MAXDEG+1], int *pt_poly_pc_deg, int poly_pc[MAXDEG+1]);

void timestamp (void);
