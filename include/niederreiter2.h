//  NIEDERREITER2.H
//
//  Local Parameters:
//
//    uint64_t MAXDEG, the highest degree of polynomial
//    to be handled. 
//
//    uint64_t DIM_MAX, the maximum dimension that will be used.
//
//    uint64_t NBITS, the number of bits (not counting the sign) in a
//    fixed-point integer.
//
//    double RECIP, the multiplier which changes the
//    integers in NEXTQ into the required real values in QUASI.
//
#define MAXDEG 67U
#define DIM_MAX 20U
#define NBITS 63U

static double RECIP = 1.0 / (double)(1ULL << NBITS);

void calcc2(uint64_t dimen, uint64_t cj[DIM_MAX][NBITS]);

void calcv2(uint64_t maxv, uint64_t px_deg, uint64_t px[MAXDEG + 1], uint64_t* b_deg, uint64_t b[MAXDEG + 1], uint64_t v[]);

void niederreiter2(uint64_t dim, uint64_t* seed, double quasi[]);

double* niederreiter2_generate(uint64_t dim_num, uint64_t n, uint64_t* seed);

void plymul2(uint64_t pa_deg, uint64_t pa[MAXDEG + 1], uint64_t pb_deg, uint64_t pb[MAXDEG + 1], uint64_t* pc_deg, uint64_t pc[MAXDEG + 1]);

void timestamp();
