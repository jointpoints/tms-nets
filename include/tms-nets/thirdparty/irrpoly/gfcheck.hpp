/**
 * @file    polynomialgf.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#pragma once

#include "pipeline.hpp"
#include "gfpoly.hpp"

namespace irrpoly {

/**
 * Calculates greatest common divisor for two polynomials.
 * Originally taken from Boost library, then made some changes.
 */
[[nodiscard]]
auto gcd(gfpoly m, gfpoly n) -> gfpoly;

namespace detail {

/**
 * Quick computation of t^n for T and N being integer types.
 * Originally taken from Boost library, then made some changes.
 */
template<class T, class N>
[[nodiscard]]
auto integer_power(T t, N n) -> T {
    switch (n) {
    case 0:return static_cast<T>(1U);
    case 1:return t;
    case 2:return t * t;
    case 3:return t * t * t;
    default:T result = integer_power(t, n / 2);
        return (n & 1U) ? result * result * t : result * result;
    }
}

/**
 * Calculates derivative for given polynomial.
 */
[[nodiscard]]
auto derivative(const gfpoly &poly) -> gfpoly;

/**
 * Calculates (x^pow) % mod.
 * Reduces the amount of memory used during computation.
 * Work principle is that during column division x^n is represented as [0 ... 0 1],
 * so we can take only some of zeroes, perform division partially, then pad more
 * zeroes and continue division. As a bonus there are some cases when we'll
 * get 0 ... 0 1 (which is x^m) somewhere during division process. In this case we
 * could calculate tmp = (n - m) and skip (k * tmp) steps because after each
 * tmp division iterations we well get x^(m - k*tmp).
 * TODO: check FLINT sources for similar function and rewrite this to speed up.
 */
[[nodiscard]]
auto x_pow_mod(uintmax_t pow, const gfpoly &mod) -> gfpoly;

} // namespace detail

/**
 * This function implements Berlekamp's irreducibility test for polynomials over GF[P].
 * Before all computations common cases are checked: if deg(poly) = 1 then poly
 * is irreducible. If zero-indexed term is zero and poly is not zero then poly has
 * factor x, and so is reducible.
 * First step is computing the derivative poly', if it is zero - then polynomial
 * is a power of some other polynomial, and so is reducible.
 * Second step is calculating gcd(poly, poly'). If it is non-constant - then
 * poly has some factors common with poly', and so is reducible.
 * Third step is building Berlekamp's matrix B(m,m) and calculating it's rank.
 * Its rows consists of coefficients of polynomials x^(iP) (mod poly), 0 < i < n.
 * For more information read article "A Formalization of Berlekamp’s Factorization
 * Algorithm" by Davison, Joosten, Thiemann and Yamada.
 * Then B - I is calculated and rank(B - I) is found. If rank(B - I) == deg(poly) - 1
 * then poly is irreducible, otherwise it is reducible. For rank calculation
 * matrix is reduced to a stepwise form and number of steps is calculated,
 * this number is equal to matrix rank.
 * TODO: change the way matrix rank is calculated to make algorithm faster.
 */
[[nodiscard]]
auto is_irreducible_berlekamp(const gfpoly &poly) -> bool;

/**
 * This function implements Rabin's irreducibility test for polynomials over Galois field.
 * Alghoritm is fully described in article "Analysis of Rabin's irreducibility
 * test for polynomials over finite Fields" by Panario, Pittel, Richmond and Viola.
 * Added common case checks as in Berlekamp's test above.
 */
[[nodiscard]]
auto is_irreducible_rabin(const gfpoly &poly) -> bool;

/**
 * This function implements Ben-Or's irreducibility test for polynomials over Galois field.
 * Alghoritm's pseudocode is provided in article "Tests and constructions of
 * irreducible polynomials over finite fields" by Gao and Panario.
 * Added common case checks as in Berlekamp's test above.
 */
[[nodiscard]]
auto is_irreducible_benor(const gfpoly &poly) -> bool;

/**
 * This function performs quickest irreducibility test, defined by benchmark results.
 * TODO: implement Distinct Degree Factorization algorithm, find out cases then it is fastest.
 */
[[nodiscard]]
inline
auto is_irreducible(const gfpoly &poly) -> bool {
    switch (poly.base()) {
    case 2: return is_irreducible_berlekamp(poly);
    default: return is_irreducible_benor(poly);
    }
}

/**
 * This function implements primitivity test for polynomials over Galois field.
 * Alghoritm is fully described in article "Primitive polynomials over finite
 * fields" by Hansen and Mullen. Added common case checks as in Berlekamp's test above.
 */
[[nodiscard]]
auto is_primitive_definition(const gfpoly &poly) -> bool;

/**
 * This function performs quickest primitivity test, defined by benchmark results.
 */
[[nodiscard]]
inline
auto is_primitive(const gfpoly &val) -> bool {
    return is_irreducible(val) ? is_primitive_definition(val) : false;
}

namespace multithread {

struct check_result {
    bool irreducible;
    bool primitive;
};

/**
 * Irreducibility tests available.
 */
enum class irreducible_method {
    nil, ///< do not test
    berlekamp, ///< Berlekam's test
    rabin, ///< Rabin's test
    benor, ///< Ben-Or's test
    recommended, ///< fastest test
};

/**
 * Primitivity tests available.
 */
enum class primitive_method {
    nil, ///< не проверять
    definition, ///< проверка по определению
    recommended, ///< оптимальный алгоритм
};

using polychecker = pipeline<gfpoly, check_result>;

/**
 * Creates payload_fn for irrpoly::multithread::pipeline.
 * In case nil method is selected - the result is true.
 */
[[nodiscard]]
auto make_check_func(
    irreducible_method irr_meth, primitive_method prim_meth)
-> typename pipeline<gfpoly, check_result>::payload_fn;

} // namespace multithread

} // namespace irrpoly
