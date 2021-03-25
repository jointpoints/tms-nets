#include "../../../include/tms-nets/thirdparty/irrpoly/gfcheck.hpp"





// namespace <detail>





[[nodiscard]]
auto irrpoly::detail::derivative(const gfpoly &poly) -> gfpoly {
    if (poly.is_zero() || poly.degree() == 0) {
        return gfpoly(poly.field());
    }
    std::vector<uintmax_t> res(poly.size() - 1, 0);
    for (uintmax_t i = 1; i < poly.size(); ++i) {
        res[i - 1] = (i % poly.base()) * poly[i];
    }
    return gfpoly(poly.field(), res);
}



[[nodiscard]]
auto irrpoly::detail::x_pow_mod(uintmax_t pow, const gfpoly &mod) -> gfpoly {
    const auto n = mod.degree();
    gfpoly xn = gfpoly(mod.field(), 1) << n; // x^n
    gfpoly res(mod.field(), 1);

    uintmax_t d = 0;
    uintmax_t tmp = 0;
    for (auto m = res.degree(); pow + m >= n; m = res.degree()) {
        tmp = n - m, pow -= tmp, res <<= tmp;
        d = (res == xn) ? (d) ? (d -= pow, pow %= d, 0) : pow : d;
        res %= mod;
    }
    return res << pow;
}





// Check functions





[[nodiscard]]
auto irrpoly::gcd(gfpoly m, gfpoly n) -> gfpoly {
    if (m.is_zero() || n.is_zero()) {
        throw std::domain_error("arguments must be strictly positive");
    }
    if (m.degree() < n.degree()) {
        std::swap(m, n);
    }
    gfpoly u0 = m, u1 = gfpoly(m.field(), 1), u2 = gfpoly(m.field()),
        v0 = n, v1 = gfpoly(m.field()), v2 = gfpoly(m.field(), 1),
        w0 = gfpoly(m.field()), w1 = gfpoly(m.field()), w2 = gfpoly(m.field()),
        q = gfpoly(m.field());
    while (v0) {
        q = u0 / v0;
        w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
        u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
    }
    return u0;
}



[[nodiscard]]
auto irrpoly::is_irreducible_berlekamp(const gfpoly &poly) -> bool {
    if (poly.is_zero()) {
        return false;
    }
    const auto n = poly.degree();

    if (n == 0 || (poly[0] == 0 && n > 1)) {
        return false;
    }
    if (n == 1) {
        return true;
    }

    // builds matrix B - I and calculates it's rank
    auto berlekampMatrixRank = [](const gfpoly &val) {
        uintmax_t i = 0, j = 0, k = 0, l = 0;
        const auto n = val.degree();
        std::vector<std::vector<gfn>> B(n, std::vector<gfn>(n, gfn(val.field()))); // B = 0
        for (i = 0; i < n; ++i) {
            // B[i,*] = x ^ ip (mod val)
            const auto poly = detail::x_pow_mod(i * val.base(), val);
            for (j = 0, k = poly.degree(); j <= k; ++j) {
                B[i][j] += poly[j];
            }
            B[i][i] -= 1; // B - I
        }

        // reduces matrix to stepwise form
        bool f = false;
        gfn num(val.field());
        for (i = k = 0; i < n && k < n; ++k) {
            f = !!B[i][k];
            for (j = i + 1; j < n; ++j) {
                if (B[j][k]) {
                    if (f) {
                        num = B[j][k] / B[i][k];
                        B[j][k] = 0;
                        for (l = k + 1; l < n; ++l) {
                            B[j][l] -= B[i][l] * num;
                        }
                    } else {
                        for (l = k; l < n; ++l) {
                            std::swap(B[i][l], B[j][l]);
                        }
                        f = true;
                    }
                }
            }
            i += f;
        }
        return i;
    };

    // algorithm begins here
    auto d = detail::derivative(poly);
    return !!d && gcd(poly, d).degree() == 0 &&
        berlekampMatrixRank(poly) == poly.degree() - 1;
}



[[nodiscard]]
auto irrpoly::is_irreducible_rabin(const gfpoly &poly) -> bool {
    if (poly.is_zero()) {
        return false;
    }
    const auto n = poly.degree();

    if (n == 0 || (poly[0] == 0 && n > 1)) {
        return false;
    }
    if (n == 1) {
        return true;
    }

    // returns list of distinct prime divisors of n
    auto factorize = [](uintmax_t n) {
        std::vector<uintmax_t> list;
        const auto begin = n;
        for (uintmax_t d = 2; d * d <= n; ++d) {
            if (n % d) {
                continue;
            }
            list.emplace_back(begin / d);
            while (n % d == 0) {
                n /= d;
            }
        }
        if (n != 1) {
            list.emplace_back(begin / n);
        }
        return list;
    };

    auto P = poly.base();
    auto list = factorize(n);
    gfpoly tmp(poly.field()), x = gfpoly(poly.field(), {0, 1});
    for (auto i: list) {
        tmp = detail::x_pow_mod(detail::integer_power(P, i), poly) - x;
        if (tmp.is_zero() || gcd(poly, tmp).degree() > 0) {
            return false;
        }
    }

    tmp = detail::x_pow_mod(detail::integer_power(P, n), poly) - x;
    return tmp.is_zero();
}



[[nodiscard]]
auto irrpoly::is_irreducible_benor(const gfpoly &poly) -> bool {
    if (poly.is_zero()) {
        return false;
    }
    const auto n = poly.degree();

    if (n == 0 || (poly[0] == 0 && n > 1)) {
        return false;
    }
    if (n == 1) {
        return true;
    }

    auto P = poly.base();
    gfpoly tmp(poly.field()), x = gfpoly(poly.field(), {0, 1});
    for (uintmax_t m = n / 2, i = 1; i <= m; ++i) {
        tmp = detail::x_pow_mod(detail::integer_power(P, i), poly) - x;
        if (tmp.is_zero() || gcd(poly, tmp).degree() > 0) {
            return false;
        }
    }
    return true;
}



[[nodiscard]]
auto irrpoly::is_primitive_definition(const gfpoly &poly) -> bool {
    if (poly.is_zero()) {
        return false;
    }
    const auto n = poly.degree();

    if (n == 0 || (poly[0] == 0 && n > 1)) {
        return false;
    }
    if (n == 1 && poly[0] == 0) {
        return true;
    } // val = k * x + 0

    // this algorithm is defined only for normalized polynomials
    const auto npoly = poly / poly[n];

    // degenerate case
    auto P = npoly.base();
    if (P == 2 && npoly == gfpoly(npoly.field(), {1, 1})) {
        return false;
    }

    auto mp = gfn(npoly.field(), npoly[0]);
    mp = (n % 2) ? -mp : mp;

    // returns list of distinct prime divisors of n except 1 and n if it's prime
    auto factorize = [](uintmax_t n) {
        std::vector<uintmax_t> list;
        const auto begin = n;
        for (uintmax_t d = 2; d * d <= n; ++d) {
            if (n % d) {
                continue;
            }
            list.emplace_back(d);
            while (n % d == 0) {
                n /= d;
            }
        }
        if (n != 1 && n != begin) {
            list.emplace_back(n);
        }
        return list;
    };

    if (P > 2) {
        const auto p = P - 1;
        auto list = (p == 2) ? std::vector<uintmax_t>{2} : factorize(p);
        auto m = list.size() - 1;
        auto tmp = mp;
        for (uint32_t i = 1; i <= p; ++i, tmp *= mp) {
            if (i != p / list[m]) {
                continue;
            }
            if (tmp == 1) {
                return false;
            }
            if (m == 0) {
                break;
            }
            m -= 1;
        }
    }

    uintmax_t r = (detail::integer_power(P, n) - 1) / (P - 1);
    auto tmp = detail::x_pow_mod(r, poly) - mp;
    if (tmp) {
        return false;
    }

    auto list3 = factorize(r);
    const auto m = list3.size();
    for (size_t i = 0; i < m; ++i) {
        tmp = detail::x_pow_mod(r / list3[i], npoly);
        if (tmp.is_zero() || tmp.degree() == 0) {
            return false;
        }
    }

    return true;
}





// namespace <multithread>





[[nodiscard]]
auto irrpoly::multithread::make_check_func(
    irreducible_method irr_meth, primitive_method prim_meth)
-> typename pipeline<gfpoly, check_result>::payload_fn {
    return [=](const gfpoly &poly, std::optional<check_result> &res) {
        auto result = check_result{true, true};

        switch (irr_meth) {
        case irreducible_method::recommended:
            result.irreducible = is_irreducible(poly);
            break;
        case irreducible_method::berlekamp:
            result.irreducible = is_irreducible_berlekamp(poly);
            break;
        case irreducible_method::rabin:
            result.irreducible = is_irreducible_rabin(poly);
            break;
        case irreducible_method::benor:
            result.irreducible = is_irreducible_benor(poly);
            break;
        default:; // irreducible_method::nil
        }

        switch (prim_meth) {
        case primitive_method::recommended:
            result.primitive = result.irreducible ?
                               is_primitive(poly) : false;
            break;
        case primitive_method::definition:
            result.primitive = result.irreducible ?
                               is_primitive_definition(poly) : false;
            break;
        default:; // primitive_method::nil
        }

        res.emplace(result);
    };
}
