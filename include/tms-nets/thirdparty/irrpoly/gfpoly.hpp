/**
 * @file    gfpoly.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#pragma once

#include "gf.hpp"

#include <vector>
#include <functional>
#include <initializer_list>
#include <string>
#include <sstream>
#include <cctype>

namespace irrpoly {

/**
 * gfpoly represents a polynomial over Galois field.
 * This class is originally taken from Boost library but was significantly changed.
 * Get by index operation [i] returns polynomial term x^i.
 * Polynomial is either zero or reduced which means that leading term is non-zero.
 */
class gfpoly final {
private:
    gf m_field;
    std::vector<uintmax_t> m_data; ///< polynomial coefficients

public:
    /**
     * Generates random polynomial over provided Galois field of given degree.
     */
    static
    auto random(const gf &field, uintmax_t degree) -> gfpoly {
        static std::random_device rd;
#ifdef __LP64__
        static std::mt19937_64 gen(rd());
#else
        static std::mt19937 gen(rd());
#endif
        std::uniform_int_distribution<uint_fast64_t> dis(0, field->base() - 1);
        std::vector<uintmax_t> data;
        data.reserve(degree + 1);
        for (uintmax_t i = 0; i < degree; ++i) {
            data.push_back(dis(gen));
        }
        data.push_back(1);
        while (data[0] == 0) {
            data[0] = dis(gen);
        }
        return gfpoly(field, std::move(data));
    }

    [[nodiscard]]
    auto value() const -> const std::vector<uintmax_t>&;

private:
    /**
     * Removes leading zeroes.
     */
    auto reduce() -> gfpoly &;

public:
    explicit
    gfpoly(const gf &field);

    gfpoly(const gf &field, const std::vector<uintmax_t> &l);

    gfpoly(const gf &field, std::vector<uintmax_t> &&l);

    auto operator=(const std::vector<uintmax_t> &l) -> gfpoly &;

    gfpoly(const gf &field, std::initializer_list<uintmax_t> l);

    auto operator=(std::initializer_list<uintmax_t> l) -> gfpoly &;

    gfpoly(const gfpoly &p);

    gfpoly(gfpoly &&p);

    auto operator=(const gfpoly &p) -> gfpoly &;

    explicit
    gfpoly(const gfn &value);

    auto operator=(const gfn &value) -> gfpoly &;

    gfpoly(const gf &field, uintmax_t value);

    auto operator=(uintmax_t value) -> gfpoly &;

    [[nodiscard]]
    auto field() const -> const gf &;

    [[nodiscard]]
    auto base() const -> uintmax_t;

    [[nodiscard]]
    auto size() const -> uintmax_t;

    /**
     * Returns polynomial degree. For zero polynomial degree is undefined.
     */
    [[nodiscard]]
    auto degree() const -> uintmax_t;

    /**
     * Polynomial couldn't be mutated by index because it would be possible to
     * replace some term by gfn instance with field different to polynomial's one.
     */
    auto operator[](uintmax_t i) const -> uintmax_t;

    [[nodiscard]]
    auto is_zero() const -> bool;

    explicit operator bool() const;

    [[maybe_unused]]
    auto set_zero() -> gfpoly &;

private:
    /**
     * UNSAFE! Requires lb and rb to lay between 0 and base().
     * Returned number also lays between 0 and base().
     */
    [[nodiscard]]
    auto add(uintmax_t lb, uintmax_t rb) const -> uintmax_t;

    /**
     * UNSAFE! Requires lb and rb to lay between 0 and base().
     * Returned number also lays between 0 and base().
     */
    [[nodiscard]]
    auto sub(uintmax_t lb, uintmax_t rb) const -> uintmax_t;

    /**
     * UNSAFE! Requires lb and rb to lay between 0 and base().
     * Returned number also lays between 0 and base().
     */
    [[nodiscard]]
    auto mul(uintmax_t lb, uintmax_t rb) const -> uintmax_t;

    /**
     * UNSAFE! Requires lb and rb to lay between 0 and base().
     * Returned number also lays between 0 and base().
     */
    [[nodiscard]]
    auto div(uintmax_t lb, uintmax_t rb) const -> uintmax_t;

    /**
     * UNSAFE! Requires rb to lay between 0 and base().
     * Returned number also lays between 0 and base().
     */
    [[nodiscard]]
    auto neg(uintmax_t rb) const -> uintmax_t;

    using OP = uintmax_t (gfpoly::*)(uintmax_t, uintmax_t) const;

    auto transform(uintmax_t value, OP op) -> gfpoly &;

    auto transform(const gfn &value, OP op) -> gfpoly &;

    auto transform(const gfpoly &value, OP op) -> gfpoly &;

public:
    auto operator+=(uintmax_t value) -> gfpoly &;

    auto operator+=(const gfn &value) -> gfpoly &;

    auto operator+=(const gfpoly &value) -> gfpoly &;

    auto operator-=(uintmax_t value) -> gfpoly &;

    auto operator-=(const gfn &value) -> gfpoly &;

    auto operator-=(const gfpoly &value) -> gfpoly &;

    auto operator*=(uintmax_t value) -> gfpoly &;

    auto operator*=(const gfn &value) -> gfpoly &;

    auto operator/=(uintmax_t value) -> gfpoly &;

    auto operator/=(const gfn &value) -> gfpoly &;

    template<class U>
    auto operator%=(const U & /*value*/) -> gfpoly & {
        // we can always divide by a scalar, so there is no remainder
        return set_zero();
    }

private:
    auto multiply(const gfpoly &a, const gfpoly &b) -> gfpoly &;

public:
    auto operator*=(const gfpoly &value) -> gfpoly &;

private:
    /**
     * This is the core method of the whole class. It takes the most time during computations.
     * TODO: try implementing polynomial division using Discrete Fourier Transform.
     */
    static auto division(gfpoly u, const gfpoly &v) -> std::pair<gfpoly, gfpoly> {
        uintmax_t const m = u.size() - 1, n = v.size() - 1;
        uintmax_t k = m - n;
        gfpoly q(u.field());
        q.m_data.resize(m - n + 1, 0);

        auto division_impl = [](gfpoly *q, gfpoly *u, const gfpoly &v, auto n, auto k) {
            q->m_data[k] = q->div(u->m_data[n + k], v.m_data[n]);
            for (auto j = n + k; j > k;) {
                j--;
                u->m_data[j] = u->sub(u->m_data[j],
                                      q->mul(q->m_data[k], v.m_data[j - k]));
            }
        };

        do {
            division_impl(&q, &u, v, n, k);
        } while (k-- != 0);
        u.m_data.resize(n, 0);
        return std::make_pair(q.reduce(), u.reduce());
    }

    /**
     * Calculates a / b and a % b, returning the pair (quotient, remainder) together
     * because the same amount of computation yields both.
     * This function is not defined for division by zero: user beware.
     */
    static auto quotient_remainder(const gfpoly &dividend, const gfpoly &divisor)
    -> std::pair<gfpoly, gfpoly> {
        if (dividend.size() < divisor.size()) {
            return std::make_pair(gfpoly(dividend.field()), dividend);
        }
        return division(dividend, divisor);
    }

public:
    auto operator/=(const gfpoly &value) -> gfpoly &;

    auto operator%=(const gfpoly &value) -> gfpoly &;

    /**
     * Logically equal to operation this /= x^n. Defined only when such devision is possible.
     */
    template<typename U>
    auto operator>>=(U const &n) -> gfpoly & {
        if (n > degree() || !std::all_of(
            m_data.begin(), m_data.begin() + n,
            [](uintmax_t x) { return x == 0; })) {
            throw std::logic_error("division is impossible");
        }
        m_data.erase(m_data.begin(), m_data.begin() + n);
        return *this;
    }

    /**
     * Logically equal to operation this *= x^n.
     */
    template<typename U>
    auto operator<<=(U const &n) -> gfpoly & {
        reduce();
        m_data.insert(m_data.begin(), n, 0);
        return *this;
    }

    friend auto operator-(gfpoly /*a*/) -> gfpoly;

    friend auto operator*(const gfpoly & /*a*/, const gfpoly & /*b*/) -> gfpoly;

    friend auto operator/(const gfpoly & /*a*/, const gfpoly & /*b*/) -> gfpoly;

    friend auto operator%(const gfpoly & /*a*/, const gfpoly & /*b*/) -> gfpoly;

    friend auto operator==(const gfpoly & /*a*/, const gfpoly & /*b*/) -> bool;

    friend auto operator!=(const gfpoly & /*a*/, const gfpoly & /*b*/) -> bool;
};

inline
auto operator-(gfpoly a) -> gfpoly {
    std::transform(a.m_data.begin(), a.m_data.end(), a.m_data.begin(),
                   [&](uintmax_t x) { return a.neg(x); });
    return a.reduce();
}

inline
auto operator+(const gfpoly &a, const gfpoly &b) -> gfpoly {
    gfpoly result(a);
    result += b;
    return result;
}

inline
auto operator+(gfpoly &&a, const gfpoly &b) -> gfpoly {
    a += b;
    return a;
}

inline
auto operator+(const gfpoly &a, gfpoly &&b) -> gfpoly {
    b += a;
    return b;
}

inline
auto operator+(gfpoly &&a, gfpoly &&b) -> gfpoly {
    a += b;
    return a;
}

inline
auto operator-(const gfpoly &a, const gfpoly &b) -> gfpoly {
    gfpoly result(a);
    result -= b;
    return result;
}

inline
auto operator-(gfpoly &&a, const gfpoly &b) -> gfpoly {
    a -= b;
    return a;
}

inline
auto operator-(const gfpoly &a, gfpoly &&b) -> gfpoly {
    b -= a;
    return -b;
}

inline
auto operator-(gfpoly &&a, gfpoly &&b) -> gfpoly {
    a -= b;
    return a;
}

inline
auto operator*(const gfpoly &a, const gfpoly &b) -> gfpoly {
    gfpoly result(a.field());
    return result.multiply(a, b);
}

inline
auto operator/(const gfpoly &a, const gfpoly &b) -> gfpoly {
    return gfpoly::quotient_remainder(a, b).first;
}

inline
auto operator%(const gfpoly &a, const gfpoly &b) -> gfpoly {
    return gfpoly::quotient_remainder(a, b).second;
}

template<class U>
inline
auto operator+(gfpoly a, const U &b) -> gfpoly {
    a += b;
    return a;
}

template<class U>
inline
auto operator-(gfpoly a, const U &b) -> gfpoly {
    a -= b;
    return a;
}

template<class U>
inline
auto operator*(gfpoly a, const U &b) -> gfpoly {
    a *= b;
    return a;
}

template<class U>
inline
auto operator/(gfpoly a, const U &b) -> gfpoly {
    a /= b;
    return a;
}

template<class U>
inline
auto operator%(const gfpoly &a, const U & /*unused*/) -> gfpoly {
    return gfpoly(a.field());
}

template<class U>
inline
auto operator+(const U &a, gfpoly b) -> gfpoly {
    b += a;
    return b;
}

template<class U>
inline
auto operator-(const U &a, gfpoly b) -> gfpoly {
    b -= a;
    return -b;
}

template<class U>
inline
auto operator*(const U &a, gfpoly b) -> gfpoly {
    b *= a;
    return b;
}

inline
auto operator==(const gfpoly &a, const gfpoly &b) -> bool {
    return a.m_data == b.m_data;
}

inline
auto operator!=(const gfpoly &a, const gfpoly &b) -> bool {
    return a.m_data != b.m_data;
}

template<typename U>
inline
auto operator>>(gfpoly a, const U &b) -> gfpoly {
    a >>= b;
    return a;
}

template<typename U>
inline
auto operator<<(gfpoly a, const U &b) -> gfpoly {
    a <<= b;
    return a;
}

template<class charT, class traits>
auto operator<<(std::basic_ostream<charT, traits> &os, const gfpoly &poly)
-> std::basic_ostream<charT, traits> & {
    os << "{ ";
    for (uintmax_t i = 0; i < poly.size(); ++i) {
        if (i) {
            os << ", ";
        }
        os << poly[i];
    }
    os << " }";
    return os;
}

template<class charT, class traits>
auto operator>>(std::basic_istream<charT, traits> &is, gfpoly &poly)
-> std::basic_istream<charT, traits> & {
    charT tmp;
    uintmax_t num = 0;
    std::string str;
    std::vector<uintmax_t> vec;
    while (is.good() && is.get(tmp) && tmp != '{') {
        if (tmp != ' ' && tmp != '\n') {
            throw std::invalid_argument("wrong input");
        }
    }
    if (tmp != '{') {
        throw std::invalid_argument("wrong input");
    }
    while (is.good() && is.get(tmp) && tmp != '}') {
        if (tmp == ',' || tmp == ' ' || tmp == '\n') {
            if (!str.empty()) {
                std::stringstream(str) >> num;
                str.clear();
                vec.emplace_back(num);
            }
        } else if (std::isdigit(tmp)) {
            str += tmp;
        } else {
            throw std::invalid_argument("wrong input");
        }
    }
    if (tmp != '}') {
        throw std::invalid_argument("wrong input");
    }
    poly = gfpoly(poly.field(), vec);
    return is;
}

} // namespace irrpoly
