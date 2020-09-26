/**
 * @file    gf.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#pragma once

#include "nn.hpp"

#include <algorithm>
#include <utility>
#include <vector>
#include <cstdint>
#include <ostream>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <memory>

namespace irrpoly {

/**
 * Binary operations for two gfn instances are correctly defined only
 * when field is the same for both of them. By default this is checked
 * only in Debug configuration and no checks performed in Release to speed
 * up computations. If you are not sure in correctness of your code add
 * #define IRRPOLY_RELEASE_CHECKED before #include <irrpoly.h> to enable
 * checks for Release configuration.
 */
#if !defined(NDEBUG) || defined(IRRPOLY_RELEASE_CHECKED) // Debug or Release Checked
#define CHECK_FIELD(comparison) \
    if (!(comparison)) { \
        throw std::logic_error("field check failed"); \
    }
#else // Release
#define CHECK_FIELD(comparison)
#endif

class gfbase;

/**
 * gf type represents PRIME Galois field. It is a shared pointer that couldn't
 * contain nullptr value that allows to get rid of null checks at runtime.
 * gf instance could be constructed only then associated Galois field exists.
 * This fact is checked by calculating multiplicative inverse for every elements.
 * In PRIME field all multiplicative inverse elements must exist. The smallest
 * field you can create is GF[2], the largest is GF[65521] for 32-bit architecture
 * and GF[4294967291] for 64-bit architecture. gfn instance must always be passed
 * by reference and copied at the very last moment.
 * TODO: add support for COMPOSITE fields P^N.
 */
using gf = dropbox::oxygen::nn_shared_ptr<gfbase>;

class gfbase final {
private:
    const uintmax_t m_base; ///< field base, always could be converted to intmax_t
    std::vector<uintmax_t> m_inv; ///< multiplicative inverses for all elements

    explicit
    gfbase(uintmax_t /*base*/);

    friend auto make_gf(uintmax_t /*base*/) -> gf;

public:
    [[nodiscard]]
    auto base() const -> uintmax_t;

    /**
     * Returns multiplicative inverse for given number.
     */
    [[nodiscard]]
    auto mul_inv(uintmax_t /*val*/) const -> uintmax_t;
};

auto operator==(const gf &lb, const gf &rb) -> bool {
    return lb->base() == rb->base();
}

auto operator!=(const gf &lb, const gf &rb) -> bool {
    return lb->base() != rb->base();
}

/**
 * gfn type represents number in GF[P]. The number is always within 0 and P-1.
 */
class gfn final {
private:
    gf m_field;
    uintmax_t m_val;

public:
    /**
     * Generates a random number from  range [0, P-1].
     */
    static
    auto random(const gf &field) -> gfn {
        static std::random_device rd;
#ifdef __LP64__
        static std::mt19937_64 gen(rd());
#else
        static std::mt19937 gen(rd());
#endif
        std::uniform_int_distribution<uint_fast64_t> dis(0, field->base() - 1);
        return gfn(field, dis(gen));
    }

    [[nodiscard]]
    auto value() const -> uintmax_t {
        return m_val;
    }

    explicit
    gfn(const gf &field) : m_field(field), m_val(0) {}

    gfn(const gf &field, const uintmax_t val) :
        m_field(field), m_val(val % m_field->base()) {}

    gfn(const gfn &other) = default;

    gfn(gfn &&other) = default;

    auto operator=(const gfn &other) -> gfn & {
        if (this != &other) {
            // m_field == nullptr means that gfn instance is uninitialised
            // this happens during std::move, std::swap and inside some std::vector methods
            CHECK_FIELD(m_field == nullptr || m_field == other.m_field)
            m_field = other.m_field;
            m_val = other.m_val;
        }
        return *this;
    }

    [[nodiscard]]
    auto base() const -> uintmax_t {
        return m_field->base();
    }

    auto operator=(const uintmax_t other) -> gfn & {
        m_val = other % base();
        return *this;
    }

    [[nodiscard]]
    auto field() const -> const gf & {
        return m_field;
    }

    [[nodiscard]]
    auto operator+() const -> gfn {
        return gfn{*this};
    }

    [[nodiscard]]
    auto operator+(const gfn &other) const -> gfn {
        CHECK_FIELD(m_field == other.field())
        return gfn{m_field, m_val + other.m_val};
    }

    [[nodiscard]]
    auto operator+(const uintmax_t other) const -> gfn {
        return gfn{m_field, m_val + (other % base())};
    }

    friend
    auto operator+(uintmax_t /*other*/, const gfn & /*curr*/) -> gfn;

    auto operator+=(const gfn &other) -> gfn & {
        CHECK_FIELD(m_field == other.field())
        m_val = (m_val + other.m_val) % base();
        return *this;
    }

    auto operator+=(const uintmax_t other) -> gfn {
        m_val = (m_val + (other % base())) % base();
        return *this;
    }

    auto operator++() -> gfn & {
        m_val = (m_val + 1) % base();
        return *this;
    }

    [[nodiscard]]
    auto operator++(int) & -> gfn {
        gfn tmp{*this};
        m_val = (m_val + 1) % base();
        return tmp;
    }

    [[nodiscard]]
    auto operator-() const -> gfn {
        gfn tmp{*this};
        tmp.m_val = (base() - m_val) % base();
        return tmp;
    }

    [[nodiscard]]
    auto operator-(const gfn &other) const -> gfn {
        CHECK_FIELD(m_field == other.field())
        return gfn{m_field, base() + m_val - other.m_val};
    }

    [[nodiscard]]
    auto operator-(const uintmax_t other) const -> gfn {
        return gfn{m_field, base() + m_val - (other % base())};
    }

    friend
    auto operator-(uintmax_t /*other*/, const gfn & /*curr*/) -> gfn;

    auto operator-=(const gfn &other) -> gfn & {
        CHECK_FIELD(m_field == other.field())
        m_val = (base() + m_val - other.m_val) % base();
        return *this;
    }

    auto operator-=(const uintmax_t other) -> gfn {
        m_val = (base() + m_val - (other % base())) % base();
        return *this;
    }

    auto operator--() -> gfn & {
        m_val = (base() + m_val - 1) % base();
        return *this;
    }

    [[nodiscard]]
    auto operator--(int) & -> gfn {
        gfn tmp{*this};
        m_val = (base() + m_val - 1) % base();
        return tmp;
    }

    [[nodiscard]]
    auto operator*(const gfn &other) const -> gfn {
        CHECK_FIELD(m_field == other.field())
        return gfn{m_field, m_val * other.m_val};
    }

    [[nodiscard]]
    auto operator*(const uintmax_t other) const -> gfn {
        return gfn{m_field, m_val * (other % base())};
    }

    friend
    auto operator*(uintmax_t /*other*/, const gfn & /*curr*/) -> gfn;

    auto operator*=(const gfn &other) -> gfn & {
        CHECK_FIELD(m_field == other.field())
        m_val = (m_val * other.m_val) % base();
        return *this;
    }

    auto operator*=(const uintmax_t other) -> gfn {
        m_val = (m_val * (other % base())) % base();
        return *this;
    }

    /**
     * Returns multiplicative inverse for current gfn instance.
     */
    [[maybe_unused]] [[nodiscard]]
    auto mul_inv() -> gfn {
        return gfn{m_field, m_field->mul_inv(m_val)};
    }

    /**
     * Division is defined as multiplication by multiplicative inverse.
     */
    [[nodiscard]]
    auto operator/(const gfn &other) const -> gfn {
        CHECK_FIELD(m_field == other.field())
        switch (other.m_val) {
        case 0:throw std::invalid_argument("division by zero");
        default:return gfn{m_field, m_val * m_field->mul_inv(other.m_val)};
        }
    }

    [[nodiscard]]
    auto operator/(const uintmax_t other) const -> gfn {
        switch (other % base()) {
        case 0:throw std::invalid_argument("division by zero");
        default:
            return gfn{m_field,
                       m_val * m_field->mul_inv(other % base())};
        }
    }

    friend
    auto operator/(uintmax_t /*other*/, const gfn & /*curr*/) -> gfn;

    auto operator/=(const gfn &other) -> gfn & {
        CHECK_FIELD(m_field == other.field())
        switch (other.m_val) {
        case 0:throw std::invalid_argument("division by zero");
        default:m_val = (m_val * m_field->mul_inv(other.m_val)) % base();
            return *this;
        }
    }

    auto operator/=(const uintmax_t other) -> gfn {
        switch (other % base()) {
        case 0:throw std::invalid_argument("division by zero");
        default:m_val = (m_val * m_field->mul_inv(other % base())) % base();
            return *this;
        }
    }

    [[nodiscard]]
    auto is_zero() const -> bool {
        return 0 == m_val;
    }

    explicit operator bool() const {
        return 0 != m_val;
    }

    template<class charT, class traits>
    friend
    auto operator<<(std::basic_ostream<charT, traits> & /*os*/, const gfn & /*val*/)
    -> std::basic_ostream<charT, traits> &;

    template<class charT, class traits>
    friend
    auto operator>>(std::basic_istream<charT, traits> & /*is*/, gfn & /*val*/)
    -> std::basic_istream<charT, traits> &;
};

#define GFN_COMPARISON_OPERATORS(op) \
    inline \
    auto operator op(const gfn &l, const gfn &r) -> bool { \
        CHECK_FIELD(l.field() == r.field()) \
        return l.value() op r.value(); \
    } \
    \
    inline \
    auto operator op(const gfn &l, const uintmax_t r) -> bool { \
        return l.value() op (r % l.base()); \
    } \
    \
    inline \
    auto operator op(const uintmax_t l, const gfn &r) -> bool { \
        return (l % r.base()) op r.value(); \
    }

GFN_COMPARISON_OPERATORS(==)
GFN_COMPARISON_OPERATORS(!=)
GFN_COMPARISON_OPERATORS(<)
GFN_COMPARISON_OPERATORS(<=)
GFN_COMPARISON_OPERATORS(>)
GFN_COMPARISON_OPERATORS(>=)

#undef GFN_COMPARISON_OPERATORS

[[nodiscard]]
auto operator+(const uintmax_t other, const gfn &curr) -> gfn {
    return gfn{curr.m_field, (other % curr.base()) + curr.m_val};
}

[[nodiscard]]
auto operator-(const uintmax_t other, const gfn &curr) -> gfn {
    return gfn{curr.m_field, curr.base() + (other % curr.base()) - curr.m_val};
}

[[nodiscard]]
auto operator*(const uintmax_t other, const gfn &curr) -> gfn {
    return gfn{curr.m_field, (other % curr.base()) * curr.m_val};
}

[[nodiscard]]
auto operator/(const uintmax_t other, const gfn &curr) -> gfn {
    switch (curr.m_val) {
    case 0:throw std::invalid_argument("division by zero");
    default:
        return gfn{curr.m_field,
                   (other % curr.base()) * curr.m_field->mul_inv(curr.m_val)};
    }
}

template<class charT, class traits>
auto operator<<(std::basic_ostream<charT, traits> &os, const gfn &val)
-> std::basic_ostream<charT, traits> & {
    return os << val.m_val;
}

template<class charT, class traits>
auto operator>>(std::basic_istream<charT, traits> &is, gfn &val)
-> std::basic_istream<charT, traits> & {
    is >> val.m_val;
    val.m_val %= val.base();
    return is;
}

inline
gfbase::gfbase(const uintmax_t base) : m_base(base), m_inv(base, 0) {
    if (base == 0) {
        throw std::logic_error("empty field");
    }
    if (base == 1) {
        throw std::logic_error("field could contain only zero");
    }
    if (UINTMAX_MAX / (base - 1) < (base - 1)) {
        throw std::logic_error("too large field");
    }

    auto i_base = static_cast<intmax_t>(base);
    auto inv_calc = [](const intmax_t base, const intmax_t val) -> uintmax_t {
        intmax_t u0 = base, u1 = 1, u2 = 0,
            v0 = val, v1 = 0, v2 = 1, w0 = 0, w1 = 0, w2 = 0, q = 0;
        while (v0 > 0) {
            q = u0 / v0;
            w0 = u0 - q * v0, w1 = u1 - q * v1, w2 = u2 - q * v2;
            u0 = v0, u1 = v1, u2 = v2, v0 = w0, v1 = w1, v2 = w2;
        }
        if (u0 > 1) {
            throw std::logic_error("multiplicative inverse don't exist");
        }
        return static_cast<uintmax_t>(u2 < 0 ? (base + u2) : (u2));
    };

    m_inv[1] = 1;
    for (uintmax_t i = 2; i < m_base; ++i) {
        if (m_inv[i]) {
            continue;
        }
        m_inv[i] = inv_calc(i_base, i);
        m_inv[m_inv[i]] = i;
    }
}

[[nodiscard]]
inline
auto gfbase::base() const -> uintmax_t {
    return m_base;
}

[[nodiscard]]
inline
auto gfbase::mul_inv(const uintmax_t val) const -> uintmax_t {
    switch (val % m_base) {
    case 0:throw std::logic_error("multiplicative inverse don't exist");
    default:return m_inv[val % m_base];
    }
}

[[nodiscard]]
inline
auto make_gf(const uintmax_t base) -> gf {
    return dropbox::oxygen::nn< std::shared_ptr<gfbase> >(dropbox::oxygen::nn(
        dropbox::oxygen::i_promise_i_checked_for_null_t{}, new gfbase(base)));
}

#undef CHECK_FIELD

} // namespace irrpoly
