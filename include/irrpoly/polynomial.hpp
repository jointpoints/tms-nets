/**
 * @file    boost/math/tools/polynomial.hpp
 * @author  John Maddock 2006, Jeremy William Murphy 2015
 * @license Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * @url     https://www.boost.org/doc/libs/1_71_0/
 */

#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <vector>
#include <cassert>
#include <ostream>
#include <algorithm>
#include <initializer_list>

namespace irrpoly {

    /**
     * Реализация класса polynomial взыта из библиотеки Boost 1.71.0.
     * За комментариями к реализации обращайтесь к первоисточнику
     * https://www.boost.org/doc/libs/1_71_0/boost/math/tools/polynomial.hpp
     */
    template<typename T>
    class polynomial;

    namespace detail {

        template<typename T, typename N>
        void
        division_impl(polynomial<T> &q, polynomial<T> &u, const polynomial<T> &v, N n, N k) {
            q[k] = u[n + k] / v[n];
            for (N j = n + k; j > k;) {
                j--;
                u[j] -= q[k] * v[j - k];
            }
        }

        template<class T, class N>
        T integer_power(T t, N n) {
            switch (n) {
                case 0:
                    return static_cast<T>(1u);
                case 1:
                    return t;
                case 2:
                    return t * t;
                case 3:
                    return t * t * t;
            }
            T result = integer_power(t, n / 2);
            result *= result;
            if (n & 1)
                result *= t;
            return result;
        }

        template<typename T>
        ::std::pair<polynomial<T>, polynomial<T> >
        division(polynomial<T> u, const polynomial<T> &v) {
            assert(v.size() <= u.size());
            assert(v);
            assert(u);

            typedef typename polynomial<T>::size_type N;

            N const m = u.size() - 1, n = v.size() - 1;
            N k = m - n;
            polynomial<T> q;
            q.data().resize(m - n + 1);

            do {
                division_impl(q, u, v, n, k);
            } while (k-- != 0);
            u.data().resize(n);
            u.normalize();
            return ::std::make_pair(q, u);
        }

        struct negate {
            template<class T>
            T operator()(T const &x) const {
                return -x;
            }
        };

        struct plus {
            template<class T, class U>
            T operator()(T const &x, U const &y) const {
                return x + y;
            }
        };

        struct minus {
            template<class T, class U>
            T operator()(T const &x, U const &y) const {
                return x - y;
            }
        };

    } // namespace detail

    /**
    * Returns the zero element for multiplication of polynomials.
    */
    template<class T>
    polynomial<T> zero_element(::std::multiplies<polynomial<T> >) {
        return polynomial<T>();
    }

    template<class T>
    polynomial<T> identity_element(::std::multiplies<polynomial<T> >) {
        return polynomial<T>(T(1));
    }

    /* Calculates a / b and a % b, returning the pair (quotient, remainder) together
    * because the same amount of computation yields both.
    * This function is not defined for division by zero: user beware.
    */
    template<typename T>
    ::std::pair<polynomial<T>, polynomial<T> >
    quotient_remainder(const polynomial<T> &dividend, const polynomial<T> &divisor) {
        assert(divisor);
        if (dividend.size() < divisor.size())
            return ::std::make_pair(polynomial<T>(), dividend);
        return detail::division(dividend, divisor);
    }


    template<class T>
    class polynomial {
    public:
        typedef typename ::std::vector<T>::value_type value_type;
        typedef typename ::std::vector<T>::size_type size_type;

        polynomial() {}

        template<class U>
        polynomial(const U *data, unsigned order)
                : m_data(data, data + order + 1) {
            normalize();
        }

        template<class I>
        polynomial(I first, I last)
                : m_data(first, last) {
            normalize();
        }

        polynomial(::std::vector<T> &&p) : m_data(::std::move(p)) {
            normalize();
        }

        template<class U>
        explicit polynomial(const U &point, U * = 0) {
            if (point != U(0))
                m_data.push_back(point);
        }

        polynomial(polynomial &&p) noexcept
                : m_data(::std::move(p.m_data)) {}

        polynomial(const polynomial &p)
                : m_data(p.m_data) {}

        template<class U>
        polynomial(const polynomial<U> &p) {
            m_data.resize(p.size());
            for (unsigned i = 0; i < p.size(); ++i) {
                m_data[i] = static_cast<T>(p[i]);
            }
        }

        polynomial(::std::initializer_list<T> l) : polynomial(::std::begin(l), ::std::end(l)) {
        }

        polynomial &
        operator=(::std::initializer_list<T> l) {
            m_data.assign(::std::begin(l), ::std::end(l));
            normalize();
            return *this;
        }

        size_type size() const { return m_data.size(); }

        size_type degree() const {
            if (size() == 0)
                throw ::std::logic_error("degree() is undefined for the zero polynomial");
            return m_data.size() - 1;
        }

        value_type &operator[](size_type i) {
            return m_data[i];
        }

        const value_type &operator[](size_type i) const {
            return m_data[i];
        }

        ::std::vector<T> const &data() const {
            return m_data;
        }

        ::std::vector<T> &data() {
            return m_data;
        }

        polynomial &operator=(polynomial &&p) noexcept {
            m_data = ::std::move(p.m_data);
            return *this;
        }

        polynomial &operator=(const polynomial &p) {
            m_data = p.m_data;
            return *this;
        }

        template<class U>
        polynomial &operator+=(const U &value) {
            addition(value);
            normalize();
            return *this;
        }

        template<class U>
        polynomial &operator-=(const U &value) {
            subtraction(value);
            normalize();
            return *this;
        }

        template<class U>
        polynomial &operator*=(const U &value) {
            multiplication(value);
            normalize();
            return *this;
        }

        template<class U>
        polynomial &operator/=(const U &value) {
            division(value);
            normalize();
            return *this;
        }

        template<class U>
        polynomial &
        operator%=(const U & /*value_type*/) {
            // We can always divide by a scalar, so there is no remainder:
            this->set_zero();
            return *this;
        }

        template<class U>
        polynomial &operator+=(const polynomial<U> &value) {
            addition(value);
            normalize();
            return *this;
        }

        template<class U>
        polynomial &operator-=(const polynomial<U> &value) {
            subtraction(value);
            normalize();
            return *this;
        }

        template<typename U, typename V>
        void multiply(const polynomial<U> &a, const polynomial<V> &b) {
            if (!a || !b) {
                this->set_zero();
                return;
            }
            ::std::vector<T> prod(a.size() + b.size() - 1, T(0));
            for (typename ::std::vector<T>::size_type i = 0; i < a.size(); ++i)
                for (typename ::std::vector<T>::size_type j = 0; j < b.size(); ++j)
                    prod[i + j] += a.m_data[i] * b.m_data[j];
            m_data.swap(prod);
        }

        template<class U>
        polynomial &operator*=(const polynomial<U> &value) {
            this->multiply(*this, value);
            return *this;
        }

        template<typename U>
        polynomial &operator/=(const polynomial<U> &value) {
            *this = quotient_remainder(*this, value).first;
            return *this;
        }

        template<typename U>
        polynomial &operator%=(const polynomial<U> &value) {
            *this = quotient_remainder(*this, value).second;
            return *this;
        }

        template<typename U>
        polynomial &operator>>=(U const &n) {
            assert(n <= m_data.size());
            m_data.erase(m_data.begin(), m_data.begin() + n);
            return *this;
        }

        template<typename U>
        polynomial &operator<<=(U const &n) {
            m_data.insert(m_data.begin(), n, static_cast<T>(0));
            normalize();
            return *this;
        }

        [[nodiscard]]
        bool is_zero() const {
            return m_data.empty();
        }

        explicit operator bool() const {
            return !m_data.empty();
        }

        void set_zero() {
            m_data.clear();
        }

        void normalize() {
            m_data.erase(
                    ::std::find_if(m_data.rbegin(), m_data.rend(), [](const T &x) -> bool { return x != T(0); }).base(),
                    m_data.end());
        }

    private:
        template<class U, class R>
        polynomial &addition(const U &value, R op) {
            if (m_data.size() == 0)
                m_data.resize(1, 0);
            m_data[0] = op(m_data[0], value);
            return *this;
        }

        template<class U>
        polynomial &addition(const U &value) {
            return addition(value, detail::plus());
        }

        template<class U>
        polynomial &subtraction(const U &value) {
            return addition(value, detail::minus());
        }

        template<class U, class R>
        polynomial &addition(const polynomial<U> &value, R op) {
            if (m_data.size() < value.size())
                m_data.resize(value.size(), 0);
            for (size_type i = 0; i < value.size(); ++i)
                m_data[i] = op(m_data[i], value[i]);
            return *this;
        }

        template<class U>
        polynomial &addition(const polynomial<U> &value) {
            return addition(value, detail::plus());
        }

        template<class U>
        polynomial &subtraction(const polynomial<U> &value) {
            return addition(value, detail::minus());
        }

        template<class U>
        polynomial &multiplication(const U &value) {
            ::std::transform(m_data.begin(), m_data.end(), m_data.begin(), [&](const T &x) -> T { return x * value; });
            return *this;
        }

        template<class U>
        polynomial &division(const U &value) {
            ::std::transform(m_data.begin(), m_data.end(), m_data.begin(), [&](const T &x) -> T { return x / value; });
            return *this;
        }

        ::std::vector<T> m_data;
    };


    template<class T>
    polynomial<T> operator+(const polynomial<T> &a, const polynomial<T> &b) {
        polynomial<T> result(a);
        result += b;
        return result;
    }

    template<class T>
    polynomial<T> operator+(polynomial<T> &&a, const polynomial<T> &b) {
        a += b;
        return a;
    }

    template<class T>
    polynomial<T> operator+(const polynomial<T> &a, polynomial<T> &&b) {
        b += a;
        return b;
    }

    template<class T>
    polynomial<T> operator+(polynomial<T> &&a, polynomial<T> &&b) {
        a += b;
        return a;
    }

    template<class T>
    polynomial<T> operator-(const polynomial<T> &a, const polynomial<T> &b) {
        polynomial<T> result(a);
        result -= b;
        return result;
    }

    template<class T>
    polynomial<T> operator-(polynomial<T> &&a, const polynomial<T> &b) {
        a -= b;
        return a;
    }

    template<class T>
    polynomial<T> operator-(const polynomial<T> &a, polynomial<T> &&b) {
        b -= a;
        return -b;
    }

    template<class T>
    polynomial<T> operator-(polynomial<T> &&a, polynomial<T> &&b) {
        a -= b;
        return a;
    }

    template<class T>
    polynomial<T> operator*(const polynomial<T> &a, const polynomial<T> &b) {
        polynomial<T> result;
        result.multiply(a, b);
        return result;
    }

    template<class T>
    polynomial<T> operator/(const polynomial<T> &a, const polynomial<T> &b) {
        return quotient_remainder(a, b).first;
    }

    template<class T>
    polynomial<T> operator%(const polynomial<T> &a, const polynomial<T> &b) {
        return quotient_remainder(a, b).second;
    }

    template<class T, class U>
    polynomial<T>
    operator+(polynomial<T> a, const U &b) {
        a += b;
        return a;
    }

    template<class T, class U>
    polynomial<T>
    operator-(polynomial<T> a, const U &b) {
        a -= b;
        return a;
    }

    template<class T, class U>
    polynomial<T>
    operator*(polynomial<T> a, const U &b) {
        a *= b;
        return a;
    }

    template<class T, class U>
    polynomial<T>
    operator/(polynomial<T> a, const U &b) {
        a /= b;
        return a;
    }

    template<class T, class U>
    polynomial<T>
    operator%(const polynomial<T> &, const U &) {
        return polynomial<T>();
    }

    template<class U, class T>
    polynomial<T>
    operator+(const U &a, polynomial<T> b) {
        b += a;
        return b;
    }

    template<class U, class T>
    polynomial<T>
    operator-(const U &a, polynomial<T> b) {
        b -= a;
        return -b;
    }

    template<class U, class T>
    polynomial<T>
    operator*(const U &a, polynomial<T> b) {
        b *= a;
        return b;
    }

    template<class T>
    bool operator==(const polynomial<T> &a, const polynomial<T> &b) {
        return a.data() == b.data();
    }

    template<class T>
    bool operator!=(const polynomial<T> &a, const polynomial<T> &b) {
        return a.data() != b.data();
    }

    template<typename T, typename U>
    polynomial<T> operator>>(polynomial<T> a, const U &b) {
        a >>= b;
        return a;
    }

    template<typename T, typename U>
    polynomial<T> operator<<(polynomial<T> a, const U &b) {
        a <<= b;
        return a;
    }

    template<class T>
    polynomial<T> operator-(polynomial<T> a) {
        ::std::transform(a.data().begin(), a.data().end(), a.data().begin(), detail::negate());
        return a;
    }

    template<class charT, class traits, class T>
    ::std::basic_ostream<charT, traits> &
    operator<<(::std::basic_ostream<charT, traits> &os, const polynomial<T> &poly) {
        os << "{ ";
        for (unsigned i = 0; i < poly.size(); ++i) {
            if (i) os << ", ";
            os << poly[i];
        }
        os << " }";
        return os;
    }

} // namespace irrpoly

#endif //POLYNOMIAL_HPP
