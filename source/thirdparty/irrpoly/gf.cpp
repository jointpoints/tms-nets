#include "../../../include/tms-nets/thirdparty/irrpoly/gf.hpp"





// Comparison operators for <gf> objects





auto irrpoly::operator==(const gf &lb, const gf &rb) -> bool {
    return lb->base() == rb->base();
}



auto irrpoly::operator!=(const gf &lb, const gf &rb) -> bool {
    return lb->base() != rb->base();
}





// Class <gfn>





[[nodiscard]]
auto irrpoly::gfn::value() const -> uintmax_t {
	return m_val;
}



irrpoly::gfn::gfn(const gf &field) : m_field(field), m_val(0) {}



irrpoly::gfn::gfn(const gf &field, const uintmax_t val) :
	m_field(field), m_val(val % m_field->base()) {}



irrpoly::gfn::gfn(const gfn &other) = default;



irrpoly::gfn::gfn(gfn &&other) = default;



auto irrpoly::gfn::operator=(const gfn &other) -> gfn & {
	if (this != &other) {
		// m_field == nullptr means that gfn instance is uninitialised
		// this happens during std::move, std::swap and inside some std::vector methods
		m_field = other.m_field;
		m_val = other.m_val;
	}
	return *this;
}



[[nodiscard]]
auto irrpoly::gfn::base() const -> uintmax_t {
	return m_field->base();
}



auto irrpoly::gfn::operator=(const uintmax_t other) -> gfn & {
	m_val = other % base();
	return *this;
}



[[nodiscard]]
auto irrpoly::gfn::field() const -> const gf & {
	return m_field;
}



[[nodiscard]]
auto irrpoly::gfn::operator+() const -> gfn {
	return gfn{*this};
}



[[nodiscard]]
auto irrpoly::gfn::operator+(const gfn &other) const -> gfn {
	return gfn{m_field, m_val + other.m_val};
}



[[nodiscard]]
auto irrpoly::gfn::operator+(const uintmax_t other) const -> gfn {
	return gfn{m_field, m_val + (other % base())};
}



auto irrpoly::gfn::operator+=(const gfn &other) -> gfn & {
	m_val = (m_val + other.m_val) % base();
	return *this;
}



auto irrpoly::gfn::operator+=(const uintmax_t other) -> gfn {
	m_val = (m_val + (other % base())) % base();
	return *this;
}



auto irrpoly::gfn::operator++() -> gfn & {
	m_val = (m_val + 1) % base();
	return *this;
}



[[nodiscard]]
auto irrpoly::gfn::operator++(int) & -> gfn {
	gfn tmp{*this};
	m_val = (m_val + 1) % base();
	return tmp;
}



[[nodiscard]]
auto irrpoly::gfn::operator-() const -> gfn {
	gfn tmp{*this};
	tmp.m_val = (base() - m_val) % base();
	return tmp;
}



[[nodiscard]]
auto irrpoly::gfn::operator-(const gfn &other) const -> gfn {
	return gfn{m_field, base() + m_val - other.m_val};
}



[[nodiscard]]
auto irrpoly::gfn::operator-(const uintmax_t other) const -> gfn {
	return gfn{m_field, base() + m_val - (other % base())};
}



auto irrpoly::gfn::operator-=(const gfn &other) -> gfn & {
	m_val = (base() + m_val - other.m_val) % base();
	return *this;
}



auto irrpoly::gfn::operator-=(const uintmax_t other) -> gfn {
	m_val = (base() + m_val - (other % base())) % base();
	return *this;
}



auto irrpoly::gfn::operator--() -> gfn & {
	m_val = (base() + m_val - 1) % base();
	return *this;
}



[[nodiscard]]
auto irrpoly::gfn::operator--(int) & -> gfn {
	gfn tmp{*this};
	m_val = (base() + m_val - 1) % base();
	return tmp;
}



[[nodiscard]]
auto irrpoly::gfn::operator*(const gfn &other) const -> gfn {
	return gfn{m_field, m_val * other.m_val};
}



[[nodiscard]]
auto irrpoly::gfn::operator*(const uintmax_t other) const -> gfn {
	return gfn{m_field, m_val * (other % base())};
}



auto irrpoly::gfn::operator*=(const gfn &other) -> gfn & {
	m_val = (m_val * other.m_val) % base();
	return *this;
}



auto irrpoly::gfn::operator*=(const uintmax_t other) -> gfn {
	m_val = (m_val * (other % base())) % base();
	return *this;
}



[[maybe_unused]] [[nodiscard]]
auto irrpoly::gfn::mul_inv() -> gfn {
	return gfn{m_field, m_field->mul_inv(m_val)};
}



[[nodiscard]]
auto irrpoly::gfn::operator/(const gfn &other) const -> gfn {
	switch (other.m_val) {
	case 0:throw std::invalid_argument("division by zero");
	default:return gfn{m_field, m_val * m_field->mul_inv(other.m_val)};
	}
}



[[nodiscard]]
auto irrpoly::gfn::operator/(const uintmax_t other) const -> gfn {
	switch (other % base()) {
	case 0:throw std::invalid_argument("division by zero");
	default:
		return gfn{m_field,
					m_val * m_field->mul_inv(other % base())};
	}
}



auto irrpoly::gfn::operator/=(const gfn &other) -> gfn & {
	switch (other.m_val) {
	case 0:throw std::invalid_argument("division by zero");
	default:m_val = (m_val * m_field->mul_inv(other.m_val)) % base();
		return *this;
	}
}



auto irrpoly::gfn::operator/=(const uintmax_t other) -> gfn {
	switch (other % base()) {
	case 0:throw std::invalid_argument("division by zero");
	default:m_val = (m_val * m_field->mul_inv(other % base())) % base();
		return *this;
	}
}



[[nodiscard]]
auto irrpoly::gfn::is_zero() const -> bool {
	return 0 == m_val;
}



irrpoly::gfn::operator bool() const {
	return 0 != m_val;
}





// External operators for <gfn>





[[nodiscard]]
auto irrpoly::operator+(const uintmax_t other, const gfn &curr) -> gfn {
    return gfn{curr.m_field, (other % curr.base()) + curr.m_val};
}



[[nodiscard]]
auto irrpoly::operator-(const uintmax_t other, const gfn &curr) -> gfn {
    return gfn{curr.m_field, curr.base() + (other % curr.base()) - curr.m_val};
}



[[nodiscard]]
auto irrpoly::operator*(const uintmax_t other, const gfn &curr) -> gfn {
    return gfn{curr.m_field, (other % curr.base()) * curr.m_val};
}



[[nodiscard]]
auto irrpoly::operator/(const uintmax_t other, const gfn &curr) -> gfn {
    switch (curr.m_val) {
    case 0:throw std::invalid_argument("division by zero");
    default:
        return gfn{curr.m_field,
                   (other % curr.base()) * curr.m_field->mul_inv(curr.m_val)};
    }
}
