#include "../../../include/tms-nets/thirdparty/irrpoly/gfpoly.hpp"





// Class <gfpoly>





[[nodiscard]]
auto irrpoly::gfpoly::value() const -> const std::vector<uintmax_t>& {
	return m_data;
}



auto irrpoly::gfpoly::reduce() -> gfpoly & {
	m_data.erase(std::find_if(
		m_data.rbegin(), m_data.rend(),
		[](uintmax_t x) { return x != 0; }
	).base(), m_data.end());
	return *this;
}



irrpoly::gfpoly::gfpoly(const gf &field) : m_field(field), m_data() {}



irrpoly::gfpoly::gfpoly(const gf &field, const std::vector<uintmax_t> &l) :
	m_field(field), m_data() {
	m_data.reserve(l.size());
	for (uintmax_t v : l) {
		m_data.push_back(v % base());
	}
	reduce();
}



irrpoly::gfpoly::gfpoly(const gf &field, std::vector<uintmax_t> &&l) :
	m_field(field), m_data(l) {
	for (uintmax_t &v : m_data) {
		v %= base();
	}
	reduce();
}



auto irrpoly::gfpoly::operator=(const std::vector<uintmax_t> &l) -> gfpoly & {
	gfpoly copy(m_field, l);
	std::swap(*this, copy);
	return *this;
}



irrpoly::gfpoly::gfpoly(const gf &field, std::initializer_list<uintmax_t> l) :
	gfpoly(field, std::move(std::vector<uintmax_t>{l})) {}



auto irrpoly::gfpoly::operator=(std::initializer_list<uintmax_t> l) -> gfpoly & {
	gfpoly copy(m_field, l);
	std::swap(*this, copy);
	return *this;
}



irrpoly::gfpoly::gfpoly(const gfpoly &p) = default;



irrpoly::gfpoly::gfpoly(gfpoly &&p) = default;



auto irrpoly::gfpoly::operator=(const gfpoly &p) -> gfpoly & {
	if (this != &p) {
		// m_field == nullptr means that gfn instance is uninitialised
		// this happens during std::move, std::swap and inside some std::vector methods
		m_field = p.m_field;
		m_data = p.m_data;
	}
	return *this;
}



irrpoly::gfpoly::gfpoly(const gfn &value) :
	m_field(value.field()), m_data() {
	if (value) {
		m_data.push_back(value.value());
	}
}



auto irrpoly::gfpoly::operator=(const gfn &value) -> gfpoly & {
	// m_field == nullptr means that gfn instance is uninitialised
	// this happens during std::move, std::swap and inside some std::vector methods
	gfpoly copy(value);
	std::swap(*this, copy);
	return *this;
}



irrpoly::gfpoly::gfpoly(const gf &field, uintmax_t value) :
	m_field(field), m_data() {
	if (value % base() != 0) {
		m_data.push_back(value % base());
	}
}



auto irrpoly::gfpoly::operator=(uintmax_t value) -> gfpoly & {
	gfpoly copy(m_field, value);
	std::swap(*this, copy);
	return *this;
}



[[nodiscard]]
auto irrpoly::gfpoly::field() const -> const gf & {
	return m_field;
}



[[nodiscard]]
auto irrpoly::gfpoly::base() const -> uintmax_t {
	return m_field->base();
}



[[nodiscard]]
auto irrpoly::gfpoly::size() const -> uintmax_t {
	return m_data.size();
}



[[nodiscard]]
auto irrpoly::gfpoly::degree() const -> uintmax_t {
	if (size() == 0) {
		throw std::logic_error("degree is undefined for zero polynomial");
	}
	return m_data.size() - 1;
}



auto irrpoly::gfpoly::operator[](uintmax_t i) const -> uintmax_t {
	return m_data[i];
}



[[nodiscard]]
auto irrpoly::gfpoly::is_zero() const -> bool {
	return m_data.empty();
}



irrpoly::gfpoly::operator bool() const {
	return !is_zero();
}



[[maybe_unused]]
auto irrpoly::gfpoly::set_zero() -> gfpoly & {
	m_data.clear();
	return *this;
}



[[nodiscard]]
auto irrpoly::gfpoly::add(uintmax_t lb, uintmax_t rb) const -> uintmax_t {
	return (lb + rb) % base();
}



[[nodiscard]]
auto irrpoly::gfpoly::sub(uintmax_t lb, uintmax_t rb) const -> uintmax_t {
	return (base() + lb - rb) % base();
}



[[nodiscard]]
auto irrpoly::gfpoly::mul(uintmax_t lb, uintmax_t rb) const -> uintmax_t {
	return (lb * rb) % base();
}



[[nodiscard]]
auto irrpoly::gfpoly::div(uintmax_t lb, uintmax_t rb) const -> uintmax_t {
	switch (rb) {
	case 0:throw std::invalid_argument("division by zero");
	default:return (lb * field()->mul_inv(rb)) % base();
	}
}



[[nodiscard]]
auto irrpoly::gfpoly::neg(uintmax_t rb) const -> uintmax_t {
	return (base() - rb) % base();
}



auto irrpoly::gfpoly::transform(uintmax_t value, OP op) -> gfpoly & {
	if (m_data.empty()) {
		m_data.resize(1, 0);
	}
	m_data[0] = std::invoke(op, this, m_data[0], value);
	return reduce();
}



auto irrpoly::gfpoly::transform(const gfn &value, OP op) -> gfpoly & {
	if (m_data.empty()) {
		m_data.resize(1, 0);
	}
	m_data[0] = std::invoke(op, this, m_data[0], value.value());
	return reduce();
}



auto irrpoly::gfpoly::transform(const gfpoly &value, OP op) -> gfpoly & {
	if (m_data.size() < value.size()) {
		m_data.resize(value.size(), 0);
	}
	for (uintmax_t i = 0; i < value.size(); ++i) {
		m_data[i] = std::invoke(op, this, m_data[i], value[i]);
	}
	return reduce();
}



auto irrpoly::gfpoly::operator+=(uintmax_t value) -> gfpoly & {
	return transform(value, &gfpoly::add);
}



auto irrpoly::gfpoly::operator+=(const gfn &value) -> gfpoly & {
	return transform(value.value(), &gfpoly::add);
}



auto irrpoly::gfpoly::operator+=(const gfpoly &value) -> gfpoly & {
	return transform(value, &gfpoly::add);
}



auto irrpoly::gfpoly::operator-=(uintmax_t value) -> gfpoly & {
	return transform(value, &gfpoly::sub);
}



auto irrpoly::gfpoly::operator-=(const gfn &value) -> gfpoly & {
	return transform(value.value(), &gfpoly::sub);
}



auto irrpoly::gfpoly::operator-=(const gfpoly &value) -> gfpoly & {
	return transform(value, &gfpoly::sub);
}



auto irrpoly::gfpoly::operator*=(uintmax_t value) -> gfpoly & {
	std::transform(m_data.begin(), m_data.end(), m_data.begin(),
					[&](uintmax_t x) -> uintmax_t { return mul(x, value); });
	return reduce();
}



auto irrpoly::gfpoly::operator*=(const gfn &value) -> gfpoly & {
	std::transform(m_data.begin(), m_data.end(), m_data.begin(),
					[&](uintmax_t x) -> uintmax_t { return mul(x, value.value()); });
	return reduce();
}



auto irrpoly::gfpoly::operator/=(uintmax_t value) -> gfpoly & {
	std::transform(m_data.begin(), m_data.end(), m_data.begin(),
					[&](uintmax_t x) -> uintmax_t { return div(x, value); });
	return reduce();
}



auto irrpoly::gfpoly::operator/=(const gfn &value) -> gfpoly & {
	std::transform(m_data.begin(), m_data.end(), m_data.begin(),
					[&](uintmax_t x) -> uintmax_t { return div(x, value.value()); });
	return reduce();
}



auto irrpoly::gfpoly::multiply(const gfpoly &a, const gfpoly &b) -> gfpoly & {
	if (!a || !b) {
		return set_zero();
	}
	std::vector<uintmax_t> prod(a.size() + b.size() - 1, 0);
	for (uintmax_t i = 0; i < a.size(); ++i) {
		for (uintmax_t j = 0; j < b.size(); ++j) {
			prod[i + j] = add(prod[i + j], mul(a.m_data[i], b.m_data[j]));
		}
	}
	m_data.swap(prod);
	return reduce();
}



auto irrpoly::gfpoly::operator*=(const gfpoly &value) -> gfpoly & {
	return multiply(*this, value);
}



auto irrpoly::gfpoly::operator/=(const gfpoly &value) -> gfpoly & {
	return *this = quotient_remainder(*this, value).first;
}



auto irrpoly::gfpoly::operator%=(const gfpoly &value) -> gfpoly & {
	return *this = quotient_remainder(*this, value).second;
}
