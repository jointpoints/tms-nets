/*!
 *	\file raref.hpp
 *
 *	\author
 *		Arseny Zakharov (Russian Technological University, KMBO-01-17, Russia, 2020)
 *	\author
 *		Daria Sabirianova (Russian Technological University, KMBO-01-17, Russia, 2020)
 *	\author
 *		Sergey Kharlamov (Russian Technological University, KMBO-01-17, Russia, 2020)
 */

#ifndef _RAREF_HPP_
#define _RAREF_HPP_

#include "common.hpp"


TsTestsReturnCode find_defect(uint64_t &ro, uint64_t m, uint64_t s, std::function<std::vector<std::vector<uint>>(uint64_t const)> const &gamma_matrix_getter);

#endif