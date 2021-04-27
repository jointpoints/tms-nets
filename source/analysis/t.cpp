/**
 * @file    t.cpp
 * 
 * @brief   Contains functions for calculations of t parameter.
 *
 * @author  Arseny Zakharov (Russian Technological University, KMBO-01-17, Russia, 2020)
 *          Daria Sabirianova (Russian Technological University, KMBO-01-17, Russia, 2020)
 *          Sergey Kharlamov (Russian Technological University, KMBO-01-17, Russia, 2020)
 *          Andrei Eliseev (JointPoints), 2021
 * 
 */
#include "../../include/tms-nets/analysis/analysis.hpp"

#include <stdexcept>    // needed for exceptions
#include <algorithm>    // needed for "equal"





// Auxiliary content





typedef std::vector<bool>           RAREFVector;
typedef std::vector<RAREFVector>    RAREFMatrix;

typedef struct RAREF
{
	RAREFMatrix L;
	RAREFMatrix T;
	std::vector<size_t> p;
} RAREF;

typedef std::vector<size_t>      Composition;
typedef std::vector<Composition> Compositions;

// Debug
/*void print_matrix(RAREFMatrix const &matrix)
{
	for (auto i : matrix)
	{
		for (auto j : i)
		{
			std::cout << j << " ";
		}
		std::cout << std::endl;
	}
}

void print_RAREF(RAREF const &r, std::string const &name)
{
	std::cout << "RAREF for " << name << ":" << std::endl;
	std::cout << "L:" << std::endl;
	print_matrix(r.L);
	std::cout << "T:" << std::endl;
	print_matrix(r.T);
	std::cout << "p:" << std::endl;
	for (auto i : r.p)
	{
		std::cout << i << " ";
	}
	std::cout << std::endl
			  << std::endl;
}*/
// Debug

/*
 * Change row to absolute value of subtaction with another row
 */
inline void subtract_row(RAREFMatrix &matrix, size_t min_ind, size_t sub_ind)
{
	std::transform(matrix[min_ind].begin(), matrix[min_ind].end(),
				   matrix[sub_ind].begin(), matrix[min_ind].begin(),
				   [](bool a, bool b) { return a != b; });
}

/*
 * Swap rows in matrix
 */
inline void swap_row(RAREFMatrix &matrix, size_t row1_ind, size_t row2_ind)
{
	std::swap(matrix[row1_ind], matrix[row2_ind]);
}

RAREFMatrix ident_matrix(size_t size)
{
	RAREFMatrix res(size, RAREFVector(size, 0));
	for (size_t i = 0; i < size; i++)
	{
		res[i][i] = 1;
	}
	return res;
}

void composition_recursive(Composition compos, size_t curr_ind, Compositions &res)
{
	if (curr_ind >= compos.size() - 1)
	{
		return;
	}
	composition_recursive(compos, curr_ind + 1, res);
	if (compos[curr_ind] > 1)
	{
		compos[curr_ind]--;
		compos[curr_ind + 1]++;
		res.push_back(compos);
		composition_recursive(compos, curr_ind, res);
	}
}

/*
 * Generate compositions of q of length u
 */
Compositions generate_compositions(size_t q, size_t u)
{
	Compositions result;
	Composition compos(u, 1);
	compos[0] = q - u + 1;
	result.push_back(compos);
	composition_recursive(compos, 0, result);
	return result;
}

/*
 * Compute binomial coefficient
 */
size_t binomial_coefficient(size_t n, size_t k)
{
	size_t res = 1;

	if (k > n - k)
	{
		k = n - k;
	}

	for (size_t i = 0; i < k; i++)
	{
		res *= (n - i);
		res /= (i + 1);
	}

	return res; 
}

/*
 * Generate compositions from [from:to] of size k
 */
Compositions binomial_coefficient(size_t from, size_t to, size_t k)
{
	size_t n = to - from + 1;
	std::vector<bool> sep(n, 0);
	std::fill(sep.begin(), sep.begin() + k, 1);

	Compositions res;

	do
	{
		Composition compos;
		for (size_t i = 0; i < n; i++)
		{
			if (sep[i])
			{
				compos.push_back(i + from);
			}
		}
		res.push_back(compos);
	} while (std::prev_permutation(sep.begin(), sep.end()));
	return res;
}

/*
 * Compute RAREF for matrix C
 */
RAREF compute_RAREF(RAREFMatrix const &C, RAREFMatrix const &L)
{
	size_t q = C.size();
	if (!q)
	{
		throw std::logic_error("\nWrong C dimensions");
	}
	size_t k = C[0].size();
	if (q > k)
	{
		throw std::logic_error("\nWrong C dimensions");
	}

	RAREF res;
	res.L = L;
	res.T = C;
	res.p = {};

	for (size_t i = 0; i < q; i++)
	{
		size_t j = 0;
		for (; j < k; j++)
		{
			if (res.T[i][j] == 0)
			{
				for (size_t i2 = i + 1; i2 < q; i2++)
				{
					if (res.T[i2][j] == 1)
					{
						swap_row(res.T, i, i2);
						swap_row(res.L, i, i2);
						break;
					}
				}
			}
			if (res.T[i][j] == 1)
			{
				res.p.push_back(j);
				break;
			}
		}
		j = j == k ? k - 1 : j;
		for (size_t i2 = 0; i2 < q; i2++)
		{
			if (i2 != i && (res.T[i2][j] == 1))
			{
				subtract_row(res.T, i2, i);
				subtract_row(res.L, i2, i);
			}
		}
	}
	return res;
}

/*
 * Call compute_RAREF with identity matrix as L
 */
RAREF compute_RAREF(RAREFMatrix const &C)
{
	size_t q = C.size();
	return compute_RAREF(C, ident_matrix(q));
}

/*
 * Update calculated RAREF for matrix with difference in only one row
 */
RAREF update_RAREF(RAREFMatrix const &C, RAREFMatrix const &C2, RAREF const &src)
{
	size_t q = C.size();
	if (!q)
	{
		throw std::logic_error("\nWrong C dimensions");
	}
	size_t k = C[0].size();
	if (q > k)
	{
		throw std::logic_error("\nWrong C dimensions");
	}
	if (C2.size() != q || C2[0].size() != k)
	{
		throw std::logic_error("\nWrong C2 dimensions");
	}

	RAREF res;
	res.L = src.L;
	res.T = src.T;
	res.p = src.p;

	if (C == C2)
	{
		return res;
	}

	size_t diff_row = 0;
	for (size_t i = 0; i < q; i++)
	{
		for (size_t j = 0; j < k; j++)
		{
			if (C[i][j] != C2[i][j])
			{
				diff_row = i;
			}
		}
	}

	for (size_t i = 0; i < q; i++)
	{
		if (res.L[i][diff_row] == 1)
		{
			swap_row(res.T, i, diff_row);
			swap_row(res.L, i, diff_row);
			break;
		}
	}

	for (size_t i = 0; i < q; i++)
	{
		if (i != diff_row && res.L[i][diff_row] == 1)
		{
			subtract_row(res.T, i, diff_row);
			subtract_row(res.L, i, diff_row);
			for (size_t j = 0; j < k; j++)
			{
				if (res.T[diff_row][j] == 1)
				{
					std::remove(res.p.begin(), res.p.end(), j);
				}
			}
		}
	}

	for (size_t i = 0; i < q; i++)
	{
		if (i != diff_row)
		{
			res.L[diff_row][i] = 0;
		}
	}

	res.T[diff_row] = C2[diff_row];

	for (size_t i = 0; i < res.p.size(); i++)
	{
		size_t pivot = res.p[i];
		size_t I = 0;
		for (size_t j = 0; j < q; j++)
		{
			if (res.T[j][pivot] == 1 && j != diff_row)
			{
				I = j;
			}
		}
		if (res.T[diff_row][pivot] == 1)
		{
			subtract_row(res.T, diff_row, I);
			subtract_row(res.L, diff_row, I);
		}
	}
	return compute_RAREF(res.T, res.L);
}

/*
 * Add rows from [0:to] of src to the end of dst
 * If <reverse> == true, rows will be added in reverse order
 */
void add_rows(RAREFMatrix &dst, RAREFMatrix const &src, size_t const to, bool const reverse = false)
{
	size_t n = to + 1;
	for (size_t i = 0; i < n; i++)
	{
		size_t ind = (reverse) ? (to - i) : (i);
		dst.push_back(src[ind]);
	}
	return;
}

Composition generate_projections(size_t s, size_t u, size_t k,
                                 std::vector<RAREFMatrix> const &gen_mat, Composition const &rho)
{
	size_t qmax = u == 2 ? k : *std::min_element(rho.begin(), rho.end());
	Composition res;
	Compositions c = binomial_coefficient(0, s - 1, u);
	for (size_t i = 0; i < c.size(); i++)
	{
		size_t ro_tilda;
		size_t curr_rank = u - 1;
		for (size_t q = qmax; q >= u; q--)
		{
			bool                flag                    = false;
			Compositions        comp_arr                = generate_compositions(q, u);
			RAREFMatrix         matrix                  = {};
			std::vector<bool>   is_section_reversed(comp_arr[0].size(), false);
			RAREF               r;
			for (size_t j = 0; j < comp_arr.size(); j++)
			{
				size_t matrix_rank;
				if (j >= 1)
				{
					// Try to update RAREF faster
					Composition abs_diff(comp_arr[0].size());
					std::transform(comp_arr[j].begin(), comp_arr[j].end(), comp_arr[j - 1].begin(),
					               abs_diff.begin(), [](size_t a, size_t b){return std::max(a, b) - std::min(a, b);});
					// We study two consequent compositions. Faster computation is applicable, if:
					//     1. a new composite matrix <matrix> differs from the previous by eactly one row,
					//        which basically means that there is a pair of two adjacent sections in the previous matrix
					//        such that if the "upper section" "loses" one row and the "lower section" "adds" one new row,
					//        the previous matrix becomes a new matrix;
					//     2. the "upper section" is packed in direct order and the "lower section" is packed in reversed
					//        order.
					bool    fast_update_available       = false;
					size_t  fast_update_lower_section   = 0;
					for (size_t section_i = 0; section_i < comp_arr[0].size() - 1; ++section_i)
					{
						if (std::abs(static_cast<int>(comp_arr[j][section_i] - comp_arr[j - 1][section_i])) > 1)
						{
							fast_update_available = false;
							break;
						}
						if ((std::abs(static_cast<int>(comp_arr[j][section_i]     - comp_arr[j - 1][section_i]))     == 1)  &&
						    (std::abs(static_cast<int>(comp_arr[j][section_i + 1] - comp_arr[j - 1][section_i + 1])) == 1)  &&
						    ((!is_section_reversed[section_i] && is_section_reversed[section_i + 1]) || (!is_section_reversed[section_i] && !is_section_reversed[section_i + 1] && (comp_arr[j - 1][section_i + 1] == 1))))
						{
							if (!fast_update_available)
							{
								fast_update_available = true;
								fast_update_lower_section = section_i + 1;
							}
							else
							{
								fast_update_available = false;
								break;
							}
						}
					}
					if (std::abs(static_cast<int>(comp_arr[j][comp_arr[0].size() - 1] - comp_arr[j - 1][comp_arr[0].size() - 1])) > 1)
						fast_update_available = false;
					// If fast update is available, update faster
					if (fast_update_available)
					{
						RAREFMatrix old_matrix = matrix;
						matrix = {};

						for (size_t P = 0; P < comp_arr[0].size(); P++)
						{
							if (P != fast_update_lower_section)
							{
								add_rows(matrix, gen_mat[c[i][P]], comp_arr[j][P] - 1);
								is_section_reversed[P] = false;
							}
							else
							{
								add_rows(matrix, gen_mat[c[i][P]], comp_arr[j][P] - 1, true);
								is_section_reversed[P] = true;
							}
						}
						r = update_RAREF(old_matrix, matrix, r);
						matrix_rank = r.p.size();
					}
					// Otherwise, update slowly :(
					else
					{
						matrix = {};
						for (size_t P = 0; P < comp_arr[0].size(); P++)
						{
							add_rows(matrix, gen_mat[c[i][P]], comp_arr[j][P] - 1);
							is_section_reversed[P] = false;
						}
						r = compute_RAREF(matrix);
						matrix_rank = r.p.size();
					}
				}
				else
				{
					for (size_t P = 0; P < comp_arr[0].size(); P++)
					{
						add_rows(matrix, gen_mat[c[i][P]], comp_arr[j][P] - 1);
						is_section_reversed[P] = false;
					}
					r = compute_RAREF(matrix);
					matrix_rank = r.p.size();
				}
				if (matrix_rank < q)
				{
					flag = true;
					break;
				}
				else
				{
					curr_rank = q;
				}
			}
			if (!flag)
			{
				break;
			}
		}
		ro_tilda = curr_rank;
		res.push_back(std::min(ro_tilda, qmax));
	}
	return res;
}

/*
 * Find defect of (t,m,s)-net
 */
Composition find_rho_inner(size_t k, size_t s, size_t dmax,
                         std::vector<RAREFMatrix> const &gen_mat)
{
	Composition rho;
	if (s == 1)
	{
		rho.push_back(k);
		return rho;
	}
	for (size_t u = 2; u <= dmax; u++)
	{
		rho = generate_projections(s, u, k, gen_mat, rho);
	}
	return rho;
}

/*
 * Cast matrix of uints to matrix of bools
 */
RAREFMatrix cast_matrix(tms::GenMat const &src)
{
	RAREFMatrix dst;
	for (tms::BasicInt i = 0; i < src.size(); ++i)
	{
		RAREFVector row;
		for (tms::BasicInt j = 0; j < src.size(); ++j)
			row.push_back(src[i][j]);
		dst.push_back(row);
	}
	return dst;
}





// Main function





tms::BasicInt tms::analysis::t(DigitalNet const &net)
{
	std::vector<RAREFMatrix> genMat;
	for (BasicInt dim_i = 0; dim_i < net.get_s(); dim_i++)
	{
		genMat.push_back(cast_matrix(net.get_generating_matrix(dim_i)));
	}
	Composition rho = find_rho_inner(net.get_m(), net.get_s(), net.get_s(), genMat);
	if (rho.size() != 1)
		throw std::runtime_error("Was not able to calculate t.");
	return net.get_m() - rho[0];
}
