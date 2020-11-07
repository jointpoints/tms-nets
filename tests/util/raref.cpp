/*!
 *	\file raref.cpp
 *
 *	\author
 *		Arseny Zakharov (Russian Technological University, KMBO-01-17, Russia, 2020)
 *	\author
 *		Daria Sabirianova (Russian Technological University, KMBO-01-17, Russia, 2020)
 *	\author
 *		Sergey Kharlamov (Russian Technological University, KMBO-01-17, Russia, 2020)
 */

#include "raref.hpp"
#include <algorithm>
#include <iostream>



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

/*
 * Change row to absolute value of subtaction with another row
 */
inline void subtract_row(RAREFMatrix &matrix, size_t min_ind, size_t sub_ind)
{
	std::transform(matrix[min_ind].begin(), matrix[min_ind].end(),
				   matrix[sub_ind].begin(), matrix[min_ind].begin(),
				   [](bool a, bool b) { return !a * b + !b * a; });
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
 * Generate compositions
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

void add_rows(RAREFMatrix &dst, RAREFMatrix const &src, size_t from, size_t to)
{
	bool forward = from <= to ? true : false;
	size_t n = std::max(from, to) - std::min(from, to) + 1;
	for (size_t i = 0; i < n; i++)
	{
		size_t ind = forward ? from + i : from - i;
		dst.push_back(src[ind]);
	}
	// print_matrix(dst);
	// std::cout << "--" << std::endl;
}

Composition generate_projections(size_t s, size_t u, size_t k,
                                 std::vector<RAREFMatrix> const &gen_mat, Composition const &ro)
{
	size_t qmax = u == 2 ? k : *std::min_element(ro.begin(), ro.end());
	Composition res;
	Compositions c = binomial_coefficient(0, s - 1, u);
	for (size_t i = 0; i < c.size(); i++)
	{
		size_t ro_tilda;
		size_t curr_rank = u - 1;
		for (size_t q = qmax; q >= u; q--)
		{
			bool flag = 0;
			Compositions comp_arr = generate_compositions(q, u);
			RAREFMatrix matrix = {};
			RAREF r;
			for (size_t j = 0; j < comp_arr.size(); j++)
			{
				size_t matrix_rank;
				if (j >= 1)
				{
					Composition abs_diff(comp_arr[0].size());
					std::transform(comp_arr[j].begin(), comp_arr[j].end(), comp_arr[j - 1].begin(),
					               abs_diff.begin(), [](size_t a, size_t b){return std::max(a, b) - std::min(a, b);});
					// std::cout << "transformed " << j << std::endl;
					if (*std::max_element(abs_diff.begin(), abs_diff.end()) == 1)
					{
						// std::cout << "if " << abs_diff.size() << std::endl;
						RAREFMatrix old_matrix = matrix;
						matrix = {};
						Composition abs_diff_ones;
						for (size_t ind = 0; ind < abs_diff.size(); ind++)
						{
							if (abs_diff[ind] == 1)
							{
								abs_diff_ones.push_back(ind);
							}
						}
						size_t coeff = abs_diff_ones[1];
						// std::cout << "here" << std::endl;

						for (size_t P = 0; P < comp_arr[0].size(); P++)
						{
							if (P != coeff)
							{
								add_rows(matrix, gen_mat[c[i][P]], 0, comp_arr[j][P] - 1);
							}
							else
							{
								add_rows(matrix, gen_mat[c[i][P]], comp_arr[j][P] - 1, 0);
							}
						}
						// std::cout << "UPD" << std::endl;
						// print_matrix(old_matrix);
						// print_matrix(matrix);
						// print_RAREF(r, "r");
						r = update_RAREF(old_matrix, matrix, r);
						// std::cout << "updated" << std::endl;
						matrix_rank = r.p.size();
					}
					else
					{
						// std::cout << "else" << std::endl;
						matrix = {};
						for (size_t P = 0; P < comp_arr[0].size(); P++)
						{
							add_rows(matrix, gen_mat[c[i][P]], 0, comp_arr[j][P] - 1);
						}
						// std::cout << "CMP" << std::endl;
						r = compute_RAREF(matrix);
						matrix_rank = r.p.size();
					}
				}
				else
				{
					for (size_t P = 0; P < comp_arr[0].size(); P++)
					{
						add_rows(matrix, gen_mat[c[i][P]], 0, comp_arr[j][P] - 1);
					}
					// std::cout << "CMP2 " << j << std::endl;
					r = compute_RAREF(matrix);
					// std::cout << "computed" << std::endl;
					matrix_rank = r.p.size();
				}
				if (matrix_rank < q)
				{
					flag = 1;
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

Composition find_defect_inner(size_t k, size_t s, size_t dmax,
                         std::vector<RAREFMatrix> const &gen_mat)
{
	Composition ro;
	if (s == 1)
	{
		ro.push_back(k);
		return ro;
	}
	for (size_t u = 2; u <= dmax; u++)
	{
		ro = generate_projections(s, u, k, gen_mat, ro);
	}
	return ro;
}

void print_matrix(RAREFMatrix const &matrix)
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
}

RAREFMatrix cast_matrix(std::vector<std::vector<uint>> const &src)
{
	RAREFMatrix dst;
	for (auto i : src)
	{
		RAREFVector row;
		for (auto j : i)
		{
			row.push_back(j);
		}
		dst.push_back(row);
	}
	return dst;
}

TsTestsReturnCode find_defect(uint &ro, uint m, uint s, std::function<std::vector<std::vector<uint>>(uint const)> const &gamma_matrix_getter)
{
	try
	{
		std::vector<RAREFMatrix> genMat;
		for (uint i = 0; i < s; i++)
		{
			genMat.push_back(cast_matrix(gamma_matrix_getter(i)));
		}
		Composition defect = find_defect_inner(m, s, s, genMat);
		if (defect.size() != 1)
		{
			return TSTESTS_RETURNCODE_FAIL_GENERAL;
		}
		ro = defect[0];
	}
	catch(const std::bad_alloc &e)
	{
		return TSTESTS_RETURNCODE_FAIL_MEMORY;
	}
	catch(const std::logic_error &e)
	{
		return TSTESTS_RETURNCODE_FAIL_INPUT;
	}
	return TSTESTS_RETURNCODE_SUCCESS;
}