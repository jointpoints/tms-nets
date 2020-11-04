/*!
 *	\file raref.cpp
 *
 *	\author
 *		Arseny Zakharov (Russian Technological University, KMBO-01-17, Russia, 2020)
 */

#include "raref.hpp"
#include <algorithm>
#include <iostream>

typedef std::vector<size_t> Composition;
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

int main(int argc, char **argv)
{
	// RAREFMatrix matrix = {{1, 2, 3, 4, 5, 6},
	//                       {2, 3, 4, 5, 6, 7}};
	// print_matrix(matrix);
	// std::cout << std::endl;
	// swap_row(matrix, 0, 1);
	// print_matrix(matrix);
	// std::cout << std::endl;
	// subtract_row(matrix, 1, 0);
	// print_matrix(matrix);
	// std::cout << std::endl;

	// auto ans = generate_compositions(6, 5);
	// for (auto i : ans){
	// 	for (auto j : i){
	// 		std::cout << j << " ";
	// 	}
	// 	std::cout << std::endl;
	// }

	RAREFMatrix a = {{1, 1, 0, 0}, {1, 0, 1, 1}, {0, 1, 1, 0}};
	RAREFMatrix b = {{1, 1, 0, 0}, {1, 0, 1, 1}, {0, 1, 1, 1}};
	RAREFMatrix x = {{0, 0, 1, 1, 0, 0}, {1, 1, 1, 1, 1, 1}, {0, 1, 1, 0, 0, 1}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1}};
	RAREFMatrix y = {{1, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {1, 0, 1, 1, 1}, {1, 0, 0, 0, 0}, {1, 0, 0, 0, 1}};
	RAREFMatrix z = {{1, 1, 1, 0, 0}, {0, 1, 1, 1, 0}, {1, 0, 1, 1, 1}, {1, 0, 0, 0, 0}, {1, 0, 0, 0, 1}};

	// RAREF ar = compute_RAREF(a);
	// RAREF br = compute_RAREF(b);
	// RAREF cr = update_RAREF(a, b, ar);
	RAREF xr = compute_RAREF(x);
	RAREF yr = compute_RAREF(y);
	RAREF zr = compute_RAREF(z);
	RAREF wr = update_RAREF(y, z, yr);

	print_RAREF(xr, "x");
	print_RAREF(yr, "y");
	print_RAREF(zr, "z");
	print_RAREF(wr, "w");
}