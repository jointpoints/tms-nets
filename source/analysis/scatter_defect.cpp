/**
 * \file
 *       scatter_defect.cpp
 *
 * \author
 *       Andrei Eliseev (JointPoints), 2021
 */
#include "../../include/tms-nets/analysis/analysis.hpp"
#include "../../include/tms-nets/thirdparty/jacobi/jacobi.hpp"





using PCAMatrixRow = std::vector<tms::Real>;
using PCAMatrix    = std::vector<PCAMatrixRow>;





tms::Point tms::analysis::scatter_defect(tms::DigitalNet const &net)
{
	CountInt        point_count     = (1ULL << net.m());
	BasicInt        s               = net.s();
	Real            scatter         = 0;
	Point           point;
	PCAMatrix       cov_matrix;
	PCAMatrix       dummy;

	Point           result(s, 0.0);

	// 1. Prepare container of covariance matrix
	cov_matrix.resize(s);
	for (auto &row : cov_matrix)
		row.resize(s, 0);
	dummy.resize(s);
	for (auto &row : dummy)
		row.resize(s, 0);

	// 2. Update covariance matrix and scatter
	for (CountInt point_i = 0; point_i < point_count; ++point_i)
	{
		point = net.generate_point(point_i);

		for (BasicInt i = 0; i < s; ++i)
			for (BasicInt j = i; j < s; ++j)
				cov_matrix[i][j] = cov_matrix[j][i] += (point[i] - 0.5) * (point[j] - 0.5);
	}

	// 3. Retrieve eigenvalues of covariance matrix
	jacobi_public_domain::Jacobi<Real, PCAMatrixRow &, PCAMatrix &> diagonaliser(s);
	diagonaliser.Diagonalize(cov_matrix, result, dummy);

	// 4. Calculate defect
	std::for_each(result.begin(), result.end(), [&scatter](Real value){scatter += value * value;});
	std::transform(result.begin(), result.end(), result.begin(), [s, scatter, point_count](Real value){return value * value / scatter - 1.0 / s;});

	return result;
}
