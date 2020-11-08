/*!
 *	\file incremental_pca.hpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2020)
 */
#include "incremental_pca.hpp"
#include "jacobi/jacobi.hpp"
#include <algorithm>





TsTestsReturnCode init_incremental_pca(IncrementalPCAInfo *const pca_info, uint8_t const dim)
{
	try
	{
		pca_info->points_count = 0;
		pca_info->scaled_covariance.resize(dim);
		pca_info->eigenvectors.resize(dim);
		pca_info->eigenvalues.resize(dim, 0.0);
		for (uint8_t i = 0; i < dim; ++i)
		{
			pca_info->scaled_covariance[i].resize(dim, 0.0);
			pca_info->eigenvectors[i].resize(dim, 0.0);
		}
		pca_info->eigen_up_to_date = true;
	}
	catch(...)
	{
		return TSTESTS_RETURNCODE_FAIL_MEMORY;
	}
	
	return TSTESTS_RETURNCODE_SUCCESS;
}





void partial_fit(IncrementalPCAInfo *const pca_info, PCAMatrix &new_points)
{
	/*
	 * Iterate over all rows/points in (new_points)
	 */
	for (uint64_t point_i = 0; point_i < new_points.size(); ++point_i)
		for (uint8_t dim_i = 0; dim_i < pca_info->eigenvalues.size(); ++dim_i)
			/*
			 * Update covariance matrix
			 */
			for (uint8_t dim_j = dim_i; dim_j < pca_info->eigenvalues.size(); ++dim_j)
				pca_info->scaled_covariance[dim_i][dim_j] = pca_info->scaled_covariance[dim_j][dim_i] += (new_points[point_i][dim_i] - .5) * (new_points[point_i][dim_j] - .5);
	
	pca_info->points_count += new_points.size();
	pca_info->eigen_up_to_date = false;
	
	return;
}





void update_eigen(IncrementalPCAInfo *const pca_info)
{
	PCAMatrix covariance(pca_info->scaled_covariance);
	jacobi_public_domain::Jacobi<TSTESTS_COORDINATE_TYPE, PCAVector&, PCAMatrix&> eigen_calc(pca_info->eigenvalues.size());
	
	eigen_calc.Diagonalize(covariance, pca_info->eigenvalues, pca_info->eigenvectors);
	pca_info->eigen_up_to_date = true;
}

PCAMatrix get_principal_axes(IncrementalPCAInfo *const pca_info)
{
	if (!pca_info->eigen_up_to_date)
		update_eigen(pca_info);
	
	return pca_info->eigenvectors;
}

PCAVector get_components_variance(IncrementalPCAInfo *const pca_info)
{
	if (!pca_info->eigen_up_to_date)
		update_eigen(pca_info);
	
	return pca_info->eigenvalues;
}

PCAVector get_components_normalised_variance(IncrementalPCAInfo *const pca_info)
{
	PCAVector variance(get_components_variance(pca_info));
	TSTESTS_COORDINATE_TYPE sum = 0.0;
	
	std::for_each(variance.begin(), variance.end(), [&sum](TSTESTS_COORDINATE_TYPE elem){sum += elem;});
	std::transform(variance.begin(), variance.end(), variance.begin(), [&sum](TSTESTS_COORDINATE_TYPE elem){return elem / sum;});
	
	return variance;
}
