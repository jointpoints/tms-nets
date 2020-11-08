/*!
 *	\file incremental_pca.hpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2020)
 */
#ifndef _INCREMENTAL_PCA_HPP_
#define _INCREMENTAL_PCA_HPP_



#include "common.hpp"



typedef std::vector<TSTESTS_COORDINATE_TYPE>    PCAVector;
typedef std::vector<PCAVector>	                PCAMatrix;

typedef struct IncrementalPCAInfo
{
	uint64_t                        points_count;
	PCAMatrix                       scaled_covariance;
	PCAMatrix                       eigenvectors;
	PCAVector                       eigenvalues;
	bool                            eigen_up_to_date;
}
IncrementalPCAInfo;



TsTestsReturnCode init_incremental_pca              (IncrementalPCAInfo *const pca_info, uint8_t const dim);
void              partial_fit                       (IncrementalPCAInfo *const pca_info, PCAMatrix &new_points);
PCAMatrix         get_principal_axes                (IncrementalPCAInfo *const pca_info);
PCAVector         get_components_variance           (IncrementalPCAInfo *const pca_info);
PCAVector         get_components_normalised_variance(IncrementalPCAInfo *const pca_info);



#endif // _INCREMENTAL_PCA_HPP_
