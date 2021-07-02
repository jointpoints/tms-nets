/**
 * @file    analysis.hpp
 * 
 * @brief   Contains functions for digital nets analyses, like figures of merit and various statistics.
 *
 * @author  Andrei Eliseev (JointPoints), 2021
 * 
 */
#ifndef TMS_NETS_ANALYSIS_HPP
#define TMS_NETS_ANALYSIS_HPP

#include "../digital_net.hpp"





/**
 * @namespace tms::analysis
 * 
 * @brief Contains functions for digital nets analyses, like figures of merit and various statistics.
 * 
 * The precision of Quasi Monte Carlo calculations is directly related to the discrepancy of the
 * digital net (e.g., search for "Koksma–Hlawka inequality"). Discrepancy itself is computationally
 * hard to be calculated directly which makes researchers introduce different metrics conventionally
 * called <b>figures of merit</b> that help one estimate the true discrepancy.
 * 
 * Statistics features that are also introduced here might be useful to verify some basic
 * properties that any "good" digital net should obtain.
 */
namespace tms::analysis
{
	


	/// @name Figures of merit
	/// @{

	/**
	 * Calculates the precise value of \f$t\f$
	 * 
	 * Finds out the least value of \f$t\f$ for which the given digital net is a \f$(t, m, s)\f$-net
	 * in base \f$2\f$ (we remind that values of \f$m\f$ and \f$s\f$ are known at the moment of net
	 * construction).
	 * 
	 * @param   net     A digital net.
	 * 
	 * @returns Precise value of parameter \f$t\f$ for the given digital net.
	 * 
	 * @throws  invalid_argument    If at least one generating matrix of the given net is degenerate
	 *                              (determinant equals zero over \f$\mathbb{F}_2\f$).
	 * 
	 * @note Exceptions will never be encountered in this function if it used with \c Niederreiter
	 * and \c Sobol nets provided in this library. Their generating matrices are \b guaranteed to
	 * always be non-degenerate. If you discover a case when this function fails to perform on these
	 * nets, [<b>report a bug</b>](https://github.com/jointpoints/tms-nets/issues/new).
	 * 
	 * @paragraph References
	 * 1. Marion P., Godin M., L'Ecuyer P. (2020) An algorithm to compute the t-value of a digital
	 * net and of its projections. *Journal of Computational and Applied Mathematics*. **371**. 112669.
	 * https://doi.org/10.1016/j.cam.2019.112669
	 */
	BasicInt            t               (DigitalNet const &net);

	///@}



	/// @name Statistics
	/// @{
	
	/**
	 * Calculates the scatter defect
	 * 
	 * Performs a special version of model-based principal component analysis that was specially
	 * optimised for digital \f$(t, m, s)\f$-nets. Here, having a singular value decomposition
	 * \f$C = U \Sigma V^*\f$ for a covariance matrix \f$C\f$, a \b scatter is defined as a sum of
	 * squares of all singular values, <b>relative influence</b> of the \f$i\f$-th principal axis
	 * is defined as a ratio of the squared \f$i\f$-th singular value to the scatter, and a
	 * <b>scatter defect</b> of the \f$i\f$-th principal axis is defined as a difference between
	 * its relative influence and the value of \f$s^{-1}\f$.
	 * 
	 * If points of the given digital net saturate the \f$s\f$-dimensional unit cube perfectly
	 * equally in all directions, the scatter defect will be \f$\overrightarrow{0} \in
	 * \mathbb{R}^s\f$. If points are scatterred along the \f$i\f$-th principal axis more than it
	 * is expected under the perfectly equal saturation, then the \f$i\f$-th component of defect
	 * vector will be positive, otherwise, negative.
	 * 
	 * @param   net     A digital net.
	 * 
	 * @returns A \ref tms::Point the \f$i\f$-th component of which equals the scatter defect along
	 * the \f$i\f$-th principal axis.
	 * 
	 * @note Principal axes always form an orthonormal basis in the \f$s\f$-dimensional space,
	 * however, there are <i>no guarantees</i> that they will match with the basis that is used to
	 * express the coordinates of digital net points.
	 * 
	 * @paragraph References
	 * 1. Eliseev A. (2020) Properties of (t, s)-sequences and their verification. Bachelor's
	 * graduation work (in Russian).
	 */
	tms::Point          scatter_defect  (DigitalNet const &net);

	/// @}



}; // namespace tms::analysis





#endif // #ifndef TMS_NETS_ANALYSIS_HPP
