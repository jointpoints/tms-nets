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
	 * @returns Precise value of parameter \f$t\f$ the given digital net.
	 */
	BasicInt t(DigitalNet const &net);

	///@}



	/// @name Statistics
	/// @{
	
	/**
	 * 
	 */
	//...

	/// @}



}; // namespace tms::analysis





#endif // #ifndef TMS_NETS_ANALYSIS_HPP
