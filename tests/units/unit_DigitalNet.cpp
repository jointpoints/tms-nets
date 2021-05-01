/**
 * \file
 *       unit_DigitalNet.cpp
 *
 * \author
 *       Andrei Eliseev (JointPoints), 2021
 */
#include "../catch2/catch_amalgamated.hpp"
#include "../../include/tms-nets.hpp"





TEST_CASE("Validation of DigitalNet class", "[nets][DigitalNet]")
{
	std::vector<tms::GenMat> generating_matrices(3);
	REQUIRE_NOTHROW( generating_matrices[0] = tms::GenMat{tms::GenMatRow{1,0,1}, tms::GenMatRow{1,0,0}, tms::GenMatRow{0,1,0}} );
	REQUIRE_NOTHROW( generating_matrices[1] = tms::GenMat{tms::GenMatRow{0,0,1}, tms::GenMatRow{0,1,0}, tms::GenMatRow{1,0,0}} );
	REQUIRE_NOTHROW( generating_matrices[2] = tms::GenMat{tms::GenMatRow{1,1,1}, tms::GenMatRow{1,0,1}, tms::GenMatRow{0,1,1}} );
	tms::DigitalNet nondeg_net(generating_matrices);

	SECTION("Check parameter m")
	{
		REQUIRE( nondeg_net.get_m() == 3 );
	}

	SECTION("Check parameter s")
	{
		REQUIRE( nondeg_net.get_s() == 3 );
	}

	SECTION("Construction of 3-dimensional point")
	{
		tms::Point point = nondeg_net.generate_point_classical(5);
		REQUIRE( point[0] == Catch::Approx(0.25).margin(0.0001) );
		REQUIRE( point[1] == Catch::Approx(0.625).margin(0.0001) );
		REQUIRE( point[2] == Catch::Approx(0.125).margin(0.0001) );
	}
}
