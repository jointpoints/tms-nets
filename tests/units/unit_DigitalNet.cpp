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

	SECTION("Check generating matrices")
	{
		REQUIRE( nondeg_net.get_generating_matrix(0) == generating_matrices[0] );
		REQUIRE( nondeg_net.get_generating_matrix(1) == generating_matrices[1] );
		REQUIRE( nondeg_net.get_generating_matrix(2) == generating_matrices[2] );
	}

	SECTION("Check direction numbers")
	{
		tms::DirNum direction_numbers = nondeg_net.get_direction_numbers(0);
		REQUIRE( direction_numbers[0] == 6 );
		REQUIRE( direction_numbers[1] == 1 );
		REQUIRE( direction_numbers[2] == 4 );
		direction_numbers = nondeg_net.get_direction_numbers(1);
		REQUIRE( direction_numbers[0] == 1 );
		REQUIRE( direction_numbers[1] == 2 );
		REQUIRE( direction_numbers[2] == 4 );
		direction_numbers = nondeg_net.get_direction_numbers(2);
		REQUIRE( direction_numbers[0] == 6 );
		REQUIRE( direction_numbers[1] == 5 );
		REQUIRE( direction_numbers[2] == 7 );
	}

	SECTION("Construction of 3-dimensional point (classical enumeration of points)")
	{
		// classical representation of index 4 is (0 0 1)^T
		tms::Point point = nondeg_net.generate_point_classical(4);
		REQUIRE( point[0] == Catch::Approx(0.5).margin(0.0001) );
		REQUIRE( point[1] == Catch::Approx(0.5).margin(0.0001) );
		REQUIRE( point[2] == Catch::Approx(0.875).margin(0.0001) );
	}

	SECTION("Construction of 3-dimensional point (Gray code enumeration of points)")
	{
		// Gray representation of index 4 is (0 1 1)^T
		tms::Point point = nondeg_net.generate_point(4);
		REQUIRE( point[0] == Catch::Approx(0.625).margin(0.0001) );
		REQUIRE( point[1] == Catch::Approx(0.75).margin(0.0001) );
		REQUIRE( point[2] == Catch::Approx(0.25).margin(0.0001) );
	}

	SECTION("Construction of 3-dimensional point (Gray code enumeration of points, integer point)")
	{
		// Gray representation of index 4 is (0 1 1)^T
		tms::IntPoint point = nondeg_net.generate_point_int(4);
		REQUIRE( point[0] == 5 );
		REQUIRE( point[1] == 6 );
		REQUIRE( point[2] == 2 );
	}
}
