/**
 * \file
 *       unit_Niederreiter.cpp
 *
 * \author
 *       Andrei Eliseev (JointPoints), 2021
 */
#include "../catch2/catch_amalgamated.hpp"
#include "../../include/tms-nets.hpp"





TEST_CASE("Validation of Niederreiter class", "[nets][Niederreiter]")
{
	std::vector<tms::GenMat> generating_matrices(3);
	REQUIRE_NOTHROW( generating_matrices[0] = tms::GenMat{tms::GenMatRow{1,0,0}, tms::GenMatRow{0,1,0}, tms::GenMatRow{0,0,1}} );
	REQUIRE_NOTHROW( generating_matrices[1] = tms::GenMat{tms::GenMatRow{1,1,1}, tms::GenMatRow{0,1,0}, tms::GenMatRow{0,0,1}} );
	REQUIRE_NOTHROW( generating_matrices[2] = tms::GenMat{tms::GenMatRow{0,1,1}, tms::GenMatRow{1,1,0}, tms::GenMatRow{0,0,1}} );
	tms::Niederreiter nondeg_net(3, 3);

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
		REQUIRE( direction_numbers[0] == 4 );
		REQUIRE( direction_numbers[1] == 2 );
		REQUIRE( direction_numbers[2] == 1 );
		direction_numbers = nondeg_net.get_direction_numbers(1);
		REQUIRE( direction_numbers[0] == 4 );
		REQUIRE( direction_numbers[1] == 6 );
		REQUIRE( direction_numbers[2] == 5 );
		direction_numbers = nondeg_net.get_direction_numbers(2);
		REQUIRE( direction_numbers[0] == 2 );
		REQUIRE( direction_numbers[1] == 6 );
		REQUIRE( direction_numbers[2] == 5 );
	}

	SECTION("Construction of 3-dimensional point (classical enumeration of points)")
	{
		// classical representation of index 4 is (0 0 1)^T
		tms::Point point = nondeg_net.generate_point_classical(4);
		REQUIRE( point[0] == Catch::Approx(0.125).margin(0.0001) );
		REQUIRE( point[1] == Catch::Approx(0.625).margin(0.0001) );
		REQUIRE( point[2] == Catch::Approx(0.625).margin(0.0001) );
	}

	SECTION("Construction of 3-dimensional point (Gray code enumeration of points)")
	{
		// Gray representation of index 4 is (0 1 1)^T
		tms::Point point = nondeg_net.generate_point(4);
		REQUIRE( point[0] == Catch::Approx(0.375).margin(0.0001) );
		REQUIRE( point[1] == Catch::Approx(0.375).margin(0.0001) );
		REQUIRE( point[2] == Catch::Approx(0.375).margin(0.0001) );
	}

	SECTION("Construction of 3-dimensional point (Gray code enumeration of points, integer point)")
	{
		// Gray representation of index 4 is (0 1 1)^T
		tms::IntPoint point = nondeg_net.generate_point_int(4);
		REQUIRE( point[0] == 3 );
		REQUIRE( point[1] == 3 );
		REQUIRE( point[2] == 3 );
	}

	SECTION("Estimation of sequence defect")
	{
		REQUIRE( nondeg_net.get_t_estimate() == 1 );
	}
}
