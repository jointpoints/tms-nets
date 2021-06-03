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
		REQUIRE( nondeg_net.m() == 3 );
	}

	SECTION("Check parameter s")
	{
		REQUIRE( nondeg_net.s() == 3 );
	}

	SECTION("Check generating matrices")
	{
		CHECK( nondeg_net.generating_matrix(0) == generating_matrices[0] );
		CHECK( nondeg_net.generating_matrix(1) == generating_matrices[1] );
		CHECK( nondeg_net.generating_matrix(2) == generating_matrices[2] );
	}

	SECTION("Check generating numbers")
	{
		tms::GenNum generating_numbers = nondeg_net.generating_numbers(0);
		CHECK( generating_numbers[0] == 6 );
		CHECK( generating_numbers[1] == 1 );
		CHECK( generating_numbers[2] == 4 );
		generating_numbers = nondeg_net.generating_numbers(1);
		CHECK( generating_numbers[0] == 1 );
		CHECK( generating_numbers[1] == 2 );
		CHECK( generating_numbers[2] == 4 );
		generating_numbers = nondeg_net.generating_numbers(2);
		CHECK( generating_numbers[0] == 6 );
		CHECK( generating_numbers[1] == 5 );
		CHECK( generating_numbers[2] == 7 );
	}

	SECTION("Construction of 3-dimensional point (classical enumeration of points)")
	{
		// classical representation of index 4 is (0 0 1)^T
		tms::Point point = nondeg_net.generate_point_classical(4);
		CHECK( point[0] == Catch::Approx(0.5).margin(0.0001) );
		CHECK( point[1] == Catch::Approx(0.5).margin(0.0001) );
		CHECK( point[2] == Catch::Approx(0.875).margin(0.0001) );
	}

	SECTION("Construction of 3-dimensional point (Gray code enumeration of points)")
	{
		// Gray representation of index 4 is (0 1 1)^T
		tms::Point point = nondeg_net.generate_point(4);
		CHECK( point[0] == Catch::Approx(0.625).margin(0.0001) );
		CHECK( point[1] == Catch::Approx(0.75).margin(0.0001) );
		CHECK( point[2] == Catch::Approx(0.25).margin(0.0001) );
	}

	SECTION("Construction of 3-dimensional point (Gray code enumeration of points, integer point)")
	{
		// Gray representation of index 4 is (0 1 1)^T
		tms::IntPoint point = nondeg_net.generate_int_point(4);
		CHECK( point[0] == 5 );
		CHECK( point[1] == 6 );
		CHECK( point[2] == 2 );
	}
}
