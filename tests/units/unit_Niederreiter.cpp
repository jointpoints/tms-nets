/**
 * \file
 *       unit_Niederreiter.cpp
 *
 * \author
 *       Andrei Eliseev (JointPoints), 2021
 */
#include "../catch2/catch_amalgamated.hpp"
#include "../../include/tms-nets.hpp"





TEST_CASE("Validation of Niederreiter class, (m, s) constructor", "[nets][Niederreiter]")
{
	std::vector<tms::GenMat> generating_matrices(3);
	REQUIRE_NOTHROW( generating_matrices[0] = tms::GenMat{tms::GenMatRow{1,0,0}, tms::GenMatRow{0,1,0}, tms::GenMatRow{0,0,1}} );
	REQUIRE_NOTHROW( generating_matrices[1] = tms::GenMat{tms::GenMatRow{1,1,1}, tms::GenMatRow{0,1,0}, tms::GenMatRow{0,0,1}} );
	REQUIRE_NOTHROW( generating_matrices[2] = tms::GenMat{tms::GenMatRow{0,1,1}, tms::GenMatRow{1,1,0}, tms::GenMatRow{0,0,1}} );
	tms::Niederreiter nondeg_net(3, 3);

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
		CHECK( generating_numbers[0] == 4 );
		CHECK( generating_numbers[1] == 2 );
		CHECK( generating_numbers[2] == 1 );
		generating_numbers = nondeg_net.generating_numbers(1);
		CHECK( generating_numbers[0] == 4 );
		CHECK( generating_numbers[1] == 6 );
		CHECK( generating_numbers[2] == 5 );
		generating_numbers = nondeg_net.generating_numbers(2);
		CHECK( generating_numbers[0] == 2 );
		CHECK( generating_numbers[1] == 6 );
		CHECK( generating_numbers[2] == 5 );
	}

	SECTION("Construction of 3-dimensional point (classical enumeration of points)")
	{
		// classical representation of index 4 is (0 0 1)^T
		tms::Point point = nondeg_net.generate_point_classical(4);
		CHECK( point[0] == Catch::Approx(0.125).margin(0.0001) );
		CHECK( point[1] == Catch::Approx(0.625).margin(0.0001) );
		CHECK( point[2] == Catch::Approx(0.625).margin(0.0001) );
	}

	SECTION("Construction of 3-dimensional point (Gray code enumeration of points)")
	{
		// Gray representation of index 4 is (0 1 1)^T
		tms::Point point = nondeg_net.generate_point(4);
		CHECK( point[0] == Catch::Approx(0.375).margin(0.0001) );
		CHECK( point[1] == Catch::Approx(0.375).margin(0.0001) );
		CHECK( point[2] == Catch::Approx(0.375).margin(0.0001) );
	}

	SECTION("Construction of 3-dimensional point (Gray code enumeration of points, integer point)")
	{
		// Gray representation of index 4 is (0 1 1)^T
		tms::IntPoint point = nondeg_net.generate_int_point(4);
		CHECK( point[0] == 3 );
		CHECK( point[1] == 3 );
		CHECK( point[2] == 3 );
	}

	SECTION("Estimation of sequence defect")
	{
		REQUIRE( nondeg_net.t_estimate() == 1 );
	}
}



TEST_CASE("Validation of Niederreiter class, polynomial constructor", "[nets][Niederreiter]")
{
	tms::Niederreiter nondeg_net(3, 3);
	REQUIRE_THROWS( nondeg_net = tms::Niederreiter(3, {{1,1,1}, {1,0,1,1}, {1,0,0,1}}) );
	REQUIRE_NOTHROW( nondeg_net = tms::Niederreiter(3, {{1,1,1}, {1,0,1,1}, {1,1,0,1}}) );

	SECTION("Check parameter s")
	{
		REQUIRE( nondeg_net.s() == 3 );
	}

	SECTION("Construction of 3-dimensional point (classical enumeration of points)")
	{
		// classical representation of index 4 is (0 0 1)^T
		tms::Point point = nondeg_net.generate_point_classical(4);
		CHECK( point[0] == Catch::Approx(0.625).margin(0.0001) );
		CHECK( point[1] == Catch::Approx(0.875).margin(0.0001) );
		CHECK( point[2] == Catch::Approx(0.625).margin(0.0001) );
	}

	SECTION("Construction of 3-dimensional point (Gray code enumeration of points)")
	{
		// Gray representation of index 4 is (0 1 1)^T
		tms::Point point = nondeg_net.generate_point(4);
		CHECK( point[0] == Catch::Approx(0.375).margin(0.0001) );
		CHECK( point[1] == Catch::Approx(0.5).margin(0.0001) );
		CHECK( point[2] == Catch::Approx(0.875).margin(0.0001) );
	}

	SECTION("Construction of 3-dimensional point (Gray code enumeration of points, integer point)")
	{
		// Gray representation of index 4 is (0 1 1)^T
		tms::IntPoint point = nondeg_net.generate_int_point(4);
		CHECK( point[0] == 3 );
		CHECK( point[1] == 4 );
		CHECK( point[2] == 7 );
	}

	SECTION("Estimation of sequence defect")
	{
		REQUIRE( nondeg_net.t_estimate() == 5 );
	}
}
