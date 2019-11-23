#include <cstdlib>
#include <iostream>
#include <fstream>

#include "niederreiter2.hpp"



int main(void)
{
	int dim     = 2;
	int expo    = 15;
	
	using GeneratorType = sequences::Niederreiter<uint64_t, 63>;

	//Create generator of (t,s) sequences with s == dim and consequently implicitly defined t.
	GeneratorType generator(dim);

	//Generate 20 (not 1 << 20!!) points of (t,s)-sequence, beginning with 10 sequence number
	std::vector<sequences::Point> points = generator.get_points_real(0, 20);

	//Load another 20 (not 1 << 20!!!) points of (t,s)-sequence, beginning with 100 sequence number into points-vector
	generator.load_points_real(points, 100);


	std::vector<GeneratorType::IntPoint> points_int;
	//ERROR: because points_int are empty
	//	generator.load_points_int(points_int, 10);

	points_int = std::vector<GeneratorType::IntPoint>(20);
	//ERROR: because points_int containing empty points
	//	generator.load_points_int(points_int, 10);

	points_int = std::vector<GeneratorType::IntPoint>(20, GeneratorType::IntPoint(dim));
	//OK: 20 consequent points of (t,s)-sequence will be loaded into points_int
	generator.load_points_int(points_int, 10);
	
	return 0;
}
