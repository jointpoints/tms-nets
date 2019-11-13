#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "niederreiter2.hpp"

#define DEFAULT_DIM_VALUE      	10ULL
#define DEFAULT_EXPO_VALUE      3ULL



int main(void)
{
	int dim     = DEFAULT_DIM_VALUE;
	int expo    = DEFAULT_EXPO_VALUE;
	
	//Create generator of (t,s) sequences with s == dim and consequently implicitly defined t.
	Nied2Generator generator(dim);
	
	//Generate 200 (not 1 << 200!!) points of (t,s)-sequence, beginning with 10 sequence number
	std::vector<Nied2Generator::RealPoint> points = generator.get_points_real(10, 200);

	//Load another 200 (not 1 << 200!!!) points of (t,s)-sequence, beginning with 100 sequence number into points-vector
	generator.load_points_real(points, 100);


	std::vector<Nied2Generator::IntPoint> points_int;
	//ERROR: because points_int are empty
	//generator.load_points_int(points_int, 10);

	points_int = std::vector<Nied2Generator::IntPoint>(20);
	//ERROR: because points_int containing empty points
	//generator.load_points_int(points_int, 10);

	points_int = std::vector<Nied2Generator::IntPoint>(20, Nied2Generator::IntPoint(dim));
	//OK: 20 consequent points of (t,s)-sequence will be loaded into points_int
	generator.load_points_int(points_int, 10);
	
	//kak-to tak
	
	return 0;
}
