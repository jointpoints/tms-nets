#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "niederreiter2.h"

int main();
void test01();



int main (void)
{
	timestamp();
	cout << "\n";
	cout << "NIEDERREITER2_PRB\n";
	cout << "	C++ version\n";
	cout << "Test the NIEDERREITER2 library.\n";

	test01();
	
	cout << "\n";
	cout << "NIEDERREITER2_PRB\n";
	cout << "	Normal end of execution.\n";

	cout << "\n";
	timestamp();

	return 0;
}



void test01 (void)
{
#	define DIM 2
	int dim_num;
	int i;
	int j;
	double r[DIM];
	int seed;
	
	seed = 0;
	
	for(i = 0; i < 16; ++i)
	{
		niederreiter2(DIM, &seed, r);
		
		for (j = 0; j < DIM; j++)
		{
			cout << setw(10) << r[j] << " ";
		}
		cout << "\n";
	}
	
	return;
#	undef DIM_MAX
}
