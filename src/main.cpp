#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "niederreiter2.h"

void generateNiederreiter2(int, uint64_t, char const *);

#define DEFAULT_PATH_VALUE 		"generated_sequence.txt"
#define DEFAULT_DIM_VALUE		3ULL
#define DEFAULT_EXPO_VALUE		10ULL


int main(void)
{
	int			dim		 = DEFAULT_DIM_VALUE;
	int			expo	 = DEFAULT_EXPO_VALUE;
	char const *filePath = DEFAULT_PATH_VALUE;
	
	timestamp();
	std::cout << "\n\n\n";
	
	generateNiederreiter2(dim, 1ULL << expo, filePath);
	
	std::cout << "\n\n\nNormal end of execution.\n";
	timestamp();
	
	return 0;
}



void generateNiederreiter2(int dim, uint64_t amount, char const *filePath)
{
	uint64_t i;
	uint64_t j;
	Real r[dim];
	uint64_t seed;
	std::ofstream outFile(filePath);
	if ( !outFile.good() )
	{
		std::cout << "\nERROR: Bad file path.\n";
	}
	
	seed = 0;
	
	for (i = 0; i < amount; ++i)
	{
		generate_next_nied2_real(dim, &seed, r);
		for (j = 0; j < dim; j++)
		{
			outFile << r[j] << " ";
		}
		outFile << '\n';
	}
	outFile.close();
	
	return;
}
