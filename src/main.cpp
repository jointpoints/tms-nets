#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "niederreiter2.h"
#include "tests.h"

int main();
void generateNiederreiter2(uint64_t, uint64_t);



int main(void)
{
	/*
	 * Change these two variables to variate amount of points and dimension of space
	 */
	uint64_t	dim		(5);
	uint64_t	amount	(1 << 16);
	
	timestamp();
	std::cout << "\n\n\n";
	
	/*
	 * Uncomment lines below to perform specified actions
	 */
	//generateNiederreiter2(dim, amount);
	//checkDefinition(dim, amount);
	
	std::cout << "\n\n\nNormal end of execution.\n";
	timestamp();

	return 0;
}



void generateNiederreiter2(uint64_t dim, uint64_t amount)
{
	int dim_num;
	int i;
	int j;
	int r[dim];
	int seed;
	std::ofstream outFile(".\\producedNets\\generated_sequence.txt");
	
	seed = 0;
	
	for(i = 0; i < amount; ++i)
	{
		niederreiter2(dim, &seed, r);
		
		for (j = 0; j < dim; j++)
		{
			outFile << r[j] << " ";
		}
		outFile << '\n';
	}
	outFile.close();
	
	return;
}
