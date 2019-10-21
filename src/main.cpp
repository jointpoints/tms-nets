#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "niederreiter2.h"

int main();
void generateNiederreiter2(uint64_t, uint64_t);



int main(void)
{
	uint64_t	dim		(4);
	uint64_t	amount	(1 << 4);
	
	timestamp();
	std::cout << "\n\n\n";

	generateNiederreiter2(dim, amount);
	
	std::cout << "\n\n\nNormal end of execution.\n";
	timestamp();

	return 0;
}



void generateNiederreiter2(uint64_t dim, uint64_t amount)
{
	int dim_num;
	int i;
	int j;
	double r[dim];
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
