/*!
 *	\file main.cpp
 *
 *	\author
 *		Michail Lazukov (Russian Technological University, KMBO-03-16, Russia, 2019),
 *		Avetik Vardanyan (Russian Technological University, KMBO-03-16, Russia, 2019)
 */



/*
 * Use the following line to control verbosity
 * Possible levels of verbosity: 1, 2, 3, 4
 */
#define VERBOSITY_LEVEL 2



#include "tests.hpp"



int main(void)
{
	using MyGeneratorType = sequences::Niederreiter<uint32_t, 32>;
	
	int cur   = 0;
	int s     = 0;
	int m     = 0;
	int s_max = 7;
	int m_max = 16;
	srand(time(0));
	while (cur < 100)
	{
	    cur += 1;
	    std::cout << "TEST #" << cur << "/100\n";
		s = rand() % s_max + 1;
	    MyGeneratorType my_generator(s);
        m = rand() % (m_max - my_generator.get_t() + 1) + my_generator.get_t();
		niederreiter_check_definition(&my_generator, m);
		niederreiter_check_uniqueness(&my_generator, m);
		std::cout << "----------------------------------\n\n";
	}

	return 0;
}
