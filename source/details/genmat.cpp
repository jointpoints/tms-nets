#include "../../include/tms-nets/details/genmat.hpp"





tms::GenMat
tms::genmat::eye(BasicInt const size)
{
	GenMat gen_mat(size);
	
	for (BasicInt i = 0; i < size; ++i)
	{
		gen_mat[i][i] = 1;
	}
	
	return gen_mat;
}

tms::GenMat
tms::genmat::make_phi(GenMat const &a, GenMat const &b)
{
	if ( a.size() == b.size() )
	{
		return b.inverse() * a;
	}
	else
	{
		throw std::length_error("\n");
	}
}

tms::GenMat
tms::genmat::make_psi(GenMat const &a, GenMat const &b)
{
	if ( a.size() == b.size() )
	{
		return a * b.inverse();
	}
	else
	{
		throw std::length_error("\n");
	}
}
