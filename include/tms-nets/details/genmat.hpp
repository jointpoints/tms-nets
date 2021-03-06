#ifndef TMS_NETS_GENMAT_HPP
#define TMS_NETS_GENMAT_HPP

#include "common.hpp"


namespace tms::genmat
{	

	GenMat eye(BasicInt const size);
	
	GenMat make_phi(GenMat const &a, GenMat const &b);
	
	GenMat make_psi(GenMat const &a, GenMat const &b);

	

	
	
	
	GenMat
	eye(BasicInt const size)
	{
		GenMat gen_mat(size);
		
		for (BasicInt i = 0; i < size; ++i)
		{
			gen_mat[i][i] = 1;
		}
		
		return gen_mat;
	}
	
	GenMat
	make_phi(GenMat const &a, GenMat const &b)
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
	
	GenMat
	make_psi(GenMat const &a, GenMat const &b)
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

}


#endif
