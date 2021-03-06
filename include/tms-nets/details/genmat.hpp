#ifndef TMS_NETS_GENMAT_HPP
#define TMS_NETS_GENMAT_HPP

#include "common.hpp"


namespace tms::genmat
{	

	GenMat eye(BasicInt const size);
	
	GenMat make_phi(GenMat const &a, GenMat const &b);
	
	GenMat make_psi(GenMat const &a, GenMat const &b);

}


#endif
