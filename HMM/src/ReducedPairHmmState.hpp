/*
 * ReducedPairHmmState.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: root
 */

#ifndef REDUCEDPAIRHMMSTATE_HPP_
#define REDUCEDPAIRHMMSTATE_HPP_

#include "DpReducedMatrix.hpp"
#include "PairHmmStateBase.hpp"

namespace EBC
{

class ReducedPairHmmState: public DpReducedMatrix, public PairHmmStateBase
{
public:
	ReducedPairHmmState(unsigned int xSize, unsigned int ySize);
	virtual ~ReducedPairHmmState();
};

} /* namespace EBC */

#endif /* REDUCEDPAIRHMMSTATE_HPP_ */
