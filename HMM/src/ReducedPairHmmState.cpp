/*
 * ReducedPairHmmState.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: root
 */

#include "ReducedPairHmmState.hpp"

namespace EBC
{

ReducedPairHmmState::ReducedPairHmmState(unsigned int xSize, unsigned int ySize)
	: DpReducedMatrix(xSize,ySize)
{

}

ReducedPairHmmState::~ReducedPairHmmState()
{

}

} /* namespace EBC */
