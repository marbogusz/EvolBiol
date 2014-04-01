/*
 * ReducedPairHmmDeleteState.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: root
 */

#include "ReducedPairHmmDeleteState.hpp"

namespace EBC
{

ReducedPairHmmDeleteState::ReducedPairHmmDeleteState(unsigned int xSize, unsigned int ySize)
	: DpReducedMatrix(xSize,ySize)
{

}

ReducedPairHmmDeleteState::~ReducedPairHmmDeleteState()
{

}

void ReducedPairHmmDeleteState::initializeData()
{
	previousRow = buffer[0];
	currentRow = buffer[1];
	currentRow[0] = previousRow[0] = minVal;
}

} /* namespace EBC */
