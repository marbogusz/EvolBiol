/*
 * ReducedPairHmmMatchState.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: root
 */

#include "ReducedPairHmmMatchState.hpp"

namespace EBC
{

ReducedPairHmmMatchState::ReducedPairHmmMatchState(unsigned int xSize, unsigned int ySize)
	: DpReducedMatrix(xSize,ySize)
{

}

ReducedPairHmmMatchState::~ReducedPairHmmMatchState()
{

}

void ReducedPairHmmMatchState::initializeData()
{
	previousRow = buffer[0];
	currentRow = buffer[1];

	for (unsigned int i=0; i< xSize; i++)
	{
		currentRow[i] = previousRow[i] = minVal;
	}
	currentRow[0] = 0;
}

} /* namespace EBC */
