/*
 * ReducedPairHmmInsertState.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: root
 */

#include "ReducedPairHmmInsertState.hpp"

namespace EBC
{

ReducedPairHmmInsertState::ReducedPairHmmInsertState(unsigned int xSize, unsigned int ySize)
	: DpReducedMatrix(xSize,ySize)
{

}

ReducedPairHmmInsertState::~ReducedPairHmmInsertState()
{

}

void ReducedPairHmmInsertState::initializeData()
{
	previousRow = buffer[0];
	currentRow = buffer[1];

	for (unsigned int i=0; i< ySize; i++)
	{
		currentRow[i] = previousRow[i] = minVal;
	}
}

} /* namespace EBC */
