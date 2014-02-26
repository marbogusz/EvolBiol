/*
 * PairHMMstate.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#include "PairHmmState.hpp"

namespace EBC
{

PairHmmState::PairHmmState(unsigned int xSize, unsigned int ySize):
		DpMatrix(xSize,ySize)
{
	// TODO Auto-generated constructor stub

}

PairHmmState::~PairHmmState()
{
	// TODO Auto-generated destructor stub
}

void PairHmmState::addTransitionProbability(PairHmmState* state, double value)
{
	this->transitionProbs.insert(std::make_pair(state,value));
}

void PairHmmState::setTransitionProbability(PairHmmState* state, double value)
{
	this->transitionProbs[state]  = value;
}

double PairHmmState::getTransitionProbability(PairHmmState* state)
{
	return this->transitionProbs[state];
}

double& PairHmmState::operator ()(unsigned int row, unsigned int column)
{
	return this->matrixData[row][column].score;
}

} /* namespace EBC */


