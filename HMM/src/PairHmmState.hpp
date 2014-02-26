/*
 * PairHMMstate.h
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#ifndef PAIRHMMSTATE_H_
#define PAIRHMMSTATE_H_

#include "DpMatrix.hpp"
#include "Dictionary.hpp"

#include <map>
#include <algorithm>


namespace EBC
{

class PairHmmState: public EBC::DpMatrix
{
protected:

	//index differences depending on the state
	//(diagonal (-1,-1) insertion (-1,0), deletion(0,-1);
	int xDifference;
	int yDifference;

	//map of transitions TO the current state

	//TODO - maybe a pointer TO a model class member...
	std::map<PairHmmState*, double> transitionProbs;


public:
	double& operator ()(unsigned int row, unsigned int column);

	PairHmmState(unsigned int xSize, unsigned int ySize);

	void addTransitionProbability(PairHmmState* state, double value);

	void setTransitionProbability(PairHmmState* state, double value);

	double getTransitionProbability(PairHmmState* state);

	virtual ~PairHmmState();
};

} /* namespace EBC */
#endif /* PAIRHMMSTATE_H_ */
