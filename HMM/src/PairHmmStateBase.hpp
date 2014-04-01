/*
 * PairHMMstate.h
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#ifndef PAIRHMMSTATEBASE_H_
#define PAIRHMMSTATEBASE_H_

#include "Dictionary.hpp"

#include <map>
#include <algorithm>
#include <cmath>


namespace EBC
{

class PairHmmStateBase
{
protected:
	//map of transitions TO the current state

	//TODO - maybe a pointer TO a model class member...
	std::map<PairHmmStateBase*, double> transitionProbs;

public:

	void addTransitionProbabilityFrom(PairHmmStateBase* state, double value);

	void setTransitionProbability(PairHmmStateBase* state, double value);

	double getTransitionProbabilityFrom(PairHmmStateBase* state);

};

} /* namespace EBC */
#endif /* PAIRHMMSTATEBASE_H_ */
