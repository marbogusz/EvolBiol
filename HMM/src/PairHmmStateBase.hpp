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
	//std::map<PairHmmStateBase*, double> transitionProbs;

	double transFromMatch;
	double transFromInsert;
	double transFromDelete;

public:

	//void addTransitionProbabilityFrom(PairHmmStateBase* state, double value);

	inline void setTransitionProbabilityFromMatch(double value)
	{
		transFromMatch = value;
	}

	inline void setTransitionProbabilityFromInsert(double value)
	{
		transFromInsert = value;
	}

	inline void setTransitionProbabilityFromDelete(double value)
	{
		transFromDelete = value;
	}

	inline double getTransitionProbabilityFromMatch()
	{
		return transFromMatch;
	}

	inline double getTransitionProbabilityFromInsert()
	{
		return transFromInsert;
	}
	inline double getTransitionProbabilityFromDelete()
	{
		return transFromDelete;
	}


};

} /* namespace EBC */
#endif /* PAIRHMMSTATEBASE_H_ */
