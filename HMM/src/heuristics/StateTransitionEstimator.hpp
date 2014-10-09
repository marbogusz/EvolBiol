/*
 * StateTransitionEstimator.hpp
 *
 *  Created on: Oct 6, 2014
 *      Author: root
 */

#ifndef STATETRANSITIONESTIMATOR_HPP_
#define STATETRANSITIONESTIMATOR_HPP_

#include "core/IOptimizable.hpp"
#include "core/Optimizer.hpp"

#include "heuristics/StateTransitionML.hpp"

namespace EBC
{

class StateTransitionEstimator : public IOptimizable
{
protected:

	Optimizer* bfgs;

	OptimizedModelParameters* modelParams;

	vector<StateTransitionML*> stmSamples;

	IndelModel* indelModel;

	Maths* maths;


public:
	StateTransitionEstimator(Definitions::OptimizationType ot);

	void addPair(vector<SequenceElement>& s1,vector<SequenceElement>& s2, double time);

	double runIteration();

	void optimize();

	virtual ~StateTransitionEstimator();

	OptimizedModelParameters* getModelParams()
	{
		return modelParams;
	}

	IndelModel* getIndelModel()
	{
		return indelModel;
	}
};

} /* namespace EBC */

#endif /* STATETRANSITIONESTIMATOR_HPP_ */
