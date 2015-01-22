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

	unsigned char gapCharacter


public:
	StateTransitionEstimator(IndelModel* im, Definitions::OptimizationType ot, unsigned int pairCategories, unsigned char);

	void addTime(double time, unsigned int triplet, unsigned int pr);

	void addPair(vector<unsigned char>* s1, vector<unsigned char>* s2, unsigned int triplet, unsigned int pr, double weight);

	double runIteration();

	void optimize();

	void clean();

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
