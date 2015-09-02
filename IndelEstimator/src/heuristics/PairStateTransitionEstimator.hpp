/*
 * StateTransitionEstimator.hpp
 *
 *  Created on: Oct 6, 2014
 *      Author: root
 */

#ifndef PAIRSTATETRANSITIONESTIMATOR_HPP_
#define PAIRSTATETRANSITIONESTIMATOR_HPP_

#include "core/IOptimizable.hpp"
#include "core/Optimizer.hpp"

#include "heuristics/StateTransitionML.hpp"

namespace EBC
{

class PairStateTransitionEstimator : public IOptimizable
{
protected:

	Optimizer* bfgs;

	OptimizedModelParameters* modelParams;

	vector<StateTransitionML*> stmSamples;

	IndelModel* indelModel;

	Maths* maths;

	unsigned char gapCharacter;

	double maxTime;

public:
	PairStateTransitionEstimator(IndelModel* im, unsigned char gapchar);

	void addPair(vector<unsigned char>* s1, vector<unsigned char>* s2, double time);

	void addPair(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, double time);

	double runIteration();

	void optimize();

	void clean();

	virtual ~PairStateTransitionEstimator();

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
