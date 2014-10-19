/*
 * StateTransitionEstimator.cpp
 *
 *  Created on: Oct 6, 2014
 *      Author: root
 */

#include "heuristics/StateTransitionEstimator.hpp"
#include "models/NegativeBinomialGapModel.hpp"

namespace EBC
{

StateTransitionEstimator::StateTransitionEstimator(Definitions::OptimizationType ot)
{
	DEBUG("Starting State Transition Estimator");
	// TODO Auto-generated constructor stub
	indelModel = new NegativeBinomialGapModel();
	maths = new Maths();

	modelParams = new OptimizedModelParameters(NULL, indelModel,0, 0, false,
	true, false, false, maths);

	bfgs = new Optimizer(modelParams, this,ot);
}

double StateTransitionEstimator::runIteration()
{
	double result = 0;

	//modelParams->outputParameters();
	indelModel->setParameters(modelParams->getIndelParameters());
	//go through matrices
	for(auto tm : stmSamples)
	{
		//set parameter for every sample
		result += tm->getLnL();
		//calculate and add to result
	}
	return result * -1.0;
}

void StateTransitionEstimator::addPair(vector<SequenceElement>& s1,
		vector<SequenceElement>& s2, double time)
{
	StateTransitionML* nmat = new StateTransitionML(indelModel,s1,s2,time);
	this->stmSamples.push_back(nmat);

}

void StateTransitionEstimator::optimize()
{
	bfgs->optimize();
	//modelParams->outputParameters();
}

StateTransitionEstimator::~StateTransitionEstimator()
{
	delete bfgs;
	delete modelParams;
    delete maths;
}

} /* namespace EBC */
