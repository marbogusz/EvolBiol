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

StateTransitionEstimator::StateTransitionEstimator(Definitions::OptimizationType ot, unsigned int pc) : stmSamples(pc)
{
	DEBUG("Starting State Transition Estimator");
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

void StateTransitionEstimator::addTime(double time, unsigned int triplet, unsigned int pr)
{
	stmSamples[2*triplet+pr] = new StateTransitionML(indelModel, time);
}

void StateTransitionEstimator::addPair(vector<SequenceElement>& s1,
		vector<SequenceElement>& s2, unsigned int triplet, unsigned int pr, double weight)
{
	DUMP("add pair for triplet " << triplet << " and pair no " << pr << " with weight " << weight);
	stmSamples[2*triplet+pr]->addSample(s1,s2, weight);
}

void StateTransitionEstimator::optimize()
{
	bfgs->optimize();
	INFO("StateTransitionEstimator results:");
	modelParams->logParameters();
	//modelParams->outputParameters();
}

StateTransitionEstimator::~StateTransitionEstimator()
{
	delete bfgs;
	delete modelParams;
    delete maths;
    for(auto tm : stmSamples)
    {
    	delete tm;
    }
    delete indelModel;
}

} /* namespace EBC */
