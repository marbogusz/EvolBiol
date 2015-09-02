/*
 * StateTransitionEstimator.cpp
 *
 *  Created on: Oct 6, 2014
 *      Author: root
 */

#include <heuristics/PairStateTransitionEstimator.hpp>
#include "models/NegativeBinomialGapModel.hpp"

namespace EBC
{

PairStateTransitionEstimator::PairStateTransitionEstimator(IndelModel* im, unsigned char gc) :
		indelModel(im), gapCharacter(gc)
{
	DEBUG("Starting State Transition Estimator");
	//indelModel = new NegativeBinomialGapModel();
	maths = new Maths();

	modelParams = new OptimizedModelParameters(NULL, indelModel,0, 0, false,
	true, false, false, maths);
	modelParams->useIndelModelInitialParameters();

	bfgs = new Optimizer(modelParams, this,Definitions::OptimizationType::BFGS);

	maxTime = Definitions::almostZero;
}

double PairStateTransitionEstimator::runIteration()
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

void PairStateTransitionEstimator::addPair(vector<unsigned char>* s1,
		vector<unsigned char>* s2, double time)
{
	if (time > maxTime)
			maxTime = time;
	StateTransitionML* stml = new StateTransitionML(indelModel, time, gapCharacter, false);
	stml->addSample(s1,s2);
	stmSamples.push_back(stml);
}

void PairStateTransitionEstimator::addPair(vector<SequenceElement*>* s1,
		vector<SequenceElement*>* s2, double time)
{
	if (time > maxTime)
			maxTime = time;
	StateTransitionML* stml = new StateTransitionML(indelModel, time, gapCharacter, false);
	stml->addSample(s1,s2);
	stmSamples.push_back(stml);
}

void PairStateTransitionEstimator::optimize()
{
	modelParams->boundLambdaBasedOnDivergence(maxTime);
	bfgs->optimize();
	indelModel->setParameters(modelParams->getIndelParameters());
	INFO("StateTransitionEstimator results:");
	modelParams->logParameters();
	//modelParams->outputParameters();
}

void PairStateTransitionEstimator::clean()
{
    for(auto tm : stmSamples)
    {
    	delete tm;
    }
}

PairStateTransitionEstimator::~PairStateTransitionEstimator()
{
	delete bfgs;
	delete modelParams;
    delete maths;
    //for(auto tm : stmSamples)
    //{
    //	delete tm;
    //}
    //delete indelModel;
}

} /* namespace EBC */
