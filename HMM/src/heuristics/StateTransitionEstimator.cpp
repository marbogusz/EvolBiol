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

StateTransitionEstimator::StateTransitionEstimator(IndelModel* im, Definitions::OptimizationType ot, unsigned int pc, unsigned char gc, bool useEq) :
		indelModel(im), stmSamples(pc), gapCharacter(gc), useStateEq(useEq)
{
	DEBUG("Starting State Transition Estimator");
	//indelModel = new NegativeBinomialGapModel();
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

	DEBUG("State Transition Estimator add time for triplet " << triplet << "\t pair " << pr << "\ttime " << time);
	stmSamples[2*triplet+pr] = new StateTransitionML(indelModel, time, gapCharacter, useStateEq);
}

void StateTransitionEstimator::addPair(vector<unsigned char>* s1,
		vector<unsigned char>* s2, unsigned int triplet, unsigned int pr)
{
	DUMP("State Transition Estimator add pair for triplet " << triplet << " and pair no " << pr );
	stmSamples[2*triplet+pr]->addSample(s1,s2);
}

void StateTransitionEstimator::optimize()
{
	bfgs->optimize();
	indelModel->setParameters(modelParams->getIndelParameters());
	INFO("StateTransitionEstimator results:");
	modelParams->logParameters();
	//modelParams->outputParameters();
}

void StateTransitionEstimator::clean()
{
    for(auto tm : stmSamples)
    {
    	delete tm;
    }
}

StateTransitionEstimator::~StateTransitionEstimator()
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
