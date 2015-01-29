/*
 * HMMEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "sampling/HMMEstimator.hpp"

#include <chrono>
#include <array>
#include <map>

using namespace std;

namespace EBC
{

HMMEstimator::HMMEstimator(Sequences* inputSeqs, Definitions::ModelType model ,
		Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, bool estimateAlpha) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories),
				gtree(new GuideTree(inputSeqs)), tst(*gtree)
{
	DEBUG("HMM Estimator starting");
	DEBUG("About to sample some triplets");
	DEBUG("Sampling triplets of sequences for gamma shape parameter estimation");
	
	maths = new Maths();
	dict = inputSequences->getDictionary();
	tal = new TripletAligner (inputSequences, gtree->getDistanceMatrix());

	tripletIdxs = tst.sampleFromTree();

	tripletIdxsSize = tripletIdxs.size();

	//2 pairs for every triplet!
	sampleWorkers.reserve(tripletIdxsSize*2);

	//FIXME - release the memory!!!! - delete pair objects and vectors (arrays) of SeqEls


    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();



	if (model == Definitions::ModelType::GTR)
	{
		substModel = new GTRModel(dict, maths,gammaRateCategories);
		DUMP("SME: Creating new GTR model");
	}
	else if (model == Definitions::ModelType::HKY85)
	{
		substModel = new HKY85Model(dict, maths,gammaRateCategories);
		DUMP("SME: Creating new HKY model");
	}
	else if (model == Definitions::ModelType::LG)
	{
			substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
			DUMP("SME: Creating new LG model");
	}

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());

	indelModel = new NegativeBinomialGapModel();

	modelParams = new OptimizedModelParameters(substModel, indelModel, 0, 0, true, true, estimateAlpha, false, maths);

	bfgs = new Optimizer(modelParams, this, Definitions::OptimizationType::BFGS);

	this->calculateInitialPairs(model);
	this->optimise();
	//this->runIteration();

	end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;

    INFO("Model Estimator elapsed time: " << elapsed_seconds.count() << " seconds");

	//substModel->summarize();

}

void HMMEstimator::calculateInitialPairs(Definitions::ModelType model)
{
	DEBUG("HMM Estimator - create initial pairs");
	//amino acid mode
	bool aaMode = false;
	double tmpd;

	double initAlpha = 0.75;
	double initKappa = 2.0;
	double initLambda = 0.02;
	double initEpsilon = 0.6;
	//k-mers tend to underestimate the distances;
	double initTimeModifier = 1.5;

	//if rateCat =  1 alpha does not matter.
	substModel->setAlpha(initAlpha);
	modelParams->setAlpha(initAlpha);

	if (model == Definitions::ModelType::GTR)
	{
		substModel->setParameters({1.0,1.0/initKappa,1.0/initKappa,1.0/initKappa,1.0/initKappa});
		modelParams->setUserSubstParams({1.0,1.0/initKappa,1.0/initKappa,1.0/initKappa,1.0/initKappa});
	}
	else if (model == Definitions::ModelType::HKY85)
	{
		substModel->setParameters({initKappa});
		modelParams->setUserSubstParams({initKappa});
	}
	//no else - no need to set anything for AAs

	modelParams->setUserIndelParams({initLambda, initEpsilon});

	substModel->calculateModel();

	indelModel->setParameters({initLambda,initEpsilon});

	for (int i = 0; i < tripletIdxs.size(); i++)
	{
		sampleWorkers.emplace_back(inputSequences->getSequencesAt(tripletIdxs[i][0]), inputSequences->getSequencesAt(tripletIdxs[i][1]),
				substModel, indelModel, gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][0],tripletIdxs[i][1]) * initTimeModifier);
		sampleWorkers.emplace_back(inputSequences->getSequencesAt(tripletIdxs[i][1]), inputSequences->getSequencesAt(tripletIdxs[i][2]),
				substModel, indelModel, gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][1],tripletIdxs[i][2]) * initTimeModifier);
	}
}

double EBC::HMMEstimator::runIteration() {
	double lnl = 0;

	substModel->setAlpha(modelParams->getAlpha());
	substModel->setParameters(modelParams->getSubstParameters());
	substModel->calculateModel();
	indelModel->setParameters(modelParams->getIndelParameters());

	for (auto &worker : sampleWorkers)
	{
		lnl += worker.optimiseDivergenceTime();
	}
	return lnl;

}

void EBC::HMMEstimator::optimise() {
	bfgs->optimize();
}


HMMEstimator::~HMMEstimator()
{
    delete maths;
    delete gtree;
    delete tal;
}

vector<double> HMMEstimator::getSubstitutionParameters()
{
	return this->modelParams->getSubstParameters();
}

vector<double> HMMEstimator::getIndelParameters()
{
	return this->modelParams->getIndelParameters();
}


double HMMEstimator::getAlpha()
{
	return this->modelParams->getAlpha();
}

} /* namespace EBC */


