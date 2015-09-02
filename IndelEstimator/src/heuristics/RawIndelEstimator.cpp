/*
 * ModelEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include <heuristics/RawIndelEstimator.hpp>
#include "core/Dictionary.hpp"
#include <chrono>
#include <array>
#include <map>

using namespace std;

namespace EBC
{

RawIndelEstimator::RawIndelEstimator(Sequences* inputSeqs, Definitions::ModelType model , DistanceMatrix* dm) :
				inputSequences(inputSeqs), distMat(dm), pst(dm), gammaRateCategories(1) // no heterogeniety
{

	DEBUG("Sampling pairs of sequences for indel estimation");

	maths = new Maths();
	dict = inputSequences->getDictionary();

	pairIdxs = pst.sampleFromTree();

	pairIdxsSize =pairIdxs.size();

	DEBUG(pairIdxsSize << " pairs sampled");

	if (pairIdxsSize == 0)
		throw ProgramException("No triplets selected for model estimation");

	fwdHMMs.resize(pairIdxsSize);
	pairwiseTimes.resize(pairIdxsSize);
	bands.resize(pairIdxsSize);

	maxTime = Definitions::almostZero;

	calculateInitials(model);

	modelParams = new OptimizedModelParameters(NULL, indelModel,0, 0, false,
	true, false, false, maths);
	modelParams->useIndelModelInitialParameters();
	bfgs = new Optimizer(modelParams, this,Definitions::OptimizationType::BFGS);


    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    modelParams->boundLambdaBasedOnDivergence(maxTime);
    bfgs->optimize();
    indelModel->setParameters(modelParams->getIndelParameters());

	end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;

    INFO("Model Estimator elapsed time: " << elapsed_seconds.count() << " seconds");
}

double RawIndelEstimator::runIteration()
{
	double lnl = 0;

	for (unsigned int i = 0; i<pairIdxsSize; i++){
		//set indel params!
		indelModel->setParameters(modelParams->getIndelParameters());
		DEBUG("Indel params " << modelParams->getIndelParameters()[0] << " " << modelParams->getIndelParameters()[1]);
		fwdHMMs[i]->setDivergenceTimeAndCalculateModels(pairwiseTimes[i]);
		lnl += fwdHMMs[i]->runAlgorithm();
	}
	return lnl;
}


void RawIndelEstimator::calculateInitials(Definitions::ModelType model)
{

	bool aaMode = false;
	double tmpd;

	double initAlpha = 0.75;
	double initKappa = 2.0;
	double initLambda = 0.025;
	double initEpsilon = 0.5;

	if (model == Definitions::ModelType::HKY85){
			//FIXME - make model idiotproof by checking if parameters are set;
			DEBUG("Setting HKY85");
			this->substModel = new HKY85Model(dict, maths,gammaRateCategories);
			substModel->setParameters({initKappa});
	}
	else if (model == Definitions::ModelType::GTR){
			DEBUG("Setting GTR");
			this->substModel = new GTRModel(dict, maths,gammaRateCategories);
			substModel->setParameters({1.0,1.0/initKappa,1.0/initKappa,1.0/initKappa,1.0/initKappa});
	}
	else if (model == Definitions::ModelType::LG){
			aaMode = true;
			substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
	}


	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
	substModel->calculateModel();

	indelModel = new NegativeBinomialGapModel();
	indelModel->setParameters({initLambda,initEpsilon});

	for (int i = 0; i < pairIdxsSize; i++)
	{

		vector<SequenceElement*>* seq1 = inputSequences->getSequencesAt(pairIdxs[i][0]);
		vector<SequenceElement*>* seq2 = inputSequences->getSequencesAt(pairIdxs[i][1]);
		pairwiseTimes[i] = distMat->getDistance(pairIdxs[i][0],pairIdxs[i][1]);
		if (pairwiseTimes[i] > maxTime)
			maxTime = pairwiseTimes[i];
		bands[i] = new Band(seq1->size(),seq2->size(),0.1);

		DEBUG("Forward HMM for sequence id " << pairIdxs[i][0] << " and " << pairIdxs[i][1] << " and distance " << pairwiseTimes[i]);
		fwdHMMs[i] = new ForwardPairHMM(seq1,seq2, substModel, indelModel, Definitions::DpMatrixType::Full, bands[i],true);
	}
}




RawIndelEstimator::~RawIndelEstimator()
{

	for (int i =0; i < pairIdxsSize; i++)
	{
		delete bands[i];
		delete fwdHMMs[i];

	}

//FIXME - clean up!

    delete maths;
}

} /* namespace EBC */
