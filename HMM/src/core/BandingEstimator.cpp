/*
 * BandingEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "core/BandingEstimator.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"
#include "hmm/DpMatrixFull.hpp"

namespace EBC
{

BandingEstimator::BandingEstimator(Definitions::AlgorithmType at, Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params, Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, GuideTree* gt) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories), pairCount(inputSequences->getPairCount()),
				hmms(pairCount), bands(pairCount), divergenceTimes(pairCount)
{
	//Banding estimator means banding enabled!

	FileLogger::DebugLogger() << "Starting Banding Estimator" << "\n";
	maths = new Maths();
	dict = inputSequences->getDictionary();

	//Helper models
	if (model == Definitions::ModelType::GTR)
	{
		substModel = new GTRModel(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::HKY85)
	{
		substModel = new HKY85Model(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::LG)
	{
			substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
	}

	indelModel = new NegativeBinomialGapModel();

	estimateSubstitutionParams = false;
	estimateIndelParams = false;
	estimateAlpha = false;
	/*
	estimateSubstitutionParams = subst_params.size() != substModel->getParamsNumber();
	estimateIndelParams = indel_params.size() == 0;
	//Do not estimate alpha here. Alpha needs to be estimated beforehand
	this->estimateAlpha = false;

	DEBUG("Pairwise banding model estimator starting");
	DEBUG("Estimate substitution parameters set to : " << estimateSubstitutionParams << " Estimate indel parameters set to : " << estimateIndelParams);
	DEBUG("Estimate alpha set to : " << estimateAlpha << " , rate categories " << gammaRateCategories << " , alpha : " << alpha);

	FileLogger::DebugLogger() << "Estimate substitution parameters set to : " << estimateSubstitutionParams << " Estimate indel parameters set to : " << estimateIndelParams << "\n";
	FileLogger::DebugLogger() << "Estimate alpha set to : " << estimateAlpha << " , rate categories " << gammaRateCategories << " , alpha : " << alpha << "\n";
	 */
	//pairwise comparison mode
	modelParams = new OptimizedModelParameters(substModel, indelModel,2, 1, estimateSubstitutionParams,
			estimateIndelParams, estimateAlpha, true, maths);

	if(!estimateIndelParams)
		modelParams->setUserIndelParams(indel_params);
	if(!estimateSubstitutionParams)
		modelParams->setUserSubstParams(subst_params);
	modelParams->setAlpha(alpha);

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
	if (estimateSubstitutionParams == false)
	{
		//set parameters and calculate the model
		substModel->setAlpha(modelParams->getAlpha());
		substModel->setParameters(modelParams->getSubstParameters());
		substModel->calculateModel();
	}

	if (estimateIndelParams == false)
	{
			//set parameters and calculate the model
		indelModel->setParameters(modelParams->getIndelParameters());
	}

	//let's assume that we have all the parameters estimated
	//need to get times!

	EvolutionaryPairHMM *hmm;
//non banden probs
//	vector<EvolutionaryPairHMM*> hmmsNB(pairCount);


	for(unsigned int i =0; i<pairCount; i++)
		{
			std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);
			DEBUG("Running band calculator for sequence " << idxs.first << " and " << idxs.second);
			FileLogger::InfoLogger() << "Running band calculator for sequence " << idxs.first << " and " << idxs.second << "\n";
			BandCalculator* bc = new BandCalculator(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					substModel, indelModel, gt->getDistanceMatrix()->getDistance(idxs.first,idxs.second));
			bands[i] = bc->getBand();
			if (at == Definitions::AlgorithmType::Viterbi)
			{
				hmm = hmms[i] = new ViterbiPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					substModel, indelModel, Definitions::DpMatrixType::Full, bands[i]);
			}
			else if (at == Definitions::AlgorithmType::Forward)
			{
				hmm = hmms[i] = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
						substModel, indelModel, Definitions::DpMatrixType::Full, bands[i]);
				//hmmsNB[i] = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
				//						substModel, indelModel, Definitions::DpMatrixType::Full, NULL);
			}
			else
			{
				throw HmmException("Wrong algorithm type - use either Forward or viterbi\n");
			}
			delete bc;
		}

	bfgs = new Optimizer(modelParams, NULL, ot);
	//bfgs->optimize();
	//this->modelParams->logParameters();

}

BandingEstimator::~BandingEstimator()
{
	// TODO Auto-generated destructor stub
	for(auto hmm : hmms)
		delete hmm;
	delete bfgs;
	delete modelParams;
    delete maths;
    for (auto bnd : bands)
    	delete bnd;
    delete indelModel;
    delete substModel;
}

void BandingEstimator::optimizePairByPair()
{
	EvolutionaryPairHMM* hmm;
	PairHmmCalculationWrapper* wrapper = new PairHmmCalculationWrapper();
	double result;

	for(unsigned int i =0; i<pairCount; i++)
	{
		FileLogger::DebugLogger() << "Optimizing distance for pair #" << i << '\n';
		hmm = hmms[i];
		//hmm->setDivergenceTime(modelParams->getDivergenceTime(0)); //zero as there's only one pair!
		wrapper->setTargetHMM(hmm);
		wrapper->setModelParameters(modelParams);
		bfgs->setTarget(wrapper);
		result = bfgs->optimize() * -1.0;
		FileLogger::DebugLogger() << "Resulting likelihood " << result << "\n";
		if (result <= Definitions::minMatrixLikelihood)
		{
			FileLogger::ErrorLogger() << "Optimization failed for pair #" << i << " Zero probability FWD" << '\n';
			bands[i]->output();
			dynamic_cast<DpMatrixFull*>(hmm->M->getDpMatrix())->outputValues(0);
			dynamic_cast<DpMatrixFull*>(hmm->X->getDpMatrix())->outputValues(0);
			dynamic_cast<DpMatrixFull*>(hmm->Y->getDpMatrix())->outputValues(0);
		}
		this->divergenceTimes[i] = modelParams->getDivergenceTime(0);
	}

	FileLogger::DebugLogger() << this->divergenceTimes;

}


double BandingEstimator::runIteration()
{
	double result = 0;
	double tmp;
	EvolutionaryPairHMM* hmm;

	if (estimateSubstitutionParams == true)
	{
			//set parameters and calculate the model
		if(this->estimateAlpha)
		{
			substModel->setAlpha(modelParams->getAlpha());
		}
		substModel->setParameters(modelParams->getSubstParameters());
		substModel->calculateModel();
	}

	if (estimateIndelParams == true)
	{
		//set parameters and calculate the model
		indelModel->setParameters(modelParams->getIndelParameters());
	}

	DEBUGN("DBG pair likelihoods : ");
	for(unsigned int i =0; i<pairCount; i++)
	{
		hmm = hmms[i];
		//FIXME - individual indel models!!! or gap opening extension class aggregator!
		// the following is not thread safe for indels!!!
		hmm->setDivergenceTime(modelParams->getDivergenceTime(i));
		//indelModel->setTime(modelParams->getDivergenceTime(i));
		//indelModel->calculate();
		tmp = hmm->runAlgorithm();
		DEBUGN("\t\t" << tmp << " " << modelParams->getDivergenceTime(i));
		FileLogger::DebugLogger() << tmp << " " << modelParams->getDivergenceTime(i) << '\n';
		result += tmp;
		//modelParams->outputParameters();
	}
	DEBUGN("\n");
	//cerr << " lnl " << result << endl;
	return result;
}

void BandingEstimator::outputDistanceMatrix(stringstream& ss)
{
	unsigned int count, pairCount;
	count = this->inputSequences->getSequenceCount();
	pairCount = this->inputSequences->getPairCount();

	ss << "\t" << this->inputSequences->getSequenceCount() << endl;

	for (unsigned int i = 0; i< count; i++)
	{
		ss << "S" << i << " ";
		for(unsigned int j=0; j< count; j++)
		{
			ss << this->modelParams->getDistanceBetween(i,j) << " ";
		}
		ss << endl;
	}
}

} /* namespace EBC */

