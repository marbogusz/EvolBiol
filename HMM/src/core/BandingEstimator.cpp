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
				inputSequences(inputSeqs), gammaRateCategories(rateCategories), pairCount(inputSequences->getPairCount()), hmms(pairCount), bands(pairCount)
{
	//Banding estimator means banding enabled!

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

	estimateSubstitutionParams = subst_params.size() != substModel->getParamsNumber();
	estimateIndelParams = indel_params.size() == 0;
	//Do not estimate alpha here. Alpha needs to be estimated beforehand
	this->estimateAlpha = false;

	DEBUG("Pairwise banding model estimator starting");
	DEBUG("Estimate substitution parameters set to : " << estimateSubstitutionParams << " Estimate indel parameters set to : " << estimateIndelParams);
	DEBUG("Estimate alpha set to : " << estimateAlpha << " , rate categories " << gammaRateCategories << " , alpha : " << alpha);

	modelParams = new OptimizedModelParameters(substModel, indelModel,inputSequences->getSequenceCount(), pairCount, estimateSubstitutionParams,
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
/*
	hmms[0]->setDivergenceTime(0.5);
	hmmsNB[0]->setDivergenceTime(0.5);

	hmms[0]->runAlgorithm();
	hmmsNB[0]->runAlgorithm();

	DEBUG("Match Banded");
	dynamic_cast<DpMatrixFull*>(hmms[0]->M->getDpMatrix())->outputValues(0);
	DEBUG("Match UNbanded");
	dynamic_cast<DpMatrixFull*>(hmmsNB[0]->M->getDpMatrix())->outputValues(0);

	DEBUG("Insert Banded");
	dynamic_cast<DpMatrixFull*>(hmms[0]->X->getDpMatrix())->outputValues(0);
	DEBUG("Insert UNbanded");
	dynamic_cast<DpMatrixFull*>(hmmsNB[0]->X->getDpMatrix())->outputValues(0);
	DEBUG("Delete Banded");
	dynamic_cast<DpMatrixFull*>(hmms[0]->Y->getDpMatrix())->outputValues(0);
	DEBUG("Delete UNbanded");
	dynamic_cast<DpMatrixFull*>(hmmsNB[0]->Y->getDpMatrix())->outputValues(0);
*/

	bfgs = new Optimizer(modelParams, this, ot);
	bfgs->optimize();
	this->modelParams->outputParameters();

/*
	for(unsigned int i =0; i<pairCount; i++)
	{
		std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);

		if (at == Definitions::AlgorithmType::Viterbi)
		{
			hmm = hmms[i] = new ViterbiPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
				false, substModel, indelModel, 0, Definitions::DpMatrixType::Limited);
		}
		else if (at == Definitions::AlgorithmType::Forward)
		{
			hmm = hmms[i] = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					false, substModel, indelModel, 0, Definitions::DpMatrixType::Limited);
		}
		else
		{
			throw HmmException("Wrong algorithm type - use either Forward or viterbi\n");
		}
	}

	bfgs = new Optimizer(modelParams, this, ot);
	bfgs->optimize();
*/
}

BandingEstimator::~BandingEstimator()
{
	// TODO Auto-generated destructor stub
	delete bfgs;
	delete modelParams;
    delete maths;
    for (auto bnd : bands)
    	delete bnd;
}

double BandingEstimator::runIteration()
{
	double result = 0;
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

	for(unsigned int i =0; i<pairCount; i++)
	{
		hmm = hmms[i];
		//FIXME - individual indel models!!! or gap opening extension class aggregator!
		// the following is not thread safe for indels!!!
		hmm->setDivergenceTime(modelParams->getDivergenceTime(i));
		//indelModel->setTime(modelParams->getDivergenceTime(i));
		//indelModel->calculate();
		result += hmm->runAlgorithm();
		//modelParams->outputParameters();
	}
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
