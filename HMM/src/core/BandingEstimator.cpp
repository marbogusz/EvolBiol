/*
 * BandingEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "core/BandingEstimator.hpp"
#include "core/PhylogeneticTree.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"
#include "hmm/DpMatrixFull.hpp"
#include "extras/LikelihoodSurfacePlotter.hpp"
#include <chrono>
#include <ctime>

namespace EBC
{

BandingEstimator::BandingEstimator(Definitions::AlgorithmType at, Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params, Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, GuideTree* g) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories), pairCount(inputSequences->getPairCount()),
				/*hmms(pairCount), bands(pairCount),*/ divergenceTimes(pairCount), algorithm(at), gt(g)
{
	//Banding estimator means banding enabled!

	DEBUG("Starting Banding Estimator");
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

	modelParams->boundDivergenceBasedOnLambda(indel_params[0]);

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

/*
	for(unsigned int i =0; i< pairCount; i++)
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
*/
	numopt = new BrentOptimizer(modelParams, NULL);
	//bfgs->optimize();
	//this->modelParams->logParameters();

}

BandingEstimator::~BandingEstimator()
{
	//for(auto hmm : hmms)
	//	delete hmm;
	delete numopt;
	delete modelParams;
    delete maths;
    //for (auto bnd : bands)
    //	delete bnd;
    delete indelModel;
    delete substModel;
}

void BandingEstimator::optimizePairByPair()
{
	EvolutionaryPairHMM* hmm;
	Band* band;
	DistanceMatrix* dm = gt->getDistanceMatrix();
	PairHmmCalculationWrapper* wrapper = new PairHmmCalculationWrapper();
	double result;
	chrono::time_point<chrono::system_clock> start, end;
	chrono::duration<double> elapsed_seconds;

	ifstream myfile ("tree.sim_1");
	string newick;
	if (myfile.is_open())
	{
		getline (myfile,newick);
		myfile.close();
	}

	PhylogeneticTree ptree(this->inputSequences);
	ptree.fromNewick(newick);




	for(unsigned int i =0; i< pairCount; i++)
	{
		DEBUG("Optimizing distance for pair #" << i);
		std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);
		INFO("Running pairwise calculator for sequence id " << idxs.first << " and " << idxs.second
				<< " ,number " << i+1 <<" out of " << pairCount << " pairs" );
		BandCalculator* bc = new BandCalculator(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
				substModel, indelModel, gt->getDistanceMatrix()->getDistance(idxs.first,idxs.second));
		band = bc->getBand();
		if (algorithm == Definitions::AlgorithmType::Viterbi)
		{
			DEBUG("Creating Viterbi algorithm to optimize the pairwise divergence time...");
			hmm = new ViterbiPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					substModel, indelModel, Definitions::DpMatrixType::Full, band);
		}
		else if (algorithm == Definitions::AlgorithmType::Forward)
		{
			DEBUG("Creating forward algorithm to optimize the pairwise divergence time...");
			hmm = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					substModel, indelModel, Definitions::DpMatrixType::Full, band);
		}

		//hmm->setDivergenceTimeAndCalculateModels(modelParams->getDivergenceTime(0)); //zero as there's only one pair!

		//LikelihoodSurfacePlotter lsp;
		//lsp.setTargetHMM(hmm);
		//lsp.getLikelihoodSurface();


		wrapper->setTargetHMM(hmm);
		DUMP("Set model parameter in the hmm...");
		wrapper->setModelParameters(modelParams);
		modelParams->setUserDivergenceParams({bc->getClosestDistance()});
		//modelParams->setUserDivergenceParams({ptree.distanceById(idxs.first,idxs.second)});
		numopt->setTarget(wrapper);
		numopt->setAccuracy(bc->getBrentAccuracy());

		start = chrono::system_clock::now();

		result = numopt->optimize() * -1.0;

		end = chrono::system_clock::now();
		elapsed_seconds = end-start;

		INFO("Estimated Divergence Time " << modelParams->getDivergenceTime(0));
		INFO("Real divergence time " << ptree.distanceById(idxs.first,idxs.second));
	    INFO("Computation time " << elapsed_seconds.count() << "s\n");


		DEBUG("Likelihood after pairwise optimization: " << result);
		if (result <= (Definitions::minMatrixLikelihood /2.0))
		{
			DEBUG("Optimization failed for pair #" << i << " Zero probability FWD");
			band->output();
			dynamic_cast<DpMatrixFull*>(hmm->M->getDpMatrix())->outputValuesWithBands(band->getMatchBand() ,band->getInsertBand(),band->getDeleteBand(),'|', '-');
			dynamic_cast<DpMatrixFull*>(hmm->X->getDpMatrix())->outputValuesWithBands(band->getInsertBand(),band->getMatchBand() ,band->getDeleteBand(),'\\', '-');
			dynamic_cast<DpMatrixFull*>(hmm->Y->getDpMatrix())->outputValuesWithBands(band->getDeleteBand(),band->getMatchBand() ,band->getInsertBand(),'\\', '|');
		}
		this->divergenceTimes[i] = modelParams->getDivergenceTime(0);

		delete band;
		delete bc;
		delete hmm;
	}

	DEBUG("Optimized divergence times:");
	DEBUG(this->divergenceTimes);
}


double BandingEstimator::runIteration()
{
	double result = 0;
	double tmp;
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

