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

	EvolutionaryPairHMM *hmm;

	numopt = new BrentOptimizer(modelParams, NULL);

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
	//	BandCalculator* bc = new BandCalculator(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
	//					substModel, indelModel, ptree.distanceById(idxs.first,idxs.second));
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
		numopt->setBounds(bc->getLeftBound(), bc->getRightBound() < 0 ? modelParams->divergenceBound : bc->getRightBound());


		start = chrono::system_clock::now();

		result = numopt->optimize() * -1.0;

		end = chrono::system_clock::now();
		elapsed_seconds = end-start;

		INFO("Estimated Divergence Time " << modelParams->getDivergenceTime(0));
		INFO("Real divergence time " << ptree.distanceById(idxs.first,idxs.second));
		double divdelta =  (modelParams->getDivergenceTime(0) - ptree.distanceById(idxs.first,idxs.second))/ptree.distanceById(idxs.first,idxs.second);

		if(std::abs(divdelta) < 0.025){
			INFO("Divergence delta % \e[38;5;46m " << divdelta*100 );
		}
		else if(std::abs(divdelta) < 0.05){
			INFO("Divergence delta % \e[38;5;76m " << divdelta*100 );
		}
		else if(std::abs(divdelta) < 0.1){
			INFO("Divergence delta % \e[38;5;106m " << divdelta*100 );
		}
		else if(std::abs(divdelta) < 0.25){
			INFO("Divergence delta % \e[38;5;136m " << divdelta*100 );
		}
		else if(std::abs(divdelta) < 0.5){
			INFO("Divergence delta % \e[38;5;166m " << divdelta*100 );
		}
		else{
			INFO("Divergence delta % \e[38;5;196m " << divdelta*100 );
		}



	    //INFO("Computation time " << elapsed_seconds.count() << "s\n");


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

