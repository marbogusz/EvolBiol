/*
 * TripletModelEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "heuristics/TripletModelEstimator.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"


namespace EBC
{


TripletModelEstimator::BFGS::BFGS(TripletModelEstimator* enclosing, Definitions::OptimizationType ot) : optimizationType(ot)
{
	parent = enclosing;
	paramsCount = parent->modelParams->optParamCount();
	this->initParams.set_size(paramsCount);
	this->lowerBounds.set_size(paramsCount);
	this->upperBounds.set_size(paramsCount);

	parent->modelParams->toDlibVector(initParams,lowerBounds,upperBounds);

	//cerr << "DLIB optimizer init with " << paramsCount << " parameters" << endl;
}

TripletModelEstimator::BFGS::~BFGS()
{
}

double TripletModelEstimator::BFGS::objectiveFunction(const column_vector& bfgsParameters)
{
	this->parent->modelParams->fromDlibVector(bfgsParameters);
	return parent->runIteration();
}


const column_vector TripletModelEstimator::BFGS::objectiveFunctionDerivative(const column_vector& bfgsParameters)
{
	column_vector results(this->paramsCount);
	return results;
}


void TripletModelEstimator::BFGS::optimize()
{
	using std::placeholders::_1;
	std::function<double(const column_vector&)> f_objective= std::bind( &TripletModelEstimator::BFGS::objectiveFunction, this, _1 );
	double likelihood;

	switch(optimizationType)
	{
		case Definitions::OptimizationType::BFGS:
		{
			likelihood = dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
					dlib::objective_delta_stop_strategy(1e-8),
					f_objective,
					derivative(f_objective),
					initParams,
					lowerBounds,
					upperBounds);
			break;
		}
		case Definitions::OptimizationType::BOBYQA:
		{
			likelihood = dlib::find_min_bobyqa(f_objective, initParams, parent->modelParams->optParamCount()+4,
					lowerBounds,upperBounds, 0.05, 1e-7, 20000 );
			break;
		}
	}
	this->parent->modelParams->fromDlibVector(initParams);
	parent->modelParams->outputParameters();
	cout  << likelihood << "\n";

}


TripletModelEstimator::TripletModelEstimator(Sequences* inputSeqs, Definitions::ModelType model ,
		Definitions::OptimizationType ot,
		unsigned int rateCategories, double alpha, bool estimateAlpha) : inputSequences(inputSeqs), gammaRateCategories(rateCategories)
{
	maths = new Maths();
	dict = inputSequences->getDictionary();

	//Helper models
	//FIXME - get some static definitions or sth!!
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


	estimateSubstitutionParams = true;
	this->estimateAlpha = estimateAlpha;

	modelParams = new OptimizedModelParameters(substModel, indelModel,inputSequences->getSequenceCount(), pairCount, estimateSubstitutionParams,
			estimateIndelParams, estimateAlpha, userTime < 0, maths);

	if(!estimateIndelParams)
		modelParams->setUserIndelParams(indel_params);
	if(!estimateSubstitutionParams)
		modelParams->setUserSubstParams(subst_params);
	modelParams->setAlpha(alpha);

	if (userTime > 0)
	{
		//cerr << "User time " << userTime << endl;
		vector<double> times(pairCount);
		for (auto it = times.begin(); it < times.end(); it++)
		{
			*it = 2.0*(userTime/pairCount);
		}
		modelParams->setUserDivergenceParams(times);
	}


	bandFactor = bandPercentage;
	bandingEnabled = banding;

	EvolutionaryPairHMM* hmm;

	for(unsigned int i =0; i<pairCount; i++)
	{
		std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);

		if (at == Definitions::AlgorithmType::Viterbi)
		{
			hmm = hmms[i] = new ViterbiPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
				dict, model, banding, bandPercentage,rateCategories, maths, Definitions::DpMatrixType::Limited);
		}
		else if (at == Definitions::AlgorithmType::Forward)
		{
			hmm = hmms[i] = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
				dict, model, banding, bandPercentage,rateCategories, maths, Definitions::DpMatrixType::Limited);
		}
		else
		{
			throw HmmException("Wrong algorithm type - use either Forward or viterbi\n");
		}
		hmm->setModelFrequencies(inputSequences->getElementFrequencies());
	}

	bfgs = new BFGS(this,ot);
	bfgs->optimize();
}

TripletModelEstimator::~TripletModelEstimator()
{
	// TODO Auto-generated destructor stub
	delete bfgs;
	delete modelParams;
	//delete Y;
	//delete X;
	//delete M;
	//delete substModel;
	//delete indelModel;
    delete maths;
}

double TripletModelEstimator::runIteration()
{
	double result = 0;
	EvolutionaryPairHMM* hmm;

	//this->modelParams->outputParameters();
	//cerr << "iteration " << endl;
	for(unsigned int i =0; i<pairCount; i++)
	{
		hmm = hmms[i];
		hmm->setModelParameters(modelParams->getIndelParameters(),modelParams->getSubstParameters(),
				modelParams->getDivergenceTime(i), modelParams->getAlpha());
		result += hmm->runAlgorithm();
	}
	//cerr << result << endl;
	return result;
}

void TripletModelEstimator::outputDistanceMatrix(stringstream& ss)
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
