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
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletAligner.hpp"

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
		Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, bool estimateAlpha) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories), substs(3),
				gtree(inputSeqs), tst(gtree)
{
	maths = new Maths();
	dict = inputSequences->getDictionary();
	//we're running triplets, sequence count is 3!
	unsigned int sequence_count = 3;



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


	//TODO - maybe providing user times would be a good idea
	modelParams = new OptimizedModelParameters(substModel, NULL,3, 3, estimateSubstitutionParams,
			false, estimateAlpha, false, maths);

	//alpha is an initial alpha!!
	modelParams->setAlpha(alpha);


	//theoretically one can set the times from DM
	/*
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
	*/

	SubstitutionModelBase* smodel;
	for(unsigned int i =0; i<sequence_count; i++)
	{
		if (model == Definitions::ModelType::GTR)
		{
			smodel =  substs[i] = new GTRModel(dict, maths,gammaRateCategories);
		}
		else if (model == Definitions::ModelType::HKY85)
		{
			smodel =  substs[i] = new HKY85Model(dict, maths,gammaRateCategories);
		}
		else if (model == Definitions::ModelType::LG)
		{
			smodel =  substs[i] = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
		}
		smodel->setObservedFrequencies(inputSequences->getElementFrequencies());

		}

	vector<array<unsigned int, 3> > tripletIdxs = tst.sampleFromDM();
	//iterate over triplets and align

	for (auto idx : tripletIdxs)
	{
		TripletAligner tal(inputSequences, idx, gtree.getDistanceMatrix());
		tripleAlignments.push_back(tal.align());
	}

	//TripletSamplingTree tst(GuideTree(inputSequences));
	//triplets = tst.sampleFromDM();

	bfgs = new BFGS(this,ot);
	bfgs->optimize();
}

TripletModelEstimator::~TripletModelEstimator()
{
	// TODO Auto-generated destructor stub
	delete bfgs;
	delete modelParams;
    delete maths;


}

double TripletModelEstimator::runIteration()
{
	double result = 0;


	for (unsigned int i = 0; i< substs.size(); i++){
		substs[i]->setAlpha(modelParams->getAlpha());
		substs[i]->setParameters(modelParams->getSubstParameters());
		substs[i]->setTime(modelParams->getDivergenceTime(i));
		substs[i]->calculatePt();
		substs[i]->calculateSitePatterns();
	}

	//get triplet alignments



	for(auto it : this->tripleAlignments )
	{
		vector<SequenceElement>& s0 = it[0];
		vector<SequenceElement>& s1 = it[1];
		vector<SequenceElement>& s2 = it[2];

		//3 sequences
		for(int rt = 0; rt < dict->getAlphabetSize(); rt++)
		{
			for(int pos = 0; pos < s0.size(); pos++ )
			{
				//TODO
			}
		}

	}


	//cerr << result << endl;
	return result;
}

} /* namespace EBC */
