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
				inputSequences(inputSeqs), gammaRateCategories(rateCategories),
				gtree(inputSeqs), tst(gtree)
{
	maths = new Maths();
	dict = inputSequences->getDictionary();
	//we're running triplets, sequence count is 3!

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



	DEBUG("Gamma rates " << gammaRateCategories);

	estimateSubstitutionParams = true;
	this->estimateAlpha = estimateAlpha;


	//TODO - maybe providing user times would be a good idea
	//modelParams = new OptimizedModelParameters(substModel, NULL,3, 3, estimateSubstitutionParams,
	//		false, estimateAlpha, true, maths);


	//hack!

	//modelParams = new OptimizedModelParameters(substModel, NULL,3, 3, false,
	//			false, false, false, maths);


	//alpha
	//vector<double> userTimes = {0.76904, 0.12799, 0.40248};
	//vector<double> userTimes = {0.76904, 0.12799, 0.40248};
	//vector<double> userParameters = {1.60233,  0.38555,  0.45816,  0.51288,  0.48568};

	//noalpha
	//vector<double> userTimes = {0.51844, 0.14299, 0.29331};
	//vector<double> userParameters = {1.44141,  0.45891,  0.53792,  0.57142,  0.55240};

	//modelParams->setUserDivergenceParams(userTimes);
	//modelParams->setUserSubstParams(userParameters);
	//modelParams->setAlpha(1.43271);


	//hack ends

	//alpha is an initial alpha!!
	//modelParams->setAlpha(alpha);

	SubstitutionModelBase* smodel;


	DEBUG("About to sample some triplets");
	vector<array<unsigned int, 3> > tripletIdxs = tst.sampleFromDM();

	patterns.resize(tripletIdxs.size());
	ptMatrices.resize(tripletIdxs.size());

	DEBUG("Creating branch-specific substitution models");

	for(unsigned int i =0; i<ptMatrices.size(); i++)
	{
		ptMatrices[i][0] = new PMatrixTriple(substModel);
		ptMatrices[i][1] = new PMatrixTriple(substModel);
		ptMatrices[i][2] = new PMatrixTriple(substModel);
	}

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies(tripletIdxs[0]));

	for (auto idx : tripletIdxs)
	{
		TripletAligner tal(inputSequences, idx, gtree.getDistanceMatrix());
		tripleAlignments.push_back(tal.align());
	}

	for(int al = 0; al < tripleAlignments.size(); al++)
	{
		for(int pos = 0; pos < tripleAlignments[al][0].size(); pos++)
		{
			//count patterns
			patterns[al][{{tripleAlignments[al][0][pos].getMatrixIndex(), tripleAlignments[al][1][pos].getMatrixIndex(),tripleAlignments[al][2][pos].getMatrixIndex()}}]++;
		}
	}

	//TripletSamplingTree tst(GuideTree(inputSequences));
	//triplets = tst.sampleFromDM();

	modelParams = new OptimizedModelParameters(substModel, NULL,3, 3*tripleAlignments.size(), estimateSubstitutionParams,
			false, estimateAlpha, true, maths);

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


	double partial1;

	SubstitutionModelBase* smodel;
	SequenceElement* seTmp;

	modelParams->outputParameters();

	substModel->setAlpha(modelParams->getAlpha());
	substModel->setParameters(modelParams->getSubstParameters());
	substModel->calculateModel();


	substModel->summarize();

	for (unsigned int i = 0; i< ptMatrices.size(); i++)
	{
		for(unsigned int j=0;j<3;j++)
		{
			DEBUG("Triplet no, node no, time: " << i << " " << j << " "<< modelParams->getDivergenceTime(3*i +j));
			ptMatrices[i][j]->setTime(modelParams->getDivergenceTime(3*i +j));
			ptMatrices[i][j]->calculate();
		}
	}

	for(int al = 0; al < patterns.size(); al++)
	{
		for(auto it : patterns[al])
		{
			partial1 = 0;
			for(int rt = 0; rt < dict->getAlphabetSize(); rt++)
			{
				partial1 += ptMatrices[al][0]->getTripleSitePattern(rt,it.first, ptMatrices[al][1],ptMatrices[al][2]);
			}  //->getPairSitePattern(it.first[0],it.first[1]) * it.second;
			result += log(partial1)* it.second;
			//cout << "Pattern " << it.first[0] << " "  << it.first[1] << " " << it.first[2]  << "\t\t\t" << it.second << "\t\t" << log(partial1) << endl;
		}
	}
	DEBUG("lnl result:" << result);
	return result * -1.0;

}

} /* namespace EBC */
