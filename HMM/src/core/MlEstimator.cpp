/*
 * MlEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "core/MlEstimator.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"

namespace EBC
{


MlEstimator::BFGS::BFGS(MlEstimator* enclosing, Definitions::OptimizationType ot) : optimizationType(ot)
{
	parent = enclosing;
	paramsCount = parent->modelParams->optParamCount();
	this->initParams.set_size(paramsCount);
	this->lowerBounds.set_size(paramsCount);
	this->upperBounds.set_size(paramsCount);

	parent->modelParams->toDlibVector(initParams,lowerBounds,upperBounds);

	cerr << "DLIB optimizer init with " << paramsCount << " parameters" << endl;
}

MlEstimator::BFGS::~BFGS()
{
}

double MlEstimator::BFGS::objectiveFunction(const column_vector& bfgsParameters)
{
	this->parent->modelParams->fromDlibVector(bfgsParameters);
	return parent->runIteration();
}


const column_vector MlEstimator::BFGS::objectiveFunctionDerivative(const column_vector& bfgsParameters)
{
	column_vector results(this->paramsCount);
	return results;
}


void MlEstimator::BFGS::optimize()
{
	using std::placeholders::_1;
	std::function<double(const column_vector&)> f_objective= std::bind( &MlEstimator::BFGS::objectiveFunction, this, _1 );
	double likelihood;

	switch(optimizationType)
	{
		case Definitions::OptimizationType::BFGS:
		{
			dlib::objective_delta_stop_strategy strt(1e-7);
			strt.be_verbose();
			likelihood = dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
					strt,
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
	parent->modelParams->logParameters();
	//parent->substs[0]->summarize();
	cout  << likelihood << "\t";

}


MlEstimator::MlEstimator(Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params,   Definitions::OptimizationType ot,
		unsigned int rateCategories, double alpha, bool estimateAlpha, double userTime, bool alignViterbi) : inputSequences(inputSeqs), gammaRateCategories(rateCategories),
		pairCount(inputSequences->getPairCount()), hmms(pairCount), ptMatrices(pairCount, nullptr), useViterbi(alignViterbi), patterns(pairCount)

{
	maths = new Maths();
	dict = inputSequences->getDictionary();

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

	estimateSubstitutionParams = true;
	estimateIndelParams = false;
	this->estimateAlpha = estimateAlpha;

	//hack here!

	useViterbi = false;

	modelParams = new OptimizedModelParameters(substModel, NULL,2, 1, false,
				false, false, true, maths);
	vector<double> userTimes = {1.30288};
	vector<double> userParameters = {1.34838,  0.27850,  0.22691,  0.57785,  0.38311};

	modelParams->setUserDivergenceParams(userTimes);
	modelParams->setUserSubstParams(userParameters);
	modelParams->setAlpha(0.7);


	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());

	for(unsigned int i =0; i<pairCount; i++)
	{
		ptMatrices[i] = new PMatrixDouble(substModel);
		std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);
		vector<SequenceElement>  s1 = inputSequences->getSequencesAt(idxs.first);
		vector<SequenceElement>  s2 = inputSequences->getSequencesAt(idxs.second);
		for(int j = 0; j< s1.size(); j++)
		{
			patterns[i][{{s1[j].getMatrixIndex(),s2[j].getMatrixIndex()}}]++;
		}
	}

	/*
	SubstitutionModelBase* smodel;
			for(unsigned int i =0; i<1; i++)
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
				//smodel->setAlpha(modelParams->getAlpha());
				//smodel->setParameters(modelParams->getSubstParameters());
				//smodel->setTime(modelParams->getDivergenceTime(i));
				//smodel->calculatePt();
			}



	modelParams = new OptimizedModelParameters(substModel, indelModel,inputSequences->getSequenceCount(), pairCount, estimateSubstitutionParams,
			estimateIndelParams, estimateAlpha, true, maths);


	//if(indel_params.size() > 0)
	//	modelParams->setUserIndelParams(indel_params);
	//if(subst_params.size() > 0)
	//	modelParams->setUserSubstParams(subst_params);
	modelParams->setAlpha(alpha);

	//FIXME - provide divergence times later
	if (userTime > 0)
	{
		cerr << "User time " << userTime << endl;
		vector<double> times(pairCount);
		vector<double> proportions(inputSequences->getSequenceCount());
		double multiplier = 0.1;
		double relativeLen  = 0;

		//do proportions of lengths
		for (auto it = proportions.begin(); it < proportions.end(); it++)
		{
			*it = multiplier;
			relativeLen += multiplier;
			multiplier += 0.05;
		}

		//actual lengths = sum of pairs, divided by relative len times actual time

		unsigned int counter =0;
		unsigned int i1;
		unsigned int i2;
		double tmptime;

		for (auto it = times.begin(); it < times.end(); it++)
		{

			i1 = inputSequences->getPairOfSequenceIndices(counter).first;
			i2 = inputSequences->getPairOfSequenceIndices(counter).second;

			tmptime = ((proportions[i1] + proportions[i2]) / relativeLen) * userTime;

			*it = tmptime;
			counter++;
		}
		modelParams->setUserDivergenceParams(times);
	}


	if (useViterbi)
	{
		ViterbiPairHMM* hmm;

		for(unsigned int i =0; i<pairCount; i++)
		{
			std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);
			hmm = hmms[i] = new ViterbiPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					dict, model, false, 1 ,rateCategories, maths, Definitions::DpMatrixType::Full);
			hmm->setModelFrequencies(inputSequences->getElementFrequencies());
			hmm->setModelParameters(modelParams->getIndelParameters(),modelParams->getSubstParameters(),
				modelParams->getDivergenceTime(i), modelParams->getAlpha());

			//FIXME - HACK!!!!!!
			hmm->runAlgorithm();
		}

	}
	else
	{
		SubstitutionModelBase* smodel;
		for(unsigned int i =0; i<pairCount; i++)
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
			//smodel->setAlpha(modelParams->getAlpha());
			//smodel->setParameters(modelParams->getSubstParameters());
			//smodel->setTime(modelParams->getDivergenceTime(i));
			//smodel->calculatePt();
		}
	}
*/

	//do sitePatterns

	bfgs = new BFGS(this,ot);
	bfgs->optimize();
}

MlEstimator::~MlEstimator()
{
	delete bfgs;
	delete modelParams;
	//delete Y;
	//delete X;
	//delete M;
	//delete substModel;
	//delete indelModel;
    delete maths;
}

double MlEstimator::runIteration()
{
	double result = 0;
	if (useViterbi)
	{
		ViterbiPairHMM* hmm;
	//

		for(unsigned int i =0; i<pairCount; i++)
		{
			hmm = hmms[i];
			hmm->setDivergenceTimeAndCalculateModels(modelParams->getDivergenceTime(i));
			result += hmm->getViterbiSubstitutionLikelihood();
		}
	}
	else
	{
		//this->modelParams->outputParameters();
		substModel->setAlpha(modelParams->getAlpha());
		substModel->setParameters(modelParams->getSubstParameters());
		substModel->calculateModel();

		substModel->summarize();

		for(unsigned int i =0; i<pairCount; i++)
		{
			//this calculates the matrix(matrices for a gamma model)
			ptMatrices[i]->setTime(modelParams->getDivergenceTime(i));
			ptMatrices[i]->calculate();

			//go through the map of patterns!
			for(auto it : patterns[i])
			{
				result += ptMatrices[i]->getPairSitePattern(it.first[0],it.first[1]) * it.second;
			}

			ptMatrices[i]->summarize();

		}
	}

	cerr << result << endl;
	return result * -1.0;
}


} /* namespace EBC */
