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
#include "core/PhylogeneticTree.hpp"

namespace EBC
{

MlEstimator::MlEstimator(Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params,   Definitions::OptimizationType ot,
		unsigned int rateCategories, double alpha, bool estimateAlpha, vector<double> userTime, bool alignViterbi) : inputSequences(inputSeqs), gammaRateCategories(rateCategories),
		pairCount(inputSequences->getPairCount()), hmms(pairCount), ptMatrices(pairCount, nullptr), useViterbi(alignViterbi), patterns(pairCount),
		optTimes(pairCount)

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

	substModel->setAlpha(alpha);

	/*
	estimateSubstitutionParams = true;
	if (subst_params.size() > 0){
		estimateSubstitutionParams = false;
	}
	*/

	estimateIndelParams = false;
	estimateSubstitutionParams = false;
	//can't estimate alpha!
	//unless done before using triplets!

	this->estimateAlpha = false;

	//hack here!

	useViterbi = false;

	/*
	modelParams = new OptimizedModelParameters(substModel, NULL,inputSequences->getSequenceCount(), pairCount, estimateSubstitutionParams,
			estimateIndelParams, estimateAlpha, true, maths);

	*/
	modelParams = new OptimizedModelParameters(substModel, NULL,2, 1, false,
				false, false, true, maths);

	if (estimateSubstitutionParams)
		modelParams->generateInitialSubstitutionParameters();
	//modelParams->generateInitialDistanceParameters();

	//modelParams->setUserDivergenceParams(userTime);

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
	if (!estimateSubstitutionParams){
		substModel->setParameters(subst_params);
		substModel->calculateModel();
	}

	//bfgs = new BFGS(this,ot);
	numopt = new BrentOptimizer(modelParams, NULL);

	ifstream myfile ("tree.sim_1");
	string newick;
	if (myfile.is_open())
	{
		getline (myfile,newick);
		myfile.close();
	}

	PhylogeneticTree ptree(this->inputSequences);
	ptree.fromNewick(newick);

	for(unsigned int i =0; i<pairCount; i++)
	{
		ptMatrices[i] = new PMatrixDouble(substModel);
		std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);
		vector<SequenceElement*>*  s1 = inputSequences->getSequencesAt(idxs.first);
		vector<SequenceElement*>*  s2 = inputSequences->getSequencesAt(idxs.second);
		for(int j = 0; j< s1->size(); j++)
		{
			patterns[i][{{(*s1)[j]->getMatrixIndex(),(*s2)[j]->getMatrixIndex()}}]++;
		}
		currentPair = i;
		modelParams->setUserDivergenceParams({userTime[i]});
		INFO("Running pairwise calculator for sequence id " << idxs.first << " and " << idxs.second);

		numopt->setTarget(this);
		numopt->setAccuracy(1e-6);
		numopt->setBounds(0.0000001, modelParams->divergenceBound);

		numopt->optimize() * -1.0;

		//bfgs->optimize();
		optTimes[i] = modelParams->getDivergenceTime(0);



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

	}
}

MlEstimator::~MlEstimator()
{
	//delete bfgs;
	delete numopt;
	delete modelParams;
    delete maths;
}

double MlEstimator::runIteration()
{
	double result = 0;
	if (estimateSubstitutionParams){
		substModel->setParameters(modelParams->getSubstParameters());
		substModel->calculateModel();
	}
	//for(unsigned int i =0; i<pairCount; i++)
	//{
		//this calculates the matrix(matrices for a gamma model)
		ptMatrices[currentPair]->setTime(modelParams->getDivergenceTime(0));
		ptMatrices[currentPair]->calculate();

		//go through the map of patterns!
		for(auto it : patterns[currentPair])
		{
			result += ptMatrices[currentPair]->getPairSitePattern(it.first[0],it.first[1]) * it.second;
		}

	//}

	//cerr << result << endl;
	return result * -1.0;
}


} /* namespace EBC */
