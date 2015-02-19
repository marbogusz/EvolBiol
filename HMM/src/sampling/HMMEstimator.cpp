/*
 * HMMEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "sampling/HMMEstimator.hpp"
#include "sampling/PathExtractor.hpp"
//#include "sampling/ExpTester.hpp"

#include <chrono>
#include <array>
#include <map>

using namespace std;

namespace EBC
{

HMMEstimator::HMMEstimator(Sequences* inputSeqs, Definitions::ModelType model ,
		Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, bool estimateAlpha) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories),
				gtree(new GuideTree(inputSeqs)), tst(*gtree)
{
	DEBUG("HMM Estimator starting");
	DEBUG("About to sample some triplets");
	DEBUG("Sampling triplets of sequences for gamma shape parameter estimation");
	
	maths = new Maths();
	dict = inputSequences->getDictionary();
	tal = new TripletAligner (inputSequences, gtree->getDistanceMatrix());

	tripletIdxs = tst.sampleFromTree();

	tripletIdxsSize = tripletIdxs.size();

	//2 pairs for every triplet!
	sampleWorkers.reserve(tripletIdxsSize*2);

	//FIXME - release the memory!!!! - delete pair objects and vectors (arrays) of SeqEls

    //chrono::time_point<chrono::system_clock> start, end;
    //start = chrono::system_clock::now();

	if (model == Definitions::ModelType::GTR)
	{
		substModel = new GTRModel(dict, maths,gammaRateCategories);
		DUMP("SME: Creating new GTR model");
	}
	else if (model == Definitions::ModelType::HKY85)
	{
		substModel = new HKY85Model(dict, maths,gammaRateCategories);
		DUMP("SME: Creating new HKY model");
	}
	else if (model == Definitions::ModelType::LG)
	{
			substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
			DUMP("SME: Creating new LG model");
	}

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());

	indelModel = new NegativeBinomialGapModel();

	modelParams = new OptimizedModelParameters(substModel, indelModel, 0, 0, false, true, estimateAlpha, false, maths);

	bfgs = new Optimizer(modelParams, this, Definitions::OptimizationType::BFGS);

	this->calculateInitialPairs(model);
	this->optimise();

	//for (auto &worker : sampleWorkers)
	//{
	//	worker.doExtraStuff();
	//}

	//FIXME - case study start;
	/*
	this->calculateInitialPairs(model);

	ForwardPairHMM phmm(inputSequences->getSequencesAt(0), inputSequences->getSequencesAt(1),
			substModel, indelModel, Definitions::DpMatrixType::Full, nullptr,true);
	phmm.setDivergenceTimeAndCalculateModels(0.5);

	double fws, vts;

	fws = phmm.runAlgorithm();

	cerr << "FWD lnl " << fws << endl;

	PathExtractor ex(inputSequences->getRawSequenceAt(0), inputSequences->getRawSequenceAt(1), inputSequences->getSequencesAt(0), inputSequences->getSequencesAt(1));
	ex.genAllSeqs(vector<SequenceElement*>(),vector<SequenceElement*>(),dict->getSequenceElement('-'),0,0,3);

	vts = ex.getSumLnl(&phmm, maths) * -1.0;

	cerr << "SUM lnl " << vts << endl;

	cerr << "Difference " << fws-vts << endl;


	//FIXME - case study end
*/
	//cout << "Times : " << endl;
	//for (auto &worker : sampleWorkers)
	//{
	//	cout << worker.getDivergence() << endl;
//	}

	//this->runIteration();
/*
    ExpTester tstr;

    chrono::time_point<chrono::system_clock> start, end, endD;
    start = chrono::system_clock::now();

    cerr << tstr.testExpFloat(100000000) << endl;

	end = chrono::system_clock::now();

    cerr << tstr.testExpDouble(100000000) << endl;

    endD = chrono::system_clock::now();


    chrono::duration<double> elapsed_seconds = end-start;
    chrono::duration<double> elapsed_secondsD = endD-end;

    cerr << "Exp tester seconds float : " << elapsed_seconds.count() << endl;
    cerr << "Exp tester seconds double : " << elapsed_secondsD.count() << endl;
*/
    //INFO("Model Estimator elapsed time: " << elapsed_seconds.count() << " seconds");

	//substModel->summarize();

}

void HMMEstimator::calculateInitialPairs(Definitions::ModelType model)
{
	DEBUG("HMM Estimator - create initial pairs");
	//amino acid mode
	bool aaMode = false;
	double tmpd;

	double initAlpha = 0.75;
	double initKappa = 2.5;
	double initLambda = 0.025;
	double initEpsilon = 0.5;
	//k-mers tend to underestimate the distances;
	double initTimeModifier = 1.5;

	//if rateCat =  1 alpha does not matter.
	substModel->setAlpha(initAlpha);
	modelParams->setAlpha(initAlpha);

	if (model == Definitions::ModelType::GTR)
	{
		substModel->setParameters({1.0,1.0/initKappa,1.0/initKappa,1.0/initKappa,1.0/initKappa});
		modelParams->setUserSubstParams({1.0,1.0/initKappa,1.0/initKappa,1.0/initKappa,1.0/initKappa});
	}
	else if (model == Definitions::ModelType::HKY85)
	{
		substModel->setParameters({initKappa});
		modelParams->setUserSubstParams({initKappa});
	}
	//no else - no need to set anything for AAs

	modelParams->setUserIndelParams({initLambda, initEpsilon});

	substModel->calculateModel();

	indelModel->setParameters({initLambda,initEpsilon});

	for (int i = 0; i < tripletIdxs.size(); i++)
	{

		sampleWorkers.emplace_back(inputSequences->getSequencesAt(tripletIdxs[i][0]), inputSequences->getSequencesAt(tripletIdxs[i][1]),
				substModel, indelModel, gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][0],tripletIdxs[i][1]) * initTimeModifier);
		sampleWorkers.emplace_back(inputSequences->getSequencesAt(tripletIdxs[i][1]), inputSequences->getSequencesAt(tripletIdxs[i][2]),
			substModel, indelModel, gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][1],tripletIdxs[i][2]) * initTimeModifier);

		//sampleWorkers.emplace_back(inputSequences->getSequencesAt(0), inputSequences->getSequencesAt(1),
		//				substModel, indelModel, gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][0],tripletIdxs[i][1]) * initTimeModifier);
		//sampleWorkers.emplace_back(inputSequences->getSequencesAt(1), inputSequences->getSequencesAt(2),
		//	substModel, indelModel, gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][1],tripletIdxs[i][2]) * initTimeModifier);
	}

}

double EBC::HMMEstimator::runIteration() {
	double lnl = 0;

	substModel->setAlpha(modelParams->getAlpha());
	substModel->setParameters(modelParams->getSubstParameters());
	substModel->calculateModel();
	indelModel->setParameters(modelParams->getIndelParameters());

	for (auto &worker : sampleWorkers)
	{
		lnl += worker.optimiseDivergenceTime();
	}

	//cerr << "lnL HMMEST " <<  lnl << "\t" << modelParams->getIndelParameters()[0] << "\t" << modelParams->getIndelParameters()[1] << endl;

	return lnl;

}

void EBC::HMMEstimator::optimise() {
	bfgs->optimize();
	//run again!

	/*
	cerr << modelParams->getIndelParameters()[0] << "\t" << modelParams->getIndelParameters()[1] << endl;

	substModel->setAlpha(modelParams->getAlpha());
	substModel->setParameters(modelParams->getSubstParameters());
	substModel->calculateModel();
	indelModel->setParameters(modelParams->getIndelParameters());

	for (auto& sw : sampleWorkers){
		//cerr << "Divergence " << sw.getDivergence() << endl;
		sw.reSample();
	}
	//DUMP("HMM estimator Second Pass ^^^^");
	bfgs->optimize();
	//cout << "Divergence begin " << sampleWorkers.begin()->getDivergence() << endl;
*/
}


HMMEstimator::~HMMEstimator()
{
    delete maths;
    delete gtree;
    delete tal;
}

vector<double> HMMEstimator::getSubstitutionParameters()
{
	return this->modelParams->getSubstParameters();
}

vector<double> HMMEstimator::getIndelParameters()
{
	return this->modelParams->getIndelParameters();
}


double HMMEstimator::getAlpha()
{
	return this->modelParams->getAlpha();
}

} /* namespace EBC */


