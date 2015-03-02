/*
 * ModelEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "core/Dictionary.hpp"
#include "heuristics/ModelEstimator.hpp"
#include <chrono>
#include <array>
#include <map>

using namespace std;

namespace EBC
{

ModelEstimator::ModelEstimator(Sequences* inputSeqs, Definitions::ModelType model ,
		Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, bool estimateAlpha) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories),
				gtree(new GuideTree(inputSeqs)), tst(*gtree), estAlpha(estimateAlpha), estIndel(true), estSubst(true)
{

	DEBUG("About to sample some triplets");
	DEBUG("Sampling triplets of sequences for gamma shape parameter estimation");
	
	maths = new Maths();
	dict = inputSequences->getDictionary();
	tal = new TripletAligner (inputSequences, gtree->getDistanceMatrix());

	tripletIdxs = tst.sampleFromTree();

	tripletIdxsSize = tripletIdxs.size();

	//FIXME - release the memory!!!! - delete pair objects and vectors (arrays) of SeqEls

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

	this->calculateInitialHMMs(model);

	//delete initial simplified model
	delete substModel;

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

	sme = new SubstitutionModelEstimator(inputSequences, substModel ,ot, rateCategories, alpha, estimateAlpha, tripletIdxs.size());
	ste = new StateTransitionEstimator(indelModel, ot, 2*tripletIdxs.size(), dict->getGapID());

	end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;

    INFO("Model Estimator elapsed time: " << elapsed_seconds.count() << " seconds");
    cerr <<  "|||||||||||| Model Estimator elapsed time: " << elapsed_seconds.count() << "s ||||||||||||||\n";

	//substModel->summarize();

}

void ModelEstimator::estimateParameters()
{
	DUMP("Model Estimator estimate parameters iteration");

	int cap;
	//cap = Definitions::pathInformativeCount < alSamplesTriplet->size() ? Definitions::pathInformativeCount : alSamplesTriplet->size();

	for(unsigned int trp = 0; trp < tripletIdxsSize; trp++)
	{
		auto it= alSamplesTriplet[trp].begin();
		cap = Definitions::pathInformativeCount < alSamplesTriplet[trp].size() ? Definitions::pathInformativeCount : alSamplesTriplet[trp].size();
		for (unsigned int i=0; i < cap; i++)
		{
			DUMP("About to add triplet with lnl " << it->first << "\ttotal lnl " << (SamplesBranch1Lnls[trp]+SamplesBranch2Lnls[trp]));
			sme->addTriplet(it->second, trp, (double)exp((it->first)-(SamplesBranch1Lnls[trp]+SamplesBranch2Lnls[trp])));
			it++;
		}
	}
	sme->optimize();

	//INDEL PART
	double tb1,tb2,tb3;

	for(unsigned int trp = 0; trp < tripletIdxsSize; trp++)
	{
		tb1 = sme->getTripletDivergence(trp,0);
		tb2 = sme->getTripletDivergence(trp,1);
		tb3 = sme->getTripletDivergence(trp,2);

		ste->addTime(tb1+tb2,trp,0);
		ste->addTime(tb3+tb2,trp,1);

		auto itPr1=alSamplesBranch1[trp].begin();
		auto itPr2=alSamplesBranch2[trp].begin();

		cap = min(alSamplesBranch1[trp].size(),alSamplesBranch2[trp].size());
		cap = Definitions::pathInformativeCount < cap ? Definitions::pathInformativeCount : cap;

		for (unsigned int i=0; i < cap; i++)
		{
			ste->addPair((itPr1->second)->first,(itPr1->second)->second,trp,0,(double)exp((itPr1->first)-(SamplesBranch1Lnls[trp])));
			ste->addPair((itPr2->second)->first,(itPr2->second)->second,trp,1,(double)exp((itPr2->first)-(SamplesBranch2Lnls[trp])));
			itPr1++;
			itPr2++;
		}
	}
	ste->optimize();

	substModel->summarize();
	indelModel->summarize();

	sme->clean();
	ste->clean();
}

void ModelEstimator::calculateInitialHMMs(Definitions::ModelType model)
{
	DEBUG("EstimateTripleAligment");
	//amino acid mode
	bool aaMode = false;
	double tmpd;

	double initAlpha = 0.75;
	double initKappa = 2.0;
	double initLambda = 0.05;
	double initEpsilon = 0.5;
	//k-mers tend to underestimate the distances;
	double initTimeModifier = 1.5;

	vector<double> timeModifiers = {0.75, 1.0, 1.5};
	vector<double> lambdas = {0.03, 0.07};
	vector<double> alphas = {0.5,1.0};

	ForwardPairHMM *f1, *f2;

	if (model == Definitions::ModelType::GTR || model == Definitions::ModelType::HKY85)
	{
			//FIXME - make model idiotproof by checking if parameters are set;
			DEBUG("Setting HKY85");
			this->substModel = new HKY85Model(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::LG)
	{
			aaMode = true;
			substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
	}


	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
	//alpha setting will have no effect if we're dealing with 1 rate category
	substModel->setAlpha(initAlpha);
	if (!aaMode)
		substModel->setParameters({initKappa});
	//substModel->calculateModel();

	indelModel = new NegativeBinomialGapModel();
	indelModel->setParameters({initLambda,initEpsilon});

	vector<pair<Band*, Band*> > bandPairs(tripletIdxs.size());

	vector<array<vector<SequenceElement*>*,3> > seqsA(tripletIdxs.size());
	vector<array<double,3> > distancesA(tripletIdxs.size());

	vector<array<ForwardPairHMM*,2> > fwdHMMs;

	for (int i = 0; i < tripletIdxs.size(); i++)
	{


		seqsA[i][0] = inputSequences->getSequencesAt(tripletIdxs[i][0]);
		DUMP("Triplet " << i << " sequence 1:");
		DUMP(inputSequences->getRawSequenceAt(tripletIdxs[i][0]));
		seqsA[i][1] = inputSequences->getSequencesAt(tripletIdxs[i][1]);
		DUMP("Triplet " << i << " sequence 2:");
		DUMP(inputSequences->getRawSequenceAt(tripletIdxs[i][1]));
		seqsA[i][2] = inputSequences->getSequencesAt(tripletIdxs[i][2]);
		DUMP("Triplet " << i << " sequence 3:");
		DUMP(inputSequences->getRawSequenceAt(tripletIdxs[i][2]));

		unsigned int len1 = seqsA[i][0]->size();
		unsigned int len2 = seqsA[i][1]->size();
		unsigned int len3 = seqsA[i][2]->size();

		//0-1
		tmpd = gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][0],tripletIdxs[i][1]);
		DUMP("Triplet " << i << " guide distance between seq 1 and 2 " << tmpd);
		distancesA[i][0] = tmpd;
		//1-2
		tmpd = gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][1],tripletIdxs[i][2]);
		DUMP("Triplet " << i << " guide distance between seq 2 and 3 " << tmpd);
		distancesA[i][1] = tmpd;
		//0-2
		tmpd = gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][0],tripletIdxs[i][2]);
		DUMP("Triplet " << i << " guide distance between seq 1 and 3 " << tmpd);
		distancesA[i][2] = tmpd;
		//bandPairs[i] = make_pair(new Band(len1,len2),new Band(len2,len3));
		bandPairs[i] = make_pair(nullptr,nullptr);

		fwdHMMs[i][0] = new ForwardPairHMM(seqsA[i][0],seqsA[i][1], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first,true);
		fwdHMMs[i][1] = new ForwardPairHMM(seqsA[i][1],seqsA[i][2], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].second,true);
	}
	double bestLnl = Definitions::minMatrixLikelihood;
	double currentLnl;
	double bestA, bestL, bestTm;

	for (auto tm : timeModifiers){
		for(auto l : lambdas){
			indelModel->setParameters({l,initEpsilon});
			for(auto a : alphas){
				substModel->setAlpha(a);
				if (estAlpha)
					substModel->calculateModel();
				currentLnl = 0;
				for (int i = 0; i < tripletIdxs.size(); i++){
					f1 = fwdHMMs[i][0];
					f2 = fwdHMMs[i][1];

					f1->setDivergenceTimeAndCalculateModels(distancesA[i][0]*tm);
					f2->setDivergenceTimeAndCalculateModels(distancesA[i][1]*tm);

					currentLnl += (f1->runAlgorithm() + f2->runAlgorithm()) * -1.0;
					if (currentLnl > bestLnl){
						bestLnl = currentLnl;
						bestA = a;
						bestL = l;
						bestTm = tm;

					}
				}

			}
		}
	}
	substModel->setAlpha(bestA);
	if (estAlpha)
		substModel->calculateModel();
	indelModel->setParameters({bestL,initEpsilon});

	//Fwd + bwd + MPD
	for (int i = 0; i < tripletIdxs.size(); i++){
		f1 = fwdHMMs[i][0];
		f2 = fwdHMMs[i][1];
		f1->setDivergenceTimeAndCalculateModels(distancesA[i][0]*bestTm);
		f2->setDivergenceTimeAndCalculateModels(distancesA[i][1]*bestTm);

		BackwardPairHMM b1(seqsA[i][0],seqsA[i][1], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first);
		BackwardPairHMM b2(seqsA[i][1],seqsA[i][2], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].second);

		b1.setDivergenceTimeAndCalculateModels(distancesA[i][0]*bestTm);
		b2.setDivergenceTimeAndCalculateModels(distancesA[i][1]*bestTm);

		b1.runAlgorithm();
		b2.runAlgorithm();

		b1.calculatePosteriors(f1);
		b2.calculatePosteriors(f2);

		b1.calculateMaximumPosteriorMatrix();
		b2.calculateMaximumPosteriorMatrix();



		delete f1;
		delete f2;

	}


	//align triplets based on  posteriors


}




ModelEstimator::~ModelEstimator()
{
	for (auto pr : samplingHMMs)
	{
		delete pr.first;
		delete pr.second;
	}

	this->alSamplesBranch1.resize(tripletIdxs.size());
	this->alSamplesBranch2.resize(tripletIdxs.size());
	this->alSamplesTriplet.resize(tripletIdxs.size());

//FIXME - clean up!

    delete maths;
    delete sme;
    delete ste;
    delete gtree;
    delete tal;

}

vector<double> ModelEstimator::getSubstitutionParameters()
{
	return this->sme->getModelParams()->getSubstParameters();
}

vector<double> ModelEstimator::getIndelParameters()
{
	return this->ste->getModelParams()->getIndelParameters();
}

double ModelEstimator::getAlpha()
{
	return this->sme->getModelParams()->getAlpha();
}

} /* namespace EBC */
