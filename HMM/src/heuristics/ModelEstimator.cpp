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
				gtree(new GuideTree(inputSeqs)), tst(*gtree)
{
	DEBUG("About to sample some triplets");
	DEBUG("Sampling triplets of sequences for gamma shape parameter estimation");
	
	maths = new Maths();
	dict = inputSequences->getDictionary();
	tal = new TripletAligner (inputSequences, gtree->getDistanceMatrix());

	tripletIdxs = tst.sampleFromTree();

	this->alSamplesBranch1.resize(tripletIdxs.size());
	this->alSamplesBranch2.resize(tripletIdxs.size());

	this->alSamplesTriplet.resize(tripletIdxs.size());


	this->SamplesBranch1Lnls.resize(tripletIdxs.size());
	this->SamplesBranch2Lnls.resize(tripletIdxs.size());
	this->SamplesTripletLnls.resize(tripletIdxs.size());

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

	this->calculateInitialHMMs(model);

	this
	//1 pair for now
	ste = new StateTransitionEstimator(ot, 1);

	ste->addTime(1.0,0,0);

	auto it=alSamples.rbegin();

	int cap = Definitions::pathInformativeCount < alSamples.size() ? Definitions::pathInformativeCount : alSamples.size();

	for (unsigned int i=0; i < cap; i++)
	{
		DUMP("lnl: " << it->first << "\t total lnl: " << totalSampleLnl);
		ste->addPair((it->second).first,(it->second).second,0,0,(double)exp((it->first)-(this->totalSampleLnl)));
		it++;
	}


	ste->optimize();

	DEBUG("Re-estimating model parameters");
	//make another pass

	indelModel =  ste->getIndelModel();
	//indelModel->summarize();

/*
	sme = new SubstitutionModelEstimator(inputSeqs, model ,ot, rateCategories, alpha, estimateAlpha, tripletIdxs.size());

	for(unsigned int al = 0; al < tripletIdxs.size(); al++)
	{
		for (unsigned int i=0; i < Definitions::modelEstimatorPathSamples; i++)
			sme->addTriplet(sampleTripleAlignment(al),al);
	}
	sme->optimize();

	double tb1,tb2,tb3;

	ste = new StateTransitionEstimator(ot, tripletIdxs.size()*2);

	for(int al = 0; al < tripletIdxs.size(); al++)
	{
		tb1 = sme->getModelParams()->getDivergenceTime(al*3);
		tb2 = sme->getModelParams()->getDivergenceTime((al*3)+1);
		tb3 = sme->getModelParams()->getDivergenceTime((al*3)+2);

		ste->addTime(tb1+tb2,al,0);
		ste->addTime(tb2+tb3,al,1);

		DEBUG("Branch divergence used by STE : 1:" << tb1 << " 2: " << tb2 << " 3: " << tb3);
		for (unsigned int i=0; i < Definitions::modelEstimatorPathSamples; i++)
		{
			ste->addPair(smpldPairs[al][i][0],smpldPairs[al][i][1],al,0);
			ste->addPair(smpldPairs[al][i][2],smpldPairs[al][i][3],al,1);
		}
	}

	ste->optimize();

	DEBUG("Re-estimating model parameters");
	//make another pass

	indelModel =  ste->getIndelModel();
	substModel =  sme->getSubstModel();

	end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;

    INFO("Model Estimator elapsed time: " << elapsed_seconds.count() << " seconds");
    cerr <<  "|||||||||||| Model Estimator elapsed time: " << elapsed_seconds.count() << "s ||||||||||||||\n";

	//substModel->summarize();
	//indelModel->summarize();
	//substModel->summarize();
	//we have new alignments!
	//re-estimate
	//do Viterbi using the estimates
	//construct triplets
	 * 
	 */
}

void ModelEstimator::estimateParameters()
{
	sme = new SubstitutionModelEstimator(inputSeqs, model ,ot, rateCategories, alpha, estimateAlpha, tripletIdxs.size());

	for(unsigned int al = 0; al < tripletIdxs.size(); al++)
	{
		for (unsigned int i=0; i < Definitions::modelEstimatorPathSamples; i++)
			sme->addTriplet(sampleTripleAlignment(al),al);
	}
	sme->optimize();

	double tb1,tb2,tb3;
}

void ModelEstimator::calculateInitialHMMs(Definitions::ModelType model)
{
	DEBUG("EstimateTripleAligment");
	//amino acid mode
	bool aaMode = false;

	double initAlpha = 0.75;
	double initKappa = 2.5;
	double initLambda = 0.02;
	double initEpsilon = 0.5;
	//k-mers tend to underestimate the distances;
	double initTimeModifier = 1.5;

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
	substModel->calculateModel();

	indelModel = new NegativeBinomialGapModel();

	vector<pair<Band*, Band*> > bandPairs(tripletIdxs.size());
	vector<array<vector<SequenceElement>,3> > seqsA(tripletIdxs.size());
	vector<array<double,3> > distancesA(tripletIdxs.size());

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

		unsigned int len1 = seqsA[i][0].size();
		unsigned int len2 = seqsA[i][1].size();
		unsigned int len3 = seqsA[i][2].size();

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


		f1 = new ForwardPairHMM(seqsA[i][0],seqsA[i][1], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first,true);
		f2 = new ForwardPairHMM(seqsA[i][1],seqsA[i][2], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].second,true);

		this->samplingHMMs.push_back(make_pair(f1,f2));

	}
}



void ModelEstimator::sampleAlignments()
{
	//FIXME  - check if maps are empty, zero if not
	//Same with likelihood values

	int sampleCount = Definitions::pathSampleCount;
	int analysisCount = Definitions::pathInformativeCount;
	double totalBranchSampleLnl = 0;
	double totalPairSampleLnl = 0;
	double lnl;
	int ctr;
	pair<string, string> smplPr;
	pair<double, pair<vector<SequenceElement>, vector<SequenceElement> > > pr;
	totalSampleLnl = Definitions::minMatrixLikelihood;

	//do the first sample

	for (int i = 0; i < samplingHMMs.size(); i++)
	{
		for(ctr = 1; ctr < sampleCount; ctr++){
			auto branch1Sample = samplingHMMs[i].first->sampleAlignment();
			auto branch2Sample = samplingHMMs[i].second->sampleAlignment();

			alSamplesBranch1[i].insert(branch1Sample);
			alSamplesBranch1[i].insert(branch2Sample);

			SamplesBranch1Lnls[i] += branch1Sample.first;
			SamplesBranch2Lnls[i] += branch2Sample.first;

			auto tripletSample = tal->align(branch1Sample.second, branch2Sample.second)
			alSamplesTriplet[i].insert(make_pair(branch1Sample.first+branch2Sample.first, tripletSample));
		}
	}
}

ModelEstimator::~ModelEstimator()
{
	for (auto pr : samplingHMMs)
	{
		delete pr.first;
		delete pr.second;
	}
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
