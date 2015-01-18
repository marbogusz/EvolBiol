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
	int cap;

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
	this->sampleAlignments();

	//SUBSTITUTION PART
	sme = new SubstitutionModelEstimator(inputSequences, model ,ot, rateCategories, alpha, estimateAlpha, tripletIdxs.size());

	auto it= alSamplesTriplet.rbegin();
	cap = Definitions::pathInformativeCount < alSamplesTriplet.size() ? Definitions::pathInformativeCount : alSamplesTriplet.size();

	for(unsigned int trp = 0; trp < tripletIdxs.size(); trp++)
	{
		auto it= alSamplesTriplet[trp].rbegin();
		cap = Definitions::pathInformativeCount < alSamplesTriplet[trp].size() ? Definitions::pathInformativeCount : alSamplesTriplet[trp].size();
		for (unsigned int i=0; i < cap; i++)
		{
			sme->addTriplet(it->second, trp, (double)exp((it->first)-(SamplesBranch1Lnls[trp]+SamplesBranch2Lnls[trp])));
			it++;
		}
	}
	sme->optimize();

	//INDEL PART
	ste = new StateTransitionEstimator(ot, 2*tripletIdxs.size());
	double tb1,tb2,tb3;

	for(unsigned int trp = 0; trp < tripletIdxs.size(); trp++)
	{
		tb1 = sme->getTripletDivergence(trp,0);
		tb2 = sme->getTripletDivergence(trp,1);
		tb3 = sme->getTripletDivergence(trp,2);

		ste->addTime(tb1+tb2,trp,0);
		ste->addTime(tb3+tb2,trp,1);

		auto itPr1=alSamplesBranch1[trp].rbegin();
		auto itPr2=alSamplesBranch2[trp].rbegin();

		cap = min(alSamplesBranch1[trp].size(),alSamplesBranch2[trp].size());
		cap = Definitions::pathInformativeCount < cap ? Definitions::pathInformativeCount : cap;

		for (unsigned int i=0; i < cap; i++)
		{
			ste->addPair((itPr1->second).first,(itPr1->second).second,trp,0,(double)exp((itPr1->first)-(SamplesBranch1Lnls[trp])));
			ste->addPair((itPr2->second).first,(itPr2->second).second,trp,1,(double)exp((itPr2->first)-(SamplesBranch2Lnls[trp])));
			itPr1++;
			itPr2++;
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

	substModel->summarize();
	indelModel->summarize();
}

void ModelEstimator::estimateParameters()
{

}

void ModelEstimator::calculateInitialHMMs(Definitions::ModelType model)
{
	DEBUG("EstimateTripleAligment");
	//amino acid mode
	bool aaMode = false;
	double tmpd;

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
	int cap;
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

	for (unsigned int i = 0; i < samplingHMMs.size(); i++)
	{
		for(ctr = 1; ctr < sampleCount; ctr++){
			auto branch1Sample = samplingHMMs[i].first->sampleAlignment(dict);
			auto branch2Sample = samplingHMMs[i].second->sampleAlignment(dict);

			alSamplesBranch1[i].insert(branch1Sample);
			alSamplesBranch2[i].insert(branch2Sample);

			SamplesBranch1Lnls[i] += branch1Sample.first;
			SamplesBranch2Lnls[i] += branch2Sample.first;
		}
/*
		auto itPr1=alSamplesBranch1[i].rbegin();
		auto itPr2=alSamplesBranch2[i].rbegin();

		cap = min(alSamplesBranch1[i].size(),alSamplesBranch2[i].size());
		cap = Definitions::pathInformativeCount < cap ? Definitions::pathInformativeCount : cap;

		for (unsigned int i=0; i < cap; i++)
		{
			auto tripletSample = tal->align(itPr1->second, itPr2->second);
			alSamplesTriplet[i].insert(make_pair(itPr1->first+itPr2->first, tripletSample));
			itPr1++;
			itPr2++;
		}
*/
		map<double, pair<string, string> > alSamples1;
		map<double, pair<string, string> > alSamples2;
		for(ctr = 0; ctr < 10000; ctr++){
			auto pr1 = samplingHMMs[i].first->sampleAlignment(inputSequences->getRawSequenceAt(tripletIdxs[i][0]), inputSequences->getRawSequenceAt(tripletIdxs[i][1]));
			auto v1= dict->translate(pr1.first);
			auto v2= dict->translate(pr1.second);
			double tlnl = samplingHMMs[i].first->getAlignmentLikelihood(v1,v2);
			alSamples1.insert(make_pair(tlnl, pr1));
		}

		DUMP("Top 100 alignments String ver");
			auto its=alSamples1.rbegin();
			int cap = 100 < alSamples1.size() ? 100 : alSamples1.size();

			for (int c=0;c<cap; c++)
			{
				stringstream ss1;
				DUMP(its->first);
				DUMP(its->second.first);
				DUMP(its->second.second);
				its++;
			}

		DUMP("Top 100 alignments SeqEL ver");
		auto it=alSamplesBranch1[i].rbegin();
		cap = 100 < alSamplesBranch1[i].size() ? 100 : alSamplesBranch1[i].size();

		for (int c=0;c<cap; c++)
		{
			stringstream ss1;
			DUMP(it->first);
			for (auto el : it->second.first)
				ss1 << el.getSymbol();
			DUMP(ss1.str());
			stringstream ss2;
			for (auto el : it->second.second)
				ss2 << el.getSymbol();
			DUMP(ss2.str());
			it++;
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
