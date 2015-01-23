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

	tripletIdxsSize = tripletIdxs.size();

	//FIXME - release the memory!!!! - delete pair objects and vectors (arrays) of SeqEls

	this->alSamplesBranch1.resize(tripletIdxsSize);
	this->alSamplesBranch2.resize(tripletIdxsSize);
	this->alSamplesTriplet.resize(tripletIdxsSize);

	this->SamplesBranch1Lnls.resize(tripletIdxsSize);
	this->SamplesBranch2Lnls.resize(tripletIdxsSize);
	this->SamplesTripletLnls.resize(tripletIdxsSize);

	for (unsigned int i = 0; i < tripletIdxsSize; i++){
		alSamplesBranch1[i].reserve(Definitions::pathSampleCount);
		alSamplesBranch2[i].reserve(Definitions::pathSampleCount);
		alSamplesTriplet[i].reserve(Definitions::pathInformativeCount);
	}

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

	this->calculateInitialHMMs(model);
	this->sampleAlignments();

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

	estimateParameters();
	rescoreSamples();
	estimateParameters();
	rescoreSamples();
	estimateParameters();
	rescoreSamples();
	estimateParameters();




	//Re-score



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

void ModelEstimator::rescoreSamples()
{
	DUMP("Model Estimator rescoring samples");

	double totalSampleLnlB1;
	double totalSampleLnlB2;

	double lnlB1, lnlB2;

	double tb1,tb2,tb3;

	//re-calculate transitions etc.
	for (unsigned int trpIdx = 0; trpIdx < tripletIdxs.size(); trpIdx++)
	{
		tb1 = sme->getTripletDivergence(trpIdx,0);
		tb2 = sme->getTripletDivergence(trpIdx,1);
		tb3 = sme->getTripletDivergence(trpIdx,2);
		samplingHMMs[trpIdx].first->setDivergenceTimeAndCalculateModels(tb1+tb2);
		samplingHMMs[trpIdx].first->setDivergenceTimeAndCalculateModels(tb3+tb2);
	}

	for (unsigned int trpIdx = 0; trpIdx < tripletIdxs.size(); trpIdx++)
	{
		totalSampleLnlB1 = totalSampleLnlB2 = Definitions::minMatrixLikelihood;

		auto itPr1=alSamplesBranch1[trpIdx].begin();
		auto itPr2=alSamplesBranch2[trpIdx].begin();
		auto itTrp=alSamplesTriplet[trpIdx].begin();

		for (unsigned int pos = 0; pos < alSamplesTriplet[trpIdx].size(); pos ++)
		{
			lnlB1 = samplingHMMs[trpIdx].first->getAlignmentLikelihood(itPr1->second->first, itPr1->second->second,dict);
			lnlB2 = samplingHMMs[trpIdx].second->getAlignmentLikelihood(itPr2->second->first, itPr2->second->second,dict);

			itPr1->first = lnlB1;
			itPr2->first = lnlB2;
			itTrp->first = lnlB1+lnlB2;


			totalSampleLnlB1 = maths->logSum(totalSampleLnlB1, lnlB1);
			totalSampleLnlB2 = maths->logSum(totalSampleLnlB2, lnlB2);
			itPr1++;
			itPr2++;
			itTrp++;
		}

		SamplesBranch1Lnls[trpIdx] = totalSampleLnlB1;
		SamplesBranch2Lnls[trpIdx] = totalSampleLnlB2;

		std::sort(alSamplesBranch1[trpIdx].begin(),alSamplesBranch1[trpIdx].end(),
				[](const std::pair<double,pair<vector<unsigned char>*, vector<unsigned char>*>* > &left,
						const std::pair<double,pair<vector<unsigned char>*, vector<unsigned char>*>* > &right)
			{
				return left.first >= right.first;
			}
		);
		std::sort(alSamplesBranch2[trpIdx].begin(),alSamplesBranch2[trpIdx].end(),
				[](const std::pair<double,pair<vector<unsigned char>*, vector<unsigned char>*>* > &left,
						const std::pair<double,pair<vector<unsigned char>*, vector<unsigned char>*>* > &right)
				{
					return left.first >= right.first;
				}
		);

		std::sort(alSamplesTriplet[trpIdx].begin(),alSamplesTriplet[trpIdx].end(),
						[](const std::pair<double, array<vector<unsigned char>*,3>* > &left,
								const std::pair<double, array<vector<unsigned char>*,3>* > &right)
						{
							return left.first >= right.first;
						}
				);
	}
}


void ModelEstimator::sampleAlignments()
{
	DEBUG("\n****MODEL ESTIMATOR ALIGNMNENT SAMPLER****\n");
	//FIXME  - check if maps are empty, zero if not
	//Same with likelihood values
	int cap;
	int sampleCount = Definitions::pathSampleCount;
	int analysisCount = Definitions::pathInformativeCount;
	double totalSampleLnlB1;
	double totalSampleLnlB2;
	double lnlB1, lnlB2;
	int ctr;

	//do the first sample

	for (unsigned int trpIdx = 0; trpIdx < samplingHMMs.size(); trpIdx++)
	{
		totalSampleLnlB1 = totalSampleLnlB2 = Definitions::minMatrixLikelihood;
		for(ctr = 1; ctr < sampleCount; ctr++){
			auto branch1Sample = samplingHMMs[trpIdx].first->sampleAlignment(dict,lnlB1);
			auto branch2Sample = samplingHMMs[trpIdx].second->sampleAlignment(dict,lnlB2);

			alSamplesBranch1[trpIdx].push_back(make_pair(lnlB1, branch1Sample));
			alSamplesBranch2[trpIdx].push_back(make_pair(lnlB2, branch2Sample));

			totalSampleLnlB1 = maths->logSum(totalSampleLnlB1, lnlB1);
			totalSampleLnlB2 = maths->logSum(totalSampleLnlB2, lnlB2);
		}

		SamplesBranch1Lnls[trpIdx] = totalSampleLnlB1;
		SamplesBranch2Lnls[trpIdx] = totalSampleLnlB2;

		std::sort(alSamplesBranch1[trpIdx].begin(),alSamplesBranch1[trpIdx].end(),
			[](const std::pair<double,pair<vector<unsigned char>*, vector<unsigned char>*>* > &left,
					const std::pair<double,pair<vector<unsigned char>*, vector<unsigned char>*>* > &right)
			{
				return (left.first >= right.first);
			}
		);
		std::sort(alSamplesBranch2[trpIdx].begin(),alSamplesBranch2[trpIdx].end(),
				[](const std::pair<double,pair<vector<unsigned char>*, vector<unsigned char>*>* > &left,
						const std::pair<double,pair<vector<unsigned char>*, vector<unsigned char>*>* > &right)
				{
					return (left.first >= right.first);
				}
		);


		DUMP("****** Sampler totalSampleLnlB1\t" << totalSampleLnlB1 << " totalSampleLnlB2\t" << totalSampleLnlB2  );

		auto itPr1=alSamplesBranch1[trpIdx].begin();
		auto itPr2=alSamplesBranch2[trpIdx].begin();

		cap = min(alSamplesBranch1[trpIdx].size(),alSamplesBranch2[trpIdx].size());
		cap = Definitions::pathInformativeCount < cap ? Definitions::pathInformativeCount : cap;

		DUMP("****** Sampler " << cap << " triple alignments will be added"  );

		for (unsigned int i=0; i < cap; i++)
		{
			//auto tripletSample = tal->align(itPr1->second, itPr2->second);
			DUMP("****** Sampler adding triple alignment with lnl " << (itPr1->first+itPr2->first)
					<< "\tBranch 1 wt " << exp(itPr1->first - SamplesBranch1Lnls[trpIdx] ) << "\tBranch 2 wt " << exp(itPr2->first - SamplesBranch2Lnls[trpIdx])
					<< "\tTotal wt " << exp(itPr1->first+itPr2->first-SamplesBranch1Lnls[trpIdx]-SamplesBranch2Lnls[trpIdx]));
			alSamplesTriplet[trpIdx].push_back(make_pair((itPr1->first+itPr2->first), tal->align(itPr1->second, itPr2->second)));
			itPr1++;
			itPr2++;
		}

		/*
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
		*/
	}
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
	indelModel->setParameters({initLambda,initEpsilon});

	vector<pair<Band*, Band*> > bandPairs(tripletIdxs.size());
	vector<array<vector<SequenceElement*>*,3> > seqsA(tripletIdxs.size());
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


		f1 = new ForwardPairHMM(seqsA[i][0],seqsA[i][1], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first,true);
		f2 = new ForwardPairHMM(seqsA[i][1],seqsA[i][2], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].second,true);

		f1->setDivergenceTimeAndCalculateModels(distancesA[i][0]*initTimeModifier);
		f2->setDivergenceTimeAndCalculateModels(distancesA[i][1]*initTimeModifier);

		f1->runAlgorithm();
		f2->runAlgorithm();

		this->samplingHMMs.push_back(make_pair(f1,f2));

	}
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
