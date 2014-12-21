/*
 * ModelEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "heuristics/ModelEstimator.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/BackwardPairHMM.hpp"
#include "hmm/DpMatrixFull.hpp"
#include <chrono>

namespace EBC
{

ModelEstimator::ModelEstimator(Sequences* inputSeqs, Definitions::ModelType model ,
		Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, bool estimateAlpha) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories),
				gtree(new GuideTree(inputSeqs)), tst(*gtree)
{
	DEBUG("About to sample some triplets");
	DEBUG("Sampling triplets osf sequences for gamma shape parameter estimation");
	
	maths = new Maths();
	dict = inputSequences->getDictionary();

	tal = new TripletAligner (inputSequences, gtree->getDistanceMatrix());

	tripletIdxs = tst.sampleFromTree();

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

	this->estimateTripleAlignment(model);

	sme = new SubstitutionModelEstimator(inputSeqs, model ,ot, rateCategories, alpha, estimateAlpha, tripletIdxs.size());

	for(int al = 0; al < tripleAlignments.size(); al++)
	{
		sme->addTriplet(tripleAlignments[al]);
	}

	sme->optimize();

	double tb1,tb2,tb3;
	
	tb1 = sme->getModelParams()->getDivergenceTime(0);
	tb2 = sme->getModelParams()->getDivergenceTime(1);
	tb3 = sme->getModelParams()->getDivergenceTime(2);
	
	int cap;
	//2 pairs
	
	ste = new StateTransitionEstimator(ot, 2);
	
	ste->addTime(tb1+tb2,0,0);
	ste->addTime(tb3+tb2,0,1);
	
	
	auto it=alSamples1.rbegin();
	cap = Definitions::pathInformativeCount < alSamples1.size() ? Definitions::pathInformativeCount : alSamples1.size();
	for (unsigned int i=0; i < cap; i++)
	{
		//DUMP("lnl: " << it->first << "\t total lnl: " << totalSampleLnl);
		ste->addPair((it->second).first,(it->second).second,0,0,(double)exp((it->first)-(this->totalSampleLnl1)));
		it++;
	}
	
	auto it2=alSamples2.rbegin();
	cap = Definitions::pathInformativeCount < alSamples2.size() ? Definitions::pathInformativeCount : alSamples2.size();
	for (unsigned int i=0; i < cap; i++)
	{
		//DUMP("lnl: " << it->first << "\t total lnl: " << totalSampleLnl);
		ste->addPair((it2->second).first,(it2->second).second,0,1,(double)exp((it2->first)-(this->totalSampleLnl2)));
		it2++;
	}
	ste->optimize();
	
	alSamples1.clear();
	alSamples2.clear();
	

	DEBUG("Re-estimating model parameters");
	//make another pass
	tripleAlignments.clear();
	pairAlignments.clear();


	indelModel =  ste->getIndelModel();
	substModel =  sme->getSubstModel();
	
	//indelModel->summarize();


	for (int idx = 0; idx < tripletIdxs.size(); idx++)
	{

//		DEBUG("First Viterbi Pair " << idx);
		vphmm = new ViterbiPairHMM(inputSeqs->getSequencesAt(tripletIdxs[idx][0]), inputSeqs->getSequencesAt(tripletIdxs[idx][1]), substModel, indelModel);
		tb1 = sme->getModelParams()->getDivergenceTime(idx*3);
		tb2 = sme->getModelParams()->getDivergenceTime((idx*3)+1);
		tb3 = sme->getModelParams()->getDivergenceTime((idx*3)+2);

		vphmm->setDivergenceTime(tb1+tb2);
		//vphmm->summarize();
		vphmm->runAlgorithm();
		auto p1 =  vphmm->getAlignment(inputSeqs->getRawSequenceAt(tripletIdxs[idx][0]), inputSeqs->getRawSequenceAt(tripletIdxs[idx][1]));
		delete vphmm;
//		DEBUG("Second Viterbi Pair " << idx);
		vphmm = new ViterbiPairHMM(inputSeqs->getSequencesAt(tripletIdxs[idx][1]), inputSeqs->getSequencesAt(tripletIdxs[idx][2]), substModel, indelModel);
		vphmm->setDivergenceTime(tb2+tb3);
		vphmm->runAlgorithm();
		auto p2 =  vphmm->getAlignment(inputSeqs->getRawSequenceAt(tripletIdxs[idx][1]), inputSeqs->getRawSequenceAt(tripletIdxs[idx][2]));
		delete vphmm;
		tripleAlignments.push_back(tal->align(p1,p2));
		pairAlignments.push_back({{dict->translate(p1.first),dict->translate(p1.second),dict->translate(p2.first),dict->translate(p2.second)}});
		//tripleAlignments.push_back(tal.align());
	}

	//delete sme;
	//delete ste;



	sme = new SubstitutionModelEstimator(inputSeqs, model ,ot, rateCategories, alpha, estimateAlpha, tripletIdxs.size());

	for(int al = 0; al < tripleAlignments.size(); al++)
	{
		sme->addTriplet(tripleAlignments[al]);
	}

	sme->optimize();
	
	tb1 = sme->getModelParams()->getDivergenceTime(0);
	tb2 = sme->getModelParams()->getDivergenceTime(1);
	tb3 = sme->getModelParams()->getDivergenceTime(2);

	substModel =  sme->getSubstModel();
	indelModel =  ste->getIndelModel();

	delete fphmm1;
	delete fphmm2;

	DEBUG("*******    Forward  ********");

	fphmm1 = new ForwardPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][0]), inputSequences->getSequencesAt(tripletIdxs[0][1]),substModel, indelModel,Definitions::DpMatrixType::Full, nullptr);
	fphmm1->setDivergenceTime(tb1+tb2);
	fphmm1->runAlgorithm();

	sampleAlignments(fphmm1,alSamples1,totalSampleLnl1);

	fphmm2 = new ForwardPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][1]), inputSequences->getSequencesAt(tripletIdxs[0][2]),substModel, indelModel,Definitions::DpMatrixType::Full, nullptr);
	fphmm2->setDivergenceTime(tb3+tb2);
	fphmm2->runAlgorithm();

	sampleAlignments(fphmm2,alSamples2,totalSampleLnl2);

	
	ste = new StateTransitionEstimator(ot, 2);
	
	ste->addTime(tb1+tb2,0,0);
	ste->addTime(tb3+tb2,0,1);
	


	it=alSamples1.rbegin();
	cap = Definitions::pathInformativeCount < alSamples1.size() ? Definitions::pathInformativeCount : alSamples1.size();
	for (unsigned int i=0; i < cap; i++)
	{
		//DUMP("lnl: " << it->first << "\t total lnl: " << totalSampleLnl);
		ste->addPair((it->second).first,(it->second).second,0,0,(double)exp((it->first)-(this->totalSampleLnl1)));
		it++;
	}
	
	it2=alSamples2.rbegin();
	cap = Definitions::pathInformativeCount < alSamples2.size() ? Definitions::pathInformativeCount : alSamples2.size();
	for (unsigned int i=0; i < cap; i++)
	{
		//DUMP("lnl: " << it->first << "\t total lnl: " << totalSampleLnl);
		ste->addPair((it2->second).first,(it2->second).second,0,1,(double)exp((it2->first)-(this->totalSampleLnl2)));
		it2++;
	}
	ste->optimize();



	//end time measurments here!
	end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;

    cerr <<  "|||||||||||| elapsed time: " << elapsed_seconds.count() << "s ||||||||||||||\n";


	indelModel =  ste->getIndelModel();
	substModel =  sme->getSubstModel();

	//substModel->summarize();
	//indelModel->summarize();
	//substModel->summarize();
	//we have new alignments!
	//re-estimate
	//do Viterbi using the estimates
	//construct triplets
}

void ModelEstimator::sampleAlignments(ForwardPairHMM* hmm, map<double, pair<vector<SequenceElement>, vector<SequenceElement> > >& alSamples, double& totalSampleLnl)
{



	int sampleCount = Definitions::pathSampleCount;
	int analysisCount = Definitions::pathInformativeCount;
	totalSampleLnl = 0;
	double lnl;
	int ctr;
	pair<string, string> smplPr;
	pair<double, pair<vector<SequenceElement>, vector<SequenceElement> > > pr;
	totalSampleLnl = Definitions::minMatrixLikelihood;

	//do the first sample



	for(ctr = 1; ctr < sampleCount; ctr++){
		smplPr = hmm->sampleAlignment(inputSequences->getRawSequenceAt(tripletIdxs[0][0]), inputSequences->getRawSequenceAt(tripletIdxs[0][1]));
		pr = make_pair(0.0, make_pair(dict->translate(smplPr.first), dict->translate(smplPr.second)));



		lnl = hmm->getAlignmentLikelihood(pr.second.first, pr.second.second);
		totalSampleLnl = maths->logSum(totalSampleLnl, lnl);
		pr.first = lnl;
		//cerr << smplPr.first << endl;
		//cerr << smplPr.second << endl;
		//cerr << lnl << endl;
		alSamples.insert(pr);

	}
}

void ModelEstimator::estimateTripleAlignment(Definitions::ModelType model)
{
	DEBUG("EstimateTripleAligment");
	if (model == Definitions::ModelType::GTR || model == Definitions::ModelType::HKY85)
	{
			//FIXME - make model idiotproof by checking if parameters are set;
			DEBUG("Setting HKY85");
			this->substModel = new HKY85Model(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::LG)
	{
			substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
	}

	DEBUG("Setting Frequencies");
	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
	substModel->setParameters({2.5});
	substModel->setAlpha(0.75);
	substModel->calculateModel();
	DEBUG("Calculated models");

	pair<string, string> vp1;
	pair<string, string> vp2;
	pair<string, string> fp1;
	pair<string, string> fp2;
	double lnlp1, lnlp2;
	double tb1, tb2, tmp;

	double l,e,t;


	indelModel = new NegativeBinomialGapModel();
	//FIXME - hardcodes
	//indelModel->setParameters({0.05, 0.5});

	ViterbiPairHMM* vphmm1;
	ViterbiPairHMM* vphmm2;

		for (int idx = 0; idx < tripletIdxs.size(); idx++)
		{
			lnlp1 = std::numeric_limits<double>::max();
			//lnlp2 = std::numeric_limits<double>::max();

			vphmm1 = new ViterbiPairHMM(inputSequences->getSequencesAt(tripletIdxs[idx][0]), inputSequences->getSequencesAt(tripletIdxs[idx][1]),substModel, indelModel);
			vphmm2 = new ViterbiPairHMM(inputSequences->getSequencesAt(tripletIdxs[idx][1]), inputSequences->getSequencesAt(tripletIdxs[idx][2]),substModel, indelModel);

			tb1 = gtree->getDistanceMatrix()->getDistance(tripletIdxs[idx][0],tripletIdxs[idx][1]);
		    tb2 = gtree->getDistanceMatrix()->getDistance(tripletIdxs[idx][1],tripletIdxs[idx][2]);
			for (auto lambda : {0.02, 0.06})
				for (auto epsilon : {0.25, 0.5}){
					indelModel->setParameters({lambda, epsilon});
					for (auto time : {0.4,1.0,1.6})
					{
						//DEBUG("First Viterbi Pair " << idx << "time multiplier" << time);
						vphmm1->setDivergenceTime(tb1*time);
						vphmm1->runAlgorithm();

						vphmm2->setDivergenceTime(tb2*time);
						vphmm2->runAlgorithm();

						tmp = vphmm1->getViterbiSubstitutionLikelihood() + vphmm2->getViterbiSubstitutionLikelihood();
						DEBUG("Fwd likelihood sum for lambda " << lambda << " epsilon " << epsilon << " time mult " << time << " : " << tmp << "\t\t" << lnlp1);
						if(tmp < lnlp1)
						{
							lnlp1=tmp;
							l = lambda;
							t = time;
							e = epsilon;
							vp1 =  vphmm1->getAlignment(inputSequences->getRawSequenceAt(tripletIdxs[idx][0]), inputSequences->getRawSequenceAt(tripletIdxs[idx][1]));
							vp2 =  vphmm2->getAlignment(inputSequences->getRawSequenceAt(tripletIdxs[idx][1]), inputSequences->getRawSequenceAt(tripletIdxs[idx][2]));
						}
					}
				}
			tripleAlignments.push_back(tal->align(vp1,vp2));
			pairAlignments.push_back({{dict->translate(vp1.first),dict->translate(vp1.second),dict->translate(vp2.first),dict->translate(vp2.second)}});
			delete vphmm1;
			delete vphmm2;
		}
		indelModel->setParameters({l, e});

		tb1 = gtree->getDistanceMatrix()->getDistance(tripletIdxs[0][0],tripletIdxs[0][1]);
	    tb2 = gtree->getDistanceMatrix()->getDistance(tripletIdxs[0][1],tripletIdxs[0][2]);

		DEBUG("*******    Forward  ********");

		fphmm1 = new ForwardPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][0]), inputSequences->getSequencesAt(tripletIdxs[0][1]),substModel, indelModel,Definitions::DpMatrixType::Full, nullptr);
		fphmm1->setDivergenceTime(t*tb1);
		fphmm1->runAlgorithm();

		sampleAlignments(fphmm1,alSamples1,totalSampleLnl1);

		fphmm2 = new ForwardPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][1]), inputSequences->getSequencesAt(tripletIdxs[0][2]),substModel, indelModel,Definitions::DpMatrixType::Full, nullptr);
		fphmm2->setDivergenceTime(t*tb2);
		fphmm2->runAlgorithm();

		sampleAlignments(fphmm2,alSamples2,totalSampleLnl2);
		
		//delete indelModel;
		//delete substModel;
}



ModelEstimator::~ModelEstimator()
{
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
