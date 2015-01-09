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
#include <map>

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

	ste = new StateTransitionEstimator(ot);

	for(int al = 0; al < tripleAlignments.size(); al++)
	{

		tb1 = sme->getModelParams()->getDivergenceTime(al*3);
		tb2 = sme->getModelParams()->getDivergenceTime((al*3)+1);
		tb3 = sme->getModelParams()->getDivergenceTime((al*3)+2);

		DEBUG(tb1 << " " << tb2 << " " << tb3);
		ste->addPair(tripleAlignments[al][0],tripleAlignments[al][1],tb1+tb2);
		ste->addPair(tripleAlignments[al][1],tripleAlignments[al][2],tb2+tb3);
		//ste->addPair(pairAlignments[al][0],pairAlignments[al][1],tb1+tb2);
		//ste->addPair(pairAlignments[al][2],pairAlignments[al][3],tb2+tb3);
	}
	ste->optimize();

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

	delete sme;
	delete ste;



	sme = new SubstitutionModelEstimator(inputSeqs, model ,ot, rateCategories, alpha, estimateAlpha, tripletIdxs.size());

	for(int al = 0; al < tripleAlignments.size(); al++)
	{
		sme->addTriplet(tripleAlignments[al]);
	}

	sme->optimize();

	ste = new StateTransitionEstimator(ot);

	for(int al = 0; al < tripleAlignments.size(); al++)
	{
		tb1 = sme->getModelParams()->getDivergenceTime(al*3);
		tb2 = sme->getModelParams()->getDivergenceTime((al*3)+1);
		tb3 = sme->getModelParams()->getDivergenceTime((al*3)+2);
		ste->addPair(tripleAlignments[al][0],tripleAlignments[al][1],tb1+tb2);
		ste->addPair(tripleAlignments[al][1],tripleAlignments[al][2],tb2+tb3);
		//ste->addPair(pairAlignments[al][0],pairAlignments[al][1],tb1+tb2);
		//ste->addPair(pairAlignments[al][2],pairAlignments[al][3],tb2+tb3);
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
	substModel->setParameters({2.0});
	substModel->setAlpha(99);
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

	tripletIdxs[0][0] = 0;
	tripletIdxs[0][1] = 1;
	tripletIdxs[0][2] = 2;

		for (int idx = 0; idx < tripletIdxs.size(); idx++)
		{
			lnlp1 = std::numeric_limits<double>::max();
			//lnlp2 = std::numeric_limits<double>::max();

			vphmm1 = new ViterbiPairHMM(inputSequences->getSequencesAt(tripletIdxs[idx][0]), inputSequences->getSequencesAt(tripletIdxs[idx][1]),substModel, indelModel);
			vphmm2 = new ViterbiPairHMM(inputSequences->getSequencesAt(tripletIdxs[idx][1]), inputSequences->getSequencesAt(tripletIdxs[idx][2]),substModel, indelModel);

			tb1 = gtree->getDistanceMatrix()->getDistance(tripletIdxs[idx][0],tripletIdxs[idx][1]);
		    tb2 = gtree->getDistanceMatrix()->getDistance(tripletIdxs[idx][1],tripletIdxs[idx][2]);
			for (auto lambda : {0.03})
				for (auto epsilon : {0.8}){
					indelModel->setParameters({lambda, epsilon});
					for (auto time : {0.4})
					{
						//DEBUG("First Viterbi Pair " << idx << "time multiplier" << time);
						vphmm1->setDivergenceTime(time);
						vphmm1->runAlgorithm();

						vphmm2->setDivergenceTime(time);
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

			tb1 = 0.4;//gtree->getDistanceMatrix()->getDistance(tripletIdxs[0][0],tripletIdxs[0][1]);
		    tb2 = 0.4;//gtree->getDistanceMatrix()->getDistance(tripletIdxs[0][1],tripletIdxs[0][2]);
			indelModel->setParameters({l, e});

			//DEBUG("###########################MODEL PARAMETERS##################");
			//DEBUG("HKY85 kappa " << 2.5 << " discrete gamma " << this->gammaRateCategories <<  "categories, alpha = 0.75");
			//DEBUG("Indel parameters : lambda: " << l << " geometric gap opening (epsilon) " << e );
			//DEBUG("Divergence time : " << (tb1*t));

			substModel->summarize();


			DEBUG("*******VITERBI********");
			vphmm1 = new ViterbiPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][0]), inputSequences->getSequencesAt(tripletIdxs[0][1]),substModel, indelModel);
			vphmm1->setDivergenceTime(tb1);
			vphmm1->runAlgorithm();


			//vphmm2 = new ViterbiPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][1]), inputSequences->getSequencesAt(tripletIdxs[0][2]),substModel, indelModel);
			//vphmm2->setDivergenceTime(tb2*t);
			//vphmm2->runAlgorithm();
			DEBUG("Viterbi alignment NO EQUILIBRIUMS");
			vp1 =  vphmm1->getBestAlignment(inputSequences->getRawSequenceAt(tripletIdxs[0][0]), inputSequences->getRawSequenceAt(tripletIdxs[0][1]));
			DEBUG(vp1.first);
			DEBUG(vp1.second);
			//DUMP("Second Pair");
			//vp2 =  vphmm2->getBestAlignment(inputSequences->getRawSequenceAt(tripletIdxs[0][1]), inputSequences->getRawSequenceAt(tripletIdxs[0][2]));
			//DUMP(vp2.first);
			//DUMP(vp2.second);
			//dynamic_cast<DpMatrixFull*>(vphmm1->M->getDpMatrix())->outputValues(0);
			//dynamic_cast<DpMatrixFull*>(vphmm1->X->getDpMatrix())->outputValues(0);
			//dynamic_cast<DpMatrixFull*>(vphmm1->Y->getDpMatrix())->outputValues(0);
			DUMP("*******VITERBI lnL NO Equilibruims********");
			DUMP(vphmm1->getAlignmentLikelihood(dict->translate(vp1.first), dict->translate(vp1.second)));
			/*
			DUMP("*******VITERBI Lnl2********");
			DUMP(vphmm2->getAlignmentLikelihood(dict->translate(vp2.first), dict->translate(vp2.second)));
			*/

			vphmm1 = new ViterbiPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][0]), inputSequences->getSequencesAt(tripletIdxs[0][1]),substModel, indelModel,Definitions::DpMatrixType::Full, nullptr, true);
			vphmm1->setDivergenceTime(tb1);
			vphmm1->runAlgorithm();
			//vphmm2 = new ViterbiPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][1]), inputSequences->getSequencesAt(tripletIdxs[0][2]),substModel, indelModel);
			//vphmm2->setDivergenceTime(tb2*t);
			//vphmm2->runAlgorithm();
			DEBUG("Viterbi alignment WITH EQUILIBRIUMS");
			vp1 =  vphmm1->getBestAlignment(inputSequences->getRawSequenceAt(tripletIdxs[0][0]), inputSequences->getRawSequenceAt(tripletIdxs[0][1]));
			DEBUG(vp1.first);
			DEBUG(vp1.second);
			//DUMP("Second Pair");
			//vp2 =  vphmm2->getBestAlignment(inputSequences->getRawSequenceAt(tripletIdxs[0][1]), inputSequences->getRawSequenceAt(tripletIdxs[0][2]));
			//DUMP(vp2.first);
			//DUMP(vp2.second);
			//dynamic_cast<DpMatrixFull*>(vphmm1->M->getDpMatrix())->outputValues(0);
			//dynamic_cast<DpMatrixFull*>(vphmm1->X->getDpMatrix())->outputValues(0);
			//dynamic_cast<DpMatrixFull*>(vphmm1->Y->getDpMatrix())->outputValues(0);
			DUMP("*******VITERBI lnL WITH EQUILIBRIUMS********");
			DUMP(vphmm1->getAlignmentLikelihood(dict->translate(vp1.first), dict->translate(vp1.second)));

			//DUMP("*******VITERBI Lnl2********");
			//DUMP(vphmm2->getAlignmentLikelihood(dict->translate(vp2.first), dict->translate(vp2.second)));

			ForwardPairHMM* fphmm1;
			ForwardPairHMM* fphmm2;

			DEBUG("*******    Forward  ********");
			fphmm1 = new ForwardPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][0]), inputSequences->getSequencesAt(tripletIdxs[0][1]),substModel, indelModel,Definitions::DpMatrixType::Full, nullptr, true);
			fphmm1->setDivergenceTime(tb1);
			fphmm1->runAlgorithm();
			//fphmm2 = new ForwardPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][1]), inputSequences->getSequencesAt(tripletIdxs[0][2]),substModel, indelModel,Definitions::DpMatrixType::Full, nullptr);
			//fphmm2->setDivergenceTime(tb2*t);
			//fphmm2->runAlgorithm();

			DUMP("Forward best alignment");
			fp1 = fphmm1->getBestAlignment(inputSequences->getRawSequenceAt(tripletIdxs[0][0]), inputSequences->getRawSequenceAt(tripletIdxs[0][1]));

			DEBUG(fp1.first);
			DEBUG(fp1.second);

			//DUMP("Second Pair");
			//fp2 = fphmm2->getBestAlignment(inputSequences->getRawSequenceAt(tripletIdxs[0][1]), inputSequences->getRawSequenceAt(tripletIdxs[0][2]));

			//DUMP(fp2.first);
			//DUMP(fp2.second);
/*
			DUMP("*******FWD Lnl1********");
			DUMP(fphmm1->getAlignmentLikelihood(dict->translate(fp1.first), dict->translate(fp1.second)));
			DUMP("*******FWD Lnl2********");
			DUMP(fphmm2->getAlignmentLikelihood(dict->translate(fp2.first), dict->translate(fp2.second)));
*/
			BackwardPairHMM* bphmm1;
			//BackwardPairHMM* bphmm2;

			DEBUG("*******    Backward  ********");
			bphmm1 = new BackwardPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][0]), inputSequences->getSequencesAt(tripletIdxs[0][1]),substModel, indelModel,Definitions::DpMatrixType::Full, nullptr);
			bphmm1->setDivergenceTime(tb1);
			bphmm1->runAlgorithm();

			//DUMP("#####Match BW########");
			//dynamic_cast<DpMatrixFull*>(bphmm1->M->getDpMatrix())->outputValues(0);
			//DUMP("#####Insert BW########");
			//dynamic_cast<DpMatrixFull*>(bphmm1->X->getDpMatrix())->outputValues(0);
			//DUMP("#####Delete BW########");
			//dynamic_cast<DpMatrixFull*>(bphmm1->Y->getDpMatrix())->outputValues(0);

			//bphmm2 = new BackwardPairHMM(inputSequences->getSequencesAt(tripletIdxs[0][1]), inputSequences->getSequencesAt(tripletIdxs[0][2]),substModel, indelModel,Definitions::DpMatrixType::Full, nullptr);
			//bphmm2->setDivergenceTime(tb2*t);
			//bphmm2->runAlgorithm();

			DEBUG("*******   Calculate posteriors   ********");
			bphmm1->calculatePosteriors(fphmm1);
			//bphmm2->calculatePosteriors(fphmm2);

			double vlnl, flnl;

			DEBUG("*******VITERBI Lnl1********");
			vlnl = fphmm1->getAlignmentLikelihood(dict->translate(vp1.first), dict->translate(vp1.second));
			vlnl -=0.0001;
			DEBUG(vlnl);
			DEBUG("*******FWD Lnl1********");
			flnl = fphmm1->getAlignmentLikelihood(dict->translate(fp1.first), dict->translate(fp1.second));
			DEBUG(flnl);



			//DUMP("**************************************************");
			//DUMP("***********Forward Samples 1st pair***************");





			int s1len = (inputSequences->getRawSequenceAt(tripletIdxs[0][0])).size()+1;
			int s2len = (inputSequences->getRawSequenceAt(tripletIdxs[0][1])).size()+1;

			vector<vector<double> > posteriors (s1len , std::vector<double> ( s2len, 0.00000001 ));

			double tlnl;

			vector<double> freqs(50,0.0);

			map<double, pair<string, string> > alSamples;

			int ctr;
			int total = 0;
			for(ctr = 0; ctr < 1000; ctr++){
				auto pr1 = fphmm1->sampleAlignment(inputSequences->getRawSequenceAt(tripletIdxs[0][0]), inputSequences->getRawSequenceAt(tripletIdxs[0][1]));
				tlnl = bphmm1->getAlignmentLikelihood(dict->translate(pr1.first), dict->translate(pr1.second),false, posteriors);

				alSamples.insert(make_pair(tlnl, pr1));

				for(int cats=0; cats<50;cats++)
				{
					if (tlnl > vlnl-((cats+1)*0.5)){
						freqs[cats] +=1;
						total ++;
						break;
					}
				}

			}



			stringstream sstr;

			sstr << endl;
			double lv;


			for(unsigned int i=0; i < s1len; i++)
			{
				for(unsigned int j=0; j < s2len; j++)
				{
					lv = log(posteriors[i][j]/ctr);
					if (lv > -5.0)
						sstr << (int)round(lv*-1.0);
					else sstr << ".";

				}
				sstr << endl;
			}
			DUMP("&&&&&&&&COUNTED POSTERIORS FOR MATCH&&&&&&&&&&&&&&&&&&");
			DUMP(sstr.str());


			for(int ct =0; ct < freqs.size(); ct++)
				INFO("lnl up to " << ((ct+1)*-0.5) << "\t\t " << (freqs[ct]/ctr));

			DEBUG("Top 100 alignments");
			auto it=alSamples.rbegin();
			int cap = 100 < alSamples.size() ? 100 : alSamples.size();

			for (int c=0;c<cap; c++)
			{
				DEBUG(it->first);
				DEBUG(it->second.first);
				DEBUG(it->second.second);
				it++;
			}


		}
		delete indelModel;
		delete substModel;
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
