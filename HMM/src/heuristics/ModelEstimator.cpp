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

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

	this->estimateTripleAlignment(model);

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

void ModelEstimator::estimateTripleAlignment(Definitions::ModelType model)
{
	DEBUG("EstimateTripleAligment");
	//amino acid mode
	bool aaMode = false;

	EvolutionaryPairHMM *f1, *f2;

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

	vector<double> alphas = {0.75};
	vector<double> kappas = {2.5};
	vector<double> lambdas ={0.02, 0.05};
	vector<double> epsilons = {0.3, 0.6};
	vector<double> timeMult = {1.0, 2.0};

	double tmpd;

	int aBest, kBest, lBest, eBest, tBest;

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());

	double lnlBest = std::numeric_limits<double>::max();

	indelModel = new NegativeBinomialGapModel();

	vector<pair<Band*, Band*> > bandPairs(tripletIdxs.size());
	vector<array<vector<SequenceElement>,3> > seqsA(tripletIdxs.size());
	vector<array<double,3> > distancesA(tripletIdxs.size());

	for (int i = 0; i < tripletIdxs.size(); i++)
	{
		//tripletIdxs[idx][0] //first  index
		//tripletIdxs[idx][1] //second  index
		//tripletIdxs[idx][2] //third  index

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


		f1 = new ForwardPairHMM(seqsA[i][0],seqsA[i][1], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first);
		f2 = new ForwardPairHMM(seqsA[i][1],seqsA[i][2], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].second);

	}

	for(int a=0; a < alphas.size(); a++)
	{
		substModel->setAlpha(alphas[a]);
		for(int k=0; k < kappas.size(); k++)
		{
			substModel->setParameters({kappas[k]});
			substModel->calculateModel();
			for(int e=0; e < epsilons.size(); e++)
				for(int l=0; l < lambdas.size(); l++)
				{
					indelModel->setParameters({lambdas[l], epsilons[e]});
					for(int t=0; t < timeMult.size(); t++)
					{
						double lnl = 0;
						for(int h=0; h < tripletIdxs.size(); h++)
						{
							DUMP("Triplet " << h << " forward calculation for alpha " << alphas[a] << " k " << kappas[k] << " l " << lambdas[l] << " e " << epsilons[e] << " t "<< timeMult[t]);
							//posteriorHmms[h].first->setDivergenceTime(distancesA[h][0]*timeMult[t]);
							//posteriorHmms[h].second->setDivergenceTime(distancesA[h][1]*timeMult[t]);
							//lnl += posteriorHmms[h].first->runAlgorithm() + posteriorHmms[h].second->runAlgorithm();

							/*
							BackwardPairHMM* bw1 = new BackwardPairHMM(seqsA[h][0],seqsA[h][1], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[h].first);
							bw1->setDivergenceTime(distancesA[h][0]*timeMult[t]);
							bw1->runAlgorithm();
							DUMP(" POSTERIORS for the first pair alpha " << alphas[a]  << " lambda " << lambdas[l] << " epsilon " << epsilons[e] << " timeM "<< timeMult[t]);
							bw1->calculatePosteriors(dynamic_cast<ForwardPairHMM*>(hmmsA[h].first));
							delete bw1;
							*/
						}
						if (lnl < lnlBest)
						{
							lnlBest = lnl;
							aBest = a;
							kBest = k;
							eBest = e;
							lBest = l;
							tBest = t;
						}
					}
				}
		}
	}

	DUMP("Best values a " << alphas[aBest] << " k " << kappas[kBest] << " l " << lambdas[lBest] << " e " << epsilons[eBest] << " t "<< timeMult[tBest]);
	//found the best combination
	//Run fwd+bwd to get posteriors!

	substModel->setAlpha(alphas[aBest]);
	substModel->setParameters({kappas[kBest]});
	substModel->calculateModel();
	indelModel->setParameters({lambdas[lBest], epsilons[eBest]});

	delete indelModel;
	delete substModel;
	
	

}
void ModelEstimator::sampleAlignments(ForwardPairHMM* hmm)
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
