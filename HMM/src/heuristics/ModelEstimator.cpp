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

	posteriorHmms.resize(tripletIdxs.size());
	smpldPairs.resize(tripletIdxs.size());

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

	this->estimateTripleAlignment(model);

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


		(posteriorHmms[i]).first = f1;
		(posteriorHmms[i]).second = f2;
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
						for(int h=0; h < posteriorHmms.size(); h++)
						{
							DUMP("Triplet " << h << " forward calculation for alpha " << alphas[a] << " k " << kappas[k] << " l " << lambdas[l] << " e " << epsilons[e] << " t "<< timeMult[t]);
							posteriorHmms[h].first->setDivergenceTime(distancesA[h][0]*timeMult[t]);
							posteriorHmms[h].second->setDivergenceTime(distancesA[h][1]*timeMult[t]);
							lnl += posteriorHmms[h].first->runAlgorithm() + posteriorHmms[h].second->runAlgorithm();

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

	for(int i=0; i < posteriorHmms.size(); i++)
	{
		DUMP("Triplet " << i << " re-calculating forward likelihoods for the best parameters");
		posteriorHmms[i].first->setDivergenceTime(distancesA[i][0]*timeMult[tBest]);
		posteriorHmms[i].second->setDivergenceTime(distancesA[i][1]*timeMult[tBest]);
		//forward probs
		posteriorHmms[i].first->runAlgorithm();
		posteriorHmms[i].second->runAlgorithm();

		ViterbiPairHMM* v1 = new ViterbiPairHMM(seqsA[i][0],seqsA[i][1], substModel, indelModel, Definitions::DpMatrixType::Full);

		ViterbiPairHMM* v2 = new ViterbiPairHMM(seqsA[i][1],seqsA[i][2], substModel, indelModel, Definitions::DpMatrixType::Full);

		v1->setDivergenceTime(distancesA[i][0]*timeMult[tBest]);
		v2->setDivergenceTime(distancesA[i][1]*timeMult[tBest]);

		v1->runAlgorithm();
		v2->runAlgorithm();

		auto a1 = v1->getAlignment(inputSequences->getRawSequenceAt(tripletIdxs[i][0]), inputSequences->getRawSequenceAt(tripletIdxs[i][1]));
		auto a2 = v2->getAlignment(inputSequences->getRawSequenceAt(tripletIdxs[i][1]), inputSequences->getRawSequenceAt(tripletIdxs[i][2]));

		DUMP("Extra triplet viterbi : ");
		auto al = tal->align(a1,a2);

		/*
		//now backward!

		DUMP("Triplet " << i <<  " model Estimator First backward calculation");
		BackwardPairHMM* bw1 = new BackwardPairHMM(seqsA[i][0],seqsA[i][1], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first);
		bw1->setDivergenceTime(distancesA[i][0]*timeMult[tBest]);
		bw1->runAlgorithm();
		DUMP("Triplet " << i <<  " model Estimator Second backward calculation");
		BackwardPairHMM* bw2 = new BackwardPairHMM(seqsA[i][1],seqsA[i][2], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first);
		bw2->setDivergenceTime(distancesA[i][1]*timeMult[tBest]);
		bw2->runAlgorithm();

		DUMP("Triplet " << i <<  " model Estimator First Pair Posteriors");
		bw1->calculatePosteriors(dynamic_cast<ForwardPairHMM*>(posteriorHmms[i].first));
		DUMP("Triplet " << i <<  " model Estimator Second Pair Posteriors");
		bw2->calculatePosteriors(dynamic_cast<ForwardPairHMM*>(posteriorHmms[i].second));
		 */

		//delete posteriorHmms[i].first;
		//delete posteriorHmms[i].second;

		//posteriorHmms[i].first = bw1;
		//posteriorHmms[i].second = bw2;

		//ERROR("Ready to sample");

		//auto al = bw1->sampleAlignment(inputSequences->getRawSequenceAt(tripletIdxs[i][0]), inputSequences->getRawSequenceAt(tripletIdxs[i][1]));

		//DUMP("alignment ");
		//cerr << al->first;
		//cerr << al->second;

	}
	delete indelModel;
	delete substModel;

}
array<vector<SequenceElement>, 3> ModelEstimator::sampleTripleAlignment(unsigned int triplet)
{
	DUMP("Sample triple alignment for triplet " << triplet);
	Dictionary* di = inputSequences->getDictionary();

	pair<string,string> p1 = dynamic_cast<ForwardPairHMM*>(posteriorHmms[triplet].first)->sampleAlignment(inputSequences->getRawSequenceAt(tripletIdxs[triplet][0]), inputSequences->getRawSequenceAt(tripletIdxs[triplet][1]));
	pair<string,string> p2 = dynamic_cast<ForwardPairHMM*>(posteriorHmms[triplet].second)->sampleAlignment(inputSequences->getRawSequenceAt(tripletIdxs[triplet][1]), inputSequences->getRawSequenceAt(tripletIdxs[triplet][2]));


	DUMP("p1.1 " << p1.first);
	DUMP("p1.2 " << p1.second);
	DUMP("p2.1 " << p2.first);
	DUMP("p2.2 " << p2.second);

	smpldPairs[triplet].push_back({di->translate(p1.first), di->translate(p1.second),
		di->translate(p2.first), di->translate(p2.second)});

	//TODO - do only 1 translation - pass translated seqs to tal
	return tal->align(p1,p2);
}


ModelEstimator::~ModelEstimator()
{
    delete maths;
    delete sme;
    delete ste;
    delete gtree;
    delete tal;

    for (auto p : posteriorHmms ){
    	delete p.first;
    	delete p.second;
    }
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
