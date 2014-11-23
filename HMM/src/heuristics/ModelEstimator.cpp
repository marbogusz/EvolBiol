/*
 * ModelEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

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
	DEBUG("Sampling triplets osf sequences for gamma shape parameter estimation");
	
	maths = new Maths();
	dict = inputSequences->getDictionary();

	tal = new TripletAligner (inputSequences, gtree->getDistanceMatrix());

	tripletIdxs = tst.sampleFromTree();

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

	this->estimateTripleAlignment(model);

	throw HmmException("Test quits");

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
	}
	ste->optimize();

	DEBUG("Re-estimating model parameters");
	//make another pass
	tripleAlignments.clear();


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
	}
	ste->optimize();

	//end time measurments here!
	end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;

    INFO("Model Estimator elapsed time: " << elapsed_seconds.count() << " seconds");
    cerr <<  "|||||||||||| Model Estimator elapsed time: " << elapsed_seconds.count() << "s ||||||||||||||\n";


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

	//vector<double> alphas = {0.5, 1.0, 3.0};
	//vector<double> kappas = {1.5, 3};
	//vector<double> lambdas ={0.02, 0.05};
	//vector<double> epsilons = {0.3, 0.6};

	vector<double> alphas = {3.0};
	vector<double> kappas = {3};
	vector<double> lambdas ={0.02};
	vector<double> epsilons = {0.6};
	//time multipliers
	vector<double> timeMult = {0.4, 1.0};

	int aBest, kBest, lBest, eBest, tBest;

	DEBUG("Setting Frequencies");
	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());


	double lnlBest = std::numeric_limits<double>::max();


	indelModel = new NegativeBinomialGapModel();
	//FIXME - hardcodes
	indelModel->setParameters({0.05, 0.5});

	vector<pair<Band*, Band*> > bandPairs(tripletIdxs.size());
	vector<array<vector<SequenceElement>,3> > seqsA(tripletIdxs.size());
	vector<array<double,3> > distancesA(tripletIdxs.size());
	vector<pair<EvolutionaryPairHMM*, EvolutionaryPairHMM*>> hmmsA(tripletIdxs.size());

	for (int i = 0; i < tripletIdxs.size(); i++)
	{
		//tripletIdxs[idx][0] //first  index
		//tripletIdxs[idx][1] //second  index
		//tripletIdxs[idx][2] //third  index

		seqsA[i][0] = inputSequences->getSequencesAt(tripletIdxs[i][0]),
		seqsA[i][1] = inputSequences->getSequencesAt(tripletIdxs[i][1]);
		seqsA[i][2] = inputSequences->getSequencesAt(tripletIdxs[i][2]);

		unsigned int len1 = seqsA[i][0].size();
		unsigned int len2 = seqsA[i][1].size();
		unsigned int len3 = seqsA[i][2].size();

		//0-1
		distancesA[i][0] = gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][0],tripletIdxs[i][1]);
		//1-2
		distancesA[i][1] = gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][1],tripletIdxs[i][2]);
		//0-2
		distancesA[i][2] = gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][0],tripletIdxs[i][2]);


//bandPairs[i] = make_pair(new Band(len1,len2),new Band(len2,len3));
		bandPairs[i] = make_pair(nullptr,nullptr);

		f1 = new ForwardPairHMM(seqsA[i][0],seqsA[i][1], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first);
		f2 = new ForwardPairHMM(seqsA[i][1],seqsA[i][2], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].second);

		hmmsA[i] = make_pair(f1,f2);
	}

	for(int a=0; a < alphas.size(); a++)
	{
		substModel->setAlpha(alphas[a]);
		for(int k=0; k < 1; k++)
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
						for(int h=0; h < hmmsA.size(); h++)
						{
							hmmsA[h].first->setDivergenceTime(distancesA[h][0]*timeMult[t]);
							hmmsA[h].second->setDivergenceTime(distancesA[h][1]*timeMult[t]);
							lnl += hmmsA[h].first->runAlgorithm() + hmmsA[h].second->runAlgorithm();
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

	for(int i=0; i < hmmsA.size(); i++)
	{
		hmmsA[i].first->setDivergenceTime(distancesA[i][0]*timeMult[tBest]);
		hmmsA[i].second->setDivergenceTime(distancesA[i][1]*timeMult[tBest]);
		//forward probs
		hmmsA[i].first->runAlgorithm();
		hmmsA[i].second->runAlgorithm();
		//now backward!

		DUMP("Model Estimator First bwd calc");
		BackwardPairHMM* bw1 = new BackwardPairHMM(seqsA[i][0],seqsA[i][1], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first);
		bw1->setDivergenceTime(distancesA[i][0]*timeMult[tBest]);
		bw1->runAlgorithm();
		DUMP("Model Estimator Second bwd calc");
		BackwardPairHMM* bw2 = new BackwardPairHMM(seqsA[i][1],seqsA[i][2], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first);
		bw2->setDivergenceTime(distancesA[i][1]*timeMult[tBest]);
		bw2->runAlgorithm();

		DUMP("Model Estimator First Pair Posteriors");
		bw1->calculatePosteriors(dynamic_cast<ForwardPairHMM*>(hmmsA[i].first));
		DUMP("Model Estimator Second Pair Posteriors");
		bw2->calculatePosteriors(dynamic_cast<ForwardPairHMM*>(hmmsA[i].second));

		delete hmmsA[i].first;
		delete hmmsA[i].second;

		hmmsA[i].first = bw1;
		hmmsA[i].second = bw2;

		ERROR("Ready to sample");

		auto al = bw1->sampleAlignment(inputSequences->getRawSequenceAt(tripletIdxs[i][0]), inputSequences->getRawSequenceAt(tripletIdxs[i][1]));

		//DUMP("alignment ");
		cerr << al.first;
		cerr << al.second;

	}
	delete indelModel;
	delete substModel;
}



/*
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

	pair<string, string> p1;
	pair<string, string> p2;
	double lnlp1, lnlp2;
	double tb1, tb2, tmp;


	indelModel = new NegativeBinomialGapModel();
	//FIXME - hardcodes
	indelModel->setParameters({0.05, 0.5});


		for (int idx = 0; idx < tripletIdxs.size(); idx++)
		{
			lnlp1 = std::numeric_limits<double>::max();
			lnlp2 = std::numeric_limits<double>::max();
			for (auto time : {0.25,0.5,0.75,1.0})
			{
				DEBUG("First Viterbi Pair " << idx << " time " << time);

				vphmm = new ViterbiPairHMM(inputSequences->getSequencesAt(tripletIdxs[idx][0]), inputSequences->getSequencesAt(tripletIdxs[idx][1]),substModel, indelModel);
				tb1 = time; //sme->getModelParams()->getDivergenceTime(idx*3);
				tb2 = time; //sme->getModelParams()->getDivergenceTime((idx*3)+1);
				vphmm->setDivergenceTime(tb1);
				//vphmm->summarize();
				vphmm->runAlgorithm();
				tmp = vphmm->getViterbiSubstitutionLikelihood();
				DEBUG("lnlp1 " << tmp << "\t\t" << lnlp1);
				if(tmp < lnlp1)
				{
					lnlp1=tmp;
					DEBUG("Setting pair 1 with time " << time);
					p1 =  vphmm->getAlignment(inputSequences->getRawSequenceAt(tripletIdxs[idx][0]), inputSequences->getRawSequenceAt(tripletIdxs[idx][1]));
				}
				delete vphmm;
				DEBUG("Second Viterbi Pair " << idx << " time " << time);
	//			DEBUG("Second Viterbi Pair " << idx);
				vphmm = new ViterbiPairHMM(inputSequences->getSequencesAt(tripletIdxs[idx][1]), inputSequences->getSequencesAt(tripletIdxs[idx][2]),substModel, indelModel);
				vphmm->setDivergenceTime(tb2);
				vphmm->runAlgorithm();
				tmp = vphmm->getViterbiSubstitutionLikelihood();
				DEBUG("lnlp2 " << tmp << "\t\t"<< lnlp2);
				if(tmp < lnlp2)
				{
					lnlp2=tmp;
					DEBUG("Setting pair 2 with time " << time);
					p2 =  vphmm->getAlignment(inputSequences->getRawSequenceAt(tripletIdxs[idx][1]), inputSequences->getRawSequenceAt(tripletIdxs[idx][2]));
				}
				delete vphmm;
				//tripleAlignments.push_back(tal.align());
			}
			tripleAlignments.push_back(tal->align(p1,p2));
		}
		delete indelModel;
		delete substModel;
}

*/

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
