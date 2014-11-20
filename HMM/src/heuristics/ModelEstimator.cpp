/*
 * ModelEstimator.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#include "heuristics/ModelEstimator.hpp"
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
/*
void ModelEstimator::estimateTripleAlignment(Definitions::ModelType model)
{
	DEBUG("EstimateTripleAligment");
	//amino acid mode
	bool aaMode = false;

	Band* band1;
	Band* band2;

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



	vector<double> alphas = {0.5, 1.0, 3.0};
	vector<double> kappas = {1.5, 3};
	vector<double> lambdas ={0.02, 0.05};
	vector<double> epsilons = {0.3, 0.6};
	vector<double> times = {0.2, 0.6, 1.0, 1.4};

	DEBUG("Setting Frequencies");
	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());


	//substModel->setParameters({2.5});
	//substModel->setAlpha(0.75);
	//substModel->calculateModel();
	//DEBUG("Calculated models");



	pair<string, string> p1;
	pair<string, string> p2;
	double lnlp1, lnlp2;
	double tb1, tb2, tmp;


	indelModel = new NegativeBinomialGapModel();
	//FIXME - hardcodes
	indelModel->setParameters({0.05, 0.5});

	for (int idx = 0; idx < tripletIdxs.size(); idx++)
	{
		//tripletIdxs[idx][0] //first  index
		//tripletIdxs[idx][1] //second  index
		//tripletIdxs[idx][2] //third  index
		auto seq1 = inputSequences->getSequencesAt(tripletIdxs[idx][0]);
		auto seq2 = inputSequences->getSequencesAt(tripletIdxs[idx][1]);
		auto seq3 = inputSequences->getSequencesAt(tripletIdxs[idx][2]);

		unsigned int len1 = seq1.size();
		unsigned int len2 = seq2.size();
		unsigned int len3 = seq3.size();

		band1 = new Band(len1,len2);
		band2 = new Band(len2,len3);

		f1 = new ForwardPairHMM(seq1,seq2, substModel, indelModel, band1);
		f2 = new ForwardPairHMM(seq2,seq3, substModel, indelModel, band2);

		for(int a=0; a < alphas.size(); a++)
		{
			substModel->setAlpha(alphas[a]);
			for(int k=0; aaMode ? 1 : k < kappas.size(); k++)
			{
				substModel->setParameters({kappas[k]});
				substModel->calculateModel();
				for(int e=0; e < epsilons.size(); e++)
					for(int l=0; l < lambdas.size(); l++)
					{
						indelModel->setParameters({lambdas[l], epsilons[e]});
						for(int t1=0; t1 < times.size(); t1++)
						{
							f1->setDivergenceTime(times[t1]);
						}
						for(int t2=0; t2 < times.size(); t2++)
						{
							f2->setDivergenceTime(times[t2]);
						}
					}
			}
		}
		delete band2;
		delete band1
	}




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
