//============================================================================
// Name        : HMM.cpp
// Author      : Marcin Bogusz
// Version     :
// Copyright   :
//============================================================================

#include "core/CommandReader.hpp"
#include "core/Sequences.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "core/HmmException.hpp"
#include "core/PairwiseEstimator.hpp"
#include "core/BandingEstimator.hpp"
#include "core/MlEstimator.hpp"
#include "core/BioNJ.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "heuristics/GuideTree.hpp"
#include "heuristics/GotohAlgorithm.hpp"
#include "heuristics/TripletAligner.hpp"
#include "heuristics/ModelEstimator.hpp"
#include <array>

using namespace std;
using namespace EBC;

int main(int argc, char ** argv) {

	//Set output Precision to 2
	//FIXME - should normally be set to >= 6
	cout << fixed << setprecision(2);
	cerr << fixed << setprecision(2);

	try
	{
		CommandReader* cmdReader = new CommandReader(argc, argv);
		DEBUG("Get parser" << endl);

		IParser* parser = cmdReader->getParser();
		DEBUG("Creating alignment");

		Sequences* inputSeqs = new Sequences(parser, cmdReader->getSequenceType(),cmdReader->isFixedAlignment());
		if (cmdReader->isMLE())
		{
			//MlEstimator* me = new MlEstimator(inputSeqs, cmdReader->getModelType() ,cmdReader->getIndelParams(),
			//					cmdReader->getSubstParams(), cmdReader->getOptimizationType(), cmdReader->getCategories(),
			//					cmdReader->getAlpha(), cmdReader->estimateAlpha(),cmdReader->getDistance(), cmdReader->isFixedAlignment() == false);
			//DEBUG("Creating TripletModelEstimator");

			ModelEstimator* tme = new ModelEstimator(inputSeqs, cmdReader->getModelType(),
					cmdReader->getOptimizationType(), cmdReader->getCategories(), cmdReader->getAlpha(),
					cmdReader->estimateAlpha());


			//tme->getModelParameters();

		}

		else
		{
			/*
			ForwardPairHMM* fwdHMM = new ForwardPairHMM(inputSeqs, cmdReader->getModelType() ,
					cmdReader->getIndelParams(),cmdReader->getSubstParams(), cmdReader->getOptimizationType(),
					cmdReader->getBanding(), cmdReader->getBandFactor(), cmdReader->getDistance(),
					cmdReader->getCategories(), cmdReader->getAlpha(), cmdReader->estimateAlpha());



			ForwardPairHMM* fwdHMM = new ForwardPairHMM(inputSeqs->getSequencesAt(0), inputSeqs->getSequencesAt(1),
					inputSeqs->getDictionary(), cmdReader->getModelType() , cmdReader->getBanding(), cmdReader->getBandFactor(),
					cmdReader->getCategories(), cmdReader->getAlpha(), new Maths());


			fwdHMM->setModelFrequencies(inputSeqs->getElementFrequencies());
			fwdHMM->setModelParameters(cmdReader->getIndelParams(),cmdReader->getSubstParams(), cmdReader->getDistance(),0);
			fwdHMM->runForwardAlgorithm();
			*/

			ModelEstimator* tme = new ModelEstimator(inputSeqs, cmdReader->getModelType(),
					cmdReader->getOptimizationType(), cmdReader->getCategories(), cmdReader->getAlpha(),
					cmdReader->estimateAlpha());


			vector<double> indelParams;
			vector<double> substParams;
			double alpha;

			substParams = tme->getSubstitutionParameters();
			indelParams = tme->getIndelParameters();
			if(cmdReader->estimateAlpha())
				alpha = tme->getAlpha();


			//FIXME - hardcoding substitution parameters and alpha to come from the estimator
			BandingEstimator* be = new BandingEstimator(cmdReader->getAlgorithmType(), inputSeqs, cmdReader->getModelType() ,indelParams,
					substParams, cmdReader->getOptimizationType(), cmdReader->getCategories(),alpha, tme->getGuideTree());

			//DEBUG ("Running bionj");

			//change bionj init here!
			//BioNJ nj(inputSeqs->getSequenceCount(), be->getOptimizedTimes());
			//DEBUG("Final tree : " << nj.calculate());
			//cout << nj.calculate() << endl;

		}

		//ForwardPairHMM* epHMM = new ForwardPairHMM(inputSeqs);

		//ViterbiPairHMM* epHMM = new ViterbiPairHMM(inputSeqs);
		//epHMM->runViterbiAlgorithm();
		//epHMM->runForwardAlgorithm();
		//epHMM->getResults();
		//delete epHMM;
		//ForwardPairHMM* fwdHMM = new ForwardPairHMM(inputSeqs,true);
		//fwdHMM->summarize();
	}
	catch(HmmException& pe)
	{
		cerr << pe.what();
	}
	catch(exception &ex)
	{
		cerr << ex.what();
	}
	return 0;
}
