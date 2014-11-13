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

#include "core/FileLogger.hpp"

using namespace std;
using namespace EBC;

int main(int argc, char ** argv) {

	//Set output Precision to 2
	//FIXME - should normally be set to >= 6
	cout << fixed << setprecision(1);
	cerr << fixed << setprecision(1);

	try
	{

		//FIXME - nothing happens when the model does not get specified!

		CommandReader* cmdReader = new CommandReader(argc, argv);
		ofstream treefile;
		FileLogger::start(cmdReader->getLoggingLevel(), (string(cmdReader->getInputFileName()).append(".hmm.log")));

		treefile.open((string(cmdReader->getInputFileName()).append(".hmm.tree")).c_str(),ios::out);

		IParser* parser = cmdReader->getParser();

		FileLogger::InfoLogger() << "Reading sequences" << "\n";
		FileLogger::DebugLogger() << "Creating alignment" << "\n";

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
			delete tme;

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

			FileLogger::InfoLogger() << "Creating Model Parameters heuristics" << "\n";
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
			be->optimizePairByPair();

			//DEBUG ("Running bionj");

			//change bionj init here!
			BioNJ nj(inputSeqs->getSequenceCount(), be->getOptimizedTimes());
			//DEBUG("Final tree : " << nj.calculate());
			treefile << nj.calculate() << endl;


			treefile.close();
			delete be;
			delete tme;

		}


		delete inputSeqs;
		delete parser;
		delete cmdReader;

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

	FileLogger::stop();
	return 0;
}
