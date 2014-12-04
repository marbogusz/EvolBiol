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
#include <chrono>
#include <ctime>

#include "core/FileLogger.hpp"

using namespace std;
using namespace EBC;

int main(int argc, char ** argv) {

	//Set output Precision to 2
	//FIXME - should normally be set to >= 6
	//cout << fixed << setprecision(4);
	//cerr << fixed << setprecision(4);

	try
	{
		//Get some time statistics
	    chrono::time_point<chrono::system_clock> start, end;
	    start = chrono::system_clock::now();

		//FIXME - nothing happens when the model does not get specified!

		CommandReader* cmdReader = new CommandReader(argc, argv);
		ofstream treefile;

		FileLogger::start(cmdReader->getLoggingLevel(), (string(cmdReader->getInputFileName()).append(".hmm.log")));



		IParser* parser = cmdReader->getParser();

		//FileLogger::DebugLogger().setCerr();
		//FileLogger::DumpLogger().setCerr();
		//FileLogger::InfoLogger().setCerr();

		INFO("Reading input sequences...");
		DEBUG("Creating alignment object...");

		Sequences* inputSeqs = new Sequences(parser, cmdReader->getSequenceType(),cmdReader->isFixedAlignment());
		if (cmdReader->isMLE())
		{
			//tme->getModelParameters();
			INFO("Creating Model Parameters heuristics...");
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

			cerr << indelParams[0] << "\t" << indelParams[1] << "\n";
			delete tme;

		}

		else
		{

			treefile.open((string(cmdReader->getInputFileName()).append(".hmm.tree")).c_str(),ios::out);
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


			INFO("Creating Model Parameters heuristics...");
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


			delete tme;
			//FileLogger::Logger() << "True indel paramteres     : ";
			//FileLogger::Logger() << cmdReader->getIndelParams() << '\n';
			FileLogger::Logger() << "Estimated indel paramteres: ";
			FileLogger::Logger() << indelParams << '\n';
			//FileLogger::Logger() << "True substitution paramteres     : ";
			//FileLogger::Logger() << cmdReader->getSubstParams();
			FileLogger::Logger() << "Estimated substitution paramteres: ";
			FileLogger::Logger() << substParams;
			//FileLogger::Logger() << "True alpha      : " << cmdReader->getAlpha() << "\n";
			FileLogger::Logger() << "Estimated alpha : " << alpha << "\n";


			//FIXME - hardcoding substitution parameters and alpha to come from the estimator
			BandingEstimator* be = new BandingEstimator(cmdReader->getAlgorithmType(), inputSeqs, cmdReader->getModelType() ,indelParams,
					substParams, cmdReader->getOptimizationType(), cmdReader->getCategories(),alpha, tme->getGuideTree());
			be->optimizePairByPair();

			INFO ("Running BioNJ");

			//change bionj init here!
			BioNJ nj(inputSeqs->getSequenceCount(), be->getOptimizedTimes(), inputSeqs);
			//DEBUG("Final tree : " << nj.calculate());
			treefile << nj.calculate() << endl;


			treefile.close();
			delete be;

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


		end = chrono::system_clock::now();
	    chrono::duration<double> elapsed_seconds = end-start;
	    std::time_t end_time = chrono::system_clock::to_time_t(end);

	    INFO("Finished computation at " << std::ctime(&end_time) << " elapsed time: " << elapsed_seconds.count() << "s\n");

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
