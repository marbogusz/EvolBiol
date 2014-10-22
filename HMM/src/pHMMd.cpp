//============================================================================
// Name        : HMM.cpp
// Author      : Marcin Bogusz
// Version     :
// Copyright   :
//============================================================================

#include "core/CommandReader.hpp"
#include "core/Sequences.hpp"
#include "hmm/BasicViterbi.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "core/HmmException.hpp"
#include "core/PairwiseEstimator.hpp"
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

	//Set output Precision to 6
	cout << fixed << setprecision(8);
	cerr << fixed << setprecision(8);

	try
	{
		CommandReader* cmdReader = new CommandReader(argc, argv);
		DEBUG("Get parser" << endl);

		IParser* parser = cmdReader->getParser();
		DEBUG("Creating alignment");

		Sequences* inputSeqs = new Sequences(parser, cmdReader->getSequenceType(),cmdReader->isFixedAlignment());
		Sequences* inputSeqsTrue = new Sequences(cmdReader->getParserTrueAlignment(), cmdReader->getSequenceType(),true);

		if (cmdReader->isMLE())
		{

			ModelEstimator* tme = new ModelEstimator(inputSeqs, inputSeqsTrue, cmdReader->getModelType(),
					cmdReader->getOptimizationType(), cmdReader->getCategories(), cmdReader->getAlpha(),
					cmdReader->estimateAlpha());

		}

		else
		{

			ModelEstimator* tme = new ModelEstimator(inputSeqs, inputSeqsTrue, cmdReader->getModelType(),
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
			PairwiseEstimator* pe = new PairwiseEstimator(cmdReader->getAlgorithmType(), inputSeqs, cmdReader->getModelType() ,indelParams/*cmdReader->getIndelParams()*/,
					/*cmdReader->getSubstParams()*/ substParams, cmdReader->getOptimizationType(), cmdReader->getBanding(), cmdReader->getBandFactor(),
					cmdReader->getCategories(), /*cmdReader->getAlpha()*/ alpha, /*cmdReader->estimateAlpha()*/ false,cmdReader->getDistance());



			DEBUG ("Running bionj");

			//change bionj init here!
			BioNJ nj(inputSeqs->getSequenceCount(), pe->getOptimizedTimes());
			//DEBUG("Final tree : " << nj.calculate());
			cout << nj.calculate() << endl;



			delete tme;
			delete pe;
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
