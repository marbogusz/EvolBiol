//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

#include "core/CommandReader.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/BandingEstimator.hpp"
#include "core/MlEstimator.hpp"
#include "core/BioNJ.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "heuristics/ModelEstimator.hpp"
#include <array>
#include <chrono>
#include <ctime>

#include "core/FileLogger.hpp"

#include "core/OptimizedModelParameters.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/ViterbiPairHMM.hpp"

#include "sampling/HMMEstimator.hpp";

using namespace std;
using namespace EBC;

int main(int argc, char ** argv) {


	cout << fixed << setprecision(8);
	cerr << fixed << setprecision(8);

	//cout << Definitions::notice;

	try
	{
		//Get some time statistics
	    chrono::time_point<chrono::system_clock> start, end;
	    start = chrono::system_clock::now();

		CommandReader* cmdReader = new CommandReader(argc, argv);
		ofstream treefile;
		ofstream distfile;

		FileLogger::start(cmdReader->getLoggingLevel(), (string(cmdReader->getInputFileName()).append(Definitions::logExt)));



		IParser* parser = cmdReader->getParser();

		//Remove gaps if the user provides a MSA file
		bool removeGaps = true;

		//FileLogger::DebugLogger().setCerr();
		//FileLogger::DumpLogger().setCerr();
		//FileLogger::InfoLogger().setCerr();

		INFO("Reading input sequences...");
		DEBUG("Creating alignment object...");

		Sequences* inputSeqs = new Sequences(parser, cmdReader->getSequenceType(),removeGaps);

		INFO("Creating Model Parameters heuristics...");

		//cout << "Estimating evolutionary model parameters..." << endl;

		//ModelEstimator* tme = new ModelEstimator(inputSeqs, cmdReader->getModelType(),
		//		cmdReader->getOptimizationType(), cmdReader->getCategories(), cmdReader->getAlpha(),
		//		cmdReader->estimateAlpha());

		vector<double> indelParams;
		vector<double> substParams;
		double alpha = cmdReader->getAlpha();

		substParams = {};//tme->getSubstitutionParameters();
		indelParams = {0.01,0.5};//tme->getIndelParameters();

		alpha = 2;

		if(cmdReader->estimateAlpha()){
			alpha = 2; //tme->getAlpha();
		}

/*
		try{
			substParams = cmdReader->getSubstParams();
		}
		//do nothing - if exception, this means no user-specified params
		catch(HmmException& pe){
			substParams = tme->getSubstitutionParameters();
		}

		try{
			indelParams = cmdReader->getIndelParams();
		}
		//do nothing - if exception, this means no user-specified params
		catch(HmmException& pe){
			indelParams = tme->getIndelParameters();
		}
*/
		//cout << "Estimating pairwise distances..." << endl;



		//get epsilon, lambda and times from zorro code


		BandingEstimator* be = new BandingEstimator(Definitions::AlgorithmType::Forward, inputSeqs, cmdReader->getModelType() ,indelParams,
				substParams, cmdReader->getOptimizationType(), cmdReader->getCategories(),alpha, NULL);
		be->optimizePairByPair();


		//auto distances = be->getOptimizedTimes();
		//auto seqCount =  inputSeqs->getSequenceCount();

		delete be;

		//delete tme;


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
		ERROR(pe.what());
	}
	catch(exception &ex)
	{
		ERROR(ex.what());
	}

	FileLogger::stop();
	return 0;
}
