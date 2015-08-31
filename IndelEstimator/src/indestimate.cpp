//============================================================================
// Name        : HMM.cpp
// Author      : Marcin Bogusz
// Version     :
// Copyright   :
//============================================================================

#include <core/ProgramException.hpp>
#include "core/CommandReader.hpp"
#include "core/Sequences.hpp"
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

using namespace std;
using namespace EBC;

int main(int argc, char ** argv) {

	cout << fixed << setprecision(6);
	cerr << fixed << setprecision(6);

	bool rawMode = false;

	try
	{
		//Get some time statistics
	    chrono::time_point<chrono::system_clock> start, end;
	    start = chrono::system_clock::now();

	    CommandReader* cmdReader = new CommandReader(argc, argv);
	    FileLogger::start(cmdReader->getLoggingLevel(), (string(cmdReader->getInputFileName()).append(".indel.log")));

		IParser* parser = cmdReader->getParser();
		rawMode = cmdReader->isRawInput();


		INFO("Reading input sequences...");

		//while working in raw mode indels will be removed if a MSA was provided!
		Sequences* inputSeqs = new Sequences(parser, cmdReader->getSequenceType(),!rawMode);
		INFO("Reading tree file...");
		ifstream treeFile (cmdReader->getTreeFileName());
		string newick;
		if (treeFile.is_open())
		{
			//FIXME - parse newick file in case the string is not located in the first line
			getline (treeFile,newick);
			treeFile.close();
		}
		else{
			throw ProgramExcepton("Can't open the tree file");
		}
		PhylogeneticTree* ptree = new PhylogeneticTree(this->inputSeqs);
		ptree->fromNewick(newick);


		IndelModel* indelModel = new NegativeBinomialGapModel();


		if(rawMode){

			SubstitutionModelBase* substModel;

			INFO("Creating Model Parameters heuristics...");
			/*
			ModelEstimator* tme = new ModelEstimator(inputSeqs, cmdReader->getModelType(),
					cmdReader->getOptimizationType(), cmdReader->getCategories(), cmdReader->getAlpha(),
					cmdReader->estimateAlpha());

			vector<double> indelParams;
			vector<double> substParams;
			double alpha = 100;

			substParams = tme->getSubstitutionParameters();
			indelParams = tme->getIndelParameters();
			if(cmdReader->estimateAlpha())
				alpha = tme->getAlpha();

			//FIXME - hardcoding substitution parameters and alpha to come from the estimator
			BandingEstimator* be = new BandingEstimator(cmdReader->getAlgorithmType(), inputSeqs, cmdReader->getModelType() ,indelParams,
					substParams, cmdReader->getOptimizationType(), cmdReader->getCategories(),alpha, tme->getGuideTree());
			be->optimizePairByPair();

			INFO("Indel parameters");
			INFO(indelParams);


			delete be;
			delete tme;
			*/
		}
		else{
			//Alignment mode

			//need a set of pairs
			//getPairs (with ids)

			int i = 1;
			std::pair<unsigned int, unsigned int> idxs = inputSeqs->getPairOfSequenceIndices(i);


			//starting values;
			double initLambda = 0.025;
			double initEpsilon = 0.5;

			indelModel->setParameters({initLambda,initEpsilon});

			StateTransitionEstimator ste = new StateTransitionEstimator(indelModel, ot, 2*tripletIdxsSize, dict->getGapID(),false);


			for(unsigned int i =0; i< inputSeqs->getPairCount(); i++)
			{
				std::pair<unsigned int, unsigned int> idxs = inputSeqs->getPairOfSequenceIndices(i);

				ste->addTime(tb1+tb2,trp,0);
				ste->addTime(tb3+tb2,trp,1);

				ste->addPair(pairAlignments[trp][0],pairAlignments[trp][1],trp,0);
				ste->addPair(pairAlignments[trp][2],pairAlignments[trp][3],trp,1);

			}
			ste->optimize();

			//need distances for those pairs!


		}


		delete ptree;
	    delete inputSeqs;
		delete parser;
		delete cmdReader;



		end = chrono::system_clock::now();
	    chrono::duration<double> elapsed_seconds = end-start;
	    std::time_t end_time = chrono::system_clock::to_time_t(end);

	    INFO("Finished computation at " << std::ctime(&end_time) << " elapsed time: " << elapsed_seconds.count() << "s\n");

	}
	catch(ProgramException& pe)
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
