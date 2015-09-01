//============================================================================
// Name        : HMM.cpp
// Author      : Marcin Bogusz
// Version     :
// Copyright   :
//============================================================================


#include <heuristics/PairSamplingTree.hpp>
#include <heuristics/RawIndelEstimator.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>
#include <vector>
#include <chrono>
#include <ctime>

#include "core/ProgramException.hpp"
#include "core/CommandReader.hpp"
#include "core/Sequences.hpp"
#include "core/FileLogger.hpp"
#include "core/OptimizedModelParameters.hpp"
#include "core/PhylogeneticTree.hpp"

#include "heuristics/StateTransitionEstimator.hpp"

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
			throw ProgramException("Can't open the tree file");
		}
		PhylogeneticTree* ptree = new PhylogeneticTree(inputSeqs);
		ptree->fromNewick(newick);


		IndelModel* indelModel = new NegativeBinomialGapModel();


		if(rawMode){

			SubstitutionModelBase* substModel;

			INFO("Creating Model Parameters heuristics...");

			PairSamplingTree tst(ptree->getDistanceMatrix());
			vector<array<unsigned int, 3> > tripletIdxs = tst.sampleFromTree();
			unsigned int tripletIdxsSize = tripletIdxs.size();

			ModelEstimator* tme = new ModelEstimator(inputSeqs, cmdReader->getModelType(), ptree->getDistanceMatrix(),
					cmdReader->getCategories(), cmdReader->getAlpha(),
					cmdReader->estimateAlpha());

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

			StateTransitionEstimator* ste = new StateTransitionEstimator(indelModel, inputSeqs->getDictionary()->getGapID());


			for(unsigned int i =0; i< inputSeqs->getPairCount(); i++)
			{
				std::pair<unsigned int, unsigned int> idxs = inputSeqs->getPairOfSequenceIndices(i);
				ste->addPair(inputSeqs->getSequencesAt(idxs.first), inputSeqs->getSequencesAt(idxs.second), ptree->distanceById(idxs.first,idxs.second));
			}
			ste->optimize();
			auto params = ste->getModelParams()->getIndelParameters();

			//need distances for those pairs!
			cout << "Indel rate " << params[0] << endl;
			cout << "Indel length " << params[1] << endl;
			delete ste;


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
