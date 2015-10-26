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

#include "models/CodonModel.hpp"


using namespace std;
using namespace EBC;

int main(int argc, char ** argv) {

	//Set output Precision to 2
	//FIXME - should normally be set to >= 6
	cout << fixed << setprecision(8);
	cerr << fixed << setprecision(8);

	try
	{
		//Get some time statistics
	    chrono::time_point<chrono::system_clock> start, end;
	    start = chrono::system_clock::now();

		//FIXME - nothing happens when the model does not get specified!

		CommandReader* cmdReader = new CommandReader(argc, argv);
		ofstream treefile;

		FileLogger::start(cmdReader->getLoggingLevel(), (string(cmdReader->getInputFileName()).append(Definitions::logExt)));

		IParser* parser = cmdReader->getParser();

		//FileLogger::DebugLogger().setCerr();
		//FileLogger::DumpLogger().setCerr();
		FileLogger::InfoLogger().setCerr();

		INFO("Reading input sequences...");
		DEBUG("Creating alignment object...");

		Sequences* inputSeqs = new Sequences(parser, cmdReader->getSequenceType(),cmdReader->isFixedAlignment());

		//No indel or subst params!
		Optimizer* bfgs;
		Dictionary* dict;
		SubstitutionModelBase* substModel;
		IndelModel* indelModel;
		Maths* maths;
		Definitions::AlgorithmType algorithm;
		OptimizedModelParameters* modelParams;

		maths = new Maths();
		dict = inputSeqs->getDictionary();

		//no gamma with codon models
		substModel = new CodonModel(dict, maths,1);
		indelModel = new NegativeBinomialGapModel();

		substModel->setObservedFrequencies(inputSeqs->getElementFrequencies());

		modelParams = new OptimizedModelParameters(substModel, indelModel,2, 1, true,
							true, false, true, maths);

		modelParams->generateInitialDistanceParameters();
		modelParams->generateInitialIndelParameters();
		modelParams->generateInitialSubstitutionParameters();

		EvolutionaryPairHMM *hmm;

		bfgs = new Optimizer(modelParams, NULL, cmdReader->getOptimizationType());

		PairHmmCalculationWrapper* wrapper = new PairHmmCalculationWrapper();

		std::pair<unsigned int, unsigned int> idxs = inputSeqs->getPairOfSequenceIndices(0);

		unsigned int len1, len2;

		len1 = inputSeqs->getSequencesAt(idxs.first)->size();
		len2 = inputSeqs->getSequencesAt(idxs.second)->size();

		//BAND it ?
		hmm = new ForwardPairHMM(inputSeqs->getSequencesAt(idxs.first), inputSeqs->getSequencesAt(idxs.second),
							substModel, indelModel, Definitions::DpMatrixType::Full, new Band(len1,len2,0.25));

		wrapper->setTargetHMM(hmm);
		wrapper->setIndelModel(indelModel);
		wrapper->setSubstModel(substModel);
		wrapper->setModelParameters(modelParams);

		bfgs->setTarget(wrapper);
		bfgs->optimize();

		INFO("Divergence time " << modelParams->getDivergenceTime(0));
		INFO(modelParams->getSubstParameters());
		INFO(modelParams->getIndelParameters());

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
