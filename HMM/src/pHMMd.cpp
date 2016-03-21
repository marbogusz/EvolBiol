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
#include <thread>

#include "core/FileLogger.hpp"

#include "core/OptimizedModelParameters.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/ViterbiPairHMM.hpp"

#include "models/CodonModel.hpp"


using namespace std;
using namespace EBC;

int main(int argc, char ** argv) {

	cerr << "Starting dm\n";
	//Set output Precision to 2
	//FIXME - should normally be set to >= 6
	cout << fixed << setprecision(6);
	cerr << fixed << setprecision(6);


	std::this_thread::sleep_for(std::chrono::milliseconds(30000));
	cerr << " ELLO\n";
	cout << " Should not be here!\n";
	return 0;

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

		//indelModel->setParameters({Definitions::almostZero,0.5});

		substModel->setObservedFrequencies(inputSeqs->getElementFrequencies());

		modelParams = new OptimizedModelParameters(substModel, indelModel,2, 1, true,
							true, false, true, maths);

		//modelParams->generateInitialDistanceParameters();
		//modelParams->generateInitialIndelParameters();
		//modelParams->generateInitialSubstitutionParameters();
		modelParams->setUserIndelParams({0.01,0.5});
		modelParams->setUserDivergenceParams({0.2});
		modelParams->setUserSubstParams({2.0, 0.1});

		//indelModel->setParameters({0.00001,0.00001});
		//substModel->setParameters({16.4, 0.00352});
		//substModel->calculateModel();


		EvolutionaryPairHMM *hmm;

		bfgs = new Optimizer(modelParams, NULL, cmdReader->getOptimizationType());

		PairHmmCalculationWrapper* wrapper = new PairHmmCalculationWrapper();



		for (unsigned int pi = 0; pi<inputSeqs->getPairCount(); pi++){

			std::pair<unsigned int, unsigned int> idxs = inputSeqs->getPairOfSequenceIndices(pi);

			unsigned int len1, len2;

			len1 = inputSeqs->getSequencesAt(idxs.first)->size();
			len2 = inputSeqs->getSequencesAt(idxs.second)->size();

			Band* band = new Band(len1,len2,0.3);


/*
			auto ip =  cmdReader->getIndelParams();

			substModel->setParameters(cmdReader->getSubstParams());
			substModel->calculateModel();
			//BAND it ?


			hmm = new ForwardPairHMM(inputSeqs->getSequencesAt(idxs.first), inputSeqs->getSequencesAt(idxs.second),
							substModel, indelModel, Definitions::DpMatrixType::Full, band);


			double lam = Definitions::almostZero;
			double tim = 0.1;//Definitions::almostZero;
			double lnl = 0;

			cout << "Time\tLambda\tLnL" << endl;
			while(tim < 1.5){
				lam = Definitions::almostZero;
				while(lam < 0.1){
					indelModel->setParameters({lam,ip[1]});
					hmm->setDivergenceTimeAndCalculateModels(tim);
					lnl = hmm->runAlgorithm() * -1.0;
					lam += (0.0025);
					//lam += (Definitions::almostZero)*4;
					cout << tim << "\t" << lam << "\t" << lnl << endl;
				}
				tim += 0.025;
			}
*/


			hmm = new ForwardPairHMM(inputSeqs->getSequencesAt(idxs.first), inputSeqs->getSequencesAt(idxs.second),
					substModel, indelModel, Definitions::DpMatrixType::Full, band);

			wrapper->setTargetHMM(hmm);
			wrapper->setIndelModel(indelModel);
			wrapper->setSubstModel(substModel);
			wrapper->setModelParameters(modelParams);

			bfgs->setTarget(wrapper);
			bfgs->optimize();

			//detect if we're close to the bounds

			double lambda, divergence;
			bool runAgain;

			do{

				runAgain = false;
				lambda = modelParams->getIndelParameters()[0];
				divergence  = modelParams->getDivergenceTime(0);

				//cerr << " L " << lambda << " D " << divergence << endl;


				if(lambda > (Definitions::lambdaHiBound * 0.975)){
					runAgain = true;
					//cerr << " Lambda big \n";
					Definitions::lambdaHiBound = Definitions::lambdaHiBound * 2.0;
					Definitions::divergenceBound = Definitions::divergenceBound / 2.0;

				}	//check if we're close to the band
				else if(divergence > (Definitions::divergenceBound * 0.975)){
					runAgain = true;
					//cerr << " Divergence big \n";
					Definitions::lambdaHiBound = Definitions::lambdaHiBound / 2.0;
					Definitions::divergenceBound = Definitions::divergenceBound * 2.0;
				}

				if(runAgain){
					//cerr << "Run again...\n";
					//cerr << "New bounds " << Definitions::lambdaHiBound << " " << Definitions::divergenceBound << endl;
					indelModel->resetBounds();
					modelParams->resetBounds();
					bfgs->optimize();
				}

			}
			while(runAgain == true);


/*
			auto vh = new ViterbiPairHMM(inputSeqs->getSequencesAt(idxs.first), inputSeqs->getSequencesAt(idxs.second),
					substModel, indelModel, Definitions::DpMatrixType::Full, NULL);

			vh->setDivergenceTimeAndCalculateModels(modelParams->getDivergenceTime(0));
			vh->runAlgorithm();
			auto al = vh->getStringAlignment();

			cout << al.first << endl;
			cout << al.second << endl;

			delete vh;
			*/
			cout << modelParams->getSubstParameters()[0] << "\t" << modelParams->getSubstParameters()[1] <<
			    "\t" <<  modelParams->getIndelParameters()[0] << "\t" <<  modelParams->getIndelParameters()[1] <<
				"\t" << modelParams->getDivergenceTime(0) << "\t" << inputSeqs->getSequenceName(idxs.first) << "\t" << inputSeqs->getSequenceName(idxs.second) << "\t" << pi << "\n";

			//INFO(inputSeqs->getSequenceName(idxs.first) << " " << inputSeqs->getSequenceName(idxs.second));

			//INFO("Divergence time " << modelParams->getDivergenceTime(0));
			//INFO(modelParams->getSubstParameters());
			//INFO(modelParams->getIndelParameters());



			delete hmm;
			delete band;
		}

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
