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
		ofstream distfile;

		FileLogger::start(cmdReader->getLoggingLevel(), (string(cmdReader->getInputFileName()).append(Definitions::logExt)));



		IParser* parser = cmdReader->getParser();

		//FileLogger::DebugLogger().setCerr();
		//FileLogger::DumpLogger().setCerr();
		//FileLogger::InfoLogger().setCerr();

		INFO("Reading input sequences...");
		DEBUG("Creating alignment object...");

		Sequences* inputSeqs = new Sequences(parser, cmdReader->getSequenceType(),cmdReader->isFixedAlignment());
		if (cmdReader->isFdist() || cmdReader->isVdist())
		{
			vector<double> indelParams;
			vector<double> substParams;
			substParams = cmdReader->getSubstParams();
			indelParams = cmdReader->getIndelParams();
			double alpha = cmdReader->getAlpha();

			Optimizer* bfgs;
			Dictionary* dict;
			SubstitutionModelBase* substModel;
			IndelModel* indelModel;
			Maths* maths;
			Definitions::AlgorithmType algorithm;
			OptimizedModelParameters* modelParams;

			maths = new Maths();
			dict = inputSeqs->getDictionary();

			if (cmdReader->getModelType() == Definitions::ModelType::GTR)
			{
				substModel = new GTRModel(dict, maths,cmdReader->getCategories());
			}
			else if (cmdReader->getModelType() == Definitions::ModelType::HKY85)
			{
				substModel = new HKY85Model(dict, maths,cmdReader->getCategories());
			}
			else if (cmdReader->getModelType() == Definitions::ModelType::LG)
			{
					substModel = new AminoacidSubstitutionModel(dict, maths,cmdReader->getCategories(),Definitions::aaLgModel);
			}

			indelModel = new NegativeBinomialGapModel();

			modelParams = new OptimizedModelParameters(substModel, indelModel,2, 1, false,
					false, false, true, maths);

			modelParams->generateInitialDistanceParameters();


			//modelParams->setUserIndelParams(indelParams);
			//modelParams->setUserSubstParams(substParams);

			substModel->setObservedFrequencies(inputSeqs->getElementFrequencies());
			substModel->setAlpha(alpha);
			substModel->setParameters(substParams);
			substModel->calculateModel();


			indelModel->setParameters(indelParams);

			EvolutionaryPairHMM *hmm;

			bfgs = new Optimizer(modelParams, NULL, cmdReader->getOptimizationType());

			PairHmmCalculationWrapper* wrapper = new PairHmmCalculationWrapper();

			std::pair<unsigned int, unsigned int> idxs = inputSeqs->getPairOfSequenceIndices(0);

			if (cmdReader->isVdist())
			{
				DEBUG("Creating Viterbi algorithm to optimize the pairwise divergence time...");
				hmm = new ViterbiPairHMM(inputSeqs->getSequencesAt(idxs.first), inputSeqs->getSequencesAt(idxs.second),
					substModel, indelModel, Definitions::DpMatrixType::Full, nullptr);
			}
			else if (cmdReader->isFdist())
			{
				DEBUG("Creating forward algorithm to optimize the pairwise divergence time...");
				hmm = new ForwardPairHMM(inputSeqs->getSequencesAt(idxs.first), inputSeqs->getSequencesAt(idxs.second),
					substModel, indelModel, Definitions::DpMatrixType::Full, nullptr);
			}

			//hmm->setDivergenceTime(modelParams->getDivergenceTime(0)); //zero as there's only one pair!
			wrapper->setTargetHMM(hmm);
			wrapper->setModelParameters(modelParams);
			bfgs->setTarget(wrapper);
			bfgs->optimize();
			cout << modelParams->getDivergenceTime(0);
		}

		else if (cmdReader->isMLE())
		{

			//DISTANCE-BASED estimation using an alignment

			vector<double> indelParams;
			vector<double> substParams;
			double alpha;

			substParams = cmdReader->getSubstParams();
			indelParams = cmdReader->getIndelParams();
			alpha = cmdReader->getAlpha();

			ModelEstimator* tme;
/*
			if (substParams.size() == 0){
				tme = new ModelEstimator(inputSeqs, cmdReader->getModelType(),
								cmdReader->getOptimizationType(), cmdReader->getCategories(), cmdReader->getAlpha(),
								false);

				substParams = tme->getSubstitutionParameters();

			}
*/

			GuideTree gt(inputSeqs);


			MlEstimator mle(inputSeqs, cmdReader->getModelType(), indelParams, substParams,
							cmdReader->getOptimizationType(), cmdReader->getCategories(), alpha,
							false, gt.getDistances(), false);


			//cout << tme.getOptimizedTimes()[0];

			BioNJ nj(inputSeqs->getSequenceCount(), mle.getOptimizedTimes(), inputSeqs);
			DEBUG("Final tree : " << nj.calculate());
			string treeStr = nj.calculate();

			treefile.open((string(cmdReader->getInputFileName()).append(".nj.tree")).c_str(),ios::out);
			treefile << treeStr << endl;
			treefile.close();

		}

		else
		{



			INFO("Creating Model Parameters heuristics...");
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



			//FileLogger::Logger() << "True indel paramteres     : ";
			//FileLogger::Logger() << cmdReader->getIndelParams() << '\n';
			//FileLogger::Logger() << "Estimated indel paramteres: ";
			//FileLogger::Logger() << indelParams << '\n';
			//FileLogger::Logger() << "True substitution paramteres     : ";
			//FileLogger::Logger() << cmdReader->getSubstParams();
			//FileLogger::Logger() << "Estimated substitution paramteres: ";
			//FileLogger::Logger() << substParams;
			//FileLogger::Logger() << "True alpha      : " << cmdReader->getAlpha() << "\n";
			//FileLogger::Logger() << "Estimated alpha : " << alpha << "\n";

			//FIXME - hardcoding substitution parameters and alpha to come from the estimator
			BandingEstimator* be = new BandingEstimator(cmdReader->getAlgorithmType(), inputSeqs, cmdReader->getModelType() ,indelParams,
					substParams, cmdReader->getOptimizationType(), cmdReader->getCategories(),alpha, tme->getGuideTree());
			be->optimizePairByPair();


			auto distances = be->getOptimizedTimes();
			auto seqCount =  inputSeqs->getSequenceCount();

			//output distance matrix
			distfile.open((string(cmdReader->getInputFileName()).append(Definitions::distMatExt)).c_str(),ios::out);
			distfile << inputSeqs->getSequenceCount() << endl;
			for (unsigned int seqId = 0; seqId < seqCount; seqId++){
				distfile << inputSeqs->getSequenceName(seqId) << "        ";
				for(unsigned int j = 0; j<seqId; j++)
				{

					distfile << " " << distances[(seqId - j - 1) + (j*seqCount) - (((1+j)/2.0)*(j*1.0))];
				}
				distfile << endl;
			}
			distfile.close();


			DEBUG ("Running BioNJ");
			//change bionj init here!
			BioNJ nj(inputSeqs->getSequenceCount(), be->getOptimizedTimes(), inputSeqs);
			//DEBUG("Final tree : " << nj.calculate());
			string treeStr = nj.calculate();


			INFO("Indel parameters");
			INFO(indelParams);
			INFO("Substitution parameters");
			INFO(substParams);
			INFO("Gamma parameters (alpha and rate categories)");
			INFO(alpha << '\t' << cmdReader->getCategories());
			INFO("Newick tree");
			INFO(treeStr);


			treefile.open((string(cmdReader->getInputFileName()).append(Definitions::treeExt)).c_str(),ios::out);
			treefile << treeStr << endl;
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
