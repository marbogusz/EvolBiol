//============================================================================
// Name        : HMM.cpp
// Author      : Marcin Bogusz
// Version     :
// Copyright   :
//============================================================================

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
	cout << fixed << setprecision(6);
	cerr << fixed << setprecision(6);

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
		if (cmdReader->isFdist() || cmdReader->isVdist())
		{
			vector<double> indelParams;
			vector<double> substParams;
			substParams = cmdReader->getSubstParams();
			indelParams = cmdReader->getIndelParams();

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
				substModel = new GTRModel(dict, maths,1);
			}
			else if (cmdReader->getModelType() == Definitions::ModelType::HKY85)
			{
				substModel = new HKY85Model(dict, maths,1);
			}
			else if (cmdReader->getModelType() == Definitions::ModelType::LG)
			{
					substModel = new AminoacidSubstitutionModel(dict, maths,1,Definitions::aaLgModel);
			}

			indelModel = new NegativeBinomialGapModel();

			modelParams = new OptimizedModelParameters(substModel, indelModel,2, 1, false,
					false, false, true, maths);


			modelParams->setUserIndelParams(indelParams);
			modelParams->setUserSubstParams(substParams);

			substModel->setObservedFrequencies(inputSeqs->getElementFrequencies());
			substModel->setParameters(modelParams->getSubstParameters());
			substModel->calculateModel();


			indelModel->setParameters(modelParams->getIndelParameters());

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

			vector<double> indelParams;
			vector<double> substParams;
			double alpha = alpha = cmdReader->getAlpha();

			double dist =  cmdReader->getDistance();
			if (dist < 0)
				dist = 1.0;

			substParams = cmdReader->getSubstParams();
			indelParams = cmdReader->getIndelParams();

			//tme->getModelParameters();

			//cerr << "Alpha " << alpha << endl;
			//cerr << "Rate cat " << cmdReader->getCategories() << endl;
			INFO("Creating Model Parameters heuristics...");

			ModelEstimator* tme = new ModelEstimator(inputSeqs, cmdReader->getModelType(),
					cmdReader->getOptimizationType(), cmdReader->getCategories(), cmdReader->getAlpha(),
					cmdReader->estimateAlpha());

			//tme->getModelParameters();
			substParams = tme->getSubstitutionParameters();
			indelParams = tme->getIndelParameters();
			if(cmdReader->estimateAlpha())
				alpha = tme->getAlpha();

			cout << alpha << '\t' << indelParams[0] << '\t' << indelParams[1];
			for (auto param : substParams)
				cout  << '\t' << param;

			delete tme;


/*
			HMMEstimator* tme = new HMMEstimator(inputSeqs, cmdReader->getModelType(),
					cmdReader->getOptimizationType(), cmdReader->getCategories(), alpha,
					cmdReader->estimateAlpha(), substParams, indelParams, dist);

			substParams = tme->getSubstitutionParameters();
			indelParams = tme->getIndelParameters();
			//if(cmdReader->estimateAlpha())
			//	alpha = tme->getAlpha();

			//cout << "Final indel params" << endl;
			//cout << (1.0-exp(indelParams[0]*-1.0*dist)) << "\t" << indelParams[1] << "\n";
			//cout << "Final subst params" << endl;
			//cout << substParams[0] << endl;

			delete tme;
*/
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


			treefile.open((string(cmdReader->getInputFileName()).append(".hmm.tree")).c_str(),ios::out);
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
