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
#include "heuristics/TripletModelEstimator.hpp"
#include <array>

using namespace std;
using namespace EBC;

int main(int argc, char ** argv) {

	//Set output Precision to 6
	cout << fixed << setprecision(6);
	cerr << fixed << setprecision(6);

	try
	{
		CommandReader* cmdReader = new CommandReader(argc, argv);
		DEBUG("Get parser" << endl);

		IParser* parser = cmdReader->getParser();
		DEBUG("Creating alignment");

		Sequences* inputSeqs = new Sequences(parser, cmdReader->getSequenceType(),cmdReader->isFixedAlignment());
		if (cmdReader->isMLE())
		{


			/*
			MlEstimator* me = new MlEstimator(inputSeqs, cmdReader->getModelType() ,cmdReader->getIndelParams(),
								cmdReader->getSubstParams(), cmdReader->getOptimizationType(), cmdReader->getCategories(),
								cmdReader->getAlpha(), cmdReader->estimateAlpha(),cmdReader->getDistance(), cmdReader->isFixedAlignment() == false);

			cout << cmdReader->getAlpha() << '\t' << cmdReader->getDistance() << '\n';
			 */

			//GuideTree gt(inputSeqs);
			//gt.constructTree();
			//auto samples = gt.getSampledTriplets();
			//array<unsigned int, 3> arr1 = samples[0];

			//GotohAlgorithm* algo = new GotohAlgorithm(inputSeqs->getRawSequenceAt(0),inputSeqs->getRawSequenceAt(1));
			//algo->run();
			//std::pair<string, string> al = algo->getAlignment();

			//DEBUG(al.first);
			//DEBUG(al.second);

			//TripletAligner tal(inputSeqs, {0,1,2});
			//tal.align();

			DEBUG("Creating TripletModelEstimator");

			TripletModelEstimator* tme = new TripletModelEstimator(inputSeqs, cmdReader->getModelType(),
					cmdReader->getOptimizationType(), cmdReader->getCategories(), cmdReader->getAlpha(),
					cmdReader->estimateAlpha());

			//tme->getModelParameters();

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

			cout << cmdReader->getAlpha() << '\t' << cmdReader->getDistance() << '\t';

			PairwiseEstimator* pe = new PairwiseEstimator(cmdReader->getAlgorithmType(), inputSeqs, cmdReader->getModelType() ,cmdReader->getIndelParams(),
					cmdReader->getSubstParams(), cmdReader->getOptimizationType(), cmdReader->getBanding(), cmdReader->getBandFactor(),
					cmdReader->getCategories(), cmdReader->getAlpha(), cmdReader->estimateAlpha(),cmdReader->getDistance());


			//string distFile = cmdReader->getInputFileName();
			//distFile.insert(0,"distances_");
			//distFile.replace(distFile.end()-3,distFile.end(),"phy");
			//stringstream ss;
			//pe->outputResults(ss);

			//ofstream of(distFile);
			//of << ss.str();
			//of.close();



			//cerr << "running bionj\n";

			//change bionj init here!
			//BioNJ nj(inputSeqs->getSequenceCount(), pe->getOptimizedTimes());
			//nj.calculate();

			//scerr << cmdReader->getInputFileName() << endl;

			//double* estimatedParams = fwdHMM->getMlParameters();

			//string vitFile = cmdReader->getInputFileName();

			//for (unsigned int i = 0; i< fwdHMM->getTotalParameters(); i++)
			//{
			//	cout << estimatedParams[i] << "\t";
			//}
			//cout <<  vitFile << endl;
/*
			if(cmdReader->isOutputViterbiAlignment())
			{

				string vitFile = cmdReader->getInputFileName();
				vitFile.insert(0,"viterbi_");
				vitFile.replace(vitFile.end()-3,vitFile.end(),"fas");
				vector<double> a,b;
				stringstream ss;
				BasicViterbi* bv = new BasicViterbi(inputSeqs, cmdReader->getModelType(),a,0,b, cmdReader->getCategories(), cmdReader->getAlpha(), estimatedParams);
				bv->runViterbiAlgorithm();
				bv->getResults(ss);

				ofstream of(vitFile);
				of << ss.str();
				of.close();

				delete bv;

			}

			delete fwdHMM;
*/
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
