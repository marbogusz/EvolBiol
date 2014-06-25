//============================================================================
// Name        : HMM.cpp
// Author      : Marcin Bogusz
// Version     :
// Copyright   :
//============================================================================

#include "CommandReader.hpp"
#include "Sequences.hpp"
#include "BasicViterbi.hpp"
#include "ForwardPairHMM.hpp"
#include "HmmException.hpp"
#include "PairwiseEstimator.hpp"
#include "BioNJ.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

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

		Sequences* inputSeqs = new Sequences(parser, cmdReader->getSequenceType());
		DEBUG("Creating the HMM");

		if(cmdReader->isForward())
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


			PairwiseEstimator* pe = new PairwiseEstimator(inputSeqs, cmdReader->getModelType() ,cmdReader->getIndelParams(),
					cmdReader->getSubstParams(), cmdReader->getOptimizationType(), cmdReader->getBanding(), cmdReader->getBandFactor(),
					cmdReader->getCategories(), cmdReader->getAlpha(), cmdReader->estimateAlpha());


			//string distFile = cmdReader->getInputFileName();
			//distFile.insert(0,"distances_");
			//distFile.replace(distFile.end()-3,distFile.end(),"phy");
			//stringstream ss;
			//pe->outputResults(ss);

			//ofstream of(distFile);
			//of << ss.str();
			//of.close();

			cerr << "running bionj\n";

			BioNJ nj(inputSeqs->getSequenceCount(), pe->getOptimizedTimes());
			nj.calculate();

			cout << cmdReader->getInputFileName() << endl;

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

		else if (cmdReader->isViterbi())
		{

		}

		else
		{
			throw HmmException("specify either V or F option");
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
