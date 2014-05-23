//============================================================================
// Name        : HMM.cpp
// Author      : Marcin Bogusz
// Version     :
// Copyright   :
//============================================================================

#include "CommandReader.hpp"
#include "Sequences.hpp"
#include "ViterbiPairHMM.hpp"
#include "ForwardPairHMM.hpp"
#include "ParseException.hpp"
#include <iostream>
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

		Sequences* inputSeqs = new Sequences(parser);
		DEBUG("Creating the HMM");

		if(cmdReader->isForward())
		{
			//Check model

			//check evol model parameters provided

			//check indel parameters

			//ForwardPairHMM* fwdHMM = new ForwardPairHMM(inputSeqs,);

		}

		else if (cmdReader->isViterbi())
		{

		}

		else
		{
			throw ProgramException("specify either V or F option");
		}
		//ForwardPairHMM* epHMM = new ForwardPairHMM(inputSeqs);

		//ViterbiPairHMM* epHMM = new ViterbiPairHMM(inputSeqs);
		//epHMM->runViterbiAlgorithm();
		//epHMM->runForwardAlgorithm();
		//epHMM->getResults();
		//delete epHMM;

		cout << "FORWARD : ";
		//ForwardPairHMM* fwdHMM = new ForwardPairHMM(inputSeqs,true);
		//fwdHMM->summarize();
	}
	catch(ProgramException& pe)
	{
		cerr << pe.what();
	}
	catch(exception &ex)
	{
		cerr << ex.what();
	}
	return 0;
}
