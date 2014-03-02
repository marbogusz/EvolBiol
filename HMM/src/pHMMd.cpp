//============================================================================
// Name        : HMM.cpp
// Author      : Marcin Bogusz
// Version     :
// Copyright   :
//============================================================================

#include "CommandReader.hpp"
#include "Sequences.hpp"
//#include "ViterbiPairHMM.hpp"
#include "ForwardPairHMM.hpp"
#include "ParseException.hpp"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace EBC;

int main(int argc, char ** argv) {
	cout << fixed << setprecision(7);
	cerr << fixed << setprecision(7);

	cout << "Pair-HMM Distance estimation" << endl; // prints ML

	try
	{
		CommandReader* reader = new CommandReader(argc, argv);
		DEBUG("Get parser" << endl);
		IParser* parser = reader->getParser();
		DEBUG("Creating alignment");
		Sequences* inputSeqs = new Sequences(parser);
		DEBUG("Creating the HMM");
		ForwardPairHMM* epHMM = new ForwardPairHMM(inputSeqs);
		//epHMM->runViterbiAlgorithm();
		//epHMM->getResults();
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
