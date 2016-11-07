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

class Derivator
{
public:

	double uTime;
	vector<double> uIndelParams;
	vector<double> uSubstParams;

	vector<double> tempIndelParams;
	vector<double> tempSubstParams;
	vector<double> tempIndelParamsD;
	vector<double> tempSubstParamsD;


	SubstitutionModelBase* substModel;
	IndelModel* indelModel;
	ForwardPairHMM *hmm;

	double lDelta = 0.000001;
	double wDelta = 0.000001;

	void resetIndelParams(){
		tempIndelParams = uIndelParams;
		tempIndelParamsD = uIndelParams;
	}
	void resetSubstParams(){
		tempSubstParams = uSubstParams;
		tempSubstParamsD = uSubstParams;
	}

	void init(){
		for(int i = 0; i <=1; i++){
					if(uIndelParams[i] == 0.0){
						uIndelParams[i] = Definitions::almostZero;
					}
					if(uSubstParams[i] == 0.0){
						uSubstParams[i] = Definitions::almostZero;
					}
				}

				resetIndelParams();
				resetSubstParams();
	}

	double getGradient(){
		init();
		auto fpl0 = getDerivativeL(uIndelParams[0]);
		indelModel->setParameters(uIndelParams);
		auto fpw0 = getDerivativeW(uSubstParams[1]);
		cout  << fpl0 << '\t' << fpw0;
	}

	double getHessian(){
		init();

/*
		cout << uTime << endl;
		cout << uIndelParams[0] << endl;
		cout << uIndelParams[1] << endl;
		cout << uSubstParams[0] << endl;
		cout << uSubstParams[1] << endl;
*/

		auto fpl0 = getDerivativeL(uIndelParams[0]);
		auto fpld = getDerivativeL(uIndelParams[0] + lDelta);


		auto h11 = (fpld  - fpl0) / lDelta;


		tempIndelParams[0] = tempIndelParams[0] + lDelta;
		indelModel->setParameters(tempIndelParams);

		auto h21  = (getDerivativeW(uSubstParams[1]) - fpl0) /wDelta;

		resetIndelParams();
		indelModel->setParameters(uIndelParams);

		auto fpw0 = getDerivativeW(uSubstParams[1]);
		auto fpwd = getDerivativeW(uSubstParams[1] + wDelta);

		auto h22 = (fpwd  - fpw0) / wDelta;

		tempSubstParams[1] = tempSubstParams[1] + wDelta;
		substModel->setParameters(tempSubstParams);
		substModel->calculateModel();

		auto h12 = (getDerivativeL(uIndelParams[0]) - fpw0) /lDelta;

		cout << h11 << '\t' << h12 << '\t' << h21 << '\t' << h22;

		//cout << fpl0 << endl;
		//cout << fpld << endl;

		//cout << "MARCIN 11 " << h11 << endl;
		//cout << "MARCIN 21 " << h21 << endl;
		//cout << "MARCIN 22 " << h22 << endl;
		//cout << "MARCIN 12 " << h12 << endl;
	}


	double getDerivativeL(double point){
		tempIndelParams[0] = point;
		tempIndelParamsD[0] = point + lDelta;

		indelModel->setParameters(tempIndelParams);
		hmm->setDivergenceTimeAndCalculateModels(uTime);

		auto dp = -hmm->runAlgorithm();

		//cout << "lnl L point " << dp << endl;

		indelModel->setParameters(tempIndelParamsD);
		hmm->setDivergenceTimeAndCalculateModels(uTime);

		auto dpd = -hmm->runAlgorithm();

		//cout << "lnl L delta " << dpd << endl;
		resetIndelParams();
		return (dpd - dp)/lDelta;
	}

	double getDerivativeW(double point){
		tempSubstParams[1] = point;
		tempSubstParamsD[1] = point + wDelta;

		substModel->setParameters(tempSubstParams);
		substModel->calculateModel();
		hmm->setDivergenceTimeAndCalculateModels(uTime);

		auto dp = -hmm->runAlgorithm();

		substModel->setParameters(tempSubstParamsD);
		substModel->calculateModel();

		hmm->setDivergenceTimeAndCalculateModels(uTime);

		auto dpd = -hmm->runAlgorithm();
		resetSubstParams();
		return (dpd - dp)/wDelta;
	}

/*
			auto f_L_W = hmm->setDivergenceTimeAndCalculateModels(uTime);

			auto tmpSubstParams = uSubstParams;
			auto tmpIndelParams = uIndelParams;

			tmpIndelParams[0] = tmpIndelParams[0] + lDelta;
			indelModel->setParameters(tmpIndelParams);

			auto f_Ld_W = hmm->setDivergenceTimeAndCalculateModels(uTime);

			auto d1l = (f_Ld_W - f_L_W) / lDelta;

			tmpIndelParams[0] = tmpIndelParams[0] + lDelta;
			indelModel->setParameters(tmpIndelParams);

			auto f_L2d_W = hmm->setDivergenceTimeAndCalculateModels(uTime);

			//2nd lambda derivative
			auto d2l = (f_L2d_W - f_Ld_W) / lDelta;

			tmpIndelParams = uIndelParams;
			indelModel->setParameters(tmpIndelParams);

			tmpSubstParams[1] = tmpSubstParams[1] + wDelta;
*/

};

int main(int argc, char ** argv) {

	//cerr << "Starting dm\n";
	//Set output Precision to 2
	//FIXME - should normally be set to >= 6
	cout << fixed << setprecision(6);
	cerr << fixed << setprecision(6);


	//std::this_thread::sleep_for(std::chrono::milliseconds(30000));
	//cerr << " ELLO\n";
	//cout << " Should not be here!\n";
	//return 0;

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
		//FileLogger::InfoLogger().setCerr();

		INFO("Reading input sequences...");
		DEBUG("Creating alignment object...");

		Sequences* inputSeqs = new Sequences(parser, cmdReader->getSequenceType(),cmdReader->isFixedAlignment());

		//No indel or subst params!
		Dictionary* dict;
		SubstitutionModelBase* substModel;
		IndelModel* indelModel;
		Maths* maths;
		Definitions::AlgorithmType algorithm;

		maths = new Maths();
		dict = inputSeqs->getDictionary();

		double lDelta = 0.000001;
		double wDelta = 0.0001;

		//no gamma with codon models
		substModel = new CodonModel(dict, maths,1);
		indelModel = new NegativeBinomialGapModel();

		substModel->setObservedFrequencies(inputSeqs->getElementFrequencies());


		double uTime = cmdReader->getDistance();
		auto uIndelParams = cmdReader->getIndelParams();
		auto uSubstParams = cmdReader->getSubstParams();
		double lnl0 = 0.0;

		indelModel->setParameters(uIndelParams);
		substModel->setParameters(uSubstParams);
		substModel->calculateModel();



		ForwardPairHMM *hmm;

		std::pair<unsigned int, unsigned int> idxs = inputSeqs->getPairOfSequenceIndices(0);

		Band* band = NULL;

		hmm = new ForwardPairHMM(inputSeqs->getSequencesAt(idxs.first), inputSeqs->getSequencesAt(idxs.second),
									substModel, indelModel, Definitions::DpMatrixType::Full, band);



		Derivator dv;

		dv.hmm = hmm;
		dv.substModel = substModel;
		dv.indelModel = indelModel;
		dv.uTime = uTime;
		dv.uIndelParams = uIndelParams;
		dv.uSubstParams = uSubstParams;

		//dv.getHessian();
		dv.getGradient();


		delete hmm;
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









