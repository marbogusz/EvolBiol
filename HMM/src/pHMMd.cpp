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

	//cerr << "Starting dm\n";
	//Set output Precision to 2
	//FIXME - should normally be set to >= 6
	cout << fixed << setprecision(8);
	cerr << fixed << setprecision(8);


	//std::this_thread::sleep_for(std::chrono::milliseconds(30000));
	//cerr << " ELLO\n";
	//cout << " Should not be here!\n";
	//return 0;

	try
	{
		//Get some time statistics
	    chrono::time_point<chrono::system_clock> start, end;
	    start = chrono::system_clock::now();

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

		//no gamma with codon models
		substModel = new CodonModel(dict, maths,1);
		indelModel = new NegativeBinomialGapModel();

		//indelModel->setParameters({Definitions::almostZero,0.5});

		substModel->setObservedFrequencies(inputSeqs->getElementFrequencies());


		double uTime = cmdReader->getDistance();
		auto uIndelParams = cmdReader->getIndelParams();
		auto uSubstParams = cmdReader->getSubstParams();
		double lnl0 = 0.0;

		indelModel->setParameters(uIndelParams);
		substModel->setParameters(uSubstParams);
		substModel->calculateModel();


		ForwardPairHMM *fhmm;
		BackwardPairHMM *bhmm;

		double totlnl;

		std::pair<unsigned int, unsigned int> idxs = inputSeqs->getPairOfSequenceIndices(0);

		Band* band = nullptr;

		fhmm = new ForwardPairHMM(inputSeqs->getSequencesAt(idxs.first), inputSeqs->getSequencesAt(idxs.second),
									substModel, indelModel, Definitions::DpMatrixType::Full, band);
		fhmm->setDivergenceTimeAndCalculateModels(uTime);

		totlnl = fhmm->runAlgorithm();

		bhmm = new BackwardPairHMM(inputSeqs->getSequencesAt(idxs.first), inputSeqs->getSequencesAt(idxs.second),
											substModel, indelModel, Definitions::DpMatrixType::Full, band);
		bhmm->setDivergenceTimeAndCalculateModels(uTime);

		bhmm->runAlgorithm();

		bhmm->calculatePosteriors(fhmm);

		bhmm->calculateMaximumPosteriorMatrix();

		//auto mp1 = bhmm->getMPAlignment();


		//DUMP("MPD aligment for sequence id " << idxs.first << " and " << idxs.second);
		//DUMP(mp1.first);
		//DUMP(mp1.second);

		//cout << inputSequences->getSequenceName(idxs.first) << "\t" <<  inputSequences->getSequenceName(idxs.second) << "\t" << bhmm->getAlignmentPosteriors(inputSequences->getAlignmentsAt(idxs.first), inputSequences->getAlignmentsAt(idxs.second));


		//cout << mp1.first << endl;
		//cout << mp1.second << endl;


		vector<SequenceElement*>* ssq1;
		vector<SequenceElement*>* ssq2;

/*
		cout << "FWD   lnL\t" << totlnl << endl;

		cout << "MPD" << endl;
		cout << "MPD" << endl;


*/

		ofstream mpdfile;

		mpdfile << fixed << setprecision(8);

		mpdfile.open ((string(cmdReader->getInputFileName()).append(".hmm.msa.samples")));

		std::map<double,pair<string, pair<vector<double>*, pair<vector<SequenceElement*>*, vector<SequenceElement*>*> >*> >  smap;

		pair<vector<double>*, pair<vector<SequenceElement*>*, vector<SequenceElement*>*> > combo1  = bhmm->getMPDPosteriors(dict);

		ssq1 = combo1.second.first;
		ssq2 = combo1.second.second;

		smap.insert(std::make_pair(-1.0 * bhmm->getAlignmentLikelihood(ssq1,ssq2),std::make_pair("MPD alignment" ,&combo1)));



		//cout << "MPD  lnL\t" << bhmm->getAlignmentLikelihood(ssq1,ssq2) << endl;
/*
		mpdfile << "MPD" <<endl;
		mpdfile << bhmm->getAlignmentLikelihood(ssq1,ssq2) << endl;

		for(auto val : *ssq1){
			mpdfile << val->getSymbol();
		}
		mpdfile << endl;
		for(auto val : *ssq2){
			mpdfile << val->getSymbol();
		}
		mpdfile << endl;

		for(auto val : *(combo1.first)){
			mpdfile << val << "\t";
		}
		mpdfile << endl;
*/

		pair<vector<double>*, pair<vector<SequenceElement*>*, vector<SequenceElement*>*> > combo2 = fhmm->getBestAlignmentFromForward(dict,bhmm);
		ssq1 = combo2.second.first;
		ssq2 = combo2.second.second;

		smap.insert(std::make_pair(-1.0 * bhmm->getAlignmentLikelihood(ssq1,ssq2),std::make_pair("BEST FWD alignment" ,&combo2)));

/*
		mpdfile << "Best fwd" << endl;
		mpdfile << bhmm->getAlignmentLikelihood(ssq1,ssq2) << endl;
*/
/*
		cout << "BEST FWD" << endl;
		cout << "BEST FWD" << endl;

		cout << "BEST fwd lnL\t" << bhmm->getAlignmentLikelihood(ssq1,ssq2) << endl;
*/
/*
		for(auto val : *ssq1){
			mpdfile << val->getSymbol();
		}
			mpdfile << endl;
		for(auto val : *ssq2){
			mpdfile << val->getSymbol();
		}
			mpdfile << endl;

		for(auto val : *(combo2.first)){
			mpdfile << val << "\t";
		}
*/
/*
		cout << "SAMPLES" << endl;
		cout << "SAMPLES" << endl;
*/




/*
		for(int s=0; s<200; s++){
			//cout << "SAMPLE " << s << endl;
			pair<vector<double>*, pair<vector<SequenceElement*>*, vector<SequenceElement*>*> > combo = fhmm->sampleAlignmentFromForward(dict,bhmm);
			ssq1 = combo.second.first;
			ssq2 = combo.second.second;

			string sname = "sample " + std::to_string(s+1);

			smap.insert(std::make_pair(-1.0 *  bhmm->getAlignmentLikelihood(ssq1,ssq2),std::make_pair(sname,&combo)));

		}
*/
		int ctr = 1;

		for(auto el : smap){
			if(ctr > 20){
				break;
			}
			mpdfile << el.second.first << endl;
			mpdfile << (-1.0 * el.first) << endl;

			ssq1 = (*(el.second.second)).second.first;
			ssq2 = (*(el.second.second)).second.second;

			for(auto val : *ssq1){
					mpdfile << val->getSymbol();
			}
					mpdfile << endl;
			for(auto val : *ssq2){
					mpdfile << val->getSymbol();
			}
			mpdfile << endl;

			for(auto val : *((*(el.second.second)).first)){
						mpdfile << val << "\t";
					}
			mpdfile << endl;
			ctr++;
		}
/*
		cout << "sample lnL\t" << bhmm->getAlignmentLikelihood(ssq1,ssq2) << endl;

							for(auto val : *ssq1){
									cout << val->getSymbol();
							}
								cout << endl;
							for(auto val : *ssq2){
								cout << val->getSymbol();
							}
							cout << endl;

							for(auto val : *(combo.first)){
								cout << val << "\t";
							}
							cout << endl;
							*/

		mpdfile.close();

		delete fhmm;
		delete bhmm;



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
