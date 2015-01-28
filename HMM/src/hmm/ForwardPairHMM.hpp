/*
 * ForwardPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef FORWARDPAIRHMM_HPP_
#define FORWARDPAIRHMM_HPP_

#include "hmm/EvolutionaryPairHMM.hpp"
#include "sampling/HMMPathSample.hpp"

namespace EBC
{

class ForwardPairHMM: public EBC::EvolutionaryPairHMM
{
friend class BackwardPairHMM;
protected:

	vector<double> userIndelParameters;
	vector<double> userSubstParameters;

public:
	ForwardPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2,
			SubstitutionModelBase* smdl, IndelModel* imdl,
			Definitions::DpMatrixType mt, Band* bandObj = nullptr, bool useEquilibriumProbabilities = false);

	virtual ~ForwardPairHMM();

	double runAlgorithm();

	pair<string,string> getBestAlignment(string&a, string& b);

	pair<string,string> sampleAlignment(string&a, string& b);

	//pair< double, pair<vector<SequenceElement*>*, vector<SequenceElement*>* >* > sampleAlignment(Dictionary* dictionary);

	pair<vector<unsigned char>*, vector<unsigned char>* >* sampleAlignment(Dictionary* dictionary, double& lnl);

	void sampleAlignment(HMMPathSample& sample);

	double calculateSampleLnL(HMMPathSample& sample);


};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
