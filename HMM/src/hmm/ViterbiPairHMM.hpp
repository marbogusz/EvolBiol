/*
 * ViterbiPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef VITERBIFULLPAIRHMM_HPP_
#define VITERBIFULLPAIRHMM_HPP_

#include "hmm/EvolutionaryPairHMM.hpp"
#include "sampling/HMMPathSample.hpp"

namespace EBC
{

class ViterbiPairHMM: public EBC::EvolutionaryPairHMM
{

protected:

	vector<double> userIndelParameters;
	vector<double> userSubstParameters;

	vector<std::pair<unsigned int, unsigned int> > alignment;

	double getMax(double m, double x, double y, unsigned int i, unsigned int j, PairwiseHmmStateBase* state);

public:
	ViterbiPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2,
			SubstitutionModelBase* smdl, IndelModel* imdl,
			Definitions::DpMatrixType mt = Definitions::DpMatrixType::Full, Band* bandObj = nullptr, bool useEquilibriumFreqs = false);

	virtual ~ViterbiPairHMM();

	double runAlgorithm();

	pair<string, string> getAlignment(string&a, string& b);

	void getAlignment(HMMPathSample& sample);

	double getViterbiSubstitutionLikelihood();

	pair<string,string> getBestAlignment(string&a, string& b);
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
