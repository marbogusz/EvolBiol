/*
 * ViterbiPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef VITERBIFULLPAIRHMM_HPP_
#define VITERBIFULLPAIRHMM_HPP_

#include "hmm/EvolutionaryPairHMM.hpp"

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
	ViterbiPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, bool banding,
			SubstitutionModelBase* smdl, IndelModel* imdl, unsigned int bandPercentage, Definitions::DpMatrixType mt);

	virtual ~ViterbiPairHMM();

	double runAlgorithm();

	double getViterbiSubstitutionLikelihood();

};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */