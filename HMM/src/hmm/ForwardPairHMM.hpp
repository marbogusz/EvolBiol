/*
 * ForwardPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef FORWARDPAIRHMM_HPP_
#define FORWARDPAIRHMM_HPP_

#include "hmm/EvolutionaryPairHMM.hpp"

namespace EBC
{

class ForwardPairHMM: public EBC::EvolutionaryPairHMM
{

protected:

	vector<double> userIndelParameters;
	vector<double> userSubstParameters;

public:
	ForwardPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, bool banding,
			SubstitutionModelBase* smdl, IndelModel* imdl, unsigned int bandPercentage, Definitions::DpMatrixType mt);

	virtual ~ForwardPairHMM();

	double runAlgorithm();

};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
