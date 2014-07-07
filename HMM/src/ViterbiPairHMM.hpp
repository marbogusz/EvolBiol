/*
 * ViterbiPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef VITERBIFULLPAIRHMM_HPP_
#define VITERBIFULLPAIRHMM_HPP_

#include "EvolutionaryPairHMM.hpp"

namespace EBC
{

class ViterbiPairHMM: public EBC::EvolutionaryPairHMM
{

protected:

	vector<double> userIndelParameters;
	vector<double> userSubstParameters;

	double getMax(double m, double x, double y, unsigned int i, unsigned int j, PairwiseHmmStateBase* state);

public:
	ViterbiPairHMM(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict,  Definitions::ModelType model, bool banding,
			unsigned int bandPercentage, unsigned int rateCategories, Maths*);

	virtual ~ViterbiPairHMM();

	double runViterbiAlgorithm();

};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
