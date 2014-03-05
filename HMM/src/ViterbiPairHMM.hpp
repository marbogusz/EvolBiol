/*
 * ViterbiPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef VITERBIPAIRHMM_HPP_
#define VITERBIPAIRHMM_HPP_

#include "EvolutionaryPairHMM.hpp"

namespace EBC
{

class ViterbiPairHMM: public EBC::EvolutionaryPairHMM
{
protected:

	double getMax(double m, double x, double y, unsigned int i, unsigned int j, PairHmmState* state);

	virtual void initializeModels();


public:
	ViterbiPairHMM(Sequences* inputSeqs);

	virtual ~ViterbiPairHMM();

	void runViterbiAlgorithm();

	void getResults();
};

} /* namespace EBC */
#endif /* VITERBIPAIRHMM_HPP_ */
