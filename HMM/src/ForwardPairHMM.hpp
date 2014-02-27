/*
 * ForwardPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef FORWARDPAIRHMM_HPP_
#define FORWARDPAIRHMM_HPP_

#include "EvolutionaryPairHMM.hpp"

namespace EBC
{

class ForwardPairHMM: public EBC::EvolutionaryPairHMM
{
protected:

	virtual void initializeModels();

public:
	ForwardPairHMM(Sequences* inputSeqs);

	virtual ~ForwardPairHMM();

	void runForwardAlgorithm();
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
