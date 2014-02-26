/*
 * PairHmmMatchState.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef PAIRHMMMATCHSTATE_HPP_
#define PAIRHMMMATCHSTATE_HPP_

#include "PairHmmState.hpp"

namespace EBC
{

class PairHmmMatchState: public EBC::PairHmmState
{
protected:

	void initializeData();

	double gapOpen;
	double gapExtension;

public:
	PairHmmMatchState(unsigned int x, unsigned int y, double gapOpening, double gapExtension);


	virtual ~PairHmmMatchState();
};

} /* namespace EBC */
#endif /* PAIRHMMMATCHSTATE_HPP_ */
