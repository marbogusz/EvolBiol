/*
 * StateTransitionEstimator.hpp
 *
 *  Created on: Oct 6, 2014
 *      Author: root
 */

#ifndef STATETRANSITIONESTIMATOR_HPP_
#define STATETRANSITIONESTIMATOR_HPP_

#include "core/IOptimizable.hpp"
#include "core/Optimizer.hpp"

namespace EBC
{

class StateTransitionEstimator : public IOptimizable
{
protected:

	Optimizer* bfgs;

public:
	StateTransitionEstimator();

	void runIteration();

	virtual ~StateTransitionEstimator();
};

} /* namespace EBC */

#endif /* STATETRANSITIONESTIMATOR_HPP_ */
