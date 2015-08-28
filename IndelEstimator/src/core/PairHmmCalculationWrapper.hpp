/*
 * PairHmmCalcWrapper.hpp
 *
 *  Created on: Nov 12, 2014
 *      Author: marcin
 */

#ifndef CORE_PAIRHMMCALCULATIONWRAPPER_HPP_
#define CORE_PAIRHMMCALCULATIONWRAPPER_HPP_

#include "core/IOptimizable.hpp"
#include "hmm/EvolutionaryPairHMM.hpp"
#include "core/OptimizedModelParameters.hpp"

namespace EBC
{

class PairHmmCalculationWrapper : public IOptimizable
{
private:
	EvolutionaryPairHMM* phmm;
	OptimizedModelParameters* modelParams;

public:
	PairHmmCalculationWrapper();

	virtual ~PairHmmCalculationWrapper();

	double runIteration();

	void setTargetHMM(EvolutionaryPairHMM* hmm);
	void setModelParameters(OptimizedModelParameters* mp);


};

}
#endif /* CORE_PAIRHMMCALCULATIONWRAPPER_HPP_ */
