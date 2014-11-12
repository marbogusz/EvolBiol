/*
 * PairHmmCalculationWrapper.cpp
 *
 *  Created on: Nov 12, 2014
 *      Author: marcin
 */

#include <core/PairHmmCalculationWrapper.hpp>

namespace EBC
{

PairHmmCalculationWrapper::PairHmmCalculationWrapper() {
	// TODO Auto-generated constructor stub

}

PairHmmCalculationWrapper::~PairHmmCalculationWrapper() {
	// TODO Auto-generated destructor stub
}

double PairHmmCalculationWrapper::runIteration() {

	this->phmm->setDivergenceTime(modelParams->getDivergenceTime(0));
	return this->phmm->runAlgorithm();
}

void PairHmmCalculationWrapper::setTargetHMM(EvolutionaryPairHMM* hmm) {
	this->phmm = hmm;
}

void EBC::PairHmmCalculationWrapper::setModelParameters(OptimizedModelParameters* mp) {
	this->modelParams = mp;
}

}

