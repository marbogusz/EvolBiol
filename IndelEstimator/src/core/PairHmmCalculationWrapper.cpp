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
}

PairHmmCalculationWrapper::~PairHmmCalculationWrapper() {

}

double PairHmmCalculationWrapper::runIteration() {

	this->phmm->setDivergenceTimeAndCalculateModels(modelParams->getDivergenceTime(0));
	return this->phmm->runAlgorithm();
}

void PairHmmCalculationWrapper::setTargetHMM(EvolutionaryPairHMM* hmm) {
	this->phmm = hmm;
}

void EBC::PairHmmCalculationWrapper::setModelParameters(OptimizedModelParameters* mp) {
	this->modelParams = mp;
}

}

