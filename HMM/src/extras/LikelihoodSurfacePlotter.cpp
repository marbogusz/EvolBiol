/*
 * LikelihoodSurfacePlotter.cpp
 *
 *  Created on: Aug 3, 2015
 *      Author: marcin
 */

#include "extras/LikelihoodSurfacePlotter.hpp"

namespace EBC {

LikelihoodSurfacePlotter::LikelihoodSurfacePlotter() {
	// TODO Auto-generated constructor stub

}

void LikelihoodSurfacePlotter::setTargetHMM(EvolutionaryPairHMM* hmm) {
	this->phmm = hmm;
}

void LikelihoodSurfacePlotter::getLikelihoodSurface()
{
	FileLogger::Logger()  << "Time\tLnl\n";
	double step = 0.05;
	double maxDiv = 50;
	double time = step;
	double result;
	while (time < maxDiv){
		phmm->setDivergenceTimeAndCalculateModels(time);
		result = phmm->runAlgorithm();
		if (!std::isnan(result))
			FileLogger::Logger() << time << "\t" << result*-1.0 << "\n";
		time +=step;
	}
}

} /* namespace EBC */
