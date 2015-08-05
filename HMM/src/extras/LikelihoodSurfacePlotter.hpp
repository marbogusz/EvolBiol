/*
 * LikelihoodSurfacePlotter.hpp
 *
 *  Created on: Aug 3, 2015
 *      Author: marcin
 */

#ifndef EXTRAS_LIKELIHOODSURFACEPLOTTER_HPP_
#define EXTRAS_LIKELIHOODSURFACEPLOTTER_HPP_

#include "hmm/EvolutionaryPairHMM.hpp"

namespace EBC {

class LikelihoodSurfacePlotter
{
private:
	EvolutionaryPairHMM* phmm;

public:
	LikelihoodSurfacePlotter();

	void setTargetHMM(EvolutionaryPairHMM* hmm);

	void getLikelihoodSurface();
};

} /* namespace EBC */

#endif /* EXTRAS_LIKELIHOODSURFACEPLOTTER_HPP_ */
