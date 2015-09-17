//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

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
