/*
 * BandCalculator.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: root
 */

#include <heuristics/BandCalculator.hpp>

namespace EBC
{

BandCalculator::BandCalculator(vector<SequenceElement>& s1, vector<SequenceElement>& s2, SubstitutionModelBase* sm, IndelModel* im, double divergenceTime) :
		fwd(3), seq1(s1), seq2(s2), substModel(sm), indelModel(im), time(divergenceTime)
{
	// TODO Auto-generated constructor stub
	this->ptMatrix =  new PMatrixDouble(substModel);
	this->trProbs = new TransitionProbabilities(indelModel);
	array<double,3> multipliers = {0.5,1,1.5};
	unsigned int best = 0;
	double tmpRes = std::numeric_limits<double>::max();
	double lnl;

	for(unsigned int i = 0; i < fwd.size(); i++)
	{

		DEBUG(i << " with divergence time " << divergenceTime);
		fwd[i] = new ForwardPairHMM(seq1,seq2, false, substModel,indelModel, 0, Definitions::DpMatrixType::Full);
		fwd[i]->setDivergenceTime(time*multipliers[i]);
		lnl = fwd[i]->runAlgorithm();
		DEBUG(i << " with divergence time " << divergenceTime << " lnl " << lnl);
		if(lnl < tmpRes)
		{
			best = i;
			tmpRes = lnl;
		}
		//bwd[i] = new BackwardPairHMM(seq1,seq2, false, substModel,indelModel, 0, Definitions::DpMatrixType::Full);
		//bwd[i]->setDivergenceTime(time*multipliers[i]);
		//bwd[i]->runAlgorithm();
	}

	DEBUG("Best " << best);

	bwd =  new BackwardPairHMM(seq1,seq2, false, substModel,indelModel, 0, Definitions::DpMatrixType::Full);
	DEBUG("BWD Set time");
	bwd->setDivergenceTime(time*multipliers[best]);
	DEBUG("BWD Run");
	bwd->runAlgorithm();

	//combine fwd and bwd metrics into one!
	band = new Band(s2.size());
	bwd->calculatePosteriors(fwd[best]);
	band->processPosteriorProbabilities(bwd);
}

BandCalculator::~BandCalculator()
{
	for(auto ptr : fwd)
		delete ptr;
	delete bwd;
}

} /* namespace EBC */
