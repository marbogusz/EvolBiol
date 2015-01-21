/*
 * StateTransitionMatrix.cpp
 *
 *  Created on: Oct 6, 2014
 *      Author: root
 */

#include <heuristics/StateTransitionML.hpp>
#include <cmath>

namespace EBC
{

void StateTransitionML::addSample(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, double weight)
{
	Definitions::StateId previousState;

		if ((*s1)[0]->isIsGap())
			previousState = Definitions::StateId::Delete;
		else if((*s2)[0]->isIsGap())
			previousState = Definitions::StateId::Insert;
		else
			previousState = Definitions::StateId::Match;



		for(int pos = 1; pos < s1->size(); pos++)
		{
			if ((*s1)[pos]->isIsGap())
			{
				//Delete
				if(previousState == Definitions::StateId::Match)
					counts[Definitions::StateId::Match][Definitions::StateId::Delete] += weight;
				if(previousState == Definitions::StateId::Insert)
					counts[Definitions::StateId::Insert][Definitions::StateId::Delete] += weight;
				if(previousState == Definitions::StateId::Delete)
					counts[Definitions::StateId::Delete][Definitions::StateId::Delete] += weight;

				previousState = Definitions::StateId::Delete;
			}
			else if((*s2)[pos]->isIsGap())
			{
				//Insert
				if(previousState == Definitions::StateId::Match)
					counts[Definitions::StateId::Match][Definitions::StateId::Insert] += weight;
				if(previousState == Definitions::StateId::Insert)
					counts[Definitions::StateId::Insert][Definitions::StateId::Insert] += weight;
				if(previousState == Definitions::StateId::Delete)
					counts[Definitions::StateId::Delete][Definitions::StateId::Insert] += weight;
				previousState = Definitions::StateId::Insert;
			}
			else
			{
				//Match
				if(previousState == Definitions::StateId::Match)
					counts[Definitions::StateId::Match][Definitions::StateId::Match] += weight;
				if(previousState == Definitions::StateId::Insert)
					counts[Definitions::StateId::Insert][Definitions::StateId::Match] += weight;
				if(previousState == Definitions::StateId::Delete)
					counts[Definitions::StateId::Delete][Definitions::StateId::Match] += weight;

				previousState = Definitions::StateId::Match;
			}
		}

}

StateTransitionML::StateTransitionML(IndelModel* im, double tme)
{
	this->tpb = new TransitionProbabilities(im);
	this->time = tme;
	tpb->setTime(time);

	for(int i = 0; i < Definitions::stateCount; i++)
			for(int j = 0; j<Definitions::stateCount; j++)
				counts[i][j] = 0.0;

	DEBUG("State Transition ML for time " << time);
}

StateTransitionML::~StateTransitionML()
{
	delete tpb;
}

void StateTransitionML::calculatePIs()
{
	//get the matrix and calculate
	pis[Definitions::StateId::Delete] = ((1.0-md[0][0])+(md[0][1]*(1.0-md[0][0]+md[1][0])/(md[1][1]-1.0-md[0][1])))/(((md[0][1]-md[2][1])*(1.0-md[0][0]+md[1][0])/(md[1][1]-1.0-md[0][1]))+md[2][0]-md[0][0]+1);
	pis[Definitions::StateId::Insert] = ((pis[Definitions::StateId::Delete]*(md[0][1]-md[2][1]))-md[0][1])/(md[1][1]-1.0-md[0][1]);
	pis[Definitions::StateId::Match] = 1.0 -pis[Definitions::StateId::Insert] - pis[Definitions::StateId::Delete];

}

void StateTransitionML::calculateParameters()
{
	//Match, insert, delete
	md[0][0] = 1.0-2*g;
	md[1][1] = e+((1.0-e)*g);
	md[2][2] = e+((1.0-e)*g);
	md[0][1] = g;
	md[0][2] = g;

	md[1][0] = (1.0-e)*(1-2*g);
	md[2][0] = (1.0-e)*(1-2*g);

	md[2][1] = (1.0-e)*g;
	md[1][2] = (1.0-e)*g;
}


double StateTransitionML::getLnL()
{
	double lnl  = 0;
	//set matrix based on gap probs
	tpb->calculate();
	e = tpb->getGapExtension();
	g = tpb->getGapOpening();

	calculateParameters();
	calculatePIs();

	for(int i = 0; i < Definitions::stateCount; i++)
			for(int j = 0; j<Definitions::stateCount; j++)
			{
				//go through site patterns
				lnl += counts[i][j] * log(pis[i]*md[i][j]);
			}
	if (std::isnan(lnl))
	{
		ERROR("NAN extension " << e);
		ERROR("NAN opening " << g);
		ERROR("ERROR - EXITING WITHOUT DOING CALCLULATIONS due to wrong extension/opening probabilities : " << e << " " << g);
		//FIXME - remove this!!!!!
		//FIXME - fix this estimation, don't exit like this - only for test purposes
		exit(0);
	}

	return lnl;
	//return likelihood
}

} /* namespace EBC */
