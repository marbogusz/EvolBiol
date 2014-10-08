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


StateTransitionML::StateTransitionML(IndelModel* im, vector<SequenceElement>&s1, vector<SequenceElement>& s2, double tme)
{
	this->tpb = new TransitionProbabilities(im);
	this->time = tme;
	tpb->setTime(time);

	DEBUG("State Transition ML for time " << time);

	Definitions::StateId previousState;

	for(int i = 0; i < Definitions::stateCount; i++)
		for(int j = 0; j<Definitions::stateCount; j++)
			counts[i][j] = 0;

	if (s1[0].isIsGap())
		previousState = Definitions::StateId::Delete;
	else if(s2[0].isIsGap())
		previousState = Definitions::StateId::Insert;
	else
		previousState = Definitions::StateId::Match;



	for(int pos = 1; pos < s1.size(); pos++)
	{
		if (s1[pos].isIsGap())
		{
			//Delete
			if(previousState == Definitions::StateId::Match)
				counts[Definitions::StateId::Match][Definitions::StateId::Delete]++;
			if(previousState == Definitions::StateId::Insert)
				counts[Definitions::StateId::Insert][Definitions::StateId::Delete]++;
			if(previousState == Definitions::StateId::Delete)
				counts[Definitions::StateId::Delete][Definitions::StateId::Delete]++;

			previousState = Definitions::StateId::Delete;
		}
		else if(s2[pos].isIsGap())
		{
			//Insert
			if(previousState == Definitions::StateId::Match)
				counts[Definitions::StateId::Match][Definitions::StateId::Insert]++;
			if(previousState == Definitions::StateId::Insert)
				counts[Definitions::StateId::Insert][Definitions::StateId::Insert]++;
			if(previousState == Definitions::StateId::Delete)
				counts[Definitions::StateId::Delete][Definitions::StateId::Insert]++;
			previousState = Definitions::StateId::Insert;
		}
		else
		{
			//Match
			if(previousState == Definitions::StateId::Match)
				counts[Definitions::StateId::Match][Definitions::StateId::Match]++;
			if(previousState == Definitions::StateId::Insert)
				counts[Definitions::StateId::Insert][Definitions::StateId::Match]++;
			if(previousState == Definitions::StateId::Delete)
				counts[Definitions::StateId::Delete][Definitions::StateId::Match]++;

			previousState = Definitions::StateId::Match;
		}
	}

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
	return lnl;
	//return likelihood
}

} /* namespace EBC */
