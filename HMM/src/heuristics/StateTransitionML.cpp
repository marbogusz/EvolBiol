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

void StateTransitionML::addSample(vector<unsigned char>* s1, vector<unsigned char>* s2)
{
	Definitions::StateId previousState;

	if ((*s1)[0] == isGap)
		previousState = Definitions::StateId::Delete;
	else if((*s2)[0] == isGap)
		previousState = Definitions::StateId::Insert;
	else
		previousState = Definitions::StateId::Match;

	firstState = previousState;



	for(int pos = 1; pos < s1->size(); pos++)
	{
		if ((*s1)[pos] == isGap)
		{
			//Delete
			if(previousState == Definitions::StateId::Match)
				counts[Definitions::StateId::Match][Definitions::StateId::Delete] += 1;
			if(previousState == Definitions::StateId::Insert)
				counts[Definitions::StateId::Insert][Definitions::StateId::Delete] += 1;
			if(previousState == Definitions::StateId::Delete)
				counts[Definitions::StateId::Delete][Definitions::StateId::Delete] += 1;

			previousState = Definitions::StateId::Delete;
		}
		else if((*s2)[pos] == isGap)
		{
			//Insert
			if(previousState == Definitions::StateId::Match)
				counts[Definitions::StateId::Match][Definitions::StateId::Insert] += 1;
			if(previousState == Definitions::StateId::Insert)
				counts[Definitions::StateId::Insert][Definitions::StateId::Insert] += 1;
			if(previousState == Definitions::StateId::Delete)
				counts[Definitions::StateId::Delete][Definitions::StateId::Insert] += 1;
			previousState = Definitions::StateId::Insert;
		}
		else
		{
			//Match
			if(previousState == Definitions::StateId::Match)
				counts[Definitions::StateId::Match][Definitions::StateId::Match] += 1;
			if(previousState == Definitions::StateId::Insert)
				counts[Definitions::StateId::Insert][Definitions::StateId::Match] += 1;
			if(previousState == Definitions::StateId::Delete)
				counts[Definitions::StateId::Delete][Definitions::StateId::Match] += 1;

			previousState = Definitions::StateId::Match;
		}
	}

}

StateTransitionML::StateTransitionML(IndelModel* im, double tme, unsigned char gap, bool stEq) : isGap(gap), useStateEq(stEq)
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

	if (g < 0.0){
		ERROR("NAN opening " << g);
		g = Definitions::almostZero;
	}

	calculateParameters();

	calculatePIs();

	for(int i = 0; i < Definitions::stateCount; i++)
			for(int j = 0; j<Definitions::stateCount; j++)
			{
				if(counts[i][j] == 0)
					continue;
				//go through site patterns
				if (this->useStateEq)
					lnl += counts[i][j] * log(pis[i]*md[i][j]);
				else
					lnl += counts[i][j] * log(md[i][j]);
			}
	if(!useStateEq){
		lnl += log((md[0][firstState] * pis[0]) + (md[1][firstState] * pis[1]) + (md[2][firstState] * pis[2]));
	}

	//THIS SHOULD NEVER HAPPEN
	if (std::isnan(lnl))
	{
		ERROR("NAN extension " << e);
		ERROR("NAN opening " << g);
		ERROR("ERROR - EXITING WITHOUT DOING CALCLULATIONS due to wrong extension/opening probabilities : " << e << " " << g);
		//FIXME - remove this!!!!!
		//FIXME - fix this estimation, don't exit like this - only for test purposes
		cerr << "ERROR - EXITING WITHOUT DOING CALCLULATIONS due to wrong extension/opening probabilities : " << e << " " << g << endl;
		exit(0);
	}

	return lnl;
	//return likelihood
}

} /* namespace EBC */
