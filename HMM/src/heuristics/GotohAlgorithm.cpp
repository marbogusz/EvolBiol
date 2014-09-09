/*
 * GotohAlgorithm.cpp
 *
 *  Created on: Sep 23, 2013
 *      Author: mbogusz
 */

#include "heuristics/GotohAlgorithm.hpp"

using namespace std;

namespace EBC
{

GotohAlgorithm::~GotohAlgorithm()
{
	// TODO Auto-generated destructor stub
}

GotohAlgorithm::GotohAlgorithm(string& a, string& b)
{
	H = NULL;
	V = NULL;
	S = NULL;

	srand((unsigned)time(0));

	setSequences(a,b);
}

void GotohAlgorithm::setDistance(double distance)
{
	if (scores != NULL)
		delete scores;
	else scores = new GotohScoringMatrix(distance);
}

void GotohAlgorithm::setSequences(string& a, string& b)
{
	this->seq_a = a;
	this->seq_b = b;

	seqALen = seq_a.length();
	seqBLen = seq_b.length();

	xSize = seqALen + 1;
	ySize = seqBLen + 1;

	//reset matrices!

	scores->setStringPair(seq_a, seq_b);

	if(H != NULL)
		delete H;
	if(V != NULL)
			delete V;
	if(S != NULL)
			delete S;


	H = new GotohHMatrix(xSize, ySize, scores);
	V = new GotohVMatrix(xSize, ySize, scores);
	S = new GotohSMatrix(xSize, ySize, scores);

	H->setSMatrix(S);
	V->setSMatrix(S);
}

void GotohAlgorithm::run()
{
	processMatrices();

	DEBUG("S score" << S->scoreAt(xSize-1, ySize-1));
	DEBUG("V score" << V->scoreAt(xSize-1, ySize-1));
	DEBUG("H score" << H->scoreAt(xSize-1, ySize-1));
}

void GotohAlgorithm::processMatrices()
{
	for(int i=1; i<xSize; i++)
	{
		for(int j=1; j<ySize; j++)
		{
			S->evaluate(i,j);
		}
	}
}

}

std::pair<string, string> EBC::GotohAlgorithm::getAlignment()
{
	return this->S->getAlignment(seq_a, seq_b);
}
