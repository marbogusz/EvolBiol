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

GotohAlgorithm::GotohAlgorithm(string& a, string& b) :
		seq_a(a), seq_b(b)
{
	cout << "Gotoh initializing ..." << endl;
	seqALen = seq_a.length();
	seqBLen = seq_b.length();

	cout << a << endl;
	cout << b << endl;

	xSize = seqALen + 1;
	ySize = seqBLen + 1;

	srand((unsigned)time(0));

	scores = new GotohScoringMatrix(0.5, -1.0, -0.2);
	scores->setStringPair(seq_a, seq_b);

	H = new GotohHMatrix(xSize, ySize, scores);
	V = new GotohVMatrix(xSize, ySize, scores);
	S = new GotohSMatrix(xSize, ySize, scores);

	H->setSMatrix(S);
	V->setSMatrix(S);
	//S gets the matrices set automatically
}

void GotohAlgorithm::setSequences(string& a, string& b)
{
	this->seq_a = a;
	this->seq_b = b;



	seqALen = seq_a.length();
	seqBLen = seq_b.length();
}

void GotohAlgorithm::run()
{
	processMatrices();
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
