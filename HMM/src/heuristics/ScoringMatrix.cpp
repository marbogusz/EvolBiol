/*
 * ScoringMatrix.cpp
 *
 *  Created on: Sep 25, 2013
 *      Author: mbogusz
 */

#include "heuristics/ScoringMatrix.hpp"
#include <cmath>

using namespace std;

namespace EBC
{

ScoringMatrix::ScoringMatrix()
{
	gapPenalty = Definitions::defaultGapPenalty;
}

ScoringMatrix::~ScoringMatrix()
{
	for(int i=0; i< matrixSize; i++)
			delete[] scores[i];
	delete[] scores;
}

ScoringMatrix::ScoringMatrix(unsigned int matrixSize, double distance, Dictionary* dict) : matrixSize(matrixSize), dict(dict)
{
	scores = new double*[matrixSize];
	for(int i=0; i< matrixSize; i++)
		scores[i] = new double[matrixSize];

	if(matrixSize == Definitions::nucleotideCount)
		scoresFromDistanceJC(distance);
	else if (matrixSize == Definitions::aminoacidCount)
	{
		scoresFromBLOSUM();
	}
}

void ScoringMatrix::scoresFromBLOSUM()
{
	for(int i=0; i< matrixSize; i++)
		std::copy(Definitions::blosum62[i], Definitions::blosum62[i]+matrixSize,scores[i]);
}

void ScoringMatrix::scoresFromDistanceJC(double distance)
{
	//log odds score
	//FIXME - assumption of a nucletides!
	//Change so it works with aminoacids

	double exP = exp(4*distance/-3.0);
	double exI = exp(distance * -1.0 * Definitions::initialLambda);
	double p0 = 1+3*exP;
	double p1 = 1-exP;
	double pI = 1-exI;
	p0 = log(p0);
	p1 = log(p1);
	pI = log(pI);

	for(int i=0; i<4;i++)
		for(int j=0;j<4; j++)
		{
			if (i==j)
				scores[i][j] = p0;
			else
				scores[i][j] = p1;
		}
	gapPenalty = pI;

	DEBUG("Scoring matrix match and mismatch scores : " << p0 << ", " << p1);
	DEBUG("Scoring matrix gap penalty : " << pI);

}

double ScoringMatrix::getScore(char& a, char& b)
{
	return scores[dict->getSymbolIndex(a)][dict->getSymbolIndex(b)];

}

double ScoringMatrix::getAlignmentScore(pair<string, string>* alignment)
{
	//FIXME - remove completely
	double score = 0;
	//for (int i = 0; i < alignment->first.size(); i++)
	//{
	//	score += getScore(alignment->first[i], alignment->second[i]);
	//}
	return score;
}

} /* namespace EBC */


