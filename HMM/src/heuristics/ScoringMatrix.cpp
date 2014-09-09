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



ScoringMatrix::ScoringMatrix(double distance)
{
	scoresFromDistanceJC(distance);
}

void ScoringMatrix::scoresFromDistanceJC(double distance)
{
	//log odds score
	//FIXME - assumption of a nucletides!
	//Change so it works with aminoacids

	double exP = exp(4*distance/-3.0);
		double exI = exp(distance/-3.0);
		double p0 = 1+3*exP;
		double p1 = 1-exP;
		double pI = 1-exI;
		p0 = log(4*p0);
		p1 = log(4*p1);
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
}

double ScoringMatrix::getScore(char& a, char& b)
{
	int i,j;
	switch(a)
	{
		case 'A':
		case 'a':  	i = 0;
					break;
		case 'C':
		case 'c':  	i = 1;
					break;
		case 'G':
		case 'g':  	i = 2;
					break;
		case 'T':
		case 't':  	i = 3;
					break;
		default :  	return gapPenalty;
	}
	switch(b)
	{
		case 'A':
		case 'a':  j = 0;  break;

		case 'C':
		case 'c':  j = 1;  break;

		case 'G':
		case 'g':  j = 2;  break;

		case 'T':
		case 't':  j = 3;  break;

		default :  return gapPenalty;
	}

	return scores[i][j];

}

double ScoringMatrix::getAlignmentScore(pair<string, string>* alignment)
{
	double score = 0;
	for (int i = 0; i < alignment->first.size(); i++)
	{
		score += getScore(alignment->first[i], alignment->second[i]);
	}
	return score;
}

} /* namespace EBC */


