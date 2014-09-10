/*
 * EvaluationMatrix.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#include "heuristics/EvaluationMatrix.hpp"
#include <algorithm>

void EBC::EvaluationMatrix::allocateData()
{
	matrixData = new TraceStep*[xSize];
	for(int i=0; i<xSize; i++)
	{
		matrixData[i] = new TraceStep[ySize];
	}
}

EBC::EvaluationMatrix::EvaluationMatrix(unsigned int xS, unsigned int yS, GotohScoringMatrix* sMat) :
		xSize(xS), ySize(yS), scoring(sMat)
{
	maxVal = std::max(xSize,ySize) * scoring->getScoreByIndex(0,0);
	minVal = std::max(xSize,ySize) * std::min(scoring->getScoreByIndex(0,1), scoring->getGapOpening());
	allocateData();
}

void EBC::EvaluationMatrix::setValueScore(unsigned int x, unsigned int y, double value)
{
	//TODO - check bounds
	matrixData[x][y].score = value;
}

double EBC::EvaluationMatrix::scoreAt(unsigned int i, unsigned int j)
{
	return matrixData[i][j].score;
}

void EBC::EvaluationMatrix::traceback(string& seq_a, string& seq_b, std::pair<string,string>* alignment)
{
	unsigned int i = xSize-1;
	unsigned int j = ySize-1;

	EvaluationMatrix* currentMat = this;

	while(i>0 && j >0)
	{
		if (currentMat->matrixData[i][j].diag)
		{
			alignment->first += seq_a[i-1];
			alignment->second += seq_b[j-1];
			currentMat = currentMat->matrixData[i][j].src;
			i--;
			j--;
		}
		else if (currentMat->matrixData[i][j].hor)
		{
			alignment->second += seq_b[j-1];
			alignment->first += '-';
			currentMat = currentMat->matrixData[i][j].src;
			j--;
		}
		//vert
		else
		{
			alignment->first += seq_a[i-1];
			alignment->second += '-';
			currentMat = currentMat->matrixData[i][j].src;
			i--;
		}

	}
}

std::pair<string, string> EBC::EvaluationMatrix::getAlignment(string& seq_a,
		string& seq_b)
{
	string tmp1, tmp2;
	tmp1.reserve(seq_a.size() + (int)seq_a.size() * 0.2);
	tmp2.reserve(seq_b.size() + (int)seq_b.size() * 0.2);
	unsigned int i = xSize-1;
	unsigned int j = ySize-1;


	EvaluationMatrix* currentMat = this;

	while(i>0 || j >0)
	{
		if (currentMat->matrixData[i][j].diag)
		{
			tmp1 += seq_a[i-1];
			tmp2 += seq_b[j-1];
			currentMat = currentMat->matrixData[i][j].src;
			i--;
			j--;
		}
		else if (currentMat->matrixData[i][j].hor)
		{
			tmp2 += seq_b[j-1];
			tmp1 += '-';
			currentMat = currentMat->matrixData[i][j].src;
			j--;
		}
			//vert
		else
		{
			tmp1 += seq_a[i-1];
			tmp2 += '-';
			currentMat = currentMat->matrixData[i][j].src;
			i--;
		}

	}

	std::reverse(tmp1.begin(), tmp1.end());
	std::reverse(tmp2.begin(), tmp2.end());
	std::pair<string, string> alignment = std::make_pair(tmp1,tmp2);

	return alignment;

}
