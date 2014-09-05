/*
 * GotohVMatrix.cpp
 *
 *  Created on: Feb 13, 2014
 *      Author: root
 */

#include "heuristics/GotohVMatrix.hpp"
#include "heuristics/GotohSMatrix.hpp"

namespace EBC
{

	void GotohVMatrix::initializeData()
	{
		for (int j = 0; j< this->ySize; j++)
		{
			this->matrixData[0][j].score = minVal;
		}
	}

	GotohVMatrix::GotohVMatrix(unsigned int xSize, unsigned int ySize, GotohScoringMatrix* sMat) :
			EvaluationMatrix(xSize, ySize, sMat)
	{
		initializeData();
	}

	void GotohVMatrix::setSMatrix(GotohSMatrix* SM)
	{
		sMatrix = SM;
		sMatrix->setVMatrix(this);
	}

	double GotohVMatrix::evaluate(unsigned int i, unsigned int j)
		{
			double vScore = matrixData[i-1][j].score + scoring->getGapExtension();
			double mScore = sMatrix->scoreAt(i-1,j) + scoring->getGapOpening() +
					scoring->getGapExtension();

			TraceStep* trace = &matrixData[i][j];

			if(mScore >= vScore)
			{
				trace->diag = true;
				trace->score = mScore;
				trace->src = sMatrix;
				return mScore;
			}
			else
			{
				trace->vert = true;
				trace->score = vScore;
				trace->src = this;
				return vScore;
			}
		}

} /* namespace EBC */


