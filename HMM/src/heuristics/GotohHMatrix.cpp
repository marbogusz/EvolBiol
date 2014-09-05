/*
 * GotohHMatrix.cpp
 *
 *  Created on: Feb 13, 2014
 *      Author: root
 */

#include "heuristics/GotohHMatrix.hpp"
#include "heuristics/GotohSMatrix.hpp"

namespace EBC
{

	void GotohHMatrix::initializeData()
	{
		for (unsigned int i = 0; i< this->xSize; i++)
		{
			this->matrixData[i][0].score = minVal;
		}
	}

	GotohHMatrix::GotohHMatrix(unsigned int xSize, unsigned int ySize, GotohScoringMatrix* sMat) :
			EvaluationMatrix(xSize, ySize, sMat)
	{
		initializeData();
	}

	void GotohHMatrix::setSMatrix(GotohSMatrix* sm)
	{
		sMatrix = sm;
		sMatrix->setHMatrix(this);
	}

	double GotohHMatrix::evaluate(unsigned int i, unsigned int j)
		{
			double hScore = matrixData[i][j-1].score + scoring->getGapExtension();
			double mScore = sMatrix->scoreAt(i,j-1) + scoring->getGapOpening() +
					scoring->getGapExtension();

			TraceStep* trace = &matrixData[i][j];

			if(mScore >= hScore)
			{
				trace->diag = true;
				trace->score = mScore;
				trace->src = sMatrix;
				return mScore;
			}
			else
			{
				trace->hor = true;
				trace->score = hScore;
				trace->src = this;
				return hScore;
			}
		}
} /* namespace EBC */
