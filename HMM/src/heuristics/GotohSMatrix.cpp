/*
 * GotohSMatrix.cpp
 *
 *  Created on: Feb 13, 2014
 *      Author: root
 */

#include "heuristics/GotohSMatrix.hpp"
#include "heuristics/GotohVMatrix.hpp"
#include "heuristics/GotohHMatrix.hpp"

namespace EBC
{
	void GotohSMatrix::initializeData()
	{
		this->matrixData[0][0].score = 0;

		for (int i = 1; i< this->xSize; i++)
		{
			this->matrixData[i][0].score = scoring->getGapOpening() + i * scoring->getGapExtension();
		}

		for (int j = 1; j< this->ySize; j++)
		{
			this->matrixData[0][j].score = scoring->getGapOpening() + j * scoring->getGapExtension();
		}
	}

	void GotohSMatrix::setHMatrix(GotohHMatrix* hm)
	{
		hMatrix = hm;
	}

	void GotohSMatrix::setVMatrix(GotohVMatrix* vm)
	{
		vMatrix = vm;
	}

	void GotohSMatrix::traceback()
	{
		/*
		EvaluationMatrix* currentMatrix = this;
		int i = this->xSize;
		int j = this->ySize;

		while (i !=0 || j!=0)
		{
			if (this->matrixData[i][j].diag)
			{
				cerr << "diag" << endl;
			}
			else if (this->matrixData[i][j].hor)
			{
				currMatrix = hMatrix;
				cerr << ho
			}
		}
		*/
	}

	double GotohSMatrix::evaluate(unsigned int i, unsigned int j)
	{
		double vScore = vMatrix->evaluate(i,j);
		double hScore = hMatrix->evaluate(i,j);
		double mScore = matrixData[i-1][j-1].score +
				scoring->getScoreByIndex(i-1,j-1);

		TraceStep* trace = &matrixData[i][j];


		if(mScore >= vScore && mScore >= hScore)
		{
			trace->diag = true;
			trace->score = mScore;
			trace->src = this;
			return mScore;
		}
		else if (vScore >= hScore)
		{
			trace->vert = true;
			trace->score = vScore;
			trace->src = vMatrix;
			return vScore;
		}
		else
		{
			trace->hor = true;
			trace->score = hScore;
			trace->src = hMatrix;
			return hScore;
		}
	}

	GotohSMatrix::GotohSMatrix(unsigned int xSize, unsigned int ySize, GotohScoringMatrix* sMat) :
			EvaluationMatrix(xSize, ySize, sMat)
	{

		initializeData();
	}
} /* namespace EBC */
