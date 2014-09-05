/*
 * EvaluationMatrix.h
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#ifndef EVALUATIONMATRIX_H_
#define EVALUATIONMATRIX_H_

#include "heuristics/GotohScoringMatrix.hpp"
#include <limits>

namespace EBC
{



class EvaluationMatrix
{

protected:
	struct TraceStep
	{
		bool vert;
		bool hor;
		bool diag;
		double score;
		EvaluationMatrix* src;

		TraceStep() : vert(false), hor(false), diag(false),
				score(std::numeric_limits<double>::min()), src(NULL){}
	};

private:

	void allocateData();

protected:

	unsigned int xSize, ySize;


	double minVal;
	double maxVal;

	GotohScoringMatrix* scoring;

	virtual void initializeData() = 0;

public:

	TraceStep ** matrixData;

	EvaluationMatrix(unsigned int xSize, unsigned int ySize, GotohScoringMatrix* sMat);

	void setValueScore(unsigned int x,unsigned int y, double value);

	unsigned int pathsAt(unsigned int i, unsigned int j);

	double scoreAt(unsigned int i, unsigned int j);

	virtual double evaluate(unsigned int i, unsigned int j) = 0;

	void traceback(string& seq_a, string& seq_b, std::pair<string,string>* alignment);


};

} /* namespace EBC */
#endif /* EVALUATIONMATRIX_H_ */
