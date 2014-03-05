/*
 * DpMatrix.h
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#ifndef EVALUATIONMATRIX_H_
#define EVALUATIONMATRIX_H_

#include <limits>
#include <iostream>

using namespace std;

namespace EBC
{

class DpMatrix
{
	struct TraceStep
	{
		//3 traceback elements
		bool vert;
		bool hor;
		bool diag;
		//score/probability
		double score;
		DpMatrix* src;


		TraceStep() : vert(false), hor(false),
			diag(false), score(0), src(NULL) {}
	};


private:

	void allocateData();

protected:

	unsigned int xSize, ySize;

	double minVal;
	double maxVal;

	virtual void initializeData() = 0;

	virtual void setWholeRow(unsigned int row, double value);

	virtual void setWholeCol(unsigned int col, double value);


public:

	TraceStep ** matrixData;

	void traceback(string& seqA, string& seqB, std::pair<string,string>* alignment);

	DpMatrix(unsigned int xSize, unsigned int ySize);

	virtual ~DpMatrix();

	void setValue(unsigned int x,unsigned int y, double value);

	double valueAt(unsigned int i, unsigned int j);

	void setDiagonalAt(unsigned int i, unsigned int j);

	void setHorizontalAt(unsigned int i, unsigned int j);

	void setVerticalAt(unsigned int i, unsigned int j);

	void setSrc(unsigned int i, unsigned int j, DpMatrix*);

	void outputTrace(unsigned int);

	void outputValues(unsigned int);


};

} /* namespace EBC */
#endif /* EVALUATIONMATRIX_H_ */
