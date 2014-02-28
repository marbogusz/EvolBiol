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

struct TraceStep
{
	//3 traceback elements
	bool vert;
	bool hor;
	bool diag;
	//score/probability
	double score;

	TraceStep() : vert(false), hor(false), 
		diag(false), score(0) {}
};

class DpMatrix
{

private:

	void allocateData();

protected:

	unsigned int xSize, ySize;

	double minVal;
	double maxVal;

	virtual void initializeData() = 0;

public:

	TraceStep ** matrixData;

	void traceback(string& seqA, string& seqB, std::pair<string,string>* alignment);

	DpMatrix(unsigned int xSize, unsigned int ySize);

	void setValue(unsigned int x,unsigned int y, double value);

	double valueAt(unsigned int i, unsigned int j);

	void setDiagonalAt(unsigned int i, unsigned int j);

	void setHorizontalAt(unsigned int i, unsigned int j);

	void setVerticalAt(unsigned int i, unsigned int j);

	void outputTrace();


};

} /* namespace EBC */
#endif /* EVALUATIONMATRIX_H_ */
