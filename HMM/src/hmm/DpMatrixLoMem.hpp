/*
 * DpMatrixLoMem.h
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#ifndef DPLOMEMMATRIX_H_
#define DPLOMEMMATRIX_H_

#include <limits>
#include <iostream>
#include "hmm/DpMatrixBase.hpp"

using namespace std;

namespace EBC
{

class DpMatrixLoMem : public DpMatrixBase
{

protected:

	void allocateData();
	//Two lines of data!
	double* buffer[2];

	//First and second row pointers
	double* previousRow;
	double* currentRow;

	unsigned int currentRowIndex;
	unsigned int nextRowIndex;
	int previousRowIndex;

	void setValue(unsigned int col, double value);

	double valueAtColumn(unsigned int col);

	double valueAtLeft(unsigned int col);

	double valueAtTop(unsigned int col);

	double valueAtDiagonal(unsigned int col);

	inline void nextRow();

	void clear();

public:

	void setValue(unsigned int x,unsigned int y, double value);

	double valueAt(unsigned int i, unsigned int j);

	void setSrc(unsigned int i, unsigned int j, DpMatrixBase*);

	void setDiagonalAt(unsigned int i, unsigned int j);

	void setHorizontalAt(unsigned int i, unsigned int j);

	void setVerticalAt(unsigned int i, unsigned int j);

	void setWholeRow(unsigned int row, double value);

	void setWholeCol(unsigned int col, double value);

	void traceback(string& seq_a, string& seq_b, std::pair<string,string>* alignment) {}

	void tracebackRaw(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict, vector<std::pair<unsigned int, unsigned int> >&) {}


	DpMatrixLoMem(unsigned int xSize, unsigned int ySize);

	virtual ~DpMatrixLoMem();




};

} /* namespace EBC */
#endif /* DPREDUCEDMATRIX_H_ */
