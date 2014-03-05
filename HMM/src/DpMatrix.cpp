/*
 * DpMatrix.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#include "DpMatrix.hpp"
#include <iostream>

using namespace std;

void EBC::DpMatrix::allocateData()
{
	matrixData = new TraceStep*[xSize];
	for(int i=0; i<xSize; i++)
	{
		matrixData[i] = new TraceStep[ySize];
	}
}

EBC::DpMatrix::~DpMatrix()
{
	for(int i=0; i<xSize; i++)
	{
		delete[] matrixData[i];
	}
	delete matrixData;
}

EBC::DpMatrix::DpMatrix(unsigned int xS, unsigned int yS) :
		xSize(xS), ySize(yS)
{
	maxVal = std::numeric_limits<double>::max();
	minVal = -10000;//std::numeric_limits<double>::min();
	allocateData();
}

void EBC::DpMatrix::setValue(unsigned int x, unsigned int y, double value)
{
	//TODO - check bounds
	matrixData[x][y].score = value;
}

void EBC::DpMatrix::setWholeRow(unsigned int row, double value)
{
	for (int i=0; i<ySize; i++)
	{
		matrixData[row][i].score = value;
		matrixData[row][i].hor = true;
		matrixData[row][i].src = this;
	}
}

void EBC::DpMatrix::setWholeCol(unsigned int col, double value)
{
	for (int i=0; i<xSize; i++)
	{
		matrixData[i][col].score = value;
		matrixData[i][col].vert = true;
		matrixData[i][col].src = this;
	}
}


double EBC::DpMatrix::valueAt(unsigned int i, unsigned int j)
{
	return matrixData[i][j].score;
}

void EBC::DpMatrix::outputTrace(unsigned int bound=0)
{
	unsigned int xl = bound !=0 ? bound : xSize;
	unsigned int yl = bound !=0 ? bound : ySize;

	for(int i=0; i < xl; i++)
	{
		for(int j=0; j < yl; j++)
		{
			TraceStep& ts = matrixData[i][j];
					if (ts.vert &&  !ts.hor && !ts.diag)
						cout << "|";
					else if (ts.vert && ts.hor && !ts.diag)
						cout << "J";
					else if (ts.vert && ts.hor && ts.diag)
						cout << "*";
					else if (!ts.vert && ts.hor && !ts.diag)
						cout << "-";
					else if (!ts.vert && !ts.hor && ts.diag)
						cout << "\\";
					else if (ts.vert && ! ts.hor && ts.diag)
						cout << "V";
					else if (!ts.vert && ts.hor && ts.diag)
						cout << "7";
					else cout << ".";
		}
		cout << endl;
	}
}

void EBC::DpMatrix::outputValues(unsigned int bound=0)
{
	unsigned int xl = bound !=0 ? bound : xSize;
	unsigned int yl = bound !=0 ? bound : ySize;

	for(int i=0; i < xl; i++)
	{
		for(int j=0; j < yl; j++)
		{
			TraceStep& ts = matrixData[i][j];


			cout << (ts.score) << "\t";
		}
		cout << endl;
	}
}



void EBC::DpMatrix::setDiagonalAt(unsigned int i, unsigned int j)
{
	matrixData[i][j].diag = true;
}

void EBC::DpMatrix::setHorizontalAt(unsigned int i, unsigned int j)
{
	matrixData[i][j].hor = true;
}

void EBC::DpMatrix::setVerticalAt(unsigned int i, unsigned int j)
{
	matrixData[i][j].vert = true;
}

void EBC::DpMatrix::setSrc(unsigned int i, unsigned int j, DpMatrix* src)
{
	matrixData[i][j].src = src;
}

void EBC::DpMatrix::traceback(string& seq_a, string& seq_b, std::pair<string,string>* alignment)
{
	unsigned int i = xSize-1;
	unsigned int j = ySize-1;

	DpMatrix* currentMat = this;

	while(i>0 || j >0)
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













