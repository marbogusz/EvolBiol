/*
 * DpMatrixFull.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: root
 */

#include "hmm/DpMatrixFull.hpp"
#include <iostream>
#include <cmath>

using namespace std;

void EBC::DpMatrixFull::allocateData()
{
	matrixData = new TraceStep*[xSize];
	for(unsigned int i=0; i<xSize; i++)
	{
		matrixData[i] = new TraceStep[ySize];
	}
}

EBC::DpMatrixFull::~DpMatrixFull()
{
	for(unsigned int i=0; i<xSize; i++)
	{
		delete[] matrixData[i];
	}
	delete[] matrixData;
}

EBC::DpMatrixFull::DpMatrixFull(unsigned int xS, unsigned int yS) : DpMatrixBase(xS,yS)
{
	this->allocateData();
}

void EBC::DpMatrixFull::setValue(unsigned int x, unsigned int y, double value)
{
	//TODO - check bounds
	matrixData[x][y].score = value;
}

void EBC::DpMatrixFull::setWholeRow(unsigned int row, double value)
{
	for (unsigned int i=0; i<ySize; i++)
	{
		matrixData[row][i].score = value;
		matrixData[row][i].hor = true;
		matrixData[row][i].src = this;
	}
}

void EBC::DpMatrixFull::setWholeCol(unsigned int col, double value)
{
	for (unsigned int i=0; i<xSize; i++)
	{
		matrixData[i][col].score = value;
		matrixData[i][col].vert = true;
		matrixData[i][col].src = this;
	}
}


double EBC::DpMatrixFull::valueAt(unsigned int i, unsigned int j)
{
	return matrixData[i][j].score;
}

void EBC::DpMatrixFull::outputTrace(unsigned int bound=0)
{
	unsigned int xl = bound !=0 ? bound : xSize;
	unsigned int yl = bound !=0 ? bound : ySize;

	for(unsigned int i=0; i < xl; i++)
	{
		for(unsigned int j=0; j < yl; j++)
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

void EBC::DpMatrixFull::outputValues(unsigned int bound=0)
{
	unsigned int xl = bound !=0 ? bound : xSize;
	unsigned int yl = bound !=0 ? bound : ySize;

	for(unsigned int i=0; i < xl; i++)
	{
		for(unsigned int j=0; j < yl; j++)
		{
			TraceStep& ts = matrixData[i][j];


			cout << (ts.score) << "\t";
		}
		cout << endl;
	}
}



void EBC::DpMatrixFull::setDiagonalAt(unsigned int i, unsigned int j)
{
	matrixData[i][j].diag = true;
}

void EBC::DpMatrixFull::setHorizontalAt(unsigned int i, unsigned int j)
{
	matrixData[i][j].hor = true;
}

void EBC::DpMatrixFull::setVerticalAt(unsigned int i, unsigned int j)
{
	matrixData[i][j].vert = true;
}

void EBC::DpMatrixFull::setSrc(unsigned int i, unsigned int j, DpMatrixBase* src)
{
	matrixData[i][j].src = static_cast<DpMatrixFull*>(src);
}

void EBC::DpMatrixFull::tracebackRaw(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict, vector<std::pair<unsigned int, unsigned int> >& al)
{
		unsigned int i = xSize-1;
		unsigned int j = ySize-1;

		string a1 = "";
		string a2 = "";

		DpMatrixFull* currentMat = this;

		while(i>0 || j >0)
		{
			if (currentMat->matrixData[i][j].diag)
			{
				al.push_back(std::make_pair(s1[i-1].getMatrixIndex(), s2[j-1].getMatrixIndex()));
				currentMat = currentMat->matrixData[i][j].src;

				a1 += dict->getSymbolAt(s1[i-1].getMatrixIndex());
				a2 += dict->getSymbolAt(s2[j-1].getMatrixIndex());

				i--;
				j--;
			}
			else if (currentMat->matrixData[i][j].hor)
			{
				al.push_back(std::make_pair(dict->getAlphabetSize(), s2[j-1].getMatrixIndex()));
				currentMat = currentMat->matrixData[i][j].src;

				a1 += '-';
				a2 += dict->getSymbolAt(s2[j-1].getMatrixIndex());
				j--;
			}
			//vert
			else
			{
				al.push_back(std::make_pair( s1[i-1].getMatrixIndex(),dict->getAlphabetSize()));
				currentMat = currentMat->matrixData[i][j].src;
				a1 += dict->getSymbolAt(s1[i-1].getMatrixIndex());
				a2 += '-';
				i--;
			}

		}

		//std::cerr << a1 << endl;
		//std::cerr << a2 << endl << endl;

}

void EBC::DpMatrixFull::traceback(string& seq_a, string& seq_b, std::pair<string,string>* alignment)
{
	unsigned int i = xSize-1;
	unsigned int j = ySize-1;

	DpMatrixFull* currentMat = this;

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













