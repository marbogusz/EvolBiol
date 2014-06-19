/*
 * PairwiseHMMstate.h
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#ifndef PAIRWISEHMMSTATEBASE_H_
#define PAIRWISEHMMSTATEBASE_H_

#include "Dictionary.hpp"
#include "DpMatrixBase.hpp"

#include <map>
#include <algorithm>
#include <cmath>


namespace EBC
{

class PairwiseHmmStateBase
{
protected:

	double transFromMatch;
	double transFromInsert;
	double transFromDelete;

	DpMatrixBase* dpMatrix;



public:

	virtual void initializeData()=0;

	inline double getValueAt(unsigned int row, unsigned int column)
	{
		return this->dpMatrix->valueAt(row,column);
	}

	inline void setValueAt(unsigned int row, unsigned int column, double data)
	{
		this->dpMatrix->setValue(row,column,data);
	}

	virtual ~PairwiseHmmStateBase()
	{
		delete this->dpMatrix;
	}

	void setSourceMatrixPtr(unsigned int row, unsigned int column, PairwiseHmmStateBase* ptr)
	{
		dpMatrix->setSrc(row,column, ptr->getDpMatrix());
	}

	DpMatrixBase* getDpMatrix()
	{
		return this->dpMatrix;
	}

	virtual void setDirection(unsigned int, unsigned int)=0;

	//void addTransitionProbabilityFrom(PairHmmStateBase* state, double value);

	inline void setTransitionProbabilityFromMatch(double value)
	{
		transFromMatch = value;
	}

	inline void setTransitionProbabilityFromInsert(double value)
	{
		transFromInsert = value;
	}

	inline void setTransitionProbabilityFromDelete(double value)
	{
		transFromDelete = value;
	}

	inline double getTransitionProbabilityFromMatch()
	{
		return transFromMatch;
	}

	inline double getTransitionProbabilityFromInsert()
	{
		return transFromInsert;
	}
	inline double getTransitionProbabilityFromDelete()
	{
		return transFromDelete;
	}

	void traceback(string& seqA, string& seqB, std::pair<string,string>* alignment)
	{
		this->dpMatrix->traceback(seqA,seqB,alignment);
	}


};

} /* namespace EBC */
#endif /* PAIRWISEHMMSTATEBASE_H_ */