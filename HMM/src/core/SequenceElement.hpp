/*
 * SequenceElement.h
 *
 *  Created on: Jan 14, 2014
 *      Author: root
 */

#ifndef SEQUENCEELEMENT_H_
#define SEQUENCEELEMENT_H_

#include <vector>
#include <string>
using namespace std;


namespace EBC
{

class SequenceElement
{
//FIXME - this is potentially slow - rework

protected:
	bool isGap;
	unsigned char matrixIndex;
	string symbol;
	//vector<short> alternativeIndexes;
public:
	SequenceElement(bool, short, short*, string smbl);

	inline bool isIsGap() const
	{
		return isGap;
	}

	inline short getMatrixIndex()
	{
		return matrixIndex;
	}

	inline string& getSymbol()
	{
		return symbol;
	}
};

} /* namespace EBC */
#endif /* SEQUENCEELEMENT_H_ */
