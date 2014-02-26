/*
 * SequenceElement.h
 *
 *  Created on: Jan 14, 2014
 *      Author: root
 */

#ifndef SEQUENCEELEMENT_H_
#define SEQUENCEELEMENT_H_

#include <vector>
using namespace std;


namespace EBC
{

class SequenceElement
{
protected:
	bool isGap;
	short matrixIndex;
	vector<short> alternativeIndexes;
public:
	SequenceElement(bool, short, short*);

	bool isIsGap() const
	{
		return isGap;
	}

	short getMatrixIndex() const
	{
		return matrixIndex;
	}
};

} /* namespace EBC */
#endif /* SEQUENCEELEMENT_H_ */
