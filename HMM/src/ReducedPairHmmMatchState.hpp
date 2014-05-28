/*
 * ReducedPairHmmMatchState.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: root
 */

#ifndef REDUCEDPAIRHMMMATCHSTATE_HPP_
#define REDUCEDPAIRHMMMATCHSTATE_HPP_

#include "DpReducedMatrix.hpp"
#include "PairHmmStateBase.hpp"
#include <iomanip>

namespace EBC
{

class ReducedPairHmmMatchState: public DpReducedMatrix, public PairHmmStateBase
{
public:
	ReducedPairHmmMatchState(unsigned int xSize, unsigned int ySize);

	virtual ~ReducedPairHmmMatchState();

	void initializeData();

	inline void outputRow()
	{
		cout << fixed << setprecision(1);
		for (int i=0; i < ySize; i++ )
			if (currentRow[i] < -9999.9)
			{
				cout << "-10k" << "\t";
			}
			else
				cout << currentRow[i] << "\t";
		cout << endl;
		cout << fixed << setprecision(6);
	}

	inline void nextRow()
	{
		//DEBUG("MC");
		//DEBUGV(currentRow, ySize);
		double* tmp = previousRow;
		previousRow = currentRow;
		currentRow = tmp;
		currentRow[0] = minVal;
		std::fill(currentRow,currentRow+ySize, -10000.0);
	}

};

} /* namespace EBC */

#endif /* REDUCEDPAIRHMMMATCHSTATE_HPP_ */
