/*
 * IOptimizable.h
 *
 *  Created on: Oct 6, 2014
 *      Author: root
 */

#ifndef IOPTIMIZABLE_H_
#define IOPTIMIZABLE_H_

namespace EBC
{

class IOptimizable
{
public:

	virtual double runIteration() = 0;
};

} /* namespace EBC */

#endif /* IOPTIMIZABLE_H_ */
