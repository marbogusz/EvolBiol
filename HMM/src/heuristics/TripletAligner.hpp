/*
 * TripletAligner.hpp
 *
 *  Created on: Sep 8, 2014
 *      Author: root
 */

#ifndef TRIPLETALIGNER_HPP_
#define TRIPLETALIGNER_HPP_

#include "core/Sequences.hpp"
#include "heuristics/GotohAlgorithm.hpp"
#include "core/SequenceElement.hpp"
#include <array>

using namespace std;

namespace EBC
{

class TripletAligner
{
private:

	array<SequenceElement, 3> tripletAlignment;

	Sequences* inputSeqs;
	unsigned int s1, s2,s3;

	void align();


public:
	TripletAligner(Sequences* inputSeq, array<unsigned int, 3>& triplet);

	const array<SequenceElement, 3>& getTripletAlignment() const
	{
		return tripletAlignment;
	}
};

} /* namespace EBC */

#endif /* TRIPLETALIGNER_HPP_ */
