/*
 * Dictionary.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef NUCDICTIONARY_H_
#define NUCDICTIONARY_H_

#include <vector>
#include <string>
#include <map>
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"
#include "core/Dictionary.hpp"

using namespace std;

namespace EBC
{
	class NucleotideDictionary : public Dictionary
	{
	public:
		NucleotideDictionary(Definitions::FrequencyScheme fs);
	};

}


#endif /* NUCDICTIONARY_H_ */
