/*
 * ParseException.cpp
 *
 *  Created on: Sep 23, 2013
 *      Author: mbogusz
 */

#include <core/ProgramException.hpp>
#include "core/Definitions.hpp"

namespace EBC {

ProgramException::ProgramException()
{
	// General message
	msg = "Parse Exception";
}

ProgramException::~ProgramException() throw()
{
}

ProgramException::ProgramException(string message) : msg(message)
{
	ERROR(message);
}

const char* ProgramException::what() const throw ()
{
	return this->msg.c_str();
}

} /* namespace EBC */
