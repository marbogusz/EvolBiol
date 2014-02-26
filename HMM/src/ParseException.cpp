/*
 * ParseException.cpp
 *
 *  Created on: Sep 23, 2013
 *      Author: mbogusz
 */

#include "ParseException.hpp"

namespace EBC {

ProgramException::ProgramException()
{
	// TODO Auto-generated constructor stub
	// General message
	msg = "Parse Exception";
}

ProgramException::~ProgramException() throw()
{
	// TODO Auto-generated destructor stub
}

ProgramException::ProgramException(string message) : msg(message)
{
}

const char* ProgramException::what() const throw ()
{
	return this->msg.c_str();
}

} /* namespace EBC */
