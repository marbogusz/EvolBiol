/*
 * ParseException.h
 *
 *  Created on: Sep 23, 2013
 *      Author: mbogusz
 */

#ifndef PARSEEXCEPTION_H_
#define PARSEEXCEPTION_H_

#include <exception>
#include <string>

using namespace std;

namespace EBC {

class ProgramException: public std::exception {

private:
	string msg;

public:
	ProgramException();
	ProgramException(string message);
	virtual ~ProgramException() throw();
	virtual const char* what() const throw();

};

} /* namespace EBC */
#endif /* PARSEEXCEPTION_H_ */
