/*
 * FileLogger.cpp

 *
 *  Created on: Nov 11, 2014
 *      Author: marcin
 */

#include "core/FileLogger.hpp"

namespace EBC
{

std::ofstream FileLogger::logFile;

FileLogger FileLogger::errL;
FileLogger FileLogger::wrnL;
FileLogger FileLogger::dbgL;
FileLogger FileLogger::dmpL;
FileLogger FileLogger::infL;

FileLogger& FileLogger::DumpLogger()
{
	if(dmpL.active)
		logFile << "[DUMP]\t";
	return dmpL;
}
FileLogger& FileLogger::DebugLogger()
{
	if(dbgL.active)
		logFile << "[DEBUG]\t";
	return dbgL;
}
FileLogger& FileLogger::ErrorLogger()
{
	if(errL.active)
		logFile << "[ERROR]\t";
	return errL;
}
FileLogger& FileLogger::InfoLogger()
{	if(infL.active)
		logFile << "[INFO]\t";
	return infL;
}
FileLogger& FileLogger::WarningLogger()
{	if(wrnL.active)
		logFile << "[WARN]\t";
	return wrnL;
}

}



