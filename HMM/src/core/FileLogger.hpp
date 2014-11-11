#ifndef FILE_LOGGER_HPP
#define FILE_LOGGER_HPP


#include <fstream>
#include <string>
#include <vector>


namespace EBC
{


class FileLogger
{

//FIXME - create a nice log with static members
public:
	enum class logT { L_ERR, L_WARN, L_INF, L_DBG };

	static FileLogger& getLogger();

	void start(const std::string& fname)
	{
		//FIXME - file creation check
		instance.logFile.open (fname.c_str());
	}
	void stop()
	{
	    if (instance.logFile.is_open())
	    {
	        instance.logFile.close();
	    }
	}

	friend FileLogger &operator << (FileLogger &logger, const logT l_type) {

		switch (l_type) {
			case logT::L_ERR:
				logger.logFile << "[ERROR]: ";
				break;

			case logT::L_WARN:
				logger.logFile << "[WARNING]: ";
				break;
			default:
				logger.logFile << "[INFO]: ";
				break;
		}

		return logger;
	}

	template <typename T>
	friend FileLogger &operator << (FileLogger &logger, const T& param) {

		logger.logFile << param;
		return logger;
	}

	template <typename T>
	friend FileLogger &operator << (FileLogger &logger, const std::vector<T>& v)
	{
		for(unsigned int i = 0; i < v.size(); i++)
		{
			logger.logFile << v[i] << "\t";
		}
		if (v.size() != 0)
			logger.logFile << std::endl;
	}

private:
    FileLogger() {}
    FileLogger(const FileLogger& logger) {}
    FileLogger& operator = (const FileLogger& logger) {}
	static FileLogger instance;
	std::ofstream logFile;

};

}


#endif // FILELOGGER_HPP
