#include "std.h"
#include "Logger.h"
#include <string>
using namespace SPVLog;

//#ifndef UNICODE
//typedef std::string String 
//#else
//typedef std::wstring String 
//#endif
namespace
{
	std::wstring GetFolderPath(const std::wstring& Filepath)
	//std::wstring GetFolderPath(const TCHAR Filepath)
	{
		std::wstring::size_type pos = Filepath.find_last_of(L"\\");
		if(std::wstring::npos == pos || 0 == pos)
		{
			return L"";
		}
		return Filepath.substr(0, pos);
	}

	std::wstring GetExeFolderPath()
	{
		static const int BUF_SIZE = 1024;
		TCHAR path[BUF_SIZE];
		::GetModuleFileName(NULL, path, sizeof(path));
		return GetFolderPath(path);
	}

	std::wstring GetDefaultLogFilepath()
	{
		return GetExeFolderPath() + L"\\spv.log";
	}

	void widen(const std::string &src, std::wstring &dest)
	{
		wchar_t *wcs = new wchar_t[src.length() + 1];
		mbstowcs(wcs, src.c_str(), src.length() + 1);
		dest = wcs;
		delete [] wcs;
	}
}

CSPVLogger* CSPVLogger::m_Instance = NULL;

CSPVLogger::CSPVLogger()
	: m_LogFilepath(L"")
{
	m_LogFilepath = GetDefaultLogFilepath();
}
CSPVLogger::~CSPVLogger()
{
}

CSPVLogger* CSPVLogger::GetInstance()
{
	if(!m_Instance)
	{
		m_Instance = new CSPVLogger();
	}
	return m_Instance;
}

bool CSPVLogger::Write(const CSPVAbstractLog& Formatter) const
{
	time_t rawtime;
	time ( &rawtime );
	struct tm* timeinfo = localtime( &rawtime );

	std::ofstream OP(m_LogFilepath.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
	std::ostream& Out = OP ? OP : std::cout;
	
	if(OP)
	{
		OP << "Time: " << asctime(timeinfo) << Formatter.Out() << "****************************************************************************\n";
		OP.close();
	}
	else
	{
		std::wstring Output;
		widen(Formatter.Out(), Output);
		std::wcout << L"error : fail to open log, " << m_LogFilepath << L"\n";
		std::wcout << L"  Time: " << asctime(timeinfo) << Output << L"****************************************************************************\n";
	}
	return true;
}

bool CSPVLogger::SetFilepath(const std::wstring& Filepath)
{
	m_LogFilepath = Filepath;
	return true;
}

CSPVAbstractLog::CSPVAbstractLog(LogType Type)
{}
CSPVAbstractLog::~CSPVAbstractLog()
{}

CSPVSystemInfoLog::~CSPVSystemInfoLog()
{}

std::string CSPVSystemInfoLog::Out() const
{
	std::stringstream SystemInformation;

	SYSTEM_INFO SI;
	GetSystemInfo(&SI);
	SystemInformation	<<"No Of Processor: "			<< SI.dwNumberOfProcessors		<<"\n"
						<<"Processor Type: "			<< SI.dwProcessorType			<<"\n"
						<<"Processor Architecture: "	<< SI.wProcessorArchitecture	<<"\n";
	return SystemInformation.str();
}

CSPVDeviceInfoLog::~CSPVDeviceInfoLog()
{}
std::string CSPVDeviceInfoLog::Out() const
{
	std::stringstream deviceInfo;
	int dev,driverVersion = 0,runtimeVersion = 0, deviceCount = 0;
	cudaDeviceProp deviceProp;

	cudaGetDeviceCount(&deviceCount);
	cudaDriverGetVersion(&driverVersion);
	deviceInfo  	<< "Device Count: "     << deviceCount									<< "\n"
					<< "Device Version: "   << driverVersion/1000							<< "\n\n";


	for(int dev = 0 ; dev < deviceCount ; dev++)
	{
		cudaGetDeviceProperties(&deviceProp, dev);
		cudaRuntimeGetVersion(&runtimeVersion);

		deviceInfo  << "Device Name: "		<< deviceProp.name								<< "\n"
					<< "Global Memory: "	<< (float)deviceProp.totalGlobalMem/1048576.0f	<< "\n"
					<< "Cuda Cores: "		<< ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) * deviceProp.multiProcessorCount << "\n"
					<< "Max Thread/Block: " << deviceProp.maxThreadsPerBlock					<< "\n\n";
	}
	
	return deviceInfo.str();
}

CSPVUserLogFormatter::~CSPVUserLogFormatter()
{}

CSPVUserLogFormatter::CSPVUserLogFormatter(
					CSPVAbstractLog::LogType Type,
					Integer GPUID,
					const std::string& msg_description,
					const std::string& file_name,
					Integer line_number)
					:m_Type(Type),
					 m_GPUID(GPUID),
					 m_MsgDescription(msg_description),
					 m_FileName(file_name),
					 m_LineNumber(line_number)
{
}

std::string CSPVUserLogFormatter::Out() const
{
	std::stringstream userLogMessage;
	userLogMessage
		<< "Type : " << m_Type << ", "
		<< "GPUID : " << m_GPUID << ", "
		<< "File Name : " << m_FileName << ", "
		<< "Line Number : " << m_LineNumber << ","
		<< "Message : " << m_MsgDescription <<"\n";
	return userLogMessage.str();
}

//#include "std.h"
//#include "Logger.h"
//#include <string>
//using namespace SPVLog;
//
//namespace
//{
//	std::wstring GetFolderPath(const std::wstring& Filepath)
//	{
//		std::wstring::size_type pos = Filepath.find_last_of(L"\\");
//		if(std::wstring::npos == pos || 0 == pos)
//		{
//			return L"";
//		}
//		return Filepath.substr(0, pos);
//	}
//
//	std::wstring GetExeFolderPath()
//	{
//		static const int BUF_SIZE = 1024;
//		TCHAR path[BUF_SIZE];
//		::GetModuleFileName(NULL, path, sizeof(path));
//		return GetFolderPath(path);
//	}
//
//	std::wstring GetDefaultLogFilepath()
//	{
//		return GetExeFolderPath() + L"\\spv.log";
//	}
//
//	void widen(const std::string &src, std::wstring &dest)
//	{
//		wchar_t *wcs = new wchar_t[src.length() + 1];
//		mbstowcs(wcs, src.c_str(), src.length() + 1);
//		dest = wcs;
//		delete [] wcs;
//	}
//}
//
//CSPVLogger* CSPVLogger::m_Instance = NULL;
//
//CSPVLogger::CSPVLogger()
//	: m_LogFilepath(L"")
//{
//	m_LogFilepath = GetDefaultLogFilepath();
//}
//CSPVLogger::~CSPVLogger()
//{
//}
//
//CSPVLogger* CSPVLogger::GetInstance()
//{
//	if(!m_Instance)
//	{
//		m_Instance = new CSPVLogger();
//	}
//	return m_Instance;
//}
//
//bool CSPVLogger::Write(const CSPVAbstractLog& Formatter) const
//{
//	time_t rawtime;
//	time ( &rawtime );
//	struct tm* timeinfo = localtime( &rawtime );
//
//	std::ofstream OP(m_LogFilepath.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
//	std::ostream& Out = OP ? OP : std::cout;
//	
//	if(OP)
//	{
//		OP << "Time: " << asctime(timeinfo) << Formatter.Out() << "****************************************************************************\n";
//		OP.close();
//	}
//	else
//	{
//		std::wstring Output;
//		widen(Formatter.Out(), Output);
//		std::wcout << L"error : fail to open log, " << m_LogFilepath << L"\n";
//		std::wcout << L"  Time: " << asctime(timeinfo) << Output << L"****************************************************************************\n";
//	}
//	return true;
//}
//
//bool CSPVLogger::SetFilepath(const std::wstring& Filepath)
//{
//	m_LogFilepath = Filepath;
//	return true;
//}
//
//CSPVAbstractLog::CSPVAbstractLog(LogType Type)
//{}
//CSPVAbstractLog::~CSPVAbstractLog()
//{}
//
//CSPVSystemInfoLog::~CSPVSystemInfoLog()
//{}
//
//std::string CSPVSystemInfoLog::Out() const
//{
//	std::stringstream SystemInformation;
//
//	SYSTEM_INFO SI;
//	GetSystemInfo(&SI);
//	SystemInformation	<<"No Of Processor: "			<< SI.dwNumberOfProcessors		<<"\n"
//						<<"Processor Type: "			<< SI.dwProcessorType			<<"\n"
//						<<"Processor Architecture: "	<< SI.wProcessorArchitecture	<<"\n";
//	return SystemInformation.str();
//}
//
//CSPVDeviceInfoLog::~CSPVDeviceInfoLog()
//{}
//std::string CSPVDeviceInfoLog::Out() const
//{
//	std::stringstream deviceInfo;
//	int dev,driverVersion = 0,runtimeVersion = 0, deviceCount = 0;
//	cudaDeviceProp deviceProp;
//
//	cudaGetDeviceCount(&deviceCount);
//	cudaDriverGetVersion(&driverVersion);
//	deviceInfo  	<< "Device Count: "     << deviceCount									<< "\n"
//					<< "Device Version: "   << driverVersion/1000							<< "\n\n";
//
//
//	for(int dev = 0 ; dev < deviceCount ; dev++)
//	{
//		cudaGetDeviceProperties(&deviceProp, dev);
//		cudaRuntimeGetVersion(&runtimeVersion);
//
//		deviceInfo  << "Device Name: "		<< deviceProp.name								<< "\n"
//					<< "Global Memory: "	<< (float)deviceProp.totalGlobalMem/1048576.0f	<< "\n"
//					<< "Cuda Cores: "		<< ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) * deviceProp.multiProcessorCount << "\n"
//					<< "Max Thread/Block: " << deviceProp.maxThreadsPerBlock					<< "\n\n";
//	}
//	
//	return deviceInfo.str();
//}
//
//CSPVUserLogFormatter::~CSPVUserLogFormatter()
//{}
//
//CSPVUserLogFormatter::CSPVUserLogFormatter(
//					CSPVAbstractLog::LogType Type,
//					Integer GPUID,
//					const std::string& msg_description,
//					const std::string& file_name,
//					Integer line_number)
//					:m_Type(Type),
//					 m_GPUID(GPUID),
//					 m_MsgDescription(msg_description),
//					 m_FileName(file_name),
//					 m_LineNumber(line_number)
//{
//}
//
//std::string CSPVUserLogFormatter::Out() const
//{
//	std::stringstream userLogMessage;
//	userLogMessage
//		<< "Type : " << m_Type << ", "
//		<< "GPUID : " << m_GPUID << ", "
//		<< "File Name : " << m_FileName << ", "
//		<< "Line Number : " << m_LineNumber << ","
//		<< "Message : " << m_MsgDescription <<"\n";
//	return userLogMessage.str();
//}
