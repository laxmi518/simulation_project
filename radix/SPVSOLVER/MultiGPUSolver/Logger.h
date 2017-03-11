#pragma once

namespace SPVLog
{
	class CSPVAbstractLog;

	class CSPVLogger
	{
		CSPVLogger();
		static CSPVLogger* m_Instance;

	private:
		std::wstring m_LogFilepath;

	public:
		static CSPVLogger* GetInstance();
		~CSPVLogger();

		bool SetFilepath(const std::wstring& Filepath);
		bool Write(const CSPVAbstractLog& Formatter) const;
	};

	class CSPVAbstractLog
	{
	public:
		enum LogType
		{
			T_ERR,
			T_WARNING,
			T_INFO,
		};
		static std::string Str(LogType Type);

	protected:
		LogType m_Type;

	public:
		CSPVAbstractLog(LogType Type = T_INFO);
		virtual ~CSPVAbstractLog();
		virtual std::string Out() const = 0;

	public:
		LogType GetType() const { return m_Type; };
	};

	class CSPVSystemInfoLog : public CSPVAbstractLog
	{
	public:
		CSPVSystemInfoLog()
			: CSPVAbstractLog(CSPVAbstractLog::T_INFO)
		{};
		virtual ~CSPVSystemInfoLog();
		virtual std::string Out() const;
	};

	class CSPVDeviceInfoLog : public CSPVAbstractLog
	{
	public:
		CSPVDeviceInfoLog()
			: CSPVAbstractLog(CSPVAbstractLog::T_INFO)
		{};
		virtual ~CSPVDeviceInfoLog();
		virtual std::string Out() const;
	};

	class CSPVUserLogFormatter : public CSPVAbstractLog
	{
		Integer m_GPUID;
		std::string m_MsgDescription;
		std::string m_FileName;
		std::string m_FunctionName;
		Integer m_LineNumber;
		LogType m_Type ;

	public:
		CSPVUserLogFormatter(
			CSPVAbstractLog::LogType Type = CSPVAbstractLog::T_INFO,
			Integer GPUID = 0,
			const std::string& msg_description = "",
			const std::string& file_name = "",
			Integer line_number = 0);
		virtual ~CSPVUserLogFormatter();

		virtual std::string Out() const;
	};
}

#define CHAGE_LOG_FILE(Filepath) \
{ \
	SPVLog::CSPVLogger* Logger = SPVLog::CSPVLogger::GetInstance(); \
	if(Logger) \
	{ \
		Logger->SetFilepath(Filepath); \
	} \
	else \
	{ \
		std::cout << "error : fail to change log file path, " << Filepath << "\n"; \
	} \
} \

#define SYSTEM_INFO_LOG() \
{ \
	SPVLog::CSPVLogger* Logger = SPVLog::CSPVLogger::GetInstance(); \
	if(Logger) \
	{ \
		Logger->Write(SPVLog::CSPVSystemInfoLog()); \
	} \
	else \
	{ \
		SPVLog::CSPVSystemInfoLog Log; \
		std::cout << Log.Out() << "\n"; \
	} \
} \

#define DEVICE_INFO_LOG() \
{ \
	SPVLog::CSPVLogger* Logger = SPVLog::CSPVLogger::GetInstance(); \
	if(Logger) \
	{ \
		Logger->Write(SPVLog::CSPVDeviceInfoLog()); \
	} \
	else \
	{ \
		SPVLog::CSPVDeviceInfoLog Log; \
		std::cout << Log.Out() << "\n"; \
	} \
} \

#define USER_LOG(LogType, GpuId, Msg) \
{ \
	SPVLog::CSPVLogger* Logger = SPVLog::CSPVLogger::GetInstance(); \
	if(Logger) \
	{ \
		SPVLog::CSPVUserLogFormatter UserLog(LogType, GpuId, Msg, __FILE__, __LINE__); \
		Logger->Write(UserLog); \
	} \
	else \
	{ \
		SPVLog::CSPVSystemInfoLog Log; \
		std::cout << Log.Out() << "\n"; \
	} \
} \

#define ERR_TYPE SPVLog::CSPVAbstractLog::T_ERR

inline int ConvertSMVer2Cores(int major, int minor)
{
	// Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
	typedef struct {
		int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
		int Cores;
	} sSMtoCores;

	sSMtoCores nGpuArchCoresPerSM[] = 
	{ { 0x10,  8 }, // Tesla Generation (SM 1.0) G80 class
	  { 0x11,  8 }, // Tesla Generation (SM 1.1) G8x class
	  { 0x12,  8 }, // Tesla Generation (SM 1.2) G9x class
	  { 0x13,  8 }, // Tesla Generation (SM 1.3) GT200 class
	  { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
	  { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
	  {   -1, -1 }
	};

	int index = 0;
	while (nGpuArchCoresPerSM[index].SM != -1) {
		if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor) ) {
			return nGpuArchCoresPerSM[index].Cores;
		}
		index++;
	}
	//printf("MapSMtoCores SM %d.%d is undefined (please update to the latest SDK)!\n", major, minor);
	return -1;
}