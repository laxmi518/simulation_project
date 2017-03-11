#pragma once
class Timer;
class CETProfile
{
private:
	Timer* m_pTimer;
	char* m_pFunctionName;
	static const char* m_pFilePath;
	double m_AverageTime;
	int m_TotalStep;
	
public:
	CETProfile();
	CETProfile(char* pFunctionName);										
	virtual ~CETProfile(void);
	static void setFilePath(const char * pFilePath);
	void writeTime(void);
	void CETProfile::DisplayTime(int Step);
	void CETProfile::SetAverageTime(int Step);
};
