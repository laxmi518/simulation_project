#include "ETProfile.h"
#include "Timer.h"
#include <fstream>
#include <iostream>

const char *CETProfile::m_pFilePath=".\\times.txt";
CETProfile::CETProfile()
:m_AverageTime(0.0)
,m_TotalStep(0)
{
	m_pTimer=new Timer();
	m_pTimer->start();
	m_pFunctionName = NULL;
}
CETProfile::CETProfile(char* pFunctionName)
:m_AverageTime(0.0)
,m_TotalStep(0)
{
	m_pTimer=new Timer();
	m_pTimer->start();
	m_pFunctionName=pFunctionName;	
}

CETProfile::~CETProfile(void)
{
	m_pTimer->stop();
	writeTime();
	delete[] m_pTimer;
}

void CETProfile::setFilePath(const char* pFilePath)
{
	m_pFilePath=pFilePath;
	std::fstream timerFile(m_pFilePath,std::ios::out);
	if(timerFile.is_open())
	{
		timerFile.close();
	}
}

void CETProfile::writeTime(void)
{
	std::fstream timerFile(m_pFilePath,std::ios::app);
	if(timerFile.is_open())
	{
		timerFile<<"Total Setp : "<<m_TotalStep<<std::endl;
		long Days,Hours,Minutes,Second,totalSecond;
		totalSecond = static_cast<long>(m_pTimer->getElapsedTimeInSec());
		timerFile<<"Total Elapsed Time: "<<totalSecond<<" Sec"<<std::endl;
		timerFile<<"Average Time: "<<m_AverageTime<<" Sec per step"<<std::endl;

		Minutes = totalSecond / 60;
		Second = totalSecond % 60;
		Hours = Minutes / 60;
		Minutes = Minutes % 60;
		Days = Hours / 24;
		Hours = Hours % 24;
		timerFile<< "Time D/H/M/S : " <<Days << " Days"<< Hours << " Hours" <<Minutes << " Minutes" << Second << " Seconds"<<std::endl;

		timerFile.close();
	}
	else
	{
		std::cout<<"Error writing in file";
	}
}
void CETProfile::DisplayTime(int Step)
{
	std::cout<<"STEP:"<<Step<<" :: TIME:"<<m_pTimer->getElapsedTimeInSec()<<" SEC"<<std::endl;
}
void CETProfile::SetAverageTime(int Step)
{
	m_TotalStep = Step;
	if(Step <= 0)
	{
		m_AverageTime = 0.0;
		return;
	}
	m_AverageTime = m_pTimer->getElapsedTimeInSec() / Step;	
}
