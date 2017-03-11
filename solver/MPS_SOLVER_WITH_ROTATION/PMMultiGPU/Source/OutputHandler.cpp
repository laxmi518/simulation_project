#include "stdafx.h"
#include "OutputHandler.h"
#include "IO.h"

bool COutputHandler::m_IsAscii = true;
COutputHandler::COutputHandler(void)
{
}

COutputHandler::~COutputHandler(void)
{
}

//This is a thread
//must be executed as a thread
DWORD WINAPI COutputHandler::threadProc(LPVOID value)
{		//cout<<"i am printing";
	COutputFrame* OF = (COutputFrame*)value;
	OF->Output(COutputHandler::m_IsAscii);	
	return 0;
}

CCTStatusType COutputHandler::Output(COutputFrame* OutputFrame, bool ISAscii)
{
	m_IsAscii = ISAscii;
	if(!OutputFrame)
	{
		return CCT_PARAMERR;
	}
	OutputFrame->m_ThreadHandle = CreateThread(NULL,0,threadProc,OutputFrame,0,NULL);
	if(OutputFrame->m_ThreadHandle==NULL)
	{
		return CCT_ETCERR;
	}
	return CCT_NOERR;
}
CCTStatusType COutputFrame::Output(bool IsAscii)
{
	CCTStatusType Status =  IO::Save(m_FileName, m_ParticleNum, m_OutputParticles,m_ModelPosition,m_CurrentTime,IsAscii,
									m_aOutputParticleID,m_aOutputParticlePosition,m_aOutputParticleVelocity,m_aOutputParticlePressure,
									m_aOutputParticleDensity,m_aOutputParticleTemperature,	m_aOutputParticleKineticViscosity, 
									m_aOutputParticleSolidPhaseRate, m_aOutputParticleType,m_RotateAngle);
	m_Empty = true;
	CCT_ERROR_CHECK(Status);	
	return CCT_NOERR;
}
