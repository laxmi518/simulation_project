#include "stdafx.h"
#include "DtFrame.h"

DtFrame::DtFrame(void)
{
}

DtFrame::~DtFrame(void)
{
}
CCTStatusType DtFrame::Load(const std::string& FileName)
{
	std::ifstream DtFile(FileName.c_str());
	if(!DtFile.is_open())
	{
		return CCT_FILEERR;
	}
	while(DtFile.is_open() && DtFile.good())
	{
		Integer Step;
		DeltaTime Dt;
		DtFile >> Step;
		DtFile >> Dt.Dt;
		DtFile >> Dt.CourantNumber;
		DtFile >> Dt.InflowStep;
		DtFile >> Dt.OutputStep;
		if(DtFile.good())
		{
			m_Step.push_back(Step);
			m_Dt.push_back(Dt);
		}
	}
	DtFile.close();
	return CCT_NOERR;
}
Integer DtFrame::GetStep(Integer NextItr) const
{
	if((NextItr>= 0)  && m_Dt.size() > NextItr)
	{
		const std::vector<Integer>::const_iterator StepItr = m_Step.begin() + NextItr;
		return *StepItr;
	}
	return -1;
}

const DeltaTime* DtFrame::GetDtTime(Integer NextItr) const
{
	if((NextItr>= 0) && m_Dt.size() > NextItr)
	{
		const std::vector<DeltaTime>::const_iterator DtItr = m_Dt.begin() + NextItr;
		return &(*DtItr);
		
	}
	return NULL;
}
