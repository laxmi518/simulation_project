#pragma once

class DtFrame
{
public:
	DtFrame(void);
public:
	~DtFrame(void);
private:
	std::vector<Integer> m_Step;
	std::vector<DeltaTime> m_Dt;
public:
	CCTStatusType Load(const std::string& FileName);
	Integer GetStep(Integer NextItr) const;
	const DeltaTime* GetDtTime(Integer NextItr) const;
	
};
