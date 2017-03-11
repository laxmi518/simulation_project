#include "DeviceSelector.h"
#include "OptionParser.h"
#include <boost/algorithm/string.hpp>
using namespace SystemUtility;

CDevicePriorityListMaker::CDevicePriorityListMaker(const COptionParser& Parser)
: m_RefParser(Parser)
{
}

CDevicePriorityListMaker::~CDevicePriorityListMaker()
{
}
void CDevicePriorityListMaker::SetDelimiter(const CDevicePriorityListDelimiter& delim)
{
	m_Delimiter = delim;
}
void CDevicePriorityListMaker::Make(const KeyValueArgs& Args, std::vector<int>& PriorityList) const
{
	PriorityList.clear();
	if(!Args.Has("d"))
	{
		return;
	}
	std::vector<std::string> Tokens;
	boost::algorithm::split(Tokens, Args["d"], boost::is_any_of(m_Delimiter.ToString()));
	if(Tokens.empty() || Tokens.front().empty())
	{
		return;
	}
	for(size_t i = 0; i < Tokens.size(); ++i)
	{
		PriorityList.push_back(atoi(Tokens[i].c_str()));
	}
}
std::string CDevicePriorityListMaker::DisplayVersion(const KeyValueArgs& Args) const
{
	//PriorityList.clear();
	if(!Args.Has("v"))
	{
		return "";
	}
	return "Current Version is 2.1";
	/*std::vector<std::string> Tokens;
	boost::algorithm::split(Tokens, Args["d"], boost::is_any_of(m_Delimiter.ToString()));
	if(Tokens.empty() || Tokens.front().empty())
	{
		return;
	}
	for(size_t i = 0; i < Tokens.size(); ++i)
	{
		PriorityList.push_back(atoi(Tokens[i].c_str()));
	}*/
}

CDeviceSelector::CDeviceSelector()
: m_IsInitialized(false)
{
}

CDeviceSelector::~CDeviceSelector()
{
}

bool CDeviceSelector::Initialize(int NumDevices)
{
	m_NumDevices = NumDevices;
	if(1 > m_NumDevices)
	{
		return false;
	}
	m_IsInitialized = false;
	m_DeviceId.clear();
	for(int i = 0; i < m_NumDevices; ++i)
	{
		m_DeviceId.push_back(i);
	}
	m_IsInitialized = true;
	return m_IsInitialized;
}

bool CDeviceSelector::IsInitialized() const
{
	return m_IsInitialized;
}

bool CDeviceSelector::SetPriority(std::vector<int>& Devices)
{
	if(!IsInitialized())
	{
		return false;
	}
	std::vector<int> DevicePriority(m_NumDevices, 0);
	std::vector<bool> AlreadyAssigned(m_NumDevices, false);

	int Count = 0;
	for(size_t i = 0; i < Devices.size(); ++i)
	{
		if(Devices[i] < 0 || Devices[i] >= m_NumDevices)
		{
			continue;
		}
		AlreadyAssigned[Devices[i]] = true;
		DevicePriority[Count] = Devices[i];
		++Count;
		if(m_NumDevices <= Count)
		{
			break;
		}
	}

	for(size_t i = 0; i < AlreadyAssigned.size(); ++i)
	{
		if(AlreadyAssigned[i])
		{
			continue;
		}
		AlreadyAssigned[i] = true;
		DevicePriority[Count] = m_DeviceId[i];
		++Count;
	}
	m_DeviceId = DevicePriority;
	return true;
}

bool CDeviceSelector::Empty() const
{
	return m_DeviceId.empty();
}

size_t CDeviceSelector::Size() const
{
	return m_DeviceId.size();
}

int CDeviceSelector::operator[](size_t Index) const
{
	return m_DeviceId[Index];
}