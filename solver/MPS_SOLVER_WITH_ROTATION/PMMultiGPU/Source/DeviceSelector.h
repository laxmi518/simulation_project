#pragma once
#include <vector>

namespace SystemUtility
{

class COptionParser;
struct KeyValueArgs;


class CDevicePriorityListDelimiter
{
	std::string m_Delimiter;

public:
	CDevicePriorityListDelimiter(const std::string& Delimiter = ",")
		: m_Delimiter(Delimiter)
	{};

	CDevicePriorityListDelimiter(const CDevicePriorityListDelimiter& obj)
		: m_Delimiter(obj.m_Delimiter)
	{};

	virtual ~CDevicePriorityListDelimiter()
	{};

	CDevicePriorityListDelimiter& operator=(const CDevicePriorityListDelimiter& obj)
	{
		m_Delimiter = obj.m_Delimiter;
		return *this;
	};

	const std::string& ToString() const
	{
		return m_Delimiter;
	};
};


class CDevicePriorityListMaker
{
	const COptionParser& m_RefParser;
	CDevicePriorityListDelimiter m_Delimiter;

public:
	CDevicePriorityListMaker(const COptionParser& Parser);
	virtual ~CDevicePriorityListMaker();

	void SetDelimiter(const CDevicePriorityListDelimiter& delim);
	void Make(const KeyValueArgs& Args, std::vector<int>& PriorityList) const;

	std::string DisplayVersion(const KeyValueArgs& Args) const;
};


class CDeviceSelector
{
	bool m_IsInitialized;
	int m_NumDevices;
	std::vector<int> m_DeviceId;

public:
	CDeviceSelector();
	virtual ~CDeviceSelector();

	bool Initialize(int NumDevices);
	bool IsInitialized() const;
	bool SetPriority(std::vector<int>& Devices);
	bool Empty() const;
	size_t Size() const;
	int operator[](size_t Index) const;
};

}