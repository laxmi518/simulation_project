#include "MacAddress.h"
#include <iostream>
#include <iomanip>
#include <sstream>

#include <Iphlpapi.h>
#include <Assert.h>
#pragma comment(lib, "iphlpapi.lib")


namespace
{
	template<typename INTEGER>
	INTEGER convert(BYTE value)
	{
		int converted = value & 0xff;
		return converted;
	}
}

MacAddressDisplayType::MacAddressDisplayType()
{
	Initialize();
}

MacAddressDisplayType::MacAddressDisplayType(const MacAddressDisplayType& obj)
	: nType(obj.nType)
	, aType(obj.aType)
	, hType(obj.hType)
{
}

MacAddressDisplayType::~MacAddressDisplayType()
{
}

void MacAddressDisplayType::Initialize()
{
	nType = HEX;
	aType = UPPER;
	hType = ON;
}

MacAddressDisplayType& MacAddressDisplayType::operator=(const MacAddressDisplayType& obj)
{
	nType = obj.nType;
	hType = obj.hType;
	aType = obj.aType;
	return *this;
}

CMacAddress::CMacAddress()
{
	Initialize();
}

CMacAddress::~CMacAddress()
{
}

void CMacAddress::Initialize()
{
	isLoaded = false;
	loadedAddressLength = 0;
	memset(address, 0, MAX_ADAPTER_ADDRESS_LENGTH);
	displayType.Initialize();
}

bool CMacAddress::Load()
{
	Initialize();
	IP_ADAPTER_INFO AdapterInfo[16];		// Allocate information for up to 16 NICs
	DWORD dwBufLen = sizeof(AdapterInfo);	// Save memory size of buffer
	DWORD dwStatus = GetAdaptersInfo(		// Call GetAdapterInfo
							AdapterInfo,	// [out] buffer to receive data
							&dwBufLen);		// [in] size of receive data buffer

	if(ERROR_SUCCESS != dwStatus)
	{
		return isLoaded;	// must be false because of Initialize()
	}
	isLoaded = true; 
	PIP_ADAPTER_INFO pAdapterInfo = AdapterInfo; // Contains pointer to
	loadedAddressLength = pAdapterInfo->AddressLength;
	for(size_t i = 0; i < MAX_ADAPTER_ADDRESS_LENGTH; ++i)
	{
		address[i] = pAdapterInfo->Address[i];
	}
	return isLoaded;
}
bool CMacAddress::LoadAllMacAddress()
{
	Initialize();
	int ValidNoOfMACAddress = 0;
	IP_ADAPTER_INFO AdapterInfo[16];			// Allocate information for up to 16 NICs
	DWORD dwBufLen = sizeof(AdapterInfo);		// Save the memory size of buffer

	DWORD dwStatus = GetAdaptersInfo(			// Call GetAdapterInfo
		AdapterInfo,							// [out] buffer to receive data
		&dwBufLen);								// [in] size of receive data buffer
	assert(dwStatus == ERROR_SUCCESS);			// Verify return value is valid, no buffer overflow

	if(ERROR_SUCCESS != dwStatus)
		return false;
	isLoaded = true;
	PIP_ADAPTER_INFO pAdapterInfo = AdapterInfo;// Contains pointer to current adapter info
	loadedAddressLength = pAdapterInfo->AddressLength;
	do {
		for(size_t i = 0; i < MAX_ADAPTER_ADDRESS_LENGTH; ++i)
		{
			address[i] = pAdapterInfo->Address[i];
		}
		m_ValidMacAddress.push_back(Address());
		pAdapterInfo = pAdapterInfo->Next;		// Progress through linked list

		printf("MAC Address: %02X-%02X-%02X-%02X-%02X-%02X\n", 
		address[0], address[1], address[2], address[3], address[4], address[5]);
	}
	while(pAdapterInfo);						// Terminate if last adapter
	return true;
}
std::vector<std::string> CMacAddress::GetAllMac()
{
	return m_ValidMacAddress;
}

std::string CMacAddress::Address() const
{
	if(!isLoaded)
	{
		return "";
	}
	std::stringstream ss;
	
	if(displayType.nType == MacAddressDisplayType::HEX)
	{
		if(MacAddressDisplayType::LOWER == displayType.aType)
		{
			ss << std::hex;
		}
		else
		{
			ss << std::hex << std::uppercase;
		}
		ss << std::setw(2) << std::setfill('0');
	}
	else
	{
		ss << std::dec;
	}
	ss  << convert<int>(address[0]);
	for(size_t i = 1; i < loadedAddressLength; ++i)
	{
		if(MacAddressDisplayType::ON == displayType.hType)
		{
			ss << "-";
		}
		//ss << convert<int>(address[i]);
		int num = convert<int>(address[i]);
		if(10 <= num && num <= 15)
		{
				ss << '0' << num;
		}
		else
		{
				ss << num;
		}
	}
	return ss.str();

	//New Change
	
}