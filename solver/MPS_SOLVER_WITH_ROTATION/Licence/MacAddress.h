#pragma once
#include <Windows.h>
#include <iphlpapi.h>
#include <string>
#include <vector>

struct MacAddressDisplayType
{
	enum NumType
	{
		DEC = 0,
		HEX,
	};

	enum AlphabetType
	{
		LOWER = 0,
		UPPER,
	};

	enum HyphenType
	{
		OFF = 0,
		ON,
	};

	NumType nType;
	AlphabetType aType;
	HyphenType hType;

	MacAddressDisplayType();
	MacAddressDisplayType(const MacAddressDisplayType& obj);
	~MacAddressDisplayType();

	void Initialize();
	MacAddressDisplayType& operator=(const MacAddressDisplayType& obj);
};

class CMacAddress
{
	std::vector<std::string> m_ValidMacAddress;
	bool isLoaded;
	size_t loadedAddressLength;
	BYTE address[MAX_ADAPTER_ADDRESS_LENGTH];
	MacAddressDisplayType displayType;

public:
	CMacAddress();
	virtual ~CMacAddress();

	void Initialize();
	bool Load();	///< load mac address of this computer
	bool LoadAllMacAddress();
	
public:
	bool IsLoaded() const { return isLoaded; };
	std::string Address() const;						///< get mac address as string

public:
	void SetDisplayType(const MacAddressDisplayType& type) { displayType = type; };
	std::vector<std::string> GetAllMac();
};