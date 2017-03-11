// ClientTest.cpp : Defines the entry point for the console application.
//
//#include "CSmartReader.h"
#include <sstream>
#include <iostream>
#include "utils.h"
#include <vector>
#include "resource.h"
#include "Licence.h"
#include "MacAddress.h"
#include "../PMMultiGPU/Source/File.h"
#include <Windows.h>
#include <Iphlpapi.h>
#include <Assert.h>
#pragma comment(lib, "iphlpapi.lib")

#define LICENSE_OFF

namespace
{
	/**
	* get license file path that must be in same folder where
	* executable is.
	*/
	std::wstring GetLicenseFilePath()
	{
		std::wstring ExecutableFolder = File::GetExeFolderPath();
		std::wstring LicenseFilePath = ExecutableFolder + L"\\License.dat";
		return LicenseFilePath;
	}
}

CLicenceValidate::CLicenceValidate()
{
}

CLicenceValidate::~CLicenceValidate()
{
}

bool CLicenceValidate::LicenceCheck()
{
#ifdef LICENSE_OFF
	return true;
#else
    std::string publicKey;
    std::vector<char> vecPublic;
    utils::SaveResToVector(L"PUBLIC_KEY", IDR_FILE_PUBLIC_KEY, &vecPublic);
    if (vecPublic.empty())
    {
		std::cout << "No license file\nPress Any key..." << std::endl;
        std::cin.get();
        return false;
    }
    publicKey.assign(vecPublic.begin(), vecPublic.end());
	std::string str =""; 

	CMacAddress macAddress;

	if(!macAddress.Load())
	{
		std::cout << "error: get MAC address\n";
		return false;
	}
	str = macAddress.Address();
	std::stringstream trimmer;
	trimmer << str;
	trimmer >> str;

	try
	{
		std::vector<char> vec;
		utils::LoadFileToVector(GetLicenseFilePath(), vec);
		if (utils::RsaVerifyVector(publicKey, str, vec))
		{
			std::cout << "License is valid. " << std::endl;
		}
		else
		{
			std::cout << "License is invalid. Program is being closed" << std::endl;
			return false;
		}
	}
	catch(const std::logic_error& ex)
	{
		std::cout << "You do not have license.dat file installed. Please put it in program dir." << std::endl;
		return false;
	}
	return true;
#endif
}