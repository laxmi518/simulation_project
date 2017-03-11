// GetMACUuid.cpp : Defines the entry point for the console application.
//
// Author:	Khalid Shaikh [Shake@ShakeNet.com]
// Date:	April 5th, 2002
//
// This program fetches the MAC address of the localhost by creating a UUID
// and obtaining the IP address through that
//
// Prewindows 2000 one should replace the function below 
//			UuidCreateSequential with UuidCreate
// Microsoft decided that it was not wise to put the MAC address in the UUID
// hence changed the algorithm in Windows 2000 and later to not include it.
//
// Supported in Windows NT/2000/XP
// Supported in Windows 95/98/Me
//
// Supports single NIC card.
/*
#include <iostream>
*/
#include <stdio.h>
#include <Windows.h>
#include <rpc.h>
#include <rpcdce.h>
#pragma comment(lib, "rpcrt4.lib")

namespace MACADDRESSGETTER
{
	// Prints the MAC address stored in a 6 byte array to stdout
	void PrintMACaddress(unsigned char MACData[])
	{
		printf("MAC Address: %02X-%02X-%02X-%02X-%02X-%02X\n", 
			MACData[0], MACData[1], MACData[2], MACData[3], MACData[4], MACData[5]);
	}
	// Fetches the MAC address and prints it
	void GetMACaddress(void)
	{
		unsigned char MACData[6];

		UUID uuid;
		UuidCreateSequential( &uuid );				// Ask OS to create UUID

		for (int i=2; i<8; i++)						// Bytes 2 through 7 inclusive are MAC address
			MACData[i - 2] = uuid.Data4[i];

		PrintMACaddress(MACData);					// Print MAC address
	}
}

/*
int _tmain(int argc, _TCHAR* argv[])
{
	GetMACaddress();							// Obtain MAC address of adapters

	return 0;
}
*/
