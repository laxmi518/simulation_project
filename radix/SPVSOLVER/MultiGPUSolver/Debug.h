#pragma once
#include "std.h"

namespace Debug
{
	template<typename T>
	bool Dump(const std::string& Filename, size_t Num, T* DeviceArr)
	{
		if(0 == Num || NULL == DeviceArr)
		{
			return false;
		}
		std::ofstream Out(Filename.c_str());
		if(!Out)
		{
			return false;
		}

		size_t Size = Num * sizeof(T);
		T* HostArr = (T*)malloc(Size);
		cudaMemcpy(HostArr, DeviceArr, Size, cudaMemcpyDeviceToHost);
		for(size_t i = 0; i < Num; ++i)
		{
			Out << i << "," << HostArr[i] << "\n";
		}
		Out.close();
		free(HostArr);
		return true;
	}
	template<typename T>
	bool DumpSTLDistMagintude(const std::string& FileName, size_t Num, T* DeviceArr)
	{
		if(0==Num || NULL == DeviceArr)
		{
			return false;
		}
		std::ofstream Out(FileName.c_str());
		if(!Out)
		{
			return false;
		}
		size_t Size = Num * sizeof(T);
		T* HostArr = (T*)malloc(Size);
		cudaMemcpy(HostArr,DeviceArr,Size, cudaMemcpyDeviceToHost);
		for(size_t i = 0 ; i < Num ; ++i)
		{
			Out << i << "," << HostArr[i].Magnitude <<"\n";
		}
		Out.close();
		free(HostArr);
		return true;
	}
	//bool Dump(const std::string& Filename, size_t Num)
	//{
	//	std::ofstream Out(Filename.c_str());
	//	if(!Out)
	//	{
	//		return false;
	//	}

	//	size_t Size = Num * sizeof(T);
	//	T* HostArr = (T*)malloc(Size);
	//	cudaMemcpy(HostArr, DeviceArr, Size, cudaMemcpyDeviceToHost);
	//	for(size_t i = 0; i < Num; ++i)
	//	{
	//		Out << i << "," << HostArr[i] << "\n";
	//	}
	//	Out.close();
	//	free(HostArr);
	//	return true;
	//}
}
