//#ifdef _DEBUG
//#include <vld.h>
//#endif
#include "stdafx.h"
#include <cuda_runtime_api.h>
#include "Licence.h"
#include "MultiGPUHandler.h"

int main(int argc, char** argv)
{
	/*Validate Licence Starts*/
	CLicenceValidate LicenceCheck;
	bool Validate = false;
	Validate = LicenceCheck.LicenceCheck();
	if(!Validate)
		return -1;
	/*Validate Licence Ends*/


	int deviceCount = 0;
	if (cudaGetDeviceCount(&deviceCount) != cudaSuccess) {
		printf("cudaGetDeviceCount FAILED CUDA Driver and Runtime version may be mismatched.\n");
		printf("\nFAILED\n");
		return -1;
	}
	
	// This function call returns 0 if there are no CUDA capable devices.
    if (deviceCount == 0)
	{
		printf("There is no device supporting CUDA\n");
		return -1;
	}
	if(deviceCount == 1)
	{
		int dev = 0;
		cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
		// This function call returns 9999 for both major & minor fields, if no CUDA capable devices are present
		if (deviceProp.major == 9999 && deviceProp.minor == 9999)
		{
                printf("There is no device supporting CUDA.\n");
				return -1;
		}
	}  
	printf("%d Devices Found\n",deviceCount);

	CMultiGPUHandler* pCalculationManager = NULL;
	cublasInit();
	
	pCalculationManager = new CMultiGPUHandler();	
	
	CCTStatusType StatusType = pCalculationManager->Run(argc, argv);
	if(CCT_NOERR != StatusType)
	{
		delete pCalculationManager;
		cublasShutdown();
		return -1;
	}
	cublasShutdown();
	delete pCalculationManager;

	return 0;
}