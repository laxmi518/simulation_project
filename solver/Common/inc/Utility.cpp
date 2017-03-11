#include "STdAfx.h"
#include <sys\types.h>
#include <sys\stat.h> 
#include "Utility.h"
#include <Gl\gl.h>
#include <cuda_runtime.h>

bool IsGLError()
{
	if(glGetError() != GL_NO_ERROR)
	{
		return false;
	}
	return true;
	}
bool ISCUDAError()
{
	cudaError status;
	status = cudaGetLastError();
	if(status != cudaSuccess)
	{
		std::string Error(cudaGetErrorString(status));
		return false;
	}
	return true;
}
bool ISCUDAError(cudaError_t Status)
{
	if(Status != cudaSuccess)
	{		
		return false;
	}
	return true;
}
//CCTStatusType CudaSafeCall(cudaError_t Status)
//{
//	if(cudaSuccess != Status)
//	{
//		printf(cudaGetErrorString(Status));
//		return CCT_CUDAERR;	
//	}
//	return CCT_NOERR;
//}
void Normalize(double& A, double& B, double& C)
{
	double Mag = sqrt(A * A + B * B + C * C);
	if(Mag > 0.0)
	{
		double Res = 1 / Mag;
		
		A *= Res;
		B *= Res;
		C *= Res;
	}
}

void Swap(float& a, float&b)
{
	float temp = a;
	a = b;
	b = temp;
}

__int64 SizeOfFile64( const char * szFileName ) 
{ 
  struct __stat64 fileStat; 
  int err = _stat64( szFileName, &fileStat ); 
  if (0 != err) return 0; 
  return fileStat.st_size; 
}
int SizeOfFile( const char * szFileName ) 
{ 
  struct stat fileStat; 
  int err = stat( szFileName, &fileStat ); 
  if (0 != err) return 0; 
  return fileStat.st_size; 
}