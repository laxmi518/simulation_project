#include <stdio.h>
#include "DataType.h"
#include "DeviceProcess.h"
#include "DeviceProcess_Kernel.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/random.h>
#include <thrust/generate.h>
#include <thrust/detail/type_traits.h>
#include "IO.h"
#include "Configuration.h"

#define USECDUASTREAM

//Round a / b to nearest higher integer value	
unsigned int iDivUp(unsigned int a, unsigned int b){
	return (a % b != 0) ? (a / b + 1) : (a / b);	
}
class CThreadScaler
{
private:
	Integer Dg;
	Integer Db;
public:
	CThreadScaler(Integer NumThreads)
	{
		Db = min ( BLOCK_MAX_DIM, NumThreads);
		if(Db > 0)
		{
			Dg = iDivUp(NumThreads, Db);
		}else
		{
			Dg = 0;
		}
	}
	Integer Grids()
	{
		return Dg;
	}
	Integer Blocks()
	{
		return Db;
	}
};
extern "C"
{
	/*inline void check_cuda_errors(const char *filename, const int line_number)
	{
		#ifdef DEBUG
		  cudaThreadSynchronize();
		  cudaError_t error = cudaGetLastError();
		  if(error != cudaSuccess)
		  {
			printf("CUDA error at %s:%i: %s\n", filename, line_number, cudaGetErrorString(error));
			exit(-1);
		  }
		#endif
	}*/
	CCTStatusType CalcTurbulenceViscosity(cudaStream_t &Stream,Integer ComputeParticleNum)
	{
		if(ComputeParticleNum > 0)
		{
			CCTStatusType StatusType;
			CThreadScaler TS(ComputeParticleNum);
			CalcTurbulenceViscosity_Kernel<<<TS.Grids(), TS.Blocks() ,0,Stream>>>(ComputeParticleNum);
			//std::string kernelName = "CalcTurbulenceViscosity";
			//WriteConstant(kernelName,ComputeParticleNum);
			//check_cuda_errors(__FILE__, __LINE__);
			StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}
	//Check Particle outside compute zone starts
	CCTStatusType CheckParticleOutsideComputeZone(cudaStream_t &Stream,Integer *ParticleNum)
	{
		if((*ParticleNum) > 0)
		{
			CThreadScaler TS(*ParticleNum);
			CheckParticleOutsideComputeZone_Kernel<<<TS.Grids(),TS.Blocks(),0,Stream>>>(ParticleNum);
			//check_cuda_errors(__FILE__, __LINE__);
			CCTStatusType StatusType;
			StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}
	//Check particle Outside Compute Zone Ends

	CCTStatusType CalcExplicitly(cudaStream_t &Stream,Integer ComputeParticleNum)
	{
		if((ComputeParticleNum) > 0)
		{
			CCTStatusType StatusType;
			CThreadScaler TS(ComputeParticleNum);
			CalcExplicitly_Kernel<<<TS.Grids(),TS.Blocks(),0,Stream>>>(ComputeParticleNum);
			//check_cuda_errors(__FILE__, __LINE__);
			StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}
	CCTStatusType ResetTriangleTopology(cudaStream_t &Stream,Integer CellNum, CCell* aCell)
	{
		if(CellNum > 0)
		{
			CThreadScaler TS(CellNum);
			ResetTriangleTopology_Kernel<<<TS.Grids(),TS.Blocks(),0,Stream>>>(CellNum, aCell);
			//check_cuda_errors(__FILE__, __LINE__);
			CCTStatusType StatusType;
			StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}
	CCTStatusType RegisterTriangleTopology(cudaStream_t &Stream,CTriangle* daTriangle, Integer TriangleNum, CCell* daCell, Integer CellNum)
	{
		unsigned int DbTriangle = min ( BLOCK_MAX_DIM, TriangleNum);
		unsigned int DgTriangle = iDivUp(TriangleNum, DbTriangle);
		RegisterTriangleTopology_Kernel<<<DgTriangle, DbTriangle,0,Stream>>>(daTriangle, TriangleNum, daCell, CellNum);
		//check_cuda_errors(__FILE__, __LINE__);
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType UpdateTrianglePosition(cudaStream_t &Stream,const Integer TriangleNum, CTriangle* daTriangles)
	{
		unsigned int DbTriangle = min ( BLOCK_MAX_DIM, TriangleNum);
		unsigned int DgTriangle = iDivUp(TriangleNum, DbTriangle);
		UpdateTrianglePosition_Kernel<<<DgTriangle, DbTriangle,0,Stream>>>(TriangleNum, daTriangles);
		//check_cuda_errors(__FILE__, __LINE__);
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType RotateTrianglePosition(cudaStream_t &Stream,const Integer TriangleNum,CTriangle* daTriangles, const Integer analysisStep)
	{
		unsigned int DbTriangle = min ( BLOCK_MAX_DIM, TriangleNum);
		unsigned int DgTriangle = iDivUp(TriangleNum, DbTriangle);
		RotateTrianglePosition_Kernel<<<DgTriangle, DbTriangle,0,Stream>>>(TriangleNum, daTriangles, analysisStep);
		//check_cuda_errors(__FILE__, __LINE__);
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType ResetWallPosition(cudaStream_t &Stream,const Integer TriangleNum,const Integer AnalysisStep,const CTriangle* daTriangles)
	{
		unsigned int DbTriangle = min ( BLOCK_MAX_DIM, TriangleNum);
		unsigned int DgTriangle = iDivUp(TriangleNum, DbTriangle);
		ResetWallPosition_Kernel<<<DgTriangle, DbTriangle,0,Stream>>>(TriangleNum,AnalysisStep, daTriangles);
		//check_cuda_errors(__FILE__, __LINE__);
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}	

	CCTStatusType CalcSTLDistance(cudaStream_t &Stream, Integer ComputeParticleNum)
	{
		if(ComputeParticleNum > 0)
		{
			CThreadScaler TS(ComputeParticleNum);
			CalcSTLDistance_Kernel<<<TS.Grids(),TS.Blocks(),0,Stream>>>(ComputeParticleNum);
			//check_cuda_errors(__FILE__, __LINE__);
			CCTStatusType StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}
	CCTStatusType CaculateCellIDandInitializeHash(cudaStream_t &Stream,Integer ParticleNum,Integer CellNum,int* dGridParticleHash, int* dGridParticleIndex, Scalar3* particlePosition)
	{
		if(ParticleNum > 0)
		{
			CThreadScaler TS(ParticleNum);
			// calculate grid hash
			calcHashD<<<TS.Grids(), TS.Blocks() ,0,Stream>>>(ParticleNum,dGridParticleHash,dGridParticleIndex,particlePosition);
			//check_cuda_errors(__FILE__, __LINE__);
			CCTStatusType StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}
	
	CCTStatusType reorderDataAndFindCellStart(cudaStream_t& Stream,Integer numParticles,Integer numCells , int*  gridParticleHash, int*  gridParticleIndex, int*  cellStart, int*  cellEnd)
	{
		if(numParticles > 0)
		{
			CThreadScaler TS(numParticles);
			// set all cells to empty
			cudaMemsetAsync(cellStart, 0xffffffff, numCells*sizeof(int),Stream);
			cudaMemsetAsync(cellEnd,   0xffffffff, numCells*sizeof(int),Stream);
			int smemSize = sizeof(int)*(TS.Blocks()+1);
			reorderDataAndFindCellStartD<<< TS.Grids(), TS.Blocks(), smemSize,Stream>>>(numParticles,numCells,gridParticleHash,	gridParticleIndex,cellStart, cellEnd);
			//check_cuda_errors(__FILE__, __LINE__);
			CCTStatusType StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}
	CCTStatusType SortUsingThrust(Integer MaxParticleNum, Integer * daNumberHash, Integer* daNumberIndex)
	{
		CCTStatusType  Status = CCT_NOERR;
		if(MaxParticleNum > 0)
		{
			thrust::sort_by_key(thrust::device_ptr<Integer>(daNumberHash),
								thrust::device_ptr<Integer>(daNumberHash + MaxParticleNum),
								thrust::device_ptr<Integer>(daNumberIndex));
			//check_cuda_errors(__FILE__, __LINE__);
			Status = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(Status);
		}
		return Status;
	}
	CCTStatusType StableSortUsingThrust(Integer MaxParticleNum, Integer * daNumberHash, Integer* daNumberIndex)
	{
		CCTStatusType  Status = CCT_NOERR;
		if(MaxParticleNum > 0)
		{
			thrust::stable_sort_by_key(thrust::device_ptr<Integer>(daNumberHash),
									   thrust::device_ptr<Integer>(daNumberHash + MaxParticleNum),
									   thrust::device_ptr<Integer>(daNumberIndex));
			//check_cuda_errors(__FILE__, __LINE__);
			Status = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(Status);
		}
		return Status;
	}
	CCTStatusType CalcDragEffect(cudaStream_t &Stream,Integer ComputeParticleNum)
	{
		if((ComputeParticleNum) > 0)
		{
			CThreadScaler TS(ComputeParticleNum);
			CalcDragEffect_Kernel<<<TS.Grids(),TS.Blocks(),0,Stream>>>(ComputeParticleNum);
			//check_cuda_errors(__FILE__, __LINE__);
			CCTStatusType StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}
	CCTStatusType CalcExplicitPressure(cudaStream_t &Stream,Integer ComputeParticleNum)
	{
		if((ComputeParticleNum) > 0)
		{
			CThreadScaler TS(ComputeParticleNum);
			CalcExplicitPressure_Kernel<<<TS.Grids(),TS.Blocks(),0,Stream>>>(ComputeParticleNum);
			//check_cuda_errors(__FILE__, __LINE__);
			CCTStatusType StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}
	CCTStatusType CalcExplicitPressureGradient(cudaStream_t &Stream,Integer ComputeParticleNum)
	{
		if((ComputeParticleNum) > 0)
		{
			CThreadScaler TS(ComputeParticleNum);
			CalcExplicitPressureGradient_Kernel<<<TS.Grids(),TS.Blocks(),0,Stream>>>(ComputeParticleNum);
			//check_cuda_errors(__FILE__, __LINE__);
			CCTStatusType StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}
	CCTStatusType CalcTemperatureFactor(cudaStream_t &Stream,Integer ComputeParticleNum)
	{
		if((ComputeParticleNum) > 0)
		{
			CThreadScaler TS(ComputeParticleNum);
			CalcTemperatureFactor_Kernel<<<TS.Grids(),TS.Blocks(),0,Stream>>>(ComputeParticleNum);
			//check_cuda_errors(__FILE__, __LINE__);
			CCTStatusType StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}

	CCTStatusType InitializeDeviceMemConst(CParameter Parameter,Integer ParticleNum,CTriangle * Triangles,Integer TriangleNum,CTriangleParameters * TriangleParameters,
		Integer MaxParticleNum, CDistance * STLDistance,Integer * StlID, CCell * Cell,Integer CellNum,
		Integer * CellStart,Integer * CellEnd,Integer * GridParticleIndex, CGridBox BoundingBox)
	{
		CCTStatusType Status;
		Status = CudaSafeCall(cudaMemcpyToSymbol(CONSTANT_PARAMETER, &Parameter, sizeof(CParameter)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_ParticleNum, &ParticleNum, sizeof(Integer)));
		CCT_ERROR_CHECK(Status);		
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daTriangles, &Triangles, sizeof(Triangles)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_TriangleNum, &TriangleNum, sizeof(Integer)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daTrianglesParameters, &TriangleParameters, sizeof(TriangleParameters)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_MaxParticleNum, &MaxParticleNum, sizeof(Integer)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daSTLDistance, &STLDistance, sizeof(STLDistance)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daSTLID, &StlID, sizeof(STLID)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daCell, &Cell, sizeof(Cell)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_CellNum, &CellNum, sizeof(Integer)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_cellStart, &CellStart, sizeof(CellStart)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_cellEnd, &CellEnd, sizeof(CellEnd)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_gridParticleIndex, &GridParticleIndex, sizeof(GridParticleIndex)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall( cudaMemcpyToSymbol(CONSTANT_BOUNDINGBOX, &BoundingBox, sizeof(CGridBox)) );
		CCT_ERROR_CHECK(Status);
		//check_cuda_errors(__FILE__, __LINE__);
		return CCT_NOERR;
	}
	CCTStatusType InitializeDeviceConstOutPutParticles(Integer * OutputParticleID, Scalar3 * OutputParticlePosition, Scalar3 * OutputParticleVelocity,
		Scalar * OutputParticlePressure, Scalar * OutputParticleDensity, Scalar * OutputParticleTemperature, Scalar * OutputParticleKineticViscosity,
		Scalar * OutputParticleSolidPhaseRate, ParticleType * OutputParticleType)
	{
		CCTStatusType Status;
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daOutputParticleID,&OutputParticleID, sizeof(OutputParticleID)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daOutputParticlePosition,&OutputParticlePosition, sizeof(OutputParticlePosition)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daOutputParticleVelocity,&OutputParticleVelocity, sizeof(OutputParticleVelocity)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daOutputParticlePressure,&OutputParticlePressure, sizeof(OutputParticlePressure)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daOutputParticleDensity,&OutputParticleDensity, sizeof(OutputParticleDensity)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daOutputParticleTemperature,&OutputParticleTemperature, sizeof(OutputParticleTemperature)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daOutputParticleKineticViscosity,&OutputParticleKineticViscosity, sizeof(OutputParticleKineticViscosity)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daOutputParticleSolidPhaseRate,&OutputParticleSolidPhaseRate, sizeof(OutputParticleSolidPhaseRate)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daOutputParticleType,&OutputParticleType, sizeof(OutputParticleType)));
		CCT_ERROR_CHECK(Status);
		//check_cuda_errors(__FILE__, __LINE__);
		return CCT_NOERR;
	}
	CCTStatusType InitializeDeviceConstInputParticles(Integer * InputParticleID, Scalar3 * InputParticlePosition, Scalar3 * InputParticleVelocity,
		Scalar * InputParticlePressure, Scalar * InputParticleDensity, Scalar * InputParticleTemperature, Scalar * InputParticleKineticViscosity,
		Scalar * InputParticleSolidPhaseRate, ParticleType * InputParticleType, Scalar* ParticleTurbulaceViscosity, Scalar* ParticleStrainTensorProduct)
	{
		CCTStatusType Status;
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daParticleID,&InputParticleID, sizeof(InputParticleID)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daParticlePosition,&InputParticlePosition, sizeof(InputParticlePosition)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daParticleVelocity,&InputParticleVelocity, sizeof(InputParticleVelocity)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daParticlePressure,&InputParticlePressure, sizeof(InputParticlePressure)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daParticleDensity,&InputParticleDensity, sizeof(InputParticleDensity)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daParticleTemperature,&InputParticleTemperature, sizeof(InputParticleTemperature)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daParticleKineticViscosity,&InputParticleKineticViscosity, sizeof(InputParticleKineticViscosity)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daParticleSolidPhaseRate,&InputParticleSolidPhaseRate, sizeof(InputParticleSolidPhaseRate)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daParticleType,&InputParticleType, sizeof(InputParticleType)));
		CCT_ERROR_CHECK(Status);
		//Turbulace 
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daParticleTurbulaceViscosity,&ParticleTurbulaceViscosity, sizeof(ParticleTurbulaceViscosity)));		
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daParticleStrainTensorProduct,&ParticleStrainTensorProduct, sizeof(ParticleStrainTensorProduct)));		
		CCT_ERROR_CHECK(Status);
		//check_cuda_errors(__FILE__, __LINE__);
		return CCT_NOERR;
	}
	CCTStatusType ParticleNumberToConst(Integer ParticleNum)
	{
		CCTStatusType Status;
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_ParticleNum, &ParticleNum, sizeof(Integer)));
		//check_cuda_errors(__FILE__, __LINE__);
		CCT_ERROR_CHECK(Status);
		return CCT_NOERR;
	}
	CCTStatusType DragParametersToConst(DragParameter *InputDragParameter,Scalar3 *InputDragAcc,Scalar* InputDragTemperature ,Integer DragTriangleNum,Integer * MagnifierCount,CDragTriangle * DragTriangles)
	{
		CCTStatusType Status;
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daSTLDragParameter,&InputDragParameter, sizeof(InputDragParameter)));
		CCT_ERROR_CHECK(Status);

		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daDragAcc,&InputDragAcc, sizeof(InputDragAcc)));
		CCT_ERROR_CHECK(Status);

		//For Drag Temperature
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daDragTemperature,&InputDragTemperature, sizeof(InputDragTemperature)));
		CCT_ERROR_CHECK(Status);

		Status = CudaSafeCall(cudaMemcpyToSymbol(c_DragTriangleNum, &DragTriangleNum, sizeof(DragTriangleNum)));
		CCT_ERROR_CHECK(Status);
		
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daMagnifierCount, &MagnifierCount, sizeof(MagnifierCount)));
		CCT_ERROR_CHECK(Status);

		Status = CudaSafeCall(cudaMemcpyToSymbol(c_daDragTriangles, &DragTriangles, sizeof(DragTriangles)));
		CCT_ERROR_CHECK(Status);
		//check_cuda_errors(__FILE__, __LINE__);
		return CCT_NOERR;
	}

}