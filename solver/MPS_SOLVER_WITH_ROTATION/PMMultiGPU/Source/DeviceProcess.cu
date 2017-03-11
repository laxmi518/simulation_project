#include <stdio.h>
#include "..\\..\\..\\Common\\inc\\DataType.h"
#include "DeviceProcess.h"
#include "DeviceProcess_Kernel.h"

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
		Dg = iDivUp(NumThreads, Db);
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
	CCTStatusType RegisterParticleTopology(cudaStream_t &Stream,Integer ParticleNum)
	{
		CThreadScaler TS(ParticleNum);
		RegisterParticleTopology_Kernel<<<TS.Grids(), TS.Blocks() ,0,Stream>>>();
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType ResetParticleTopology(cudaStream_t &Stream,Integer CellNum, CCell* aCell, Integer* ParticleHash, Integer ParticleNum)
	{
		cudaMemset(ParticleHash, -1, ParticleNum * sizeof(Integer));
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		StatusType = ResetCellParticleTopology(Stream,CellNum, aCell);
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType ResetCellParticleTopology(cudaStream_t &Stream,Integer CellNum, CCell* aCell)
	{
		CThreadScaler TS(CellNum);
		ResetParticleTopology_Kernel<<<TS.Grids(),TS.Blocks() ,0,Stream>>>(CellNum, aCell);
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType CalcExplicitly(cudaStream_t &Stream,Integer ParticleNum,Scalar InnerAcceleration)
	{
		CThreadScaler TS(ParticleNum);
		CalcExplicitly_Kernel<<<TS.Grids(),TS.Blocks() ,0,Stream>>>(InnerAcceleration);
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType ResetTriangleTopology(cudaStream_t &Stream,Integer CellNum, CCell* aCell)
	{
		CThreadScaler TS(CellNum);
		ResetTriangleTopology_Kernel<<<TS.Grids(),TS.Blocks() ,0,Stream>>>(CellNum, aCell);
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType RegisterTriangleTopology(cudaStream_t &Stream,CTriangle* daTriangle, Integer TriangleNum, CCell* daCell, Integer CellNum)
	{
		unsigned int DbTriangle = min ( BLOCK_MAX_DIM, TriangleNum);
		unsigned int DgTriangle = iDivUp(TriangleNum, DbTriangle);
		RegisterTriangleTopology_Kernel<<<DgTriangle, DbTriangle ,0,Stream>>>(daTriangle, TriangleNum, daCell, CellNum);
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType UpdateTrianglePosition(cudaStream_t &Stream,const Integer TriangleNum, CTriangle* daTriangles)
	{
		unsigned int DbTriangle = min ( BLOCK_MAX_DIM, TriangleNum);
		unsigned int DgTriangle = iDivUp(TriangleNum, DbTriangle);
		UpdateTrianglePosition_Kernel<<<DgTriangle, DbTriangle, 0, Stream>>>(TriangleNum, daTriangles);
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType CalcTemperatureFactor(cudaStream_t &Stream,Integer ParticleNum)
	{
		CThreadScaler TS(ParticleNum);
		CalcTemperatureFactor_Kernel<<<TS.Grids(),TS.Blocks(), 0, Stream>>>();
		CCTStatusType StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType CalcSTLDistance(cudaStream_t &Stream,Integer ParticleNum)
	{
		CThreadScaler TS(ParticleNum);
		CalcSTLDistance_Kernel<<<TS.Grids(),TS.Blocks(), 0, Stream>>>();
		CCTStatusType StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	
	CCTStatusType RelocateParticleData(cudaStream_t &Stream,Integer ParticleNum)
	{
		CThreadScaler TS(ParticleNum);
		RelocateParticleData_Kernel<<<TS.Grids(),TS.Blocks(), 0, Stream>>>();
		CCTStatusType StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}	
	CCTStatusType CalcImplicitExplicitly(cudaStream_t &Stream,const Integer ParticleNum)
	{
		CThreadScaler TS(ParticleNum);
		CalcImplicitExplicitly_Kernel<<<TS.Grids(), TS.Blocks(), 0, Stream>>>();
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	CCTStatusType CalcPressureExplicit(cudaStream_t &Stream,const Integer ParticleNum)
	{
		CThreadScaler TS(ParticleNum);
		CalcPressureExplicit_Kernel<<<TS.Grids(),TS.Blocks(), 0, Stream>>>();
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}	
	CCTStatusType RotateTrianglePosition(cudaStream_t &Stream,const Integer TriangleNum, const Integer analysisStep)
	{
		unsigned int DbTriangle = min ( BLOCK_MAX_DIM, TriangleNum);
		unsigned int DgTriangle = iDivUp(TriangleNum, DbTriangle);
		RotateTrianglePosition_Kernel<<<DgTriangle, DbTriangle, 0, Stream>>>(TriangleNum, analysisStep);
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
		return CCT_NOERR;
	}
	//CCTStatusType ResetTrianglePosition(cudaStream_t &Stream,const Integer TriangleNum, CTriangle* deviceTriangle, const CTriangle* OriginalTriangle)			
	//{			
	//	unsigned int DbTriangle = min ( BLOCK_MAX_DIM, TriangleNum);		
	//	unsigned int DgTriangle = iDivUp(TriangleNum, DbTriangle);		
	//	ResetTrianglePosition_Kernel<<<DgTriangle, DbTriangle, 0, Stream>>>(TriangleNum, deviceTriangle ,OriginalTriangle);		
	//	CCTStatusType StatusType;		
	//	StatusType = CudaSafeCall(cudaGetLastError());		
	//	CCT_ERROR_CHECK(StatusType);		
	//	return CCT_NOERR;		
	//}
	CCTStatusType ResetWallPosition(cudaStream_t &Stream,const Integer TriangleNum,const Integer AnalysisStep,const CTriangle* daTriangles)
	{
		unsigned int DbTriangle = min ( BLOCK_MAX_DIM, TriangleNum);
		unsigned int DgTriangle = iDivUp(TriangleNum, DbTriangle);
		ResetWallPosition_Kernel<<<DgTriangle, DbTriangle,0,Stream>>>(TriangleNum,AnalysisStep, daTriangles);
		CCTStatusType StatusType;
		StatusType = CudaSafeCall(cudaGetLastError());
		CCT_ERROR_CHECK(StatusType);
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
		return CCT_NOERR;
	}
	CCTStatusType InitializeDeviceConstInputParticles(Integer * InputParticleID, Scalar3 * InputParticlePosition, Scalar3 * InputParticleVelocity,
		Scalar * InputParticlePressure, Scalar * InputParticleDensity, Scalar * InputParticleTemperature, Scalar * InputParticleKineticViscosity,
		Scalar * InputParticleSolidPhaseRate, ParticleType * InputParticleType)
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
		return CCT_NOERR;
	}
	CCTStatusType InitializeDeviceMemConst(CParameter Parameter,Integer ParticleNum,CTriangle * Triangles,Integer TriangleNum,CTriangleParameters * TriangleParameters,
		Integer MaxParticleNum, CDistance * STLDistance,Integer * StlID, CCell * Cell,Integer CellNum,
		Integer BucketNum,/*Integer * GridParticleIndex, CGridParams GridParams,*/
		CGridBox BoundingBox,Integer * ParticleHash)
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
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_BucketNum, &BucketNum, sizeof(Integer)));
		CCT_ERROR_CHECK(Status);
		/*
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_cellStart, &CellStart, sizeof(CellStart)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_cellEnd, &CellEnd, sizeof(CellEnd)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_gridParticleIndex, &GridParticleIndex, sizeof(GridParticleIndex)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(CONSTANT_GRIDPARAMS, &GridParams, sizeof(CGridParams)));
		CCT_ERROR_CHECK(Status);
		*/
		Status = CudaSafeCall( cudaMemcpyToSymbol(CONSTANT_BOUNDINGBOX, &BoundingBox, sizeof(CGridBox)) );
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_dParticleHash, &ParticleHash, sizeof(ParticleHash)));
		CCT_ERROR_CHECK(Status);

		return CCT_NOERR;
	}
	CCTStatusType ParticleNumberToConst(Integer ParticleNum)
	{
		CCTStatusType Status;
		Status = CudaSafeCall(cudaMemcpyToSymbol(c_ParticleNum, &ParticleNum, sizeof(Integer)));
		CCT_ERROR_CHECK(Status);
		return CCT_NOERR;
	}
	//Drag@Rajan
	CCTStatusType DragParametersToConst(DragParameter *InputDragParameter,Scalar3 *InputDragAcc,Scalar* InputDragTemperature ,Integer DragTriangleNum)
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
		return CCT_NOERR;
	}

	CCTStatusType CalcDragEffect(cudaStream_t &Stream, Integer ComputeParticleNum)
	{
		if(ComputeParticleNum > 0)
		{
			CThreadScaler TS(ComputeParticleNum);

			CalcDragEffect_Kernel<<<TS.Grids(),TS.Blocks(),0,Stream>>>();
			
			CCTStatusType StatusType = CudaSafeCall(cudaGetLastError());
			CCT_ERROR_CHECK(StatusType);
		}
		return CCT_NOERR;
	}
}
