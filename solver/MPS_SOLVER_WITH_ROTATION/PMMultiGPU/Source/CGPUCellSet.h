#pragma once
#define MAXNEIGHBORGRIDS 26
#include "InnerPressureModel.h"
class CMultiGPUHandler;

#define INDIVIDUALTYPE

class CCellSet
{
	
public:
	CCellSet();
	~CCellSet();
	CParameter m_Parameter;
	CGridBox m_BoundingBox;
	
	CCell* m_haCell;
	CCell* m_daCell;
	Integer m_CellNum;
	Integer m_ComputeCellNum;
	Integer m_GhostCellNum;

	CParticle* m_haParticle; //Host Particle Array
	CParticle* m_daParticle; //Devide particle Array
	//CParticle* m_daOutputParticle;
	Integer m_ParticleNum; //Number of current Particles
	Integer m_MaxParticleNum;//Maximum allocated size

	CParticle* m_haBucket; // Bucket particle to fill later if necessary
	CParticle* m_daBucket; // Device particle to fill later if necessary
	Integer m_BucketNum; // Number of Bucket Particles
	

	CTriangle * m_haTriangles;// Cell and Grid of the STLTriangles
	CTriangle * m_daTriangles;// Cell and Grid of the STLTriangles
	CTriangle * m_daOriginalTriangles;// Storing Original Cell and Grid of the STLTriangles
	Integer m_TriangleNum;

	CTriangleParameters * m_haTriangleParameters;

	//Drag@Rajan
	DragParameter	*m_haSTLDragPrameter;
	DragParameter	*m_daSTLDragPrameter;
	Integer m_DragTriangleNum;


	CTriangleParameters * m_daTriangleParameters;
	Integer m_TriangleModelNumber;

	CDragTriangle * m_haDragTriangles;// Cell and Grid of the STLTriangles
	CDragTriangle * m_daDragTriangles;// Cell and Grid of the STLTriangles
	
	Scalar3* m_daDragAcc;
	//Scalar3* m_haDragAcc;
	Scalar* m_daDragTemperature;
	//Scalar* m_haDragTemperature;

	Scalar3* m_daSubgridScaleStress;
	
	
	CDistance *m_daSTLDistance;// Holds the distance between particle and triangle
	Integer* m_daSTLID;// Holds the ID of the neighbor triangle

	

	Integer* m_dParticleHash;// Stores Next ID's of the particle in the cell
	Integer* m_hParticleHash;

	Scalar* m_dX;
	Scalar* m_hX;

	Scalar* m_dX0;
	Scalar* m_hX0;

	Scalar* m_dMInv;
	Scalar* m_hMInv;

	Scalar* m_dAii;
	Scalar* m_hAii;

	Scalar* m_dB;
	Scalar* m_hB;

	Scalar* m_dEps;
	Scalar* m_hEps;

	Scalar* m_dReduceData;
	Scalar* m_hReduceData;

	Scalar m_ParticleMass;


	CMultiGPUHandler* m_pOutputHandler;
	CInnerPressureModel* m_pInnerPressureModel;

	Integer m_OldInNumber;
	Integer m_NewInNumber;
	Integer m_OutPutNum;

	cudaStream_t m_CudaStream; 

#ifdef INDIVIDUALTYPE
	//@Rajan Individual Components
	Integer			*m_haParticleID;
	Scalar3			*m_haParticlePosition;
	Scalar3			*m_haParticleVelocity;
	Scalar			*m_haParticlePressure;
	Scalar			*m_haParticleDensity;
	Scalar			*m_haParticleTemperature;
	Scalar			*m_haParticleKineticViscosity;
	Scalar			*m_haParticleSolidPhaseRate;
	ParticleType	*m_haParticleType;

	Integer			*m_daParticleID;
	Scalar3			*m_daParticlePosition;
	Scalar3			*m_daParticleVelocity;
	Scalar			*m_daParticlePressure;
	Scalar			*m_daParticleDensity;
	Scalar			*m_daParticleTemperature;
	Scalar			*m_daParticleKineticViscosity;
	Scalar			*m_daParticleSolidPhaseRate;
	ParticleType	*m_daParticleType;

	
	Integer			*m_daOutputParticleID;
	Scalar3			*m_daOutputParticlePosition;
	Scalar3			*m_daOutputParticleVelocity;
	Scalar			*m_daOutputParticlePressure;
	Scalar			*m_daOutputParticleDensity;
	Scalar			*m_daOutputParticleTemperature;
	Scalar			*m_daOutputParticleKineticViscosity;
	Scalar			*m_daOutputParticleSolidPhaseRate;  //14 22 29 33 36 40 38 
	ParticleType	*m_daOutputParticleType;

	Integer			*m_haBucketID;
	Scalar3			*m_haBucketPosition;
	Scalar3			*m_haBucketVelocity;
	Scalar			*m_haBucketPressure;
	Scalar			*m_haBucketDensity;
	Scalar			*m_haBucketTemperature;
	Scalar			*m_haBucketKineticViscosity;
	Scalar			*m_haBucketSolidPhaseRate;
	ParticleType	*m_haBucketType;
/*
	Integer			*m_daBucketID;
	Scalar3			*m_daBucketPosition;
	Scalar3			*m_daBucketVelocity;
	Scalar			*m_daBucketPressure;
	Scalar			*m_daBucketDensity;
	Scalar			*m_daBucketTemperature;
	Scalar			*m_daBucketKineticViscosity;
	Scalar			*m_daBucketSolidPhaseRate;
	ParticleType	*m_daBucketType;
*/
#endif
	
private:
	// バウンディングボックスのオフセット
	void Offset(const Scalar Offset);
	Integer3 m_GridID;
	CCellSet* m_aNeighborGrid[MAXNEIGHBORGRIDS];
	Integer m_NeighborGridNum;
public:
	void SetGridID(Integer x, Integer y, Integer z);
	void SetNeighborGrids(CCellSet* aGrid, Integer3 Section);
	void SetComputeZone(const CBox* const ComputeZone, CParameter* pParameter, const CBox* const ParentZone);
	void SetStartStep(const Integer StartStep);
	CCTStatusType InitializeHostMemory(Integer MaxParticleNum);
	void UninitializeHostMemory();
	//bool RegisterParticleTopology(const CParticle * const Particle);
	//bool RegisterTriangleTopology(const CTriangle * const Triangle, const Integer TID);
	void SetTriangles(CTriangle* Triangles, Integer TriangleNum);
	void SetDragTriangles(CDragTriangle* DragTriangles, Integer DragTriangleNum,CInnerPressureModel* InnerPressureModel);
	CParticle* GetGhostBuffer();	
	CCTStatusType SetParticleData(CParticle* haParticle, Integer ParticleNum);
	CCTStatusType SetBucketData(CParticle* aBucket, Integer BucketNum);

public:
	CCTStatusType InitializeDeviceMemory();
	CCTStatusType UnInitializeDeviceMemory();
	CCTStatusType InitializeDeviceConstantMemory();
	CCTStatusType ParticleDeviceToHost();
	CCTStatusType RelocateAllParticleTopology();
	CCTStatusType RegisterParticleTopology();
	CCTStatusType RegisterTriangleTopology();
	CCTStatusType RegisterDragTriangleTopology();
	CCTStatusType CalculateSTLDistance();
	CCTStatusType SetComputeType(ComputeType Type);

private:
	// Thread Management
	//! handle to the thread threadProc
	HANDLE m_ThreadHandle;
	//! handle to the mutex
	HANDLE m_MutexHandle;// I think that now it is gone unnecessary but it is used in ServerHandler constructor
	// Thread id
	DWORD m_ThreadID;
	HANDLE m_Output; 
	HANDLE m_Ghost;
	HANDLE m_Buffer;
	HANDLE m_OutputBuffer;
	HANDLE m_AddParticles;
	CRITICAL_SECTION m_CriticalSection;
public:
	CCTStatusType WaitForOutput();
	CCTStatusType WaitForGhost();

	CCTStatusType SetBuffer();
	CCTStatusType ResetBuffer();
	CCTStatusType WaitForBuffer();
	CCTStatusType SetOutputBuffer();
	CCTStatusType WaitForOutputBuffer();
	CCTStatusType SetAddedParticlesEvent();

private:
	void run(void);
	
	CCTStatusType CalculatePressure(Integer ParticleNum, Scalar* daEps,Scalar* dReduceData,Scalar* haEps, Integer LoopMax, Scalar ConvergenceConstant);
	CCTStatusType CudaSafeCall(cudaError_t Status);
	CCTStatusType Calculate();
	CCTStatusType CalculateTwoPhase();
	CCTStatusType Compute(Integer i);

	CCTStatusType CalculatePressureCG();
	CCTStatusType CalculatePressureNJPCG();

	CCTStatusType CheckAndUpdateDt(Integer CurrentStep);

	Integer ParticlesBeyondComputeZone();
	Integer ParticlesInsideComputeZone();
	Scalar HeightFromDensity();

public:
	void WaitToStop(void);
	CCTStatusType ResetOutputEvent();
	CCTStatusType ResetGhostEvent();
	bool IsEnd();
	void RegisterTestParticleTopology_Kernel();
	CCTStatusType MInverseR();
	CCTStatusType SolveVelocity();

	CCTStatusType SetOutputEvent();

private:
	Integer m_DeviceID;
	bool m_End;
	bool m_Err;
	Integer m_StartStep;

	Integer m_NextItr;

public:
	//Newly Added
	CCTStatusType CalcSTLDistance();
	CCTStatusType CalcExplicitly();
	CCTStatusType CalcPressureExplicit();
	CCTStatusType CalcPressureExplicitGradient();
	CCTStatusType UpdateTrianglePosition();
	CCTStatusType RotateTrianglePosition(Integer i);
	CCTStatusType AddToOutputBuffer(CParticle* ParticleBuffer, Integer & ParticleNum);
	CCTStatusType AddToOutputBuffer(bool Output,Integer BufferSize, Integer* ParticleBufferID, Scalar3* ParticleBufferPosition,
									Scalar3* ParticleBufferVelcity, Scalar* ParticleBufferPressure,
									Scalar* ParticleBufferDensity, Scalar* ParticleBufferTemperature, Scalar* ParticleBufferKineticViscosity,
									Scalar* ParticleBufferSolidPhaseRate, ParticleType* ParticleBufferType);

	CCTStatusType SetParticleData(Integer BufferSize,
									  Integer* ParticleBufferID, Scalar3* ParticleBufferPosition,
									  Scalar3* ParticleBufferVelcity, Scalar* ParticleBufferPressure,
									  Scalar* ParticleBufferDensity, Scalar* ParticleBufferTemperature, Scalar* ParticleBufferKineticViscosity,
									  Scalar* ParticleBufferSolidPhaseRate, ParticleType* ParticleBufferType);

	CCTStatusType ResetWallPosition(Integer AnalysisStep);
	CCTStatusType ParticleHostToDevice();
	CCTStatusType ParticleNumberToContantMemory();
	CCTStatusType WaitForAddedParticles();

	CCTStatusType CalcTemperatureFactor();

	CCTStatusType InitializeGPU(Integer deviceID);

	CCTStatusType StreamSynchronize();

	void SetTrianglesParameters(CTriangleParameters * TriangleParameter, Integer TriangleModelNumber);

	void SetDragParameters(DragParameter * dragParameter, Integer dragTriangleNumber);

	//friend class CMultiGPUHandler;

	//Drag@Rajan
	CCTStatusType CalculateDragEffect();
};
#define VALIDATE_CELL(CellID, CellNum) if(CellID < 0 || CellID >= CellNum) return CCT_CALCERR ;