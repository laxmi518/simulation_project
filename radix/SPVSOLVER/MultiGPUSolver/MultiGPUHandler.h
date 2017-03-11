#pragma once
#include "DeviceSelector.h"
class CCellSet;
class CConfiguration;
class COutputFrame;
class DtFrame;
class CMultiGPUHandler 
{
public:
	CMultiGPUHandler(void);
public:
	~CMultiGPUHandler(void);
private:
	CGridBox m_CalculationBox;
	CBox m_ComputeZone;
	Integer m_GridNum;
	CCellSet* m_aGrid;
	Scalar* m_aResults;

	//STARTS for individual Output particles----------------------E&T Nepal August 2011--------------------------
	Integer			*m_aOutputParticleID;
	Scalar3			*m_aOutputParticlePosition;
	Scalar3			*m_aOutputParticleVelocity;
	Scalar			*m_aOutputParticlePressure;
	Scalar			*m_aOutputParticleDensity;
	Scalar			*m_aOutputParticleTemperature;
	Scalar			*m_aOutputParticleKineticViscosity;
	Scalar			*m_aOutputParticleSolidPhaseRate;
	ParticleType	*m_aOutputParticleType;

	//ENDS for individual Output particles----------------------E&T Nepal August 2011--------------------------

	Integer m_ParticleNum;
	Integer m_MaxParticleNum;
		
	HANDLE m_Output; 
	CParameter* m_pParameter;
	CConfiguration* m_pConfiguration;

	CParameterCoefficients* m_pParameterCoefficient;

	CTriangle * m_haTriangles;// Cell and Grid of the STLTriangles
	Integer m_TriangleNum;

	//Triangle Parameters
	CTriangleParameters * m_haTriangleParameters;
	Integer m_TriangleModelNumber;

	CDragTriangle* m_haDragTriangles;
	Integer m_DragTriangleNum;	
	DragParameter	*m_haSTLDragPrameter;		//added by Arpan 20130128

	//For Multiple Bucket
	MultipleBucket * m_MultiBucket;
	Integer m_NumberOfBuckets;
	std::vector<Integer> m_TrianglePerModel;
	std::vector<Scalar3> m_ModelPosition;

	COutputFrame* m_aOutputFrame;// To store the Output
	Integer m_OutputFrameNum;
	COutputFrame* m_pOutputFrame;

	Integer m_StartStep;

	Integer m_NextItr;

public:
	void SetComputeZone(CBox* ComputeZone);
	void SetComputeZone(CParticle* aParticle, Integer ParticleNum, CTriangle* aTriangle, Integer TriangleNum);
	
	CCTStatusType CalculateParameterCoefficient();
	
	CCTStatusType Initialize(Integer ParticleNum,Integer * haParticleID, Scalar3* haParticlePosition, Scalar3* haParticleVelocity,
										Scalar* haParticlePressure, Scalar* haParticleDensity, Scalar* haParticleTemperature,
										Scalar* haParticleKineticViscosity, Scalar* haParticleSolidPhaseRate, ParticleType* haParticleType);
	CCTStatusType Calculate();
	CCTStatusType AddToOutputBuffer(Integer BufferSize,
									Integer* ParticleBufferID, Scalar3* ParticleBufferPosition,
									Scalar3* ParticleBufferVelcity, Scalar* ParticleBufferPressure,
									Scalar* ParticleBufferDensity, Scalar* ParticleBufferTemperature, Scalar* ParticleBufferKineticViscosity,
									Scalar* ParticleBufferSolidPhaseRate, ParticleType* ParticleBufferType);
	CCTStatusType WaitForOutputBuffer();
	
private:
	void Check(const CTriangle* pTriangle);
	void Check(const CParticle* const pParticle);
	CCTStatusType InitializeGrids();	
	CCTStatusType InitializeHostMemory();
	void InitializeComputeZone(Integer3 Section, CParameter* pParameter);
	//void InitializeParticleData(CParticle* aParticle, Integer ParticleNum);
	void InitializeParticleData(Integer ParticleNum,
										Integer * haParticleID, Scalar3* haParticlePosition, Scalar3* haParticleVelocity,
										Scalar* haParticlePressure, Scalar* haParticleDensity, Scalar* haParticleTemperature,
										Scalar* haParticleKineticViscosity, Scalar* haParticleSolidPhaseRate, ParticleType* haParticleType);
	void InitializeTriangleData();
	void InitializeDragTriangleData();
	//CCTStatusType InitializeOutputMemory();
	CCTStatusType WaitForOutput();	
	CCTStatusType Output(Integer Step, std::string OutputType);
	void OutputComputeZone();
	CCTStatusType ResetOutput();

	CCTStatusType InitializeOutputBuffer();
	void InitializeStartStep();
	void WaitForOutputFrame();
	CCTStatusType WaitForOutputThreads();

	void SetModelID(const Integer Model);

public:
	CCTStatusType Run(int argc, char** argv);

private:
	void SetModelPosition(Integer Model, const Scalar3* const Position );

	CCTStatusType CalculateLambdaValueLocal();
	CCTStatusType AddBucketParticles(const CParticle * bucket, const Integer bucketNumber);
	void InsertBucket(const CParticle * bucket, const Integer bucketNumber);
	CCTStatusType SetAddedBucket();
	CCTStatusType Synchronize();

	DtFrame* m_pDt;
	CCTStatusType CheckAndUpdateDt(Integer CurrentStep);
public:
	const DtFrame* GetDtFrame() const;
	
	void InitializeCells();
	//CCTStatusType InitializeGPU();
	CCTStatusType RegisterParticleTopology(Integer step);
	CCTStatusType RegisterTriangleTopology();
	CCTStatusType CalcExplicitly();
	CCTStatusType MoveTriangles();

	CCTStatusType RestWallPositions(Integer AnalysisStep); //For Bucket

	CCTStatusType ParticleDeviceToHost();
	CCTStatusType ParticleHostToDevice(); //FOr Bucket

	CCTStatusType ParticleNumberToContantMemory();

	CCTStatusType CalculateSTLDistance();
	void InitializeParameters();
	
	/////////////////////////////////////
	// system setting
	/////////////////////////////////////
	SystemUtility::CDeviceSelector m_DeviceSelector;

	////////Calculate Explicit//////////////////
	CCTStatusType CalcExplicitPressure();
	CCTStatusType CalcExplicitPressureGradient();
	CCTStatusType CalcTemperatureFactor();
	CCTStatusType CalcDragEffect();
	
	CCTStatusType LoadMultiBucket();
	CCTStatusType RotateTrianglePosition(const Integer analysisStep);
	CCTStatusType GetDataOnMultiGPUHandler(bool Output);

	//Turbulance
	CCTStatusType CalcTurbulenceViscosity();
	
	CCTStatusType ConstantsDeviceToHost(Integer Step);
	CCTStatusType TransferRequiredCheck(CCTStatusType Error,Integer Step);
	CCTStatusType SetOutputPath();
	public:
		bool m_IsAnyTriangleLinearMovable;
		bool m_IsAnyTriangleRotate;
		bool m_IsReTriangleRegisterNecessary;
};
