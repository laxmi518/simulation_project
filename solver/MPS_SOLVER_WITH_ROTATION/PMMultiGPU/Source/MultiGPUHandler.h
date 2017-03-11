#pragma once
#include "DeviceSelector.h"
class CCellSet;
class CConfiguration;
class COutputFrame;
class DtFrame;
class CInnerPressureModel;
class CMultiGPUHandler 
{
public:
	CMultiGPUHandler(void);
public:
	~CMultiGPUHandler(void);
private:
	CBox m_ComputeZone;
	Integer m_GridNum;
	CCellSet* m_aGrid;

	CParticle* m_OutputParticles;
	Integer m_ParticleNum;
	Integer m_MaxParticleNum;

	CParticle* m_haBucket; // Bucket particle to fill later if necessary
	Integer m_BucketNum;

	HANDLE m_Output; 
	CParameter* m_pParameter;
	CConfiguration* m_pConfiguration;

	CTriangle * m_haTriangles;// Cell and Grid of the STLTriangles
	//CTriangle * m_daTriangles;// Cell and Grid of the STLTriangles
	Integer m_TriangleNum;

	//Triangle Parameters
	CTriangleParameters * m_haTriangleParameters;
	Integer m_TriangleModelNumber;

	CDragTriangle* m_haDragTriangles;
	
	CInnerPressureModel* m_pInnerPressrueModel;

	//For Multiple Bucket
	MultipleBucket * m_MultiBucket;
	Integer m_NumberOfBuckets;
	//For MultiBucket Ends

	//Integer m_TrianglePerModel[10];// Stores the number of Triangles of each model
	//Integer m_ModelNum;// Number of Models Loaded
	std::vector<Integer> m_TrianglePerModel;
	std::vector<Scalar3> m_ModelPosition;

	COutputFrame* m_aOutputFrame;// To store the Output
	Integer m_OutputFrameNum;
	COutputFrame* m_pOutputFrame;

	Scalar m_StartStep;

	Integer m_NextItr;

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

	//Drag @Rajan
	DragParameter	*m_haSTLDragPrameter;		//added by Arpan 20130128
	Integer m_DragTriangleNum;
public:
	void SetComputeZone(CBox* ComputeZone);
	void SetComputeZone(CParticle* aParticle, Integer ParticleNum, CTriangle* aTriangle, Integer TriangleNum);
	CCTStatusType Initialize(CParticle* aParticle,Integer ParticleNum);
	CCTStatusType AddToOutputBuffer(CParticle* ParticleBuffer, Integer BufferSize);
	/*
	AddToOutputBuffer(Integer BufferSize,		  Integer* ParticleBufferID, Scalar3* ParticleBufferPosition,
												  Scalar3* ParticleBufferVelcity, Scalar* ParticleBufferPressure,
												  Scalar* ParticleBufferDensity, Scalar* ParticleBufferTemperature, Scalar* ParticleBufferKineticViscosity,
												  Scalar* ParticleBufferSolidPhaseRate, ParticleType* ParticleBufferType);

*/
	CCTStatusType WaitForOutputBuffer();
	//CCTStatusType InitializeBucket(CParticle* aBucket,Integer BucketNum);
private:
	void Check(const CTriangle* pTriangle);
	void Check(const CParticle* const pParticle);
	CCTStatusType InitializeGrids();	
	void InitializeNeighbour(Integer3 Section);
	void InitializeComputeZone(Integer3 Section, CParameter* pParameter);
	CCTStatusType InitializeHostMemory();
	void InitializeParticleData(CParticle* aParticle, Integer ParticleNum);
	void InitializeBucketData();
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
	//void SetModelVelocity(Integer Model, const Scalar3* const Velocity );
	//void SetModelTemperature(const Integer Model, const Scalar Temperature);
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

public:
	//void SetParamForSTLRotation(Integer Model, const Scalar3* CenterOfRotation,const Scalar3* const SecondPointOfRotation,const Scalar AngleOfRotation, const bool resetToOriginalState);

	CCTStatusType CalculateNew();

	CCTStatusType LoadMultiBucket();

	CCTStatusType BucketAddition(Integer AnalysisStep, Integer &AddedTImes, bool & IsBucketAdded);

	CCTStatusType GetDataOnMultiGPUHandler(bool Output);

	CCTStatusType ParticleDeviceToHost();

public:
	/////////////////////////////////////
	// system setting
	/////////////////////////////////////
	SystemUtility::CDeviceSelector m_DeviceSelector;
};
