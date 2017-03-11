#pragma once
#define MAXNEIGHBORGRIDS 26
class CMultiGPUHandler;
class CCellSet
{

public:
	CCellSet();
	~CCellSet();
	std::string m_Outpath;
	CParameter m_Parameter;
	CParameterCoefficients m_ParameterCoefficient;
	CGridBox m_BoundingBox;

	int*  m_dGridParticleHash; // grid hash value for each particle
	int*  m_dGridParticleIndex;// particle index for each particle
	int*  m_dCellStart;        // index of start of each cell in sorted list
	int*  m_dCellEnd;          // index of end of cell

	int* m_haCellStart;
	int* m_haCellEnd;

	int   m_gridSortBits;
	int* m_hGridParticleHash;
	int* m_hGridParticleIndex;

	CCell* m_haCell;
	CCell* m_daCell;
	Integer m_CellNum;

	Integer m_ParticleNum; //Number of current Particles
	Integer m_MaxParticleNum;//Maximum allocated size

	
	CTriangle * m_haTriangles;// Cell and Grid of the STLTriangles
	CTriangle * m_daTriangles;// Cell and Grid of the STLTriangles
	CTriangle * m_daOriginalTriangles;// Cell and Grid of the STLTriangles
	Integer m_TriangleNum;

	
	CTriangleParameters * m_haTriangleParameters;
	CTriangleParameters * m_daTriangleParameters;

	DragParameter	*m_haSTLDragPrameter;
	DragParameter	*m_daSTLDragPrameter;
	Integer m_DragTriangleNum;

	Integer m_TriangleModelNumber;

	CDragTriangle * m_haDragTriangles;// Cell and Grid of the STLTriangles
	CDragTriangle * m_daDragTriangles;// Cell and Grid of the STLTriangles

	Scalar3* m_daDragAcc;
	Scalar3* m_haDragAcc;
	Scalar* m_daDragTemperature;
	Scalar* m_haDragTemperature;
	
	CDistance *m_daSTLDistance;// Holds the distance between particle and triangle
	Integer* m_daSTLID;// Holds the ID of the neighbor triangle	

	CDistance *m_haSTLDistance;
	Integer* m_haSTLID;

	
	Integer			*m_haParticleID;
	Scalar3			*m_haParticlePosition;
	Scalar3			*m_haParticleVelocity;
	Scalar			*m_haParticlePressure;
	Scalar			*m_haParticleDensity;
	Scalar			*m_haParticleTemperature;
	Scalar			*m_haParticleKineticViscosity;
	Scalar			*m_haParticleSolidPhaseRate;
	ParticleType	*m_haParticleType;

	//TurbulanceF
	Scalar			*m_haParticleTurbulaceViscosity;
	Scalar			*m_haParticleStrainTensorProduct;

	Integer			*m_haMagnifierCount;

	Integer			*m_daParticleID;
	Scalar3			*m_daParticlePosition;
	Scalar3			*m_daParticleVelocity;
	Scalar			*m_daParticlePressure;
	Scalar			*m_daParticleDensity;
	Scalar			*m_daParticleTemperature;
	Scalar			*m_daParticleKineticViscosity;
	Scalar			*m_daParticleSolidPhaseRate;
	ParticleType	*m_daParticleType;
	//Turbulance
	Scalar			*m_daParticleTurbulaceViscosity;
	Scalar			*m_daParticleStrainTensorProduct;
	Integer			*m_daMagnifierCount;////AA

	Integer			*m_daOutputParticleID;
	Scalar3			*m_daOutputParticlePosition;
	Scalar3			*m_daOutputParticleVelocity;
	Scalar			*m_daOutputParticlePressure;
	Scalar			*m_daOutputParticleDensity;
	Scalar			*m_daOutputParticleTemperature;
	Scalar			*m_daOutputParticleKineticViscosity;
	Scalar			*m_daOutputParticleSolidPhaseRate;
	ParticleType	*m_daOutputParticleType;

	cudaStream_t m_CudaStream;
	cudaStream_t m_CudaStream1;
	cublasHandle_t m_CublasHandle;
	Integer3 m_GridID;	
public:
	void SetGridID(Integer x, Integer y, Integer z);
	void SetNeighborGrids(CCellSet* aGrid, Integer3 Section);
	void SetComputeZone(const CBox* const ComputeZone, CParameter* pParameter/*,CParameterCoefficients* pParameterCoefficient*/,  const CBox* const ParentZone);
	void SetStartStep(const Integer StartStep);
	CCTStatusType InitializeHostMemory(Integer MaxParticleNum);
	void UninitializeHostMemory();
	
	void SetTriangles(CTriangle* Triangles, Integer TriangleNum);
	void SetDragTriangles(CDragTriangle* DragTriangles, Integer DragTriangleNum/*,CInnerPressureModel* InnerPressureModel*/);

	CCTStatusType SetParticleData(Integer ParticleNum,
		Integer * haParticleID, Scalar3* haParticlePosition, Scalar3* haParticleVelocity,
		Scalar* haParticlePressure, Scalar* haParticleDensity, Scalar* haParticleTemperature,
		Scalar* haParticleKineticViscosity, Scalar* haParticleSolidPhaseRate, ParticleType* haParticleType);
private:
	CCTStatusType InitializeDeviceMemory();
	CCTStatusType UnInitializeDeviceMemory();
	CCTStatusType InitializeDeviceConstantMemory();
public:
	CCTStatusType ParticleDeviceToHost();
	CCTStatusType ParticleHostToDevice();
	CCTStatusType ConstantsDeviceToHost(Integer Step);
public:
	CCTStatusType CalculateSTLDistance();
public:
	CCTStatusType RegisterParticleTopology(Integer step);
	CCTStatusType RegisterTriangleTopology();

	CCTStatusType CudaSafeCall(cudaError_t Status);
	CCTStatusType ResetWallPosition(const Integer AnalysisStep);
	CCTStatusType ParticleNumberToContantMemory();

	
private:
	Integer m_DeviceID;
	Integer m_StartStep;
	HANDLE m_AddParticles;
	Integer m_ID;
public:

	CCTStatusType CheckParticleOutsideComputeZone();

	void InitializeCells(int DeviceId, int GridId,int GridNum, CGridBox& TotalGridBox);
	void InitializeParameter(CParameter* pParameter, CParameterCoefficients* pParameterCoeff);
	CCTStatusType InitializeGPU(Integer deviceID);
	CCTStatusType FinalizeGPU();
	CCTStatusType MoveTriangles();

	//FOr Rotation
	CCTStatusType RotateTrianglePosition(const Integer analysisStep);

	void AddToBuffer(bool Output, Integer &BufferSize,Integer * ParticleBufferID,
		Scalar3* ParticleBufferPosition,
		Scalar3* ParticleBufferVelcity, Scalar* ParticleBufferPressure,
		Scalar* ParticleBufferDensity, Scalar* ParticleBufferTemperature, Scalar* ParticleBufferKineticViscosity,
		Scalar* ParticleBufferSolidPhaseRate, ParticleType* ParticleBufferType);
	CCTStatusType StreamSynchronize();

	//Core Calculation
	CCTStatusType CalculateDragEffect();
	CCTStatusType CalcExplicitly();
	CCTStatusType CalcExplicitPressure();
	CCTStatusType CalcExplicitPressureGradient();
	CCTStatusType CalcTemperatureFactor();
	
	//Turbulance
	CCTStatusType CalcTurbulenceViscosity();
	Integer getCellNum();
	CCTStatusType SetCudaDevice();
	void SetTrianglesParameters(CTriangleParameters * TriangleParameter, Integer TriangleModelNumber);
	void SetDragParameters(DragParameter * dragParameter, Integer dragTriangleNumber);

	//void RegisterTriangle(const int TriangleNumber,			const CTriangle * Triangle,
	//							const CGridBox *const BoundingBox,	const CGridParams* const GridParams,
	//							Integer ID,							CCell* daCell,
	//							Integer CellNum);

	CCTStatusType SetAddedParticlesEvent();

	static __inline__ __host__ __device__ bool IsInclude(const CGridBox* const BoundingBox, const Scalar3 * const Position,const Scalar Tolerance)
	{	
		if(Position->x < (BoundingBox->m_ComputeZone.m_MinBound.x - Tolerance))
		{
			return false;
		}
		if(Position->y < (BoundingBox->m_ComputeZone.m_MinBound.y - Tolerance))
		{
			return false;
		}
		if(Position->z < (BoundingBox->m_ComputeZone.m_MinBound.z - Tolerance))
		{
			return false;
		}
		if(Position->x >= (BoundingBox->m_ComputeZone.m_MaxBound.x + Tolerance)) 
		{
			return false;
		}
		if(Position->y >= (BoundingBox->m_ComputeZone.m_MaxBound.y + Tolerance))
		{
			return false;
		}
		if(Position->z >= (BoundingBox->m_ComputeZone.m_MaxBound.z + Tolerance))
		{
			return false;
		}
		return true;
	}

	CCTStatusType SetOutPath(std::string& Path);
	//Test Only
	//std::stringstream filenameAfter;
	
};
