#include "stdafx.h"
#include "CGPUCellSet.h"

#include "Utility_Inline.h"
#include "DeviceProcess.h"

#include "MultiGPUHandler.h"
#include "DtFrame.h"

#define INDIVIDUALTYPE
#define USECDUASTREAM

//http://faq.programmerworld.net/programming/win32-multithreading-and-synchronization.html

CCellSet::CCellSet()
:m_haCell(NULL), m_daCell(NULL), m_CellNum(0)
,m_ComputeCellNum(NULL)
,m_GhostCellNum(NULL)

,m_haParticle(NULL) 
,m_daParticle(NULL)
//,m_daOutputParticle(NULL)
,m_ParticleNum(0)
,m_MaxParticleNum(0)

,m_haBucket(NULL)
,m_daBucket(NULL)
,m_BucketNum(0)

,m_haTriangles(NULL)
,m_daTriangles(NULL)
,m_daOriginalTriangles(NULL)  //For storing original Triangles in Device Memory
,m_TriangleNum(0)

,m_haDragTriangles(NULL)
,m_daDragTriangles(NULL)
,m_DragTriangleNum(0)
,m_daDragAcc(NULL)
,m_daDragTemperature(NULL)

,m_daSubgridScaleStress(NULL)

,m_daSTLDistance(NULL)
,m_daSTLID(NULL)

,m_dParticleHash(NULL)// Stores Next ID's of the particle in the cell
,m_hParticleHash(NULL)

,m_dEps(NULL)
,m_hEps(NULL)

,m_dMInv(NULL)
,m_hMInv(NULL)
,m_dAii(NULL)
,m_hAii(NULL)

,m_ThreadHandle(NULL)
,m_MutexHandle(NULL)
,m_ThreadID(-1)
,m_End(true)
,m_NeighborGridNum(0)
,m_Output(NULL)
,m_Buffer(NULL)
,m_OutputBuffer(NULL)
,m_AddParticles(NULL)
,m_pOutputHandler(NULL)
,m_Ghost(NULL)
,m_dReduceData(NULL)
,m_hReduceData(NULL)
,m_Err(false)
,m_ParticleMass(0)
,m_StartStep(0)
,m_NextItr(0)
,m_hX(NULL)
,m_hX0(NULL)
,m_hB(NULL)
,m_pInnerPressureModel(NULL)
,m_OldInNumber(0)
,m_NewInNumber(0)
,m_OutPutNum(0)

#ifdef INDIVIDUALTYPE
,m_haParticleID(0)
,m_haParticlePosition(0)
,m_haParticleVelocity(0)
,m_haParticlePressure(0)
,m_haParticleDensity(0)
,m_haParticleTemperature(0)
,m_haParticleKineticViscosity(0)
,m_haParticleSolidPhaseRate(0)
,m_haParticleType(NULL)

,m_daParticleID(0)
,m_daParticlePosition(0)
,m_daParticleVelocity(0)
,m_daParticlePressure(0)
,m_daParticleDensity(0)
,m_daParticleTemperature(0)
,m_daParticleKineticViscosity(0)
,m_daParticleSolidPhaseRate(0)
,m_daParticleType(NULL)

,m_daOutputParticleID(0)
,m_daOutputParticlePosition(0)
,m_daOutputParticleVelocity(0)
,m_daOutputParticlePressure(0)
,m_daOutputParticleDensity(0)
,m_daOutputParticleTemperature(0)
,m_daOutputParticleKineticViscosity(0)
,m_daOutputParticleSolidPhaseRate(0)
,m_daOutputParticleType(NULL)

,m_haBucketID(0)
,m_haBucketPosition(0)
,m_haBucketVelocity(0)
,m_haBucketPressure(0)
,m_haBucketDensity(0)
,m_haBucketTemperature(0)
,m_haBucketKineticViscosity(0)
,m_haBucketSolidPhaseRate(0)
,m_haBucketType(NULL)

/*
,m_daBucketID(0)
,m_daBucketPosition(0)
,m_daBucketVelocity(0)
,m_daBucketPressure(0)
,m_daBucketDensity(0)
,m_daBucketTemperature(0)
,m_daBucketKineticViscosity(0)
,m_daBucketSolidPhaseRate(0)
,m_daBucketType(NULL)
*/
,m_CudaStream(NULL)

#endif
,m_haTriangleParameters(NULL)
,m_daTriangleParameters(NULL)
,m_TriangleModelNumber(0)

,m_haSTLDragPrameter(NULL)
,m_daSTLDragPrameter(NULL)

,m_DeviceID(0)
{
	m_BoundingBox.m_BufferedZone.m_MaxBound.x = FLT_MAX;
	m_BoundingBox.m_BufferedZone.m_MaxBound.y = FLT_MAX;
	m_BoundingBox.m_BufferedZone.m_MaxBound.z = FLT_MAX;

	m_BoundingBox.m_BufferedZone.m_MinBound.x = -FLT_MAX;
	m_BoundingBox.m_BufferedZone.m_MinBound.y = -FLT_MAX;
	m_BoundingBox.m_BufferedZone.m_MinBound.z = -FLT_MAX;

	m_BoundingBox.m_ComputeZone.m_MaxBound = make_Scalar3(0.0,0.0,0.0);
	m_BoundingBox.m_ComputeZone.m_MinBound = make_Scalar3(0.0,0.0,0.0);
}

CCellSet::~CCellSet()
{
	UninitializeHostMemory();
	UnInitializeDeviceMemory();
}

void CCellSet::SetGridID(Integer x, Integer y, Integer z)
{
	m_GridID.x = x;
	m_GridID.y = y;
	m_GridID.z = z;
}
void CCellSet::SetNeighborGrids(CCellSet* aGrid, Integer3 Section)
{
	Integer NID = 0;
	for(Integer i = -1 ; i <=1; ++i)
	{
		for(Integer j = -1; j <=1; ++j)
		{
			for(Integer k = -1; k <=1; ++k)
			{
				if( i == 0 && j == 0 && k == 0)
				{
					continue;
				}
				bool pass = true;
				Integer X = (Integer)m_GridID.x + i;
				if( X < 0 || X >= Section.x)
				{
					pass = false;
				}
				Integer Y = (Integer)m_GridID.y + j;
				if( Y < 0 || Y >= Section.y)
				{
					pass = false;
				}
				Integer Z = (Integer)m_GridID.z + k;
				if( Z < 0 || Z >= Section.z)
				{
					pass = false;
				}				
				if(pass)
				{
					Integer GID = (Integer)((Section.y * Section.z) * X + Section.z * Y + Z);
					m_aNeighborGrid[NID++] = &aGrid[GID];
				}				
			}
		}
	}
	m_NeighborGridNum = NID;
}
void CCellSet::SetComputeZone(const CBox* const ComputeZone, CParameter* pParameter, const CBox* const ParentZone)
{
	//m_pParameter = pParameter;
	memcpy((void*)&m_Parameter, (void*)pParameter,sizeof(CParameter));

	m_ParticleMass = m_Parameter.Density * pow(m_Parameter.InitialParticleDistance, 3);

	Scalar3 MaxBound = ComputeZone->m_MaxBound;
	Scalar3 MinBound = ComputeZone->m_MinBound;
	m_BoundingBox.m_ComputeZone.m_MaxBound = MaxBound;
	m_BoundingBox.m_ComputeZone.m_MinBound = MinBound;

	m_BoundingBox.m_ParentZone.m_MaxBound = ParentZone->m_MaxBound;
	m_BoundingBox.m_ParentZone.m_MinBound = ParentZone->m_MinBound;

	Scalar MaxVal;
	//#ifdef _DEBUG
	//	MaxVal = std::max(pParameter->LaplacianInfluenceRadiusCoefficient,pParameter->GradientInfluenceRadiusCoefficient);
	//#else
	MaxVal =  max(pParameter->LaplacianInfluenceRadiusCoefficient,pParameter->GradientInfluenceRadiusCoefficient);
	//#endif
	if(m_Parameter.bUpdateInnerPressure)
	{
		MaxVal = max(pParameter->SurfaceInfluenceRadiusCoefficient, MaxVal);
	}

	Scalar MaxInfluenceRadius = pParameter->InitialParticleDistance * MaxVal;

	/*m_BoundingBox.m_ExcludeZone.m_MaxBound = make_Scalar3(MaxBound.x - 2 * MaxInfluenceRadius, MaxBound.y - 2 * MaxInfluenceRadius, MaxBound.z - 2 * MaxInfluenceRadius);;
	m_BoundingBox.m_ExcludeZone.m_MinBound = make_Scalar3(MinBound.x + 2 * MaxInfluenceRadius, MinBound.y + 2 * MaxInfluenceRadius, MinBound.z + 2 * MaxInfluenceRadius);*/

	m_BoundingBox.m_PressureZone.m_MaxBound = make_Scalar3(MaxBound.x + MaxInfluenceRadius, MaxBound.y + MaxInfluenceRadius, MaxBound.z + MaxInfluenceRadius);
	m_BoundingBox.m_PressureZone.m_MinBound = make_Scalar3(MinBound.x - MaxInfluenceRadius, MinBound.y - MaxInfluenceRadius, MinBound.z - MaxInfluenceRadius);

	m_BoundingBox .m_BufferedZone.m_MaxBound = make_Scalar3(MaxBound.x + 2 * MaxInfluenceRadius, MaxBound.y + 2 * MaxInfluenceRadius, MaxBound.z + 2 * MaxInfluenceRadius);
	m_BoundingBox .m_BufferedZone.m_MinBound = make_Scalar3(MinBound.x - 2 * MaxInfluenceRadius, MinBound.y - 2 * MaxInfluenceRadius, MinBound.z - 2 * MaxInfluenceRadius);
	m_BoundingBox .m_CellSize = MaxInfluenceRadius;

	Scalar3 GridSize = make_Scalar3(
		(MaxBound.x - MinBound.x) / MaxInfluenceRadius + 4,
		(MaxBound.y - MinBound.y) / MaxInfluenceRadius + 4,
		(MaxBound.z - MinBound.z) / MaxInfluenceRadius + 4
		);
	m_BoundingBox .m_GridSize = make_Integer3(ceil(GridSize.x), ceil(GridSize.y), ceil(GridSize.z));
}

void CCellSet::SetStartStep(const Integer StartStep)
{
	m_StartStep = StartStep;
}

CCTStatusType CCellSet::InitializeHostMemory(Integer MaxParticleNum)
{
	CCTStatusType Status = CCT_NOERR;
	m_MaxParticleNum = MaxParticleNum;
	m_CellNum = m_BoundingBox .m_GridSize.x * m_BoundingBox .m_GridSize.y * m_BoundingBox .m_GridSize.z;
	m_haCell = new CCell[m_CellNum];
	if(!m_haCell)
	{
		return CCT_MEMERR;
	}
	for(Integer i = 0; i < m_CellNum; ++i)
	{
		m_haCell[i].m_TriangleNum = 0;
		m_haCell[i].m_HashID = -1;
		m_haCell[i].m_ParticleID = -1;
		for(Integer j = 0; j < CELLMAXTRIANGLENUM; ++j)
		{
			m_haCell[i].m_TriangleID[j] = -1;
		}
	}
	m_ComputeCellNum = (m_BoundingBox.m_GridSize.x - 2 * OFFSET )* (m_BoundingBox.m_GridSize.y - 2 * OFFSET)* (m_BoundingBox.m_GridSize.z - 2 * OFFSET);
	Integer MidCellNum = (m_BoundingBox.m_GridSize.x - 2 * OFFSET - 2 * GHOST_OFFSET )* (m_BoundingBox.m_GridSize.y - 2 * OFFSET - 2 * GHOST_OFFSET)* (m_BoundingBox.m_GridSize.z - 2 * OFFSET - 2 * GHOST_OFFSET);//(m_BoundingBox.m_GridSize.x - 2 * OFFSET - 2 * GHOST_OFFSET
	m_GhostCellNum = m_ComputeCellNum - MidCellNum;	

	m_hEps = new Scalar[m_MaxParticleNum];
	if(!m_hEps)
	{
		return CCT_MEMERR;
	}	
	/*
	Status = CudaSafeCall(cudaMallocHost((void**)&m_haParticle, m_MaxParticleNum * sizeof(CParticle)));
	CCT_ERROR_CHECK(Status);
	*/
	/*
	m_haParticle = new CParticle[m_MaxParticleNum];
	if(!m_haParticle)
	{
		return CCT_MEMERR;
	}
	*/
	m_hReduceData = new Scalar[BLOCK_MAX_DIM];
	if(!m_hReduceData)
	{
		return CCT_MEMERR;
	}
	m_hParticleHash = new Integer[m_MaxParticleNum];
	if(!m_hParticleHash)
	{
		return CCT_MEMERR;
	}
#ifdef INDIVIDUALTYPE	
	//@Rajan Individual Components
	m_haParticleID				= new Integer[m_MaxParticleNum];
	m_haParticlePosition		= new Scalar3[m_MaxParticleNum];
	m_haParticleVelocity		= new Scalar3[m_MaxParticleNum];
	m_haParticlePressure		= new Scalar[m_MaxParticleNum];
	m_haParticleDensity			= new Scalar[m_MaxParticleNum];
	m_haParticleTemperature		= new Scalar[m_MaxParticleNum];
	m_haParticleKineticViscosity= new Scalar[m_MaxParticleNum];
	m_haParticleSolidPhaseRate	= new Scalar[m_MaxParticleNum];
	m_haParticleType			= new ParticleType[m_MaxParticleNum];
#endif

	return CCT_NOERR;
}
void CCellSet::UninitializeHostMemory()
{
	std::cout <<"Deinitializing Host Memory\n";
	if(m_haCell)
	{
		delete(m_haCell);
	}
	//if(m_haParticle)
	//{
	//	delete(m_haParticle);
	//}
	if(m_hEps)
	{
		delete(m_hEps);
	}
	if(m_hX)
	{
		delete(m_hX);
	}
	if(m_hX0)
	{
		delete(m_hX0);
	}
	if(m_hB)
	{
		delete(m_hB);
	}
	if(m_hParticleHash)
	{
		delete(m_hParticleHash);
	}
	if(m_hReduceData)
	{
		delete(m_hReduceData);
	}
	if(m_hMInv)
	{
		delete (m_hMInv);
	}
#ifdef INDIVIDUALTYPE
	//@Rajan STARTS for individual particles
	if(m_haParticleID)
	{
		delete[] m_haParticleID;
		m_haParticleID = NULL;
	}
	if(m_haParticlePosition)
	{
		delete[] m_haParticlePosition;
		m_haParticlePosition = NULL;
	}
	if(m_haParticleVelocity)
	{
		delete[] m_haParticleVelocity;
		m_haParticleVelocity  = NULL;
	}
	if(m_haParticlePressure)
	{
		delete[] m_haParticlePressure;
		m_haParticlePressure = NULL;
	}
	if(m_haParticleDensity)
	{
		delete[] m_haParticleDensity;
		m_haParticleDensity = NULL;
	}
	if(m_haParticleTemperature)
	{
		delete[] m_haParticleTemperature;
		m_haParticleTemperature = NULL;
	}
	if(m_haParticleKineticViscosity)
	{
		delete[] m_haParticleKineticViscosity;
		m_haParticleKineticViscosity = NULL;
	}
	if(m_haParticleSolidPhaseRate)
	{
		delete[] m_haParticleSolidPhaseRate;
		m_haParticleSolidPhaseRate = NULL;
	}
	if(m_haParticleType)
	{
		delete[] m_haParticleType;
		m_haParticleType = NULL;
	}
	//@Rajan STARTS for individual particles
#endif
}
void CCellSet::SetTriangles(CTriangle* Triangles, Integer TriangleNum)
{
	m_haTriangles = Triangles;
	m_TriangleNum = TriangleNum;
}
void CCellSet::SetTrianglesParameters(CTriangleParameters * TriangleParameter, Integer TriangleModelNumber)
{
	m_haTriangleParameters = TriangleParameter;
	m_TriangleModelNumber = TriangleModelNumber;
}
void CCellSet::SetDragParameters(DragParameter * dragParameter, Integer dragTriangleNumber)
{
	m_haSTLDragPrameter = dragParameter;
	m_DragTriangleNum = dragTriangleNumber;
}
void CCellSet::SetDragTriangles(CDragTriangle* DragTriangles, Integer DragTriangleNum, CInnerPressureModel* InnerPressureModel)
{
	m_haDragTriangles = DragTriangles;
	m_DragTriangleNum = DragTriangleNum;
	m_pInnerPressureModel = InnerPressureModel;
}
void CCellSet::WaitToStop(void)
{	
	WaitForSingleObject(m_ThreadHandle,INFINITE);
}
CCTStatusType CCellSet::InitializeGPU(Integer deviceID)
{
	CCTStatusType Status;
	m_DeviceID = deviceID;

	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaStreamCreate(&m_CudaStream));
	CCT_ERROR_CHECK(Status);

	Status = this->InitializeDeviceMemory();
	CCT_ERROR_CHECK(Status);

	Status = this->InitializeDeviceConstantMemory();
	CCT_ERROR_CHECK(Status);

	std::cout<<"Max Number from CGPUCellSet is :- "<<m_MaxParticleNum<<std::endl;

	return CCT_NOERR;
}
CCTStatusType CCellSet::StreamSynchronize()
{
	CCTStatusType Status = CCT_NOERR;

	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaStreamSynchronize(m_CudaStream));
	CCT_ERROR_CHECK(Status);

	/*Status = CudaSafeCall(cudaThreadSynchronize());
	CCT_ERROR_CHECK(Status);*/

	return Status;
}


CCTStatusType CCellSet::InitializeDeviceMemory()
{
	CCTStatusType Status;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaMalloc((void**)&m_daCell,m_CellNum * sizeof(CCell)));
	CCT_ERROR_CHECK(Status);

	/*Status = CudaSafeCall(cudaMalloc((void**)&m_daOutputParticle,m_MaxParticleNum * sizeof(CParticle)));
	CCT_ERROR_CHECK(Status);*/

	Status = CudaSafeCall(cudaMalloc((void**)&m_daSTLDistance,m_MaxParticleNum * sizeof(CDistance)));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaMalloc((void**)&m_daSTLID,m_MaxParticleNum * sizeof(Integer)));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaMalloc((void**)&m_dParticleHash,m_MaxParticleNum * sizeof(Integer)));
	CCT_ERROR_CHECK(Status);

	if(m_TriangleNum > 0)
	{
		Status = CudaSafeCall(cudaMalloc((void**)&m_daTriangles, m_TriangleNum * sizeof(CTriangle)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpy(m_daTriangles,m_haTriangles, m_TriangleNum * sizeof(CTriangle), cudaMemcpyHostToDevice));
		CCT_ERROR_CHECK(Status);

		Status = CudaSafeCall(cudaMalloc((void**)&m_daOriginalTriangles, m_TriangleNum * sizeof(CTriangle)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpy(m_daOriginalTriangles,m_haTriangles, m_TriangleNum * sizeof(CTriangle), cudaMemcpyHostToDevice));
		CCT_ERROR_CHECK(Status);
		
		//SetTriangleParameters
		Status = CudaSafeCall(cudaMalloc((void**)&m_daTriangleParameters,m_TriangleModelNumber * sizeof(CTriangleParameters)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpy(m_daTriangleParameters,m_haTriangleParameters, m_TriangleModelNumber * sizeof(CTriangleParameters), cudaMemcpyHostToDevice));
		CCT_ERROR_CHECK(Status);
	}
	
	Status = CudaSafeCall(cudaMemcpy(m_daCell, m_haCell,m_CellNum * sizeof(CCell),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	
	Status = CudaSafeCall(cudaMemset(m_daSTLID,-1,m_MaxParticleNum * sizeof(Integer)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemset(m_dParticleHash,-1,m_MaxParticleNum * sizeof(Integer)));
	CCT_ERROR_CHECK(Status);
	/*Status = CudaSafeCall(cudaMemset(m_daOutputParticle,0,m_MaxParticleNum * sizeof(CParticle)));
	CCT_ERROR_CHECK(Status);*/


	//Drag@Rajan Starts
	if(m_DragTriangleNum > 0)
	{
		Status = CudaSafeCall(cudaMalloc((void**)&m_daDragAcc , m_MaxParticleNum * sizeof(Scalar3)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemset(m_daDragAcc , 0 , m_MaxParticleNum * sizeof(Scalar3)));
		CCT_ERROR_CHECK(Status);

		//For drag Temperature
		Status = CudaSafeCall(cudaMalloc((void**)&m_daDragTemperature, m_MaxParticleNum * sizeof(Scalar)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemset(m_daDragTemperature , 0 , m_MaxParticleNum * sizeof(Scalar)));
		CCT_ERROR_CHECK(Status);		

		Status = CudaSafeCall(cudaMalloc((void**)&m_daSTLDragPrameter,m_DragTriangleNum * sizeof(DragParameter)));
		CCT_ERROR_CHECK(Status);
		
		Status = CudaSafeCall(cudaMemcpy(m_daSTLDragPrameter , m_haSTLDragPrameter , m_DragTriangleNum * sizeof(DragParameter) , cudaMemcpyHostToDevice));
		CCT_ERROR_CHECK(Status);	
		
		/*
		for(int i = 0 ; i < m_DragTriangleNum ; i++)
		{
			Scalar * CoefficientMember = m_haSTLDragPrameter[i].CoeffOfMember;
			Status = CudaSafeCall(cudaMalloc((void**)&m_haSTLDragPrameter[i].CoeffOfMember,		sizeof(Scalar) * (m_haSTLDragPrameter[i].RankofPolynomialEquation + 1)));
			CCT_ERROR_CHECK(Status);

			Status = CudaSafeCall(cudaMemcpy(&m_daSTLDragPrameter[i],&m_haSTLDragPrameter[i],sizeof(DragParameter),cudaMemcpyHostToDevice));
			CCT_ERROR_CHECK(Status);

			Status = CudaSafeCall(cudaMemcpy(m_haSTLDragPrameter[i].CoeffOfMember,	CoefficientMember,	sizeof(Scalar) * (m_haSTLDragPrameter[i].RankofPolynomialEquation + 1),	cudaMemcpyHostToDevice));
			CCT_ERROR_CHECK(Status);
			m_haSTLDragPrameter[i].CoeffOfMember = CoefficientMember;		
		}
		*/
	}
	//Drag@Rajan Ends
	
#ifdef INDIVIDUALTYPE
	//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
	Status = CudaSafeCall(cudaMalloc((void**)&m_daParticleID,m_MaxParticleNum * sizeof(Integer)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daParticlePosition,m_MaxParticleNum * sizeof(Scalar3)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daParticleVelocity,m_MaxParticleNum * sizeof(Scalar3)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daParticlePressure,m_MaxParticleNum * sizeof(Scalar)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daParticleDensity,m_MaxParticleNum * sizeof(Scalar)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daParticleTemperature,m_MaxParticleNum * sizeof(Scalar)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daParticleKineticViscosity,m_MaxParticleNum * sizeof(Scalar)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daParticleSolidPhaseRate,m_MaxParticleNum * sizeof(Scalar)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daParticleType,m_MaxParticleNum * sizeof(ParticleType)));
	CCT_ERROR_CHECK(Status);
	//ENDS for individual particles----------------------E&T Nepal August 2011--------------------------

	//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
	Status = CudaSafeCall(cudaMalloc((void**)&m_daOutputParticleID,m_MaxParticleNum * sizeof(Integer)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daOutputParticlePosition,m_MaxParticleNum * sizeof(Scalar3)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daOutputParticleVelocity,m_MaxParticleNum * sizeof(Scalar3)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daOutputParticlePressure,m_MaxParticleNum * sizeof(Scalar)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daOutputParticleDensity,m_MaxParticleNum * sizeof(Scalar)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daOutputParticleTemperature,m_MaxParticleNum * sizeof(Scalar)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daOutputParticleKineticViscosity,m_MaxParticleNum * sizeof(Scalar)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daOutputParticleSolidPhaseRate,m_MaxParticleNum * sizeof(Scalar)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daOutputParticleType,m_MaxParticleNum * sizeof(ParticleType)));
	CCT_ERROR_CHECK(Status);
	//ENDS for individual particles----------------------E&T Nepal August 2011--------------------------

	//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
	/*
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleID,m_haParticleID,m_ParticleNum * sizeof(Integer),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticlePosition,m_haParticlePosition,m_ParticleNum * sizeof(Scalar3),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleVelocity,m_haParticleVelocity,m_ParticleNum * sizeof(Scalar3),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticlePressure,m_haParticlePressure,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleDensity,m_haParticleDensity,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleTemperature,m_haParticleTemperature,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleKineticViscosity,m_haParticleKineticViscosity,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleSolidPhaseRate,m_haParticleSolidPhaseRate,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleType,m_haParticleType,m_ParticleNum * sizeof(ParticleType),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	*/

	
	Status = CudaSafeCall(cudaMemcpy(m_daParticleID,m_haParticleID,m_ParticleNum * sizeof(Integer),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daParticlePosition,m_haParticlePosition,m_ParticleNum * sizeof(Scalar3),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daParticleVelocity,m_haParticleVelocity,m_ParticleNum * sizeof(Scalar3),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daParticlePressure,m_haParticlePressure,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daParticleDensity,m_haParticleDensity,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daParticleTemperature,m_haParticleTemperature,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daParticleKineticViscosity,m_haParticleKineticViscosity,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daParticleSolidPhaseRate,m_haParticleSolidPhaseRate,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daParticleType,m_haParticleType,m_ParticleNum * sizeof(ParticleType),cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);

	this->ParticleDeviceToHost();
	
	//ENDS for individual particles----------------------E&T Nepal August 2011--------------------------
#endif
	return CCT_NOERR;
}
CCTStatusType CCellSet::UnInitializeDeviceMemory()
{
	std::cout<<"DEINITIALIZE DEVIVE MEMORY\n";
	if(m_daCell)
	{
		cudaFree(m_daCell);
	}
	if(m_daParticle)
	{
		cudaFree(m_daParticle);
	}
	if(m_daBucket)
	{
		cudaFree(m_daBucket);
	}
	if(m_dEps)
	{
		cudaFree(m_dEps);
	}
	if(m_dX)
	{
		cudaFree(m_dX);
	}
	if(m_dX0)
	{
		cudaFree(m_dX0);
	}
	if(m_dB)
	{
		cudaFree(m_dB);
	}

	if(m_dMInv)
	{
		cudaFree(m_dMInv);
	}
	if(m_dAii)
	{
		cudaFree(m_dAii);
	}
	if(m_daSTLDistance)
	{
		cudaFree(m_daSTLDistance);
	}

	if(m_daSTLID)
	{
		cudaFree(m_daSTLID);
	}

	if(m_dParticleHash)
	{
		cudaFree(m_dParticleHash);
	}
	if(m_dReduceData)
	{
		cudaFree(m_dReduceData);
	}

	if(m_daDragTriangles)
	{
		cudaFree(m_daDragTriangles);
	}
	if(m_daDragAcc)
	{
		cudaFree(m_daDragAcc);
	}
	if(m_daDragTemperature)
	{
		cudaFree(m_daDragTemperature);
	}
	if(m_daParticleID)
	{
		cudaFree(m_daParticleID);
	}
	if(m_daParticlePosition)
	{
		cudaFree(m_daParticlePosition);
	}
	if(m_daParticleVelocity)
	{
		cudaFree(m_daParticleVelocity);
	}
	if(m_daParticleDensity)
	{
		cudaFree(m_daParticleDensity);
	}
	if(m_daParticlePressure)
	{
		cudaFree(m_daParticlePressure);
	}
	if(m_daParticleTemperature)
	{
		cudaFree(m_daParticleTemperature);
	}
	if(m_daParticleKineticViscosity)
	{
		cudaFree(m_daParticleKineticViscosity);
	}
	if(m_daParticleSolidPhaseRate)
	{
		cudaFree(m_daParticleSolidPhaseRate);
	}
	if(m_daParticleType)
	{
		cudaFree(m_daParticleType);
	}
	if(m_daOutputParticleID)
	{
		cudaFree(m_daOutputParticleID);
	}
	if(m_daOutputParticlePosition)
	{
		cudaFree(m_daOutputParticlePosition);
	}
	if(m_daOutputParticleVelocity)
	{
		cudaFree(m_daOutputParticleVelocity);
	}
	if(m_daOutputParticleDensity)
	{
		cudaFree(m_daOutputParticleDensity);
	}
	if(m_daParticlePressure)
	{
		cudaFree(m_daParticlePressure);
	}
	if(m_daOutputParticleTemperature)
	{
		cudaFree(m_daOutputParticleTemperature);
	}
	if(m_daOutputParticleKineticViscosity)
	{
		cudaFree(m_daOutputParticleKineticViscosity);
	}
	if(m_daOutputParticleSolidPhaseRate)
	{
		cudaFree(m_daOutputParticleSolidPhaseRate);
	}
	if(m_daOutputParticleType)
	{
		cudaFree(m_daOutputParticleType);
	}

	m_ParticleNum = 0;
	m_MaxParticleNum = 0;
	m_CellNum = 0;
	m_BucketNum = 0;
	return CCT_NOERR;
}
CCTStatusType CCellSet::InitializeDeviceConstantMemory()
{
	CCTStatusType Status;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = InitializeDeviceConstOutPutParticles(m_daOutputParticleID, m_daOutputParticlePosition, m_daOutputParticleVelocity, m_daOutputParticlePressure, m_daOutputParticleDensity, m_daOutputParticleTemperature,
												  m_daOutputParticleKineticViscosity, m_daOutputParticleSolidPhaseRate, m_daOutputParticleType);
	CCT_ERROR_CHECK(Status);
	
	Status = InitializeDeviceConstInputParticles(m_daParticleID, m_daParticlePosition, m_daParticleVelocity, m_daParticlePressure, m_daParticleDensity, m_daParticleTemperature, m_daParticleKineticViscosity,
												 m_daParticleSolidPhaseRate, m_daParticleType);
	CCT_ERROR_CHECK(Status);

	Status = InitializeDeviceMemConst(m_Parameter, m_ParticleNum, m_daTriangles, m_TriangleNum, m_daTriangleParameters,  m_MaxParticleNum, m_daSTLDistance, m_daSTLID,
									  m_daCell, m_CellNum, m_BucketNum, m_BoundingBox, m_dParticleHash);
	CCT_ERROR_CHECK(Status);

	Status = DragParametersToConst(m_daSTLDragPrameter, m_daDragAcc, m_daDragTemperature, m_DragTriangleNum);

	return CCT_NOERR;
}
CCTStatusType CCellSet::ParticleDeviceToHost()
{
		
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);	

	//std::cout<<"Particle Device to Host Started\t";

#ifdef USECDUASTREAM
	/*Status = StreamSynchronize();
	CCT_ERROR_CHECK(Status);
	*/

	Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleID			,m_daOutputParticleID,m_ParticleNum * sizeof(Integer), cudaMemcpyDeviceToHost, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_haParticlePosition		,m_daOutputParticlePosition ,m_ParticleNum * sizeof(Scalar3), cudaMemcpyDeviceToHost, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleVelocity		,m_daOutputParticleVelocity ,m_ParticleNum * sizeof(Scalar3), cudaMemcpyDeviceToHost, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_haParticlePressure		,m_daOutputParticlePressure ,m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleDensity		,m_daOutputParticleDensity ,m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleTemperature	,m_daOutputParticleTemperature ,m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleKineticViscosity,m_daOutputParticleKineticViscosity ,m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleSolidPhaseRate,m_daOutputParticleSolidPhaseRate ,m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleType			,m_daOutputParticleType ,m_ParticleNum * sizeof(ParticleType), cudaMemcpyDeviceToHost, m_CudaStream));
	CCT_ERROR_CHECK(Status);

	Status = StreamSynchronize();
	CCT_ERROR_CHECK(Status);
#else
	Status = CudaSafeCall(cudaMemcpy(m_haParticleID,m_daOutputParticleID,m_ParticleNum * sizeof(Integer), cudaMemcpyDeviceToHost));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_haParticlePosition,m_daOutputParticlePosition,m_ParticleNum * sizeof(Scalar3), cudaMemcpyDeviceToHost));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_haParticleVelocity,m_daOutputParticleVelocity,m_ParticleNum * sizeof(Scalar3), cudaMemcpyDeviceToHost));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_haParticlePressure,m_daOutputParticlePressure,m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_haParticleDensity,m_daOutputParticleDensity,m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_haParticleTemperature,m_daOutputParticleTemperature,m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_haParticleKineticViscosity,m_daOutputParticleKineticViscosity,m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_haParticleSolidPhaseRate,m_daOutputParticleSolidPhaseRate,m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_haParticleType,m_daOutputParticleType,m_ParticleNum * sizeof(ParticleType), cudaMemcpyDeviceToHost));
	CCT_ERROR_CHECK(Status);
#endif
	//std::cout<<"Particle Device to Host Ends\n";
	return CCT_NOERR;
}
CCTStatusType CCellSet::ParticleHostToDevice()
{
	CCTStatusType Status = CCT_NOERR;
	
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);	

	//std::cout<< "Particle Host to Device Starts\t";
		
#ifdef USECDUASTREAM
	/*Status = StreamSynchronize();
	CCT_ERROR_CHECK(Status);
	*/

	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleID			, m_haParticleID				,m_ParticleNum * sizeof(Integer)	,cudaMemcpyHostToDevice, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticlePosition		, m_haParticlePosition			,m_ParticleNum * sizeof(Scalar3)	,cudaMemcpyHostToDevice, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleVelocity		, m_haParticleVelocity			,m_ParticleNum * sizeof(Scalar3)	,cudaMemcpyHostToDevice, m_CudaStream));
	CCT_ERROR_CHECK(Status);		
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticlePressure		, m_haParticlePressure			,m_ParticleNum * sizeof(Scalar)		,cudaMemcpyHostToDevice, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleDensity		, m_haParticleDensity			,m_ParticleNum * sizeof(Scalar)		,cudaMemcpyHostToDevice, m_CudaStream));
	CCT_ERROR_CHECK(Status);				
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleTemperature	, m_haParticleTemperature		,m_ParticleNum * sizeof(Scalar)		,cudaMemcpyHostToDevice, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleKineticViscosity, m_haParticleKineticViscosity,m_ParticleNum * sizeof(Scalar)		,cudaMemcpyHostToDevice, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleSolidPhaseRate, m_haParticleSolidPhaseRate	,m_ParticleNum * sizeof(Scalar)		,cudaMemcpyHostToDevice, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleType			, m_haParticleType				,m_ParticleNum * sizeof(ParticleType),cudaMemcpyHostToDevice, m_CudaStream));
	CCT_ERROR_CHECK(Status);
	
	Status = StreamSynchronize();
	CCT_ERROR_CHECK(Status);
#else
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleID			, m_haParticleID				,m_ParticleNum * sizeof(Integer), cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticlePosition		, m_haParticlePosition			,m_ParticleNum * sizeof(Scalar3), cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleVelocity		, m_haParticleVelocity			,m_ParticleNum * sizeof(Scalar3), cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticlePressure		, m_haParticlePressure			,m_ParticleNum * sizeof(Scalar), cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleDensity		, m_haParticleDensity			,m_ParticleNum * sizeof(Scalar), cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleTemperature	, m_haParticleTemperature		,m_ParticleNum * sizeof(Scalar), cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleKineticViscosity,m_haParticleKineticViscosity  ,m_ParticleNum * sizeof(Scalar), cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleSolidPhaseRate	, m_haParticleSolidPhaseRate	,m_ParticleNum * sizeof(Scalar), cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpy(m_daOutputParticleType			,m_haParticleType				,m_ParticleNum * sizeof(ParticleType), cudaMemcpyHostToDevice));
	CCT_ERROR_CHECK(Status);

#endif
	//std::cout<<"Particle Host to Device Ends \n";
	return CCT_NOERR;
}
CCTStatusType CCellSet::SetParticleData(CParticle* haParticle, Integer ParticleNum)
{
	m_ParticleNum = ParticleNum;

#ifdef INDIVIDUALTYPE
	if(m_ParticleNum > 0)
	{
		for(int i = 0 ; i < m_ParticleNum ; ++i)
		{
			m_haParticleID[i]		= haParticle[i].ID;
			m_haParticlePosition[i] = haParticle[i].Position;
			m_haParticleVelocity[i] = haParticle[i].Velocity;
			m_haParticlePressure[i] = haParticle[i].Pressure;
			m_haParticleDensity[i]  = haParticle[i].Density;
			m_haParticleTemperature[i]		= haParticle[i].Temperature;
			m_haParticleKineticViscosity[i] = haParticle[i].KineticViscosity;
			m_haParticleSolidPhaseRate[i]	= haParticle[i].SolidPhaseRate;
			m_haParticleType[i]				= haParticle[i].Type;		
		}
		//@Rajan Ends for individual particles
	}
#endif

	//memcpy(m_haParticle, haParticle, ParticleNum * sizeof(CParticle));
	return CCT_NOERR;
}
CCTStatusType CCellSet::SetParticleData(Integer BufferSize,
									  Integer* ParticleBufferID, Scalar3* ParticleBufferPosition,
									  Scalar3* ParticleBufferVelcity, Scalar* ParticleBufferPressure,
									  Scalar* ParticleBufferDensity, Scalar* ParticleBufferTemperature, Scalar* ParticleBufferKineticViscosity,
									  Scalar* ParticleBufferSolidPhaseRate, ParticleType* ParticleBufferType)
{
	m_ParticleNum = BufferSize;

#ifdef INDIVIDUALTYPE
	if(m_ParticleNum > 0)
	{
		for(int i = 0 ; i < m_ParticleNum ; ++i)
		{
			/*
			m_haParticleID[i]				= ParticleBufferID[i];
			m_haParticlePosition[i]			= ParticleBufferPosition[i];
			m_haParticleVelocity[i]			= ParticleBufferVelcity[i];
			m_haParticlePressure[i]			= ParticleBufferPressure[i];
			m_haParticleDensity[i]			= ParticleBufferDensity[i];
			m_haParticleTemperature[i]		= ParticleBufferTemperature[i];
			m_haParticleKineticViscosity[i] = ParticleBufferKineticViscosity[i];
			m_haParticleSolidPhaseRate[i]	= ParticleBufferSolidPhaseRate[i];
			m_haParticleType[i]				= ParticleBufferType[i];
			*/

			memcpy(&m_haParticleID[i],					&ParticleBufferID[i],sizeof(Integer));
			memcpy(&m_haParticlePosition[i],			&ParticleBufferPosition[i],sizeof(Scalar3));
			memcpy(&m_haParticleVelocity[i],			&ParticleBufferVelcity[i],sizeof(Scalar3));
			memcpy(&m_haParticlePressure[i],			&ParticleBufferPressure[i],sizeof(Scalar));
			memcpy(&m_haParticleDensity[i],				&ParticleBufferDensity[i],sizeof(Scalar));
			memcpy(&m_haParticleTemperature[i],			&ParticleBufferTemperature[i],sizeof(Scalar));
			memcpy(&m_haParticleKineticViscosity[i],	&ParticleBufferKineticViscosity[i],sizeof(Scalar));
			memcpy(&m_haParticleSolidPhaseRate[i],		&ParticleBufferSolidPhaseRate[i],sizeof(Scalar));
			memcpy(&m_haParticleType[i],				&ParticleBufferType[i],sizeof(ParticleType));
		}	
	}
#endif

	//memcpy(m_haParticle, haParticle, ParticleNum * sizeof(CParticle));
	return CCT_NOERR;
}
CCTStatusType CCellSet::SetBucketData(CParticle* aBucket, Integer BucketNum)
{
	m_haBucket = aBucket;
	m_BucketNum = BucketNum;

	#ifdef INDIVIDUALTYPE
	if(m_BucketNum > 0)
	{
		//@Rajan Individual Components
		m_haBucketID				= new Integer[m_BucketNum];
		m_haBucketPosition			= new Scalar3[m_BucketNum];
		m_haBucketVelocity			= new Scalar3[m_BucketNum];
		m_haBucketPressure			= new Scalar[m_BucketNum];
		m_haBucketDensity			= new Scalar[m_BucketNum];
		m_haBucketTemperature		= new Scalar[m_BucketNum];
		m_haBucketKineticViscosity	= new Scalar[m_BucketNum];
		m_haBucketSolidPhaseRate	= new Scalar[m_BucketNum];
		m_haBucketType				= new ParticleType[m_BucketNum];
	
		for(Integer i = 0 ; i < m_BucketNum ; ++i)
		{
			m_haBucketID[i]					= aBucket[i].ID;
			m_haBucketPosition[i]			= aBucket[i].Position;
			m_haBucketVelocity[i]			= aBucket[i].Velocity;
			m_haBucketPressure[i]			= aBucket[i].Pressure;
			m_haBucketDensity[i]			= aBucket[i].Density;
			m_haBucketTemperature[i]		= aBucket[i].Temperature;
			m_haBucketKineticViscosity[i]	= aBucket[i].KineticViscosity;
			m_haBucketSolidPhaseRate[i]		= aBucket[i].SolidPhaseRate;
			m_haBucketType[i]				= aBucket[i].Type;
		}
	}
#endif	
	return CCT_NOERR;
}
CCTStatusType CCellSet::RelocateAllParticleTopology()
{
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = RelocateParticleData(m_CudaStream, m_ParticleNum);
	CCT_ERROR_CHECK(Status);	
	return CCT_NOERR;
}
CCTStatusType CCellSet::RegisterParticleTopology()
{
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = ResetParticleTopology(m_CudaStream, m_CellNum, m_daCell, m_dParticleHash, m_ParticleNum);
	CCT_ERROR_CHECK(Status);
	Status = ::RegisterParticleTopology(m_CudaStream, m_ParticleNum);
	CCT_ERROR_CHECK(Status);		
	return CCT_NOERR;
}
CCTStatusType CCellSet::RegisterTriangleTopology()
{
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = ResetTriangleTopology(m_CudaStream, m_CellNum, m_daCell);
	CCT_ERROR_CHECK(Status);
	Status = ::RegisterTriangleTopology(m_CudaStream, m_daTriangles,m_TriangleNum,m_daCell,m_CellNum);
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
}
CCTStatusType CCellSet::CalculateSTLDistance()
{
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaMemset(m_daSTLID, -1, m_ParticleNum * sizeof(Integer)));
	CCT_ERROR_CHECK(Status);
	Status = ::CalcSTLDistance(m_CudaStream, m_ParticleNum);
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
}
CCTStatusType CCellSet::WaitForOutput()
{
	if(m_Err)
	{
		return CCT_CUDAERR;
	}
	if(!m_End)
	{
		WaitForSingleObject(m_Output,INFINITE);
		ResetEvent(m_Output);
		return CCT_NOERR;
	}
	else
	{
		return CCT_NOERR;
	}
}
CCTStatusType CCellSet::WaitForGhost()
{
	if(!m_End)
	{
		WaitForSingleObject(m_Ghost,INFINITE);
		return CCT_NOERR;
	}
	else
	{
		return CCT_ETCERR;
	}

}
CParticle* CCellSet::GetGhostBuffer()
{
	return m_haParticle;
}
CCTStatusType CCellSet::SetOutputEvent()
{
	SetEvent(m_Output);
	return CCT_NOERR;
}
CCTStatusType CCellSet::ResetOutputEvent()
{
	ResetEvent(m_Output);
	return CCT_NOERR;
}
CCTStatusType CCellSet::ResetGhostEvent()
{
	ResetEvent(m_Ghost);
	return CCT_NOERR;
}
bool CCellSet::IsEnd()
{
	return m_End;
}
CCTStatusType CCellSet::CudaSafeCall(cudaError_t Status)
{
	if(cudaSuccess != Status)
	{
		printf("Device : %d\ :",m_DeviceID);
		printf(cudaGetErrorString(Status));	
		m_Err = true;
		return CCT_CUDAERR;	
	}
	return CCT_NOERR;
}
CCTStatusType CCellSet::SetBuffer()
{
	bool evt = SetEvent(m_Buffer);	
	return CCT_NOERR;
}
CCTStatusType CCellSet::ResetBuffer()
{
	ResetEvent(m_Buffer);
	return CCT_NOERR;
}
CCTStatusType CCellSet::WaitForBuffer()
{
	if(!m_End)
	{
		WaitForSingleObject(m_Buffer,INFINITE);
		return CCT_NOERR;
	}
	else
	{
		return CCT_ETCERR;
	}
	return CCT_NOERR;
}

CCTStatusType CCellSet::SetOutputBuffer()
{
	bool evt = SetEvent(m_OutputBuffer);	
	return CCT_NOERR;
}
CCTStatusType CCellSet::WaitForOutputBuffer()
{
	if(!m_End)
	{
		WaitForSingleObject(m_OutputBuffer,INFINITE);
		ResetEvent(m_OutputBuffer);
		return CCT_NOERR;
	}
	else
	{
		return CCT_ETCERR;
	}

	return CCT_NOERR;
}

CCTStatusType CCellSet::ResetWallPosition(Integer AnalysisStep)
{
	CCTStatusType Status;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	
	Status = ::ResetWallPosition(m_CudaStream,m_TriangleNum,AnalysisStep,m_daOriginalTriangles);
	CCT_ERROR_CHECK(Status);	
	return CCT_NOERR;
}
CCTStatusType CCellSet::WaitForAddedParticles()
{
	if(!m_End)
	{
		WaitForSingleObject(m_AddParticles,INFINITE);	
		ResetEvent(m_AddParticles);
		return CCT_NOERR;
	}
	else
	{
		return CCT_ETCERR;
	}
}
CCTStatusType CCellSet::SetAddedParticlesEvent()
{
	SetEvent(m_AddParticles);
	return CCT_NOERR;
}
CCTStatusType CCellSet::SetComputeType(ComputeType Type)
{
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaMemcpyToSymbol(COMPUTETYPE, &Type, sizeof(ComputeType)));
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
}
CCTStatusType CCellSet::CheckAndUpdateDt(Integer CurrentStep)
{
	const DtFrame* Dt = m_pOutputHandler->GetDtFrame();
	if(Dt)
	{
		Integer DtStep = Dt->GetStep(m_NextItr);
		if((DtStep >= 0) && CurrentStep >= DtStep)
		{
			const DeltaTime* DtTime = Dt->GetDtTime(m_NextItr);
			if(DtTime)
			{
				m_Parameter.Dt = DtTime->Dt;
				m_Parameter.CourantNumber = DtTime->CourantNumber;
				m_Parameter.ParticleGenerationStep = DtTime->InflowStep;
				m_Parameter.OutputStep = DtTime->OutputStep;
				++m_NextItr;

				CCTStatusType Status;
				Status = CudaSafeCall( cudaMemcpyToSymbol(PARAMETER, (void*)&m_Parameter, sizeof(CParameter)) );
				CCT_ERROR_CHECK(Status);
			}
		}
	}
	return CCT_NOERR;
}
CCTStatusType CCellSet::CalcSTLDistance()
{
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaMemset(m_daSTLID, -1, m_ParticleNum * sizeof(Integer)));
	CCT_ERROR_CHECK(Status);
	Status = ::CalcSTLDistance(m_CudaStream, m_ParticleNum);
	CCT_ERROR_CHECK(Status);
	return Status;
}
CCTStatusType CCellSet::CalcExplicitly()
{
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Scalar InnerAcceleration = 0.0;
	Status = ::CalcExplicitly(m_CudaStream,m_ParticleNum,InnerAcceleration);
	CCT_ERROR_CHECK(Status);
	return Status;
}
CCTStatusType CCellSet::CalcPressureExplicit()
{
	CCTStatusType Status = CCT_NOERR;
	
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = ::CalcPressureExplicit(m_CudaStream, m_ParticleNum);
	CCT_ERROR_CHECK(Status);
	return Status;
}
CCTStatusType CCellSet::CalcTemperatureFactor()
{
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = ::CalcTemperatureFactor(m_CudaStream, m_ParticleNum);
	CCT_ERROR_CHECK(Status);
	return Status;
}
CCTStatusType CCellSet::CalcPressureExplicitGradient()
{
	CCTStatusType Status = CCT_NOERR;
	
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = ::CalcImplicitExplicitly(m_CudaStream, m_ParticleNum);
	CCT_ERROR_CHECK(Status);
	return Status;
}
CCTStatusType CCellSet::UpdateTrianglePosition()
{
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = ::UpdateTrianglePosition(m_CudaStream, m_TriangleNum,m_daTriangles);
	CCT_ERROR_CHECK(Status);
	return Status;
}
CCTStatusType CCellSet::RotateTrianglePosition(Integer i)
{
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = ::RotateTrianglePosition(m_CudaStream, m_TriangleNum,i);
	CCT_ERROR_CHECK(Status);
	return Status;
}

CCTStatusType CCellSet::AddToOutputBuffer(CParticle* ParticleBuffer, Integer & ParticleNum)
{
	for(Integer i = 0;i < m_ParticleNum; ++ i)
	{
		Integer ID = m_haParticle[i].ID;
		if(ID < 0 || ID >= m_ParticleNum)
		{
			continue;
		}
		memcpy(&ParticleBuffer[ID], &m_haParticle[i],sizeof(CParticle));
		//memcpy(&m_OutputParticles[ID], &ParticleBuffer[i], sizeof(CParticle));
	}
	return CCT_NOERR;
}

CCTStatusType CCellSet::AddToOutputBuffer( bool Output, Integer BufferSize,
												  Integer* ParticleBufferID, Scalar3* ParticleBufferPosition,
												  Scalar3* ParticleBufferVelcity, Scalar* ParticleBufferPressure,
												  Scalar* ParticleBufferDensity, Scalar* ParticleBufferTemperature, Scalar* ParticleBufferKineticViscosity,
												  Scalar* ParticleBufferSolidPhaseRate, ParticleType* ParticleBufferType)
{	
	
	for(Integer i = 0;i < BufferSize; ++ i)
	{
		memcpy(&ParticleBufferID[i],				&m_haParticleID[i],sizeof(Integer));		

		//Deleting the partiles outside of the Compute Zone.
		if(Output && (!IsInclude(&m_BoundingBox.m_ComputeZone,&m_haParticlePosition[i],m_Parameter.Tolerance)))
		{
			ParticleBufferPosition[i].x = 0;
			ParticleBufferPosition[i].y = 0;
			ParticleBufferPosition[i].z = 0;			

			ParticleBufferVelcity[i].x = 0;
			ParticleBufferVelcity[i].y = 0;
			ParticleBufferVelcity[i].z = 0;

			ParticleBufferPressure[i]		= 0;
			ParticleBufferDensity[i]		= 0;
			ParticleBufferTemperature[i]	= 0;
		}
		else
		{
			memcpy(&ParticleBufferPosition[i],			&m_haParticlePosition[i],sizeof(Scalar3));
			memcpy(&ParticleBufferVelcity[i],			&m_haParticleVelocity[i],sizeof(Scalar3));

			memcpy(&ParticleBufferPressure[i],			&m_haParticlePressure[i],sizeof(Scalar));
			memcpy(&ParticleBufferDensity[i],			&m_haParticleDensity[i],sizeof(Scalar));
			memcpy(&ParticleBufferTemperature[i],		&m_haParticleTemperature[i],sizeof(Scalar));
		}
		memcpy(&ParticleBufferKineticViscosity[i],	&m_haParticleKineticViscosity[i],sizeof(Scalar));
		memcpy(&ParticleBufferSolidPhaseRate[i],	&m_haParticleSolidPhaseRate[i],sizeof(Scalar));
		memcpy(&ParticleBufferType[i],				&m_haParticleType[i],sizeof(ParticleType));
	}
	m_ParticleNum = BufferSize;
	return CCT_NOERR;
}
CCTStatusType CCellSet::ParticleNumberToContantMemory()
{
	CCTStatusType Status = CCT_NOERR;

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = ::ParticleNumberToConst(m_ParticleNum);
	CCT_ERROR_CHECK(Status);

	//Status = CudaSafeCall(cudaMemcpyToSymbol(PARTICLENUM, &m_ParticleNum, sizeof(Integer)));
	//CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
}
//Drag@Rajan
CCTStatusType CCellSet::CalculateDragEffect()
{
	CCTStatusType Status = CCT_NOERR;
	if(m_ParticleNum >0)
	{
		//Set Device ID
		Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
		CCT_ERROR_CHECK(Status);

		Status = ::CalcDragEffect(m_CudaStream, m_ParticleNum);
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
