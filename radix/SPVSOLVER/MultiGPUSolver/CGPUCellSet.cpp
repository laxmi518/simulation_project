#include "std.h"
#include "Debug.h"
#include "CGPUCellSet.h"

#include "Utility_Inline.h"
#include "DeviceProcess.h"

#include "MultiGPUHandler.h"
#include "DtFrame.h"
#include "IO.h"
#include "Configuration.h"
//http://faq.programmerworld.net/programming/win32-multithreading-and-synchronization.html

CCellSet::CCellSet()
:m_haCell(NULL), m_daCell(NULL), m_CellNum(0)

,m_gridSortBits(18)  // increase this for larger grids

,m_dGridParticleHash(0)
,m_dGridParticleIndex(0)
,m_hGridParticleHash(0)
,m_hGridParticleIndex(0)

,m_ParticleNum(0)
,m_MaxParticleNum(0)

,m_haTriangles(NULL)
,m_daTriangles(NULL)
,m_daOriginalTriangles(NULL)
,m_TriangleNum(0)

,m_haDragTriangles(NULL)
,m_daDragTriangles(NULL)
,m_DragTriangleNum(0)
,m_daDragAcc(NULL)
,m_haDragAcc(NULL)
,m_daDragTemperature(NULL)
,m_haDragTemperature(NULL)
,m_haSTLDragPrameter(NULL)
,m_daSTLDragPrameter(NULL)
,m_daSTLDistance(NULL)
,m_daSTLID(NULL)

,m_StartStep(0)
,m_DeviceID(0)

,m_haParticleID(0)
,m_haParticlePosition(0)
,m_haParticleVelocity(0)
,m_haParticlePressure(0)
,m_haParticleDensity(0)
,m_haParticleTemperature(0)
,m_haParticleKineticViscosity(0)
,m_haParticleSolidPhaseRate(0)
,m_haParticleType(NULL)
,m_haMagnifierCount(NULL)
,m_haParticleTurbulaceViscosity(NULL)
,m_haParticleStrainTensorProduct(NULL)
,m_daParticleTurbulaceViscosity(0)
,m_daParticleStrainTensorProduct(0)
,m_daMagnifierCount(NULL)//Added by Ambika

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
,m_haTriangleParameters(NULL)
,m_daTriangleParameters(NULL)
,m_TriangleModelNumber(0)
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
	cudaDeviceReset();
}
void CCellSet::SetGridID(Integer x, Integer y, Integer z)
{
	m_GridID.x = x;
	m_GridID.y = y;
	m_GridID.z = z;
}
void CCellSet::SetComputeZone(const CBox* const ComputeZone, CParameter* pParameter/*, CParameterCoefficients* pParameterCoefficient*/, const CBox* const ParentZone)
{
	//m_pParameter = pParameter;
	memcpy((void*)&m_Parameter, (void*)pParameter,sizeof(CParameter));
	//memcpy((void*)&m_ParameterCoefficient,(void*)pParameterCoefficient,sizeof(CParameterCoefficients));

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
	/*if(m_Parameter.bUpdateInnerPressure)
	{
		MaxVal = max(pParameter->SurfaceInfluenceRadiusCoefficient, MaxVal);
	}*/

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
	m_BoundingBox.m_GridSize = make_Integer3((Integer)ceil(GridSize.x), (Integer)ceil(GridSize.y), (Integer)ceil(GridSize.z));
}

void CCellSet::SetStartStep(const Integer StartStep)
{
	m_StartStep = StartStep;
}

CCTStatusType CCellSet::InitializeHostMemory(Integer MaxParticleNum)
{
	//CCTStatusType Status = CCT_NOERR;
	//m_MaxParticleNum = MaxParticleNum;
	//m_CellNum = m_BoundingBox .m_GridSize.x * m_BoundingBox .m_GridSize.y * m_BoundingBox .m_GridSize.z;
	//printf("\nTotal Cells : %d\n",m_CellNum);
	//m_haCell = new CCell[m_CellNum];
	//if(!m_haCell)
	//{
	//	return CCT_MEMERR;
	//}
	//m_haCellStart = new int[m_CellNum];
	//if(!m_haCellStart)
	//{
	//	return CCT_MEMERR;
	//}
	//m_haCellEnd = new int[m_CellNum];
	//if(!m_haCellEnd)
	//{
	//	return CCT_MEMERR;
	//}
	//for(Integer i = 0; i < m_CellNum; ++i)
	//{
	//	m_haCell[i].m_TriangleNum = 0;
	//	for(Integer j = 0; j < CELLMAXTRIANGLENUM; ++j)
	//	{
	//		m_haCell[i].m_TriangleID[j] = -1;
	//	}
	//}
	////////////////////For Radix Sort Starts/////////////////////////////////////////
	//m_hGridParticleHash = new int[m_MaxParticleNum];
	//if(!m_hGridParticleHash)
	//{
	//	return CCT_MEMERR;
	//}
	//m_hGridParticleIndex = new int[m_MaxParticleNum];
	//if(!m_hGridParticleIndex)
	//{
	//	return CCT_MEMERR;
	//}
	//if(!m_hGridParticleIndex)
	//{
	//	return CCT_MEMERR;
	//}
	//m_haParticleID					= new Integer[m_MaxParticleNum];
	//m_haParticlePosition			= new Scalar3[m_MaxParticleNum];
	//m_haParticleVelocity			= new Scalar3[m_MaxParticleNum];
	//m_haParticlePressure			= new Scalar[m_MaxParticleNum];
	//m_haParticleDensity				= new Scalar[m_MaxParticleNum];
	//m_haParticleTemperature			= new Scalar[m_MaxParticleNum];
	//m_haParticleKineticViscosity	= new Scalar[m_MaxParticleNum];
	//m_haParticleSolidPhaseRate		= new Scalar[m_MaxParticleNum];
	//m_haParticleType				= new ParticleType[m_MaxParticleNum];

	//m_haSTLDistance = new CDistance[m_MaxParticleNum];
	//m_haSTLID = new Integer[m_MaxParticleNum];
	//m_haMagnifierCount = new Integer[m_MaxParticleNum];
	//m_haDragTemperature = new Scalar[m_MaxParticleNum];
	//m_haDragAcc = new Scalar3[m_MaxParticleNum];
	//m_haParticleTurbulaceViscosity = new Scalar[m_MaxParticleNum];
	//m_haParticleStrainTensorProduct = new Scalar[m_MaxParticleNum];

	//return Status;

	CCTStatusType Status = CCT_NOERR;
	m_MaxParticleNum = MaxParticleNum;
	m_CellNum = m_BoundingBox .m_GridSize.x * m_BoundingBox .m_GridSize.y * m_BoundingBox .m_GridSize.z;
	printf("\nTotal Cells : %d\n",m_CellNum);
	//m_haCell = new CCell[m_CellNum];
	 size_t size =  m_CellNum * sizeof(CCell);   //changed laxmi kadariya
     
	cudaMallocHost((void **) &m_haCell,size);
	if(!m_haCell)
	{
		return CCT_MEMERR;
	}
	//m_haCellStart = new int[m_CellNum];
	cudaMallocHost((void **) &m_haCellStart, m_CellNum * sizeof(int));

	if(!m_haCellStart)
	{
		return CCT_MEMERR;
	}
	//m_haCellEnd = new int[m_CellNum];
	cudaMallocHost((void **) &m_haCellEnd, m_CellNum * sizeof(int));

	
	if(!m_haCellEnd)
	{
		return CCT_MEMERR;
	}
	for(Integer i = 0; i < m_CellNum; ++i)
	{
		m_haCell[i].m_TriangleNum = 0;
		for(Integer j = 0; j < CELLMAXTRIANGLENUM; ++j)
		{
			m_haCell[i].m_TriangleID[j] = -1;
		}
	}
	//////////////////For Radix Sort Starts/////////////////////////////////////////
	//m_hGridParticleHash = new int[m_MaxParticleNum];
      cudaMallocHost((void **) &m_hGridParticleHash, m_MaxParticleNum * sizeof(int)); //laxmi kadariya
	if(!m_hGridParticleHash)
	{
		return CCT_MEMERR;
	}
	//m_hGridParticleIndex = new int[m_MaxParticleNum];
	 cudaMallocHost((void **) &m_hGridParticleIndex, m_MaxParticleNum * sizeof(int)); //laxmi kadariya
	if(!m_hGridParticleIndex)
	{
		return CCT_MEMERR;
	}
	if(!m_hGridParticleIndex)
	{
		return CCT_MEMERR;
	}
	//m_haParticleID					= new Integer[m_MaxParticleNum];
     cudaMallocHost((void **) &m_haParticleID, m_MaxParticleNum * sizeof(Integer)); //laxmi kadariya
	//m_haParticlePosition			= new Scalar3[m_MaxParticleNum];
	  cudaMallocHost((void **) &m_haParticlePosition, m_MaxParticleNum * sizeof(Scalar3)); //laxmi kadariya
	//m_haParticleVelocity			= new Scalar3[m_MaxParticleNum];
cudaMallocHost((void **) &m_haParticleVelocity, m_MaxParticleNum * sizeof(Scalar3)); //laxmi kadariya
	//m_haParticlePressure			= new Scalar[m_MaxParticleNum];
cudaMallocHost((void **) &m_haParticlePressure, m_MaxParticleNum * sizeof(Scalar)); //laxmi kadariya
	//m_haParticleDensity				= new Scalar[m_MaxParticleNum];
cudaMallocHost((void **) &m_haParticleDensity, m_MaxParticleNum * sizeof(Scalar)); //laxmi kadariya
	//m_haParticleTemperature			= new Scalar[m_MaxParticleNum];
cudaMallocHost((void **) &m_haParticleTemperature, m_MaxParticleNum * sizeof(Scalar)); //laxmi kadariya
	//m_haParticleKineticViscosity	= new Scalar[m_MaxParticleNum];
cudaMallocHost((void **) &m_haParticleKineticViscosity, m_MaxParticleNum * sizeof(Scalar)); //laxmi kadariya
	//m_haParticleSolidPhaseRate		= new Scalar[m_MaxParticleNum];
cudaMallocHost((void **) &m_haParticleSolidPhaseRate, m_MaxParticleNum * sizeof(Scalar)); //laxmi kadariya
	//m_haParticleType				= new ParticleType[m_MaxParticleNum];
cudaMallocHost((void **) &m_haParticleType, m_MaxParticleNum * sizeof(ParticleType)); //laxmi kadariya

	//m_haSTLDistance = new CDistance[m_MaxParticleNum];
cudaMallocHost((void **) &m_haSTLDistance, m_MaxParticleNum * sizeof(CDistance)); //laxmi kadariya
	//m_haSTLID = new Integer[m_MaxParticleNum];
cudaMallocHost((void **) &m_haSTLID, m_MaxParticleNum * sizeof(Integer)); //laxmi kadariya
	//m_haMagnifierCount = new Integer[m_MaxParticleNum];
cudaMallocHost((void **) &m_haMagnifierCount, m_MaxParticleNum * sizeof(Integer)); //laxmi kadariya
	//m_haDragTemperature = new Scalar[m_MaxParticleNum];
cudaMallocHost((void **) &m_haDragTemperature, m_MaxParticleNum * sizeof(Scalar)); //laxmi kadariya
	//m_haDragAcc = new Scalar3[m_MaxParticleNum];
cudaMallocHost((void **) &m_haDragAcc, m_MaxParticleNum * sizeof(Scalar3)); //laxmi kadariya
	//m_haParticleTurbulaceViscosity = new Scalar[m_MaxParticleNum];
cudaMallocHost((void **) &m_haParticleTurbulaceViscosity, m_MaxParticleNum * sizeof(Scalar)); //laxmi kadariya
	//m_haParticleStrainTensorProduct = new Scalar[m_MaxParticleNum];
cudaMallocHost((void **) &m_haParticleStrainTensorProduct, m_MaxParticleNum * sizeof(Scalar)); //laxmi kadariya

	return Status;
}

void CCellSet::UninitializeHostMemory()
{
	if(m_haCell)
	{
		//delete(m_haCell); original
		cudaFreeHost(m_haParticleID);
	}
	if(m_haParticleID)
	{
		cudaFreeHost(m_haParticleID);
	}
	if(m_haParticlePosition)
	{
		cudaFreeHost(m_haParticlePosition);
	}
	if(m_haParticleVelocity)
	{
		cudaFreeHost(m_haParticleVelocity);
	}
	if(m_haParticlePressure)
	{
		cudaFreeHost(m_haParticlePressure);
	}
	if(m_haParticleDensity)
	{
		cudaFreeHost(m_haParticleDensity);
	}
	if(m_haParticleTemperature)
	{
		cudaFreeHost(m_haParticleTemperature);
	}
	if(m_haParticleKineticViscosity)
	{
		cudaFreeHost(m_haParticleKineticViscosity);
	}
	if(m_haParticleSolidPhaseRate)
	{
		cudaFreeHost(m_haParticleSolidPhaseRate);
	}
	if(m_haParticleType)
	{
		cudaFreeHost(m_haParticleType);
	}
	//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
}
//
//CCTStatusType CCellSet::InitializeHostMemory(Integer MaxParticleNum)
//{
//	CCTStatusType Status = CCT_NOERR;
//	m_MaxParticleNum = MaxParticleNum;
//	m_CellNum = m_BoundingBox .m_GridSize.x * m_BoundingBox .m_GridSize.y * m_BoundingBox .m_GridSize.z;
//	printf("\nTotal Cells : %d\n",m_CellNum);
//	m_haCell = new CCell[m_CellNum];
//	if(!m_haCell)
//	{
//		return CCT_MEMERR;
//	}
//	m_haCellStart = new int[m_CellNum];
//	if(!m_haCellStart)
//	{
//		return CCT_MEMERR;
//	}
//	m_haCellEnd = new int[m_CellNum];
//	if(!m_haCellEnd)
//	{
//		return CCT_MEMERR;
//	}
//	for(Integer i = 0; i < m_CellNum; ++i)
//	{
//		m_haCell[i].m_TriangleNum = 0;
//		for(Integer j = 0; j < CELLMAXTRIANGLENUM; ++j)
//		{
//			m_haCell[i].m_TriangleID[j] = -1;
//		}
//	}
//	//////////////////For Radix Sort Starts/////////////////////////////////////////
//	m_hGridParticleHash = new int[m_MaxParticleNum];
//	if(!m_hGridParticleHash)
//	{
//		return CCT_MEMERR;
//	}
//	m_hGridParticleIndex = new int[m_MaxParticleNum];
//	if(!m_hGridParticleIndex)
//	{
//		return CCT_MEMERR;
//	}
//	if(!m_hGridParticleIndex)
//	{
//		return CCT_MEMERR;
//	}
//	m_haParticleID					= new Integer[m_MaxParticleNum];
//	m_haParticlePosition			= new Scalar3[m_MaxParticleNum];
//	m_haParticleVelocity			= new Scalar3[m_MaxParticleNum];
//	m_haParticlePressure			= new Scalar[m_MaxParticleNum];
//	m_haParticleDensity				= new Scalar[m_MaxParticleNum];
//	m_haParticleTemperature			= new Scalar[m_MaxParticleNum];
//	m_haParticleKineticViscosity	= new Scalar[m_MaxParticleNum];
//	m_haParticleSolidPhaseRate		= new Scalar[m_MaxParticleNum];
//	m_haParticleType				= new ParticleType[m_MaxParticleNum];
//
//	m_haSTLDistance = new CDistance[m_MaxParticleNum];
//	m_haSTLID = new Integer[m_MaxParticleNum];
//	m_haMagnifierCount = new Integer[m_MaxParticleNum];
//	m_haDragTemperature = new Scalar[m_MaxParticleNum];
//	m_haDragAcc = new Scalar3[m_MaxParticleNum];
//	m_haParticleTurbulaceViscosity = new Scalar[m_MaxParticleNum];
//	m_haParticleStrainTensorProduct = new Scalar[m_MaxParticleNum];
//
//	return Status;
//}
//void CCellSet::UninitializeHostMemory()
//{
//	if(m_haCell)
//	{
//		delete(m_haCell);
//	}
//	if(m_haParticleID)
//	{
//		delete(m_haParticleID);
//	}
//	if(m_haParticlePosition)
//	{
//		delete(m_haParticlePosition);
//	}
//	if(m_haParticleVelocity)
//	{
//		delete(m_haParticleVelocity);
//	}
//	if(m_haParticlePressure)
//	{
//		delete(m_haParticlePressure);
//	}
//	if(m_haParticleDensity)
//	{
//		delete(m_haParticleDensity);
//	}
//	if(m_haParticleTemperature)
//	{
//		delete(m_haParticleTemperature);
//	}
//	if(m_haParticleKineticViscosity)
//	{
//		delete(m_haParticleKineticViscosity);
//	}
//	if(m_haParticleSolidPhaseRate)
//	{
//		delete(m_haParticleSolidPhaseRate);
//	}
//	if(m_haParticleType)
//	{
//		delete(m_haParticleType);
//	}
//	//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
//}
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
void CCellSet::SetDragTriangles(CDragTriangle* DragTriangles, Integer DragTriangleNum/*, CInnerPressureModel* InnerPressureModel*/)
{
	m_haDragTriangles = DragTriangles;
	m_DragTriangleNum = DragTriangleNum;
	//m_pInnerPressureModel = InnerPressureModel;
}
CCTStatusType CCellSet::InitializeDeviceMemory()
{
	//************Radix Sort
	CCTStatusType Status;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_dGridParticleHash,m_MaxParticleNum * sizeof(int)));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaMalloc((void**)&m_dGridParticleIndex,m_MaxParticleNum * sizeof(int)));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaMalloc((void**)&m_dCellStart,m_CellNum * sizeof(int)));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaMalloc((void**)&m_dCellEnd,m_CellNum * sizeof(int)));
	CCT_ERROR_CHECK(Status);
	//***********Radix Sort

	//	CCTStatusType Status;
	Status = CudaSafeCall(cudaMalloc((void**)&m_daCell,m_CellNum * sizeof(CCell)));
	CCT_ERROR_CHECK(Status);


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

	Status = CudaSafeCall(cudaMalloc((void**)&m_daSTLDistance,m_MaxParticleNum * sizeof(CDistance)));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMalloc((void**)&m_daSTLID,m_MaxParticleNum * sizeof(Integer)));
	CCT_ERROR_CHECK(Status);

	if(m_TriangleNum > 0)
	{
		Status = CudaSafeCall(cudaMalloc((void**)&m_daTriangles, m_TriangleNum * sizeof(CTriangle)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpy(m_daTriangles,m_haTriangles, m_TriangleNum * sizeof(CTriangle), cudaMemcpyHostToDevice));
		CCT_ERROR_CHECK(Status);
		//For the Backup of OriginalTriangle
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
		if(m_Parameter.bTurbulance)
	{
		Status = CudaSafeCall(cudaMalloc((void**)&m_daParticleTurbulaceViscosity, m_MaxParticleNum * sizeof(Scalar)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemset(m_daParticleTurbulaceViscosity,0,m_MaxParticleNum * sizeof(Scalar)));
		CCT_ERROR_CHECK(Status);

		Status = CudaSafeCall(cudaMalloc((void**)&m_daParticleStrainTensorProduct, m_MaxParticleNum * sizeof(Scalar)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemset(m_daParticleStrainTensorProduct,0,m_MaxParticleNum * sizeof(Scalar)));
		CCT_ERROR_CHECK(Status);
	}
	Status = CudaSafeCall(cudaMemcpyAsync(m_daCell, m_haCell,m_CellNum * sizeof(CCell),cudaMemcpyHostToDevice,m_CudaStream));
	CCT_ERROR_CHECK(Status);
	//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleID,m_haParticleID,m_ParticleNum * sizeof(Integer),cudaMemcpyHostToDevice,m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticlePosition,m_haParticlePosition,m_ParticleNum * sizeof(Scalar3),cudaMemcpyHostToDevice,m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleVelocity,m_haParticleVelocity,m_ParticleNum * sizeof(Scalar3),cudaMemcpyHostToDevice,m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticlePressure,m_haParticlePressure,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice,m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleDensity,m_haParticleDensity,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice,m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleTemperature,m_haParticleTemperature,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice,m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleKineticViscosity,m_haParticleKineticViscosity,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice,m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleSolidPhaseRate,m_haParticleSolidPhaseRate,m_ParticleNum * sizeof(Scalar),cudaMemcpyHostToDevice,m_CudaStream));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleType,m_haParticleType,m_ParticleNum * sizeof(ParticleType),cudaMemcpyHostToDevice,m_CudaStream));
	CCT_ERROR_CHECK(Status);
	//ENDS for individual particles----------------------E&T Nepal August 2011--------------------------

	Status = CudaSafeCall(cudaMemsetAsync(m_daSTLID,-1,m_MaxParticleNum * sizeof(Integer),m_CudaStream));
	CCT_ERROR_CHECK(Status);

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

		Status = CudaSafeCall(cudaMalloc((void**)&m_daMagnifierCount, m_MaxParticleNum * sizeof(Integer)));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemset(m_daMagnifierCount , 0 , m_MaxParticleNum * sizeof(Integer)));
		CCT_ERROR_CHECK(Status);	
	}
	return CCT_NOERR;
}
CCTStatusType CCellSet::UnInitializeDeviceMemory()
{
	CCTStatusType Status = CCT_NOERR;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	if(m_daCell)
	{
		Status = CudaSafeCall(cudaFree(m_daCell));
		CCT_ERROR_CHECK(Status);
	}
	if(m_daSTLDistance)
	{
		Status = CudaSafeCall(cudaFree(m_daSTLDistance));
		CCT_ERROR_CHECK(Status);
	}

	if(m_daSTLID)
	{
		Status = CudaSafeCall(cudaFree(m_daSTLID));
		CCT_ERROR_CHECK(Status);
	}
	if(m_daParticleID)
	{
		Status = CudaSafeCall(cudaFree(m_daParticleID));
		CCT_ERROR_CHECK(Status);
	}
	if(m_daParticlePosition)
	{
		Status = CudaSafeCall(cudaFree(m_daParticlePosition));
		CCT_ERROR_CHECK(Status);
	}
	if(m_daParticleVelocity)
	{
		Status = CudaSafeCall(cudaFree(m_daParticleVelocity));
		CCT_ERROR_CHECK(Status);
	}
	if(m_daParticlePressure)
	{
		Status = CudaSafeCall(cudaFree(m_daParticlePressure));
		CCT_ERROR_CHECK(Status);
	}
	if(m_daParticleDensity)
	{
		Status = CudaSafeCall(cudaFree(m_daParticleDensity));
	}
	if(m_daParticleTemperature)
	{
		Status = CudaSafeCall(cudaFree(m_daParticleTemperature));
		CCT_ERROR_CHECK(Status);
	}
	if(m_daParticleKineticViscosity)
	{
		Status = CudaSafeCall(cudaFree(m_daParticleKineticViscosity));
		CCT_ERROR_CHECK(Status);
	}
	if(m_daParticleSolidPhaseRate)
	{
		Status = CudaSafeCall(cudaFree(m_daParticleSolidPhaseRate));
		CCT_ERROR_CHECK(Status);
	}
	if(m_daParticleType)
	{
		Status = CudaSafeCall(cudaFree(m_daParticleType));
		CCT_ERROR_CHECK(Status);
	}
	//prashant
	if(m_daParticleTurbulaceViscosity)
	{
		cudaFree(m_daParticleTurbulaceViscosity);
	}
	if(m_daParticleStrainTensorProduct)
	{
		cudaFree(m_daParticleStrainTensorProduct);
	}
	//ENDS for individual particles----------------------E&T Nepal August 2011--------------------------
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
	if(m_dGridParticleHash)
	{
		Status = CudaSafeCall(cudaFree(m_dGridParticleHash));
		CCT_ERROR_CHECK(Status);
	}
	if(m_dGridParticleIndex)
	{
		Status = CudaSafeCall(cudaFree(m_dGridParticleIndex));
		CCT_ERROR_CHECK(Status);
	}
	if(m_CublasHandle)
	{
		cublasStatus_t cs = cublasDestroy(m_CublasHandle);
		if(cs != CUBLAS_STATUS_SUCCESS)
		{
			return CCT_CUDAERR;
		}
	}
	Status = CudaSafeCall(cudaStreamDestroy(m_CudaStream));
	CCT_ERROR_CHECK(Status);

	//laxmi kadariya
	Status = CudaSafeCall(cudaStreamDestroy(m_CudaStream1));
	CCT_ERROR_CHECK(Status);

	m_ParticleNum = 0;
	m_MaxParticleNum = 0;
	m_CellNum = 0;
	return CCT_NOERR;
}
CCTStatusType CCellSet::InitializeDeviceConstantMemory()
{
	
	CCTStatusType Status = CCT_NOERR;
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = InitializeDeviceConstOutPutParticles(m_daOutputParticleID, m_daOutputParticlePosition, m_daOutputParticleVelocity, m_daOutputParticlePressure, m_daOutputParticleDensity,
												  m_daOutputParticleTemperature, m_daOutputParticleKineticViscosity, m_daOutputParticleSolidPhaseRate, m_daOutputParticleType);
	CCT_ERROR_CHECK(Status);
	
	Status = InitializeDeviceConstInputParticles(m_daParticleID, m_daParticlePosition, m_daParticleVelocity, m_daParticlePressure, m_daParticleDensity, m_daParticleTemperature,
												 m_daParticleKineticViscosity, m_daParticleSolidPhaseRate, m_daParticleType,m_daParticleTurbulaceViscosity, m_daParticleStrainTensorProduct);
	CCT_ERROR_CHECK(Status);
	
	Status = InitializeDeviceMemConst(m_Parameter,m_ParticleNum,m_daTriangles,m_TriangleNum,m_daTriangleParameters,  m_MaxParticleNum,m_daSTLDistance,m_daSTLID,
									  m_daCell,m_CellNum,m_dCellStart,m_dCellEnd,m_dGridParticleIndex,m_BoundingBox);
	CCT_ERROR_CHECK(Status);

	Status = DragParametersToConst(m_daSTLDragPrameter, m_daDragAcc, m_daDragTemperature, m_DragTriangleNum,m_daMagnifierCount,m_daDragTriangles);
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
}
CCTStatusType CCellSet::ParticleDeviceToHost()
{
	CCTStatusType Status = CCT_NOERR;
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);	

	if(m_ParticleNum > 0)
	{
		//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
		//Status = SetCudaDevice();
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleID,		m_daOutputParticleID,			m_ParticleNum * sizeof(Integer), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticlePosition,	m_daOutputParticlePosition,		m_ParticleNum * sizeof(Scalar3), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleVelocity,	m_daOutputParticleVelocity,		m_ParticleNum * sizeof(Scalar3), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticlePressure,	m_daOutputParticlePressure,		m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleDensity,	m_daOutputParticleDensity,		m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);

		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleTemperature,		m_daOutputParticleTemperature,			m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleKineticViscosity,	m_daOutputParticleKineticViscosity,		m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleSolidPhaseRate,	m_daOutputParticleSolidPhaseRate,		m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleType,				m_daOutputParticleType,					m_ParticleNum * sizeof(ParticleType), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);

		Status = StreamSynchronize();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CCellSet::ConstantsDeviceToHost(Integer Step)
{
	CCTStatusType Status = CCT_NOERR;
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);	
	bool inCodition =false;
	if(m_ParticleNum > 0)
	{
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleTurbulaceViscosity,	m_daParticleTurbulaceViscosity,		m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);

		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleStrainTensorProduct,	m_daParticleStrainTensorProduct,		m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haSTLDistance,		m_daSTLDistance,			m_ParticleNum * sizeof(CDistance), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haSTLID,	m_daSTLID,		m_ParticleNum * sizeof(Integer), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);

		Status = CudaSafeCall(cudaMemcpyAsync(m_haMagnifierCount,	m_daMagnifierCount,		m_ParticleNum * sizeof(Integer), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);

		Status = CudaSafeCall(cudaMemcpyAsync(m_haDragTemperature,	m_daDragTemperature,		m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);

		Status = CudaSafeCall(cudaMemcpyAsync(m_haDragAcc,	m_daDragAcc,		m_ParticleNum * sizeof(Scalar3), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);

		inCodition = true;
	}
	if(m_DragTriangleNum >0)
	{
		Status = CudaSafeCall(cudaMemcpyAsync(m_haSTLDragPrameter,	m_daSTLDragPrameter, m_DragTriangleNum * sizeof(DragParameter), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		inCodition = true;
	}
	if(m_CellNum > 0)
	{
		Status = CudaSafeCall(cudaMemcpyAsync(m_haCell,	m_daCell, m_CellNum * sizeof(CCell), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		inCodition = true;

	}
	if(m_TriangleModelNumber >0)
	{
		Status = CudaSafeCall(cudaMemcpyAsync(m_haTriangleParameters,	m_daTriangleParameters, m_TriangleModelNumber * sizeof(CTriangleParameters), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		inCodition = true;

	}
	if(m_TriangleNum > 0)
	{
		Status = CudaSafeCall(cudaMemcpyAsync(m_haTriangles,	m_daTriangles, m_TriangleNum * sizeof(CTriangle), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		inCodition = true;
	}
	if(inCodition)
	{
		Status = StreamSynchronize();
		CCT_ERROR_CHECK(Status);
	}
	std::stringstream Filename;
	Filename << m_Outpath<< "\\Constants"<<"_Step_"<<Step<<".txt";
	
	IO::SaveConstantData(m_ParticleNum,m_TriangleNum, m_DragTriangleNum, m_CellNum,m_TriangleModelNumber,m_haSTLDistance,m_haSTLID,m_haMagnifierCount,m_haDragTemperature,m_haDragAcc,m_haParticleTurbulaceViscosity,m_haParticleStrainTensorProduct,m_haSTLDragPrameter,m_haCell,m_haTriangleParameters,m_haTriangles,Filename.str(),Step);
	
	if(m_ParticleNum > 0)
	{
		//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
		//Status = SetCudaDevice();
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleID,		m_daOutputParticleID,			m_ParticleNum * sizeof(Integer), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticlePosition,	m_daOutputParticlePosition,		m_ParticleNum * sizeof(Scalar3), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleVelocity,	m_daOutputParticleVelocity,		m_ParticleNum * sizeof(Scalar3), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticlePressure,	m_daOutputParticlePressure,		m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleDensity,	m_daOutputParticleDensity,		m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);

		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleTemperature,		m_daOutputParticleTemperature,			m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleKineticViscosity,	m_daOutputParticleKineticViscosity,		m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleSolidPhaseRate,	m_daOutputParticleSolidPhaseRate,		m_ParticleNum * sizeof(Scalar), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_haParticleType,				m_daOutputParticleType,					m_ParticleNum * sizeof(ParticleType), cudaMemcpyDeviceToHost, m_CudaStream));
		CCT_ERROR_CHECK(Status);

		Status = StreamSynchronize();
		CCT_ERROR_CHECK(Status);
	}
	//std::stringstream ParticleFileName;
	//ParticleFileName << m_Outpath<< "\\ConstantsOutputParticle"<<"_Step_"<<Step<<".dat";
	//IO::SaveConstantOutputParticleData(ParticleFileName.str(),m_ParticleNum,m_haParticleID, m_haParticlePosition,m_haParticleVelocity,m_haParticlePressure,m_haParticleDensity,m_haParticleTemperature,m_haParticleKineticViscosity,m_haParticleSolidPhaseRate,m_haParticleType);
	return CCT_NOERR;
}
CCTStatusType CCellSet::ParticleHostToDevice()
{
	CCTStatusType Status = CCT_NOERR;
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	if(m_ParticleNum > 0)
	{
		//Status = SetCudaDevice();
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleID				, m_haParticleID				,m_ParticleNum * sizeof(Integer)	,cudaMemcpyHostToDevice, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticlePosition		, m_haParticlePosition			,m_ParticleNum * sizeof(Scalar3)	,cudaMemcpyHostToDevice, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleVelocity		, m_haParticleVelocity			,m_ParticleNum * sizeof(Scalar3)	,cudaMemcpyHostToDevice, m_CudaStream));
		CCT_ERROR_CHECK(Status);		
		Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticlePressure		, m_haParticlePressure			,m_ParticleNum * sizeof(Scalar)		,cudaMemcpyHostToDevice, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleDensity			, m_haParticleDensity			,m_ParticleNum * sizeof(Scalar)		,cudaMemcpyHostToDevice, m_CudaStream));
		CCT_ERROR_CHECK(Status);				
		Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleTemperature		, m_haParticleTemperature		,m_ParticleNum * sizeof(Scalar)		,cudaMemcpyHostToDevice, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleKineticViscosity, m_haParticleKineticViscosity	,m_ParticleNum * sizeof(Scalar)		,cudaMemcpyHostToDevice, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleSolidPhaseRate	, m_haParticleSolidPhaseRate	,m_ParticleNum * sizeof(Scalar)		,cudaMemcpyHostToDevice, m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemcpyAsync(m_daOutputParticleType			, m_haParticleType				,m_ParticleNum * sizeof(ParticleType),cudaMemcpyHostToDevice, m_CudaStream));
		CCT_ERROR_CHECK(Status);

		Status = StreamSynchronize();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}


CCTStatusType CCellSet::SetAddedParticlesEvent()
{
	SetEvent(m_AddParticles);
	return CCT_NOERR;
}
CCTStatusType CCellSet::SetParticleData(Integer ParticleNum,
										Integer * haParticleID, Scalar3* haParticlePosition, Scalar3* haParticleVelocity,
										Scalar* haParticlePressure, Scalar* haParticleDensity, Scalar* haParticleTemperature,
										Scalar* haParticleKineticViscosity, Scalar* haParticleSolidPhaseRate, ParticleType* haParticleType)
{
	m_ParticleNum = ParticleNum;
	if(m_ParticleNum > 0)
	{
		for(int i = 0 ; i < m_ParticleNum ; ++i)
		{
			memcpy(&m_haParticleID[i],					&haParticleID[i],sizeof(Integer));
			memcpy(&m_haParticlePosition[i],			&haParticlePosition[i],sizeof(Scalar3));
			memcpy(&m_haParticleVelocity[i],			&haParticleVelocity[i],sizeof(Scalar3));
			memcpy(&m_haParticlePressure[i],			&haParticlePressure[i],sizeof(Scalar));
			memcpy(&m_haParticleDensity[i],				&haParticleDensity[i],sizeof(Scalar));
			memcpy(&m_haParticleTemperature[i],			&haParticleTemperature[i],sizeof(Scalar));
			memcpy(&m_haParticleKineticViscosity[i],	&haParticleKineticViscosity[i],sizeof(Scalar));
			memcpy(&m_haParticleSolidPhaseRate[i],		&haParticleSolidPhaseRate[i],sizeof(Scalar));
			memcpy(&m_haParticleType[i],				&haParticleType[i],sizeof(ParticleType));
		}	
	}
	return CCT_NOERR;
}
CCTStatusType CCellSet::RegisterParticleTopology(Integer step)
{
	CCTStatusType Status = CCT_NOERR;
	/*Status = SetCudaDevice();*/

	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	Status = StreamSynchronize();
	CCT_ERROR_CHECK(Status);
	//Partilce Deletion that goes outside of compute zone
	Status = CheckParticleOutsideComputeZone();
	CCT_ERROR_CHECK(Status);
		
	// calculate grid hash //particleSystem.cpp L260
	Status = ::CaculateCellIDandInitializeHash( m_CudaStream,
												m_ParticleNum,
												m_CellNum,
												m_dGridParticleHash,
												m_dGridParticleIndex,
												m_daOutputParticlePosition
												);
	CCT_ERROR_CHECK(Status);
	
	Status = StreamSynchronize();
	CCT_ERROR_CHECK(Status);

	Status = ::SortUsingThrust(   m_ParticleNum,
										m_dGridParticleHash,
										m_dGridParticleIndex);
	CCT_ERROR_CHECK(Status);
	
	Status = StreamSynchronize();
	CCT_ERROR_CHECK(Status);

	// reorder particle arrays into sorted order and
    // find start and end of each cell
	Status =  reorderDataAndFindCellStart(	m_CudaStream,
											m_ParticleNum,
											m_CellNum,
											m_dGridParticleHash,
											m_dGridParticleIndex,
											m_dCellStart,
											m_dCellEnd
											);
	CCT_ERROR_CHECK(Status);
	
	Status = StreamSynchronize();
	CCT_ERROR_CHECK(Status);

	return CCT_NOERR;
}
CCTStatusType CCellSet::RegisterTriangleTopology()
{
	CCTStatusType Status = CCT_NOERR;
	//Status = SetCudaDevice();
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	//CCT_ERROR_CHECK(Status);
	Status = ResetTriangleTopology(m_CudaStream,m_CellNum, m_daCell);
	CCT_ERROR_CHECK(Status);
	Status = ::RegisterTriangleTopology(m_CudaStream,m_daTriangles,m_TriangleNum,m_daCell,m_CellNum);
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
}
CCTStatusType CCellSet::CalculateSTLDistance()
{
	CCTStatusType Status = CCT_NOERR;
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	if(m_ParticleNum >0)
	{
		//Status = SetCudaDevice();
		CCT_ERROR_CHECK(Status);
		Status = CudaSafeCall(cudaMemsetAsync(m_daSTLID, -1, m_ParticleNum * sizeof(Integer),m_CudaStream));
		CCT_ERROR_CHECK(Status);
		Status = CalcSTLDistance(m_CudaStream,m_ParticleNum);
		CCT_ERROR_CHECK(Status);

#ifdef _DEBUG
		std::stringstream fileName;
		fileName << "test\\STLDist_" <<m_ID <<"_.csv";
		Debug::DumpSTLDistMagintude(fileName.str(),m_ParticleNum,m_daSTLDistance);
#endif
	}
	return CCT_NOERR;
}
CCTStatusType CCellSet::CudaSafeCall(cudaError_t Status)
{
		
	if(Status == cudaErrorInvalidValue )
	{
		printf("cudaErrorInvalidValue");
	}
	if(cudaSuccess != Status)
	{
		printf("Device : %d\ :",m_DeviceID);
		//printf("Device : %d :",m_ID);
		printf(cudaGetErrorString(Status));	
		return CCT_CUDAERR;	
	}
	return CCT_NOERR;
}
CCTStatusType CCellSet::ResetWallPosition(const Integer AnalysisStep)
{
	CCTStatusType Status = CCT_NOERR;
/*	Status = SetCudaDevice();
	CCT_ERROR_CHECK(Status)*/;
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	Status = ::ResetWallPosition(m_CudaStream,m_TriangleNum,AnalysisStep,m_daOriginalTriangles);
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
}

void CCellSet::InitializeCells(int DeviceId, int GridId,int GridNum, CGridBox& TotalGridBox) //Max and Min region of Box ,TotalGridBox gives the compute zone.
{
	m_ID = DeviceId;
	m_BoundingBox = TotalGridBox;
	Integer3 MyGridSize = TotalGridBox.m_GridSize;
	MyGridSize.x = TotalGridBox.m_GridSize.x / GridNum;
	m_BoundingBox.m_GridSize = MyGridSize;

	Integer TotalCells = TotalGridBox.m_GridSize.x * TotalGridBox.m_GridSize.y * TotalGridBox.m_GridSize.z;
	m_CellNum = TotalCells;
}
Integer CCellSet::getCellNum()
{
	return m_CellNum;
}
CCTStatusType CCellSet::FinalizeGPU()
{
	CCTStatusType Status = CCT_NOERR;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	Status = UnInitializeDeviceMemory();
	CCT_ERROR_CHECK(Status);
	return Status;
}
CCTStatusType CCellSet::InitializeGPU(Integer deviceID)
{
	/*CCTStatusType Status = CCT_NOERR;
	Status = SetCudaDevice();
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaStreamCreate(&m_CudaStream));
	CCT_ERROR_CHECK(Status);*/
	CCTStatusType Status;
	m_DeviceID = deviceID;

	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = CudaSafeCall(cudaStreamCreate(&m_CudaStream));
	CCT_ERROR_CHECK(Status);
   //Laxmi kadariya
	Status = CudaSafeCall(cudaStreamCreate(&m_CudaStream1));
	CCT_ERROR_CHECK(Status);
	cublasStatus_t cs = cublasCreate(&m_CublasHandle);
	if(cs != CUBLAS_STATUS_SUCCESS)
	{
		return CCT_CUDAERR;
	}
	cs = cublasSetStream(m_CublasHandle, m_CudaStream); 
	if(cs != CUBLAS_STATUS_SUCCESS)
	{
		return CCT_CUDAERR;
	}
	Status = InitializeDeviceMemory();
	CCT_ERROR_CHECK(Status);
	Status = InitializeDeviceConstantMemory();
	CCT_ERROR_CHECK(Status);

	std::cout<<"Max Number from CGPUCellSet is :- "<<m_MaxParticleNum<<std::endl;

	return Status;
}
CCTStatusType CCellSet::CheckParticleOutsideComputeZone()
{
	CCTStatusType Status = CCT_NOERR;
	/*Status = SetCudaDevice();
	CCT_ERROR_CHECK(Status);*/
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	if(m_ParticleNum > 0)
	{
		Status = StreamSynchronize();
		CCT_ERROR_CHECK(Status);

		Status = ParticleDeviceToHost();
		CCT_ERROR_CHECK(Status);

		for(int index = 0 ; index < m_ParticleNum ; index++)
		{
			if(!IsInclude(&m_BoundingBox,&m_haParticlePosition[index],m_Parameter.Tolerance))
			{
				m_haParticlePosition[index]			= m_haParticlePosition[m_ParticleNum - 1];
				m_haParticleVelocity[index]			= m_haParticleVelocity[m_ParticleNum - 1];

				m_haParticlePressure[index]			= m_haParticlePressure[m_ParticleNum - 1];
				m_haParticleDensity[index]			= m_haParticleDensity[m_ParticleNum - 1];
				m_haParticleTemperature[index]		= m_haParticleTemperature[m_ParticleNum - 1];
				m_haParticleKineticViscosity[index]	= m_haParticleKineticViscosity[m_ParticleNum - 1];
				m_haParticleSolidPhaseRate[index]	= m_haParticleSolidPhaseRate[m_ParticleNum - 1];
				m_haParticleType[index]				= m_haParticleType[m_ParticleNum - 1];
				m_haParticleID[index]				= m_haParticleID[m_ParticleNum - 1];
				m_ParticleNum--;
			}
			else
			{
				continue;
			}
		}
		Status = ParticleHostToDevice();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}

CCTStatusType CCellSet::CalcExplicitly()
{
	CCTStatusType Status = CCT_NOERR;
	//Set Device ID
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	if(m_ParticleNum > 0)
	{
			Status = ::CalcExplicitly(m_CudaStream,m_ParticleNum);
			CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CCellSet::MoveTriangles()
{
	CCTStatusType Status = CCT_NOERR;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	Status = ::UpdateTrianglePosition(m_CudaStream,m_TriangleNum,m_daTriangles);
	CCT_ERROR_CHECK(Status);

//	//Display Triangle
//#if 1
//	cudaMemcpy(m_haTriangles,m_daTriangles, m_TriangleNum * sizeof(CTriangle), cudaMemcpyDeviceToHost);
//	std::stringstream filename;
//	filename<<"test\\TriangleInfo.dat";
//	std::ofstream of((filename.str().c_str()));
//	if(of)
//	{		
//		of << m_TriangleNum * 3<<"\n";
//		for(int i = 0 ; i < m_ParticleNum ; i++)
//		{
//			of<< m_haTriangles[i].Vertex[0].x <<" " << m_haTriangles[i].Vertex[0].y <<" " << m_haTriangles[i].Vertex[0].z << " 0 0 0 0 0 0 0 0 0"<< "\n";
//			of<< m_haTriangles[i].Vertex[1].x <<" " << m_haTriangles[i].Vertex[1].y <<" " << m_haTriangles[i].Vertex[1].z << " 0 0 0 0 0 0 0 0 0"<< "\n";
//			of<< m_haTriangles[i].Vertex[2].x <<" " << m_haTriangles[i].Vertex[2].y <<" " << m_haTriangles[i].Vertex[2].z << " 0 0 0 0 0 0 0 0 0"<< "\n";
//		}	
//		of<<"0\n0\n0 0 0 0";
//	}
//	of.close();
//	
//#endif
	return CCT_NOERR;
}
CCTStatusType CCellSet::RotateTrianglePosition(const Integer analysisStep)
{
	CCTStatusType Status = CCT_NOERR;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	Status = ::RotateTrianglePosition(m_CudaStream,m_TriangleNum,m_daTriangles,analysisStep);
	CCT_ERROR_CHECK(Status);

	return CCT_NOERR;
}
void CCellSet::AddToBuffer( bool Output, Integer &BufferSize, Integer * ParticleBufferID,
								 Scalar3* ParticleBufferPosition,
								 Scalar3* ParticleBufferVelcity, Scalar* ParticleBufferPressure,
								 Scalar* ParticleBufferDensity, Scalar* ParticleBufferTemperature, Scalar* ParticleBufferKineticViscosity,
								 Scalar* ParticleBufferSolidPhaseRate, ParticleType* ParticleBufferType)
{
	//for(int i = m_OutputParticleStart ; i < m_OutputParticleStart + m_ComputeParticleNum; ++i)
	BufferSize = m_ParticleNum;
	for(int i = 0 ; i < m_ParticleNum ; ++i)
	{
		int id = m_haParticleID[i];
		//ParticleBufferID[id] = id;

		if(!IsInclude(&m_BoundingBox,&m_haParticlePosition[i],m_Parameter.Tolerance) || id < 0)
		{
			if(Output)
			{
				ParticleBufferID[i] = -1;

				ParticleBufferPosition[i]			= make_Scalar3(0,0,0);
				ParticleBufferVelcity[i]			= make_Scalar3(0,0,0);
				ParticleBufferPressure[i]			= 0;
				ParticleBufferDensity[i]			= 0;
				ParticleBufferTemperature[i]		= 0;
				ParticleBufferKineticViscosity[i]	= 0;
				ParticleBufferSolidPhaseRate[i]		= 0;
				ParticleBufferType[i]				= TYPE_HEAVY_WEIGHT;
			}
			else
			{
				ParticleBufferID[i] = -1;

				ParticleBufferPosition[i]			= m_haParticlePosition[i];
				ParticleBufferVelcity[i]			= m_haParticleVelocity[i];
				ParticleBufferPressure[i]			= m_haParticlePressure[i];
				ParticleBufferDensity[i]			= m_haParticleDensity[i];
				ParticleBufferTemperature[i]		= m_haParticleTemperature[i];
				ParticleBufferKineticViscosity[i]	= m_haParticleKineticViscosity[i];
				ParticleBufferSolidPhaseRate[i]		= m_haParticleSolidPhaseRate[i];
				ParticleBufferType[i]				= m_haParticleType[i];
			}
		}
		else
		{
			ParticleBufferID[i] = id;

			ParticleBufferPosition[i]			= m_haParticlePosition[i];
			ParticleBufferVelcity[i]			= m_haParticleVelocity[i];
			ParticleBufferPressure[i]			= m_haParticlePressure[i];
			ParticleBufferDensity[i]			= m_haParticleDensity[i];
			ParticleBufferTemperature[i]		= m_haParticleTemperature[i];
			ParticleBufferKineticViscosity[i]	= m_haParticleKineticViscosity[i];
			ParticleBufferSolidPhaseRate[i]		= m_haParticleSolidPhaseRate[i]; //id;
			ParticleBufferType[i]				= m_haParticleType[i];
			
		}
#if 0 
		
		if(id < 0)
		{
			continue;
		}
		if(!IsInclude(&m_BoundingBox,&m_haParticlePosition[i],m_Parameter.Tolerance))
		{
			if(Output)
			{
				ParticleBufferID[id] = -1;

				ParticleBufferPosition[id]			= make_Scalar3(0,0,0);
				ParticleBufferVelcity[id]			= make_Scalar3(0,0,0);
				ParticleBufferPressure[id]			= 0;
				ParticleBufferDensity[id]			= 0;
				ParticleBufferTemperature[id]		= 0;
				ParticleBufferKineticViscosity[id]	= 0;
				ParticleBufferSolidPhaseRate[id]	= 0;
				ParticleBufferType[id]				= TYPE_HEAVY_WEIGHT;
			}
			else
			{
				ParticleBufferID[id] = -1;

				ParticleBufferPosition[id]			= m_haParticlePosition[i];
				ParticleBufferVelcity[id]			= m_haParticleVelocity[i];
				ParticleBufferPressure[id]			= m_haParticlePressure[i];
				ParticleBufferDensity[id]			= m_haParticleDensity[i];
				//ParticleBufferTemperature[id]		= m_haParticleTemperature[i];
				ParticleBufferTemperature[id]	= ParticleBufferID[id]; // Test;
				ParticleBufferKineticViscosity[id]	= m_haParticleKineticViscosity[i];
				ParticleBufferSolidPhaseRate[id]	= m_haParticleSolidPhaseRate[i]; //id;
				ParticleBufferType[id]				= ELASTICPARTICLE; // m_haParticleType[i]; Test
			}
		}
		else
		{
			ParticleBufferID[id] = id;

			ParticleBufferPosition[id]			= m_haParticlePosition[i];
			ParticleBufferVelcity[id]			= m_haParticleVelocity[i];
			ParticleBufferPressure[id]			= m_haParticlePressure[i];
			ParticleBufferDensity[id]			= m_haParticleDensity[i];
			ParticleBufferTemperature[id]		= m_haParticleTemperature[i];
			ParticleBufferKineticViscosity[id]	= m_haParticleKineticViscosity[i];
			ParticleBufferSolidPhaseRate[id]	= m_haParticleSolidPhaseRate[i]; //id;
			ParticleBufferType[id]				= m_haParticleType[i];

			// Test
			if(Output)
			{
				if(m_haParticlePosition[i].x > 1790)
				{
					ParticleBufferTemperature[id] = -1000;
				}
			}
		}
#endif

#if 0
		if(Output && (!IsInclude(&m_BoundingBox,&m_haParticlePosition[i],m_Parameter.Tolerance) || id < 0))
		{
			ParticleBufferID[id] = -1;

			ParticleBufferPosition[id]			= make_Scalar3(0,0,0);
			ParticleBufferVelcity[id]			= make_Scalar3(0,0,0);
			ParticleBufferPressure[id]			= 0;
			ParticleBufferDensity[id]			= 0;
			ParticleBufferTemperature[id]		= 0;
			ParticleBufferKineticViscosity[id]	= 0;
			ParticleBufferSolidPhaseRate[id]	= 0;
			ParticleBufferType[id]				= TYPE_HEAVY_WEIGHT;
		}
		else if((!IsInclude(&m_BoundingBox,&m_haParticlePosition[i],m_Parameter.Tolerance) || id < 0))
		{
			ParticleBufferID[id] = -1;

			ParticleBufferPosition[id]			= m_haParticlePosition[i];
			ParticleBufferVelcity[id]			= m_haParticleVelocity[i];
			ParticleBufferPressure[id]			= m_haParticlePressure[i];
			ParticleBufferDensity[id]			= m_haParticleDensity[i];
			ParticleBufferTemperature[id]		= m_haParticleTemperature[i];
			ParticleBufferKineticViscosity[id]	= m_haParticleKineticViscosity[i];
			ParticleBufferSolidPhaseRate[id]	= m_haParticleSolidPhaseRate[i]; //id;
			ParticleBufferType[id]				= ELASTICPARTICLE; // m_haParticleType[i]; Test
		}
		else	
		{
			ParticleBufferID[id] = id;

			ParticleBufferPosition[id]			= m_haParticlePosition[i];
			ParticleBufferVelcity[id]			= m_haParticleVelocity[i];
			ParticleBufferPressure[id]			= m_haParticlePressure[i];
			ParticleBufferDensity[id]			= m_haParticleDensity[i];
			ParticleBufferTemperature[id]		= m_haParticleTemperature[i];
			ParticleBufferKineticViscosity[id]	= m_haParticleKineticViscosity[i];
			ParticleBufferSolidPhaseRate[id]	= m_haParticleSolidPhaseRate[i]; //id;
			ParticleBufferType[id]				= m_haParticleType[i];
		}
#endif
	}
}
CCTStatusType CCellSet::ParticleNumberToContantMemory()
{
	CCTStatusType Status = CCT_NOERR;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);

	Status = ::ParticleNumberToConst(m_ParticleNum);
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
}
CCTStatusType CCellSet::StreamSynchronize()
{
	CCTStatusType Status = CCT_NOERR;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
//	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	Status = CudaSafeCall(cudaStreamSynchronize(m_CudaStream));
	CCT_ERROR_CHECK(Status);
	return Status;
}
void CCellSet::InitializeParameter(CParameter* pParameter, CParameterCoefficients* pParameterCoeff)
{
	m_Parameter = *pParameter;
	m_ParameterCoefficient = *pParameterCoeff;
}
CCTStatusType CCellSet::CalcExplicitPressure()
{
	CCTStatusType Status = CCT_NOERR;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	if(m_ParticleNum >0 )
	{
		Status = ::CalcExplicitPressure(m_CudaStream, m_ParticleNum);
		CCT_ERROR_CHECK(Status);
	}
	return Status;
}
CCTStatusType CCellSet::CalculateDragEffect()
{
	CCTStatusType Status = CCT_NOERR;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	if(m_ParticleNum >0 )
	{			
		Status = ::CalcDragEffect(m_CudaStream,m_ParticleNum);
		CCT_ERROR_CHECK(Status);		
	}
	return Status;
}
CCTStatusType CCellSet::CalcExplicitPressureGradient()
{
	CCTStatusType Status = CCT_NOERR;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	if(m_ParticleNum >0 )
	{			
		Status = ::CalcExplicitPressureGradient(m_CudaStream,m_ParticleNum);
		CCT_ERROR_CHECK(Status);		
	}
	return Status;
}
CCTStatusType CCellSet::CalcTemperatureFactor()
{
	CCTStatusType Status = CCT_NOERR;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	if(m_ParticleNum > 0)
	{
		Status = ::CalcTemperatureFactor(m_CudaStream,m_ParticleNum);
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CCellSet::CalcTurbulenceViscosity()
{
	CCTStatusType Status = CCT_NOERR;
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	if(m_ParticleNum > 0)
	{
		Status = ::CalcTurbulenceViscosity(m_CudaStream,m_ParticleNum);
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CCellSet::SetCudaDevice()
{
	CCTStatusType Status = CCT_NOERR;
//	Status = CudaSafeCall(cudaSetDevice(m_ID));
	Status = CudaSafeCall(cudaSetDevice(m_DeviceID));
	CCT_ERROR_CHECK(Status);
	return Status;
}

//void CCellSet::RegisterTriangle(const int TriangleNumber,			const CTriangle * Triangle,
//								const CGridBox *const BoundingBox,	const CGridParams* const GridParams,
//								Integer ID,							CCell* daCell,
//								Integer CellNum)
//{
//	for(int i= 0; i < TriangleNumber;++i)
//	{
//		Scalar3 A = Triangle[i].Vertex[0]; // Vert0;
//		Scalar3 B = Triangle[i].Vertex[1]; // Vert1;
//		Scalar3 C = Triangle[i].Vertex[2]; // Vert2;
//		Scalar3 MidPoint = make_Scalar3((A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3, (A.z + B.z + C.z) / 3);
//		const Scalar CellSize = BoundingBox->m_CellSize ;
//		Scalar MagBC; 
//		Scalar MagAB; 
//		Scalar MagCA; 
//		do{
//			Scalar3 AB = B - A;
//			Scalar3 BA = A - B;
//			Scalar3 BC = C - B;
//			Scalar3 CB = B - C;
//			Scalar3 CA = A - C;
//			Scalar3 AC = C - A;
//
//			MagBC = Magnitude(&BC);
//			MagAB = Magnitude(&AB);
//			MagCA = Magnitude(&CA);
//
//			if(MagAB > 0.0)
//			{
//				RegisterLine(1,BoundingBox,GridParams,A, make_Scalar3(AB.x / MagAB, AB.y / MagAB, AB.z / MagAB), MagAB, ID, daCell, CellNum);
//			}		
//			if(MagBC > 0.0)
//			{
//				RegisterLine(1,BoundingBox,GridParams,B, make_Scalar3(BC.x / MagBC, BC.y / MagBC, BC.z / MagBC), MagBC,ID, daCell, CellNum);
//			}		
//			if(MagCA > 0.0)
//			{
//				RegisterLine(1,BoundingBox,GridParams,C, make_Scalar3(CA.x / MagCA, CA.y / MagCA, CA.z / MagCA) , MagCA,  ID , daCell, CellNum);
//			}		
//			Scalar3 AM = AddVector(&AB, &AC);
//			Scalar3 BM = AddVector(&BA, &BC);
//			Scalar3 CM = AddVector(&CA, &CB);
//			Scalar MagAM = Magnitude(&AM);
//			Scalar MagBM = Magnitude(&BM);
//			Scalar MagCM = Magnitude(&CM);
//			if(MagAM > 0.0)
//			{
//				AM.x /= MagAM; AM.y /= MagAM; AM.z /= MagAM;
//				A.x += AM.x * CellSize;
//				A.y += AM.y * CellSize;
//				A.z += AM.z * CellSize;
//			}
//			if(MagBM > 0.0)
//			{
//				BM.x /= MagBM; BM.y /= MagBM; BM.z /= MagBM;
//				B.x += BM.x * CellSize;
//				B.y += BM.y * CellSize;
//				B.z += BM.z * CellSize;
//			}
//			if(MagCM > 0.0)
//			{
//				CM.x /= MagCM; CM.y /= MagCM; CM.z /= MagCM;
//				C.x += CM.x * CellSize;
//				C.y += CM.y * CellSize;
//				C.z += CM.z * CellSize;
//			}		
//
//		}while(MagBC >= CellSize && MagAB >= CellSize && MagCA >= CellSize);
//	}
//}
//

CCTStatusType CCellSet::SetOutPath(std::string& OutputPath )
{
	m_Outpath = OutputPath;
	return CCT_NOERR;
}