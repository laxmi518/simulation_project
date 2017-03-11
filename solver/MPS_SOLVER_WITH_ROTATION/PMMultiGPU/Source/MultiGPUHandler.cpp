#include "stdafx.h"
#include "MultiGPUHandler.h"
#include "CGPUCellSet.h"
#include "Utility_Inline.h"
#include "IO.h"
#include "Configuration.h"
#include "OutputHandler.h"
#include "ETProfile.h"
#include "DtFrame.h"
#include "InnerPressureModel.h"
#include "OptionParser.h"

#define MULTIBUCKET

CMultiGPUHandler::CMultiGPUHandler(void)
:m_GridNum(0)
,m_aGrid(NULL)
,m_OutputParticles(NULL)
,m_MaxParticleNum(0)
,m_Output(NULL)
,m_pParameter(NULL)
,m_pConfiguration(NULL)
,m_haBucket(NULL)
,m_BucketNum(0)
,m_ParticleNum(0)
,m_pOutputFrame(NULL)
,m_aOutputFrame(NULL)
,m_OutputFrameNum(0)
,m_haDragTriangles(NULL)
,m_DragTriangleNum(0)
,m_StartStep(0)
,m_pDt(NULL)
,m_NextItr(0)
,m_pInnerPressrueModel(NULL)

//STARTS for individual Output particles
,m_aOutputParticleID(0)
,m_aOutputParticlePosition(0)
,m_aOutputParticleVelocity(0)
,m_aOutputParticlePressure(0)
,m_aOutputParticleDensity(0)
,m_aOutputParticleTemperature(0)
,m_aOutputParticleKineticViscosity(0)
,m_aOutputParticleSolidPhaseRate(0)
,m_aOutputParticleType(NULL)
//ENDS for individual Output particles

//MultiBucket Starts
,m_MultiBucket(NULL)
,m_NumberOfBuckets(0)
//MultiBucket Ends

,m_haTriangleParameters(NULL)
,m_TriangleModelNumber(0)
{
	m_ComputeZone.m_MaxBound = make_Scalar3(0.0,0.0,0.0);
	m_ComputeZone.m_MinBound = make_Scalar3(0.0,0.0,0.0);
}

CMultiGPUHandler::~CMultiGPUHandler(void)
{
	if(m_pConfiguration)
	{
		delete m_pConfiguration;
		m_pConfiguration = NULL;
	}
	if(m_aOutputFrame)
	{
		delete[] m_aOutputFrame;
		m_aOutputFrame = NULL;
	}
	if(m_haTriangles)
	{
		delete m_haTriangles;
		m_haTriangles = NULL;
	}
	if(m_haDragTriangles)
	{
		delete m_haDragTriangles;
		m_haDragTriangles = NULL;
	}	
	if(m_haBucket)
	{
		delete m_haBucket;
		m_haBucket = NULL;
	}
	if(m_pParameter)
	{
		delete m_pParameter;
		m_pParameter = NULL;
	}
	if(m_pDt)
	{
		delete m_pDt;
	}
	if(m_aGrid)
	{
		delete []m_aGrid;
		m_aGrid = NULL;
	}
	if(m_pInnerPressrueModel)
	{
		delete m_pInnerPressrueModel;
	}
}
void CMultiGPUHandler::SetComputeZone(CBox* ComputeZone)
{
	m_ComputeZone.m_MaxBound = ComputeZone->m_MaxBound;
	m_ComputeZone.m_MinBound = ComputeZone->m_MinBound;
}
void CMultiGPUHandler::SetComputeZone(CParticle* aParticle, Integer ParticleNum, CTriangle* aTriangle, Integer TriangleNum)
{
	m_ComputeZone.m_MaxBound.x = -FLT_MAX;
	m_ComputeZone.m_MaxBound.y = -FLT_MAX;
	m_ComputeZone.m_MaxBound.z = -FLT_MAX;

	m_ComputeZone.m_MinBound.x = FLT_MAX;
	m_ComputeZone.m_MinBound.y = FLT_MAX;
	m_ComputeZone.m_MinBound.z = FLT_MAX;
	if(aParticle)
	{
		for(Integer i = 0; i < ParticleNum; ++i)
		{
			Check(&aParticle[i]);
		}
	}
	if(aTriangle)
	{
		for(Integer i = 0; i < TriangleNum; ++i)
		{
			Check(&aTriangle[i]);
		}
	}
	Scalar MaxVal;
//#ifdef _DEBUG
//	MaxVal = std::max(m_pParameter->LaplacianInfluenceRadiusCoefficient,m_pParameter->GradientInfluenceRadiusCoefficient);
//#else
	MaxVal =  max(m_pParameter->LaplacianInfluenceRadiusCoefficient,m_pParameter->GradientInfluenceRadiusCoefficient);
//#endif

	Scalar OffsetValue = m_pParameter->InitialParticleDistance * MaxVal;

	m_ComputeZone.m_MaxBound.x += OffsetValue;
	m_ComputeZone.m_MaxBound.y += OffsetValue;
	m_ComputeZone.m_MaxBound.z += OffsetValue;

	m_ComputeZone.m_MinBound.x -= OffsetValue;
	m_ComputeZone.m_MinBound.y -= OffsetValue;
	m_ComputeZone.m_MinBound.z -= OffsetValue;
}
// check to find the minumun and maximum bound of the triangle
void CMultiGPUHandler::Check(const CTriangle* pTriangle)
{
	for(Integer i = 0; i < 3; ++i)
	{
		if(pTriangle->Vertex[i].x <= m_ComputeZone.m_MinBound.x)
		{
			m_ComputeZone.m_MinBound.x = pTriangle->Vertex[i].x;
		}
		if(pTriangle->Vertex[i].y <= m_ComputeZone.m_MinBound.y)
		{
			m_ComputeZone.m_MinBound.y = pTriangle->Vertex[i].y;
		}
		if(pTriangle->Vertex[i].z <= m_ComputeZone.m_MinBound.z)
		{
			m_ComputeZone.m_MinBound.z = pTriangle->Vertex[i].z;
		}
		if(pTriangle->Vertex[i].x >= m_ComputeZone.m_MaxBound.x)
		{
			m_ComputeZone.m_MaxBound.x = pTriangle->Vertex[i].x;
		}
		if(pTriangle->Vertex[i].y >= m_ComputeZone.m_MaxBound.y)
		{
			m_ComputeZone.m_MaxBound.y = pTriangle->Vertex[i].y;
		}
		if(pTriangle->Vertex[i].z >= m_ComputeZone.m_MaxBound.z)
		{
			m_ComputeZone.m_MaxBound.z = pTriangle->Vertex[i].z;
		}
	}
}
// check to find the minumun and maximum bound of the particle
void CMultiGPUHandler::Check(const CParticle* const pParticle)
{	
	if(pParticle->Position.x <= m_ComputeZone.m_MinBound.x)
	{
		m_ComputeZone.m_MinBound.x = pParticle->Position.x;
	}
	if(pParticle->Position.y <= m_ComputeZone.m_MinBound.y)
	{
		m_ComputeZone.m_MinBound.y = pParticle->Position.y;
	}
	if(pParticle->Position.z <= m_ComputeZone.m_MinBound.z)
	{
		m_ComputeZone.m_MinBound.z = pParticle->Position.z;
	}
	if(pParticle->Position.x >= m_ComputeZone.m_MaxBound.x)
	{
		m_ComputeZone.m_MaxBound.x = pParticle->Position.x;
	}
	if(pParticle->Position.y >= m_ComputeZone.m_MaxBound.y)
	{
		m_ComputeZone.m_MaxBound.y = pParticle->Position.y;
	}
	if(pParticle->Position.z >= m_ComputeZone.m_MaxBound.z)
	{
		m_ComputeZone.m_MaxBound.z = pParticle->Position.z;
	}
}
CCTStatusType CMultiGPUHandler::Initialize(CParticle* aParticle,Integer ParticleNum)
{	
	if(!aParticle)
	{
		printf("ERR: Null parameter");
		return CCT_PARAMERR;
	}	
	//Scalar InfluenceRadius = pParameter->InitialParticleDistance * pParameter->LaplacianInfluenceRadiusCoefficient > pParameter->GradientInfluenceRadiusCoefficient ? pParameter->LaplacianInfluenceRadiusCoefficient : pParameter->GradientInfluenceRadiusCoefficient;
	CCTStatusType Status = CCT_NOERR;
	Status =  InitializeGrids();
	CCT_ERROR_CHECK(Status);
	InitializeParticleData(aParticle,ParticleNum);
	if(m_haBucket)
	{
		InitializeBucketData();
	}
	InitializeTriangleData();
	//Status = InitializeOutputMemory();
	Status = InitializeOutputBuffer();
	CCT_ERROR_CHECK(Status);

	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::InitializeGrids()
{
	Integer3 Section = m_pConfiguration->GetSection();
	m_GridNum = Section.x * Section.y * Section.z;
	m_aGrid = new CCellSet[m_GridNum];
	//m_aGrid = new CCellSet;
	if(!m_aGrid)
	{
		printf("ERR: STL file load");("ERR: Memory");
		return CCT_MEMERR;
	}
	InitializeComputeZone(Section, m_pParameter);
	InitializeNeighbour(Section);
	CCTStatusType Status;
	Status = InitializeHostMemory();
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
}
void CMultiGPUHandler::InitializeComputeZone(Integer3 Section, CParameter* pParameter)
{
	for(Integer i = 0; i < Section.x; ++i)
	{
		for(Integer j = 0; j < Section.y; ++j)
		{
			for(Integer k = 0; k < Section.z; ++k)
			{				
				Integer GridID = (Section.y * Section.z) * i + Section.z * j + k;
				m_aGrid[GridID].SetGridID(i,j,k);				
				CBox SectionComputeZone = GetSectionComputeZone(m_ComputeZone, Section,i,j,k);
				m_aGrid[GridID].SetComputeZone(&SectionComputeZone, pParameter,&m_ComputeZone);
			}
		}
	}
}
void CMultiGPUHandler::InitializeNeighbour(Integer3 Section)
{
	for(Integer i = 0;i < m_GridNum;++i)
	{
		m_aGrid[i].SetNeighborGrids(m_aGrid,Section);
	}
}
CCTStatusType CMultiGPUHandler::InitializeHostMemory()
{
	CCTStatusType Status;
	for(Integer i = 0; i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].InitializeHostMemory(m_MaxParticleNum);
		if(Status != CCT_NOERR)
		{
			for(Integer j = i; j >= 0; --j)
			{
				m_aGrid[j].UninitializeHostMemory();
				return Status;
			}
		}
	}
	return CCT_NOERR;
}
void CMultiGPUHandler::InitializeParticleData(CParticle* aParticle, Integer ParticleNum)
{
	m_ParticleNum = ParticleNum;
	for(Integer j = 0; j < m_GridNum;++j) 
	{
		m_aGrid[j].SetParticleData(aParticle, ParticleNum);
	}
}
void CMultiGPUHandler::InitializeBucketData()
{
	if(m_haBucket && m_BucketNum > 0)
	{
		for(Integer j = 0; j < m_GridNum;++j) 
		{
			m_aGrid[j].SetBucketData(m_haBucket, m_BucketNum);
		}
	}
}
void CMultiGPUHandler::InitializeTriangleData()
{
	for(Integer j = 0; j < m_GridNum;++j) 
	{
		m_aGrid[j].SetTriangles(m_haTriangles, m_TriangleNum);
		m_aGrid[j].SetTrianglesParameters(m_haTriangleParameters,m_TriangleModelNumber);

		//Drag @Rajan
		m_aGrid[j].SetDragParameters(m_haSTLDragPrameter, m_DragTriangleNum);
	}
}
void CMultiGPUHandler::InitializeDragTriangleData()
{
	for(Integer j = 0; j < m_GridNum;++j) 
	{
		m_aGrid[j].SetDragTriangles(m_haDragTriangles, m_DragTriangleNum,m_pInnerPressrueModel);
	}
}
void CMultiGPUHandler::InitializeStartStep()
{
	for(Integer j = 0; j < m_GridNum;++j) 
	{
		m_aGrid[j].SetStartStep(m_StartStep);
	}
}
CCTStatusType CMultiGPUHandler::CalculateNew()
{
	CCTStatusType Status = CCT_NOERR;
		
	Status = CudaSafeCall(cudaGetLastError());
	CCT_ERROR_CHECK(Status);

	Integer deviceID = m_DeviceSelector[0];

	std::cout<<"Max Number from CMultiGPUHandler is :- "<<m_MaxParticleNum<<std::endl;
	Status = m_aGrid->InitializeGPU(deviceID);
	CCT_ERROR_CHECK(Status);

	m_Output = CreateEvent(NULL, TRUE, FALSE, NULL);
	if(m_Output == NULL)
	{
		return CCT_ETCERR;
	}	
	std::string TimerPath = m_pConfiguration->GetOutputPath();
	TimerPath.append("Times.txt");
	CETProfile::setFilePath(TimerPath.c_str());
	CETProfile Ep;
	Integer FileNo = m_StartStep / m_pConfiguration->GetOutputInterval();
	Integer AddedTImes = 0;
	Integer WallStep = 0;
	Integer i = m_StartStep;
	Integer STLRotationSteps = 0;

	for(; i <= m_pParameter->AnalysisStep; ++i)
	{
		Status = CheckAndUpdateDt(i);
		CCT_ERROR_CHECK(Status);

		Status = m_aGrid->RegisterParticleTopology();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Register Particle\n";

		Status = m_aGrid->RegisterTriangleTopology();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Register Triangle\n";

		Status = m_aGrid->CalcSTLDistance();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : STL Distance\n";

		if(m_DragTriangleNum > 0)		
		{
			Status = m_aGrid->CalculateDragEffect();
			CCT_ERROR_CHECK(Status);
		}
		
		if( true == m_pParameter->bUpdateTemperature && 0 == i % m_pParameter->ThermalStep)
		{
			Status = m_aGrid->CalcTemperatureFactor();
			CCT_ERROR_CHECK(Status);
			std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Temperature Factor\n";
		}

		Status = m_aGrid->CalcExplicitly();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Calculate Explicit\n";

		Status = m_aGrid->RelocateAllParticleTopology();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Relocate Particle\n";

		Status = m_aGrid->RegisterParticleTopology();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Register Particle\n";

		//Explicit output
		if((i % m_pConfiguration->GetOutputInterval()) == 0 && m_pParameter->bExplicitOut)
		{
			Status = m_aGrid->ParticleDeviceToHost();
			WaitForOutputFrame();			
			
			m_pOutputFrame->SetModelPosition(i, m_pParameter->Dt, m_pConfiguration->GetModelVelocity(),m_pConfiguration->GetModelPosition(),
											m_pConfiguration->GetModelResetInterval(),	m_pConfiguration->GetModelResetToOriginal());

			Status = this->Output(FileNo, "Explicit");
			CCT_ERROR_CHECK(Status);
		}

		Status = m_aGrid->CalcSTLDistance();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : STL Distance\n";

		Status = m_aGrid->CalcPressureExplicit();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Pressure Explicit \n";

		Status = m_aGrid->CalcPressureExplicitGradient();
		CCT_ERROR_CHECK(Status);		
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Calculate Gradient\n";

		Status = m_aGrid->RelocateAllParticleTopology();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Relocate Particle\n";
		
		//Setup STL Rotaion
		if( i!=0 )
		{
			m_pOutputFrame->SetRotationAngle(i,m_pConfiguration->GetRotationAngle(), m_pConfiguration->GetRotationStartTime(),m_pConfiguration->GetTimeToReachRPM() ,m_pParameter);
		}
		//Output Steps Starts
		if( 0 == (i % m_pConfiguration->GetOutputInterval()))
		{
			Status = this->Synchronize();
			CCT_ERROR_CHECK(Status);

			Status = this->ParticleDeviceToHost();
			CCT_ERROR_CHECK(Status);

			WaitForOutputFrame();			
			
			m_pOutputFrame->SetModelPosition(i, m_pParameter->Dt, m_pConfiguration->GetModelVelocity(),m_pConfiguration->GetModelPosition(),
											m_pConfiguration->GetModelResetInterval(),	m_pConfiguration->GetModelResetToOriginal());
			
			Status = this->Synchronize();
			CCT_ERROR_CHECK(Status);

			Status = this->Output(FileNo, "Implicit");
			CCT_ERROR_CHECK(Status);

			FileNo++;
		}
		else
		{
			Status = this->Synchronize();
			CCT_ERROR_CHECK(Status);
		}

		//Bucket Addition Starts
		bool IsBucketAdded = false;
		if(m_NumberOfBuckets > 0)
		{
			Status = this->BucketAddition(i,AddedTImes, IsBucketAdded);
			CCT_ERROR_CHECK(Status);		
		}
		if(IsBucketAdded)
		{
			Status = m_aGrid->RelocateAllParticleTopology();
			CCT_ERROR_CHECK(Status);
		}
		else
		{
			Status =  m_aGrid->UpdateTrianglePosition();
			CCT_ERROR_CHECK(Status);
			++WallStep;
			std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Update Triangle\n";
		}

		Status = m_aGrid->RotateTrianglePosition(i);
		CCT_ERROR_CHECK(Status);
		STLRotationSteps++;
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Rotate Trianlge\n";
				
		std::cout <<"**ATDev: "<< deviceID <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : REGISTER IMPLICIT :"<<i<<std::endl;
	}
	WaitForOutputThreads();
	Ep.SetAverageTime(i);
	Ep.DisplayTime(i);
	return Status;
}
CCTStatusType CMultiGPUHandler::ParticleDeviceToHost()
{
	CCTStatusType Status = CCT_NOERR;
	Status =  m_aGrid->ParticleDeviceToHost();
	CCT_ERROR_CHECK(Status);
	return Status;
}
CCTStatusType CMultiGPUHandler::WaitForOutput()
{
	CCTStatusType Status = CCT_NOERR;
	SetEvent(m_Output);	
	for(Integer i = 0; i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].WaitForOutput();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::AddToOutputBuffer(CParticle* ParticleBuffer, Integer BufferSize)
{
	for(Integer i = 0;i < BufferSize; ++ i)
	{
		Integer ID = ParticleBuffer[i].ID;
		if(ID < 0 || ID >= m_ParticleNum)
		{
			continue;
		}
		memcpy(&m_OutputParticles[ID], &ParticleBuffer[i], sizeof(CParticle));
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::WaitForOutputBuffer()
{
	WaitForSingleObject(m_Output,INFINITE);
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::Output(Integer Step, std::string OutputType)
{
	CCTStatusType Status = CCT_NOERR;	
	Status = this->GetDataOnMultiGPUHandler(true);
	CCT_ERROR_CHECK(Status);
	
	std::stringstream SaveTo;
	SaveTo << m_pConfiguration->GetOutputFileName() << "_" << OutputType << "_" << m_pConfiguration->GetOutputInterval() << "_" ;
	if(Step < 10)
	{
		SaveTo<<"0"<<Step<<".dat";
	}
	else
	{
		SaveTo<<Step<<".dat";
	}
	Scalar AnalysisTime = m_pParameter->Dt * m_pConfiguration->GetOutputInterval() * Step;
	m_pOutputFrame->SetOutputParams(m_ParticleNum, AnalysisTime,SaveTo.str());
	Status = COutputHandler::Output(m_pOutputFrame, m_pConfiguration->IsAscii());
	CCT_ERROR_CHECK(Status);
	std::cout<<"("<<Step<<") "<<OutputType<<" Out"<<std::endl;
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::GetDataOnMultiGPUHandler(bool Output)
{
	CCTStatusType Status = CCT_NOERR;	
	Status = m_aGrid->AddToOutputBuffer(Output, m_ParticleNum,		m_aOutputParticleID,	m_aOutputParticlePosition,
								m_aOutputParticleVelocity,m_aOutputParticlePressure,m_aOutputParticleDensity,
								m_aOutputParticleTemperature,m_aOutputParticleKineticViscosity,
								m_aOutputParticleSolidPhaseRate,m_aOutputParticleType);
	
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::ResetOutput()
{
	for(Integer i = 0; i < m_GridNum; ++i)
	{
		m_aGrid[i].ResetOutputEvent();
	}
	return CCT_NOERR;
}

CCTStatusType CMultiGPUHandler::Run(int argc, char** argv)
{
	
	//Device Selector Starts
	// parse command line arguments
	SystemUtility::COptionParser Parser;
	SystemUtility::PlaneArgs PlaneArg;
	SystemUtility::KeyValueArgs KeyValueArg;
	Parser.Parse(argc, argv, PlaneArg, KeyValueArg);
	//Device Selector Ends

	CCTStatusType Status = CCT_NOERR;
	if(argc < 2)
	{
		return CCT_PARAMERR;
	}
	int GPUNum;
	cudaGetDeviceCount((int*)&GPUNum);
	if(GPUNum < 1)
	{
		printf("ERR: No Suitable GPU");
		return CCT_CUDAERR;
	}

	if(!m_DeviceSelector.Initialize(GPUNum) || m_DeviceSelector.Empty())
	{
		printf("ERR: No Suitable GPU");
		return CCT_CUDAERR;
	}
	std::vector<int> DevicePriorityList;
	SystemUtility::CDevicePriorityListMaker PriorityListMaker(Parser);
	PriorityListMaker.Make(KeyValueArg, DevicePriorityList);
	
	//Display Version
	std::cout << PriorityListMaker.DisplayVersion(KeyValueArg);

	m_DeviceSelector.SetPriority(DevicePriorityList);


	m_pConfiguration = new CConfiguration();
	if(!m_pConfiguration)
	{
		printf("ERR: Out of Host Memory");
		return CCT_MEMERR;
	}
	//Status =  m_pConfiguration->ReadConfigurationFile(argv[1]);
	std::cout<<"INFO :- READ CONFIGURATION STARTS\n";
	Status = m_pConfiguration->ReadConfigurationFileWithDrag(argv[1]);
	CCT_ERROR_CHECK(Status);
	std::cout<<"INFO :- READ CONFIGURATION ENDS\n";

	Integer3 Section = m_pConfiguration->GetSection();
	if(Section.x * Section.y * Section.z == 0)
	{
		printf("ERR: Section Division Zero");
		return CCT_PARAMERR;
	}
	if(GPUNum < Section.x * Section.y * Section.z)
	{
		printf("ERR: Sections greater than no of GPU");
		return CCT_PARAMERR;
	}	

	std::cout<<"INFO :- LOAD PARAMETER FILE STARTS\n";
	Status = IO::Load(m_pConfiguration->GetParameterFile(), m_pParameter);
	CCT_ERROR_CHECK(Status);
	std::cout<<"INFO :- LOAD PARAMETER FILE ENDS\n";

	// Movable Wall model STL triangles load
	std::cout<<"INFO :- LOAD STL MODELS FILE STARTS\n";
	Status = IO::Load(m_pConfiguration->GetModelFiles(), m_TrianglePerModel,m_TriangleNum, m_haTriangles);
	CCT_ERROR_CHECK(Status);
	printf("Loaded %d Triangles\n",m_TriangleNum);
	std::cout<<"INFO :- LOAD STL MODELS FILE ENDS\n";

	//Initialize the TriangleParameters
	m_TriangleModelNumber = m_TrianglePerModel.size();
	if(m_TriangleModelNumber > 0)
	{
		m_haTriangleParameters = new CTriangleParameters[m_TriangleModelNumber];

		//m_haSTLDragPrameter		= new DragParameter[m_TriangleModelNumber];
	}

	//Drag@Rajan
	m_DragTriangleNum = m_pConfiguration->GetNumberOfDragSTL();
	if(m_DragTriangleNum > 0 )
	{
		m_haSTLDragPrameter = new DragParameter[m_DragTriangleNum];
	}
	for(unsigned int i = 0 ; i < m_DragTriangleNum ; i++)
	{
		m_haSTLDragPrameter[i] = m_pConfiguration->GetSTLDragParameter(i);
	}
	
	for(unsigned int i = 0; i < m_TrianglePerModel.size(); ++i)
	{
		Scalar3 Velocity	= m_pConfiguration->GetModelVelocity(i);
		Scalar Temperature	= m_pConfiguration->GetModelTemperature(i);

		std::cout <<"Temperature is : " << Temperature <<std::endl;

		Scalar3 Position	= m_pConfiguration->GetModelPosition(i);

		bool ResetToOriginal	= m_pConfiguration->GetModelResetToOriginal(i);
		Integer ResetInterval	= m_pConfiguration->GetModelResetInterval(i);

		//Paramter Setting for STL Rotation Starts

		bool RotatingOrNot			= m_pConfiguration->GetTrianlgeRotatingOrNot(i);
		//Scalar centerDegree		= m_pConfiguration->GetAngleDegree(i);
		Integer RotationInRPM		= m_pConfiguration->GetRotationRPM(i);
			
		//Calculating angle in Degree from Given RPM
		Scalar RotationInDegree  =	RotationInRPM * m_pParameter->Dt * 6;  // Angle of Rotation = (RPM /60) * dt * 360

		m_pConfiguration->SetRotationInAngle(RotationInDegree);

		Scalar StartTimeOfRotation		= m_pConfiguration->GetRotationStartTime(i);
		Scalar TimeToReachRPM			= m_pConfiguration->GetTimeToReachRPM(i);
		Scalar3 centerOfRotation		= m_pConfiguration->GetCenterOfRotation(i);
		Scalar3 secondPointOfRotation	= m_pConfiguration->GetOtherPointOfRotation(i);

		//For Wall Parameter in Configuration Starts
		Scalar WallFrictionCoefficient	= m_pConfiguration->GetWallFrictionCoefficient(i);
		Scalar WallThermalRestivity		= m_pConfiguration->GetWallThermalRestivity(i);
		Scalar WallThermalConductivity	= m_pConfiguration->GetWallThermalConductivity(i);
		//For Wall Parameter in Configuration Ends

		//Drag Test
		bool TestDrag = m_pConfiguration->GetModelISDrag(i);

		SetModelID(i);
		SetModelPosition(i,&Position);

		m_haTriangleParameters[i].Velocity = Velocity;
		m_haTriangleParameters[i].Temperature = Temperature;
		m_haTriangleParameters[i].LinearMovableOrNot = ResetToOriginal;
		m_haTriangleParameters[i].LinearResetInterval = ResetInterval;
		m_haTriangleParameters[i].RotatingOrNot = RotatingOrNot;
		m_haTriangleParameters[i].RotationStartTime = StartTimeOfRotation;
		m_haTriangleParameters[i].TimeToReachRPM = TimeToReachRPM;
		m_haTriangleParameters[i].CenterOfRotation = centerOfRotation;
		m_haTriangleParameters[i].SecondPointOfRotation = secondPointOfRotation;
		m_haTriangleParameters[i].AngleOfRotationInDegree = RotationInDegree;

		////For Wall Parameter in Configuration Starts
		m_haTriangleParameters[i].WallFrictionCoefficient = WallFrictionCoefficient;
		m_haTriangleParameters[i].WallThermalResistivity  = WallThermalRestivity;
		m_haTriangleParameters[i].WallThermalConductivity = WallThermalConductivity;

		std::cout <<"WallFriction		: "<<WallFrictionCoefficient<<"\n";
		std::cout <<"WallResistivity	: "<<WallThermalRestivity<<"\n";
		std::cout <<"WallConductivity	: "<<WallThermalConductivity<<"\n\n";

		//Drage @Rajan
		m_haTriangleParameters[i].isDrag = TestDrag;
		m_haTriangleParameters[i].DragTriangleID = m_pConfiguration->GetModelDragID(i);		
	}

	// (Reading of particles)
	CParticle* InitialParticles;
	Integer ParticleNum;
	std::cout<<"INFO :- LOAD PARTICLE FILE STARTS\n";
	Status = IO::Load(m_pConfiguration->GetParticleFile(), ParticleNum, InitialParticles);
	CCT_ERROR_CHECK(Status);
	printf("Loaded %d Input Particles\n",ParticleNum);
	std::cout<<"INFO :- LOAD PARTICLE FILE ENDS\n";
	
#ifdef MULTIBUCKET
	Status = LoadMultiBucket();
#else
	//starts Bucket Load E&T Nepal
	if(!m_pConfiguration->GetBucketFile().empty())
	{
		Status = IO::Load(m_pConfiguration->GetBucketFile(), m_BucketNum ,m_haBucket);
		CCT_ERROR_CHECK(Status);
		printf("Loaded %d Bucket Particles\n",m_BucketNum);
	}
#endif

	//Calculation of Maximum number of Buckets
#ifdef MULTIBUCKET
	int totalBucketNumber = 0;
	for(int i = 0 ; i < m_NumberOfBuckets; ++i)
	{
		totalBucketNumber += m_MultiBucket[i].TotalAdditionStep * m_MultiBucket[i].BucketNumber;
	}

	m_MaxParticleNum = ParticleNum + totalBucketNumber;
	std::cout <<"Maximum number of Particles is : "<< m_MaxParticleNum <<std::endl;
#else
	m_MaxParticleNum = ParticleNum + m_BucketNum * m_pParameter->TotalGenerationStep;
#endif

	//// Allocate Host Memory
	//m_haParticle = new CParticle[m_MaxParticleNum];
	//memcpy(m_haParticle, InitialParticles, m_ParticleNum * sizeof(CParticle));	
	//delete InitialParticles;

	CalculateLambdaValueLocal();

	CBox * ComputeZone = m_pConfiguration->GetComputeZone();
	if(ComputeZone)
	{
		this->SetComputeZone(ComputeZone);
	}
	else
	{
		this->SetComputeZone(InitialParticles, ParticleNum, m_haTriangles, m_TriangleNum);
	}
	Status = Initialize(InitialParticles,ParticleNum);
	CCT_ERROR_CHECK(Status);

	delete InitialParticles;

	//if(argc >= 3)
	if(PlaneArg.size() >= 3)
	{
		std::string FileName = argv[2];
		if(!FileName.empty())
		{
			std::string InnerPressureParameterFile;
			Status = IO::Load(FileName,m_haDragTriangles,m_DragTriangleNum,InnerPressureParameterFile);
			CCT_ERROR_CHECK(Status);
			std::cout<<"Loaded :"<< m_DragTriangleNum << "Drag Triangles"<<std::endl;

			m_pInnerPressrueModel = new CInnerPressureModel(m_pParameter);
			Status = m_pInnerPressrueModel->Load(InnerPressureParameterFile);
			CCT_ERROR_CHECK(Status);
			m_pInnerPressrueModel->SetOutputPath(m_pConfiguration->GetOutputPath());

			InitializeDragTriangleData();
		}
	}
	//if(argc >= 4)
	if(PlaneArg.size() >= 4)
	{
		std::string FileName = argv[3];
		if(!FileName.empty())
		{
			m_pDt = new DtFrame();
			Status =  m_pDt->Load(FileName);
			if(Status != CCT_NOERR)
			{
				delete m_pDt;
				m_pDt = NULL;
				return Status;
			}
		}
	}
	m_StartStep = m_pConfiguration->GetStartStep();
	InitializeStartStep();	
	OutputComputeZone();

	//Save The used Parameter File After LambdaValue Calcuation
	std::stringstream UsedParameterFileName;
	UsedParameterFileName << m_pConfiguration->GetOutputPath() << "\\UsedParameters.txt";

	IO::saveParameterWithName(UsedParameterFileName.str(),				m_pParameter,								m_TriangleModelNumber,
								m_pConfiguration->GetCenterOfRotation(),m_pConfiguration->GetOtherPointOfRotation(),m_pConfiguration->GetRotationRPM(),
								m_pConfiguration->GetWallFrictionCoefficient(),m_pConfiguration->GetWallThermalConductivity(),m_pConfiguration->GetWallThermalRestivity(),								
								m_haSTLDragPrameter,m_DragTriangleNum);
	//Save The Used Parameter File Ends


	//Status = this->Calculate();
	Status = this->CalculateNew();
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
}
void CMultiGPUHandler::SetModelID(const Integer Model)
{
	Integer Start = 0;
	if(Model == 0)
	{
		Start = 0;
	}
	else
	{
		//Start = m_TrianglePerModel[Model-1];
		for(int j=0;j<Model;++j)
		{
			Start += m_TrianglePerModel[j];
		}
	}
	Integer TriangleNum = m_TrianglePerModel[Model];
	for(Integer i = Start; i < Start + TriangleNum; ++i)
	{
		m_haTriangles[i].ModelIDNumber = Model;
	}
}
void CMultiGPUHandler::SetModelPosition(Integer Model, const Scalar3* const Position )
{
	Integer Start = 0;
	if(Model == 0)
	{
		Start = 0;
	}
	else
	{
		for(int j=0;j<Model;++j)
		{
			Start += m_TrianglePerModel[j];
		}
	}
	Integer TriangleNum = m_TrianglePerModel[Model];
	for(Integer i = Start; i < Start + TriangleNum; ++i)
	{
		m_haTriangles[i].Vertex[0] = m_haTriangles[i].Vertex[0] + *Position;
		m_haTriangles[i].Vertex[1] = m_haTriangles[i].Vertex[1] + *Position;
		m_haTriangles[i].Vertex[2] = m_haTriangles[i].Vertex[2] + *Position;
	}
}
void CMultiGPUHandler::WaitForOutputFrame()
{
	Integer i = 0;
	while(true)
	{
		m_pOutputFrame = m_aOutputFrame[i].GetOutputFrame();
		if(m_pOutputFrame)
		{
			m_aOutputParticleID					= m_pOutputFrame->m_aOutputParticleID;
			m_aOutputParticlePosition			= m_pOutputFrame->m_aOutputParticlePosition;
			m_aOutputParticleVelocity			= m_pOutputFrame->m_aOutputParticleVelocity;
			m_aOutputParticlePressure			= m_pOutputFrame->m_aOutputParticlePressure;
			m_aOutputParticleDensity			= m_pOutputFrame->m_aOutputParticleDensity;
			m_aOutputParticleTemperature		= m_pOutputFrame->m_aOutputParticleTemperature;
			m_aOutputParticleKineticViscosity	= m_pOutputFrame->m_aOutputParticleKineticViscosity;
			m_aOutputParticleSolidPhaseRate		= m_pOutputFrame->m_aOutputParticleSolidPhaseRate;
			m_aOutputParticleType				= m_pOutputFrame->m_aOutputParticleType;
			break;
		}
		++i;
		if(i >= m_OutputFrameNum)
		{
			i = 0;
		}
	}
}
CCTStatusType CMultiGPUHandler::InitializeOutputBuffer()
{
	SYSTEM_INFO SI;
	GetSystemInfo(&SI);
	//m_OutputFrameNum = SI.dwNumberOfProcessors;
	m_OutputFrameNum = 1;
	m_aOutputFrame = new COutputFrame[m_OutputFrameNum];
	if(!m_aOutputFrame)
	{
		return CCT_ETCERR;
	}
	for(Integer i = 0 ;i < m_OutputFrameNum; ++i)
	{
		CCTStatusType Status =  m_aOutputFrame[i].CreateOutputFrame(m_MaxParticleNum,m_pConfiguration->GetModelVelocity().size());
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::WaitForOutputThreads()
{
	for(Integer i = 0; i< m_OutputFrameNum; ++i)
	{
		WaitForSingleObject(m_aOutputFrame[i].m_ThreadHandle,INFINITE);
		CloseHandle(m_aOutputFrame[i].m_ThreadHandle);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::AddBucketParticles(const CParticle * bucket, const Integer bucketNumber)
{
	InsertBucket(bucket,bucketNumber);
	CCTStatusType Status = SetAddedBucket();
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
}
void CMultiGPUHandler::InsertBucket(const CParticle * bucket, const Integer bucketNumber)
{
	Integer j = 0;
	Integer i = 0;
	for(; i < m_ParticleNum; ++i)
	{
		if(m_aOutputParticleID[i] < 0)
		{
			m_aOutputParticleID[i]					= i;
			m_aOutputParticlePosition[i]			= bucket[j].Position;
			m_aOutputParticleVelocity[i]			= bucket[j].Velocity;
			m_aOutputParticleTemperature[i]			= bucket[j].Temperature;
			m_aOutputParticleDensity[i]				= bucket[j].Density;
			m_aOutputParticleKineticViscosity[i]	= bucket[j].KineticViscosity;
			m_aOutputParticleSolidPhaseRate[i]		= bucket[j].SolidPhaseRate;
			m_aOutputParticleType[i]				= bucket[j].Type;

			j++;

			if(j >= bucketNumber )
			{
				return;
			}
		}
	}
	if(m_ParticleNum + bucketNumber - j <= m_MaxParticleNum )
	{
		for(; j < bucketNumber ; ++j,++i)
		{
			m_aOutputParticleID[i]					= i;
			m_aOutputParticlePosition[i]			= bucket[j].Position;
			m_aOutputParticleVelocity[i]			= bucket[j].Velocity;
			m_aOutputParticleTemperature[i]			= bucket[j].Temperature;
			m_aOutputParticleDensity[i]				= bucket[j].Density;
			m_aOutputParticleKineticViscosity[i]	= bucket[j].KineticViscosity;
			m_aOutputParticleSolidPhaseRate[i]		= bucket[j].SolidPhaseRate;
			m_aOutputParticleType[i]				= bucket[j].Type;
		}
		m_ParticleNum = i;
	}
	return;
}
CCTStatusType CMultiGPUHandler::SetAddedBucket()
{
	CCTStatusType Status = CCT_NOERR;
	for(Integer k = 0;k < m_GridNum; ++k)
	{
		Status = m_aGrid[k].SetParticleData(m_ParticleNum,m_aOutputParticleID,m_aOutputParticlePosition,
											m_aOutputParticleVelocity,m_aOutputParticlePressure,
											m_aOutputParticleDensity,m_aOutputParticleTemperature,
											m_aOutputParticleKineticViscosity,m_aOutputParticleSolidPhaseRate,
											m_aOutputParticleType);
		CCT_ERROR_CHECK(Status);
		Status = m_aGrid[k].SetAddedParticlesEvent();
		CCT_ERROR_CHECK(Status);
	}
	m_pOutputFrame->SetEmpty(true);
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::Synchronize()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ; i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].StreamSynchronize();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;

	/*
	CCTStatusType Status = CCT_NOERR;
	WaitForOutputFrame();	
	for(Integer i = 0;i < m_GridNum; ++i)
	{
		m_aGrid[i].SetOutputBuffer();
	}	
	Status =  WaitForOutput();
	CCT_ERROR_CHECK(Status);	
	
	return CCT_NOERR;
	*/
}
void CMultiGPUHandler::OutputComputeZone()
{
	std::string Path = m_pConfiguration->GetOutputPath();
	Path.append("ComputeZone.txt");
	std::ofstream file(Path.c_str());
	if(file.is_open())
	{
		file << m_ComputeZone.m_MinBound.x << " " << m_ComputeZone.m_MinBound.y << " " << m_ComputeZone.m_MinBound.z<< " ";
		file << m_ComputeZone.m_MaxBound.x << " " << m_ComputeZone.m_MaxBound.y << " " << m_ComputeZone.m_MaxBound.z<< " ";
	}
}
const DtFrame* CMultiGPUHandler::GetDtFrame() const
{
	return m_pDt;
}

CCTStatusType CMultiGPUHandler::CheckAndUpdateDt(Integer CurrentStep)
{
	if(m_pDt)
	{
		Integer DtStep = m_pDt->GetStep(m_NextItr);
		if((DtStep >= 0) && CurrentStep >= DtStep)
		{
			const DeltaTime* DtTime = m_pDt->GetDtTime(m_NextItr);
			if(DtTime)
			{
				m_pParameter->Dt = DtTime->Dt;
				m_pParameter->CourantNumber = DtTime->CourantNumber;
				m_pParameter->ParticleGenerationStep = DtTime->InflowStep;
				m_pParameter->OutputStep = DtTime->OutputStep;
				++m_NextItr;
			}
		}
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::LoadMultiBucket()
{
	CCTStatusType Status = CCT_NOERR;
	if(!m_pConfiguration->GetMultiBucketFiles().empty())
	{
		m_NumberOfBuckets = m_pConfiguration->GetMultiBucketFiles().size();
		m_MultiBucket = new MultipleBucket[m_NumberOfBuckets];
		for(int j = 0 ; j < m_NumberOfBuckets ; ++j)
		{
			Status = IO::Load(m_pConfiguration->GetMultiBucketFiles()[j], m_MultiBucket[j].BucketNumber ,m_MultiBucket[j].Bucket);
			m_MultiBucket[j].AdditionInterval	= m_pConfiguration->GetMultiAdditionInterval()[j];
			m_MultiBucket[j].TotalAdditionStep	= m_pConfiguration->GetMultiTotalAdditionStep()[j];
			m_MultiBucket[j].AddedTimes = 0;
		}
	}
	return Status;
}
CCTStatusType CMultiGPUHandler::BucketAddition(Integer AnalysisStep, Integer & AddedTImes, bool & IsBucketAdded)
{
	CCTStatusType Status = CCT_NOERR;
	for(Integer p = 0 ; p < m_NumberOfBuckets ; ++p)
	{
		if((m_MultiBucket[p].AddedTimes < m_MultiBucket[p].TotalAdditionStep) && (0 == ( (AnalysisStep + 1) % m_MultiBucket[p].AdditionInterval)))
		{
			IsBucketAdded = true;

			Status = this->Synchronize();
			CCT_ERROR_CHECK(Status);

			Integer OldNum = m_ParticleNum;	
			//std::cout<< "::::::::::::::::Bucket Addition Particle Device to Host\n";
			Status = m_aGrid->ParticleDeviceToHost();
			CCT_ERROR_CHECK(Status);

			WaitForOutputFrame();

			//std::cout<< "::::::::::::::::Bucket Addition Get Data on MultiGPUHaldler\n";
			Status = this->GetDataOnMultiGPUHandler(false);
			CCT_ERROR_CHECK(Status);

			//std::cout<< "::::::::::::::::Bucket Addition Add Bucket Particles\n";
			Status = AddBucketParticles(m_MultiBucket[p].Bucket,m_MultiBucket[p].BucketNumber);
			CCT_ERROR_CHECK(Status);

			//std::cout<< "::::::::::::::::Bucket Addition Reset Wall Position\n";
			Status = m_aGrid->ResetWallPosition(AnalysisStep);
			CCT_ERROR_CHECK(Status);
			
			//std::cout<< "::::::::::::::::Bucket Addition Particle Host to Device\n";
			Status = m_aGrid->ParticleHostToDevice(); //Need to copy into output Particle
			CCT_ERROR_CHECK(Status);

			//std::cout<< "::::::::::::::::Bucket Addition particle Number to Const Memory\n";
			Status = m_aGrid->ParticleNumberToContantMemory();
			CCT_ERROR_CHECK(Status);

			++m_MultiBucket[p].AddedTimes;
			++AddedTImes;

			std::cout<<"** ATDevice : "<< "New Particles("<<OldNum<<" To "<< m_ParticleNum<<" ) Added"<< AddedTImes << " Times" <<std::endl <<std::endl;
		}
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::CalculateLambdaValueLocal()
{	
	const Scalar I_0 = m_pParameter->InitialParticleDistance;
	Scalar a = 0.0, b = 0.0;
		
	if(0<I_0<=0.0001)
	{
		a = 5.02608E-06;
		b = 3.82946E-10;
	}
	
	else if(0.0001<I_0<=0.0005)
	{
		a = 0.000124918;
		b = 1.81984E-08;
	}
	
	else if(0.0005<I_0<=0.001)
	{
		a = 0.0017886;
		b = -4.84707E-07;
	}
	
	else if(0.001<I_0<=0.002)
	{
		a = 0.00373277;
		b = -2.15486E-06;
	}
	
	else if(0.002<I_0<=0.003)
	{
		a = 0.0066746;
		b = -6.72384E-06;
	}
	
	else if(0.003<I_0<=0.004)
	{
		a = 0.0107452;
		b = -1.76386E-05;
	}
	
	else if(0.004<I_0<=0.005)
	{
		a = 0.0141993;
		b = -3.08052E-05;
	}
	
	else if(0.005<I_0<=0.006)
	{
		a = 0.0164975;
		b = -4.14937E-05;
	}
	
	else if(0.006<I_0<=0.007)
	{
		a = 0.0181782;
		b = -4.95251E-05;
	}
	
	else if(0.007<I_0<=0.008)
	{
		a = 0.0205885;
		b = -6.28759E-05;
	}
	
	else if(0.008<I_0<=0.009)
	{
		a = 0.0222761;
		b = -7.17244E-05;
	}
	
	else if(0.009<I_0<=0.01)
	{
		a = 0.0231972;
		b = -7.34117E-05;
	}
	
	else if(0.01<I_0<=0.02)
	{
		a = 0.0369119;
		b = -0.000210081;
	}
	
	else if(0.02<I_0<=0.03)
	{
		a = 0.0667457;
		b = -0.000672376;
	}
	
	else if(0.03<I_0<=0.04)
	{
		a = 0.106105;
		b = -0.00171939;
	}
	
	else if(0.04<I_0<=0.05)
	{
		a = 0.141998;
		b = -0.00308074;
	}
	
	else if(0.05<I_0<=0.06)
	{
		a = 0.164981;
		b = -0.00414971;
	}
	
	else if(0.05<I_0<=0.06)
	{
		a = 0.164981;
		b = -0.00414971;
	}
	
	else if(0.06<I_0<=0.07)
	{
		a = 0.178483;
		b = -0.00474469;
	}
	
	else if(0.07<I_0<=0.08)
	{
		a = 0.197539;
		b = -0.00567825;
	}
	
	else if(0.08<I_0<=0.09)
	{
		a = 0.208407;
		b = -0.00598107;
	}
	
	else if(0.09<I_0<=0.1)
	{
		a = 0.211437;
		b = -0.00543137;
	}
	
	else if(0.01<I_0<=0.5)
	{
		a = 0.369118;
		b = -0.021008;
	}
	
	else if(0.5<I_0<=1)
	{
		a = 1.78861;
		b = -0.48471;
	}
	
	else if(1<I_0<=10)
	{
		a = 0.423015;
		b = 4.1213;
	}
	
	else if(10<I_0<=1000)
	{
		a = 0.0138964;
		b = 745.762;
	}
	
	m_pParameter->LambdaValue = a * I_0 + b;

	return CCT_NOERR;
}
