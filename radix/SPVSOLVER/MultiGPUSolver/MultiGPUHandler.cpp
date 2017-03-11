#include "std.h"
#include "MultiGPUHandler.h"
//#include "radixsort.h"
#include "CGPUCellSet.h"
#include "Utility_Inline.h"
#include "IO.h"
#include "Configuration.h"
#include "OutputHandler.h"
#include "ETProfile.h"
#include "DtFrame.h"
#include "OptionParser.h"
#include "Logger.h"

#define MULTIBUCKET

CMultiGPUHandler::CMultiGPUHandler(void)
:m_GridNum(0)
,m_aGrid(NULL)
,m_aResults(NULL)
,m_MaxParticleNum(0)
,m_Output(NULL)
,m_pParameter(NULL)

,m_pParameterCoefficient(NULL)

,m_pConfiguration(NULL)
,m_ParticleNum(0)
,m_pOutputFrame(NULL)
,m_aOutputFrame(NULL)
,m_OutputFrameNum(0)
,m_haDragTriangles(NULL)
,m_DragTriangleNum(0)
,m_StartStep(0)
,m_pDt(NULL)
,m_NextItr(0)

,m_aOutputParticleID(0)
,m_aOutputParticlePosition(0)
,m_aOutputParticleVelocity(0)
,m_aOutputParticlePressure(0)
,m_aOutputParticleDensity(0)
,m_aOutputParticleTemperature(0)
,m_aOutputParticleKineticViscosity(0)
,m_aOutputParticleSolidPhaseRate(0)
,m_aOutputParticleType(NULL)

,m_MultiBucket(NULL)
,m_NumberOfBuckets(0)

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
	if(m_aResults)
	{
		delete m_aResults;
		m_aResults = NULL;
	}
	if(m_MultiBucket)
	{
		delete m_MultiBucket;
		m_MultiBucket = NULL;
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

	MaxVal =  max(m_pParameter->LaplacianInfluenceRadiusCoefficient,m_pParameter->GradientInfluenceRadiusCoefficient);
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
CCTStatusType CMultiGPUHandler::Initialize(Integer ParticleNum,
										   Integer * haParticleID, Scalar3* haParticlePosition, Scalar3* haParticleVelocity,
										Scalar* haParticlePressure, Scalar* haParticleDensity, Scalar* haParticleTemperature,
										Scalar* haParticleKineticViscosity, Scalar* haParticleSolidPhaseRate, ParticleType* haParticleType)
{	
	CalculateParameterCoefficient();
	CCTStatusType Status = CCT_NOERR;
	Status =  InitializeGrids();
	CCT_ERROR_CHECK(Status);
	InitializeParameters();
	InitializeParticleData(ParticleNum,
						haParticleID, haParticlePosition, haParticleVelocity,
						haParticlePressure,haParticleDensity, haParticleTemperature,
						haParticleKineticViscosity, haParticleSolidPhaseRate, haParticleType);
	
	InitializeTriangleData();
	Status = InitializeOutputBuffer();
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

CCTStatusType CMultiGPUHandler::InitializeGrids()
{
	Integer3 Section = m_pConfiguration->GetSection();
	m_GridNum = Section.x * Section.y * Section.z;
	m_aGrid = new CCellSet[m_GridNum];
	m_aResults = new Scalar[m_GridNum];
	if(!m_aGrid)
	{
		printf("ERR: STL file load");("ERR: Memory");
		return CCT_MEMERR;
	}
	InitializeComputeZone(Section, m_pParameter);
	/*Scalar MaxVal;
	MaxVal =  max(m_pParameter->LaplacianInfluenceRadiusCoefficient,m_pParameter->GradientInfluenceRadiusCoefficient);
	if(m_pParameter->bUpdateInnerPressure)
	{
		MaxVal = max(m_pParameter->SurfaceInfluenceRadiusCoefficient, MaxVal);
	}
	Scalar MaxInfluenceRadius = m_pParameter->InitialParticleDistance * MaxVal;
	Scalar3 Size = m_ComputeZone.m_MaxBound - m_ComputeZone.m_MinBound ;
	Scalar3 GridSize = Size / MaxInfluenceRadius;
	m_CalculationBox.m_ComputeZone = m_ComputeZone;
	m_CalculationBox.m_BufferedZone = m_ComputeZone;
	m_CalculationBox.m_CellSize = MaxInfluenceRadius;
	m_CalculationBox.m_GridSize =  make_Integer3((Integer)ceil(GridSize.x), (Integer)ceil(GridSize.y), (Integer)ceil(GridSize.z));
	Integer Remainder = m_CalculationBox.m_GridSize.x % m_GridNum;
	if(Remainder != 0)
	{
		m_CalculationBox.m_GridSize.x += Remainder;
	}
	InitializeCells();*/
	CCTStatusType Status;
	Status = InitializeHostMemory();
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
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
void CMultiGPUHandler::InitializeParticleData(Integer ParticleNum,
										Integer * haParticleID, Scalar3* haParticlePosition, Scalar3* haParticleVelocity,
										Scalar* haParticlePressure, Scalar* haParticleDensity, Scalar* haParticleTemperature,
										Scalar* haParticleKineticViscosity, Scalar* haParticleSolidPhaseRate, ParticleType* haParticleType)
{
	m_ParticleNum = ParticleNum;
	for(Integer j = 0; j < m_GridNum;++j) 
	{
		m_aGrid[j].SetParticleData(ParticleNum,
									haParticleID,		haParticlePosition,		haParticleVelocity,			haParticlePressure, 
									haParticleDensity,	haParticleTemperature,	haParticleKineticViscosity,	haParticleSolidPhaseRate,
									haParticleType);
	}
}
void CMultiGPUHandler::InitializeTriangleData()
{
	for(Integer j = 0; j < m_GridNum;++j) 
	{
		m_aGrid[j].SetTriangles(m_haTriangles, m_TriangleNum);
		m_aGrid[j].SetTrianglesParameters(m_haTriangleParameters,m_TriangleModelNumber);
		
		m_aGrid[j].SetDragParameters(m_haSTLDragPrameter, m_DragTriangleNum);
	}
}
void CMultiGPUHandler::InitializeDragTriangleData()
{
	for(Integer j = 0; j < m_GridNum;++j) 
	{
		m_aGrid[j].SetDragTriangles(m_haDragTriangles, m_DragTriangleNum/*,m_pInnerPressrueModel*/);
	}
}
void CMultiGPUHandler::InitializeStartStep()
{
	for(Integer j = 0; j < m_GridNum;++j) 
	{
		m_aGrid[j].SetStartStep(m_StartStep);
	}
}

CCTStatusType CMultiGPUHandler::AddToOutputBuffer(Integer BufferSize,
												  Integer* ParticleBufferID, Scalar3* ParticleBufferPosition,
												  Scalar3* ParticleBufferVelcity, Scalar* ParticleBufferPressure,
												  Scalar* ParticleBufferDensity, Scalar* ParticleBufferTemperature, Scalar* ParticleBufferKineticViscosity,
												  Scalar* ParticleBufferSolidPhaseRate, ParticleType* ParticleBufferType)
{	
	for(Integer i = 0;i < BufferSize; ++ i)
	{		
		m_aOutputParticlePosition[i]			=	ParticleBufferPosition[i];
		m_aOutputParticleVelocity[i]			=	ParticleBufferVelcity[i];
		m_aOutputParticlePressure[i]			=	ParticleBufferPressure[i];
		m_aOutputParticleDensity[i]				=	ParticleBufferDensity[i];	
		m_aOutputParticleKineticViscosity[i]	=	ParticleBufferKineticViscosity[i];
		m_aOutputParticleType[i]				=	ParticleBufferType[i];
		m_aOutputParticleTemperature[i]			=	ParticleBufferTemperature[i];
	}
	m_ParticleNum = BufferSize;
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::WaitForOutputBuffer()
{
	WaitForSingleObject(m_Output,INFINITE);
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::Output(Integer Step, std::string OutputType){
	CCTStatusType Status = CCT_NOERR;

	for(int i = 0 ;i < m_GridNum; ++i)
	{
		m_aGrid[i].AddToBuffer(true, m_ParticleNum,					m_aOutputParticleID, 	m_aOutputParticlePosition,		m_aOutputParticleVelocity,
								m_aOutputParticlePressure,			m_aOutputParticleDensity,		m_aOutputParticleTemperature,
								m_aOutputParticleKineticViscosity,	m_aOutputParticleSolidPhaseRate,m_aOutputParticleType);
	}
	std::stringstream SaveTo;
	SaveTo<<m_pConfiguration->GetOutputFileName()<<"_"<<OutputType<<"_"<<m_pConfiguration->GetOutputInterval()<<"_";
	if(Step < 10)
	{
		SaveTo<<"0"<<Step<<".dat";
	}		
	else
	{
		SaveTo<<Step<<".dat";
	}
	Integer outPutInterval	= m_pConfiguration->GetOutputInterval();
	Scalar AnalysisTime		= m_pParameter->Dt * m_pConfiguration->GetOutputInterval() * Step;
	m_pOutputFrame->SetOutputParams(m_ParticleNum, AnalysisTime,SaveTo.str());
	Status	=	COutputHandler::Output(m_pOutputFrame, m_pConfiguration->IsAscii());
	CCT_ERROR_CHECK(Status);
	std::cout<<"("<<Step<<") "<<OutputType<<" Out"<<std::endl;
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::Run(int argc, char** argv)
{
	// parse command line arguments
	SystemUtility::COptionParser Parser;
	SystemUtility::PlaneArgs PlaneArg;
	SystemUtility::KeyValueArgs KeyValueArg;
	Parser.Parse(argc, argv, PlaneArg, KeyValueArg);

	CCTStatusType StatusType = CCT_NOERR;
	if(PlaneArg.size() < 2)
	{
		return CCT_PARAMERR;
	}
	int GPUNum;
	cudaGetDeviceCount((int*)&GPUNum);
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
	//StatusType =  m_pConfiguration->ReadConfigurationFile(argv[1]);
	std::cout<<"INFO :- READ CONFIGURATION STARTS\n";
	StatusType =  m_pConfiguration->ReadConfigurationFileWithDrag(argv[1]);
	CCT_ERROR_CHECK(StatusType);
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
	StatusType = IO::Load(m_pConfiguration->GetParameterFile(), m_pParameter);
	CCT_ERROR_CHECK(StatusType);	
	std::cout<<"INFO :- LOAD PARAMETER FILE ENDS\n";

	// Movable Wall model STL triangles load
	std::cout<<"INFO :- LOAD STL MODELS FILE STARTS\n";
	StatusType = IO::Load(m_pConfiguration->GetModelFiles(), m_TrianglePerModel,m_TriangleNum, m_haTriangles);
	CCT_ERROR_CHECK(StatusType);
	printf("Loaded %d Triangles\n",m_TriangleNum);
	std::cout <<"Number of Models is :"<< m_TrianglePerModel.size()<<"\n";
	std::cout<<"INFO :- LOAD STL MODELS FILE ENDS\n";
	
	//Initialize the TriangleParameters
	m_TriangleModelNumber = m_TrianglePerModel.size();
	if(m_TriangleModelNumber > 0)
	{
		m_haTriangleParameters = new CTriangleParameters[m_TriangleModelNumber];
	}
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
		Scalar3 Velocity = m_pConfiguration->GetModelVelocity(i);
		Scalar Temperature = m_pConfiguration->GetModelTemperature(i);

		std::cout <<"Temperature is : " << Temperature <<std::endl;

		Scalar3 Position = m_pConfiguration->GetModelPosition(i);

		bool ResetToOriginal = m_pConfiguration->GetModelResetToOriginal(i);
		Integer ResetInterval = m_pConfiguration->GetModelResetInterval(i);

		//Paramter Setting for STL Rotation Starts

		bool RotatingOrNot				= m_pConfiguration->GetTrianlgeRotatingOrNot(i);
		//Scalar centerDegree			= m_pConfiguration->GetAngleDegree(i);
		Integer RotationInRPM			= m_pConfiguration->GetRotationRPM(i);
			
		//Calculating angle in Degree from Given RPM
		Scalar RotationInDegree  =	RotationInRPM * m_pParameter->Dt * 6;  // Angle of Rotation = (RPM /60) * dt * 360

		bool NegativePressureCheck = m_pConfiguration->GetNegativePressureCheck(i);

		m_pConfiguration->SetRotationInAngle(RotationInDegree);

		Scalar StartTimeOfRotation		= m_pConfiguration->GetRotationStartTime(i);
		Scalar TimeToReachRPM			= m_pConfiguration->GetTimeToReachRPM(i);
		Scalar3 centerOfRotation		= m_pConfiguration->GetCenterOfRotation(i);
		Scalar3 secondPointOfRotation	= m_pConfiguration->GetOtherPointOfRotation(i);
		Scalar  InflueneRegion	        = m_pConfiguration->GetInflueneRegion(i);
		Scalar3 DirectionVector		    = m_pConfiguration->GetDirectionVector(i);
		Scalar  PressureGas				= m_pConfiguration->GetPressureGas(i);

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
		m_haTriangleParameters[i].InfluenceRegion = InflueneRegion;
		m_haTriangleParameters[i].DirectionVector = DirectionVector;
		m_haTriangleParameters[i].PressureGas	  = PressureGas;

		m_haTriangleParameters[i].AngleOfRotationInDegree = RotationInDegree;
		m_haTriangleParameters[i].RotationInRPM = RotationInRPM;

		m_haTriangleParameters[i].NegativePressureCheck = NegativePressureCheck;
		//If any Triangle is Movable we need to update Trianlge
		//Check if any trianlge is linear movable or not
		if(ResetToOriginal)
		{
			m_IsAnyTriangleLinearMovable = true;
			m_IsReTriangleRegisterNecessary = true;
		}
		//If any Tringle is Rotating we need to rotate and upate the Triangle
		//Check if any trinagle is Roating or not
		if(RotatingOrNot)
		{
			m_IsAnyTriangleRotate = true;
			m_IsReTriangleRegisterNecessary = true;
		}
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
	StatusType = IO::Load(m_pConfiguration->GetParticleFile(), ParticleNum, InitialParticles);
	CCT_ERROR_CHECK(StatusType);
	printf("Loaded %d Input Particles\n",ParticleNum);
	std::cout<<"INFO :- LOAD PARTICLE FILE ENDS\n";
#ifdef MULTIBUCKET
	///Read Multiple Bucket Files Starts
	StatusType = LoadMultiBucket();
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
	///Read Multiple Bucket Files Ends

	m_MaxParticleNum = ParticleNum + totalBucketNumber;
	std::cout <<"Maximum number of Particles is : "<< m_MaxParticleNum <<std::endl;
#else
	m_MaxParticleNum = ParticleNum + m_BucketNum * m_pParameter->TotalGenerationStep;
#endif

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
	Integer			*initialParticleID;
	Scalar3			*initialParticlePosition;
	Scalar3			*initialParticleVelocity;
	Scalar			*initialParticlePressure;
	Scalar			*initialParticleDensity;
	Scalar			*initialParticleTemperature;
	Scalar			*initialParticleKineticViscosity;
	Scalar			*initialParticleSolidPhaseRate;
	ParticleType	*initialParticleType;

	initialParticleID		= new Integer[ParticleNum];
	initialParticlePosition = new Scalar3[ParticleNum];
	initialParticleVelocity = new Scalar3[ParticleNum];
	initialParticlePressure = new Scalar[ParticleNum];
	initialParticleDensity	= new Scalar[ParticleNum];
	initialParticleTemperature		= new Scalar[ParticleNum];
	initialParticleKineticViscosity = new Scalar[ParticleNum];
	initialParticleSolidPhaseRate	= new Scalar[ParticleNum];
	initialParticleType				= new ParticleType[ParticleNum];
	//ENDS for individual Output particles----------------------E&T Nepal August 2011--------------------------
	
	for(Integer i = 0 ; i < ParticleNum ; i++)
	{
		initialParticleID[i]				= InitialParticles[i].ID;
		initialParticlePosition[i]			= InitialParticles[i].Position;
		initialParticleVelocity[i]			= InitialParticles[i].Velocity;
		initialParticlePressure[i]			= InitialParticles[i].Pressure;
		initialParticleDensity[i]			= InitialParticles[i].Density;
		initialParticleTemperature[i]		= InitialParticles[i].Temperature;
		initialParticleKineticViscosity[i]	= InitialParticles[i].KineticViscosity;
		initialParticleSolidPhaseRate[i]	= InitialParticles[i].SolidPhaseRate;
		initialParticleType[i]				= InitialParticles[i].Type;
	}
	//ENDS for individual Output particles----------------------E&T Nepal August 2011--------------------------

	CCTStatusType Status = Initialize(ParticleNum,
									initialParticleID,initialParticlePosition,initialParticleVelocity, initialParticlePressure,initialParticleDensity, initialParticleTemperature,
									initialParticleKineticViscosity,initialParticleSolidPhaseRate,initialParticleType);
	CCT_ERROR_CHECK(Status);

	delete InitialParticles;
	delete initialParticleID;
	delete initialParticlePosition;
	delete initialParticleVelocity;
	delete initialParticlePressure;
	delete initialParticleDensity;
	delete initialParticleTemperature;
	delete initialParticleKineticViscosity;
	delete initialParticleSolidPhaseRate;
	delete initialParticleType;

		if(PlaneArg.size() >= 3)
	{
		std::string FileName = argv[2];
		if(!FileName.empty())
		{
			std::string InnerPressureParameterFile;
			Status = IO::Load(FileName,m_haDragTriangles,m_DragTriangleNum,InnerPressureParameterFile);
			CCT_ERROR_CHECK(Status);
			std::cout<<"Loaded :"<< m_DragTriangleNum << "Drag Triangles"<<std::endl;

			//m_pInnerPressrueModel = new CInnerPressureModel(m_pParameter);
			//Status = m_pInnerPressrueModel->Load(InnerPressureParameterFile);
			CCT_ERROR_CHECK(Status);
//			m_pInnerPressrueModel->SetOutputPath(m_pConfiguration->GetOutputPath());

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
	Status = SetOutputPath();
	CCT_ERROR_CHECK(Status);

	m_StartStep = m_pConfiguration->GetStartStep();
	InitializeStartStep();	
	OutputComputeZone();
	//Save The used Parameter File After LambdaValue Calcuation
	std::stringstream UsedParameterFileName;
	UsedParameterFileName << m_pConfiguration->GetOutputPath() << "\\UsedParameters.txt";
//	IO::saveParameterWithName(UsedParameterFileName.str(), m_pParameter);
	IO::saveParameterWithName(UsedParameterFileName.str(),				m_pParameter,								m_TriangleModelNumber,
								m_pConfiguration->GetCenterOfRotation(),m_pConfiguration->GetOtherPointOfRotation(),m_pConfiguration->GetRotationRPM(),
								m_pConfiguration->GetWallFrictionCoefficient(),m_pConfiguration->GetWallThermalConductivity(),m_pConfiguration->GetWallThermalRestivity(),								
								m_haSTLDragPrameter,m_DragTriangleNum);

	//Save The Used Parameter File Ends

	///////Test For the Maximum number of Particles/////////////////////
	std::cout <<"The Maximum Number of Particles : " << m_MaxParticleNum << std::endl;
		
	///////Test For the Maximum number of Particles/////////////////////
	Status = this->Calculate();
	CCT_ERROR_CHECK(Status);
	return CCT_NOERR;
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
		//Start = m_TrianglePerModel[Model-1];
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

//CCTStatusType CMultiGPUHandler::CalculateLambdaValueLocal()
//{	
//	const Scalar I_0 = m_pParameter->InitialParticleDistance * 1E-3; //conversion into metre scale
//	Scalar a = 0.0, b = 0.0;
//		
//	if( 0.0 <I_0 && I_0 <=0.0001)
//	{
//		a = static_cast<Scalar>(5.02608E-06);
//		b = static_cast<Scalar>(3.82946E-10);
//	}
//	
//	else if(0.0001<I_0 && I_0 <=0.0005)
//	{
//		a = static_cast<Scalar>(0.000124918);
//		b = static_cast<Scalar>(1.81984E-08);
//	}
//	
//	else if(0.0005<I_0 && I_0<=0.001)
//	{
//		a = static_cast<Scalar>(0.0017886);
//		b = static_cast<Scalar>(-4.84707E-07);
//	}
//	
//	else if(0.001<I_0&& I_0<=0.002)
//	{
//		a = static_cast<Scalar>(0.00373277);
//		b = static_cast<Scalar>(-2.15486E-06);
//	}
//	
//	else if(0.002<I_0&& I_0<=0.003)
//	{
//		a =static_cast<Scalar>( 0.0066746);
//		b = static_cast<Scalar>(-6.72384E-06);
//	}
//	
//	else if(0.003<I_0&& I_0<=0.004)
//	{
//		a = static_cast<Scalar>(0.0107452);
//		b = static_cast<Scalar>(-1.76386E-05);
//	}
//	
//	else if(0.004<I_0&& I_0<=0.005)
//	{
//		a = static_cast<Scalar>(0.0141993);
//		b = static_cast<Scalar>(-3.08052E-05);
//	}
//	
//	else if(0.005<I_0&& I_0<=0.006)
//	{
//		a = static_cast<Scalar>(0.0164975);
//		b = static_cast<Scalar>(-4.14937E-05);
//	}
//	
//	else if(0.006<I_0&& I_0<=0.007)
//	{
//		a = static_cast<Scalar>(0.0181782);
//		b = static_cast<Scalar>(-4.95251E-05);
//	}
//	
//	else if(0.007<I_0&& I_0<=0.008)
//	{
//		a = static_cast<Scalar>(0.0205885);
//		b = static_cast<Scalar>(-6.28759E-05);
//	}
//	
//	else if(0.008<I_0&& I_0<=0.009)
//	{
//		a = static_cast<Scalar>(0.0222761);
//		b = static_cast<Scalar>(-7.17244E-05);
//	}
//	
//	else if(0.009<I_0&& I_0<=0.01)
//	{
//		a = static_cast<Scalar>(0.0231972);
//		b = static_cast<Scalar>(-7.34117E-05);
//	}
//	
//	else if(0.01<I_0&& I_0<=0.02)
//	{
//		a = static_cast<Scalar>(0.0369119);
//		b = static_cast<Scalar>(-0.000210081);
//	}
//	
//	else if(0.02<I_0&& I_0<=0.03)
//	{
//		a = static_cast<Scalar>(0.0667457);
//		b = static_cast<Scalar>(-0.000672376);
//	}
//	
//	else if(0.03<I_0&& I_0<=0.04)
//	{
//		a = static_cast<Scalar>(0.106105);
//		b =static_cast<Scalar>( -0.00171939);
//	}
//	
//	else if(0.04<I_0&& I_0<=0.05)
//	{
//		a = static_cast<Scalar>(0.141998);
//		b = static_cast<Scalar>(-0.00308074);
//	}
//	
//	else if(0.05<I_0&& I_0<=0.06)
//	{
//		a = static_cast<Scalar>(0.164981);
//		b = static_cast<Scalar>(-0.00414971);
//	}
//	
//	else if(0.05<I_0&& I_0<=0.06)
//	{
//		a = static_cast<Scalar>(0.164981);
//		b = static_cast<Scalar>(-0.00414971);
//	}
//	
//	else if(0.06<I_0&& I_0<=0.07)
//	{
//		a = static_cast<Scalar>(0.178483);
//		b = static_cast<Scalar>(-0.00474469);
//	}
//	
//	else if(0.07<I_0&& I_0<=0.08)
//	{
//		a = static_cast<Scalar>(0.197539);
//		b = static_cast<Scalar>(-0.00567825);
//	}
//	
//	else if(0.08<I_0&& I_0<=0.09)
//	{
//		a = static_cast<Scalar>(0.208407);
//		b = static_cast<Scalar>(-0.00598107);
//	}
//	
//	else if(0.09<I_0&& I_0<=0.1)
//	{
//		a = static_cast<Scalar>(0.211437);
//		b = static_cast<Scalar>(-0.00543137);
//	}
//	
//	else if(0.01<I_0&& I_0<=0.5)
//	{
//		a = static_cast<Scalar>(0.369118);
//		b = static_cast<Scalar>(-0.021008);
//	}
//	
//	else if(0.5<I_0&& I_0<=1)
//	{
//		a = static_cast<Scalar>(1.78861);
//		b = static_cast<Scalar>(-0.48471);
//	}
//	
//	else if(1<I_0&& I_0<=10)
//	{
//		a = static_cast<Scalar>(0.423015);
//		b = static_cast<Scalar>(4.1213);
//	}
//	
//	else if(10<I_0&& I_0<=1000)
//	{
//		a = static_cast<Scalar>(0.0138964);
//		b = static_cast<Scalar>(745.762);
//	}
//	
//	m_pParameter->LambdaValue = (a * I_0 + b) * 1E6; //conversion into mm scale
//
//	return CCT_NOERR;
//}

void CMultiGPUHandler::WaitForOutputFrame()
{
	Integer i = 0;
	while(true)
	{
		m_pOutputFrame = m_aOutputFrame[i].GetOutputFrame();
		if(m_pOutputFrame)
		{
			//STARTS for ResetOutputParticles----------------------E&T Nepal August 2011--------------------------
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
			//ENDS for Reset OutputParticles----------------------E&T Nepal August 2011--------------------------
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
	Integer ModelNum = static_cast<Integer>(m_pConfiguration->GetModelVelocity().size());
	for(Integer i = 0 ;i < m_OutputFrameNum; ++i)
	{
		CCTStatusType Status =  m_aOutputFrame[i].CreateOutputFrame(m_MaxParticleNum,ModelNum);
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
CCTStatusType CMultiGPUHandler::GetDataOnMultiGPUHandler(bool Output)
{
	CCTStatusType Status = CCT_NOERR;	
	m_aGrid->AddToBuffer(Output, m_ParticleNum,		m_aOutputParticleID,	m_aOutputParticlePosition,
								m_aOutputParticleVelocity,m_aOutputParticlePressure,m_aOutputParticleDensity,
								m_aOutputParticleTemperature,m_aOutputParticleKineticViscosity,
								m_aOutputParticleSolidPhaseRate,m_aOutputParticleType);
	
	return CCT_NOERR;
}

CCTStatusType CMultiGPUHandler::AddBucketParticles(const CParticle * bucket, const Integer bucketNumber)
{
	InsertBucket(bucket, bucketNumber);
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
		//STARTS for Add Bucket on OutputParticles----------------------E&T Nepal August 2011--------------------------
		if(m_aOutputParticleID[i] < 0 )
		{
			m_aOutputParticleID[i]					= i;
			m_aOutputParticlePosition[i]			= bucket[j].Position;
			m_aOutputParticleVelocity[i]			= bucket[j].Velocity;
			m_aOutputParticlePressure[i]			= bucket[j].Pressure;
			m_aOutputParticleDensity[i]				= bucket[j].Density;
			m_aOutputParticleTemperature[i]			= bucket[j].Temperature;
			m_aOutputParticleKineticViscosity[i]	= bucket[j].KineticViscosity;
			m_aOutputParticleSolidPhaseRate[i]		= bucket[j].SolidPhaseRate;
			m_aOutputParticleType[i]				= bucket[j].Type;
			j++;
			if(j >= bucketNumber )
			{
				return;
			}			
		}
		//Ends for Add Bucket on OutputParticles----------------------E&T Nepal August 2011--------------------------
	}
	if(m_ParticleNum + bucketNumber - j <= m_MaxParticleNum )
	{
		
		//STARTS for Add Bucket on OutputParticles----------------------E&T Nepal August 2011--------------------------
		for(;j<bucketNumber;++j,++i)
		{
			m_aOutputParticleID[i]				= i;
			m_aOutputParticlePosition[i]		= bucket[j].Position;
			m_aOutputParticleVelocity[i]		= bucket[j].Velocity;
			m_aOutputParticlePressure[i]		= bucket[j].Pressure;
			m_aOutputParticleDensity[i]			= bucket[j].Density;
			m_aOutputParticleTemperature[i]		= bucket[j].Temperature;
			m_aOutputParticleKineticViscosity[i]= bucket[j].KineticViscosity;
			m_aOutputParticleSolidPhaseRate[i]  = bucket[j].SolidPhaseRate;
			m_aOutputParticleType[i]			= bucket[j].Type;
		}
		//Ends for Add Bucket on OutputParticles----------------------E&T Nepal August 2011--------------------------
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
				m_pParameter->ParticleGenerationStep = /*static_cast<Integer>*/(DtTime->InflowStep);
				m_pParameter->OutputStep = static_cast<Integer>(DtTime->OutputStep);
				++m_NextItr;
			}
		}
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler:: CalculateParameterCoefficient()
{
	m_pParameterCoefficient = new CParameterCoefficients;
	//const Scalar Re = CONSTANT_PARAMETER.LaplacianInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	m_pParameterCoefficient->ReLap = m_pParameter->LaplacianInfluenceRadiusCoefficient * m_pParameter->InitialParticleDistance;

	m_pParameterCoefficient->ReGrad = m_pParameter->GradientInfluenceRadiusCoefficient * m_pParameter->InitialParticleDistance;
	
	//const Scalar Rc = CONSTANT_PARAMETER.r_coeff * CONSTANT_PARAMETER.InitialParticleDistance;
	m_pParameterCoefficient->Rc = m_pParameter->r_coeff * m_pParameter->InitialParticleDistance;

	//const Scalar Re = CONSTANT_PARAMETER.SurfaceInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	m_pParameterCoefficient->ReSurface = m_pParameter->SurfaceInfluenceRadiusCoefficient * m_pParameter->InitialParticleDistance;

	//const Scalar MinDist = Re * 1e-3;
	m_pParameterCoefficient->MinimumDistance = static_cast<Scalar>(m_pParameterCoefficient->ReLap * 1e-3);

	//For CalcExplicitly_Kernel
	//Scalar Coefficient = ( 2.0 * CONSTANT_PARAMETER.Dimension ) / ( CONSTANT_PARAMETER.LambdaValue * CONSTANT_PARAMETER.ConstantParticleDensity);
	m_pParameterCoefficient->CoefficientCalcExplicit =  static_cast<Scalar>((2.0 * m_pParameter->Dimension ) / ( m_pParameter->LambdaValue * m_pParameter->ConstantParticleDensity));

	//For CalcTemperatureFactor_Kernel
	//const Scalar EnthalpyCoefficient = CONSTANT_PARAMETER.Dt * 2 * CONSTANT_PARAMETER.Dimension / ( CONSTANT_PARAMETER.ConstantParticleDensity * Density * CONSTANT_PARAMETER.FruidSpecificHeat);
	m_pParameterCoefficient->EnthalpyCoefficientCalcTempFactor = m_pParameter->Dt * 2 * m_pParameter->Dimension /(m_pParameter->ConstantParticleDensity * m_pParameter->FruidSpecificHeat);
	//ThermalConductivity = (2 * CONSTANT_PARAMETER.FruidThermalConductivity * CONSTANT_PARAMETER.WallThermalConductivity ) / ( CONSTANT_PARAMETER.FruidThermalConductivity + CONSTANT_PARAMETER.WallThermalConductivity);
	m_pParameterCoefficient->ThermalConductivityCalcTempFactor = (2 * m_pParameter->FruidThermalConductivity * m_pParameter->WallThermalConductivity) / (m_pParameter->FruidThermalConductivity + m_pParameter->WallThermalConductivity);

	//For Iterate_Kernel
	//Scalar CoefficientC = 2 * CONSTANT_PARAMETER.Dimension / ( CONSTANT_PARAMETER.LambdaValue * CONSTANT_PARAMETER.ConstantParticleDensity);
	m_pParameterCoefficient->CoefficientCIterate = 2 * m_pParameter->Dimension / (m_pParameter->LambdaValue * m_pParameter->ConstantParticleDensity);

	//For IntiR_KernelCG
	//Scalar CoefficientC = 2 * CONSTANT_PARAMETER.Dimension / ( CONSTANT_PARAMETER.LambdaValue * CONSTANT_PARAMETER.ConstantParticleDensity);
	m_pParameterCoefficient->CoefficientCIntiR = 2 * m_pParameter->Dimension / (m_pParameter->LambdaValue * m_pParameter->ConstantParticleDensity);

	//Scalar PressureCoefficient = (m_pParameter->LambdaValue * m_pParameter->ConstantParticleDensity * m_pParameter->Density) 
	//	/ ( 2 * m_pParameter->Dimension * m_pParameter->Dt * m_pParameter->Dt);
	m_pParameterCoefficient->PressureCoefficient =  (m_pParameter->LambdaValue * m_pParameter->ConstantParticleDensity * m_pParameter->Density) 
		/ ( 2 * m_pParameter->Dimension * m_pParameter->Dt * m_pParameter->Dt);

	return CCT_NOERR;
}

void CMultiGPUHandler::InitializeCells()
{
	if((Integer)(m_DeviceSelector.Size()) < m_GridNum)
	{
		printf("ERR: Invalid Device Number(Number of GPUs is lesser than number of grids)\n");
	}

	for(int i = 0 ;i < m_GridNum; ++i)
	{
		m_aGrid[i].InitializeCells(m_DeviceSelector[i], i ,m_GridNum,m_CalculationBox);
		printf("Total Cells : %d", m_aGrid[i].getCellNum());
	}
}
//CCTStatusType CMultiGPUHandler::InitializeGPU()
//{	
//	CCTStatusType Status = CCT_NOERR;
//	for(int i =0; i < m_GridNum; ++i)
//	{
//		Status= m_aGrid[i].InitializeGPU();
//		CCT_ERROR_CHECK(Status);
//	}
//	return Status;
//}
CCTStatusType CMultiGPUHandler::TransferRequiredCheck(CCTStatusType Error, Integer Step)
{
	if(Error!= CCT_NOERR)
		{
			//Added to check the values of Constants
			CCTStatusType TransferStatus = this->ConstantsDeviceToHost(Step);
			CCT_ERROR_CHECK(TransferStatus);
		}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::Calculate()
{
	CCTStatusType Status = CCT_NOERR;
	CCTStatusType TransferStatus =CCT_NOERR;	
	Status = CudaSafeCall(cudaGetLastError());
	CCT_ERROR_CHECK(Status);
	
	Integer deviceID = m_DeviceSelector[0];

	std::cout<<"Max Number from CMultiGPUHandler is :- "<<m_MaxParticleNum<<std::endl;
	Status = m_aGrid->InitializeGPU(deviceID);
//	Status = m_aGrid->InitializeGPU();
	CCT_ERROR_CHECK(Status);
	std::cout<<"Started " << m_GridNum << " GPU"<<std::endl;
	m_Output = CreateEvent(NULL, TRUE, FALSE, NULL);
	if(m_Output == NULL)
	{
		return CCT_ETCERR;
	}	
	std::string TimerPath = m_pConfiguration->GetOutputPath();
	TimerPath.append("Times.txt");
	CETProfile::setFilePath(TimerPath.c_str());
	CETProfile Ep;
	Integer FileNo = m_StartStep /m_pConfiguration->GetOutputInterval();
	Integer AddedTImes = 0;
	Integer WallStep = 0;
	Integer i = m_StartStep;
	Integer STLRotationSteps = 0;

	//Display CUDA VERSION STARTS
	int driverVersion = 0, runtimeVersion = 0;
	cudaDriverGetVersion(&driverVersion);
    cudaRuntimeGetVersion(&runtimeVersion);
    //Display CUDA VERSION ENDS

	std::cout<<"""""""""""""""""""""""""MULTI BUCKET SOLVER WITH ROTATION CUDA VERSION "<< driverVersion/1000 <<"."<< (driverVersion%100)/10 <<"""""""""""""""""\n\n\n";
	
	for(; i < m_pParameter->AnalysisStep; ++i)
	{
		Status = CheckAndUpdateDt(i);
		CCT_ERROR_CHECK(Status);

		Status = this->RegisterParticleTopology(i);
		CCT_ERROR_CHECK(Status);
		std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : REGISTER PARTICLE\n";
		
		Status = this->RegisterTriangleTopology();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : REGISTER TRIANGLE\n";
		
		Status = this->CalculateSTLDistance();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : CALC STL DIST\n";
		
		if(m_DragTriangleNum > 0)		
		{
			Status = this->CalcDragEffect();
			CCT_ERROR_CHECK(Status);
			std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : CALC Drag\n";

		}

		if( true == m_pParameter->bUpdateTemperature && 0 == i % m_pParameter->ThermalStep)
		{
			Status = this->CalcTemperatureFactor();
			CCT_ERROR_CHECK(Status);
			std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : CALC TEMPERATURE\n";
		}

		//For Turbulance Calculation
		if(m_pParameter->bTurbulance)
		{
			Status = this->CalcTurbulenceViscosity();
			CCT_ERROR_CHECK(Status);
			std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : Calc Turbulence Viscosity\n";//Spelling Corrected
		}

		Status = this->CalcExplicitly();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : CALC EXPLICIT\n";
		
		
		Status = this->Synchronize();
		CCT_ERROR_CHECK(Status);

		if((i % m_pConfiguration->GetOutputInterval()) == 0 && m_pParameter->bExplicitOut)
		{
			Status = Synchronize();
			CCT_ERROR_CHECK(Status);

			Status = ParticleDeviceToHost();
			CCT_ERROR_CHECK(Status);
			WaitForOutputFrame();

			m_pOutputFrame->SetModelPosition(WallStep, m_pParameter->Dt, m_pConfiguration->GetModelVelocity(),m_pConfiguration->GetModelPosition(),
											m_pConfiguration->GetModelResetInterval(),	m_pConfiguration->GetModelResetToOriginal());
			Status = Synchronize();
			CCT_ERROR_CHECK(Status);
			Status = this->Output(FileNo,"Explicit");
			CCT_ERROR_CHECK(Status);
		}
		if(m_pParameter->ConstantCheckStep==i)
		{
			CCTStatusType TransferStatus = this->ConstantsDeviceToHost(i);
			CCT_ERROR_CHECK(TransferStatus);
		}
		Status = this->RegisterParticleTopology(i);
		CCT_ERROR_CHECK(Status);
		std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : REGISTER PARTICLE\n";
				
		Status = this->CalculateSTLDistance();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : CALC STL DISTANCE\n";
		
		Status = this->CalcExplicitPressure();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : CALC PRESSURE\n";
		
		Status = this->CalcExplicitPressureGradient();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : CALC PRESSURE GRADIENT\n";
		
		//For Saving The Rotation Angle in output File
		//Setup STL Rotaion
		if( i!=0 )
		{
			m_pOutputFrame->SetRotationAngle(i,m_pConfiguration->GetRotationAngle(), m_pConfiguration->GetRotationStartTime(),m_pConfiguration->GetTimeToReachRPM() ,m_pParameter);
		}
		//Output Steps Starts
		if((i % m_pConfiguration->GetOutputInterval()) ==0)
		{
			Status = Synchronize();
			CCT_ERROR_CHECK(Status);

			Status = ParticleDeviceToHost();
			CCT_ERROR_CHECK(Status);
			std::cout <<"***---***::: OUTPUT STEP :\t" <<  m_pConfiguration->GetOutputInterval() << "\n";
			WaitForOutputFrame();
			
			m_pOutputFrame->SetModelPosition(WallStep, m_pParameter->Dt,m_pConfiguration->GetModelVelocity(),m_pConfiguration->GetModelPosition(), m_pConfiguration->GetModelResetInterval(), m_pConfiguration->GetModelResetToOriginal());
			Status = Synchronize();
			CCT_ERROR_CHECK(Status);
					
			Status = this->Output(FileNo,"Implicit");
			CCT_ERROR_CHECK(Status);
			FileNo++;
		}
		else
		{
			Status = this->Synchronize();
			CCT_ERROR_CHECK(Status);
		}

		//Multiple Bucket Addition
		for(int p = 0 ; p < m_NumberOfBuckets ; ++p)
		{
			if((m_MultiBucket[p].AddedTimes < m_MultiBucket[p].TotalAdditionStep) && (0 == ( (i+1) % m_MultiBucket[p].AdditionInterval)))
			{
				Status = Synchronize();
				CCT_ERROR_CHECK(Status);
				
				Integer OldNum = m_ParticleNum;	

				Status = this->ParticleDeviceToHost();
				CCT_ERROR_CHECK(Status);

				WaitForOutputFrame();
	
				Status = this->GetDataOnMultiGPUHandler(false);
				CCT_ERROR_CHECK(Status);

				Status = this->AddBucketParticles(m_MultiBucket[p].Bucket, m_MultiBucket[p].BucketNumber);
				ET_ERROR_BREAK(Status);

				Status = this->RestWallPositions(i);
				ET_ERROR_BREAK(Status);

				m_pOutputFrame->RestSetModelPosition(p);

				Status = this->ParticleHostToDevice();	
				ET_ERROR_BREAK(Status);

				Status = this->ParticleNumberToContantMemory();
				CCT_ERROR_CHECK(Status);

				++m_MultiBucket[p].AddedTimes;

				std::cout<<"***---***::: New Particle Num :\t" << m_ParticleNum << " (FROM "<< OldNum << ")"<<std::endl;
			}
		}		
		//Bucket Addition Ends

		Status = this->MoveTriangles();
		CCT_ERROR_CHECK(Status);
		std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : MOVE TRIANGLE\n";
		
		//Rotation Triangles Starts
		Status = this->RotateTrianglePosition(i);		
		CCT_ERROR_CHECK(Status);
		std::cout <<"**RDDev: "<< m_DeviceSelector[0] <<" :Step: " << i << " :M: "<< m_MaxParticleNum <<" :P: " << m_ParticleNum << " :T: " << m_TriangleNum <<" : ROTATE TRIANGLE\n";
		
		if(i < m_pParameter->StopWallStep)
		{
			++WallStep;
		}
		STLRotationSteps++;
		Ep.DisplayTime(i);
	}
	WaitForOutputThreads();
	Ep.SetAverageTime(i);
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::SetOutputPath()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ;i < m_GridNum; ++i)
	{
		Status =m_aGrid[i].SetOutPath(m_pConfiguration->GetOutputPath());
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::RegisterParticleTopology(Integer step)
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ;i < m_GridNum; ++i)
	{
		Status =m_aGrid[i].RegisterParticleTopology(step);
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}

CCTStatusType CMultiGPUHandler::RegisterTriangleTopology()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].RegisterTriangleTopology();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::CalcExplicitly()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].CalcExplicitly();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::MoveTriangles()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].MoveTriangles();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::RotateTrianglePosition(const Integer analysisStep)
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].RotateTrianglePosition(analysisStep);
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::RestWallPositions(Integer AnalysisStep)
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].ResetWallPosition(AnalysisStep);
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::ParticleDeviceToHost()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i =0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].ParticleDeviceToHost();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::ParticleHostToDevice()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i =0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].ParticleHostToDevice();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::ParticleNumberToContantMemory()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i =0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].ParticleNumberToContantMemory();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::CalculateSTLDistance()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i =0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].CalculateSTLDistance();
		CCT_ERROR_CHECK(Status);
	}
	return Status;
}
CCTStatusType CMultiGPUHandler::CalcDragEffect()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i =0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].CalculateDragEffect();
		CCT_ERROR_CHECK(Status);
	}
	return Status;
}
void CMultiGPUHandler::InitializeParameters()
{
	for(int i=0; i<m_GridNum;++i)
	{
		m_aGrid[i].InitializeParameter(m_pParameter,m_pParameterCoefficient);
	}
}
CCTStatusType CMultiGPUHandler::CalcExplicitPressure()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].CalcExplicitPressure();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::CalcExplicitPressureGradient()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].CalcExplicitPressureGradient();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::CalcTemperatureFactor()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].CalcTemperatureFactor();
		CCT_ERROR_CHECK(Status);
	}
	return CCT_NOERR;
}
CCTStatusType CMultiGPUHandler::CalcTurbulenceViscosity()
{
	CCTStatusType Status = CCT_NOERR;
	for(int i = 0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].CalcTurbulenceViscosity();
		CCT_ERROR_CHECK(Status);
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
			m_MultiBucket[j].AdditionInterval = m_pConfiguration->GetMultiAdditionInterval()[j];
			m_MultiBucket[j].TotalAdditionStep = m_pConfiguration->GetMultiTotalAdditionStep()[j];
			m_MultiBucket[j].AddedTimes = 0;
		}
	}
	return Status;
}
CCTStatusType CMultiGPUHandler::ConstantsDeviceToHost(Integer Step)
{
	CCTStatusType Status = CCT_NOERR;
	for(int i =0 ;i < m_GridNum; ++i)
	{
		Status = m_aGrid[i].ConstantsDeviceToHost(Step);
		CCT_ERROR_CHECK(Status);
	}
	return Status;
}
/*End*/
