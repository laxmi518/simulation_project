#pragma once
class COutputFrame
{
public:
	CParticle* m_OutputParticles;
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

	Scalar m_CurrentTime;
	std::vector<Scalar3> m_ModelPosition;

	int m_RotationTimes;

	std::vector<Scalar> m_RotateAngle;
private:
	Integer m_ParticleNum;
	bool m_Empty;
	std::string m_FileName;
	// Thread Management
	//! handle to the thread threadProc
public:
	HANDLE m_ThreadHandle;

public:
	COutputFrame()
		:m_OutputParticles(NULL)
		,m_ParticleNum(0)
		,m_Empty(true)
		,m_ThreadHandle(NULL)
		,m_CurrentTime(0.0)
		,m_RotationTimes(1)
	{
	}
	~COutputFrame()
	{
		m_FileName.clear();
		m_ModelPosition.clear();
		if(m_OutputParticles)
		{
			delete m_OutputParticles;
		}
		if(m_aOutputParticleID)
		{
			delete(m_aOutputParticleID);
		}
		if(m_aOutputParticlePosition)
		{
			delete(m_aOutputParticlePosition);
		}
		if(m_aOutputParticleVelocity)
		{
			delete(m_aOutputParticleVelocity);
		}
		if(m_aOutputParticlePressure)
		{
			delete(m_aOutputParticlePressure);
		}
		if(m_aOutputParticleDensity)
		{
			delete(m_aOutputParticleDensity);
		}
		if(m_aOutputParticleTemperature)
		{
			delete(m_aOutputParticleTemperature);
		}
		if(m_aOutputParticleKineticViscosity)
		{
			delete(m_aOutputParticleKineticViscosity);
		}
		if(m_aOutputParticleSolidPhaseRate)
		{
			delete(m_aOutputParticleSolidPhaseRate);
		}
		if(m_aOutputParticleType)
		{
			delete(m_aOutputParticleType);
		}
	}
	COutputFrame* GetOutputFrame()
	{
		if(m_Empty)
		{
			m_Empty = false;
			//CloseHandle(m_ThreadHandle);
			return this;
		}
		return NULL;
	}
	bool IsEmpty()
	{
		return m_Empty;
	}
	void SetEmpty(bool Empty)
	{
		m_Empty = Empty;
	}
	CCTStatusType CreateOutputFrame(Integer MaxParticleNum,Integer ModelNum)
	{
		CCTStatusType Status = CCT_NOERR;
		//m_OutputParticles = new CParticle[MaxParticleNum];
#if 0
			cudaMallocHost((void**)&m_aOutputParticleID,				MaxParticleNum * sizeof(Integer));
			cudaMallocHost((void**)&m_aOutputParticlePosition,			MaxParticleNum * sizeof(Scalar3));
			cudaMallocHost((void**)&m_aOutputParticleVelocity,			MaxParticleNum * sizeof(Scalar3));
			cudaMallocHost((void**)&m_aOutputParticlePressure,			MaxParticleNum * sizeof(Scalar));
			cudaMallocHost((void**)&m_aOutputParticleDensity,			MaxParticleNum * sizeof(Scalar));
			cudaMallocHost((void**)&m_aOutputParticleTemperature,		MaxParticleNum * sizeof(Scalar));
			cudaMallocHost((void**)&m_aOutputParticleKineticViscosity,	MaxParticleNum * sizeof(Scalar));
			cudaMallocHost((void**)&m_aOutputParticleSolidPhaseRate,	MaxParticleNum * sizeof(Scalar));
			cudaMallocHost((void**)&m_aOutputParticleType,				MaxParticleNum * sizeof(ParticleType));
#else
	        /*Rajan Individual Starts*/
		m_aOutputParticleID					= new Integer[MaxParticleNum];
		m_aOutputParticlePosition			= new Scalar3[MaxParticleNum];
		m_aOutputParticleVelocity			= new Scalar3[MaxParticleNum];
		m_aOutputParticlePressure			= new Scalar[MaxParticleNum];
		m_aOutputParticleDensity			= new Scalar[MaxParticleNum];
		m_aOutputParticleTemperature		= new Scalar[MaxParticleNum];
		m_aOutputParticleKineticViscosity	= new Scalar[MaxParticleNum];
		m_aOutputParticleSolidPhaseRate		= new Scalar[MaxParticleNum];
		m_aOutputParticleType				= new ParticleType[MaxParticleNum];

                /*Rajan Individual Ends */
#endif

		/*if(!m_aOutputParticleID || m_aOutputParticlePosition || m_aOutputParticleVelocity || )
		{
			return CCT_ETCERR;
		}*/
		if(!m_aOutputParticleID)
		{
			return CCT_ETCERR;
		}
		m_ModelPosition.resize(ModelNum,make_Scalar3(0,0,0));
		m_RotateAngle.resize(ModelNum,0.0);
		return Status;
	}
	void SetOutputParams(Integer ParticleNum,Scalar AnalysisTime,std::string& Path)
	{
		m_ParticleNum = ParticleNum;
		m_CurrentTime = AnalysisTime;
		m_FileName = Path;
	}
	//void SetModelPosition(Integer Step, Scalar dt, const std::vector<Scalar3>& ModelVelocity, const std::vector<Scalar3>& InitialModelPosition)
	void SetModelPosition(Integer TotalStep, Scalar dt , const std::vector<Scalar3>& ModelVelocity , const std::vector<Scalar3>& InitialModelPosition , 
							const std::vector<Integer>& GetModelResetInterval	,	 const std::vector<bool>& GetModelResetToOriginal)
	{
		Integer Step = 0;

		Scalar Lapse = 0;
		for(size_t i = 0;i < m_ModelPosition.size();++i)
		{
			if(GetModelResetToOriginal[i])
			{
				Step = TotalStep % GetModelResetInterval[i];
				Lapse = dt * (Step);
				m_ModelPosition[i].x = Lapse * ModelVelocity[i].x + InitialModelPosition[i].x;
				m_ModelPosition[i].y = Lapse * ModelVelocity[i].y + InitialModelPosition[i].y;
				m_ModelPosition[i].z = Lapse * ModelVelocity[i].z + InitialModelPosition[i].z;
			}
		}	
	}
	void SetRotationAngle(const int rotationTimes, const std::vector<Scalar>& RotationAngles,const std::vector<Scalar>& InitialTimeForZeroAngle,const std::vector<Scalar>& FinalTimeToReachAngle,const CParameter* Parameter)
	{
		//Scalar Lapse = dt * Step;
		for(Integer i = 0;i < RotationAngles.size();++i)
		{
			//To make Slope of the Rotation Starts
				Scalar angleOfRotation = 0.0;
		
				Scalar timePassed = Parameter->Dt * rotationTimes;

				int requiredStep = ( FinalTimeToReachAngle[i] - InitialTimeForZeroAngle[i]) / Parameter->Dt;
				
				if((timePassed < InitialTimeForZeroAngle[i]) || (RotationAngles[i] == 0))
				{
					angleOfRotation = 0.0;
				}
				else if((timePassed >= InitialTimeForZeroAngle[i]) && (timePassed < FinalTimeToReachAngle[i]))
				{					
					angleOfRotation = (RotationAngles[i] / (requiredStep + 1)) * (rotationTimes - (InitialTimeForZeroAngle[i] / Parameter->Dt));					
				}
				else
				{
					angleOfRotation = RotationAngles[i];
				}
			//To make Slope of the Rotation Ends
				m_RotateAngle[i] += angleOfRotation; // * rotationTimes;
		}		
	}
	CCTStatusType Output(bool IsAscii);

	void RestSetModelPosition(int modelNumberLess)
	{
		//int modelNumber = modelNumber + 1;
		////for(size_t i = 0;i < m_ModelPosition.size();++i)
		//if(modelNumber < m_ModelPosition.size())
		//{
		//	m_ModelPosition[modelNumber].x = 0;
		//	m_ModelPosition[modelNumber].y = 0;
		//	m_ModelPosition[modelNumber].z = 0;
		//}
	}
};

class COutputHandler
{
private:
	
public:
	COutputHandler(void);
	virtual ~COutputHandler(void);

private:
	/**
	*	@brief The thread procedure 
	*	
	*	@param LPVOID value is void pointer
	*/
	static DWORD WINAPI __stdcall threadProc(LPVOID value);
public:
	//Device
	static CCTStatusType Output(COutputFrame* OutputFrame, bool ISAscii);
	static bool m_IsAscii;
};
