#include "stdafx.h"
#include "IO.h"
#include "STLSerializer.h"
#include "ETLuaState.h"


CCTStatusType IO::Load(const std::string& FileName, Integer& ParticleNum, CParticle*& phParticle)
{
	std::ifstream ifs(FileName.c_str());
	if(false == ifs.is_open())
	{
		printf("ERR: Particle file load");
		return CCT_FILEERR;
	}
	
	// 粒子タイプを見て壁粒子は読み込まないようにするなどの処理が必要
	//(Wall particle particle type to see the process as necessary to read them)

	std::vector<CParticle> vParticle;
	Integer Type = 0;
	while(!ifs.eof())
	{
		CParticle Particle;
		ifs >> Particle.Position.x >> Particle.Position.y >> Particle.Position.z;
		ifs >> Particle.Velocity.x >> Particle.Velocity.y >> Particle.Velocity.z;
		ifs >> Particle.Pressure >> Particle.Density >> Particle.Temperature 
			>> Particle.KineticViscosity >> Particle.SolidPhaseRate >> Type;
		if ( 0 == Type)
		{
			Particle.Type = TYPE_HEAVY_WEIGHT;
		}
		if ( 1 == Type)
		{
			Particle.Type = TYPE_HEAVY_WEIGHT_FREE_SURFACE;
		}
		if ( 2 == Type)
		{
			Particle.Type = TYPE_LIGHT_WEIGHT;
		}
		if(!ifs.fail())
		{
			vParticle.push_back(Particle);
		}
	}

	ParticleNum = vParticle.size();
	phParticle = new CParticle[ParticleNum];
	for(Integer i = 0; i < ParticleNum; ++i)
	{
		phParticle[i].Position = vParticle[i].Position;
		phParticle[i].Velocity = vParticle[i].Velocity;
		phParticle[i].Pressure = vParticle[i].Pressure;
		//phParticle[i].MinPressure = vParticle[i].Pressure;
		phParticle[i].Density = vParticle[i].Density;
		phParticle[i].Temperature = vParticle[i].Temperature;
		phParticle[i].KineticViscosity = vParticle[i].KineticViscosity;
		phParticle[i].SolidPhaseRate = vParticle[i].SolidPhaseRate;
		phParticle[i].Type = vParticle[i].Type;
		//phParticle[i].ExplicitPosition = vParticle[i].Position;
		//phParticle[i].ExplicitVelocity = vParticle[i].Velocity;
		phParticle[i].ID = i;
	}
	return CCT_NOERR;
}

CCTStatusType IO::Load(const std::string& FileName, CParameter*& pParameter)
{
	std::ifstream ifs(FileName.c_str());
	if(false == ifs.is_open())
	{
		printf("ERR: STL file load");("ERR: STL file load");("ERR: Parameter file couldn't open");
		return CCT_FILEERR;
	}

	pParameter = new CParameter;

	ifs >> pParameter->Dimension;	//1
	ifs >> pParameter->GravityConstant.x;  //2
	ifs >> pParameter->GravityConstant.y; //3
	ifs >> pParameter->GravityConstant.z; //4

	ifs >> pParameter->Density; //5
	ifs >> pParameter->ConstantParticleDensity; //6
	ifs >> pParameter->SurfaceTensionCoefficient; //7
	ifs >> pParameter->InitialParticleDistance; //8
	ifs >> pParameter->FreeSurfaceThreshold; //9

	ifs >> pParameter->GradientInfluenceRadiusCoefficient; //10
	ifs >> pParameter->LaplacianInfluenceRadiusCoefficient;
	ifs >> pParameter->TemperatureInfluenceRadiusCoefficient;

	ifs >> pParameter->FruidInitializeialTemperature;
	ifs >> pParameter->WallInitializeialTemperature;
	ifs >> pParameter->FruidThermalResistance;
	ifs >> pParameter->WallThermalResistance; //16
	ifs >> pParameter->FruidThermalConductivity;
	ifs >> pParameter->WallThermalConductivity; //18
	ifs >> pParameter->FruidSpecificHeat;
	ifs >> pParameter->FruidDensity;

	ifs >> pParameter->Dt;
	ifs >> pParameter->AnalysisStep;
	ifs >> pParameter->ThermalStep;
	ifs >> pParameter->NegativeCalculationStep;
	ifs >> pParameter->OutputStep;

	ifs >> pParameter->LambdaID;

	ifs >> pParameter->SORJacobianRelaxationCoefficient;
	ifs >> pParameter->SORLoopMax;
	ifs >> pParameter->SORConvergenceConstant;

	std::string OnOff("");
	ifs >> OnOff;
	if(0 == OnOff.compare("ON"))
	{
		pParameter->bUpdateTemperature = true;
	}
	else
	{
		pParameter->bUpdateTemperature = false;
	}
	if(!ifs.eof())
	{
		ifs >> pParameter->ViscosityUpdateParameterA; //31
		ifs >> pParameter->ViscosityUpdateParameterB; //32
		
		/*ifs >> pParameter->FixedWallTemperature; //33
		ifs >> pParameter->MovableWallTemperature; //34
		*/

		ifs >> pParameter->MaxWallTemperature; //33
		ifs >> pParameter->MinWallTemperature; //34
	}
	if(!ifs.eof())
	{
		ifs >> pParameter->CourantNumber;
	}

	if(!ifs.eof())
	{
		ifs >> pParameter->InitialTimeForZeroAngle;      //Paramter For Rotation (Initial Time for Angle Zero)
		ifs >> pParameter->FinalTimeToReachAngle;      //Parameter for Rotation (Final Time To Reach Angle)

		/*ifs >> pParameter->CellSize.x;      
		ifs >> pParameter->CellSize.y;  */    

		ifs >> pParameter->CellSize.z;
	}	
	if(!ifs.eof())
	{
		ifs >> pParameter->MvVelocity.x;
		ifs >> pParameter->MvVelocity.y;
		ifs >> pParameter->MvVelocity.z;
	}

	if(!ifs.eof())
	{
		ifs >> pParameter->OffSet;
	}

	if(!ifs.eof())
	{
		ifs >> pParameter->ParticleGenerationStep;
		ifs >> pParameter->TotalGenerationStep;
	}
	pParameter->Alpha = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->Alpha;
	}
	pParameter->AlphaLightWeight = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->AlphaLightWeight;
	}
	pParameter->bExplicitCollisonCheck = false;
	pParameter->bExplicitCollisionOut = false;
	pParameter->bImplicitCollisonCheck = false;
	pParameter->CollisionCheckCoefficient = 0.0;
	pParameter->MinDistanceCoefficient = 0.0;
	if(!ifs.eof())
	{
		std::string OnOff("");
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bExplicitCollisonCheck = true;
		}
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bExplicitCollisionOut = true;
		}
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bImplicitCollisonCheck = true;
		}
	}
	if(!ifs.eof())
	{
		ifs >> pParameter->CollisionCheckCoefficient;
		ifs >> pParameter->MinDistanceCoefficient;
	}
	pParameter->bFreeSurfacePressureUpdate = false;
	pParameter->BasePressure = 0.0;
	pParameter->MaxPressure = 0.0;
	pParameter->Time1 = pParameter->Dt;
	pParameter->TriggerStep = 0;
	pParameter->EndStep = 0;
	pParameter->AnalysisArea = 0;
	if(!ifs.eof())
	{
		std::string OnOff("");
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bFreeSurfacePressureUpdate = true;
		}		
		ifs >> pParameter->BasePressure ;
		ifs >> pParameter->MaxPressure ;
		ifs >> pParameter->Time1;
		ifs >> pParameter->TriggerStep ;
		ifs >> pParameter->EndStep ;
		ifs >> pParameter->AnalysisArea ;
	}
	pParameter->bExplicitOut = false;
	if(!ifs.eof())
	{
		std::string OnOff("");
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bExplicitOut = true;
		}
	}
	pParameter->WallFrictionCoeffieicnt = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->WallFrictionCoeffieicnt; //60
	}

	pParameter->DragCoefficient = 0;
	if(!ifs.eof())
	{
		ifs >> pParameter->DragCoefficient;
	}
	pParameter->bTwoPhase = false;
	if(!ifs.eof())
	{
		std::string OnOff("");
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bTwoPhase = true;
		}
	}
	pParameter->StopWallStep = INT_MAX;
	if(!ifs.eof())
	{
		ifs >> pParameter->StopWallStep;
	}
	pParameter->Tolerance = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->Tolerance;
	}
	pParameter->ImplicitSurfaceThreshold = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->ImplicitSurfaceThreshold;
	}
	pParameter->bCGSolver = false;
	if(!ifs.eof())
	{
		std::string OnOff("");
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bCGSolver = true;
		}

	}

	pParameter->bUpdateInnerPressure = false;
	if(!ifs.eof())
	{
		std::string OnOff("");
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bUpdateInnerPressure = true;
		}
	}
	pParameter->SurfaceInfluenceRadiusCoefficient = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->SurfaceInfluenceRadiusCoefficient;
	}
	pParameter->SurfaceThreshold = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->SurfaceThreshold;
	}
	pParameter->GradientForce = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->GradientForce;
	}
	pParameter->StartStep = 0;
	if(!ifs.eof())
	{
		ifs >> pParameter->StartStep;
	}

	pParameter->InitialForce = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->InitialForce;
	}

	pParameter->VelocityThreshold = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->VelocityThreshold;
	}

	pParameter->bExplicitPressure = false;
	if(!ifs.eof())
	{
		std::string OnOff("");
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bExplicitPressure = true;
		}
	}
	pParameter->SonicVelocity = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->SonicVelocity;
	}
	pParameter->ArtificialPressure = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->ArtificialPressure;
	}
	pParameter->r_coeff = 0.0;
	if(!ifs.eof())
	{
		ifs >> pParameter->r_coeff;
	}
	pParameter->bImplicitVelocity = false;
	if(!ifs.eof())
	{
		std::string OnOff("");
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bImplicitVelocity = true;
		}
	}

	pParameter->bNJPCGSolver = false;
	if(!ifs.eof())
	{
		std::string OnOff("");
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bNJPCGSolver = true;
		}
	}

	if(!ifs.eof())
	{
		pParameter->PressureMax = FLT_MAX;
		if(!ifs.eof())
		{
			ifs >> pParameter->PressureMax;
		}
	}

	if(!ifs.eof())
	{
		pParameter->RLimitCoefficient = FLT_MAX;
		if(!ifs.eof())
		{
			ifs >> pParameter->RLimitCoefficient;
		}
	}

	pParameter->bTurbulance = false;
	Scalar SmagorinskyConstant = 0.1;
	Scalar FilterWidth = pParameter->InitialParticleDistance;
	if(!ifs.eof())
	{
		std::string OnOff("");
		ifs >> OnOff;
		if(0 == OnOff.compare("ON"))
		{
			pParameter->bTurbulance = true;
		}
		ifs >> pParameter->SmagorinskyConstant;
		ifs >> pParameter->FilterWidth;
	}

	//Drag permiability Constant 85
	if(!ifs.eof())
	{
		pParameter->DragPermiabilityConstant = 1;
		if(!ifs.eof())
		{
			ifs >> pParameter->DragPermiabilityConstant;
		}
	}

	ifs.close();

	//Interparticle distance information initialization
	if(2 == pParameter->Dimension)
	{
		pParameter->MaxNeighborhoodNum = 81;
	}
	else
	{
		pParameter->MaxNeighborhoodNum = 729;
	}

	return CCT_NOERR;
}

// Calculate the Normal of the Triangle
void CalcTriangleNormal(CTriangle* pTriangle)
{
	Scalar3 Vec0 = pTriangle->Vertex[1] - pTriangle->Vertex[0];
	Scalar3 Vec1 = pTriangle->Vertex[2] - pTriangle->Vertex[0];
	Scalar3 Normal;
	
	Normal.x = Vec0.y * Vec1.z - Vec1.y * Vec0.z;
	Normal.y = Vec0.z * Vec1.x - Vec1.z * Vec0.x;
	Normal.z = Vec0.x * Vec1.y - Vec1.x * Vec0.y;
	
	const Scalar Abs = sqrt_new(Normal.x * Normal.x + Normal.y * Normal.y + Normal.z * Normal.z);

	if(Abs != 0.0)
	{
		pTriangle->Normal.x = Normal.x / Abs;
		pTriangle->Normal.y = Normal.y / Abs;
		pTriangle->Normal.z = Normal.z / Abs;
	}
}


CCTStatusType IO::Load(const std::vector<std::string> &Models, std::vector<Integer>& TrianglesPerModel,Integer& TotalTriangles ,CTriangle*& pTriangle)
{
	TrianglesPerModel.resize(Models.size());
	std::string dump;
	//std::vector<CTriangle> vTriangle;
	Serializer::CTriangleContainer<float> vTriangle;
	Integer count = 0;
	for(size_t i = 0; i < Models.size(); ++i)
	{
		std::string FileName = Models[i];
		std::ifstream ifs(FileName.c_str());
		if(false == ifs.is_open())
		{
			printf("ERR: STL file load");("ERR: STL file load");
			return CCT_FILEERR;
		}		
		ifs >> dump;
		bool IsAscii = 0 == dump.compare("solid") ? true : false;
		ifs.seekg(std::ios::beg);
		if(!Serializer::CSTLSerializer::Load(FileName.c_str(), vTriangle, IsAscii))
		{
			printf("ERR: STL file load");("ERR: STL file load");
			return CCT_FILEERR;
		}
		TrianglesPerModel[i] = vTriangle.size() - count;
		count = vTriangle.size();
		ifs.close();
	}
	
	TotalTriangles = vTriangle.size();
	if(TotalTriangles > 0)
	{
		pTriangle = new CTriangle[ vTriangle.size()];
		for(Integer i = 0; i <  TotalTriangles; ++i)
		{
			pTriangle[i].Normal.x = vTriangle[i]->m_Normal.m_X;
			pTriangle[i].Normal.y = vTriangle[i]->m_Normal.m_Y;
			pTriangle[i].Normal.z = vTriangle[i]->m_Normal.m_Z;
			for(Integer j = 0; j < 3; ++j)
			{
				pTriangle[i].Vertex[j].x = vTriangle[i]->m_Coords[j].m_X;
				pTriangle[i].Vertex[j].y = vTriangle[i]->m_Coords[j].m_Y;
				pTriangle[i].Vertex[j].z = vTriangle[i]->m_Coords[j].m_Z;
				//pTriangle[i].Velocity.x  = 0;//Vector3D(0.0, 0.0, 0.0);
				//pTriangle[i].Velocity.y  = 0;
				//pTriangle[i].Velocity.z  = 0;
				//pTriangle[i].Temperature = 0.0;

			}	
			CalcTriangleNormal(&pTriangle[i]);
		}
	}
	return CCT_NOERR;
}
CCTStatusType IO::Load(const std::string& FileName, CDragTriangle*& pDragTriangles, Integer& TotalTriangleNum, std::string& InnerPressureParameterFileName)
{
	std::ifstream DragFile(FileName.c_str());
	if(!DragFile.is_open())
	{
		std::cout<< "Drag File Load failed";
		return CCT_FILEERR;
	}
	
	std::string InnerPressureParameter;
	std::vector<std::string> vDragModelFileName;
	//std::vector<std::string> vInnerPressureParamter;
	std::vector<CDragTriangle> vDragParameters;
	std::string RootPath;
	DragFile >> RootPath;
	if(RootPath[RootPath.size()-1] != '\\' )
	{
		RootPath.append("\\");
	}
	DragFile >> InnerPressureParameter;
	while(!DragFile.eof())
	{		
		std::string DragModelFileName;
		CDragTriangle DragModelParameter;
		DragFile >> DragModelFileName;
		DragFile >> DragModelParameter.DragLength;
		DragFile >> DragModelParameter.StartStep01;
		DragFile >> DragModelParameter.StartStep02;
		DragFile >> DragModelParameter.EndStep;
		DragFile >> DragModelParameter.InitialAcc.x;
		DragFile >> DragModelParameter.InitialAcc.y;
		DragFile >> DragModelParameter.InitialAcc.z;
		DragFile >> DragModelParameter.GradAcc.x;
		DragFile >> DragModelParameter.GradAcc.y;
		DragFile >> DragModelParameter.GradAcc.z;	
		if(DragFile.good())
		{
			vDragModelFileName.push_back(DragModelFileName);
			vDragParameters.push_back(DragModelParameter);
		}
	}
	Serializer::CTriangleContainer<float> vTriangle;
	std::vector<size_t> TrianglesPerModel;
	size_t count = 0;
	//Inner Pressure Calculation Starts (19 October 2011)
	/*std::string */InnerPressureParameterFileName = RootPath + InnerPressureParameter;
	/*std::ifstream InnerPressureParameterFile(InnerPressureParameterFileName.c_str());
	if(InnerPressureParameterFile.is_open())
	{
		InnerPressureParameterFile >> pParameter->TankVolume;
		InnerPressureParameterFile >> pParameter->TankLimit;
		InnerPressureParameterFile >> pParameter->Gamma;
		InnerPressureParameterFile >> pParameter->BreatherRadius;
		InnerPressureParameterFile >> pParameter->LBreather;
		InnerPressureParameterFile >> pParameter->JointRadius;
		InnerPressureParameterFile >> pParameter->PAtmospheric;
		InnerPressureParameterFile >> pParameter->JointZeta;
		InnerPressureParameterFile >> pParameter->PressureCoefficient;
		InnerPressureParameterFile >> pParameter->Vv_coef;
		InnerPressureParameterFile >> pParameter->kinematicViscosity;
	}
	InnerPressureParameterFile.close();*/
	//Inner Pressure Calculation Ends
	for( int i = 0; i < vDragModelFileName.size(); ++i)
	{
		std::string DragModelFileName = RootPath +  vDragModelFileName[i];
		std::ifstream DragModelFile(DragModelFileName.c_str());
		if(DragModelFile.is_open())
		{
			std::string dump;
			DragModelFile >> dump;
			bool IsAscii = (0 == dump.compare("solid") ? true : false);
			DragModelFile.seekg(std::ios::beg);
			DragModelFile.close();
			if(!Serializer::CSTLSerializer::Load(DragModelFileName, vTriangle, IsAscii))
			{
				std::cout<< "ERR:" << DragModelFileName << "Couldn't load";
				return CCT_FILEERR;
			}
			TrianglesPerModel.push_back(vTriangle.size() - count);
			count = vTriangle.size();			
		}
		else
		{
			std::cout<< "ERR:" << DragModelFileName << "Couldn't load";
		}
	}
	TotalTriangleNum = vTriangle.size();
	if(TotalTriangleNum > 0)
	{
		pDragTriangles = new CDragTriangle[ TotalTriangleNum];
		int tid = 0;
		for(Integer i = 0; i < TrianglesPerModel.size(); ++i)
		{
			for(Integer j = 0; j < TrianglesPerModel[i]; ++j)
			{
				pDragTriangles[tid].Normal.x = vTriangle[tid]->m_Normal.m_X;
				pDragTriangles[tid].Normal.y = vTriangle[tid]->m_Normal.m_Y;
				pDragTriangles[tid].Normal.z = vTriangle[tid]->m_Normal.m_Z;
				for(Integer k = 0; k < 3; ++k)
				{
					pDragTriangles[tid].Vertex[k].x = vTriangle[tid]->m_Coords[k].m_X;
					pDragTriangles[tid].Vertex[k].y = vTriangle[tid]->m_Coords[k].m_Y;
					pDragTriangles[tid].Vertex[k].z = vTriangle[tid]->m_Coords[k].m_Z;
					pDragTriangles[tid].DragLength = vDragParameters[i].DragLength;
					pDragTriangles[tid].StartStep01 = vDragParameters[i].StartStep01;
					pDragTriangles[tid].StartStep02 = vDragParameters[i].StartStep02;
					pDragTriangles[tid].EndStep = vDragParameters[i].EndStep;
					pDragTriangles[tid].InitialAcc = vDragParameters[i].InitialAcc;
					pDragTriangles[tid].GradAcc = vDragParameters[i].GradAcc;
				}
				////Inner Pressure Calculation Starts 19 October 2011
				pDragTriangles[tid].ID = tid;
				//Inner Pressure Calculation Ends 19 October 2011
				++tid;
			}
		}
	}
	else
	{
		std::cout<< " Zero Drag Triangles"<<std::endl;
		return CCT_ETCERR;
	}
	return CCT_NOERR;
}
CCTStatusType IO::Save(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle, const std::vector<Scalar3>& Position, Scalar CurrentTime, bool isAscii,
					   const Integer* particleID,	const Scalar3* phParticlePosition, const Scalar3* particleVelocity, const Scalar* pressure, const Scalar* density, const Scalar* temperature,
					   const Scalar* kineticViscosity,const  Scalar* solidPhaseRate, const ParticleType* type,const std::vector<Scalar>& RotationAngle)
{
	CCTStatusType StatusType = CCT_NOERR;
	if(isAscii)
	{
		StatusType = SaveAscii(FileName, ParticleNum, phParticle,Position,CurrentTime,
								particleID,	phParticlePosition, particleVelocity, pressure, density, temperature,
								kineticViscosity, solidPhaseRate, type, RotationAngle);
	}
	else
	{
		//StatusType = SaveBinary(FileName, ParticleNum, phParticle,Position,CurrentTime);
		StatusType = SaveBinary(FileName, ParticleNum, phParticle,Position,CurrentTime,
								particleID,	phParticlePosition, particleVelocity, pressure, density, temperature,
								kineticViscosity, solidPhaseRate, type, RotationAngle);
	}	
	return StatusType;
}

CCTStatusType IO::SaveAscii(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle, const std::vector<Scalar3>& Position,const Scalar CurrentTime,
							const Integer* particleID,	const Scalar3* phParticlePosition, const Scalar3* particleVelocity, const Scalar* pressure, const Scalar* density, const Scalar* temperature,
							const Scalar* kineticViscosity,const  Scalar* solidPhaseRate, const ParticleType* type,const std::vector<Scalar>& RotationAngle)
{
	std::ofstream fout(FileName.c_str());
	if(false == fout.is_open())
	{
		return CCT_FILEERR;
	}
	//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
	size_t ModelNumber = Position.size();	
	fout << ParticleNum << std::endl;
	for(Integer i = 0;i < ParticleNum; ++i)
	{
		fout << phParticlePosition[i].x << " "
			 << phParticlePosition[i].y << " "
			 << phParticlePosition[i].z << " "
			 << particleVelocity[i].x << " "
			 << particleVelocity[i].y << " "
			 << particleVelocity[i].z << " "

			 << pressure[i] << " "
			 << density[i] << " "
			 << temperature[i] << " "

			 << kineticViscosity[i] << " "

			 << solidPhaseRate[i] << " "		
			 << type[i] << "\n";
	}
	fout << CurrentTime << std::endl;
	fout << ModelNumber << std::endl;
	for(size_t i=0;i<Position.size();i++)
	{
		fout << Position[i].x << "\t" << Position[i].y << "\t" << Position[i].z << "\t" << RotationAngle[i] << std::endl;
	}
	//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------

	fout.close();
	return CCT_NOERR;
}

//CCTStatusType IO::SaveBinary(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle, const std::vector<Scalar3>& Position, Scalar CurrentTime)
CCTStatusType IO::SaveBinary(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle, const std::vector<Scalar3>& Position,const Scalar CurrentTime,
							const Integer* particleID,	const Scalar3* phParticlePosition, const Scalar3* particleVelocity, const Scalar* pressure, const Scalar* density, 
							const Scalar* temperature,	const Scalar* kineticViscosity,const  Scalar* solidPhaseRate, const ParticleType* type,const std::vector<Scalar>& RotationAngle)
{
	std::ofstream fout(FileName.c_str(),std::ofstream::binary);
	if(false == fout.is_open())
	{
		return CCT_FILEERR;
	}

	// Modified, E&T Nepal
	Integer ModelNumber = Position.size();
	
	fout.write((char*)&ParticleNum,sizeof(Integer));	

	//fout.write((char*)phParticle,sizeof(CParticle)*ParticleNum);
	for(Integer j=0;j<ParticleNum;++j)
	{
		fout.write((char*)&phParticlePosition[j],sizeof(Scalar3));
		fout.write((char*)&particleVelocity[j],sizeof(Scalar3));
		fout.write((char*)&pressure[j],sizeof(Scalar));
		fout.write((char*)&density[j],sizeof(Scalar));
		fout.write((char*)&temperature[j],sizeof(Scalar));
		fout.write((char*)&kineticViscosity[j],sizeof(Scalar));
		fout.write((char*)&solidPhaseRate[j],sizeof(Scalar));
		fout.write((char*)&type[j],sizeof(ParticleType));
	}

	fout.write((char*)&CurrentTime,sizeof(Scalar));
	
	fout.write((char*)&ModelNumber,sizeof(Integer));

	for(size_t i=0;i<Position.size();i++)
	{		
		fout.write((char*)&	Position[i],sizeof(Scalar3));
		fout.write((char*)&	RotationAngle[i],sizeof(Scalar));
		std::cout << "\n" << Position[i].x << "\t" << Position[i].y << "\t" << Position[i].z << "\t" << RotationAngle[i] <<"\n";
	}	
	// End, E&T Nepal
	
	fout.close();
   
	return CCT_NOERR;
}
// GraphR
CCTStatusType IO::Dump(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle)
{
	std::ofstream ofs(FileName.c_str());		
	if(false == ofs.is_open())
	{
		return CCT_FILEERR;
	}

	ofs << "データ形式,23" << std::endl << std::endl;
	Scalar MinX = FLT_MAX, MinY = FLT_MAX, MinZ = FLT_MAX;
	Scalar MaxX = -FLT_MAX, MaxY = -FLT_MAX, MaxZ = - FLT_MAX;

	for(Integer i = 0; i < ParticleNum; ++i)
	{	
		ofs << phParticle[i].Position.x << "," 
			<< phParticle[i].Position.y << "," 
			<< phParticle[i].Position.z << "," 
			<<"0.5," << "0" << std::endl;

		if(phParticle[i].Position.x > MaxX)
		{
			MaxX = phParticle[i].Position.x;				
		}
		if(phParticle[i].Position.x <= MinX)
		{
			MinX = phParticle[i].Position.x;				
		}
		if(phParticle[i].Position.y > MaxY)
		{
			MaxY = phParticle[i].Position.y;			
		}
		if(phParticle[i].Position.y <= MinY)
		{
			MinY = phParticle[i].Position.y;				
		}
		if(phParticle[i].Position.z > MaxZ)
		{
			MaxZ = phParticle[i].Position.z;				
		}
		if(phParticle[i].Position.z <= MinZ)
		{
			MinZ = phParticle[i].Position.z;				
		}
	}

	ofs << "\n\n[ScaleDraw]\nDrawFlg=1\nSameDataScaleFlg=0\nValueFontSize=12\nValueFontColor=0.000000 0.000000 0.000000\n";
	ofs << "LabelFontSize=12\nLabelFontColor=0.000000 0.000000 0.000000\n[ScaleX]\nScaleMin=" << MinX << "\nScaleMax=" << MaxX << "\n";
	ofs << "ScaleStep=0.100000\nScaleSubStep=0.050000\nSize=150.000000\nLabel=\nScaleColor=0.498039 0.498039 0.498039\n";
	ofs << "ScaleSubColor=0.698039 0.698039 0.698039\nValueFontPrintValueFiguret=3\nValueFontPrintValueType=0\nReverseFlg=0\n";
	ofs << "LogScaleFlg=0\nDrawSubScaleFlg=0\nAutoScaleFlg=0\n[ScaleY]\nScaleMin=" << MinY << "\nScaleMax=" << MaxY << "\n";
	ofs << "ScaleStep=0.060000\nScaleSubStep=0.030000\nSize=100.000000\nLabel=\nScaleColor=0.498039 0.498039 0.498039\n";
	ofs << "ScaleSubColor=0.698039 0.698039 0.698039\nValueFontPrintValueFiguret=3\nValueFontPrintValueType=0\n";
	ofs << "ReverseFlg=0\nLogScaleFlg=0\nDrawSubScaleFlg=0\nAutoScaleFlg=0\n";

	ofs.close();

	return CCT_NOERR;
}

CCTStatusType IO::SavePressures(const std::string& FileName, const Scalar P1, const Scalar P2 , const Scalar P3, const Scalar H_Zoint, const Scalar H_breather)
{
	//if(P2 == )
	std::ofstream ofs(FileName.c_str(),std::ios::app);		
	if(false == ofs.is_open())
	{
		return CCT_FILEERR;
	}
	ofs << P1 << "\t" << P2 << "\t" << P3 <<"\t" <<H_Zoint <<"\t" <<H_breather <<"\n";
	
	ofs.close();

	return CCT_NOERR;
}
CCTStatusType IO::SaveHeight(const std::string& FileName, const Scalar Time, const Scalar Height)
{
	//if(P2 == )
	std::ofstream ofs(FileName.c_str(),std::ios::app);		
	if(false == ofs.is_open())
	{
		return CCT_FILEERR;
	}
	ofs << Time << "\t" << Height << "\n";
	
	ofs.close();

	return CCT_NOERR;
}
void IO::saveParameterWithName(const std::string& FileName, const CParameter * pParameter,
							   Integer ModelNumber, const std::vector<Scalar3>& CenterOfRotation, const std::vector<Scalar3>& OtherPointOfRotation, const std::vector<Integer>& RotationRPM ,
							   const std::vector<Scalar>& WallFrictionCoefficient ,const std::vector<Scalar>& WallThermalConductivity ,const std::vector<Scalar>& WallThermalResistivity ,
							   const DragParameter * DragPrameter, Integer dragTriangleNum)
{
	/*std::stringstream fileName;
	fileName << "test\\Parameter.Dat";
	std::ofstream ofParameter(fileName.str().c_str());*/
	std::ofstream ofParameter(FileName.c_str());
	if(ofParameter)
	{
		ofParameter << "pParameter->Dimension \t"		<<	pParameter->Dimension <<"\n";
		ofParameter << "pParameter->Gravity Constant.x\t"<<	pParameter->GravityConstant.x <<"\n"; 
		ofParameter << "pParameter->GravityConstant.y; "	<<	pParameter->GravityConstant.y	<< "\n"; 
		ofParameter << "pParameter->GravityConstant.z;"	<<	pParameter->GravityConstant.z	<< "\n";
		ofParameter << "pParameter->Density;"	<<	pParameter->Density	<< "\n";
		ofParameter << "pParameter->ConstantParticleDensity;"	<<	pParameter->ConstantParticleDensity	<< "\n";
		ofParameter << "pParameter->SurfaceTensionCoefficient;"	<<	pParameter->SurfaceTensionCoefficient	<< "\n";
		ofParameter << "pParameter->InitialParticleDistance;"	<<	pParameter->InitialParticleDistance	<< "\n";
		ofParameter << "pParameter->FreeSurfaceThreshold;"	<<	pParameter->FreeSurfaceThreshold	<< "\n";
		ofParameter << "pParameter->GradientInfluenceRadiusCoefficient;"	<<	pParameter->GradientInfluenceRadiusCoefficient	<< "\n";
		ofParameter << "pParameter->LaplacianInfluenceRadiusCoefficient;"	<<	pParameter->LaplacianInfluenceRadiusCoefficient	<< "\n";
		ofParameter << "pParameter->TemperatureInfluenceRadiusCoefficient;"	<<	pParameter->TemperatureInfluenceRadiusCoefficient	<< "\n";
		ofParameter << "pParameter->FruidInitializeialTemperature;"	<<	pParameter->FruidInitializeialTemperature	<< "\n";
		ofParameter << "pParameter->WallInitializeialTemperature;"	<<	pParameter->WallInitializeialTemperature	<< "\n";
		ofParameter << "pParameter->FruidThermalResistance;"	<<	pParameter->FruidThermalResistance	<< "\n";
		ofParameter << "pParameter->WallThermalResistance;"	<<	pParameter->WallThermalResistance	<< "\n";
		ofParameter << "pParameter->FruidThermalConductivity;"	<<	pParameter->FruidThermalConductivity	<< "\n";
		ofParameter << "pParameter->WallThermalConductivity;"	<<	pParameter->WallThermalConductivity	<< "\n";
		ofParameter << "pParameter->FruidSpecificHeat;"	<<	pParameter->FruidSpecificHeat	<< "\n";
		ofParameter << "pParameter->FruidDensity;"	<<	pParameter->FruidDensity	<< "\n";
		ofParameter << "pParameter->Dt;"	<<	pParameter->Dt	<< "\n";
		ofParameter << "pParameter->AnalysisStep;"	<<	pParameter->AnalysisStep	<< "\n";
		ofParameter << "pParameter->ThermalStep;"	<<	pParameter->ThermalStep	<< "\n";
		ofParameter << "pParameter->NegativeCalculationStep;"	<<	pParameter->NegativeCalculationStep	<< "\n";
		ofParameter << "pParameter->OutputStep;"	<<	pParameter->OutputStep	<< "\n";
		ofParameter << "pParameter->LambdaID;"	<<	pParameter->LambdaID	<< "\n";
		ofParameter << "pParameter->SORJacobianRelaxationCoefficient;"	<<	pParameter->SORJacobianRelaxationCoefficient	<< "\n";
		ofParameter << "pParameter->SORLoopMax;"	<<	pParameter->SORLoopMax	<< "\n";
		ofParameter << "pParameter->SORConvergenceConstant;"	<<	pParameter->SORConvergenceConstant	<< "\n";
		ofParameter << "pParameter->bUpdateTemperature;"	<<	pParameter->bUpdateTemperature	<< "\n";
		ofParameter << "pParameter->ViscosityUpdateParameterA;"	<<	pParameter->ViscosityUpdateParameterA	<< "\n";
		ofParameter << "pParameter->ViscosityUpdateParameterB;"	<<	pParameter->ViscosityUpdateParameterB	<< "\n";
		ofParameter << "pParameter->MaxWallTemperature;"	<<	pParameter->MaxWallTemperature	<< "\n";
		ofParameter << "pParameter->MinWallTemperature;"	<<	pParameter->MinWallTemperature	<< "\n";
		ofParameter << "pParameter->CourantNumber;"	<<	pParameter->CourantNumber	<< "\n";
		ofParameter << "pParameter->InitialTimeForZeroAngle;"	<<	pParameter->InitialTimeForZeroAngle	<< "\n";
		ofParameter << "pParameter->FinalTimeToReachAngle;"	<<	pParameter->FinalTimeToReachAngle	<< "\n";
		ofParameter << "pParameter->CellSize.z;"	<<	pParameter->CellSize.z	<< "\n";
		ofParameter << "pParameter->MvVelocity.x;"	<<	pParameter->MvVelocity.x	<< "\n";
		ofParameter << "pParameter->MvVelocity.y;"	<<	pParameter->MvVelocity.y	<< "\n";
		ofParameter << "pParameter->MvVelocity.z;"	<<	pParameter->MvVelocity.z	<< "\n";
		ofParameter << "pParameter->OffSet;"	<<	pParameter->OffSet	<< "\n";
		ofParameter << "pParameter->ParticleGenerationStep;"	<<	pParameter->ParticleGenerationStep	<< "\n";
		ofParameter << "pParameter->TotalGenerationStep;"	<<	pParameter->TotalGenerationStep	<< "\n";
		ofParameter << "pParameter->Alpha;"	<<	pParameter->Alpha	<< "\n";
		ofParameter << "pParameter->AlphaLightWeight;"	<<	pParameter->AlphaLightWeight	<< "\n";
		ofParameter << "pParameter->bExplicitCollisonCheck = true;"	<<	pParameter->bExplicitCollisonCheck	<< "\n";
		ofParameter << "pParameter->bExplicitCollisionOut = true;"	<<	pParameter->bExplicitCollisionOut << "\n";
		ofParameter << "pParameter->bImplicitCollisonCheck = true;"	<<	pParameter->bImplicitCollisonCheck<< "\n";
		ofParameter << "pParameter->CollisionCheckCoefficient;"	<<	pParameter->CollisionCheckCoefficient	<< "\n";
		ofParameter << "pParameter->MinDistanceCoefficient;"	<<	pParameter->MinDistanceCoefficient	<< "\n";
		ofParameter << "pParameter->bFreeSurfacePressureUpdate = true;"	<<	pParameter->bFreeSurfacePressureUpdate	<< "\n";
		ofParameter << "pParameter->BasePressure ;"	<<	pParameter->BasePressure 	<< "\n";
		ofParameter << "pParameter->MaxPressure ;"	<<	pParameter->MaxPressure 	<< "\n";
		ofParameter << "pParameter->Time1;"	<<	pParameter->Time1	<< "\n";
		ofParameter << "pParameter->TriggerStep ;"	<<	pParameter->TriggerStep 	<< "\n";
		ofParameter << "pParameter->EndStep ;"	<<	pParameter->EndStep 	<< "\n";
		ofParameter << "pParameter->AnalysisArea ;"	<<	pParameter->AnalysisArea 	<< "\n";
		ofParameter << "pParameter->bExplicitOut = true;"	<<	pParameter->bExplicitOut << "\n";
		ofParameter << "pParameter->WallFrictionCoeffieicnt;"	<<	pParameter->WallFrictionCoeffieicnt	<< "\n";
		ofParameter << "pParameter->DragCoefficient;"	<<	pParameter->DragCoefficient	<< "\n";
		ofParameter << "pParameter->bTwoPhase = true;"	<<	pParameter->bTwoPhase	<< "\n";
		ofParameter << "pParameter->StopWallStep;"	<<	pParameter->StopWallStep	<< "\n";
		ofParameter << "pParameter->Tolerance;"	<<	pParameter->Tolerance	<< "\n";
		ofParameter << "pParameter->ImplicitSurfaceThreshold;"	<<	pParameter->ImplicitSurfaceThreshold	<< "\n";
		ofParameter << "pParameter->bCGSolver = true;"	<<	pParameter->bCGSolver	<< "\n";
		ofParameter << "pParameter->bUpdateInnerPressure = true;"	<<	pParameter->bUpdateInnerPressure << "\n";
		ofParameter << "pParameter->SurfaceInfluenceRadiusCoefficient;"	<<	pParameter->SurfaceInfluenceRadiusCoefficient	<< "\n";
		ofParameter << "pParameter->SurfaceThreshold;"	<<	pParameter->SurfaceThreshold	<< "\n";
		ofParameter << "pParameter->GradientForce;"	<<	pParameter->GradientForce	<< "\n";
		ofParameter << "pParameter->StartStep;"	<<	pParameter->StartStep	<< "\n";
		ofParameter << "pParameter->InitialForce;"	<<	pParameter->InitialForce	<< "\n";
		ofParameter << "pParameter->VelocityThreshold;"	<<	pParameter->VelocityThreshold	<< "\n";
		ofParameter << "pParameter->bExplicitPressure = true;"	<<	pParameter->bExplicitPressure << "\n";
		ofParameter << "pParameter->SonicVelocity;"	<<	pParameter->SonicVelocity	<< "\n";
		ofParameter << "pParameter->ArtificialPressure;"	<<	pParameter->ArtificialPressure	<< "\n";
		ofParameter << "pParameter->r_coeff;"	<<	pParameter->r_coeff	<< "\n";
		ofParameter << "pParameter->bImplicitVelocity = true;"	<<	pParameter->bImplicitVelocity << "\n";
		ofParameter << "pParameter->bNJPCGSolver = true;"	<<	pParameter->bNJPCGSolver << "\n";
		ofParameter << "pParameter->PressureMax;"	<<	pParameter->PressureMax	<< "\n";
		ofParameter << "pParameter->RLimitCoefficient;"	<<	pParameter->RLimitCoefficient	<< "\n";
		ofParameter << "pParameter->bTurbulance = true;"	<<	pParameter->bTurbulance 	<< "\n";
		ofParameter << "pParameter->SmagorinskyConstant;"	<<	pParameter->SmagorinskyConstant	<< "\n";
		ofParameter << "pParameter->FilterWidth;"	<<	pParameter->FilterWidth	<< "\n";
		ofParameter <<  "pParameter->MaxNeighborhoodNum;"	<< pParameter->MaxNeighborhoodNum <<"\n";	
		ofParameter <<  "pParameter->OutputStep;"	<< pParameter->OutputStep <<"\n";
		//Drag@Rajan
		ofParameter << "pParameter->DragPermiabilityConstant" << pParameter->DragPermiabilityConstant << "\n";
		
		for(int i = 0 ; i < ModelNumber ; ++i)
		{
			Scalar3 centerOfRotaion = CenterOfRotation[i];
			Scalar3 OtherPointOfRotaion = OtherPointOfRotation[i];
			Integer RotationinRPM = RotationRPM[i];
			Scalar AngleInDegeree = RotationinRPM * pParameter->Dt * 6;

			ofParameter << "pConfiguration->WallFricitonCoefficient ; "<< WallFrictionCoefficient[i] <<"\n";
			ofParameter << "pConfiguration->WallThermalRestivity    ; "<< WallThermalResistivity[i]  <<"\n";
			ofParameter << "pConfiguration->WallThermalConductivity ; "<< WallThermalConductivity[i] <<"\n";

			ofParameter << "pConfiguration->Center of Rotaion(X,Y,Z); ("		<< centerOfRotaion.x	<<","<< centerOfRotaion.y		<<","<<centerOfRotaion.z		<<")\n";
			ofParameter << "pConfiguration->Other Point of Rotaion(X,Y,Z); ("	<< OtherPointOfRotaion.x<<","<< OtherPointOfRotaion.y	<<","<<OtherPointOfRotaion.z	<<")\n";
			ofParameter << "pConfifuration->Rotation(RPM,Angle Degree);("		<< RotationinRPM		<<","<< AngleInDegeree			<<")\n\n";
		}
		ofParameter << "Drag Parameters Values\n";
		for(int i = 0 ; i < dragTriangleNum ; i ++)
		{
			ofParameter <<"Drag Model:ID: "<< DragPrameter[i].DragID		<<", Drag Type:"		   << DragPrameter[i].DragToken				<<", Input Token:"			<< DragPrameter[i].InputVelocityToken1				<<", Output Token:"<< DragPrameter[i].OutputVelocityToken2 <<", Influence Region:"<< DragPrameter[i].CoefficientofInflueneRegion <<"\n";
			ofParameter <<"Velcity Index:" << DragPrameter[i].VelocityIndex	<<", Velocity Exponential:"<< DragPrameter[i].VelocityExponential	<<", Velocity Coefficient:"	<< DragPrameter[i].VelocityDragPermiabilityConstant <<"\n";
			ofParameter <<"Temp Index:"	   << DragPrameter[i].TempIndex		<<", Temp Exponential:"	   << DragPrameter[i].TempExponential		<<", Temp Coefficient:"		<< DragPrameter[i].TempDragPermiabilityConstant		<<"\n";
			ofParameter <<"Constant Vx: "  << DragPrameter[i].ConstantVx	<<", Constant Vy:"		   << DragPrameter[i].ConstantVy			<<", Constant Vz:"			<< DragPrameter[i].ConstantVz						<<"\n";
			ofParameter <<"Velocity Magnitude Factor: "  << DragPrameter[i].VelocityMagnitueFactor <<"\n\n";
		}
		
		ofParameter.close();
	}
}