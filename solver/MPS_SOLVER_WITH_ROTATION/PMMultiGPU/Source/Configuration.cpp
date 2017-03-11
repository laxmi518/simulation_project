#include "stdafx.h"
#include "Configuration.h"
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/lexical_cast.hpp>

CConfiguration::CConfiguration(void)
:m_pComputeZone(NULL)
,m_IsOutputAscii(true)
,m_StartStep(0)
{
	m_Section.x = 0;
	m_Section.y = 0;
	m_Section.z = 0;

}

CConfiguration::~CConfiguration(void)
{
	DestroyAll();
}
CCTStatusType CConfiguration::ReadConfigurationFileWithDrag(const std::string& FileName)
{
	ifstream ConfigureFile;
	ConfigureFile.open(FileName.c_str());

	if(!ConfigureFile.is_open())
	{
		cout<<"Error in opening Configuration File";
		return CCT_FILEERR;
	}
	DestroyAll();

	string Line;
	stringstream LineStream;
	
	getline(ConfigureFile,m_InputPath);
	//cout<<m_InputPath<<endl;

	getline(ConfigureFile,m_ParameterFileName);
	//cout<<m_ParameterFileName<<endl;

	getline(ConfigureFile,m_ParticleFileName);
	//cout<<m_ParticleFileName<<endl;

	getline(ConfigureFile,Line);

	//LineStream.str(Line);

	while(!Line.empty())
	{
		string tempFileName;
		Integer additionInterval;
		Integer totalAdditionSteps;
		LineStream.clear(); //To clear the stream before read next line
		LineStream.str(Line);		
		LineStream >> tempFileName >> additionInterval >> totalAdditionSteps;
		m_aBucketFileName.push_back(tempFileName);
		m_aAdditionInterval.push_back(additionInterval);
		m_aTotalAdditionStep.push_back(totalAdditionSteps);
		getline(ConfigureFile,Line);
	}
	//Earlier Bucket File
	//LineStream>> m_BucketFileName >> m_AdditionIntervalTime >> m_TotalAddition;

	getline(ConfigureFile , Line);
	int CountDragID =0;
	while(!Line.empty())
	{
		string tempFileName;

		LineStream.clear(); //To clear the stream before read next line
		LineStream.str(Line);

		LineStream >> tempFileName;			//Read STL file Name
		m_aSTLFileName.push_back(tempFileName);
		
		//char IfDrag;
		std::string IfDrag;
		LineStream>>IfDrag;
		boost::to_upper(IfDrag);

		//To check where to effect values
		std::string inputToken;
		std::string outPutToken;

		if(IfDrag=="D" || IfDrag=="DE" || IfDrag == "DC" || IfDrag == "DF")
		{
			m_aSTLISDrag.push_back(true);
			m_aSTLDragID.push_back(CountDragID);
			DragParameter objDragParameter;
			objDragParameter.DragID	= CountDragID;
			//objDragParameter.DragToken= IfDrag[0];

			LineStream >> inputToken; //For Input Token
			LineStream >> outPutToken; //For Output Token

			//Convert into uppter Case
			boost::to_upper(inputToken);
			boost::to_upper(outPutToken);

			LineStream >> objDragParameter.CoefficientofInflueneRegion;

			if(IfDrag == "D" || IfDrag == "DE")    //Check whether its exponention or not
			{
				objDragParameter.DragToken	= 'E'; //Confirm its exponential Function

				if( inputToken == "V") //Check whether its Temperature or Velocity
				{
					objDragParameter.InputVelocityToken1	= 'V'; //V is used for Velocity
					objDragParameter.OutputVelocityToken2	= 'V'; //V is used for Velocity

					LineStream >> objDragParameter.VelocityIndex;
					LineStream >> objDragParameter.VelocityExponential;
					LineStream >> objDragParameter.VelocityDragPermiabilityConstant;

					//If its Velocity set Temperature Zero;
					objDragParameter.TempIndex						= 0.0;
					objDragParameter.TempExponential				= 0.0;
					objDragParameter.TempDragPermiabilityConstant	= 0.0;

					objDragParameter.VelocityMagnitueFactor = 0.0;

					std::cout << "INFO : Velocity Exponential Equation : "<< tempFileName << "\n";
				}
				else if(inputToken == "T")
				{
					objDragParameter.InputVelocityToken1	= 'T'; //T is used for Temperature
					objDragParameter.OutputVelocityToken2	= 'T'; //T is used for Temperature
					LineStream >> objDragParameter.TempIndex;
					LineStream >> objDragParameter.TempExponential;
					LineStream >> objDragParameter.TempDragPermiabilityConstant;

					//If its temperature set Velocity to zero
					objDragParameter.VelocityIndex						= 0.0;
					objDragParameter.VelocityExponential				= 0.0;
					objDragParameter.VelocityDragPermiabilityConstant	= 0.0;

					objDragParameter.VelocityMagnitueFactor = 0.0;

					std::cout << "INFO : Temperature Exponential Equation: "<< tempFileName << "\n";

				}
				else if(inputToken == "VT" || "TV")
				{
					objDragParameter.InputVelocityToken1	= 'B'; //B is used for both temperature and Velocity
					objDragParameter.OutputVelocityToken2	= 'B'; //B used for both temperature and Velocity

					LineStream >> objDragParameter.VelocityIndex;
					LineStream >> objDragParameter.VelocityExponential;
					LineStream >> objDragParameter.VelocityDragPermiabilityConstant;

					LineStream >> objDragParameter.TempIndex;
					LineStream >> objDragParameter.TempExponential;
					LineStream >> objDragParameter.TempDragPermiabilityConstant;

					std::cout << "INFO : Velocity and Temperature Exponential Equation :  "<< tempFileName << "\n";
				}
				
				objDragParameter.ConstantVx = 0.0;
				objDragParameter.ConstantVy = 0.0;
				objDragParameter.ConstantVz = 0.0;

				objDragParameter.VelocityMagnitueFactor = 0.0;
			}
			if(IfDrag == "DF")    //Check whether its exponention or not
			{
				objDragParameter.DragToken	= 'F'; //Confirm its exponential Function

				if( inputToken == "V") //Check whether its Temperature or Velocity
				{
					objDragParameter.InputVelocityToken1	= 'V'; //V is used for Velocity
					objDragParameter.OutputVelocityToken2	= 'V'; //V is used for Velocity

					LineStream >> objDragParameter.VelocityMagnitueFactor;

					//if its Constant Vecocity, set Index and exp zero;
					objDragParameter.VelocityIndex						= 0.0;
					objDragParameter.VelocityExponential				= 0.0;
					objDragParameter.VelocityDragPermiabilityConstant	= 0.0;

					//If its Constant Velocity, set Temperature Zero;
					objDragParameter.TempIndex						= 0.0;
					objDragParameter.TempExponential				= 0.0;
					objDragParameter.TempDragPermiabilityConstant	= 0.0;

					objDragParameter.ConstantVx = 0.0;
					objDragParameter.ConstantVy = 0.0;
					objDragParameter.ConstantVz = 0.0;

					std::cout << "INFO : VELOCITY MAGNITUE FACTOR : "<< tempFileName << "\n";
				}
			}
			else if(IfDrag == "DC")//If its not Exponential,than its Constant function.
			{
				objDragParameter.InputVelocityToken1	= 'V'; //V for Velocity Token
				objDragParameter.OutputVelocityToken2	= 'V'; //V for Velocity Token

				objDragParameter.DragToken= 'C';

				LineStream >> objDragParameter.ConstantVx;
				LineStream >> objDragParameter.ConstantVy;
				LineStream >> objDragParameter.ConstantVz;

				//if its Constant Vecocity, set Index and exp zero;
				objDragParameter.VelocityIndex						= 0.0;
				objDragParameter.VelocityExponential				= 0.0;
				objDragParameter.VelocityDragPermiabilityConstant	= 0.0;

				//If its Constant Velocity, set Temperature Zero;
				objDragParameter.TempIndex						= 0.0;
				objDragParameter.TempExponential				= 0.0;
				objDragParameter.TempDragPermiabilityConstant	= 0.0;

				std::cout << "INFO : Constant Velocity Equation : "<< tempFileName << "\n";
			}
			m_aSTLDragParameter.push_back(objDragParameter);
			CountDragID++;
			
			Scalar3 temp = make_Scalar3(0,0,0);		
			m_aSTLVelocity.push_back(temp);	

			Scalar Temperature=533;
			m_aSTLTemperature.push_back(Temperature);

			Scalar3 tempPos = make_Scalar3(0,0,0);
			m_aSTLPosition.push_back(tempPos);

			bool linearMovableOrNot = false;
			m_aSTLRestToOriginal.push_back(linearMovableOrNot);

			Integer ResetInterval = 0;
			m_aSTLRestInterval.push_back(ResetInterval);

			bool RotatingOrNot = false;
			m_aRotatingOrNot.push_back(RotatingOrNot);

			Integer rotationInRPM = 0;
			m_aRPM.push_back(rotationInRPM);

			Scalar startTimeOfRotation = 0;
			m_aStartTimeOfRotation.push_back(startTimeOfRotation);

			Scalar timeToReachRPM = 0;
			m_aTimeToReachRPM.push_back(timeToReachRPM);

			Scalar3 tempCenterOfRotation = make_Scalar3(0,0,0);
			m_aCenterOfRotation.push_back(tempCenterOfRotation);

			Scalar3 tempOtherPointOfRotation = make_Scalar3(0,0,0);
			m_aSecondPointOfRotation.push_back(tempOtherPointOfRotation);

			Scalar tempWallFrictionCoefficient = 0;
			m_aWallFrictionCoefficient.push_back(tempWallFrictionCoefficient);

			Scalar tempWallThermalResistivity = 0;
			m_aWallThermalRestivity. push_back(tempWallThermalResistivity);

			Scalar tempWallThermalConductivity = 0;
			m_aWallThermalConductivity.push_back(tempWallThermalConductivity);
		
		}
		else
		{
			/*DragParameter objDragParameter;
			objDragParameter.CoefficientofInflueneRegion = 0;
			objDragParameter.CoeffOfMember = new float;
			objDragParameter.CoeffOfMember[0] = 0;
			objDragParameter.DragID = 0;
			objDragParameter.DragToken = 'r';
			objDragParameter.InputVelocityToken1 = 'r';
			objDragParameter.OutputVelocityToken2 ='r';
			objDragParameter.RankofPolynomialEquation = 0;
			m_aSTLDragParameter.push_back(objDragParameter);*/

			Scalar3 temp = make_Scalar3(0,0,0);	
			temp.x = boost::lexical_cast<float>(IfDrag);
			LineStream>>temp.y>>temp.z; //Read STL Velocity
			m_aSTLVelocity.push_back(temp);	

			Scalar Temperature;
			LineStream>>Temperature;			//Read STL Temperatue
			m_aSTLTemperature.push_back(Temperature);

			Scalar3 tempPos = make_Scalar3(0,0,0);
			LineStream >> tempPos.x >> tempPos.y >> tempPos.z; //Read STL Offset
			m_aSTLPosition.push_back(tempPos);

			bool linearMovableOrNot = false;
			LineStream >> linearMovableOrNot ; //Read Either Linear Movable or Not
			m_aSTLRestToOriginal.push_back(linearMovableOrNot);

			Integer ResetInterval = 0;
			LineStream >> ResetInterval; //Read Reset interval for Movable wall
			m_aSTLRestInterval.push_back(ResetInterval);

			bool RotatingOrNot = false;
			LineStream >> RotatingOrNot ; //Read Either Rotating or Not Rotating
			m_aRotatingOrNot.push_back(RotatingOrNot);

			Integer rotationInRPM = 0;
			LineStream >> rotationInRPM;		//Read rotation  unit in RPM
			m_aRPM.push_back(rotationInRPM);

			Scalar startTimeOfRotation = 0;
			LineStream >> startTimeOfRotation; //Read Start time of Rotation
			m_aStartTimeOfRotation.push_back(startTimeOfRotation);

			Scalar timeToReachRPM = 0;
			LineStream >> timeToReachRPM;		//Read Final time to Reach RPM
			m_aTimeToReachRPM.push_back(timeToReachRPM);

			Scalar3 tempCenterOfRotation = make_Scalar3(0,0,0);
			LineStream >> tempCenterOfRotation.x >> tempCenterOfRotation.y >> tempCenterOfRotation.z; //Read Center of Rotation
			m_aCenterOfRotation.push_back(tempCenterOfRotation);

			Scalar3 tempOtherPointOfRotation = make_Scalar3(0,0,0);
			LineStream >> tempOtherPointOfRotation.x >> tempOtherPointOfRotation.y >> tempOtherPointOfRotation.z;  //Read other point of Rotation
			m_aSecondPointOfRotation.push_back(tempOtherPointOfRotation);

			Scalar tempWallFrictionCoefficient = 0;
			LineStream >> tempWallFrictionCoefficient; //Read Wall Friction Coefficient
			m_aWallFrictionCoefficient.push_back(tempWallFrictionCoefficient);

			Scalar tempWallThermalResistivity = 0;
			LineStream >> tempWallThermalResistivity; //Read Wall Thermal Restivity
			m_aWallThermalRestivity. push_back(tempWallThermalResistivity);

			Scalar tempWallThermalConductivity = 0;
			LineStream >> tempWallThermalConductivity; //Read Wall Thermal Conductivity
			m_aWallThermalConductivity.push_back(tempWallThermalConductivity);

			m_aSTLISDrag.push_back(false);
			m_aSTLDragID.push_back(-1);
		}
		getline(ConfigureFile,Line);
	}
	
	m_aNumberOfDragSTLs = CountDragID;

	getline(ConfigureFile,m_OuputPath);

	getline(ConfigureFile,m_OutputName);

	getline(ConfigureFile,Line);
	
	LineStream.clear();
	LineStream.str(Line);

	LineStream >>m_OutputInterval;
	//cout<<m_OutputInterval;
	if(ConfigureFile.is_open())
	{
		getline(ConfigureFile, Line);
		LineStream.clear();
		LineStream.str(Line);
		CBox* ComputeZone = new CBox;

		LineStream>>ComputeZone->m_MinBound.x;
		LineStream>>ComputeZone->m_MinBound.y;
		LineStream>>ComputeZone->m_MinBound.z;
		LineStream>>ComputeZone->m_MaxBound.x;
		LineStream>>ComputeZone->m_MaxBound.y;
		LineStream>>ComputeZone->m_MaxBound.z;
		// We should test wether file has the valid box or not
		bool pass = true;
		if( fabs(ComputeZone->m_MaxBound.x - ComputeZone->m_MinBound.x) == 0.0)
		{
			pass = false;
		}
		if( fabs(ComputeZone->m_MaxBound.y - ComputeZone->m_MinBound.y) == 0.0)
		{
			pass = false;
		}
		if( fabs(ComputeZone->m_MaxBound.z - ComputeZone->m_MinBound.z) == 0.0)
		{
			pass = false;
		}
		if(pass)
		{
			m_pComputeZone = ComputeZone;
		}
		else
		{
			m_pComputeZone = NULL;
		}

		getline(ConfigureFile, Line);
		LineStream.clear();
		LineStream.str(Line);
		LineStream>>m_Section.x;
		LineStream>>m_Section.y;
		LineStream>>m_Section.z;

	}
	if(ConfigureFile.is_open())
	{
		getline(ConfigureFile,Line);
		LineStream.clear();
		LineStream.str(Line);
		LineStream >> m_StartStep;
	}
	//Read whether the file writing is binary or ascii
	if(ConfigureFile.is_open())
	{
		getline(ConfigureFile,Line);
		LineStream.clear();
		LineStream.str(Line);
		std::string isAscii;
		LineStream >> isAscii;
		boost::to_upper(isAscii);
		if(isAscii.compare("BINARY") == 0)
		{
			this->m_IsOutputAscii = false;
		}
		else
			this->m_IsOutputAscii = true;
	}
	ConfigureFile.close();

	return CCT_NOERR;
}
CCTStatusType CConfiguration::ReadConfigurationFile(const std::string& FileName)
{
	ifstream ConfigureFile;
	/*ConfigureFile.open ("C:\\Users\\ar1020\\Desktop\\ETCCTPM.Config");*/
	ConfigureFile.open(FileName.c_str());

	if(!ConfigureFile.is_open())
	{
		cout<<"Error in opening Configuration File";
		return CCT_FILEERR;
	}
	DestroyAll();

	string Line;
	stringstream LineStream;
	
	getline(ConfigureFile,m_InputPath);
	//cout<<m_InputPath<<endl;

	getline(ConfigureFile,m_ParameterFileName);
	//cout<<m_ParameterFileName<<endl;

	getline(ConfigureFile,m_ParticleFileName);
	//cout<<m_ParticleFileName<<endl;

	getline(ConfigureFile,Line);

	//LineStream.str(Line);

	while(!Line.empty())
	{
		string tempFileName;
		Integer additionInterval;
		Integer totalAdditionSteps;
		LineStream.clear(); //To clear the stream before read next line
		LineStream.str(Line);		
		LineStream >> tempFileName >> additionInterval >> totalAdditionSteps;
		m_aBucketFileName.push_back(tempFileName);
		m_aAdditionInterval.push_back(additionInterval);
		m_aTotalAdditionStep.push_back(totalAdditionSteps);
		getline(ConfigureFile,Line);
	}
	//Earlier Bucket File
	//LineStream>> m_BucketFileName >> m_AdditionIntervalTime >> m_TotalAddition;

	getline(ConfigureFile , Line);

	while(!Line.empty())
	{
		string tempFileName;

		LineStream.clear(); //To clear the stream before read next line
		LineStream.str(Line);

		LineStream >> tempFileName;			//Read STL file Name
		m_aSTLFileName.push_back(tempFileName);
		
		Scalar3 temp = make_Scalar3(0,0,0);		
		LineStream>>temp.x>>temp.y>>temp.z; //Read STL Velocity
		m_aSTLVelocity.push_back(temp);	

		Scalar Temperature;
		LineStream>>Temperature;			//Read STL Temperatue
		m_aSTLTemperature.push_back(Temperature);

		Scalar3 tempPos = make_Scalar3(0,0,0);
		LineStream >> tempPos.x >> tempPos.y >> tempPos.z; //Read STL Offset
		m_aSTLPosition.push_back(tempPos);

		bool linearMovableOrNot = false;
		LineStream >> linearMovableOrNot ; //Read Either Linear Movable or Not
		m_aSTLRestToOriginal.push_back(linearMovableOrNot);

		Integer ResetInterval = 0;
		LineStream >> ResetInterval; //Read Reset interval for Movable wall
		m_aSTLRestInterval.push_back(ResetInterval);

		bool RotatingOrNot = false;
		LineStream >> RotatingOrNot ; //Read Either Rotating or Not Rotating
		m_aRotatingOrNot.push_back(RotatingOrNot);

		Integer rotationInRPM = 0;
		LineStream >> rotationInRPM;		//Read rotation  unit in RPM
		m_aRPM.push_back(rotationInRPM);

		/*Scalar angleInDegree = 0;
		LineStream >> angleInDegree; 		//Read Rotation unit in Degree
		m_aAngleOfRotationDegree.push_back(angleInDegree); */

		Scalar startTimeOfRotation = 0;
		LineStream >> startTimeOfRotation; //Read Start time of Rotation
		m_aStartTimeOfRotation.push_back(startTimeOfRotation);

		Scalar timeToReachRPM = 0;
		LineStream >> timeToReachRPM;		//Read Final time to Reach RPM
		m_aTimeToReachRPM.push_back(timeToReachRPM);

		Scalar3 tempCenterOfRotation = make_Scalar3(0,0,0);
		LineStream >> tempCenterOfRotation.x >> tempCenterOfRotation.y >> tempCenterOfRotation.z; //Read Center of Rotation
		m_aCenterOfRotation.push_back(tempCenterOfRotation);

		Scalar3 tempOtherPointOfRotation = make_Scalar3(0,0,0);
		LineStream >> tempOtherPointOfRotation.x >> tempOtherPointOfRotation.y >> tempOtherPointOfRotation.z;  //Read other point of Rotation
		m_aSecondPointOfRotation.push_back(tempOtherPointOfRotation);

		getline(ConfigureFile,Line);
	}

	getline(ConfigureFile,m_OuputPath);

	getline(ConfigureFile,m_OutputName);

	getline(ConfigureFile,Line);
	
	LineStream.clear();
	LineStream.str(Line);

	LineStream >>m_OutputInterval;
	//cout<<m_OutputInterval;
	if(ConfigureFile.is_open())
	{
		getline(ConfigureFile, Line);
		LineStream.clear();
		LineStream.str(Line);
		CBox* ComputeZone = new CBox;

		LineStream>>ComputeZone->m_MinBound.x;
		LineStream>>ComputeZone->m_MinBound.y;
		LineStream>>ComputeZone->m_MinBound.z;
		LineStream>>ComputeZone->m_MaxBound.x;
		LineStream>>ComputeZone->m_MaxBound.y;
		LineStream>>ComputeZone->m_MaxBound.z;
		// We should test wether file has the valid box or not
		bool pass = true;
		if( fabs(ComputeZone->m_MaxBound.x - ComputeZone->m_MinBound.x) == 0.0)
		{
			pass = false;
		}
		if( fabs(ComputeZone->m_MaxBound.y - ComputeZone->m_MinBound.y) == 0.0)
		{
			pass = false;
		}
		if( fabs(ComputeZone->m_MaxBound.z - ComputeZone->m_MinBound.z) == 0.0)
		{
			pass = false;
		}
		if(pass)
		{
			m_pComputeZone = ComputeZone;
		}
		else
		{
			m_pComputeZone = NULL;
		}

		getline(ConfigureFile, Line);
		LineStream.clear();
		LineStream.str(Line);
		LineStream>>m_Section.x;
		LineStream>>m_Section.y;
		LineStream>>m_Section.z;

	}
	if(ConfigureFile.is_open())
	{
		getline(ConfigureFile,Line);
		LineStream.clear();
		LineStream.str(Line);
		LineStream >> m_StartStep;
	}
	//Read whether the file writing is binary or ascii
	if(ConfigureFile.is_open())
	{
		getline(ConfigureFile,Line);
		LineStream.clear();
		LineStream.str(Line);
		std::string isAscii;
		LineStream >> isAscii;
		boost::to_upper(isAscii);
		if(isAscii.compare("BINARY") == 0)
		{
			this->m_IsOutputAscii = false;
		}
		else
			this->m_IsOutputAscii = true;
	}
	ConfigureFile.close();

	return CCT_NOERR;
}

const string CConfiguration::GetParameterFile()
{
	if(m_InputPath.empty())
	{
		return m_ParameterFileName;
	}

	string ParameterFile;
	ParameterFile.append(m_InputPath); //To append the Input Path.
	
	if(ParameterFile.find_last_of("/\\") != (ParameterFile.size()-1))  //To append the slash if required.
	{
		ParameterFile.append("\\");
	}
	ParameterFile.append(m_ParameterFileName); //To append the Parameter File Name

	return ParameterFile;	

}
const string CConfiguration::GetParticleFile()
{
	if(m_InputPath.empty())
	{
		return m_ParticleFileName;
	}

	string ParticleFile;
	ParticleFile.append(m_InputPath);
	
	if(ParticleFile.find_last_of("/\\") != (ParticleFile.size()-1))
	{
		ParticleFile.append("\\");
	}
	ParticleFile.append(m_ParticleFileName);

	return ParticleFile;
}

const string CConfiguration::GetBucketFile()
{
	if(m_InputPath.empty())
	{
		return m_BucketFileName;
	}
	if(m_BucketFileName.empty())
	{
		return m_BucketFileName;
	}

	string BucketFile;
	BucketFile.append(m_InputPath);
	
	if(BucketFile.find_last_of("/\\") != (BucketFile.size()-1))
	{
		BucketFile.append("\\");
	}
	BucketFile.append(m_BucketFileName);

	return BucketFile;	

}

const vector<string> CConfiguration::GetModelFiles() const
{
	vector<string> STLFileName;
	
	for(size_t i=0;i<m_aSTLFileName.size();i++)
	{
		if(m_InputPath.empty())
		{
			return m_aSTLFileName;
		}

		string FullPath;
		FullPath.append(m_InputPath);
		if(FullPath.find_last_of("/\\") != (FullPath.size()-1))
		{
			FullPath.append("\\");
		}
		FullPath.append(m_aSTLFileName[i]);

		STLFileName.push_back(FullPath); //push the apended STLFileName into vector.
	}

	return STLFileName;	
}
//For Multi Bucket Addtion Starts
const vector<string> CConfiguration::GetMultiBucketFiles() const
{
	vector<string> BucketFileName;
	
	for(size_t i=0;i<m_aBucketFileName.size();i++)
	{
		if(m_InputPath.empty())
		{
			return m_aBucketFileName;
		}

		string FullPath;
		FullPath.append(m_InputPath);
		if(FullPath.find_last_of("/\\") != (FullPath.size()-1))
		{
			FullPath.append("\\");
		}
		FullPath.append(m_aBucketFileName[i]);

		BucketFileName.push_back(FullPath); //push the apended STLFileName into vector.
	}

	return BucketFileName;	
}
const vector<Integer>& CConfiguration::GetMultiAdditionInterval() const
{
	return m_aAdditionInterval;
}
const vector<Integer>& CConfiguration::GetMultiTotalAdditionStep() const
{
	return m_aTotalAdditionStep;
}
//For multi Bucket addition Ends

//For Linear Motion
Scalar CConfiguration::GetAdditionIntervalTime()
{
	return m_AdditionIntervalTime;
}
Scalar CConfiguration::GetTotalAddition()
{
	return m_TotalAddition;
}
const vector<Scalar3>& CConfiguration::GetModelVelocity() const
{
	return m_aSTLVelocity;
}
const vector<Scalar3>& CConfiguration::GetModelPosition() const
{
	return m_aSTLPosition;
}
const vector<bool>& CConfiguration::GetModelResetToOriginal() const
{
	return m_aSTLRestToOriginal;
}
const vector<Integer>& CConfiguration::GetModelResetInterval() const
{
	return m_aSTLRestInterval;
}
const string CConfiguration::GetOutputFileName()
{
	if(m_OuputPath.empty())
	{
		return m_OutputName;
	}

	string OutputName;
	OutputName.append(m_OuputPath);
	
	if(OutputName.find_last_of("/\\") != (OutputName.size()-1))
	{
		OutputName.append("\\");
	}
	OutputName.append(m_OutputName);

	return OutputName;	
}
string CConfiguration::GetOutputPath() const
{
	if(m_OuputPath.empty())
	{
		return m_OuputPath;
	}
	string OutputPath;
	OutputPath.append(m_OuputPath);
	if(OutputPath.find_last_of("/\\") != (OutputPath.size()-1))
	{
		OutputPath.append("\\");
	}
	return OutputPath;
}
Integer CConfiguration::GetOutputInterval()
{
	return m_OutputInterval;
}

Scalar3 CConfiguration::GetModelVelocity(unsigned int Model) const
{
	
	if(Model < m_aSTLVelocity.size() && Model >= 0)
	{
		return m_aSTLVelocity[Model];
	}
	else
	{
		Scalar3 temp = {0,0,0};
		return temp;
	}
}
Scalar CConfiguration::GetModelTemperature(unsigned int Model) const
{
	
	if(Model < m_aSTLTemperature.size() && Model >= 0)
	{
		return m_aSTLTemperature[Model];
	}
	else
	{
		return 0.0;
	}
}
Scalar3 CConfiguration::GetModelPosition(unsigned int Model) const
{
	
	if(Model < m_aSTLPosition.size() && Model >= 0)
	{
		return m_aSTLPosition[Model];
	}
	else
	{
		Scalar3 temp = {0,0,0};
		return temp;
	}
}
bool CConfiguration::GetModelResetToOriginal(unsigned int Model) const
{	
	if(Model < m_aSTLRestToOriginal.size() && Model >= 0)
	{
		return m_aSTLRestToOriginal[Model];
	}
	else
	{
		bool temp = false;
		return temp;
	}
}
Integer CConfiguration::GetModelResetInterval(unsigned int Model) const
{	
	if(Model < m_aSTLRestInterval.size() && Model >= 0)
	{
		return m_aSTLRestInterval[Model];
	}
	else
	{
		Integer temp = 0;
		return temp;
	}
}
CBox* CConfiguration::GetComputeZone() const
{
	return m_pComputeZone;
}
void CConfiguration::DestroyAll()
{
	if(m_pComputeZone)
	{
		delete m_pComputeZone;
	}
}
Integer3 CConfiguration::GetSection() const
{
	return m_Section;
}
bool CConfiguration::IsAscii() const
{
	return m_IsOutputAscii;
}
Integer CConfiguration::GetStartStep() const
{
	return m_StartStep;
}
//For Rotating STL Starts
///For Rotating Triangle Model
bool CConfiguration::GetTrianlgeRotatingOrNot(unsigned int Model) const
{	
	if(Model < m_aRotatingOrNot.size() && Model >= 0)
	{
		return m_aRotatingOrNot[Model];
	}
	else
	{
		bool temp = false;
		return temp;
	}
}
Scalar CConfiguration::GetAngleDegree(unsigned int Model) const
{	
	if(Model < m_aAngleOfRotationDegree.size() && Model >= 0)
	{
		return m_aAngleOfRotationDegree[Model];
	}
	else
	{
		Scalar temp = 0;
		return temp;
	}
}
Integer CConfiguration::GetRotationRPM(unsigned int Model) const
{	
	if(Model < m_aRPM.size() && Model >= 0)
	{
		return m_aRPM[Model];
	}
	else
	{
		Integer temp = 0;
		return temp;
	}
}
Scalar CConfiguration::GetRotationStartTime(unsigned int Model) const
{	
	if(Model < m_aStartTimeOfRotation.size() && Model >= 0)
	{
		return m_aStartTimeOfRotation[Model];
	}
	else
	{
		Scalar temp = 0;
		return temp;
	}
}
Scalar CConfiguration::GetTimeToReachRPM(unsigned int Model) const
{	
	if(Model < m_aTimeToReachRPM.size() && Model >= 0)
	{
		return m_aTimeToReachRPM[Model];
	}
	else
	{
		Scalar temp = 0;
		return temp;
	}
}
Scalar3 CConfiguration::GetCenterOfRotation(unsigned int Model) const
{	
	if(Model < m_aCenterOfRotation.size() && Model >= 0)
	{
		return m_aCenterOfRotation[Model];
	}
	else
	{
		Scalar3 temp = {0,0,0};
		return temp;
	}
}
Scalar3 CConfiguration::GetOtherPointOfRotation(unsigned int Model) const
{	
	if(Model < m_aSecondPointOfRotation.size() && Model >= 0)
	{
		return m_aSecondPointOfRotation[Model];
	}
	else
	{
		Scalar3 temp = {0,0,0};
		return temp;
	}
}
const vector<Scalar>& CConfiguration::GetRotationAngle() const
{
	return m_aAngleOfRotationDegree;
}
const vector<Integer>& CConfiguration::GetRotationRPM() const
{
	return m_aRPM;
}
const vector<Scalar>& CConfiguration::GetRotationStartTime() const
{
	return m_aStartTimeOfRotation;
}
const vector<Scalar>& CConfiguration::GetTimeToReachRPM() const
{
	return m_aTimeToReachRPM;
}

const vector<Scalar3>& CConfiguration::GetCenterOfRotation() const
{
	return m_aCenterOfRotation;
}
const vector<Scalar3>& CConfiguration::GetOtherPointOfRotation() const
{
	return m_aSecondPointOfRotation;
}

void CConfiguration::SetRotationInAngle(const Scalar RotationDegree)
{
	m_aAngleOfRotationDegree.push_back(RotationDegree);
}
//For Rotating STL Ends

//For Drag Starts
bool CConfiguration::GetModelISDrag(unsigned int Model) const
{
	if(Model < m_aSTLISDrag.size() && Model >= 0)
	{
		return m_aSTLISDrag[Model];
	}
	else
	{
		bool temp = false;
		return temp;
	}
}
Integer CConfiguration::GetModelDragID(unsigned int Model) const
{
	if(Model < m_aSTLDragID.size() && Model >= 0)
	{
		return m_aSTLDragID[Model];
	}
	else
	{
		Integer temp = 0;
		return temp;
	}
}
DragParameter CConfiguration::GetSTLDragParameter(unsigned int Model) const
{
	if(Model < m_aSTLDragParameter.size() && Model >= 0)
	{
		return m_aSTLDragParameter[Model];
	}
	else
	{
		DragParameter temp;
		return temp;
	}
}
Integer	CConfiguration::GetNumberOfDragSTL() const
{
	return m_aNumberOfDragSTLs;
}
//For Drag Ends

Scalar CConfiguration::GetWallFrictionCoefficient(unsigned int Model) const
{	
	if(Model < m_aWallFrictionCoefficient.size() && Model >= 0)
	{
		return m_aWallFrictionCoefficient[Model];
	}
	else
	{
		Scalar temp = 0;
		return temp;
	}
}
const vector<Scalar>& CConfiguration::GetWallFrictionCoefficient() const
{
	return m_aWallFrictionCoefficient;
}
Scalar CConfiguration::GetWallThermalConductivity(unsigned int Model) const
{	
	if(Model < m_aWallThermalConductivity.size() && Model >= 0)
	{
		return m_aWallThermalConductivity[Model];
	}
	else
	{
		Scalar temp = 0;
		return temp;
	}
}
const vector<Scalar>& CConfiguration::GetWallThermalConductivity() const
{
	return m_aWallThermalConductivity;
}
Scalar CConfiguration::GetWallThermalRestivity(unsigned int Model) const
{	
	if(Model < m_aWallThermalRestivity.size() && Model >= 0)
	{
		return m_aWallThermalRestivity[Model];
	}
	else
	{
		Scalar temp = 0;
		return temp;
	}
}
const vector<Scalar>& CConfiguration::GetWallThermalRestivity() const
{
	return m_aWallThermalRestivity;
}