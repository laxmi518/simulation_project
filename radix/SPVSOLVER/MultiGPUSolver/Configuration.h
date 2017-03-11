#pragma once

//#include "Scalar3.h"
//#include "std.h"

using namespace std;

class CConfiguration
{
public:
	CConfiguration(void);
	~CConfiguration(void);

	
private:
	//Input varaible
	string m_InputPath;
	string m_ParameterFileName;
	string m_ParticleFileName;
	string m_BucketFileName;
	Scalar m_AdditionIntervalTime;
	Scalar m_TotalAddition;

	//For Linear Motion Starts
	vector<string>	m_aSTLFileName;			//STL File Name
	vector<Scalar3> m_aSTLVelocity;			//STL Velocity
	vector<Scalar>	m_aSTLTemperature;		//STL Temperatuer
	vector<Scalar3> m_aSTLPosition;			//STL Offset
	vector<bool>	m_aSTLRestToOriginal;	//Linear Movable or Not.
	vector<Integer> m_aSTLRestInterval;		//Linear Reset Interval
	//For Liner Motion Ends

	//For Multiple Bucket Addtion Starts
	vector<string>	m_aBucketFileName;		//Bucket File Name
	vector<Integer> m_aAdditionInterval;	//Bucket Addition Interval
	vector<Integer> m_aTotalAdditionStep;	//Total Bucket addition Interval
	//For Multiple Bucket Addition Ends

	///For Rotating Triangle
	vector<bool>	m_aRotatingOrNot;			//Rotating Triangle or Not
	vector<Scalar>  m_aAngleOfRotationDegree;	//Angle of Rotating in Degree
	vector<Integer> m_aRPM;						//Angle of Rotation in RPM
	vector<Scalar>	m_aStartTimeOfRotation;		//Start Time of Rotation
	vector<Scalar>	m_aTimeToReachRPM;			//Time or reach Required RPM
	vector<Scalar3> m_aCenterOfRotation;		//Center of Rotation
	vector<Scalar3> m_aSecondPointOfRotation;	//Second point of Rotation	
	vector<Scalar>  m_aCoefficientofInflueneRegion; //Influene Region
	vector<Scalar3> m_aDirectionVector;         // Direction vector
	vector<bool>	m_aNegativePressureCheck;	//Read Either Check for Negative Pressure or not
	vector<Scalar>  m_aPressureGas;				// Pressure gas for each STL

	////For Rotating Triangle Ends

	//For Temperature Analysis Starts
	vector<Scalar> m_aWallThermalRestivity;
	vector<Scalar> m_aWallThermalConductivity;
	vector<Scalar> m_aWallFrictionCoefficient;
	//For Temperature Analysis Ends


	//output variable
	string m_OuputPath;
	string m_OutputName;
	Integer m_OutputInterval;
	CBox* m_pComputeZone;
	bool m_IsOutputAscii;
	Integer3 m_Section;
	Integer m_StartStep;
	void DestroyAll();
public:
	//Function to read input file
	//CCTStatusType ReadConfigurationFile(const std::string& FileName);
	const string GetParameterFile();
	const string GetParticleFile();

	const string GetBucketFile();
	Scalar GetAdditionIntervalTime();
	Scalar GetTotalAddition();

	//Function to read STL file
	const vector<string> GetModelFiles() const;
	const vector<Scalar3>& GetModelVelocity() const ;
	const vector<Scalar3>& GetModelPosition() const;
	Scalar3 GetModelVelocity(unsigned int Model) const;
	Scalar GetModelTemperature(unsigned int Model) const;
	Scalar3 GetModelPosition(unsigned int Model) const;

	//Function to read output file
	const string GetOutputFileName();
	string GetOutputPath() const;
	Integer GetOutputInterval();
	CBox* GetComputeZone() const;
	Integer3 GetSection() const;
	bool IsAscii() const;
	Integer GetStartStep() const;

	//For multiple Bucket
	const vector<string> GetMultiBucketFiles() const;
	const vector<Integer>& GetMultiAdditionInterval() const;
	const vector<Integer>& GetMultiTotalAdditionStep() const;
	const vector<bool>& GetModelResetToOriginal() const;
	const vector<Integer>& GetModelResetInterval() const;
	bool GetModelResetToOriginal(unsigned int Model) const;
	Integer GetModelResetInterval(unsigned int Model) const;

	//For Rotating Triangle
	bool	GetTrianlgeRotatingOrNot(unsigned int Model) const;
	Scalar	GetAngleDegree(unsigned int Model) const;
	Integer GetRotationRPM(unsigned int Model) const;
	Scalar	GetRotationStartTime(unsigned int Model) const;
	Scalar	GetTimeToReachRPM(unsigned int Model) const;
	Scalar3 GetCenterOfRotation(unsigned int Model) const;
	Scalar3 GetOtherPointOfRotation(unsigned int Model) const;
	Scalar	GetInflueneRegion(unsigned int Model) const;
	Scalar3 GetDirectionVector(unsigned int Model) const;
	bool	GetNegativePressureCheck(unsigned int Model) const;
	Scalar	GetPressureGas(unsigned int Model) const;

	const vector<Scalar>& GetRotationAngle() const;
	const vector<Integer>& GetRotationRPM() const;
	const vector<Scalar>& GetRotationStartTime() const;
	const vector<Scalar>& GetTimeToReachRPM() const;
	const vector<Scalar3>& GetCenterOfRotation() const;
	const vector<Scalar3>& GetOtherPointOfRotation() const;
	const vector<Scalar>& GetInflueneRegion() const;
	const vector<Scalar3>& GetDirectionVector() const;
	const vector<bool>& GetNegativePressureCheck() const;

	void SetRotationInAngle(const Scalar RotationDegree);
	////For Rotating Triangle Ends

	//for drag 
	vector<bool>	m_aSTLISDrag;
	vector<Integer> m_aSTLDragID;
	vector<DragParameter> m_aSTLDragParameter;
	Integer			m_aNumberOfDragSTLs;
	//for drag ends

		//for Drag
	bool	GetModelISDrag(unsigned int Model) const;	//Added by Arpan
	Integer GetModelDragID(unsigned int Model) const;
	DragParameter GetSTLDragParameter(unsigned int Model) const;
	Integer GetNumberOfDragSTL() const;	//Added by Arpan 20130129

	//For Wall Parameter in Configuration File
	Scalar GetWallThermalRestivity(unsigned int Model) const;
	Scalar GetWallThermalConductivity(unsigned int Model) const;
	Scalar GetWallFrictionCoefficient(unsigned int Model) const;

	const vector<Scalar>& GetWallFrictionCoefficient() const;
	const vector<Scalar>& GetWallThermalConductivity() const;
	const vector<Scalar>& GetWallThermalRestivity() const;

	CCTStatusType ReadConfigurationFileWithDrag(const std::string& FileName);
};
