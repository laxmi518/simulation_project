#pragma once

class CInnerPressureModel
{
	//parameters

	double TankVolume;
	double TankLimit;
	double Gamma;
	//double GammaSplash;
	double BreatherRadius;
	double LBreather;
	double JointRadius;
	double PAtmospheric;
	double JointZeta;
	double PressureCoefficient;
	double Vv_coef;
	double kinematicViscosity;
	double V_gasInitial;
	double alpha;
	double P_base;
	double Rho;
	double Q_FlowRate_LiterPerMinute;	

	//For Inner Pressure Calculation
	double Qg_out;
	double QI_in;
	double Re_b;
	double Re_j;
	double V_tank_current;
	double V_gas;
	double V_j;
	double V_b;	
	double Pressure_current;
	double Pressure_old;
	double H_joint;
	double H_breather;

	double P1;
	double P2;
	double P3;

	double P1_P2;
	double P2_P3;

	bool check;
	std::string m_OutputPath;
	CParameter* m_pParameter;
public:
	Scalar CorrelatedDensity;
	CInnerPressureModel(CParameter* Parameter);
public:
	virtual ~CInnerPressureModel(void);

private:
	//For Inner Pressure Calculation
	CCTStatusType FlowCalculation(Integer particle_num_out);
	CCTStatusType CalculationParameterPreparation();
	CCTStatusType TankInernalPressureCalc();
	CCTStatusType BreatherLossCalc(Integer AnalysisStep, Integer particle_num_out);
	CCTStatusType PressureCalc(Integer AnalysisStep, Integer particle_num_out);
	
public:
	CCTStatusType Load(std::string& FileName);
	void SetOutputPath(std::string FileName);
	CCTStatusType CalcInnerPressures(Integer AnalysisStep, Integer particle_num_out,double& P1_P2, double& P2_P3);
	CCTStatusType HeightCalc(Integer AnalysisStep, Scalar Height);
};
