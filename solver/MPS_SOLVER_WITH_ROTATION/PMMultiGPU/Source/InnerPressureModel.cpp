#include "stdafx.h"
#include "InnerPressureModel.h"
#include "IO.h"

CInnerPressureModel::CInnerPressureModel(CParameter* Parameter)
:m_pParameter(Parameter)
{
	 TankVolume = 0;
	 TankLimit= 0;
	 Gamma= 0;
	 //GammaSplash = 0;
	 BreatherRadius= 0;
	 LBreather= 0;
	 JointRadius= 0;
	 PAtmospheric= 0;
	 JointZeta= 0;
	 PressureCoefficient= 0;
	 Vv_coef= 0;
	 kinematicViscosity= 0;
	 V_gasInitial= 0;
	 alpha=0;
	 P_base = 0;
	 Rho = 0;
	 Q_FlowRate_LiterPerMinute = 0;
}

CInnerPressureModel::~CInnerPressureModel(void)
{
}
//For Inner Pressure Calculation
CCTStatusType CInnerPressureModel::FlowCalculation(Integer particle_num_out)
{	
	//Gasoline into the Tank
	QI_in = particle_num_out * m_pParameter->InitialParticleDistance * m_pParameter->InitialParticleDistance * m_pParameter->InitialParticleDistance;

	//Amout of gas flowing from the tank (including vapour)
	//Qg_out = particle_num_out * m_pParameter->InitialParticleDistance * m_pParameter->InitialParticleDistance * m_pParameter->InitialParticleDistance +  Vv_coef * m_pParameter->Dt;
	P_base = alpha * ((Q_FlowRate_LiterPerMinute*1000)/(BreatherRadius*BreatherRadius));
	Qg_out = (P3 - P1) / P_base * (QI_in + Vv_coef * m_pParameter->Dt);

	if(V_tank_current < TankLimit)
	{
		Qg_out = 0;
	}

	return CCT_NOERR;	
}

//Calculation of Parameter Preparation
CCTStatusType CInnerPressureModel:: CalculationParameterPreparation()
{	
	double V_tank_pre = 0;
	
	//Breather Velocity Calculation
	double denominator = (CONSTANT_CIRCULAR * BreatherRadius * BreatherRadius);
	if(denominator != 0 )
	{
		V_b = Qg_out / denominator;
	}
	else
	{
		V_b = 0.0f;
	}

	//Breather Re Number Calculation
	if(kinematicViscosity != 0 )
	{
		Re_b = V_b * (2 * BreatherRadius) / kinematicViscosity;
	}
	else
	{
		Re_b = 0.0f;
	}

	//Breather Joint Velocity Calculation
	double denominator1 = (CONSTANT_CIRCULAR * JointRadius * JointRadius);
	if(denominator1 != 0)
	{
		V_j = Qg_out / denominator1;
	}else
	{
		V_j = 0.0;
	}

	//Calculate Number of Breather Joints(Re)
	if(kinematicViscosity != 0)
	{
		Re_j = V_j * (2 * BreatherRadius) / kinematicViscosity;
	}else
	{
		Re_j = 0.0;
	}

	//Tank Capacity Calculation
	//Tank Capacity of Previous Step
	V_tank_pre = V_tank_current;

	//Current Tank Capacity
	V_tank_current = V_tank_pre - QI_in;

	//Gas Volume Calculations
	//V_gas = V_tank_pre + Vv_coef * m_pParameter->Dt;
	V_gas = V_gas +  Vv_coef *  m_pParameter->Dt - Qg_out;

	if(V_gas < 0)
	{
		V_gas = 0;
	}

	return CCT_NOERR;
}
CCTStatusType CInnerPressureModel::TankInernalPressureCalc()
{
	//Tank Inner Pressure Calculation
	if(V_tank_current > 0 )
	{
		double ratio = V_gas / V_tank_current;
		//double CurrentGamma = (V_tank_current < TankLimit) ? GammaSplash : Gamma;
		Pressure_current = PAtmospheric * pow(ratio , Gamma);
	}else
	{
		Pressure_current = 0;
	}
	
	return CCT_NOERR;
}
CCTStatusType CInnerPressureModel::BreatherLossCalc(Integer AnalysisStep, Integer particle_num_out)
{
	//Breather joint Loss Calculations
	H_joint = JointZeta * Rho * pow(V_j,2) / 2 * static_cast<double>(9.80665);

	//Breather Loss Calculation
	//Friction Coefficient Calculation
	double lambda_f = 0.0;
	if(Re_b != 0)
	{
		lambda_f =  static_cast<double>(0.3164) / pow(Re_b , static_cast<double>(0.25));
	}

	if(BreatherRadius != 0)
	{
		H_breather = lambda_f * LBreather * Rho * (pow(V_b , 2)) / (2 * BreatherRadius) * 9.80655;
	}else
	{
		H_breather = 0;
	}
	//if( 0 == (AnalysisStep % m_pParameter->OutputStep))
	//{
		//IO::SavePressures(m_OutputPath, Pressure_current, particle_num_out, QI_in, V_gas, V_tank_current);
		//IO::SavePressures(m_OutputPath, Vv_coef*m_pParameter->Dt, Pressure_old, V_gas/V_tank_current, 0, 0);
		//IO::SavePressures(m_OutputPath, JointZeta, Rho, V_j, 0, 0);
		//IO::SavePressures(m_OutputPath, Qg_out, CONSTANT_CIRCULAR, JointRadius, 0, 0);
		//IO::SavePressures(m_OutputPath, P3, P1, P_base, particle_num_out, m_pParameter->InitialParticleDistance);
		//IO::SavePressures(m_OutputPath, P3, Vv_coef, m_pParameter->Dt, Qg_out_min, 0);		
	//}
	return CCT_NOERR;
}
CCTStatusType CInnerPressureModel::PressureCalc(Integer AnalysisStep, Integer particle_num_out)
{
	P2 = Pressure_current - H_joint - H_breather;
	P3 = Pressure_current;

	P1_P2 = (P2 - P1) * PressureCoefficient;
	P2_P3 = (P3 - P1) * PressureCoefficient;

	if( 0 == (AnalysisStep % 20))
	{
		IO::SavePressures(m_OutputPath + "InnerPressureResult.txt", AnalysisStep*m_pParameter->Dt, P2, P3-PAtmospheric, H_joint, H_breather);
	}
	return CCT_NOERR;
}
CCTStatusType CInnerPressureModel::HeightCalc(Integer AnalysisStep, Scalar Height)
{	
	if( 0 == (AnalysisStep % 20))
	{
		IO::SaveHeight(m_OutputPath + "CorrelatedHeight.dat", AnalysisStep*m_pParameter->Dt, Height);
	}
	return CCT_NOERR;
}

CCTStatusType  CInnerPressureModel::Load(std::string& FileName)
{
	//InnerPressureParameterFileName = RootPath + InnerPressureParameter;
	std::ifstream InnerPressureParameterFile(FileName.c_str());
	if(InnerPressureParameterFile.is_open())
	{
		InnerPressureParameterFile >> TankVolume;
		InnerPressureParameterFile >> TankLimit;
		InnerPressureParameterFile >> Gamma;
		//InnerPressureParameterFile >> GammaSplash;
		InnerPressureParameterFile >> BreatherRadius;
		InnerPressureParameterFile >> LBreather;
		InnerPressureParameterFile >> JointRadius;
		InnerPressureParameterFile >> PAtmospheric;
		InnerPressureParameterFile >> JointZeta;
		InnerPressureParameterFile >> PressureCoefficient;
		InnerPressureParameterFile >> Vv_coef;
		InnerPressureParameterFile >> kinematicViscosity;
		InnerPressureParameterFile >> V_gasInitial;
		InnerPressureParameterFile >> alpha;
		InnerPressureParameterFile >> Rho;
		InnerPressureParameterFile >> Q_FlowRate_LiterPerMinute;
		InnerPressureParameterFile >> CorrelatedDensity;
	}
	if(InnerPressureParameterFile.fail())
	{
		std::cout<< "Innerpressure parameter load error";
		return CCT_FILEERR;
	}
	InnerPressureParameterFile.close();

	//Initial Values
	Pressure_old = PAtmospheric;
	V_tank_current = TankVolume;
	V_gas = V_gasInitial;

	P1 = PAtmospheric;
	P2 = PAtmospheric;
	P3 = PAtmospheric;
		
	return CCT_NOERR;
}
void CInnerPressureModel::SetOutputPath(std::string FileName)
{
	m_OutputPath = FileName;
	std::ofstream OutFile(FileName.c_str());
	OutFile.close();

}
CCTStatusType CInnerPressureModel::CalcInnerPressures(Integer AnalysisStep ,Integer particle_num_out,double& P1_P2, double& P2_P3)
{
	//For Inner Pressure Calculation
	CCTStatusType Status = CCT_NOERR;
	Status = FlowCalculation(particle_num_out);
	CCT_ERROR_CHECK(Status);

	Status = CalculationParameterPreparation();
	CCT_ERROR_CHECK(Status);

	Status = TankInernalPressureCalc();
	CCT_ERROR_CHECK(Status);

	Status = BreatherLossCalc(AnalysisStep, particle_num_out);
	CCT_ERROR_CHECK(Status);

	Status = PressureCalc(AnalysisStep, particle_num_out);
	CCT_ERROR_CHECK(Status);	

	P1_P2 = this->P1_P2;
	P2_P3 = this->P2_P3;
	return CCT_NOERR;
}