#pragma once
#include <cassert>
#include <iostream>
#include <cmath>
#include <float.h>

#define SQRT_MAGIC_F 0x5f3759df
#ifdef DOUBLE
typedef double	Scalar;
typedef double3 Scalar3;
#define make_Scalar3(x,y,z) make_double3(x,y,z)
#else if
typedef float	Scalar;
typedef float3 Scalar3;
#endif
typedef int		Integer;
typedef bool	Bool;
typedef int3 Integer3;

#define M_PI 3.141592653589793238
#define SCALAR_EPSILON FLT_EPSILON
#define ET_INLINE static inline
//Maximum Block Size
#define BLOCK_MAX_DIM 256
//#define BLOCK_MAX_DIM 1024
//static __inline__ __device__ __host__  float  sqrt_datatype(const float x)
//{
//  const float xhalf = 0.5f*x;
// 
//  union // get bits for floating value
//  {
//    float x;
//    int i;
//  } u;
//  u.x = x;
//  u.i = SQRT_MAGIC_F - (u.i >> 1);  // gives initial guess y0
//  return x*u.x*(1.5f - xhalf*u.x*u.x);// Newton step, repeating increases accuracy
//
//}

static __inline__ __device__ __host__ float sqrt_datatype ( float x){
  float xhalf = 0.5f*x;
  int i = *(int*)&x;
  i = 0x5f3759df - (i>>1);
  x = *(float*)&i;
  x = x*(1.5f - xhalf*x*x);
  return 1/x;}



static __inline__ __host__ __device__ Scalar SquaredMagnitude(const Scalar3& A)
{
	return A.x * A.x + A.y * A.y + A.z * A.z;
}

static __inline__ __host__ __device__ Scalar Magnitude(const Scalar3& A)
{
	return sqrt_datatype(SquaredMagnitude(A));
}
static __inline__ __host__ __device__ Scalar3 make_Scalar3(Scalar x, Scalar y, Scalar z)
{
  Scalar3 t; t.x = x; t.y = y; t.z = z; return t;
}
static __inline__ __host__ __device__ Integer3 make_Integer3(Integer x, Integer y, Integer z)
{
  Integer3 t; t.x = x; t.y = y; t.z = z; return t;
}


//typedef double	Scalar;


static const Scalar CONSTANT_CIRCULAR = static_cast<Scalar>(3.141592653589793);
#define EPS  1e-5


enum CCTStatusType
{
	CCT_NOERR,
	CCT_FILEERR,
	CCT_PARAMERR,
	CCT_CALCERR,
	CCT_MEMERR,
	CCT_CUDAERR,
	CCT_ETCERR
};

enum ParticleType
{
	/*NOTDEFINE = -1,
	FRUID = 0,
	PRESSURED_WALL = 1, 
	WALL = 2,
	PRESSURED_MOVABLE_WALL = 3,
	MOVABLE_WALL = 4*/

	/*FLUIDPARTICLE = 100,
	ELASTICPARTICLE = 101,
	WALLPARTICLE = 102,
	DELETEDPARTICLE = 200,*/

	TYPE_HEAVY_WEIGHT = 0,
	TYPE_HEAVY_WEIGHT_FREE_SURFACE = 1,
	TYPE_LIGHT_WEIGHT = 2
};

enum ComputeType
{
	COMPUTE_OWN = 0,
	COMPUTE_HEAVY_WEIGHT = 1,
	COMPUTE_LIGHT_WEIGHT = 2
};
#define HEAVY_FLUID FRUID ;
#define LIGHT_FLUID WALL ;

enum ExpImp
{
	EXPLICIT = 0,
	IMPLICIT = 1
};
struct DragParameter
{
	Integer DragID;
	char DragToken;
	char InputVelocityToken1;
	char OutputVelocityToken2;
	Scalar CoefficientofInflueneRegion;

	
	float VelocityIndex;
	float VelocityExponential;
	float VelocityDragPermiabilityConstant;
	//float* CoeffOfMember;

	float TempIndex;
	float TempExponential;
	float TempDragPermiabilityConstant;

	float ConstantVx;
	float ConstantVy;
	float ConstantVz;

	//Constant Parameter
	float VelocityMagnitueFactor;
};
//#define PMDEBUG
struct CTriangleParameters
{
	Scalar3 Velocity;
	Scalar Temperature;

	//For MultiBucket Addition Starts
	bool	LinearMovableOrNot;
	Integer LinearResetInterval;
	//For MultiBucket Addtion Ends

	//For STL rotation Starts
	bool	RotatingOrNot;
	Scalar	AngleOfRotationInDegree;
	Scalar  RotationInRPM;
	Scalar	RotationStartTime;
	Scalar	TimeToReachRPM;
	Scalar3	CenterOfRotation;
	Scalar3	SecondPointOfRotation;
	Scalar	InfluenceRegion;
	Scalar3 DirectionVector;
	bool	NegativePressureCheck;
	Scalar	PressureGas;
	//For STL rotation Ends


	Scalar WallFrictionCoefficient;
	Scalar WallThermalResistivity;
	Scalar WallThermalConductivity;

	bool isDrag;
	Integer DragTriangleID;
};

struct CTriangle
{
	Integer ModelIDNumber;
	Scalar3 Normal;
	Scalar3 Vertex[3];
/*
	Scalar3 Velocity;
	Scalar Temperature;

	//For MultiBucket Addition Starts
	bool	LinearMovableOrNot;
	Integer LinearResetInterval;
	//For MultiBucket Addtion Ends

	//For STL rotation Starts
	bool	RotatingOrNot;
	Scalar	AngleOfRotationInDegree;
	Scalar	RotationStartTime;
	Scalar	TimeToReachRPM;
	Scalar3	CenterOfRotation;
	Scalar3	SecondPointOfRotation;
	//For STL rotation Ends
*/
};

struct CDragTriangle
{
	Scalar3 Normal;
	Scalar3 Vertex[3];
	Scalar DragLength;
	Integer StartStep01;
	Integer StartStep02;
	Integer EndStep;
	Scalar3 InitialAcc;
	Scalar3 GradAcc;
	// ID
	Integer ID;

};

struct CParticle
{
	Integer ID;
	Scalar3 Position;
	Scalar3 Velocity;
	Scalar Pressure;
	//Scalar MinPressure;
	Scalar Density;
	Scalar Temperature;
	Scalar KineticViscosity;
	Scalar SolidPhaseRate;
	ParticleType Type;
	//Integer NeighborhoodNum;

	//Near Pressure for visco-elastic
	//Scalar DummyPressure;
	//Scalar NearPressure;
	////Neardensity for visco-elastic
	//Scalar DummyDensity;
	//Scalar NearDensity;
};
	/*static Integer *m_ParticleID;
	static Scalar3 *m_ParticlePosition;
	static Scalar3 *m_ParticleVelocity;
	static Scalar *m_ParticlePressure;
	static Scalar *m_ParticleDensity;
	static Scalar *m_ParticleTemperature;
	static Scalar *m_ParticleKineticViscosity;
	static Scalar *m_ParticleSolidPhaseRate;
	static ParticleType *m_ParticleType;*/

struct MultipleBucket
{
	Integer BucketNumber;
	CParticle*  Bucket;
	Integer AdditionInterval;
	Integer TotalAdditionStep;
	Integer AddedTimes;
};

struct CDistance
{
	Scalar Magnitude;
	Scalar3 Direction;
};
struct CBox
{
	Scalar3 m_MinBound;
	Scalar3 m_MaxBound;
};
struct CGridBox
{
	CBox m_BufferedZone;
	CBox m_ComputeZone;
	CBox m_PressureZone;
	CBox m_ParentZone;
	Integer3  m_GridSize;
	Scalar m_CellSize;
};
struct CGridParams
{
	//represents the first ID of the cell
	Integer m_BufferedCellIdStart;

	//represents the last ID of the cell
	Integer m_BufferedCellIdEnd;

	//represents the first cell buffer cell excluded
	Integer m_ComputeCellIdStart;

	//represents the last cell buffer cell excluded
	Integer m_ComputeCellIdEnd;

	//Total number of Cells in a x row
	Integer m_LeftGhostCellNum;
	Integer m_LeftGhostCellEnd;

	Integer m_RightGhostCellNum;
	Integer m_RightGhostCellStart;
};
struct CGhostIndex
{
	Integer m_LeftGhostStart;
	Integer m_LeftGhostEnd;
	Integer m_RightGhostStart;
	Integer m_RightGhostEnd;
	Integer m_BufferedParticleNum;
};
#define CELLMAXTRIANGLENUM 100
struct CCell
{
	//Integer m_ParticleID;// ID of the particle included in the cell, there can be multiple particles in a cell, so there is a Hash table in CCellSet class, to store more than one particle.
	//Integer m_HashID;// ID of the Hash to be filled, this is like a temporary variable, only used for filling the particles in the cell.
	unsigned int m_TriangleNum;
	Integer m_TriangleID[CELLMAXTRIANGLENUM];	
};

// —¬‘Ì–§“x‚Æ—±Žq”–§“x‚Í“¯‚¶‚È‚Ì‚Å—víœ
struct CParameter
{
	int Dimension;	
	// 重力定数
	Scalar3 GravityConstant;

	// 粒子数密度 n0
	Scalar Density;
	
	// 粒子数密度定数 ρ0
	Scalar ConstantParticleDensity;
	
	// 表面張力係数
	Scalar SurfaceTensionCoefficient;

	// 初期粒子間距離
	Scalar InitialParticleDistance;

	// 自由表面の判定に用いる閾値
	Scalar FreeSurfaceThreshold;

	// 勾配モデル影響半径係数	
	Scalar GradientInfluenceRadiusCoefficient;

	// ラプラシアンモデル影響半径係数
	Scalar LaplacianInfluenceRadiusCoefficient;
	
	// 温度計算用影響半径係数
	Scalar TemperatureInfluenceRadiusCoefficient;

	// 温度パラメータ

	// 流体初期温度
	Scalar FruidInitializeialTemperature;
	
	// 壁初期温度
	Scalar WallInitializeialTemperature;

	// 流体熱抵抗値 R_h1
	Scalar FruidThermalResistance;

	// 壁熱抵抗値 Rh_2
	Scalar WallThermalResistance;

	// 流体熱伝導率 k_fruid
	Scalar FruidThermalConductivity;

	// 壁熱伝導率 k_wall
	Scalar WallThermalConductivity;

	// 流体比熱 Cp
	Scalar FruidSpecificHeat;

	// 流体密度→Densityと同じ
	Scalar FruidDensity;

	// 計算制御用パラメータ

	// ステップ
	Scalar Dt;

	// 総ステップ数 
	Integer AnalysisStep;

	// 熱計算ステップ
	Integer ThermalStep;
	
	// 陰計算ステップ
	Integer NegativeCalculationStep;

	// 出力ステップ
	Integer OutputStep;

	// その他

	// λ値を得る粒子番号 
	Integer LambdaID;

	// SORヤコビアンの緩和係数
	Scalar SORJacobianRelaxationCoefficient;

	// SOR法の最大繰返数
	Integer SORLoopMax;
	
	// SOR法の収束定数
	Scalar SORConvergenceConstant;
	
	//温度計算を行うかどうか
	bool bUpdateTemperature;

	// λ値
	Scalar LambdaValue;
	// Number of cells in a bounding box
	Scalar3 CellSize;
	// Velocity of the moving Grid
	Scalar3 MvVelocity;

	Scalar CourantNumber;
	Scalar ViscosityUpdateParameterA;
	Scalar ViscosityUpdateParameterB;
	Scalar MaxWallTemperature;
	Scalar MinWallTemperature;
	
	// Offset to make ta boundary distance smooth
	Scalar OffSet;

	//Modification by E&T Nepal 
	Integer ParticleGenerationStep;
	Integer TotalGenerationStep;
	//Modification by E&T Nepal Ends

	Scalar Alpha;
	Scalar AlphaLightWeight;
	Integer MaxNeighborhoodNum;

	// ON/ OFF boolen variable for collision check
	bool bExplicitCollisonCheck;
	// ON/ OFF boolen variable to output collision checked data
	bool bExplicitCollisionOut;
	// Coefficient used in collision checking
	bool bImplicitCollisonCheck;
	Scalar CollisionCheckCoefficient;
	// Coefficient to check the min distance of the particle
	Scalar MinDistanceCoefficient;
	// ON/ OFF parameter for free surface pressure update
	bool bFreeSurfacePressureUpdate;
	// ON/ OFF parameter for Explicit Output data
	
	// Parameters for free surface pressure adjustment of fluids
	Scalar BasePressure;
	Scalar MaxPressure;
	Scalar Time1;
	Integer TriggerStep;
	Integer EndStep;
	Scalar AnalysisArea;


	bool bExplicitOut;

	// Friction Coefficient of the wall
	Scalar WallFrictionCoeffieicnt;

	// Drag Coefficient
	Scalar DragCoefficient;

	// For Two Phase Particle Method
	bool bTwoPhase;

	// Stop Wall Step
	Integer StopWallStep;

	Scalar Tolerance;
	//Implicit Surface Threshold For Two Phase
	Scalar ImplicitSurfaceThreshold;

	// For implementing CG method to solve pressure
	bool bCGSolver;

	// Inner pressure update boolen parameter
	bool bUpdateInnerPressure;

	// Surface Influence Radious Coefficient for Inner Pressure
	Scalar SurfaceInfluenceRadiusCoefficient;

	//Surface Threshold For inner pressure 
	Scalar SurfaceThreshold;
	
	//Gradient Force for inner pressure
	Scalar GradientForce;
	
	//Start step for inner pressure
	Integer StartStep;

	//Inner force for inner pressure calculation
	Scalar InitialForce;

	// Velocity Threshold for inner pressure calculation
	Scalar VelocityThreshold;

	// If want to update the pressure explicitly
	bool bExplicitPressure;
	Scalar SonicVelocity;
	Scalar ArtificialPressure;
	Scalar r_coeff;

	bool bImplicitVelocity;	

	bool bNJPCGSolver;

	Scalar PressureMax;
	Scalar RLimitCoefficient;

	bool bTurbulance;
	Scalar SmagorinskyConstant;
	Scalar FilterWidth;
	Scalar TurbulanViscosityRatio;
	Scalar Cnue;
	Scalar Cepsilon;
	Scalar TurbulanceIntensity;
	//New Parameter here!!!!

	//For Negative Pressure
	bool bNegativePressure;
	//Scalar PressureGas;

	//For Rotating STL model
	Scalar InitialTimeForZeroAngle;
	Scalar FinalTimeToReachAngle;

	Integer ConstantCheckStep;

	//Scalar DragPermiabilityConstant;
};
struct CParameterCoefficients
{
	Scalar ReLap;   //ReLap
	Scalar ReGrad;
	Scalar ReSurface;
	Scalar Rc;
	Scalar MinimumDistance;

	//For CalcExplicitly_Kernel
	//Scalar Coefficient = ( 2.0 * CONSTANT_PARAMETER.Dimension ) / ( CONSTANT_PARAMETER.LambdaValue * CONSTANT_PARAMETER.ConstantParticleDensity);
	Scalar CoefficientCalcExplicit;

	//For CalcTemperatureFactor_Kernel
	//const Scalar EnthalpyCoefficient = CONSTANT_PARAMETER.Dt * 2 * CONSTANT_PARAMETER.Dimension / ( CONSTANT_PARAMETER.ConstantParticleDensity * Density * CONSTANT_PARAMETER.FruidSpecificHeat);
	Scalar EnthalpyCoefficientCalcTempFactor;
	//ThermalConductivity = (2 * CONSTANT_PARAMETER.FruidThermalConductivity * CONSTANT_PARAMETER.WallThermalConductivity ) / ( CONSTANT_PARAMETER.FruidThermalConductivity + CONSTANT_PARAMETER.WallThermalConductivity);
	Scalar ThermalConductivityCalcTempFactor;

	//For Iterate_Kernel
	//Scalar CoefficientC = 2 * CONSTANT_PARAMETER.Dimension / ( CONSTANT_PARAMETER.LambdaValue * CONSTANT_PARAMETER.ConstantParticleDensity);
	Scalar CoefficientCIterate;

	//For IntiR_KernelCG
	//Scalar CoefficientC = 2 * CONSTANT_PARAMETER.Dimension / ( CONSTANT_PARAMETER.LambdaValue * CONSTANT_PARAMETER.ConstantParticleDensity);
	Scalar CoefficientCIntiR;

	//Scalar PressureCoefficient = (m_pParameter->LambdaValue * m_pParameter->ConstantParticleDensity * m_pParameter->Density) 
	//	/ ( 2 * m_pParameter->Dimension * m_pParameter->Dt * m_pParameter->Dt);
	Scalar PressureCoefficient;
		

};


struct DeltaTime
{
	Scalar Dt;
	Scalar CourantNumber;
	Scalar InflowStep;
	Scalar OutputStep;
};
__device__ __host__ __inline__ Scalar3 operator+(const Scalar3 &V1 , const Scalar3 &V2)
{
	Scalar3 tmp;
	tmp.x = V1.x + V2.x;
	tmp.y = V1.y + V2.y;
	tmp.z = V1.z + V2.z;
	return tmp;
}
__device__ __host__ __inline__ Scalar3 operator-(const Scalar3 &V1 , const Scalar3 &V2)
{
	Scalar3 tmp;
	tmp.x = V1.x - V2.x;
	tmp.y = V1.y - V2.y;
	tmp.z = V1.z - V2.z;
	return tmp;
}

static __inline__ __host__ __device__ Scalar3 operator*(const Scalar3& V, const Scalar & S)
{
	Scalar3 tmp;
	tmp.x = V.x * S;
	tmp.y = V.y * S;
	tmp.z = V.z * S;
	return tmp;
}
static __inline__ __host__ __device__ Scalar3 operator/(const Scalar3& V, const Scalar& S)
{
	Scalar3 tmp;
	tmp.x = V.x / S;
	tmp.y = V.y / S;
	tmp.z = V.z / S;
	return tmp;
}

#define CCT_ERROR_CHECK(StatusType) if(CCT_NOERR != StatusType)	{printf("\nError in file : %s in Line : %d:: Function : %s",__FILE__,__LINE__,__FUNCTION__); return StatusType;}
#define ET_ERROR_CHECK(StatusType) if(CCT_NOERR != StatusType) {return;}
#define ET_ERROR_BREAK(StatusType) if(CCT_NOERR != StatusType) {break;}
