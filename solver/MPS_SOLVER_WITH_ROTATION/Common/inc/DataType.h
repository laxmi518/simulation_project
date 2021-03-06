#pragma once
#pragma once
#include <cassert>
#include <iostream>
#include <cmath>


typedef float	Scalar;
typedef int		Integer;
typedef bool	Bool;

typedef int Integer;


#ifdef __CUDACC__
#define ET_INLINE static __inline__ __host__ __device__

typedef float3 Scalar3;
typedef int3 Integer3;
typedef float2 Scalar2;
typedef int2 Integer2;
#else
#define ET_INLINE static inline
struct Scalar2
{
	Scalar x;
	Scalar y;
};
struct Scalar3      /*create a user defined data type Scalar3*/
{
	Scalar x;
	Scalar y; 
	Scalar z;
};
struct Integer3
{
	Integer x;
	Integer y;
	Integer z;
};
#endif


/*********************************************/
// vector operation for 2D
/*********************************************/
ET_INLINE Scalar2 operator+(const Scalar2& A, const Scalar2& B)
{
	Scalar2 ans;
	ans.x = A.x + B.x;
	ans.y = A.y + B.y;
	return ans;
}

ET_INLINE Scalar2 operator-(const Scalar2& A, const Scalar2& B)
{
	Scalar2 ans;
	ans.x = A.x - B.x;
	ans.y = A.y - B.y;
	return ans;
}

ET_INLINE Scalar2 operator*(const Scalar A, const Scalar2& B)
{
	Scalar2 ans;
	ans.x = A * B.x;
	ans.y = A * B.y;
	return ans;
}

ET_INLINE Scalar2 operator*(const Scalar2& B, const Scalar A)
{
	Scalar2 ans;
	ans.x = A * B.x;
	ans.y = A * B.y;
	return ans;
}


ET_INLINE Scalar2 operator/(const Scalar2& B, const Scalar A)
{
	Scalar2 ans;
	ans.x = B.x / A;
	ans.y = B.y / A;
	return ans;
}

ET_INLINE Scalar DotProduct(const Scalar2& A, const Scalar2& B)
{
	return A.x * B.x + A.y * B.y;
}

//inline Scalar CrossProduct(const Scalar2& A, const Scalar2& B)
//{
//	return A.x * B.y - A.y * B.x;
//}

ET_INLINE Scalar3 operator+(const Scalar3& A, const Scalar3& B)
{
	Scalar3 ans;
	ans.x = A.x + B.x;
	ans.y = A.y + B.y;
	ans.z = A.z + B.z;
	return ans;
}

ET_INLINE Scalar3 operator-(const Scalar3& A, const Scalar3& B)
{
	Scalar3 ans;
	ans.x = A.x - B.x;
	ans.y = A.y - B.y;
	ans.z = A.z - B.z;
	return ans;
}

ET_INLINE Scalar3 operator*(const Scalar A, const Scalar3& B)
{
	Scalar3 ans;
	ans.x = A * B.x;
	ans.y = A * B.y;
	ans.z = A * B.z;
	return ans;
}
ET_INLINE Scalar3 operator*(const Scalar3& B, const Scalar A)
{
	Scalar3 ans;
	ans.x = A * B.x;
	ans.y = A * B.y;
	ans.z = A * B.z;
	return ans;
}

ET_INLINE Scalar3 operator/(const Scalar3& B, const Scalar A)
{
	Scalar3 ans;
	ans.x = B.x / A;
	ans.y = B.y / A;
	ans.z = B.z / A;
	return ans;
}

ET_INLINE Scalar SquaredMagnitude(const Scalar3& A)
{
	return A.x * A.x + A.y * A.y + A.z * A.z;
}

ET_INLINE Scalar Magnitude(const Scalar3& A)
{
	return sqrt(SquaredMagnitude(A));
}

ET_INLINE void Normalize(Scalar3& A)
{
	Scalar Length = Magnitude(A);
	A.x /= Length;
	A.y /= Length;
	A.z /= Length;
}

ET_INLINE Scalar DotProduct(const Scalar3& A, const Scalar3& B)
{
	return (A.x * B.x + A.y * B.y + A.z * B.z);
}

ET_INLINE Scalar3 CrossProduct(const Scalar3& A, const Scalar3& B)
{
	Scalar3 ans;
	ans.x = (A.y * B.z - A.z * B.y);
	ans.y = (A.z * B.x - A.x * B.z);
	ans.z = (A.x * B.y - A.y * B.x);
	return ans;
}

ET_INLINE Scalar3 make_Scalar3(Scalar x, Scalar y, Scalar z)
{
  Scalar3 t; t.x = x; t.y = y; t.z = z; return t;
}
ET_INLINE Integer3 make_Integer3(Integer x, Integer y, Integer z)
{
  Integer3 t; t.x = x; t.y = y; t.z = z; return t;
}


//typedef double	Scalar;


static const Scalar CONSTANT_CIRCULAR = static_cast<Scalar>(3.141592653589793);
static const Scalar RADDEG =  static_cast<Scalar>(57.29578);
#define EPS  1e-8


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
#define BLOCK_MAX_DIM 256
//#define PMDEBUG
struct CTriangle
{
	Scalar3 Normal;
	Scalar3 Vertex[3];
	Scalar3 Velocity;
	Scalar Temperature;
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

#define CELLMAXTRIANGLENUM 100
#define CELLMAXPARTICLENUM 729
#define OFFSET 2
#define GHOST_OFFSET 3
#define PRESSURE_OFFSET 1
struct CCell
{
	Integer m_ParticleID;// ID of the particle included in the cell, there can be multiple particles in a cell, so there is a Hash table in CCellSet class, to store more than one particle.
	Integer m_HashID;// ID of the Hash to be filled, this is like a temporary variable, only used for filling the particles in the cell.
	unsigned int m_TriangleNum;
	Integer m_TriangleID[CELLMAXTRIANGLENUM];	
};

// 流体密度と粒子数密度は同じなので要削除
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
	Scalar FixedWallTemperature;
	Scalar MovableWallTemperature;
	
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
};
struct DeltaTime
{
	Scalar Dt;
	Scalar CourantNumber;
	Scalar InflowStep;
	Scalar OutputStep;
};

/*********************************************/
// data structure
/*********************************************/
struct Circle2D
{
	Scalar2 Center;
	Scalar2 LocalX;	///< always (1, 0)
	Scalar2 LocalY;	///< always (0, 1)
	Scalar2 Normal; ///< always (0, 0) because the direction is +Z
	Scalar Radius;
};

inline void Copy(const Circle2D& Src, Circle2D& Dst)
{
	Dst.Center = Src.Center;
	Dst.LocalX = Src.LocalX;
	Dst.LocalY = Src.LocalY;
	Dst.Normal = Src.Normal;
	Dst.Radius = Src.Radius;
}

inline Scalar2 CalcCoord(const Circle2D& Circle, Scalar Radian)
{
	Scalar2 Coord = Circle.Center + Circle.Radius * cos(Radian) * Circle.LocalX + Circle.Radius * sin(Radian) * Circle.LocalY;
	return Coord;
}

struct Circle3D
{
	Scalar3 Center;
	Scalar3 LocalX;
	Scalar3 LocalY;
	Scalar3 Normal;
	Scalar Radius;
};

inline void Copy(const Circle3D& Src, Circle3D& Dst)
{
	Dst.Center = Src.Center;
	Dst.LocalX = Src.LocalX;
	Dst.LocalY = Src.LocalY;
	Dst.Normal = Src.Normal;
	Dst.Radius = Src.Radius;
}

inline Scalar3 CalcLocalXFromNormal(const Scalar3& Normal)
{
	assert(!(EPS > fabs(Normal.x) && EPS > fabs(Normal.y) && EPS > fabs(Normal.z)));

	double NxAbs = fabs(Normal.x); 
	double NyAbs = fabs(Normal.y); 
	double NzAbs = fabs(Normal.z); 
	
	Scalar3 LocalX;
	if(NyAbs <= NxAbs && NyAbs <= NzAbs) 
	{
		if(NxAbs > NzAbs) 
		{ 
			LocalX.x = -Normal.z;
			LocalX.y = 0;
			LocalX.z = Normal.x;
		} 
		else 
		{ 
			LocalX.x = Normal.z;
			LocalX.y = 0;
			LocalX.z = -Normal.x;
		} 
	} 
	else if(NxAbs <= NyAbs && NxAbs <= NzAbs) 
	{ 
		if (NyAbs > NzAbs) 
		{ 
			LocalX.x = 0;
			LocalX.y = -Normal.z;
			LocalX.z = Normal.y;
		} 
		else 
		{ 
			LocalX.x = 0;
			LocalX.y = Normal.z;
			LocalX.z = -Normal.y;
		} 
	} 
	else 
	{ 
		if(NxAbs > NyAbs) 
		{ 
			LocalX.x = -Normal.y;
			LocalX.y = Normal.x;
			LocalX.z = 0;
		} 
		else 
		{ 
			LocalX.x = Normal.y;
			LocalX.y = -Normal.x;
			LocalX.z = 0;
		}
	}
	Normalize(LocalX);
	return LocalX;
}

inline Circle3D* Create(const Scalar3& Center, const Scalar3& Normal, Scalar Radius)
{
	if(EPS > fabs(Normal.x) && EPS > fabs(Normal.y) && EPS > fabs(Normal.z))
	{
		return NULL;
	}
	Circle3D* Circle = new Circle3D();
	Circle->Center = Center;
	Circle->Normal = Normal;
	Normalize(Circle->Normal);
	Circle->LocalX = CalcLocalXFromNormal(Circle->Normal);
	Circle->LocalY = CrossProduct(Circle->Normal, Circle->LocalX);
	Normalize(Circle->LocalY);
	Circle->Radius = Radius;
	return Circle;
}

inline Circle3D* Create(const Scalar3& Center, const Scalar3& LocalX, const Scalar3& Normal, Scalar Radius)
{
	if(EPS > fabs(LocalX.x) && EPS > fabs(LocalX.y) && EPS > fabs(LocalX.z))
	{
		return NULL;
	}
	if(EPS > fabs(Normal.x) && EPS > fabs(Normal.y) && EPS > fabs(Normal.z))
	{
		return NULL;
	}
	Circle3D* Circle = new Circle3D();
	Circle->Center = Center;
	Circle->LocalX = LocalX;
	Circle->LocalY = CrossProduct(Normal, LocalX);
	Circle->Normal = Normal;
	Circle->Radius = Radius;
	return Circle;
}

inline Scalar3 CalcCoord(const Circle3D& Circle, Scalar Radian)
{
	Scalar3 Coord = Circle.Center + Circle.Radius * cos(Radian) * Circle.LocalX + Circle.Radius * sin(Radian) * Circle.LocalY;
	return Coord;
}

inline bool IsValid(const Circle3D& Circle)
{
	/*if(EPS > fabs(Circle.LocalX.x) && EPS > fabs(Circle.LocalX.y) && EPS > fabs(Circle.LocalX.z))
	{
		return false;
	}
	if(EPS > fabs(Circle.LocalY.x) && EPS > fabs(Circle.LocalY.y) && EPS > fabs(Circle.LocalY.z))
	{
		return false;
	}
	if(EPS > fabs(Circle.Normal.x) && EPS > fabs(Circle.Normal.y) && EPS > fabs(Circle.Normal.z))
	{
		return false;
	}
	return EPS > Circle.Radius;*/
	return EPS < Circle.Radius;
}

inline Scalar3 CalcCircumCenter(const Scalar3& A, const Scalar3& B, const Scalar3& C)
{
	assert(EPS < SquaredMagnitude(A - B));
	assert(EPS < SquaredMagnitude(B - C));
	assert(EPS < SquaredMagnitude(C - A));

	Scalar3 AC = A - C;
	Scalar3 BC = B - C;
	Scalar TriangleArea = Magnitude(CrossProduct(AC, BC)) / 2;

	Scalar ABSquaredLength = SquaredMagnitude(A - B);
	Scalar BCSquaredLength = SquaredMagnitude(BC);
	Scalar CASquaredLength = SquaredMagnitude(C - A);

	Scalar3 CircumCenter = (
	BCSquaredLength * (CASquaredLength + ABSquaredLength - BCSquaredLength) * A +
	CASquaredLength * (ABSquaredLength + BCSquaredLength - CASquaredLength) * B +
	ABSquaredLength * (BCSquaredLength + CASquaredLength - ABSquaredLength) * C) /
	(16 * TriangleArea * TriangleArea);
	return CircumCenter;
}




//__device__ __host__ __inline__ Scalar3 operator+(const Scalar3 &V1 , const Scalar3 &V2)
//{
//	Scalar3 tmp;
//	tmp.x = V1.x + V2.x;
//	tmp.y = V1.y + V2.y;
//	tmp.z = V1.z + V2.z;
//	return tmp;
//}
//__device__ __host__ __inline__ Scalar3 operator-(const Scalar3 &V1 , const Scalar3 &V2)
//{
//	Scalar3 tmp;
//	tmp.x = V1.x - V2.x;
//	tmp.y = V1.y - V2.y;
//	tmp.z = V1.z - V2.z;
//	return tmp;
//}
//
//static __inline__ __host__ __device__ Scalar3 operator*(const Scalar3& V, const Scalar & S)
//{
//	Scalar3 tmp;
//	tmp.x = V.x * S;
//	tmp.y = V.y * S;
//	tmp.z = V.z * S;
//	return tmp;
//}
//static __inline__ __host__ __device__ Scalar3 operator/(const Scalar3& V, const Scalar& S)
//{
//	Scalar3 tmp;
//	tmp.x = V.x / S;
//	tmp.y = V.y / S;
//	tmp.z = V.z / S;
//	return tmp;
//}

#define CCT_ERROR_CHECK(StatusType) if(CCT_NOERR != StatusType)	{return StatusType;}
#define ET_ERROR_CHECK(StatusType) if(CCT_NOERR != StatusType) {return;}
#define ET_ERROR_BREAK(StatusType) if(CCT_NOERR != StatusType) {break;}
