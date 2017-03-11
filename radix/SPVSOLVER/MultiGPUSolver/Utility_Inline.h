//http://www.devx.com/tips/Tip/42853
#pragma once
#define SQRT_MAGIC_F 0x5f3759df

//static __inline__ __device__ __host__  float  sqrt_utility(const float x)
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
static __inline__ __device__ __host__ float sqrt_utility ( float x){
  float xhalf = 0.5f*x;
  int i = *(int*)&x;
  i = 0x5f3759df - (i>>1);
  x = *(float*)&i;
  x = x*(1.5f - xhalf*x*x);
  return 1/x;
}

static __inline__ __device__ __host__  Scalar CalcDistance(const Scalar3* const lhs, const Scalar3* const rhs)
{
	return sqrt_utility((lhs->x - rhs->x) * (lhs->x - rhs->x) + (lhs->y - rhs->y) * (lhs->y - rhs->y) + (lhs->z - rhs->z) * (lhs->z - rhs->z));
}
static __inline__ __device__ __host__ Scalar InnerProduct(const Scalar3* const Vector)
{
	return Vector->x * Vector->x + Vector->y * Vector->y + Vector->z * Vector->z;
}
//// Magnitude of a vector
static __inline__ __device__ __host__ Scalar Magnitude(const Scalar3* const Vector)
{
	Scalar Product =  Vector->x * Vector->x + Vector->y * Vector->y + Vector->z * Vector->z;
	return sqrt_utility(Product);
}
static __inline__ __device__ __host__ Scalar DotProduct(const Scalar3* const P1, const Scalar3* const P2)
{
	Scalar val = P1->x * P2->x + P1->y * P2->y + P1->z * P2->z;
	return val;
}
//Cross Product of vector
static __inline__ __device__ __host__ Scalar3 CrossProduct(const Scalar3& A, const Scalar3& B)
{
	Scalar3 ans;
	ans.x = (A.y * B.z - A.z * B.y);
	ans.y = (A.z * B.x - A.x * B.z);
	ans.z = (A.x * B.y - A.y * B.x);
	return ans;
}
static __inline__ __device__ __host__ Scalar SegmentRatio(const Scalar3* const P1, const Scalar3* const P2, const Scalar3* const P3)
{
	Scalar3 DenomVec = *P2 - *P1;//SubtractVector(P2, P1);

	Scalar Denom = DenomVec.x * DenomVec.x + DenomVec.y * DenomVec.y + DenomVec.z * DenomVec.z;


	//Scalar Denom = InnerProduct(&DenomVec);
	//Scalar Numer = (P3->x - P1->x) * (P2->x - P1->x) + (P3->y - P1->y) * (P2->y - P1->y) + (P3->z - P1->z) * (P2->z - P1->z);
	if(Denom == 0.0)
	{
		Denom = static_cast<Scalar>(1e-5);
	}

		//Laxmi 
	//Scalar u = Numer / Denom;

	Scalar u = (P3->x - P1->x) * (P2->x - P1->x) + (P3->y - P1->y) * (P2->y - P1->y) + (P3->z - P1->z) * (P2->z - P1->z) / Denom;
	return u;
}
//http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/

//static __inline__ __device__ __host__ Scalar3 SubtractVector(const Scalar3* const P1, const Scalar3* const P2)
//{
//	Scalar3 tmp;
//	tmp.x = P1->x - P2->x;
//	tmp.y = P1->y - P2->y;
//	tmp.z = P1->z - P2->z;
//	return tmp;
//}
static __inline__ __device__ __host__  Scalar3 CalcDistanceVector(const Scalar3* const lhs, const Scalar3* const rhs)
{
	return make_Scalar3((lhs->x - rhs->x) , (lhs->y - rhs->y),  (lhs->z - rhs->z));
}
//http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
static __inline__ __device__ __host__ Scalar3 CalcDistanceVectorLine(const Scalar3* const L0,const Scalar3* const L1, const Scalar3* const Position)
{
	Scalar u = SegmentRatio(L0, L1, Position);
	Scalar3 DistanceVector;
	if(u >= 0 && u <= 1)
	{
		Scalar3 myVec = *L1 - *L0;//SubtractVector(L1, L0);
		Scalar3 P;
		P.x = L0->x + u * myVec.x;
		P.y = L0->y + u * myVec.y;
		P.z = L0->z + u * myVec.z;
		DistanceVector = *Position - P;//SubtractVector(Position, &P);
	}
	if( u < 0)
	{
		DistanceVector = *L1 - *Position;//CalcDistanceVector(L1, Position);
	}
	if( u > 1)
	{
		DistanceVector = *L0 - *Position;//CalcDistanceVector(L0, Position);
	}
	return DistanceVector;
}
static __inline__ __device__ __host__ Scalar3 CalcNormal(const Scalar3* const P1, const Scalar3* const P2)
{
	Scalar mag = CalcDistance(P2,P1);
	Scalar3 val = make_Scalar3(0,0,0);
	if(mag > 0)
	{
		val = make_Scalar3(((P2->x - P1->x)/mag),((P2->y - P1->y)/mag),((P2->z - P1->z)/mag));
	}
	return val;
}
static __inline__ __device__ __host__ Scalar3 CalcNormalDistance(const CTriangle* const pTriangle, const Scalar3* const Position)
{
	const Scalar3 Normal = pTriangle->Normal;
	//laxmi
	//Scalar Center[3];
	//Center[0] = 0;
	//Center[1] = 0;
	//Center[2] = 0;
	///*for(Integer i = 0; i < 3; ++i)
	//{
	//	Center[0] += pTriangle->Vertex[i].x;
	//	Center[1] += pTriangle->Vertex[i].y;
	//	Center[2] += pTriangle->Vertex[i].z;
	//}*/

	//Center[0] = (pTriangle->Vertex[0].x+pTriangle->Vertex[1].x+pTriangle->Vertex[2].x);
	//	Center[1] = (pTriangle->Vertex[0].y+pTriangle->Vertex[1].y+pTriangle->Vertex[2].y);
	//	Center[2] =(pTriangle->Vertex[0].z+pTriangle->Vertex[1].z+pTriangle->Vertex[2].z);
	//Center[0] /= 3.0;
	//Center[1] /= 3.0;
	//Center[2] /= 3.0;

	//const Scalar lNumerator = Normal.x * Center[0] + Normal.y * Center[1] + Normal.z * Center[2];
	//const Scalar rNumerator = Normal.x * Position->x + Normal.y * Position->y + Normal.z * Position->z;
	//const Scalar denominator = Normal.x * Normal.x + Normal.y * Normal.y + Normal.z * Normal.z;


	Scalar3 Center;
	Center.x = 0;
	Center.y = 0;
	Center.z = 0;
	//for(Integer i = 0; i < 3; ++i)
	//{
	//	Center.x += pTriangle->Vertex[i].x;
	//	Center.y += pTriangle->Vertex[i].y;
	//	Center.z += pTriangle->Vertex[i].z;
	//}

	Center.x = (pTriangle->Vertex[0].x+pTriangle->Vertex[1].x+pTriangle->Vertex[2].x)/3.0;
	Center.y = (pTriangle->Vertex[0].y+pTriangle->Vertex[1].y+pTriangle->Vertex[2].y)/3.0;
	Center.z =(pTriangle->Vertex[0].z+pTriangle->Vertex[1].z+pTriangle->Vertex[2].z)/3.0;

	/*Center.x /= 3.0;
	Center.y /= 3.0;
	Center.z /= 3.0;*/

	/*const Scalar lNumerator = Normal.x * Center.x + Normal.y * Center.y + Normal.z * Center.z;
	const Scalar rNumerator = Normal.x * Position->x + Normal.y * Position->y + Normal.z * Position->z;
	const Scalar denominator = Normal.x * Normal.x + Normal.y * Normal.y + Normal.z * Normal.z;*/
		//const Scalar sDenom = sqrt(denominator);
	
	
	const Scalar sDenom = sqrt_utility(Normal.x * Normal.x + Normal.y * Normal.y + Normal.z * Normal.z);
	Scalar3 DistanceVector = make_Scalar3(0,0,0) ;
	if(sDenom == 0.0)
	{
		return DistanceVector;
	}
	//const Scalar Coefficient = (rNumerator - lNumerator) / sDenom;

	const Scalar Coefficient = ((Normal.x * Position->x + Normal.y * Position->y + Normal.z * Position->z) - (Normal.x * Center.x + Normal.y * Center.y + Normal.z * Center.z)) / sDenom;

	DistanceVector.x = Coefficient / sDenom * Normal.x;
	DistanceVector.y = Coefficient / sDenom * Normal.y;
	DistanceVector.z = Coefficient / sDenom * Normal.z;

	return DistanceVector;
}
//http://www.blackpawn.com/texts/pointinpoly/default.html
static __inline__ __device__ __host__ Scalar CalcDistance(const CTriangle* const pTriangle, const Scalar3* const Position,Scalar3* const DistanceVector)
{
	Scalar3 dv = CalcNormalDistance(pTriangle,Position);
	Scalar3 P = make_Scalar3(Position->x - dv.x, Position->y - dv.y, Position->z - dv.z);
	// Compute vectors
	Scalar3 V0 = pTriangle->Vertex[2] - pTriangle->Vertex[0];//SubtractVector(&pTriangle->Vertex[2], &pTriangle->Vertex[0]);
	Scalar3 V1 = pTriangle->Vertex[1] - pTriangle->Vertex[0];//SubtractVector(&pTriangle->Vertex[1], &pTriangle->Vertex[0]);
	Scalar3 V2 = P - pTriangle->Vertex[0];//SubtractVector(&P, &pTriangle->Vertex[0]);
	// Compute dot products
	Scalar dot00 = DotProduct(&V0, &V0);
	Scalar dot01 = DotProduct(&V0, &V1);
	Scalar dot02 = DotProduct(&V0, &V2);
	Scalar dot11 = DotProduct(&V1, &V1);
	Scalar dot12 = DotProduct(&V1, &V2);
	// compute barycentric coordinates

	Scalar invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	Scalar u = (dot11 * dot02 - dot01 * dot12 ) * invDenom;
	Scalar v = (dot00 * dot12 - dot01 * dot02 ) * invDenom;

	if(u >= 0 && v >= 0 && u + v > 1)
	{
		dv = CalcDistanceVectorLine(&pTriangle->Vertex[1],&pTriangle->Vertex[2], Position);
	}
	if(u < 0 )
	{
		dv = CalcDistanceVectorLine(&pTriangle->Vertex[0],&pTriangle->Vertex[1], Position);
	}
	if(v < 0 )
	{
		dv = CalcDistanceVectorLine(&pTriangle->Vertex[0],&pTriangle->Vertex[2], Position);
	}
	if(DistanceVector)
	{
		DistanceVector->x = dv.x;
		DistanceVector->y = dv.y;
		DistanceVector->z = dv.z;
	}		
	return Magnitude(&dv);	
}
static __inline__ __device__ __host__ __device__ Scalar GetWeight(const Scalar Re, const Scalar Distance)
{
	if(Distance > Re)
	{
		return 0;
	}

	if(Distance > 0.0)
	{
		return  Re / Distance - 1;
	}
	return 1e6;	
}
// addition of new weightage function (1 - r/Re) for innerpressure type calculation
static __inline__ __device__ __host__ __device__ Scalar GetSurfaceWeight(const Scalar Re, const Scalar Distance)
{
	if(Distance > Re)
	{
		return 0;
	}
	if(Distance > 0.0)
	{
		return  1 - Distance / Re;
	}
	return 1e6;
}

// Z関数の取得(勾配計算時の壁寄与分の密度）, r_e = 2.1 * l_0
static __inline__ __device__ __host__ Scalar GetZValueFunctionGradient(const Scalar InterpolatedDistance, const Scalar InitialParticleDistance)
{
	const Scalar InterpolatedPerInitial = InterpolatedDistance / InitialParticleDistance;
	Scalar a = 0.0, b = 0.0, c = 0.0;
	if(InterpolatedPerInitial <= 0.1)
	{
		a = static_cast<Scalar>(-5.0474610E+03);
		b = static_cast<Scalar>(0.1);
		c = static_cast<Scalar>(2.9725900E+01);
	}
	else if(InterpolatedPerInitial <= 0.2)
	{
		a = static_cast<Scalar>(-1.1272200E+02);
		b = static_cast<Scalar>(0.2);
		c = static_cast<Scalar>(1.8453700E+01);
	}
	else if(InterpolatedPerInitial <= 0.3)
	{
		a = static_cast<Scalar>(-4.3342000E+01);
		b = static_cast<Scalar>(0.3);
		c = static_cast<Scalar>(1.4119500E+01);
	}
	else if(InterpolatedPerInitial <= 0.4)
	{
		a = static_cast<Scalar>(-2.6328000E+01);
		b = static_cast<Scalar>(0.4);
		c = static_cast<Scalar>(1.1486700E+01);
	}
	else if(InterpolatedPerInitial <= 0.5)
	{
		a = static_cast<Scalar>(-1.9648400E+01);
		b = static_cast<Scalar>(0.5);
		c = static_cast<Scalar>(9.5218600E+00);
	}
	else if(InterpolatedPerInitial <= 0.6)
	{
		a = static_cast<Scalar>(-1.5636200E+01);
		b = static_cast<Scalar>(0.6);
		c = static_cast<Scalar>(7.9582400E+00);
	}
	else if(InterpolatedPerInitial <= 0.7)
	{
		a = static_cast<Scalar>(-1.2591500E+01);
		b = static_cast<Scalar>(0.7);
		c = static_cast<Scalar>(6.6990900E+00);
	}
	else if(InterpolatedPerInitial <= 0.8)
	{
		a = static_cast<Scalar>(-1.0988800E+01);
		b = static_cast<Scalar>(0.8);
		c = static_cast<Scalar>(5.6002100E+00);
	}
	else if(InterpolatedPerInitial <= 0.9)
	{
		a = static_cast<Scalar>(-9.0692000E+00);
		b = static_cast<Scalar>(0.9);
		c = static_cast<Scalar>(4.6932900E+00);
	}
	else if(InterpolatedPerInitial <= 1.0)
	{
		a =static_cast<Scalar>(-7.5385000E+00);
		b = static_cast<Scalar>(1);
		c = static_cast<Scalar>(3.9394400E+00);
	}
	else if(InterpolatedPerInitial <= 1.1)
	{
		a = static_cast<Scalar>(-6.9148000E+00);
		b = static_cast<Scalar>(1.1);
		c =static_cast<Scalar>(3.2479600E+00);
	}
	else if(InterpolatedPerInitial <= 1.2)
	{
		a = static_cast<Scalar>(-5.9144000E+00);
		b = static_cast<Scalar>(1.2);
		c = static_cast<Scalar>(2.6565200E+00);
	}
	else if(InterpolatedPerInitial <= 1.3)
	{
		a = static_cast<Scalar>(-5.4670000E+00);
		b = static_cast<Scalar>(1.3);
		c = static_cast<Scalar>(2.1098200E+00);
	}
	else if(InterpolatedPerInitial <= 1.4)
	{
		a = static_cast<Scalar>(-5.0626000E+00);
		b = static_cast<Scalar>(1.4);
		c = static_cast<Scalar>(1.6035600E+00);
	}
	else if(InterpolatedPerInitial <= 1.5)
	{
		a = static_cast<Scalar>(-4.6948000E+00);
		b = static_cast<Scalar>(1.5);
		c = static_cast<Scalar>(1.1340800E+00);
	}
	else if(InterpolatedPerInitial <= 1.6)
	{
		a = static_cast<Scalar>(-3.6959000E+00);
		b = static_cast<Scalar>(1.6);
		c = static_cast<Scalar>(7.6449000E-01);
	}
	else if(InterpolatedPerInitial <= 1.7)
	{
		a = static_cast<Scalar>(-2.7022800E+00);
		b = static_cast<Scalar>(1.7);
		c = static_cast<Scalar>(4.9426200E-01);
	}
	else if(InterpolatedPerInitial <= 1.8)
	{
		a = static_cast<Scalar>(-2.4819500E+00);
		b = static_cast<Scalar>(1.8);
		c = static_cast<Scalar>(2.4606700E-01);
	}
	else if(InterpolatedPerInitial <= 1.9)
	{
		a = static_cast<Scalar>(-1.4080400E+00);
		b = static_cast<Scalar>(1.9);
		c = static_cast<Scalar>(1.0526300E-01);
	}
	else if(InterpolatedPerInitial <= 2.0)
	{
		a = static_cast<Scalar>(-5.5263200E-01);
		b = static_cast<Scalar>(2);
		c = static_cast<Scalar>(4.9999800E-02);
	}
	else if(InterpolatedPerInitial <= 2.1)
	{
		// 誤差を考慮して2.15くらいまで
		a = static_cast<Scalar>(-4.9999800E-01);
		b = static_cast<Scalar>(2.1);
		c = static_cast<Scalar>(0.0000000E+00);
	}
	else
	{
		return 0.0;
	}
	Scalar Val =  a * (InterpolatedPerInitial - b) + c;
	return Val;
}
static __inline__ __device__ __host__ Scalar GetZValueFunctionGradient2D(const Scalar InterpolatedDistance, const Scalar InitialParticleDistance)
{
	const Scalar InterpolatedPerInitial = InterpolatedDistance / InitialParticleDistance;
	Scalar a = 0.0, b = 0.0, c = 0.0;
	if(InterpolatedPerInitial <= 0.01)
	{//0	<x/l_0=<	0.01	-9.996E+05	0.01	3.837E+02

		a = static_cast<Scalar>(-9.996E+05);
		b = static_cast<Scalar>(0.01);
		c = static_cast<Scalar>(3.837E+02);
	}
	else if(InterpolatedPerInitial <= 0.05)
	{//0.01	<x/l_0=<	0.05	-3.027E+02	0.05	8.101E+01

		a = static_cast<Scalar>(-3.027E+02);
		b = static_cast<Scalar>(0.05);
		c = static_cast<Scalar>(8.101E+01);
	}
	else if(InterpolatedPerInitial <= 0.1)
	{//0.05	<x/l_0=<	0.1	-3.812E+01	0.1	4.289E+01

		a = static_cast<Scalar>(3.812E+01);
		b = static_cast<Scalar>(0.1);
		c = static_cast<Scalar>(4.289E+01);
	}
	else if(InterpolatedPerInitial <= 0.2)
	{//0.1	<x/l_0=<	0.2	-1.950E+01	0.2	2.339E+01

		a = static_cast<Scalar>(-1.950E+01);
		b = static_cast<Scalar>(0.2);
		c = static_cast<Scalar>(2.339E+01);
	}
	else if(InterpolatedPerInitial <= 0.3)
	{//0.2	<x/l_0=<	0.3	-6.917E+00	0.3	1.647E+01

		a = static_cast<Scalar>(-6.917E+00);
		b = static_cast<Scalar>(0.3);
		c = static_cast<Scalar>(1.647E+01);
	}
	else if(InterpolatedPerInitial <= 0.4)
	{//0.3	<x/l_0=<	0.4	-3.779E+00	0.4	1.269E+01

		a = static_cast<Scalar>(-3.779E+00);
		b = static_cast<Scalar>(0.4);
		c = static_cast<Scalar>(1.269E+01);
	}
	else if(InterpolatedPerInitial <= 0.5)
	{//0.4	<x/l_0=<	0.5	-2.523E+00	0.5	1.017E+01

		a = static_cast<Scalar>(-2.523E+00);
		b = static_cast<Scalar>(0.5);
		c = static_cast<Scalar>(1.017E+01);
	}
	else if(InterpolatedPerInitial <= 0.6)
	{//0.5	<x/l_0=<	0.6	-1.889E+00	0.6	8.282E+00

		a = static_cast<Scalar>(-1.889E+00);
		b = static_cast<Scalar>(0.6);
		c = static_cast<Scalar>(8.282E+00);
	}
	else if(InterpolatedPerInitial <= 0.7)
	{//0.6	<x/l_0=<	0.7	-1.484E+00	0.7	6.798E+00

		a = static_cast<Scalar>(-1.484E+00);
		b = static_cast<Scalar>(0.7);
		c = static_cast<Scalar>(6.798E+00);
	}
	else if(InterpolatedPerInitial <= 0.8)
	{//0.7	<x/l_0=<	0.8	-1.214E+00	0.8	5.584E+00

		a = static_cast<Scalar>(-1.214E+00);
		b = static_cast<Scalar>( 0.8);
		c = static_cast<Scalar>(5.584E+00);
	}
	else if(InterpolatedPerInitial <= 0.9)
	{//0.8	<x/l_0=<	0.9	-9.752E-01	0.9	4.609E+00

		a = static_cast<Scalar>(-9.752E-01);
		b = static_cast<Scalar>(0.9);
		c = static_cast<Scalar>(4.609E+00);
	}
	else if(InterpolatedPerInitial <= 1.0)
	{//0.9	<x/l_0=<	1	-7.931E-01	1	3.816E+00

		a = static_cast<Scalar>(-7.931E-01);
		b = static_cast<Scalar>(1);
		c = static_cast<Scalar>(3.816E+00);
	}
	else if(InterpolatedPerInitial <= 1.1)
	{//1	<x/l_0=<	1.1	-6.940E-01	1.1	3.122E+00

		a = static_cast<Scalar>(-6.940E-01);
		b = static_cast<Scalar>(1.1);
		c = static_cast<Scalar>(3.122E+00);
	}
	else if(InterpolatedPerInitial <= 1.2)
	{//1.1	<x/l_0=<	1.2	-5.320E-01	1.2	2.590E+00

		a = static_cast<Scalar>(-5.320E-01);
		b = static_cast<Scalar>(1.2);
		c = static_cast<Scalar>(2.590E+00);
	}
	else if(InterpolatedPerInitial <= 1.3)
	{//1.2	<x/l_0=<	1.3	-4.727E-01	1.3	2.117E+00

		a = static_cast<Scalar>(-4.727E-01);
		b = static_cast<Scalar>(1.3);
		c = static_cast<Scalar>(2.117E+00);
	}
	else if(InterpolatedPerInitial <= 1.4)
	{//1.3	<x/l_0=<	1.4	-4.229E-01	1.4	1.694E+00

		a = static_cast<Scalar>(-4.229E-01);
		b = static_cast<Scalar>(1.4);
		c = static_cast<Scalar>(1.694E+00);
	}
	else if(InterpolatedPerInitial <= 1.5)
	{//1.4	<x/l_0=<	1.5	-3.806E-01	1.5	1.314E+00

		a = static_cast<Scalar>(-3.806E-01);
		b = static_cast<Scalar>(1.5);
		c = static_cast<Scalar>(1.314E+00);
	}
	else if(InterpolatedPerInitial <= 1.6)
	{//1.5	<x/l_0=<	1.6	-3.442E-01	1.6	9.693E-01

		a = static_cast<Scalar>(-3.442E-01);
		b = static_cast<Scalar>(1.6);
		c = static_cast<Scalar>(9.693E-01);
	}
	else if(InterpolatedPerInitial <= 1.7)
	{//1.6	<x/l_0=<	1.7	-3.127E-01	1.7	6.566E-01

		a = static_cast<Scalar>(-3.127E-01);
		b = static_cast<Scalar>(1.7);
		c = static_cast<Scalar>(6.566E-01);
	}
	else if(InterpolatedPerInitial <= 1.8)
	{//1.7	<x/l_0=<	1.8	-2.851E-01	1.8	3.715E-01

		a = static_cast<Scalar>(-2.851E-01);
		b = static_cast<Scalar>(1.8);
		c = static_cast<Scalar>(3.715E-01);
	}
	else if(InterpolatedPerInitial <= 1.9)
	{//1.8	<x/l_0=<	1.9	-1.820E-01	1.9	1.895E-01

		a = static_cast<Scalar>(-1.820E-01);
		b = static_cast<Scalar>(1.9);
		c = static_cast<Scalar>(1.895E-01);
	}
	else if(InterpolatedPerInitial <= 2.0)
	{//1.9	<x/l_0=<	2	-9.947E-02	2	9.000E-02

		a = static_cast<Scalar>(-9.947E-02);
		b = static_cast<Scalar>(2);
		c = static_cast<Scalar>(9.000E-02);
	}
	else if(InterpolatedPerInitial <= 2.1)
	{//2	<x/l_0=<	2.1	-9.000E-02	2.1	0.000E+00

		// 誤差を考慮して2.15くらいまで
		a =static_cast<Scalar>( -9.000E-02);
		b = static_cast<Scalar>(2.1);
		c = static_cast<Scalar>(0.0000000E+00);
	}
	else
	{
		return 0.0;
	}
	Scalar Val =  a * (InterpolatedPerInitial - b) + c;
	return Val;
}


// Z関数の取得(ラプラシアン計算時の壁寄与分の密度）
static __inline__ __device__ __host__ Scalar GetZValueFunctionLaplacian(const Scalar InterpolatedDistance, const Scalar InitialParticleDistance)
{
	const Scalar InterpolatedPerInitial = InterpolatedDistance / InitialParticleDistance;
	Scalar a = 0.0, b = 0.0, c = 0.0;
	if(InterpolatedPerInitial <= 0.1)
	{
		a =static_cast<Scalar>( -9.6379E+03);
		b = static_cast<Scalar>(0.1 );
		c = static_cast<Scalar>(1.1497E+02);
	}
	else if(InterpolatedPerInitial <= 0.2)
	{
		a =static_cast<Scalar>(-2.3979E+02);
		b= static_cast<Scalar>(0.2 );
		c = static_cast<Scalar>(9.0990E+01);
	}
	else if(InterpolatedPerInitial <= 0.3)
	{
		a = static_cast<Scalar>(-1.0777E+02);
		b = static_cast<Scalar>(0.3 );
		c = static_cast<Scalar>(8.0213E+01);
	}
	else if(InterpolatedPerInitial <= 0.4)
	{
		a =static_cast<Scalar>(-7.4088E+01); 
		b= static_cast<Scalar>(0.4 ); 
		c =	static_cast<Scalar>(7.2804E+01);

	}
	else if(InterpolatedPerInitial <= 0.5)
	{
		a =static_cast<Scalar>( -6.0259E+01);
		b =static_cast<Scalar>( 0.5 );
		c = static_cast<Scalar>(6.6778E+01);
	}
	else if(InterpolatedPerInitial <= 0.6)
	{
		a = static_cast<Scalar>(-5.2803E+01);
		b =static_cast<Scalar>( 0.6 );
		c = static_cast<Scalar>(6.1498E+01);

	}
	else if(InterpolatedPerInitial <= 0.7)
	{
		a = static_cast<Scalar>(-4.8681E+01);
		b =static_cast<Scalar>( 0.700 );
		c =static_cast<Scalar>( 5.6630E+01);
	}
	else if(InterpolatedPerInitial <= 0.8)
	{
		a =static_cast<Scalar>( -4.4630E+01);
		b =static_cast<Scalar>( 0.800 );
		c =static_cast<Scalar>( 5.2167E+01);

	}
	else if(InterpolatedPerInitial <= 0.9)
	{
		a = static_cast<Scalar>(-4.1171E+01);
		b = static_cast<Scalar>(0.900 );
		c = static_cast<Scalar>(4.8050E+01);
	}
	else if(InterpolatedPerInitial <= 1.0)
	{
		a = static_cast<Scalar>(-3.8583E+01);
		b = static_cast<Scalar>(1);
		c = static_cast<Scalar>(4.4191E+01);

	}
	else if(InterpolatedPerInitial <= 1.1)
	{
		a = static_cast<Scalar>(-3.6846E+01);
		b =static_cast<Scalar>( 1.100 );
		c = static_cast<Scalar>(4.0507E+01);
	}
	else if(InterpolatedPerInitial <= 1.2)
	{
		a =static_cast<Scalar>( -3.5476E+01); 
		b = static_cast<Scalar>(	1.2);
		c = static_cast<Scalar>(3.6959E+01);

	}
	else if(InterpolatedPerInitial <= 1.3)
	{
		a = static_cast<Scalar>(-3.4183E+01);
		b = static_cast<Scalar>(1.300 );
		c = static_cast<Scalar>(3.3541E+01);
	}
	else if(InterpolatedPerInitial <= 1.4)
	{
		a = static_cast<Scalar>(-3.1589E+01);
		b = static_cast<Scalar>(1.4);
		c =static_cast<Scalar>( 3.0382E+01);

	}
	else if(InterpolatedPerInitial <= 1.5)
	{
		a = static_cast<Scalar>(-2.9258E+01);
		b = static_cast<Scalar>(1.500 );
		c = static_cast<Scalar>(2.7456E+01);
	}
	else if(InterpolatedPerInitial <= 1.6)
	{
		a =static_cast<Scalar>( -2.7047E+01);
		b = static_cast<Scalar>(1.6);
		c = static_cast<Scalar>(2.4751E+01);
	}
	else if(InterpolatedPerInitial <= 1.7)
	{
		a = static_cast<Scalar>(-2.5690E+01);
		b = static_cast<Scalar>(1.700 );
		c = static_cast<Scalar>(2.2182E+01);
	}
	else if(InterpolatedPerInitial <= 1.8)
	{
		a = static_cast<Scalar>(-2.3305E+01);
		b = static_cast<Scalar>(1.8);
		c = static_cast<Scalar>(1.9852E+01);

	}
	else if(InterpolatedPerInitial <= 1.9)
	{
		a = static_cast<Scalar>(-2.0982E+01);
		b = static_cast<Scalar>(1.900 );
		c = static_cast<Scalar>(1.7754E+01);
	}
	else if(InterpolatedPerInitial <= 2.0)
	{
		a = static_cast<Scalar>(-1.9268E+01);
		b = static_cast<Scalar>(2);
		c = static_cast<Scalar>(1.5827E+01);

	}	
	else if(InterpolatedPerInitial <= 2.4)
	{
		a = static_cast<Scalar>(-1.6887E+01);
		b = static_cast<Scalar>(2.4);
		c = static_cast<Scalar>(9.0719E+00);

	}	
	else if(InterpolatedPerInitial <= 2.8)
	{
		a = static_cast<Scalar>(-1.1091E+01); 
		b = static_cast<Scalar>(2.8);
		c = static_cast<Scalar>(4.6355E+00);

	}
	else if(InterpolatedPerInitial <= 3.2)
	{
		a =static_cast<Scalar>( -6.5068E+00);
		b =static_cast<Scalar>( 3.2);
		c = static_cast<Scalar>(2.0327E+00);

	}	
	else if(InterpolatedPerInitial <= 3.6)
	{
		a = static_cast<Scalar>(-3.7565E+00);
		b = static_cast<Scalar>(3.6);
		c = static_cast<Scalar>(5.3011E-01);

	}	
	else if(InterpolatedPerInitial <= 4.0)
	{
		a = static_cast<Scalar>(-1.3253E+00);
		b = static_cast<Scalar>(4);
		c = static_cast<Scalar>(0.0000E+00);

	}	
	else
	{
		return 0.0;
	}
	return a * (InterpolatedPerInitial - b) + c;
}

static __inline__ __device__ __host__ Scalar GetZValueFunctionLaplacian2D(const Scalar InterpolatedDistance, const Scalar InitialParticleDistance)
{
	const Scalar InterpolatedPerInitial = InterpolatedDistance / InitialParticleDistance;
	Scalar a = 0.0, b = 0.0, c = 0.0;
	if(InterpolatedPerInitial <= 0.01)
	{//0	<x/l_0=<	0.01	-9.992E+05	0.01	7.507E+02
		a = static_cast<Scalar>(-9.992E+05);
		b =static_cast<Scalar>( 0.01);
		c = static_cast<Scalar>(7.507E+02);
	}
	else if(InterpolatedPerInitial <= 0.05)
	{//0.01	<x/l_0=<	0.05	-5.766E+02	0.05	1.741E+02
		a = static_cast<Scalar>(-5.766E+02);
		b =static_cast<Scalar>( 0.05 );
		c = static_cast<Scalar>(1.741E+02);
	}
	else if(InterpolatedPerInitial <= 0.1)
	{//0.05	<x/l_0=<	0.1	-7.280E+01	0.1	1.013E+02
		a = static_cast<Scalar>(-7.280E+01);
		b = static_cast<Scalar>(0.1 );
		c = static_cast<Scalar>(1.013E+02);
	}
	else if(InterpolatedPerInitial <= 0.2)
	{//0.1	<x/l_0=<	0.2	-3.761E+01	0.2	6.368E+01
		a =static_cast<Scalar>(-3.761E+01);
		b= static_cast<Scalar>(0.2 );
		c = static_cast<Scalar>(6.368E+01);
	}	
	else if(InterpolatedPerInitial <= 0.4)
	{//0.2	<x/l_0=<	0.4	-2.128E+01	0.4	4.240E+01
		a =static_cast<Scalar>(-2.128E+01); 
		b=static_cast<Scalar>( 0.4 ); 
		c =static_cast<Scalar>(	4.240E+01);
	}	
	else if(InterpolatedPerInitial <= 0.6)
	{//0.4	<x/l_0=<	0.6	-9.280E+00	0.6	3.312E+01
		a =static_cast<Scalar>(-9.280E+00); 
		b=static_cast<Scalar>( 0.6 ); 
		c =static_cast<Scalar>(	3.312E+01);
	}	
	else if(InterpolatedPerInitial <= 0.8)
	{//0.6	<x/l_0=<	0.8	-6.104E+00	0.8	2.701E+01
		a =static_cast<Scalar>( -6.104E+00);
		b = static_cast<Scalar>(0.8 );
		c =static_cast<Scalar>( 2.701E+01);
	}	
	else if(InterpolatedPerInitial <= 1.0)
	{//0.8	<x/l_0=<	1	-4.682E+00	1	2.233E+01
		a = static_cast<Scalar>(-4.682E+00);
		b = static_cast<Scalar>(1);
		c = static_cast<Scalar>(2.233E+01);
	}	
	else if(InterpolatedPerInitial <= 1.2)
	{//1	<x/l_0=<	1.2	-3.842E+00	1.2	1.849E+01
		a =static_cast<Scalar>( -3.842E+00); 
		b = static_cast<Scalar>(	1.2);
		c =static_cast<Scalar>( 1.849E+01);
	}	
	else if(InterpolatedPerInitial <= 1.4)
	{//1.2	<x/l_0=<	1.4	-3.260E+00	1.4	1.523E+01
		a =static_cast<Scalar>( -3.260E+00);
		b =static_cast<Scalar>( 1.4);
		c =static_cast<Scalar>( 1.523E+01);
	}	
	else if(InterpolatedPerInitial <= 1.6)
	{//1.4	<x/l_0=<	1.6	-2.763E+00	1.6	1.247E+01
		a = static_cast<Scalar>(-2.763E+00);
		b =static_cast<Scalar>( 1.6);
		c = static_cast<Scalar>(1.247E+01);
	}	
	else if(InterpolatedPerInitial <= 1.8)
	{//1.6	<x/l_0=<	1.8	-2.048E+00	1.8	1.042E+01
		a = static_cast<Scalar>(-2.048E+00);
		b = static_cast<Scalar>(1.8);
		c = static_cast<Scalar>(1.042E+01);
	}	
	else if(InterpolatedPerInitial <= 2.0)
	{//1.8	<x/l_0=<	2	-2.034E+00	2	8.385E+00
		a =static_cast<Scalar>( -2.034E+00);
		b = static_cast<Scalar>(2);
		c = static_cast<Scalar>(8.385E+00);
	}	
	else if(InterpolatedPerInitial <= 2.2)
	{//2	<x/l_0=<	2.2	-1.726E+00	2.2	6.659E+00
		a = static_cast<Scalar>(-1.726E+00);
		b = static_cast<Scalar>(2.2);
		c = static_cast<Scalar>(6.659E+00);
	}	
	else if(InterpolatedPerInitial <= 2.4)
	{//2.2	<x/l_0=<	2.4	-1.543E+00	2.4	5.116E+00
		a = static_cast<Scalar>(-1.543E+00);
		b =static_cast<Scalar>( 2.4);
		c = static_cast<Scalar>(5.116E+00);
	}	
	else if(InterpolatedPerInitial <= 2.6)
	{//2.4	<x/l_0=<	2.6	-1.306E+00	2.6	3.810E+00
		a = static_cast<Scalar>(-1.306E+00);
		b = static_cast<Scalar>(2.6);
		c = static_cast<Scalar>(3.810E+00);
	}	
	else if(InterpolatedPerInitial <= 2.8)
	{//2.6	<x/l_0=<	2.8	-1.051E+00	2.8	2.759E+00
		a = static_cast<Scalar>(-1.051E+00); 
		b = static_cast<Scalar>(2.8);
		c = static_cast<Scalar>(2.759E+00);
	}
	else if(InterpolatedPerInitial <= 3.0)
	{//2.8	<x/l_0=<	3	-8.115E-01	3	1.948E+00
		a =static_cast<Scalar>( -8.115E-01); 
		b =static_cast<Scalar>( 3);
		c = static_cast<Scalar>(1.948E+00);
	}
	else if(InterpolatedPerInitial <= 3.2)
	{//3	<x/l_0=<	3.2	-5.864E-01	3.2	1.361E+00
		a =  static_cast<Scalar>(-5.864E-01);
		b =  static_cast<Scalar>(3.2);
		c =  static_cast<Scalar>(1.361E+00);
	}	
	else if(InterpolatedPerInitial <= 3.4)
	{//3.2	<x/l_0=<	3.4	-5.298E-01	3.4	8.314E-01
		a = static_cast<Scalar>(-5.298E-01);
		b = static_cast<Scalar>(3.4);
		c = static_cast<Scalar>(8.314E-01);
	}	
	else if(InterpolatedPerInitial <= 3.6)
	{//3.4	<x/l_0=<	3.6	-3.773E-01	3.6	4.541E-01
		a =static_cast<Scalar>( -3.773E-01);
		b =static_cast<Scalar>( 3.6);
		c = static_cast<Scalar>(4.541E-01);
	}	
	else if(InterpolatedPerInitial <= 3.8)
	{//3.6	<x/l_0=<	3.8	-2.946E-01	3.8	1.594E-01
		a = static_cast<Scalar>(-2.946E-01);
		b = static_cast<Scalar>(3.8);
		c = static_cast<Scalar>(1.594E-01);
	}	
	else if(InterpolatedPerInitial <= 4.0)
	{//3.8	<x/l_0=<	4	-1.594E-01	4	0.000E+00
		a =static_cast<Scalar>( -1.594E-01);
		b =static_cast<Scalar>( 4);
		c =static_cast<Scalar>( 0.0000E+00);
	}	
	else
	{
		return 0.0;
	}
	return a * (InterpolatedPerInitial - b) + c;
}

static __inline__ __device__ __host__ Scalar GetZValueFunctionSurface(const Scalar InterpolatedDistance, const Scalar InitialParticleDistance)
{
	const Scalar InterpolatedPerInitial = InterpolatedDistance / InitialParticleDistance;
	Scalar a = 0.0, b = 0.0, c = 0.0;
	if(InterpolatedPerInitial <= 0.25)
	{//0	<r/l_0=<	0.25	-4.073872	0.25	15.144314
		a = static_cast<Scalar>(-4.073872);
		b = static_cast<Scalar>(0.25);
		c = static_cast<Scalar>(15.144314);
	}
	else if(InterpolatedPerInitial <= 0.5)
	{//0.25	<r/l_0=<	0.5	-8.203804	0.5	13.093363
		a = static_cast<Scalar>(-8.203804);
		b = static_cast<Scalar>(0.5);
		c = static_cast<Scalar>(13.093363);
	}
	else if(InterpolatedPerInitial <= 0.75)
	{//0.5	<r/l_0=<	0.75	-8.047464	0.75	11.081497
		a = static_cast<Scalar>(-8.047464);
		b = static_cast<Scalar>(0.75);
		c = static_cast<Scalar>(11.081497);
	}
	else if(InterpolatedPerInitial <= 1.0)
	{//0.75	<r/l_0=<	1	-7.628876	1	9.174278
		a = static_cast<Scalar>(-7.628876);
		b = static_cast<Scalar>(1.0);
		c = static_cast<Scalar>(9.174278);
	}
	else if(InterpolatedPerInitial <= 1.25)
	{//1	<r/l_0=<	1.25	-7.001212	1.25	7.423975

		a = static_cast<Scalar>(-7.001212);
		b = static_cast<Scalar>(1.25);
		c = static_cast<Scalar>(7.423975);
	}
	else if(InterpolatedPerInitial <= 1.5)
	{//1.25	<r/l_0=<	1.5	-6.230716	1.5	5.866296

		a = static_cast<Scalar>(-6.230716);
		b = static_cast<Scalar>(1.5);
		c = static_cast<Scalar>(5.866296);
	}
	else if(InterpolatedPerInitial <= 1.75)
	{//1.5	<r/l_0=<	1.75	-5.392776	1.75	4.518102

		a = static_cast<Scalar>(-5.392776);
		b = static_cast<Scalar>(1.75);
		c = static_cast<Scalar>(4.518102);
	}
	else if(InterpolatedPerInitial <= 2.0)
	{//1.75	<r/l_0=<	2	-4.541608	2	3.3827

		a = static_cast<Scalar>(-4.541608);
		b = static_cast<Scalar>(2.0);
		c = static_cast<Scalar>(3.3827);
	}
	else if(InterpolatedPerInitial <= 2.25)
	{//2	<r/l_0=<	2.25	-3.715672	2.25	2.453782

		a = static_cast<Scalar>(-3.715672);
		b = static_cast<Scalar>(2.25);
		c = static_cast<Scalar>(2.453782);
	}
	else if(InterpolatedPerInitial <= 2.5)
	{//2.25	<r/l_0=<	2.5	-2.94322	2.5	1.717977

		a = static_cast<Scalar>(-2.94322);
		b = static_cast<Scalar>(2.5);
		c = static_cast<Scalar>(1.717977);
	}
	else if(InterpolatedPerInitial <= 2.75)
	{//2.5	<r/l_0=<	2.75	-2.256608	2.75	1.153825

		a = static_cast<Scalar>(-2.256608);
		b = static_cast<Scalar>(2.75);
		c = static_cast<Scalar>(1.153825);
	}
	else if(InterpolatedPerInitial <= 3.0)
	{//2.75	<r/l_0=<	3	-1.66868	3	0.736655

		a = static_cast<Scalar>(-1.66868);
		b = static_cast<Scalar>(3.0);
		c = static_cast<Scalar>(0.736655);
	}
	else if(InterpolatedPerInitial <= 3.25)
	{//3	<r/l_0=<	3.25	-1.180056	3.25	0.441641

		a = static_cast<Scalar>(-1.180056);
		b = static_cast<Scalar>(3.25);
		c = static_cast<Scalar>(0.441641);
	}
	else if(InterpolatedPerInitial <= 3.5)
	{//3.25	<r/l_0=<	3.5	-0.785508	3.5	0.245264

		a = static_cast<Scalar>(-0.785508);
		b = static_cast<Scalar>(3.5);
		c = static_cast<Scalar>(0.245264);
	}
	else if(InterpolatedPerInitial <= 3.75)
	{//3.5	<r/l_0=<	3.75	-0.489308	3.75	0.122937

		a = static_cast<Scalar>(-0.489308);
		b = static_cast<Scalar>(3.75);
		c = static_cast<Scalar>(0.122937);
	}
	else if(InterpolatedPerInitial <= 4.0)
	{//3.75	<r/l_0=<	4	-0.280256	4	0.052873

		a = static_cast<Scalar>(-0.280256);
		b = static_cast<Scalar>(4.0);
		c = static_cast<Scalar>(0.052873);
	}
	else if(InterpolatedPerInitial <= 4.25)
	{//4	<r/l_0=<	4.25	-0.141756	4.25	0.017434

		a = static_cast<Scalar>(-0.141756);
		b = static_cast<Scalar>(4.25);
		c = static_cast<Scalar>(0.017434);
	}
	else if(InterpolatedPerInitial <= 4.5)
	{//4.25	<r/l_0=<	4.5	-0.055172	4.5	0.003641

		a = static_cast<Scalar>(-0.055172);
		b = static_cast<Scalar>(4.5);
		c = static_cast<Scalar>(0.003641);
	}
	else if(InterpolatedPerInitial <= 4.75)
	{//4.5	<r/l_0=<	4.75	-0.013656	4.75	0.000227

		a = static_cast<Scalar>(-0.013656);
		b = static_cast<Scalar>(4.75);
		c = static_cast<Scalar>(0.000227);
	}
	else if(InterpolatedPerInitial <= 5.0)
	{//4.75	<r/l_0=<	5	-0.000908	5	0

		a = static_cast<Scalar>(-0.000908);
		b = static_cast<Scalar>(5.0);
		c = static_cast<Scalar>(0.0);
	}	
	else
	{
		return 0.0;
	}
	Scalar Val =  a * (InterpolatedPerInitial - b) + c;
	return Val;
}

static __inline__ __device__ __host__ __device__ Scalar GetDummyWeightSquare(const Scalar Re, const Scalar Distance)
{
	Scalar q = (Distance/Re);
	if(Distance > Re)
	{
		return 0;
	}
	
	if(q < 1.0)
	{		
		return  (1-q)*(1-q);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ;
	}

	if(0 == Distance)
	{
		return 1e6;
	}
	
	return -1;
}
static __inline__ __device__ __host__ __device__ Scalar GetNearWeightCubic(const Scalar Re, const Scalar Distance)
{
	Scalar q = (Distance/Re);
	if(Distance > Re)
	{
		return 0;
	}
	
	if(q < 1.0)
	{		
		return  (1-q)*(1-q)*(1-q);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ;
	}

	if(0 == Distance)
	{
		return 1e6;
	}
	
	return -1;
}
static __inline__ __device__ __host__ Scalar GetZValueFunctionGradientDummyDensity(const Scalar InterpolatedDistance, const Scalar InitialParticleDistance)
{
	const Scalar InterpolatedPerInitial = InterpolatedDistance / InitialParticleDistance;
	Scalar a = 0.0, b = 0.0, c = 0.0;
	if(InterpolatedPerInitial <= 0.05)
	{//-1.85978	0.05	1.902776

		a = static_cast<Scalar>(-1.85978);
		b =static_cast<Scalar>( 0.05);
		c = static_cast<Scalar>(1.902776);
	}
	else if(InterpolatedPerInitial <= 0.1)
	{//-2.27896	0.1	1.788828

		a = static_cast<Scalar>(-2.27896);
		b = static_cast<Scalar>(0.1);
		c = static_cast<Scalar>(1.788828);
	}
	else if(InterpolatedPerInitial <= 0.15)
	{//-2.22854	0.15	1.677401

		a = static_cast<Scalar>(-2.22854);
		b = static_cast<Scalar>(0.15);
		c = static_cast<Scalar>(1.677401);
	}
	else if(InterpolatedPerInitial <= 0.2)
	{//-2.17656	0.2	1.568573

		a = static_cast<Scalar>(-2.17656);
		b = static_cast<Scalar>(0.2);
		c = static_cast<Scalar>(1.568573);
	}
	else if(InterpolatedPerInitial <= 0.25)
	{//-2.1219	0.25	1.462478

		a = static_cast<Scalar>(-2.1219);
		b = static_cast<Scalar>(0.25);
		c = static_cast<Scalar>(1.462478);
	}
	else if(InterpolatedPerInitial <= 0.3)
	{//-2.06388	0.3	1.359284

		a = static_cast<Scalar>(-2.06388);
		b = static_cast<Scalar>(0.3);
		c = static_cast<Scalar>(1.359284);
	}
	else if(InterpolatedPerInitial <= 0.35)
	{//-2.00216	0.35	1.259176

		a = static_cast<Scalar>(-2.00216);
		b = static_cast<Scalar>(0.35);
		c = static_cast<Scalar>(1.259176);
	}
	else if(InterpolatedPerInitial <= 0.4)
	{//-1.93652	0.4	1.16235

		a = static_cast<Scalar>(-1.93652);
		b = static_cast<Scalar>(0.4);
		c = static_cast<Scalar>(1.16235);
	}
	else if(InterpolatedPerInitial <= 0.45)
	{//-1.8672	0.45	1.06899

		a = static_cast<Scalar>(-1.8672);
		b = static_cast<Scalar>(0.45);
		c = static_cast<Scalar>(1.06899);
	}
	else if(InterpolatedPerInitial <= 0.5)
	{//-1.79454	0.5	0.979263

		a = static_cast<Scalar>(-1.79454);
		b = static_cast<Scalar>(0.5);
		c = static_cast<Scalar>(0.979263);
	}
	else if(InterpolatedPerInitial <= 0.55)
	{//-1.71922	0.55	0.893302

		a = static_cast<Scalar>(-1.71922);
		b = static_cast<Scalar>(0.55);
		c = static_cast<Scalar>(0.893302);
	}
	else if(InterpolatedPerInitial <= 0.6)
	{//-1.64176	0.6	0.811214

		a = static_cast<Scalar>(-1.64176);
		b = static_cast<Scalar>(0.6);
		c = static_cast<Scalar>(0.811214);
	}
	else if(InterpolatedPerInitial <= 0.65)
	{//-1.5613	0.65	0.733149

		a = static_cast<Scalar>(-1.5613);
		b = static_cast<Scalar>(0.65);
		c = static_cast<Scalar>(0.733149);
	}
	else if(InterpolatedPerInitial <= 0.7)
	{//-1.4783	0.7	0.659234

		a = static_cast<Scalar>(-1.4783);
		b =static_cast<Scalar>( 0.7);
		c = static_cast<Scalar>(0.659234);
	}
	else if(InterpolatedPerInitial <= 0.75)
	{//-1.39368	0.75	0.58955

		a = static_cast<Scalar>(-1.39368);
		b = static_cast<Scalar>(0.75);
		c = static_cast<Scalar>(0.58955);
	}
	else if(InterpolatedPerInitial <= 0.8)
	{//-1.30852	0.8	0.524124

		a = static_cast<Scalar>(-1.30852);
		b = static_cast<Scalar>(0.8);
		c = static_cast<Scalar>(0.524124);
	}
	else if(InterpolatedPerInitial <= 0.85)
	{//-1.2239	0.85	0.462929

		a = static_cast<Scalar>(-1.2239);
		b = static_cast<Scalar>(0.85);
		c = static_cast<Scalar>(0.462929);
	}
	else if(InterpolatedPerInitial <= 0.9)
	{//-1.14014	0.9	0.405922

		a = static_cast<Scalar>(-1.14014);
		b = static_cast<Scalar>(0.9);
		c = static_cast<Scalar>(0.405922);
	}
	else if(InterpolatedPerInitial <= 0.95)
	{//-1.0554	0.95	0.353152

		a = static_cast<Scalar>(-1.0554);
		b = static_cast<Scalar>(0.95);
		c = static_cast<Scalar>(0.353152);
	}
	else if(InterpolatedPerInitial <= 1.0)
	{//-0.97002	1	0.304651

		a = static_cast<Scalar>(-0.97002);
		b = static_cast<Scalar>(1.0);
		c = static_cast<Scalar>(0.304651);
	}
	else if(InterpolatedPerInitial <= 1.05)
	{//-0.88506	1.05	0.260398

		a = static_cast<Scalar>(-0.88506);
		b = static_cast<Scalar>(1.05);
		c = static_cast<Scalar>(0.260398);
	}
	else if(InterpolatedPerInitial <= 1.1)
	{//	-0.80146	1.1	0.220325
	
		a = static_cast<Scalar>(-0.80146);
		b = static_cast<Scalar>(1.1);
		c = static_cast<Scalar>(0.220325);
	}
	else if(InterpolatedPerInitial <= 1.15)
	{//	-0.71984	1.15	0.184333
	
		a = static_cast<Scalar>(-0.71984);
		b = static_cast<Scalar>(1.15);
		c = static_cast<Scalar>(0.184333);
	}
	else if(InterpolatedPerInitial <= 1.2)
	{// -0.64004	1.2	0.152331
		
		a = static_cast<Scalar>(-0.64004);
		b = static_cast<Scalar>(1.2);
		c =static_cast<Scalar>( 0.152331);
	}
	else if(InterpolatedPerInitial <= 1.25)
	{//	-0.56268	1.25	0.124197
	
		a =static_cast<Scalar>( -0.56268);
		b = static_cast<Scalar>(1.25);
		c = static_cast<Scalar>(0.124197);
	}
	else if(InterpolatedPerInitial <= 1.3)
	{//	-0.48848	1.3	0.099773
	
		a = static_cast<Scalar>(-0.48848);
		b = static_cast<Scalar>(1.3);
		c = static_cast<Scalar>(0.099773);
	}
	else if(InterpolatedPerInitial <= 1.35)
	{//	-0.41824	1.35	0.078861
	
		a = static_cast<Scalar>(-0.41824);
		b = static_cast<Scalar>(1.35);
		c =static_cast<Scalar>(0.078861);
	}
	else if(InterpolatedPerInitial <= 1.4)
	{//	-0.3527	1.4	0.061226
	
		a =static_cast<Scalar>( -0.3527);
		b =static_cast<Scalar>( 1.4);
		c =static_cast<Scalar>( 0.061226);
	}
	else if(InterpolatedPerInitial <= 1.45)
	{// -0.2925	1.45	0.046601
		
		a = static_cast<Scalar>(-0.2925);
		b = static_cast<Scalar>(1.45);
		c = static_cast<Scalar>(0.046601);
	}
	else if(InterpolatedPerInitial <= 1.5)
	{// -0.23824	1.5	0.034689
		
		a = static_cast<Scalar>(-0.23824);
		b = static_cast<Scalar>(1.5);
		c = static_cast<Scalar>(0.034689);
	}
	else if(InterpolatedPerInitial <= 1.55)
	{//	-0.19054	1.55	0.025162
	
		a = static_cast<Scalar>(-0.19054);
		b = static_cast<Scalar>(1.55);
		c = static_cast<Scalar>(0.025162);
	}
	else if(InterpolatedPerInitial <= 1.6)
	{//	-0.14958	1.6	0.017683
	
		a =static_cast<Scalar>( -0.14958);
		b = static_cast<Scalar>(1.6);
		c = static_cast<Scalar>(0.017683);
	}
	else if(InterpolatedPerInitial <= 1.65)
	{//	-0.11396	1.65	0.011985
	
		a = static_cast<Scalar>(-0.11396);
		b = static_cast<Scalar>(1.65);
		c = static_cast<Scalar>(0.011985);
	}
	else if(InterpolatedPerInitial <= 1.7)
	{// -0.0835	1.7	0.00781
		
		a =static_cast<Scalar>( -0.0835);
		b = static_cast<Scalar>(1.7);
		c = static_cast<Scalar>(0.00781);
	}
	else if(InterpolatedPerInitial <= 1.75)
	{//	-0.0584	1.75	0.00489
	
		a =static_cast<Scalar>( -0.0584);
		b = static_cast<Scalar>(1.75);
		c = static_cast<Scalar>(0.00489);
	}
	else if(InterpolatedPerInitial <= 1.8)
	{//	-0.0389	1.8	0.002945
	
		a = static_cast<Scalar>(-0.0389);
		b = static_cast<Scalar>(1.8);
		c = static_cast<Scalar>(0.002945);
	}
	else if(InterpolatedPerInitial <= 1.85)
	{//	-0.02516	1.85	0.001687
	
		a = static_cast<Scalar>(-0.02516);
		b = static_cast<Scalar>(1.85);
		c = static_cast<Scalar>(0.001687);
	}
	else if(InterpolatedPerInitial <= 1.9)
	{//	-0.01646	1.9	0.000864
	
		a = static_cast<Scalar>(-0.01646);
		b = static_cast<Scalar>(1.9);
		c = static_cast<Scalar>(0.000864);
	}
	else if(InterpolatedPerInitial <= 1.95)
	{// -0.01	1.95	0.000364
		
		a = static_cast<Scalar>(-0.01);
		b = static_cast<Scalar>(1.95);
		c = static_cast<Scalar>(0.000364);
	}
	else if(InterpolatedPerInitial <= 2.0)
	{//	-0.00512	2	0.000108
	
		a = static_cast<Scalar>(-0.00512);
		b = static_cast<Scalar>(2.0);
		c = static_cast<Scalar>(0.000108);
	}
	else if(InterpolatedPerInitial <= 2.05)
	{// -0.0019	2.05	0.000013
		
		a = static_cast<Scalar>(-0.0019);
		b = static_cast<Scalar>(2.05);
		c = static_cast<Scalar>(0.000013);
	}
	else if(InterpolatedPerInitial <= 2.1)
	{// -0.00026	2.1	0
		
		a = static_cast<Scalar>(-0.00026);
		b = static_cast<Scalar>(2.1);
		c = static_cast<Scalar>(0.0);
	}
	else
	{
		return 0.0;
	}
	Scalar Val =  a * (InterpolatedPerInitial - b) + c;
	return Val;
}
static __inline__ __device__ __host__ Scalar GetZValueFunctionGradientNearDensity(const Scalar InterpolatedDistance, const Scalar InitialParticleDistance)
{
	const Scalar InterpolatedPerInitial = InterpolatedDistance / InitialParticleDistance;
	Scalar a = 0.0, b = 0.0, c = 0.0;
	if(InterpolatedPerInitial <= 0.05)
	{//-2.22678	0.05	3.220183

		a = static_cast<Scalar>(-2.22678);
		b = static_cast<Scalar>(0.05);
		c = static_cast<Scalar>(3.220183);
	}
	else if(InterpolatedPerInitial <= 0.1)
	{//-2.80976	0.1	3.079695

		a = static_cast<Scalar>(-2.80976);
		b = static_cast<Scalar>(0.1);
		c = static_cast<Scalar>(3.079695);
	}
	else if(InterpolatedPerInitial <= 0.15)
	{//-2.83986	0.15	2.937702

		a = static_cast<Scalar>(-2.83986);
		b = static_cast<Scalar>(0.15);
		c = static_cast<Scalar>(2.937702);
	}
	else if(InterpolatedPerInitial <= 0.2)
	{//-2.8687	0.2	2.794267

		a = static_cast<Scalar>(-2.8687);
		b = static_cast<Scalar>(0.2);
		c = static_cast<Scalar>(2.794267);
	}
	else if(InterpolatedPerInitial <= 0.25)
	{// -2.8836	0.25	2.650087

		a = static_cast<Scalar>(-2.8836);
		b = static_cast<Scalar>(0.25);
		c = static_cast<Scalar>(2.650087);
	}
	else if(InterpolatedPerInitial <= 0.3)
	{// -2.8833	0.3	2.505922

		a = static_cast<Scalar>(-2.8833);
		b = static_cast<Scalar>(0.3);
		c = static_cast<Scalar>(2.505922);
	}
	else if(InterpolatedPerInitial <= 0.35)
	{// -2.86688	0.35	2.362578

		a = static_cast<Scalar>(-2.86688);
		b = static_cast<Scalar>(0.35);
		c = static_cast<Scalar>(2.362578);
	}
	else if(InterpolatedPerInitial <= 0.4)
	{// -2.83348	0.4	2.220904

		a = static_cast<Scalar>(-2.83348);
		b = static_cast<Scalar>(0.4);
		c = static_cast<Scalar>(2.220904);
	}
	else if(InterpolatedPerInitial <= 0.45)
	{// -2.78258	0.45	2.081775

		a = static_cast<Scalar>(-2.78258);
		b = static_cast<Scalar>(0.45);
		c = static_cast<Scalar>(2.081775);
	}
	else if(InterpolatedPerInitial <= 0.5)
	{// -2.71386	0.5	1.946082

		a = static_cast<Scalar>(-2.71386);
		b = static_cast<Scalar>(0.5);
		c = static_cast<Scalar>(1.946082);
	}
	else if(InterpolatedPerInitial <= 0.55)
	{// -2.62708	0.55	1.814728

		a = static_cast<Scalar>(-2.62708);
		b = static_cast<Scalar>(0.55);
		c = static_cast<Scalar>(1.814728);
	}
	else if(InterpolatedPerInitial <= 0.6)
	{// -2.54514	0.6	1.687471

		a = static_cast<Scalar>(-2.54514);
		b = static_cast<Scalar>(0.6);
		c = static_cast<Scalar>(1.687471);
	}
	else if(InterpolatedPerInitial <= 0.65)
	{// -2.47432	0.65	1.563755

		a = static_cast<Scalar>(-2.47432);
		b = static_cast<Scalar>(0.65);
		c = static_cast<Scalar>(1.563755);
	}
	else if(InterpolatedPerInitial <= 0.7)
	{// -2.39412	0.7	1.444049

		a = static_cast<Scalar>(-2.39412);
		b = static_cast<Scalar>(0.7);
		c = static_cast<Scalar>(1.444049);
	}
	else if(InterpolatedPerInitial <= 0.75)
	{// -2.30322	0.75	1.328888

		a = static_cast<Scalar>(-2.30322);
		b = static_cast<Scalar>(0.75);
		c = static_cast<Scalar>(1.328888);
	}
	else if(InterpolatedPerInitial <= 0.8)
	{// -2.19946	0.8	1.218915

		a = static_cast<Scalar>(-2.19946);
		b = static_cast<Scalar>(0.8);
		c = static_cast<Scalar>(1.218915);
	}
	else if(InterpolatedPerInitial <= 0.85)
	{// -2.08352	0.85	1.114739

		a = static_cast<Scalar>(-2.08352);
		b = static_cast<Scalar>(0.85);
		c = static_cast<Scalar>(1.114739);
	}
	else if(InterpolatedPerInitial <= 0.9)
	{// -1.99556	0.9	1.014961

		a = static_cast<Scalar>(-1.99556);
		b = static_cast<Scalar>(0.9);
		c = static_cast<Scalar>(1.014961);
	}
	else if(InterpolatedPerInitial <= 0.95)
	{// -1.92792	0.95	0.918565

		a = static_cast<Scalar>(-1.92792);
		b = static_cast<Scalar>(0.95);
		c = static_cast<Scalar>(0.918565);
	}
	else if(InterpolatedPerInitial <= 1.0)
	{// -1.85086	1	0.826022

		a = static_cast<Scalar>(-1.85086);
		b = static_cast<Scalar>(1.0);
		c = static_cast<Scalar>(0.826022);
	}
	else if(InterpolatedPerInitial <= 1.05)
	{// -1.76478	1.05	0.737783

		a = static_cast<Scalar>(-1.76478);
		b = static_cast<Scalar>(1.05);
		c = static_cast<Scalar>(0.737783);
	}
	else if(InterpolatedPerInitial <= 1.1)
	{//	-1.67026	1.1	0.65427
	
		a = static_cast<Scalar>(-1.67026);
		b = static_cast<Scalar>(1.1);
		c = static_cast<Scalar>(0.65427);
	}
	else if(InterpolatedPerInitial <= 1.15)
	{// -1.5791	1.15	0.575315
		
		a = static_cast<Scalar>(-1.5791);
		b = static_cast<Scalar>(1.15);
		c = static_cast<Scalar>(0.575315);
	}
	else if(InterpolatedPerInitial <= 1.2)
	{//	-1.49178	1.2	0.500726
	
		a = static_cast<Scalar>(-1.49178);
		b = static_cast<Scalar>(1.2);
		c = static_cast<Scalar>(0.500726);
	}
	else if(InterpolatedPerInitial <= 1.25)
	{// -1.3974	1.25	0.430856
		
		a = static_cast<Scalar>(-1.3974);
		b = static_cast<Scalar>(1.25);
		c = static_cast<Scalar>(0.430856);
	}
	else if(InterpolatedPerInitial <= 1.3)
	{// -1.29646	1.3	0.366033
		
		a = static_cast<Scalar>(-1.29646);
		b = static_cast<Scalar>(1.3);
		c = static_cast<Scalar>(0.366033);
	}
	else if(InterpolatedPerInitial <= 1.35)
	{//	-1.18932	1.35	0.306567
	
		a =static_cast<Scalar>( -1.18932);
		b = static_cast<Scalar>(1.35);
		c = static_cast<Scalar>(0.306567);
	}
	else if(InterpolatedPerInitial <= 1.4)
	{// -1.07642	1.4	0.252746
		
		a = static_cast<Scalar>(-1.07642);
		b = static_cast<Scalar>(1.4);
		c = static_cast<Scalar>(0.252746);
	}
	else if(InterpolatedPerInitial <= 1.45)
	{// -0.95812	1.45	0.20484
		
		a = static_cast<Scalar>(-0.95812);
		b = static_cast<Scalar>(1.45);
		c = static_cast<Scalar>(0.20484);
	}
	else if(InterpolatedPerInitial <= 1.5)
	{// -0.83474	1.5	0.163103
		
		a = static_cast<Scalar>(-0.83474);
		b = static_cast<Scalar>(1.5);
		c = static_cast<Scalar>(0.163103);
	}
	else if(InterpolatedPerInitial <= 1.55)
	{// -0.7067	1.55	0.127768
		
		a = static_cast<Scalar>(-0.7067);
		b = static_cast<Scalar>(1.55);
		c = static_cast<Scalar>(0.127768);
	}
	else if(InterpolatedPerInitial <= 1.6)
	{// -0.59698	1.6	0.097919

		a = static_cast<Scalar>( -0.59698);
		b = static_cast<Scalar>(1.6);
		c = static_cast<Scalar>(0.097919);
	}
	else if(InterpolatedPerInitial <= 1.65)
	{// -0.5119	1.65	0.072324
		
		a = static_cast<Scalar>(-0.5119);
		b = static_cast<Scalar>(1.65);
		c = static_cast<Scalar>(0.072324);
	}
	else if(InterpolatedPerInitial <= 1.7)
	{// -0.42508	1.7	0.05107
		
		a = static_cast<Scalar>(-0.42508);
		b = static_cast<Scalar>(1.7);
		c = static_cast<Scalar>(0.05107);
	}
	else if(InterpolatedPerInitial <= 1.75)
	{// -0.33652	1.75	0.034244
		
		a = static_cast<Scalar>(-0.33652);
		b = static_cast<Scalar>(1.75);
		c = static_cast<Scalar>(0.034244);
	}
	else if(InterpolatedPerInitial <= 1.8)
	{// -0.2464	1.8	0.021924
		
		a = static_cast<Scalar>(-0.2464);
		b = static_cast<Scalar>(1.8);
		c = static_cast<Scalar>(0.021924);
	}
	else if(InterpolatedPerInitial <= 1.85)
	{// -0.15504	1.85	0.014172
		
		a = static_cast<Scalar>(-0.15504);
		b = static_cast<Scalar>(1.85);
		c = static_cast<Scalar>(0.014172);
	}
	else if(InterpolatedPerInitial <= 1.9)
	{// -0.10204	1.9	0.00907
		
		a = static_cast<Scalar>(-0.10204);
		b = static_cast<Scalar>(1.9);
		c = static_cast<Scalar>(0.00907);
	}
	else if(InterpolatedPerInitial <= 1.95)
	{// -0.07936	1.95	0.005102
		
		a = static_cast<Scalar>(-0.07936);
		b = static_cast<Scalar>(1.95);
		c = static_cast<Scalar>(0.005102);
	}
	else if(InterpolatedPerInitial <= 2.0)
	{// -0.05668	2	0.002268
		
		a = static_cast<Scalar>(-0.05668);
		b = static_cast<Scalar>(2.0);
		c = static_cast<Scalar>(0.002268);
	}
	else if(InterpolatedPerInitial <= 2.05)
	{// -0.03402	2.05	0.000567
		
		a = static_cast<Scalar>(-0.03402);
		b = static_cast<Scalar>(2.05);
		c = static_cast<Scalar>(0.000567);
	}
	else if(InterpolatedPerInitial <= 2.1)
	{// -0.01134	2.1	0
		
		a = static_cast<Scalar>(-0.01134);
		b = static_cast<Scalar>(2.1);
		c = static_cast<Scalar>(0.0);
	}
	else
	{
		return 0.0;
	}
	Scalar Val =  a * (InterpolatedPerInitial - b) + c;
	return Val;
}

static __inline__ __device__ __host__ Integer GetCellID(const CGridBox* const BoundingBox, const Scalar3* const Position,Integer& i, Integer& j, Integer &k) 
{
	Scalar RelativeX = Position->x - BoundingBox->m_BufferedZone.m_MinBound.x /*+ EPS*/;
	Scalar RelativeY = Position->y - BoundingBox->m_BufferedZone.m_MinBound.y /*+ EPS*/;
	Scalar RelativeZ = Position->z - BoundingBox->m_BufferedZone.m_MinBound.z /*+ EPS*/;

	RelativeX /= BoundingBox->m_CellSize;
	RelativeY /= BoundingBox->m_CellSize;
	RelativeZ /= BoundingBox->m_CellSize;
	i = (Integer)RelativeX;
	j = (Integer)RelativeY;
	k = (Integer)RelativeZ;
	//(Section.y * Section.z) * i + Section.z * j + k;
	const Integer CellID = (BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z) * i + BoundingBox->m_GridSize.z * j + k;
	return CellID;	
}
static __inline__ __device__ __host__ Integer GetCellID(const CGridBox* const BoundingBox, const Scalar3* const Position) 
{
	Scalar RelativeX = Position->x - BoundingBox->m_BufferedZone.m_MinBound.x /*+ EPS*/;
	Scalar RelativeY = Position->y - BoundingBox->m_BufferedZone.m_MinBound.y /*+ EPS*/;
	Scalar Relativez = Position->z - BoundingBox->m_BufferedZone.m_MinBound.z /*+ EPS*/;

	RelativeX /= BoundingBox->m_CellSize;
	RelativeY /= BoundingBox->m_CellSize;
	Relativez /= BoundingBox->m_CellSize;

	Integer i = (Integer)RelativeX;
	Integer j = (Integer)RelativeY;
	Integer k = (Integer)Relativez;
	const Integer CellID = (BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z) * i + BoundingBox->m_GridSize.z * j + k;
	return CellID;	
}
static __inline__ __device__ __host__ Integer GetCellID(const CGridBox* const BoundingBox, Integer i, Integer j, Integer k) 
{	
	const Integer CellID = (BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z) * i + BoundingBox->m_GridSize.z * j + k;
	return CellID;	
}


static __inline__ __device__ __host__ Scalar3 AddVector(const Scalar3* const P1, const Scalar3* const P2)
{
	Scalar3 tmp;
	tmp.x = P1->x + P2->x;
	tmp.y = P1->y + P2->y;
	tmp.z = P1->z + P2->z;
	return tmp;
}
#ifdef __CUDACC__ 
static __inline__ __device__ void RegisterLine(const CGridBox *const BoundingBox,Scalar3 Vertex, Scalar3 Direction, Scalar Magnitude, Integer TID, CCell* daCell, Integer CellNum)
{
	Scalar TAB = 0.0;
	Scalar3 Target;
	do
	{
		Target.x = Vertex.x + TAB * Direction.x;
		Target.y = Vertex.y + TAB * Direction.y;
		Target.z = Vertex.z + TAB * Direction.z;
		Integer CID = GetCellID(BoundingBox,&Target);
		if(CID >= 0 && CID < CellNum)
		{			
			CCell* Cell = &daCell[CID];
			// Atomic compare and swap
			// if the trianglenumber is 0, 1 is set and 0 is returned
			unsigned int old = atomicCAS(&Cell->m_TriangleNum,0,1);
			if(old <= 0)
			{
				Cell->m_TriangleID[0] = TID;

			}
			else
			{
				bool HasTID = false;
				for(Integer l = 0; l < old; ++l)
				{
					if(Cell->m_TriangleID[l] == TID)
					{
						HasTID = true;
						break;
					}
				}
				if(!HasTID && old < CELLMAXTRIANGLENUM)
				{
					// Atomic increment
					// Increases the value of trianglenum by one explicitly
					unsigned int  num = atomicInc(&Cell->m_TriangleNum, CELLMAXTRIANGLENUM);
					Cell->m_TriangleID[num] = TID;
				}
			}

		}
		TAB += BoundingBox->m_CellSize ;	

	}while(TAB <= Magnitude);
}
static __inline__ __device__ void RegisterTriangle(const CGridBox *const BoundingBox,const CTriangle* const Triangle, Integer ID, CCell* daCell, Integer CellNum)
{

	Scalar3 A = Triangle->Vertex[0];
	Scalar3 B = Triangle->Vertex[1];
	Scalar3 C = Triangle->Vertex[2];
	Scalar3 MidPoint = make_Scalar3((A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3, (A.z + B.z + C.z) / 3);
	const Scalar CellSize = BoundingBox->m_CellSize ;
	Scalar MagBC; 
	Scalar MagAB; 
	Scalar MagCA; 
	do{
		Scalar3 AB = B - A;//SubtractVector(&B , &A);	
		Scalar3 BA = A - B;//SubtractVector(&A , &B);
		Scalar3 BC = C - B;//SubtractVector(&C , &B);	
		Scalar3 CB = B - C;//SubtractVector(&B , &C);	
		Scalar3 CA = A - C;//SubtractVector(&A , &C);
		Scalar3 AC = C - A;//SubtractVector(&C , &A);		

		MagBC = Magnitude(&BC);
		MagAB = Magnitude(&AB);
		MagCA = Magnitude(&CA);

		if(MagAB > 0.0)
		{
			RegisterLine(BoundingBox,A, make_Scalar3(AB.x / MagAB, AB.y / MagAB, AB.z / MagAB), MagAB, ID, daCell, CellNum);
		}		
		if(MagBC > 0.0)
		{
			RegisterLine(BoundingBox,B, make_Scalar3( BC.x / MagBC, BC.y / MagBC, BC.z / MagBC), MagBC,ID, daCell, CellNum);
		}		
		if(MagCA > 0.0)
		{
			RegisterLine(BoundingBox,C, make_Scalar3(CA.x / MagCA, CA.y / MagCA, CA.z / MagCA) , MagCA,  ID , daCell, CellNum);
		}		
		Scalar3 AM = AddVector(&AB, &AC);
		Scalar3 BM = AddVector(&BA, &BC);
		Scalar3 CM = AddVector(&CA, &CB);
		Scalar MagAM = Magnitude(&AM);
		Scalar MagBM = Magnitude(&BM);
		Scalar MagCM = Magnitude(&CM);
		if(MagAM > 0.0)
		{
			AM.x /= MagAM; AM.y /= MagAM; AM.z /= MagAM;
			A.x += AM.x * CellSize;
			A.y += AM.y * CellSize;
			A.z += AM.z * CellSize;
		}
		if(MagBM > 0.0)
		{
			BM.x /= MagBM; BM.y /= MagBM; BM.z /= MagBM;
			B.x += BM.x * CellSize;
			B.y += BM.y * CellSize;
			B.z += BM.z * CellSize;
		}
		if(MagCM > 0.0)
		{
			CM.x /= MagCM; CM.y /= MagCM; CM.z /= MagCM;
			C.x += CM.x * CellSize;
			C.y += CM.y * CellSize;
			C.z += CM.z * CellSize;
		}		

	}while(MagBC >= CellSize && MagAB >= CellSize && MagCA >= CellSize);
}
//static __inline__ __device__ void RegisterLine(const CGridBox *const BoundingBox,Scalar3 Vertex, Scalar3 Direction, Scalar Magnitude, Integer TID, CCell* daCell, Integer CellNum)
//{
//	Scalar TAB = 0.0;
//	Scalar3 Target;
//	do
//	{
//		Target.x = Vertex.x + TAB * Direction.x;
//		Target.y = Vertex.y + TAB * Direction.y;
//		Target.z = Vertex.z + TAB * Direction.z;
//		Integer CID = GetCellID(BoundingBox,&Target);
//		if(CID >= 0 && CID < CellNum)
//		{			
//			CCell* Cell = &daCell[CID];
//			// Atomic compare and swap
//			// if the trianglenumber is 0, 1 is set and 0 is returned
//			unsigned int old = atomicCAS(&Cell->m_TriangleNum,0,1);
//			if(old <= 0)
//			{
//				Cell->m_TriangleID[0] = TID;
//
//			}
//			else
//			{
//				bool HasTID = false;
//				for(Integer l = 0; l < old; ++l)
//				{
//					if(Cell->m_TriangleID[l] == TID)
//					{
//						HasTID = true;
//						break;
//					}
//				}
//				if(!HasTID && old < CELLMAXTRIANGLENUM)
//				{
//					// Atomic increment
//					// Increases the value of trianglenum by one explicitly
//					unsigned int  num = atomicInc(&Cell->m_TriangleNum, CELLMAXTRIANGLENUM);
//					Cell->m_TriangleID[num] = TID;
//				}
//			}
//
//		}
//		TAB += BoundingBox->m_CellSize ;	
//
//	}while(TAB <= Magnitude);
//}
//static __inline__ __device__ void RegisterTriangle(const CGridBox *const BoundingBox,const CTriangle* const Triangle, Integer ID, CCell* daCell, Integer CellNum)
//{
//
//	Scalar3 A = Triangle->Vertex[0];
//	Scalar3 B = Triangle->Vertex[1];
//	Scalar3 C = Triangle->Vertex[2];
//	Scalar3 MidPoint = make_Scalar3((A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3, (A.z + B.z + C.z) / 3);
//	const Scalar CellSize = BoundingBox->m_CellSize ;
//	Scalar MagBC; 
//	Scalar MagAB; 
//	Scalar MagCA; 
//	do{
//		Scalar3 AB = B - A;//SubtractVector(&B , &A);	
//		Scalar3 BA = A - B;//SubtractVector(&A , &B);
//		Scalar3 BC = C - B;//SubtractVector(&C , &B);	
//		Scalar3 CB = B - C;//SubtractVector(&B , &C);	
//		Scalar3 CA = A - C;//SubtractVector(&A , &C);
//		Scalar3 AC = C - A;//SubtractVector(&C , &A);		
//
//		MagBC = Magnitude(&BC);
//		MagAB = Magnitude(&AB);
//		MagCA = Magnitude(&CA);
//
//		if(MagAB > 0.0)
//		{
//			RegisterLine(BoundingBox,A, make_Scalar3(AB.x / MagAB, AB.y / MagAB, AB.z / MagAB), MagAB, ID, daCell, CellNum);
//		}		
//		if(MagBC > 0.0)
//		{
//			RegisterLine(BoundingBox,B, make_Scalar3( BC.x / MagBC, BC.y / MagBC, BC.z / MagBC), MagBC,ID, daCell, CellNum);
//		}		
//		if(MagCA > 0.0)
//		{
//			RegisterLine(BoundingBox,C, make_Scalar3(CA.x / MagCA, CA.y / MagCA, CA.z / MagCA) , MagCA,  ID , daCell, CellNum);
//		}		
//		Scalar3 AM = AddVector(&AB, &AC);
//		Scalar3 BM = AddVector(&BA, &BC);
//		Scalar3 CM = AddVector(&CA, &CB);
//		Scalar MagAM = Magnitude(&AM);
//		Scalar MagBM = Magnitude(&BM);
//		Scalar MagCM = Magnitude(&CM);
//		if(MagAM > 0.0)
//		{
//			AM.x /= MagAM; AM.y /= MagAM; AM.z /= MagAM;
//			A.x += AM.x * CellSize;
//			A.y += AM.y * CellSize;
//			A.z += AM.z * CellSize;
//		}
//		if(MagBM > 0.0)
//		{
//			BM.x /= MagBM; BM.y /= MagBM; BM.z /= MagBM;
//			B.x += BM.x * CellSize;
//			B.y += BM.y * CellSize;
//			B.z += BM.z * CellSize;
//		}
//		if(MagCM > 0.0)
//		{
//			CM.x /= MagCM; CM.y /= MagCM; CM.z /= MagCM;
//			C.x += CM.x * CellSize;
//			C.y += CM.y * CellSize;
//			C.z += CM.z * CellSize;
//		}		
//
//	}while(MagBC >= CellSize && MagAB >= CellSize && MagCA >= CellSize);
//}

static __inline__ __device__ void RegisterTriangle(const CGridBox *const BoundingBox,const Scalar3& Vert0, const Scalar3& Vert1, const Scalar3& Vert2, Integer ID, CCell* daCell, Integer CellNum)
{

	Scalar3 A = Vert0;
	Scalar3 B = Vert1;
	Scalar3 C = Vert2;
	Scalar3 MidPoint = make_Scalar3((A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3, (A.z + B.z + C.z) / 3);
	const Scalar CellSize = BoundingBox->m_CellSize ;
	Scalar MagBC; 
	Scalar MagAB; 
	Scalar MagCA; 
	do{
		Scalar3 AB = B - A;
		Scalar3 BA = A - B;
		Scalar3 BC = C - B;
		Scalar3 CB = B - C;
		Scalar3 CA = A - C;
		Scalar3 AC = C - A;

		MagBC = Magnitude(&BC);
		MagAB = Magnitude(&AB);
		MagCA = Magnitude(&CA);

		if(MagAB > 0.0)
		{
			RegisterLine(BoundingBox,A, make_Scalar3(AB.x / MagAB, AB.y / MagAB, AB.z / MagAB), MagAB, ID, daCell, CellNum);
		}		
		if(MagBC > 0.0)
		{
			RegisterLine(BoundingBox,B, make_Scalar3( BC.x / MagBC, BC.y / MagBC, BC.z / MagBC), MagBC,ID, daCell, CellNum);
		}		
		if(MagCA > 0.0)
		{
			RegisterLine(BoundingBox,C, make_Scalar3(CA.x / MagCA, CA.y / MagCA, CA.z / MagCA) , MagCA,  ID , daCell, CellNum);
		}		
		/*Scalar3 AM = AddVector(&AB, &AC);
		Scalar3 BM = AddVector(&BA, &BC);
		Scalar3 CM = AddVector(&CA, &CB);*/
		Scalar3 AM = AB + AC; //AddVector(&AB, &AC);
		Scalar3 BM = BA + BC;//AddVector(&BA, &BC);
		Scalar3 CM = CA + CB;//AddVector(&CA, &CB);
		Scalar MagAM = Magnitude(&AM);
		Scalar MagBM = Magnitude(&BM);
		Scalar MagCM = Magnitude(&CM);
		if(MagAM > 0.0)
		{
			AM.x /= MagAM; AM.y /= MagAM; AM.z /= MagAM;
			A.x += AM.x * CellSize;
			A.y += AM.y * CellSize;
			A.z += AM.z * CellSize;
		}
		if(MagBM > 0.0)
		{
			BM.x /= MagBM; BM.y /= MagBM; BM.z /= MagBM;
			B.x += BM.x * CellSize;
			B.y += BM.y * CellSize;
			B.z += BM.z * CellSize;
		}
		if(MagCM > 0.0)
		{
			CM.x /= MagCM; CM.y /= MagCM; CM.z /= MagCM;
			C.x += CM.x * CellSize;
			C.y += CM.y * CellSize;
			C.z += CM.z * CellSize;
		}		

	}while(MagBC >= CellSize && MagAB >= CellSize && MagCA >= CellSize);
}	
static __inline__ __host__ __device__ bool IsInclude(const CGridBox* const BoundingBox, const Scalar3 * const Position,const Scalar Tolerance)
{
#if 0
	//To check ifs its finite value or not Starts
	if((isnan(Position->x)) || (isnan(Position->y)) || (isnan(Position->z)) ) 
	{
		return false;		
	}
	if((!isfinite(Position->x)) || (!isfinite(Position->y)) || (!isfinite(Position->z)))
	{
		return false;
	}
	//To check ifs its finite value or not Ends
#endif

	if(Position->x < (BoundingBox->m_ComputeZone.m_MinBound.x - Tolerance))
	{
		return false;
	}
	if(Position->y < (BoundingBox->m_ComputeZone.m_MinBound.y - Tolerance))
	{
		return false;
	}
	if(Position->z < (BoundingBox->m_ComputeZone.m_MinBound.z - Tolerance))
	{
		return false;
	}

	if(Position->x > (BoundingBox->m_ComputeZone.m_MaxBound.x + Tolerance))
	{
		return false;
	}
	if(Position->y > (BoundingBox->m_ComputeZone.m_MaxBound.y + Tolerance))
	{
		return false;
	}
	if(Position->z > (BoundingBox->m_ComputeZone.m_MaxBound.z + Tolerance))
	{
		return false;
	}
	return true;
}
#endif
//#ifdef __CUDACC__ 
//static __inline__ __device__ void RegisterLine(const CGridBox *const BoundingBox,const CGridParams* const GridParams,Scalar3 Vertex, Scalar3 Direction, Scalar Magnitude, Integer TID, CCell* daCell, Integer CellNum)
//{
//	Scalar TAB = 0.0;
//	Scalar3 Target;
//	do
//	{
//		Target.x = Vertex.x + TAB * Direction.x;
//		Target.y = Vertex.y + TAB * Direction.y;
//		Target.z = Vertex.z + TAB * Direction.z;
//		Integer CID = GetCellID(BoundingBox,&Target) ;
//		if(CID >= 0 && CID < CellNum)
//		{			
//			CCell* Cell = &daCell[CID];
//			// Atomic compare and swap
//			// if the trianglenumber is 0, 1 is set and 0 is returned
//			unsigned int old = atomicCAS(&Cell->m_TriangleNum,0,1);
//			if(old <= 0)
//			{
//				Cell->m_TriangleID[0] = TID;
//
//			}
//			else
//			{
//				bool HasTID = false;
//				for(Integer l = 0; l < old; ++l)
//				{
//					if(Cell->m_TriangleID[l] == TID)
//					{
//						HasTID = true;
//						break;
//					}
//				}
//				if(!HasTID && old < CELLMAXTRIANGLENUM)
//				{
//					// Atomic increment
//					// Increases the value of trianglenum by one explicitly
//					unsigned int  num = atomicInc(&Cell->m_TriangleNum, CELLMAXTRIANGLENUM);
//					Cell->m_TriangleID[num] = TID;
//				}
//			}
//
//		}
//		TAB += BoundingBox->m_CellSize ;	
//
//	}while(TAB <= Magnitude);
//}
//static __inline__ __device__ void RegisterTriangle(const CGridBox *const BoundingBox,const CTriangle* const Triangle, Integer ID, CCell* daCell, Integer CellNum)
//{
//
//	Scalar3 A = Triangle->Vertex[0];
//	Scalar3 B = Triangle->Vertex[1];
//	Scalar3 C = Triangle->Vertex[2];
//	Scalar3 MidPoint = make_Scalar3((A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3, (A.z + B.z + C.z) / 3);
//	const Scalar CellSize = BoundingBox->m_CellSize ;
//	Scalar MagBC; 
//	Scalar MagAB; 
//	Scalar MagCA; 
//	do{
//		Scalar3 AB = B - A;//SubtractVector(&B , &A);	
//		Scalar3 BA = A - B;//SubtractVector(&A , &B);
//		Scalar3 BC = C - B;//SubtractVector(&C , &B);	
//		Scalar3 CB = B - C;//SubtractVector(&B , &C);	
//		Scalar3 CA = A - C;//SubtractVector(&A , &C);
//		Scalar3 AC = C - A;//SubtractVector(&C , &A);		
//
//		MagBC = Magnitude(&BC);
//		MagAB = Magnitude(&AB);
//		MagCA = Magnitude(&CA);
//
//		if(MagAB > 0.0)
//		{
//			RegisterLine(BoundingBox,A, make_Scalar3(AB.x / MagAB, AB.y / MagAB, AB.z / MagAB), MagAB, ID, daCell, CellNum);
//		}		
//		if(MagBC > 0.0)
//		{
//			RegisterLine(BoundingBox,B, make_Scalar3( BC.x / MagBC, BC.y / MagBC, BC.z / MagBC), MagBC,ID, daCell, CellNum);
//		}		
//		if(MagCA > 0.0)
//		{
//			RegisterLine(BoundingBox,C, make_Scalar3(CA.x / MagCA, CA.y / MagCA, CA.z / MagCA) , MagCA,  ID , daCell, CellNum);
//		}		
//		Scalar3 AM = AddVector(&AB, &AC);
//		Scalar3 BM = AddVector(&BA, &BC);
//		Scalar3 CM = AddVector(&CA, &CB);
//		Scalar MagAM = Magnitude(&AM);
//		Scalar MagBM = Magnitude(&BM);
//		Scalar MagCM = Magnitude(&CM);
//		if(MagAM > 0.0)
//		{
//			AM.x /= MagAM; AM.y /= MagAM; AM.z /= MagAM;
//			A.x += AM.x * CellSize;
//			A.y += AM.y * CellSize;
//			A.z += AM.z * CellSize;
//		}
//		if(MagBM > 0.0)
//		{
//			BM.x /= MagBM; BM.y /= MagBM; BM.z /= MagBM;
//			B.x += BM.x * CellSize;
//			B.y += BM.y * CellSize;
//			B.z += BM.z * CellSize;
//		}
//		if(MagCM > 0.0)
//		{
//			CM.x /= MagCM; CM.y /= MagCM; CM.z /= MagCM;
//			C.x += CM.x * CellSize;
//			C.y += CM.y * CellSize;
//			C.z += CM.z * CellSize;
//		}		
//
//	}while(MagBC >= CellSize && MagAB >= CellSize && MagCA >= CellSize);
//}
//
//
//static __inline__ __device__ void RegisterTriangle(const CGridBox *const BoundingBox,const CGridParams* const GridParams,const Scalar3& Vert0, const Scalar3& Vert1, const Scalar3& Vert2, Integer ID, CCell* daCell, Integer CellNum)
//{
//
//	Scalar3 A = Vert0;
//	Scalar3 B = Vert1;
//	Scalar3 C = Vert2;
//	Scalar3 MidPoint = make_Scalar3((A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3, (A.z + B.z + C.z) / 3);
//	const Scalar CellSize = BoundingBox->m_CellSize ;
//	Scalar MagBC; 
//	Scalar MagAB; 
//	Scalar MagCA; 
//	do{
//		Scalar3 AB = B - A;
//		Scalar3 BA = A - B;
//		Scalar3 BC = C - B;
//		Scalar3 CB = B - C;
//		Scalar3 CA = A - C;
//		Scalar3 AC = C - A;
//
//		MagBC = Magnitude(&BC);
//		MagAB = Magnitude(&AB);
//		MagCA = Magnitude(&CA);
//
//		if(MagAB > 0.0)
//		{
//			RegisterLine(BoundingBox,GridParams,A, make_Scalar3(AB.x / MagAB, AB.y / MagAB, AB.z / MagAB), MagAB, ID, daCell, CellNum);
//		}		
//		if(MagBC > 0.0)
//		{
//			RegisterLine(BoundingBox,GridParams,B, make_Scalar3( BC.x / MagBC, BC.y / MagBC, BC.z / MagBC), MagBC,ID, daCell, CellNum);
//		}		
//		if(MagCA > 0.0)
//		{
//			RegisterLine(BoundingBox,GridParams,C, make_Scalar3(CA.x / MagCA, CA.y / MagCA, CA.z / MagCA) , MagCA,  ID , daCell, CellNum);
//		}		
//		Scalar3 AM = AddVector(&AB, &AC);
//		Scalar3 BM = AddVector(&BA, &BC);
//		Scalar3 CM = AddVector(&CA, &CB);
//		Scalar MagAM = Magnitude(&AM);
//		Scalar MagBM = Magnitude(&BM);
//		Scalar MagCM = Magnitude(&CM);
//		if(MagAM > 0.0)
//		{
//			AM.x /= MagAM; AM.y /= MagAM; AM.z /= MagAM;
//			A.x += AM.x * CellSize;
//			A.y += AM.y * CellSize;
//			A.z += AM.z * CellSize;
//		}
//		if(MagBM > 0.0)
//		{
//			BM.x /= MagBM; BM.y /= MagBM; BM.z /= MagBM;
//			B.x += BM.x * CellSize;
//			B.y += BM.y * CellSize;
//			B.z += BM.z * CellSize;
//		}
//		if(MagCM > 0.0)
//		{
//			CM.x /= MagCM; CM.y /= MagCM; CM.z /= MagCM;
//			C.x += CM.x * CellSize;
//			C.y += CM.y * CellSize;
//			C.z += CM.z * CellSize;
//		}		
//
//	}while(MagBC >= CellSize && MagAB >= CellSize && MagCA >= CellSize);
//}
//static __inline__ __host__ __device__ bool IsInclude(const CGridBox* const BoundingBox, const Scalar3 * const Position,const Scalar Tolerance)
//{
//#if 0
//	//To check ifs its finite value or not Starts
//	if((isnan(Position->x)) || (isnan(Position->y)) || (isnan(Position->z)) ) 
//	{
//		return false;		
//	}
//	if((!isfinite(Position->x)) || (!isfinite(Position->y)) || (!isfinite(Position->z)))
//	{
//		return false;
//	}
//	//To check ifs its finite value or not Ends
//#endif
//
//	if(Position->x < (BoundingBox->m_ComputeZone.m_MinBound.x - Tolerance))
//	{
//		return false;
//	}
//	if(Position->y < (BoundingBox->m_ComputeZone.m_MinBound.y - Tolerance))
//	{
//		return false;
//	}
//	if(Position->z < (BoundingBox->m_ComputeZone.m_MinBound.z - Tolerance))
//	{
//		return false;
//	}
//
//	if(Position->x > (BoundingBox->m_ComputeZone.m_MaxBound.x + Tolerance))
//	{
//		return false;
//	}
//	if(Position->y > (BoundingBox->m_ComputeZone.m_MaxBound.y + Tolerance))
//	{
//		return false;
//	}
//	if(Position->z > (BoundingBox->m_ComputeZone.m_MaxBound.z + Tolerance))
//	{
//		return false;
//	}
//	return true;
//}
//#endif
static __inline__ __host__ __device__ Scalar3 GetCellPosition(const CGridBox* const BoundingBox, Integer i, Integer j, Integer k)
{	
	Scalar3 V = make_Scalar3(0.0,0.0,0.0);
	V.x = i * BoundingBox->m_CellSize + BoundingBox->m_BufferedZone.m_MinBound.x;
	V.y = j * BoundingBox->m_CellSize + BoundingBox->m_BufferedZone.m_MinBound.y;
	V.z = k * BoundingBox->m_CellSize + BoundingBox->m_BufferedZone.m_MinBound.z;
	return V;
}
//static __inline__ __host__ __device__ bool IsInclude(const CBox* const BoundingBox, const Scalar3 * const Position,const Scalar Tolerance)
//{	
//	if(Position->x < BoundingBox->m_MinBound.x - Tolerance)
//	{
//		return false;
//	}
//	if(Position->y < BoundingBox->m_MinBound.y - Tolerance)
//	{
//		return false;
//	}
//	if(Position->z < BoundingBox->m_MinBound.z - Tolerance)
//	{
//		return false;
//	}
//
//	if(Position->x >= BoundingBox->m_MaxBound.x + Tolerance) 
//	{
//		return false;
//	}
//	if(Position->y >= BoundingBox->m_MaxBound.y + Tolerance)
//	{
//		return false;
//	}
//	if(Position->z >= BoundingBox->m_MaxBound.z + Tolerance)
//	{
//		return false;
//	}
//	return true;
//}

static __inline__ __host__ __device__ CBox GetSectionComputeZone(CBox ComputeZone,Integer3 Section,Integer i, Integer j, Integer k)
{
	CBox SectionBox;
	if(Section.x == 0 || Section.y == 0 || Section.z == 0)
	{
		SectionBox.m_MinBound =  make_Scalar3(0.0,0.0,0.0);
		SectionBox.m_MaxBound =  make_Scalar3(0.0,0.0,0.0);
	}
	else
	{
		Scalar3 Size = make_Scalar3(
			(ComputeZone.m_MaxBound.x - ComputeZone.m_MinBound.x) / (Scalar) Section.x,
			(ComputeZone.m_MaxBound.y - ComputeZone.m_MinBound.y) / (Scalar) Section.y,
			(ComputeZone.m_MaxBound.z - ComputeZone.m_MinBound.z) / (Scalar) Section.z
			);

		SectionBox.m_MinBound = ComputeZone.m_MinBound +  make_Scalar3(
			Size.x * (Scalar) i,
			Size.y * (Scalar) j,
			Size.z * (Scalar) k
			);
		SectionBox.m_MaxBound = ComputeZone.m_MinBound + make_Scalar3(
			Size.x * (Scalar) (i + 1),
			Size.y * (Scalar) (j + 1),
			Size.z * (Scalar) (k + 1)
			);
	}
	return SectionBox;
}

static __inline__ __host__ __device__ void CorrectVelocity(CParameter* Parameter, Scalar3 &TargetVelocity)
{	
	/*const Scalar Constant = Parameter->InitialParticleDistance * Parameter->CourantNumber / Parameter->Dt;*/
	/*Scalar NumericConst = static_cast<Scalar>(sqrt(1.0/3.0));
	Scalar Vmax = Constant * NumericConst;

	Scalar Abs = sqrt(TargetVelocity.x * TargetVelocity.x + TargetVelocity.y * TargetVelocity.y + TargetVelocity.z * TargetVelocity.z);
	if(fabs(TargetVelocity.x) > Vmax)
	{
		TargetVelocity.x = Vmax * TargetVelocity.x / fabs(TargetVelocity.x);
	}
	if(fabs(TargetVelocity.y) > Vmax)
	{
		TargetVelocity.y = Vmax *TargetVelocity.y / fabs(TargetVelocity.y);
	}
	if(fabs(TargetVelocity.z) > Vmax)
	{
		TargetVelocity.z = Vmax * TargetVelocity.z / fabs(TargetVelocity.z);
	}	*/
	const Scalar Constant = Parameter->InitialParticleDistance * Parameter->CourantNumber / Parameter->Dt;
	Scalar Abs = sqrt_utility(TargetVelocity.x * TargetVelocity.x + TargetVelocity.y * TargetVelocity.y + TargetVelocity.z * TargetVelocity.z)/*- Parameter->SonicVelocity*/;
	if(Abs > Constant)
	{
	  TargetVelocity.x = Constant * TargetVelocity.x / Abs;
	  TargetVelocity.y = Constant * TargetVelocity.y / Abs;
	  TargetVelocity.z = Constant * TargetVelocity.z / Abs;
	}
}

static __inline__ CCTStatusType CudaSafeCall(cudaError_t Status)
{		
	if(Status == cudaErrorInvalidValue )
	{
		printf("Symbol");
	}
	if(cudaSuccess != Status)
	{
		printf(cudaGetErrorString(Status));
		/*printf("\t in Function ");
		printf(Function);*/
		return CCT_CUDAERR;	
	}
	return CCT_NOERR;
}
static __inline__ __host__ void RegisterLine(const Integer TriangleNumber, const CGridBox *const BoundingBox, const CGridParams* const GridParams, Scalar3 Vertex, 
											 Scalar3 Direction, Scalar Magnitude, Integer TID, CCell* daCell, Integer CellNum)
{
	Scalar TAB = 0.0;
	Scalar3 Target;
	do
	{
		Target.x = Vertex.x + TAB * Direction.x;
		Target.y = Vertex.y + TAB * Direction.y;
		Target.z = Vertex.z + TAB * Direction.z;
		Integer CID = GetCellID(BoundingBox,&Target) ;
		if(CID >= GridParams->m_ComputeCellIdStart && CID < GridParams->m_ComputeCellIdEnd)
		{			
			CCell* Cell = &daCell[CID  - GridParams->m_BufferedCellIdStart];
			// Atomic compare and swap
			//(old == compare ? val : old)
			// if the trianglenumber is 0, 1 is set and 0 is returned
			//unsigned int old  = atomicCAS(&Cell->m_TriangleNum,0,1);
			unsigned int old = Cell->m_TriangleNum;
			if(Cell->m_TriangleNum == 0)
			{
				Cell->m_TriangleNum = 1;
			}			
			if(old <= 0)
			{
				Cell->m_TriangleID[0] = TID;
			}
			else
			{
				bool HasTID = false;
				for(Integer l = 0; l < old; ++l)
				{
					if(Cell->m_TriangleID[l] == TID)
					{
						HasTID = true;
						break;
					}
				}
				if(!HasTID && old < CELLMAXTRIANGLENUM)
				{
					// Atomic increment
					//((old >= val) ? 0 : (old+1))
					// Increases the value of trianglenum by one explicitly
					//unsigned int  num = atomicInc(&Cell->m_TriangleNum, CELLMAXTRIANGLENUM);

					unsigned int num = Cell->m_TriangleNum;
					if(Cell->m_TriangleNum >= CELLMAXTRIANGLENUM)
					{
						Cell->m_TriangleNum = 0;
					}					
					Cell->m_TriangleID[num] = TID;
				}
			}

		}
		TAB += BoundingBox->m_CellSize ;	

	}while(TAB <= Magnitude);
}
