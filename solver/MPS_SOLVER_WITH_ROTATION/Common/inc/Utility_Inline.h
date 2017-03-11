
#pragma once
#pragma warning( disable : 4305 )

#ifdef __CUDACC__ 
static __inline__ __device__ Integer CudaGetTargetID()
{
	return blockDim.x * blockIdx.x + threadIdx.x;
}
#endif
static __inline__ __device__ __host__  Scalar CalcDistance(const Scalar3* const lhs, const Scalar3* const rhs)
{
	return sqrt((lhs->x - rhs->x) * (lhs->x - rhs->x) + (lhs->y - rhs->y) * (lhs->y - rhs->y) + (lhs->z - rhs->z) * (lhs->z - rhs->z));
}
static __inline__ __device__ __host__ Scalar InnerProduct(const Scalar3* const Vector)
{
	return Vector->x * Vector->x + Vector->y * Vector->y + Vector->z * Vector->z;
}
//// Magnitude of a vector
static __inline__ __device__ __host__ Scalar Magnitude(const Scalar3* const Vector)
{
	Scalar Product =  Vector->x * Vector->x + Vector->y * Vector->y + Vector->z * Vector->z;
	return sqrt(Product);
}
static __inline__ __device__ __host__ Scalar DotProduct(const Scalar3* const P1, const Scalar3* const P2)
{
	Scalar val = P1->x * P2->x + P1->y * P2->y + P1->z * P2->z;
	return val;
}
static __inline__ __device__ __host__ Scalar SegmentRatio(const Scalar3* const P1, const Scalar3* const P2, const Scalar3* const P3)
{
	Scalar3 DenomVec = *P2 - *P1;//SubtractVector(P2, P1);
	Scalar Denom = InnerProduct(&DenomVec);
	Scalar Numer = (P3->x - P1->x) * (P2->x - P1->x) + (P3->y - P1->y) * (P2->y - P1->y) + (P3->z - P1->z) * (P2->z - P1->z);
	if(Denom == 0.0)
	{
		Denom = 1e-5;
	}
	Scalar u = Numer / Denom;
	return u;
}
//static __inline__ __device__ __host__  Scalar3 CalcDistanceVector(const Scalar3* const lhs, const Scalar3* const rhs)
//{
//	return make_Scalar3((lhs->x - rhs->x) , (lhs->y - rhs->y),  (lhs->z - rhs->z));
//}
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
static __inline__ __device__ __host__ Scalar3 CalcNormalDistance(const CTriangle* const pTriangle, const Scalar3* const Position)
{
	const Scalar3 Normal = pTriangle->Normal;

	Scalar3 Center;
	Center.x = 0;
	Center.y = 0;
	Center.z = 0;
	for(Integer i = 0; i < 3; ++i)
	{
		Center.x += pTriangle->Vertex[i].x;
		Center.y += pTriangle->Vertex[i].y;
		Center.z += pTriangle->Vertex[i].z;
	}	
	Center.x /= 3.0;
	Center.y /= 3.0;
	Center.z /= 3.0;

	const Scalar lNumerator = Normal.x * Center.x + Normal.y * Center.y + Normal.z * Center.z;
	const Scalar rNumerator = Normal.x * Position->x + Normal.y * Position->y + Normal.z * Position->z;
	const Scalar denominator = Normal.x * Normal.x + Normal.y * Normal.y + Normal.z * Normal.z;
	const Scalar sDenom = sqrt(denominator);
	Scalar3 DistanceVector = make_Scalar3(0,0,0) ;
	if(sDenom == 0.0)
	{
		return DistanceVector;
	}
	const Scalar Coefficient = (rNumerator - lNumerator) / sDenom;

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
		a = -5.0474610E+03;
		b = 0.1;
		c = 2.9725900E+01;
	}
	else if(InterpolatedPerInitial <= 0.2)
	{
		a = -1.1272200E+02;
		b = 0.2;
		c = 1.8453700E+01;
	}
	else if(InterpolatedPerInitial <= 0.3)
	{
		a = -4.3342000E+01;
		b = 0.3;
		c = 1.4119500E+01;
	}
	else if(InterpolatedPerInitial <= 0.4)
	{
		a = -2.6328000E+01;
		b = 0.4;
		c = 1.1486700E+01;
	}
	else if(InterpolatedPerInitial <= 0.5)
	{
		a = -1.9648400E+01;
		b = 0.5;
		c = 9.5218600E+00;
	}
	else if(InterpolatedPerInitial <= 0.6)
	{
		a = -1.5636200E+01;
		b = 0.6;
		c = 7.9582400E+00;
	}
	else if(InterpolatedPerInitial <= 0.7)
	{
		a = -1.2591500E+01;
		b = 0.7;
		c = 6.6990900E+00;
	}
	else if(InterpolatedPerInitial <= 0.8)
	{
		a = -1.0988800E+01;
		b = 0.8;
		c = 5.6002100E+00;
	}
	else if(InterpolatedPerInitial <= 0.9)
	{
		a = -9.0692000E+00;
		b = 0.9;
		c = 4.6932900E+00;
	}
	else if(InterpolatedPerInitial <= 1.0)
	{
		a =-7.5385000E+00;
		b = 1;
		c = 3.9394400E+00;
	}
	else if(InterpolatedPerInitial <= 1.1)
	{
		a = -6.9148000E+00;
		b = 1.1;
		c =3.2479600E+00;
	}
	else if(InterpolatedPerInitial <= 1.2)
	{
		a = -5.9144000E+00;
		b = 1.2;
		c = 2.6565200E+00;
	}
	else if(InterpolatedPerInitial <= 1.3)
	{
		a = -5.4670000E+00;
		b = 1.3;
		c = 2.1098200E+00;
	}
	else if(InterpolatedPerInitial <= 1.4)
	{
		a = -5.0626000E+00;
		b = 1.4;
		c = 1.6035600E+00;
	}
	else if(InterpolatedPerInitial <= 1.5)
	{
		a = -4.6948000E+00;
		b = 1.5;
		c = 1.1340800E+00;
	}
	else if(InterpolatedPerInitial <= 1.6)
	{
		a = -3.6959000E+00;
		b = 1.6;
		c = 7.6449000E-01;
	}
	else if(InterpolatedPerInitial <= 1.7)
	{
		a = -2.7022800E+00;
		b = 1.7;
		c = 4.9426200E-01;
	}
	else if(InterpolatedPerInitial <= 1.8)
	{
		a = -2.4819500E+00;
		b = 1.8;
		c = 2.4606700E-01;
	}
	else if(InterpolatedPerInitial <= 1.9)
	{
		a = -1.4080400E+00;
		b = 1.9;
		c = 1.0526300E-01;
	}
	else if(InterpolatedPerInitial <= 2.0)
	{
		a = -5.5263200E-01;
		b = 2;
		c = 4.9999800E-02;
	}
	else if(InterpolatedPerInitial <= 2.1)
	{
		// 誤差を考慮して2.15くらいまで
		a = -4.9999800E-01;
		b = 2.1;
		c = 0.0000000E+00;
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

		a = -9.996E+05;
		b = 0.01;
		c = 3.837E+02;
	}
	else if(InterpolatedPerInitial <= 0.05)
	{//0.01	<x/l_0=<	0.05	-3.027E+02	0.05	8.101E+01

		a = -3.027E+02;
		b = 0.05;
		c = 8.101E+01;
	}
	else if(InterpolatedPerInitial <= 0.1)
	{//0.05	<x/l_0=<	0.1	-3.812E+01	0.1	4.289E+01

		a = 3.812E+01;
		b = 0.1;
		c = 4.289E+01;
	}
	else if(InterpolatedPerInitial <= 0.2)
	{//0.1	<x/l_0=<	0.2	-1.950E+01	0.2	2.339E+01

		a = -1.950E+01;
		b = 0.2;
		c = 2.339E+01;
	}
	else if(InterpolatedPerInitial <= 0.3)
	{//0.2	<x/l_0=<	0.3	-6.917E+00	0.3	1.647E+01

		a = -6.917E+00;
		b = 0.3;
		c = 1.647E+01;
	}
	else if(InterpolatedPerInitial <= 0.4)
	{//0.3	<x/l_0=<	0.4	-3.779E+00	0.4	1.269E+01

		a = -3.779E+00;
		b = 0.4;
		c = 1.269E+01;
	}
	else if(InterpolatedPerInitial <= 0.5)
	{//0.4	<x/l_0=<	0.5	-2.523E+00	0.5	1.017E+01

		a = -2.523E+00;
		b = 0.5;
		c = 1.017E+01;
	}
	else if(InterpolatedPerInitial <= 0.6)
	{//0.5	<x/l_0=<	0.6	-1.889E+00	0.6	8.282E+00

		a = -1.889E+00;
		b = 0.6;
		c = 8.282E+00;
	}
	else if(InterpolatedPerInitial <= 0.7)
	{//0.6	<x/l_0=<	0.7	-1.484E+00	0.7	6.798E+00

		a = -1.484E+00;
		b = 0.7;
		c = 6.798E+00;
	}
	else if(InterpolatedPerInitial <= 0.8)
	{//0.7	<x/l_0=<	0.8	-1.214E+00	0.8	5.584E+00

		a = -1.214E+00;
		b = 0.8;
		c = 5.584E+00;
	}
	else if(InterpolatedPerInitial <= 0.9)
	{//0.8	<x/l_0=<	0.9	-9.752E-01	0.9	4.609E+00

		a = -9.752E-01;
		b = 0.9;
		c = 4.609E+00;
	}
	else if(InterpolatedPerInitial <= 1.0)
	{//0.9	<x/l_0=<	1	-7.931E-01	1	3.816E+00

		a = -7.931E-01;
		b = 1;
		c = 3.816E+00;
	}
	else if(InterpolatedPerInitial <= 1.1)
	{//1	<x/l_0=<	1.1	-6.940E-01	1.1	3.122E+00

		a = -6.940E-01;
		b = 1.1;
		c = 3.122E+00;
	}
	else if(InterpolatedPerInitial <= 1.2)
	{//1.1	<x/l_0=<	1.2	-5.320E-01	1.2	2.590E+00

		a = -5.320E-01;
		b = 1.2;
		c = 2.590E+00;
	}
	else if(InterpolatedPerInitial <= 1.3)
	{//1.2	<x/l_0=<	1.3	-4.727E-01	1.3	2.117E+00

		a = -4.727E-01;
		b = 1.3;
		c = 2.117E+00;
	}
	else if(InterpolatedPerInitial <= 1.4)
	{//1.3	<x/l_0=<	1.4	-4.229E-01	1.4	1.694E+00

		a = -4.229E-01;
		b = 1.4;
		c = 1.694E+00;
	}
	else if(InterpolatedPerInitial <= 1.5)
	{//1.4	<x/l_0=<	1.5	-3.806E-01	1.5	1.314E+00

		a = -3.806E-01;
		b = 1.5;
		c = 1.314E+00;
	}
	else if(InterpolatedPerInitial <= 1.6)
	{//1.5	<x/l_0=<	1.6	-3.442E-01	1.6	9.693E-01

		a = -3.442E-01;
		b = 1.6;
		c = 9.693E-01;
	}
	else if(InterpolatedPerInitial <= 1.7)
	{//1.6	<x/l_0=<	1.7	-3.127E-01	1.7	6.566E-01

		a = -3.127E-01;
		b = 1.7;
		c = 6.566E-01;
	}
	else if(InterpolatedPerInitial <= 1.8)
	{//1.7	<x/l_0=<	1.8	-2.851E-01	1.8	3.715E-01

		a = -2.851E-01;
		b = 1.8;
		c = 3.715E-01;
	}
	else if(InterpolatedPerInitial <= 1.9)
	{//1.8	<x/l_0=<	1.9	-1.820E-01	1.9	1.895E-01

		a = -1.820E-01;
		b = 1.9;
		c = 1.895E-01;
	}
	else if(InterpolatedPerInitial <= 2.0)
	{//1.9	<x/l_0=<	2	-9.947E-02	2	9.000E-02

		a = -9.947E-02;
		b = 2;
		c = 9.000E-02;
	}
	else if(InterpolatedPerInitial <= 2.1)
	{//2	<x/l_0=<	2.1	-9.000E-02	2.1	0.000E+00

		// 誤差を考慮して2.15くらいまで
		a = -9.000E-02;
		b = 2.1;
		c = 0.0000000E+00;
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
		a = -9.6379E+03;
		b = 0.1 ;
		c = 1.1497E+02;
	}
	else if(InterpolatedPerInitial <= 0.2)
	{
		a =-2.3979E+02;
		b= 0.2 ;
		c = 9.0990E+01;
	}
	else if(InterpolatedPerInitial <= 0.3)
	{
		a = -1.0777E+02;
		b = 0.3 ;
		c = 8.0213E+01;
	}
	else if(InterpolatedPerInitial <= 0.4)
	{
		a =-7.4088E+01; 
		b= 0.4 ; 
		c =	7.2804E+01;

	}
	else if(InterpolatedPerInitial <= 0.5)
	{
		a = -6.0259E+01;
		b = 0.5 ;
		c = 6.6778E+01;
	}
	else if(InterpolatedPerInitial <= 0.6)
	{
		a = -5.2803E+01;
		b = 0.6 ;
		c = 6.1498E+01;

	}
	else if(InterpolatedPerInitial <= 0.7)
	{
		a = -4.8681E+01;
		b = 0.700 ;
		c = 5.6630E+01;
	}
	else if(InterpolatedPerInitial <= 0.8)
	{
		a = -4.4630E+01;
		b = 0.800 ;
		c = 5.2167E+01;

	}
	else if(InterpolatedPerInitial <= 0.9)
	{
		a = -4.1171E+01;
		b = 0.900 ;
		c = 4.8050E+01;
	}
	else if(InterpolatedPerInitial <= 1.0)
	{
		a = -3.8583E+01;
		b = 1;
		c = 4.4191E+01;

	}
	else if(InterpolatedPerInitial <= 1.1)
	{
		a = -3.6846E+01;
		b = 1.100 ;
		c = 4.0507E+01;
	}
	else if(InterpolatedPerInitial <= 1.2)
	{
		a = -3.5476E+01; 
		b = 	1.2;
		c = 3.6959E+01;

	}
	else if(InterpolatedPerInitial <= 1.3)
	{
		a = -3.4183E+01;
		b = 1.300 ;
		c = 3.3541E+01;
	}
	else if(InterpolatedPerInitial <= 1.4)
	{
		a = -3.1589E+01;
		b = 1.4;
		c = 3.0382E+01;

	}
	else if(InterpolatedPerInitial <= 1.5)
	{
		a = -2.9258E+01;
		b = 1.500 ;
		c = 2.7456E+01;
	}
	else if(InterpolatedPerInitial <= 1.6)
	{
		a = -2.7047E+01;
		b = 1.6;
		c = 2.4751E+01;
	}
	else if(InterpolatedPerInitial <= 1.7)
	{
		a = -2.5690E+01;
		b = 1.700 ;
		c = 2.2182E+01;
	}
	else if(InterpolatedPerInitial <= 1.8)
	{
		a = -2.3305E+01;
		b = 1.8;
		c = 1.9852E+01;

	}
	else if(InterpolatedPerInitial <= 1.9)
	{
		a = -2.0982E+01;
		b = 1.900 ;
		c = 1.7754E+01;
	}
	else if(InterpolatedPerInitial <= 2.0)
	{
		a = -1.9268E+01;
		b = 2;
		c = 1.5827E+01;

	}	
	else if(InterpolatedPerInitial <= 2.4)
	{
		a = -1.6887E+01;
		b = 2.4;
		c = 9.0719E+00;

	}	
	else if(InterpolatedPerInitial <= 2.8)
	{
		a = -1.1091E+01; 
		b = 2.8;
		c = 4.6355E+00;

	}
	else if(InterpolatedPerInitial <= 3.2)
	{
		a = -6.5068E+00;
		b = 3.2;
		c = 2.0327E+00;

	}	
	else if(InterpolatedPerInitial <= 3.6)
	{
		a = -3.7565E+00;
		b = 3.6;
		c = 5.3011E-01;

	}	
	else if(InterpolatedPerInitial <= 4.0)
	{
		a = -1.3253E+00;
		b = 4;
		c = 0.0000E+00;

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
		a = -9.992E+05;
		b = 0.01;
		c = 7.507E+02;
	}
	else if(InterpolatedPerInitial <= 0.05)
	{//0.01	<x/l_0=<	0.05	-5.766E+02	0.05	1.741E+02
		a = -5.766E+02;
		b = 0.05 ;
		c = 1.741E+02;
	}
	else if(InterpolatedPerInitial <= 0.1)
	{//0.05	<x/l_0=<	0.1	-7.280E+01	0.1	1.013E+02
		a = -7.280E+01;
		b = 0.1 ;
		c = 1.013E+02;
	}
	else if(InterpolatedPerInitial <= 0.2)
	{//0.1	<x/l_0=<	0.2	-3.761E+01	0.2	6.368E+01
		a =-3.761E+01;
		b= 0.2 ;
		c = 6.368E+01;
	}	
	else if(InterpolatedPerInitial <= 0.4)
	{//0.2	<x/l_0=<	0.4	-2.128E+01	0.4	4.240E+01
		a =-2.128E+01; 
		b= 0.4 ; 
		c =	4.240E+01;
	}	
	else if(InterpolatedPerInitial <= 0.6)
	{//0.4	<x/l_0=<	0.6	-9.280E+00	0.6	3.312E+01
		a =-9.280E+00; 
		b= 0.6 ; 
		c =	3.312E+01;
	}	
	else if(InterpolatedPerInitial <= 0.8)
	{//0.6	<x/l_0=<	0.8	-6.104E+00	0.8	2.701E+01
		a = -6.104E+00;
		b = 0.8 ;
		c = 2.701E+01;
	}	
	else if(InterpolatedPerInitial <= 1.0)
	{//0.8	<x/l_0=<	1	-4.682E+00	1	2.233E+01
		a = -4.682E+00;
		b = 1;
		c = 2.233E+01;
	}	
	else if(InterpolatedPerInitial <= 1.2)
	{//1	<x/l_0=<	1.2	-3.842E+00	1.2	1.849E+01
		a = -3.842E+00; 
		b = 	1.2;
		c = 1.849E+01;
	}	
	else if(InterpolatedPerInitial <= 1.4)
	{//1.2	<x/l_0=<	1.4	-3.260E+00	1.4	1.523E+01
		a = -3.260E+00;
		b = 1.4;
		c = 1.523E+01;
	}	
	else if(InterpolatedPerInitial <= 1.6)
	{//1.4	<x/l_0=<	1.6	-2.763E+00	1.6	1.247E+01
		a = -2.763E+00;
		b = 1.6;
		c = 1.247E+01;
	}	
	else if(InterpolatedPerInitial <= 1.8)
	{//1.6	<x/l_0=<	1.8	-2.048E+00	1.8	1.042E+01
		a = -2.048E+00;
		b = 1.8;
		c = 1.042E+01;
	}	
	else if(InterpolatedPerInitial <= 2.0)
	{//1.8	<x/l_0=<	2	-2.034E+00	2	8.385E+00
		a = -2.034E+00;
		b = 2;
		c = 8.385E+00;
	}	
	else if(InterpolatedPerInitial <= 2.2)
	{//2	<x/l_0=<	2.2	-1.726E+00	2.2	6.659E+00
		a = -1.726E+00;
		b = 2.2;
		c = 6.659E+00;
	}	
	else if(InterpolatedPerInitial <= 2.4)
	{//2.2	<x/l_0=<	2.4	-1.543E+00	2.4	5.116E+00
		a = -1.543E+00;
		b = 2.4;
		c = 5.116E+00;
	}	
	else if(InterpolatedPerInitial <= 2.6)
	{//2.4	<x/l_0=<	2.6	-1.306E+00	2.6	3.810E+00
		a = -1.306E+00;
		b = 2.6;
		c = 3.810E+00;
	}	
	else if(InterpolatedPerInitial <= 2.8)
	{//2.6	<x/l_0=<	2.8	-1.051E+00	2.8	2.759E+00
		a = -1.051E+00; 
		b = 2.8;
		c = 2.759E+00;
	}
	else if(InterpolatedPerInitial <= 3.0)
	{//2.8	<x/l_0=<	3	-8.115E-01	3	1.948E+00
		a = -8.115E-01; 
		b = 3;
		c = 1.948E+00;
	}
	else if(InterpolatedPerInitial <= 3.2)
	{//3	<x/l_0=<	3.2	-5.864E-01	3.2	1.361E+00
		a = -5.864E-01;
		b = 3.2;
		c = 1.361E+00;
	}	
	else if(InterpolatedPerInitial <= 3.4)
	{//3.2	<x/l_0=<	3.4	-5.298E-01	3.4	8.314E-01
		a = -5.298E-01;
		b = 3.4;
		c = 8.314E-01;
	}	
	else if(InterpolatedPerInitial <= 3.6)
	{//3.4	<x/l_0=<	3.6	-3.773E-01	3.6	4.541E-01
		a = -3.773E-01;
		b = 3.6;
		c = 4.541E-01;
	}	
	else if(InterpolatedPerInitial <= 3.8)
	{//3.6	<x/l_0=<	3.8	-2.946E-01	3.8	1.594E-01
		a = -2.946E-01;
		b = 3.8;
		c = 1.594E-01;
	}	
	else if(InterpolatedPerInitial <= 4.0)
	{//3.8	<x/l_0=<	4	-1.594E-01	4	0.000E+00
		a = -1.594E-01;
		b = 4;
		c = 0.0000E+00;
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
		a = -4.073872;
		b = 0.25;
		c = 15.144314;
	}
	else if(InterpolatedPerInitial <= 0.5)
	{//0.25	<r/l_0=<	0.5	-8.203804	0.5	13.093363
		a = -8.203804;
		b = 0.5;
		c = 13.093363;
	}
	else if(InterpolatedPerInitial <= 0.75)
	{//0.5	<r/l_0=<	0.75	-8.047464	0.75	11.081497
		a = -8.047464;
		b = 0.75;
		c = 11.081497;
	}
	else if(InterpolatedPerInitial <= 1.0)
	{//0.75	<r/l_0=<	1	-7.628876	1	9.174278
		a = -7.628876;
		b = 1.0;
		c = 9.174278;
	}
	else if(InterpolatedPerInitial <= 1.25)
	{//1	<r/l_0=<	1.25	-7.001212	1.25	7.423975

		a = -7.001212;
		b = 1.25;
		c = 7.423975;
	}
	else if(InterpolatedPerInitial <= 1.5)
	{//1.25	<r/l_0=<	1.5	-6.230716	1.5	5.866296

		a = -6.230716;
		b = 1.5;
		c = 5.866296;
	}
	else if(InterpolatedPerInitial <= 1.75)
	{//1.5	<r/l_0=<	1.75	-5.392776	1.75	4.518102

		a = -5.392776;
		b = 1.75;
		c = 4.518102;
	}
	else if(InterpolatedPerInitial <= 2.0)
	{//1.75	<r/l_0=<	2	-4.541608	2	3.3827

		a = -4.541608;
		b = 2.0;
		c = 3.3827;
	}
	else if(InterpolatedPerInitial <= 2.25)
	{//2	<r/l_0=<	2.25	-3.715672	2.25	2.453782

		a = -3.715672;
		b = 2.25;
		c = 2.453782;
	}
	else if(InterpolatedPerInitial <= 2.5)
	{//2.25	<r/l_0=<	2.5	-2.94322	2.5	1.717977

		a = -2.94322;
		b = 2.5;
		c = 1.717977;
	}
	else if(InterpolatedPerInitial <= 2.75)
	{//2.5	<r/l_0=<	2.75	-2.256608	2.75	1.153825

		a = -2.256608;
		b = 2.75;
		c = 1.153825;
	}
	else if(InterpolatedPerInitial <= 3.0)
	{//2.75	<r/l_0=<	3	-1.66868	3	0.736655

		a = -1.66868;
		b = 3.0;
		c = 0.736655;
	}
	else if(InterpolatedPerInitial <= 3.25)
	{//3	<r/l_0=<	3.25	-1.180056	3.25	0.441641

		a = -1.180056;
		b = 3.25;
		c = 0.441641;
	}
	else if(InterpolatedPerInitial <= 3.5)
	{//3.25	<r/l_0=<	3.5	-0.785508	3.5	0.245264

		a = -0.785508;
		b = 3.5;
		c = 0.245264;
	}
	else if(InterpolatedPerInitial <= 3.75)
	{//3.5	<r/l_0=<	3.75	-0.489308	3.75	0.122937

		a = -0.489308;
		b = 3.75;
		c = 0.122937;
	}
	else if(InterpolatedPerInitial <= 4.0)
	{//3.75	<r/l_0=<	4	-0.280256	4	0.052873

		a = -0.280256;
		b = 4.0;
		c = 0.052873;
	}
	else if(InterpolatedPerInitial <= 4.25)
	{//4	<r/l_0=<	4.25	-0.141756	4.25	0.017434

		a = -0.141756;
		b = 4.25;
		c = 0.017434;
	}
	else if(InterpolatedPerInitial <= 4.5)
	{//4.25	<r/l_0=<	4.5	-0.055172	4.5	0.003641

		a = -0.055172;
		b = 4.5;
		c = 0.003641;
	}
	else if(InterpolatedPerInitial <= 4.75)
	{//4.5	<r/l_0=<	4.75	-0.013656	4.75	0.000227

		a = -0.013656;
		b = 4.75;
		c = 0.000227;
	}
	else if(InterpolatedPerInitial <= 5.0)
	{//4.75	<r/l_0=<	5	-0.000908	5	0

		a = -0.000908;
		b = 5.0;
		c = 0.0;
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

		a = -1.85978;
		b = 0.05;
		c = 1.902776;
	}
	else if(InterpolatedPerInitial <= 0.1)
	{//-2.27896	0.1	1.788828

		a = -2.27896;
		b = 0.1;
		c = 1.788828;
	}
	else if(InterpolatedPerInitial <= 0.15)
	{//-2.22854	0.15	1.677401

		a = -2.22854;
		b = 0.15;
		c = 1.677401;
	}
	else if(InterpolatedPerInitial <= 0.2)
	{//-2.17656	0.2	1.568573

		a = -2.17656;
		b = 0.2;
		c = 1.568573;
	}
	else if(InterpolatedPerInitial <= 0.25)
	{//-2.1219	0.25	1.462478

		a = -2.1219;
		b = 0.25;
		c = 1.462478;
	}
	else if(InterpolatedPerInitial <= 0.3)
	{//-2.06388	0.3	1.359284

		a = -2.06388;
		b = 0.3;
		c = 1.359284;
	}
	else if(InterpolatedPerInitial <= 0.35)
	{//-2.00216	0.35	1.259176

		a = -2.00216;
		b = 0.35;
		c = 1.259176;
	}
	else if(InterpolatedPerInitial <= 0.4)
	{//-1.93652	0.4	1.16235

		a = -1.93652;
		b = 0.4;
		c = 1.16235;
	}
	else if(InterpolatedPerInitial <= 0.45)
	{//-1.8672	0.45	1.06899

		a = -1.8672;
		b = 0.45;
		c = 1.06899;
	}
	else if(InterpolatedPerInitial <= 0.5)
	{//-1.79454	0.5	0.979263

		a = -1.79454;
		b = 0.5;
		c = 0.979263;
	}
	else if(InterpolatedPerInitial <= 0.55)
	{//-1.71922	0.55	0.893302

		a = -1.71922;
		b = 0.55;
		c = 0.893302;
	}
	else if(InterpolatedPerInitial <= 0.6)
	{//-1.64176	0.6	0.811214

		a = -1.64176;
		b = 0.6;
		c = 0.811214;
	}
	else if(InterpolatedPerInitial <= 0.65)
	{//-1.5613	0.65	0.733149

		a = -1.5613;
		b = 0.65;
		c = 0.733149;
	}
	else if(InterpolatedPerInitial <= 0.7)
	{//-1.4783	0.7	0.659234

		a = -1.4783;
		b = 0.7;
		c = 0.659234;
	}
	else if(InterpolatedPerInitial <= 0.75)
	{//-1.39368	0.75	0.58955

		a = -1.39368;
		b = 0.75;
		c = 0.58955;
	}
	else if(InterpolatedPerInitial <= 0.8)
	{//-1.30852	0.8	0.524124

		a = -1.30852;
		b = 0.8;
		c = 0.524124;
	}
	else if(InterpolatedPerInitial <= 0.85)
	{//-1.2239	0.85	0.462929

		a = -1.2239;
		b = 0.85;
		c = 0.462929;
	}
	else if(InterpolatedPerInitial <= 0.9)
	{//-1.14014	0.9	0.405922

		a = -1.14014;
		b = 0.9;
		c = 0.405922;
	}
	else if(InterpolatedPerInitial <= 0.95)
	{//-1.0554	0.95	0.353152

		a = -1.0554;
		b = 0.95;
		c = 0.353152;
	}
	else if(InterpolatedPerInitial <= 1.0)
	{//-0.97002	1	0.304651

		a = -0.97002;
		b = 1.0;
		c = 0.304651;
	}
	else if(InterpolatedPerInitial <= 1.05)
	{//-0.88506	1.05	0.260398

		a = -0.88506;
		b = 1.05;
		c = 0.260398;
	}
	else if(InterpolatedPerInitial <= 1.1)
	{//	-0.80146	1.1	0.220325
	
		a = -0.80146;
		b = 1.1;
		c = 0.220325;
	}
	else if(InterpolatedPerInitial <= 1.15)
	{//	-0.71984	1.15	0.184333
	
		a = -0.71984;
		b = 1.15;
		c = 0.184333;
	}
	else if(InterpolatedPerInitial <= 1.2)
	{// -0.64004	1.2	0.152331
		
		a = -0.64004;
		b = 1.2;
		c = 0.152331;
	}
	else if(InterpolatedPerInitial <= 1.25)
	{//	-0.56268	1.25	0.124197
	
		a = -0.56268;
		b = 1.25;
		c = 0.124197;
	}
	else if(InterpolatedPerInitial <= 1.3)
	{//	-0.48848	1.3	0.099773
	
		a = -0.48848;
		b = 1.3;
		c = 0.099773;
	}
	else if(InterpolatedPerInitial <= 1.35)
	{//	-0.41824	1.35	0.078861
	
		a = -0.41824;
		b = 1.35;
		c =0.078861;
	}
	else if(InterpolatedPerInitial <= 1.4)
	{//	-0.3527	1.4	0.061226
	
		a = -0.3527;
		b = 1.4;
		c = 0.061226;
	}
	else if(InterpolatedPerInitial <= 1.45)
	{// -0.2925	1.45	0.046601
		
		a = -0.2925;
		b = 1.45;
		c = 0.046601;
	}
	else if(InterpolatedPerInitial <= 1.5)
	{// -0.23824	1.5	0.034689
		
		a = -0.23824;
		b = 1.5;
		c = 0.034689;
	}
	else if(InterpolatedPerInitial <= 1.55)
	{//	-0.19054	1.55	0.025162
	
		a = -0.19054;
		b = 1.55;
		c = 0.025162;
	}
	else if(InterpolatedPerInitial <= 1.6)
	{//	-0.14958	1.6	0.017683
	
		a = -0.14958;
		b = 1.6;
		c = 0.017683;
	}
	else if(InterpolatedPerInitial <= 1.65)
	{//	-0.11396	1.65	0.011985
	
		a = -0.11396;
		b = 1.65;
		c = 0.011985;
	}
	else if(InterpolatedPerInitial <= 1.7)
	{// -0.0835	1.7	0.00781
		
		a = -0.0835;
		b = 1.7;
		c = 0.00781;
	}
	else if(InterpolatedPerInitial <= 1.75)
	{//	-0.0584	1.75	0.00489
	
		a = -0.0584;
		b = 1.75;
		c = 0.00489;
	}
	else if(InterpolatedPerInitial <= 1.8)
	{//	-0.0389	1.8	0.002945
	
		a = -0.0389;
		b = 1.8;
		c = 0.002945;
	}
	else if(InterpolatedPerInitial <= 1.85)
	{//	-0.02516	1.85	0.001687
	
		a = -0.02516;
		b = 1.85;
		c = 0.001687;
	}
	else if(InterpolatedPerInitial <= 1.9)
	{//	-0.01646	1.9	0.000864
	
		a = -0.01646;
		b = 1.9;
		c = 0.000864;
	}
	else if(InterpolatedPerInitial <= 1.95)
	{// -0.01	1.95	0.000364
		
		a = -0.01;
		b = 1.95;
		c = 0.000364;
	}
	else if(InterpolatedPerInitial <= 2.0)
	{//	-0.00512	2	0.000108
	
		a = -0.00512;
		b = 2.0;
		c = 0.000108;
	}
	else if(InterpolatedPerInitial <= 2.05)
	{// -0.0019	2.05	0.000013
		
		a = -0.0019;
		b = 2.05;
		c = 0.000013;
	}
	else if(InterpolatedPerInitial <= 2.1)
	{// -0.00026	2.1	0
		
		a = -0.00026;
		b = 2.1;
		c = 0.0;
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

		a = -2.22678;
		b = 0.05;
		c = 3.220183;
	}
	else if(InterpolatedPerInitial <= 0.1)
	{//-2.80976	0.1	3.079695

		a = -2.80976;
		b = 0.1;
		c = 3.079695;
	}
	else if(InterpolatedPerInitial <= 0.15)
	{//-2.83986	0.15	2.937702

		a = -2.83986;
		b = 0.15;
		c = 2.937702;
	}
	else if(InterpolatedPerInitial <= 0.2)
	{//-2.8687	0.2	2.794267

		a = -2.8687;
		b = 0.2;
		c = 2.794267;
	}
	else if(InterpolatedPerInitial <= 0.25)
	{// -2.8836	0.25	2.650087

		a = -2.8836;
		b = 0.25;
		c = 2.650087;
	}
	else if(InterpolatedPerInitial <= 0.3)
	{// -2.8833	0.3	2.505922

		a = -2.8833;
		b = 0.3;
		c = 2.505922;
	}
	else if(InterpolatedPerInitial <= 0.35)
	{// -2.86688	0.35	2.362578

		a = -2.86688;
		b = 0.35;
		c = 2.362578;
	}
	else if(InterpolatedPerInitial <= 0.4)
	{// -2.83348	0.4	2.220904

		a = -2.83348;
		b = 0.4;
		c = 2.220904;
	}
	else if(InterpolatedPerInitial <= 0.45)
	{// -2.78258	0.45	2.081775

		a = -2.78258;
		b = 0.45;
		c = 2.081775;
	}
	else if(InterpolatedPerInitial <= 0.5)
	{// -2.71386	0.5	1.946082

		a = -2.71386;
		b = 0.5;
		c = 1.946082;
	}
	else if(InterpolatedPerInitial <= 0.55)
	{// -2.62708	0.55	1.814728

		a = -2.62708;
		b = 0.55;
		c = 1.814728;
	}
	else if(InterpolatedPerInitial <= 0.6)
	{// -2.54514	0.6	1.687471

		a = -2.54514;
		b = 0.6;
		c = 1.687471;
	}
	else if(InterpolatedPerInitial <= 0.65)
	{// -2.47432	0.65	1.563755

		a = -2.47432;
		b = 0.65;
		c = 1.563755;
	}
	else if(InterpolatedPerInitial <= 0.7)
	{// -2.39412	0.7	1.444049

		a = -2.39412;
		b = 0.7;
		c = 1.444049;
	}
	else if(InterpolatedPerInitial <= 0.75)
	{// -2.30322	0.75	1.328888

		a = -2.30322;
		b = 0.75;
		c = 1.328888;
	}
	else if(InterpolatedPerInitial <= 0.8)
	{// -2.19946	0.8	1.218915

		a = -2.19946;
		b = 0.8;
		c = 1.218915;
	}
	else if(InterpolatedPerInitial <= 0.85)
	{// -2.08352	0.85	1.114739

		a = -2.08352;
		b = 0.85;
		c = 1.114739;
	}
	else if(InterpolatedPerInitial <= 0.9)
	{// -1.99556	0.9	1.014961

		a = -1.99556;
		b = 0.9;
		c = 1.014961;
	}
	else if(InterpolatedPerInitial <= 0.95)
	{// -1.92792	0.95	0.918565

		a = -1.92792;
		b = 0.95;
		c = 0.918565;
	}
	else if(InterpolatedPerInitial <= 1.0)
	{// -1.85086	1	0.826022

		a = -1.85086;
		b = 1.0;
		c = 0.826022;
	}
	else if(InterpolatedPerInitial <= 1.05)
	{// -1.76478	1.05	0.737783

		a = -1.76478;
		b = 1.05;
		c = 0.737783;
	}
	else if(InterpolatedPerInitial <= 1.1)
	{//	-1.67026	1.1	0.65427
	
		a = -1.67026;
		b = 1.1;
		c = 0.65427;
	}
	else if(InterpolatedPerInitial <= 1.15)
	{// -1.5791	1.15	0.575315
		
		a = -1.5791;
		b = 1.15;
		c = 0.575315;
	}
	else if(InterpolatedPerInitial <= 1.2)
	{//	-1.49178	1.2	0.500726
	
		a = -1.49178;
		b = 1.2;
		c = 0.500726;
	}
	else if(InterpolatedPerInitial <= 1.25)
	{// -1.3974	1.25	0.430856
		
		a = -1.3974;
		b = 1.25;
		c = 0.430856;
	}
	else if(InterpolatedPerInitial <= 1.3)
	{// -1.29646	1.3	0.366033
		
		a = -1.29646;
		b = 1.3;
		c = 0.366033;
	}
	else if(InterpolatedPerInitial <= 1.35)
	{//	-1.18932	1.35	0.306567
	
		a = -1.18932;
		b = 1.35;
		c = 0.306567;
	}
	else if(InterpolatedPerInitial <= 1.4)
	{// -1.07642	1.4	0.252746
		
		a = -1.07642;
		b = 1.4;
		c = 0.252746;
	}
	else if(InterpolatedPerInitial <= 1.45)
	{// -0.95812	1.45	0.20484
		
		a = -0.95812;
		b = 1.45;
		c = 0.20484;
	}
	else if(InterpolatedPerInitial <= 1.5)
	{// -0.83474	1.5	0.163103
		
		a = -0.83474;
		b = 1.5;
		c = 0.163103;
	}
	else if(InterpolatedPerInitial <= 1.55)
	{// -0.7067	1.55	0.127768
		
		a = -0.7067;
		b = 1.55;
		c = 0.127768;
	}
	else if(InterpolatedPerInitial <= 1.6)
	{// -0.59698	1.6	0.097919

		a =  -0.59698;
		b = 1.6;
		c = 0.097919;
	}
	else if(InterpolatedPerInitial <= 1.65)
	{// -0.5119	1.65	0.072324
		
		a = -0.5119;
		b = 1.65;
		c = 0.072324;
	}
	else if(InterpolatedPerInitial <= 1.7)
	{// -0.42508	1.7	0.05107
		
		a = -0.42508;
		b = 1.7;
		c = 0.05107;
	}
	else if(InterpolatedPerInitial <= 1.75)
	{// -0.33652	1.75	0.034244
		
		a = -0.33652;
		b = 1.75;
		c = 0.034244;
	}
	else if(InterpolatedPerInitial <= 1.8)
	{// -0.2464	1.8	0.021924
		
		a = -0.2464;
		b = 1.8;
		c = 0.021924;
	}
	else if(InterpolatedPerInitial <= 1.85)
	{// -0.15504	1.85	0.014172
		
		a = -0.15504;
		b = 1.85;
		c = 0.014172;
	}
	else if(InterpolatedPerInitial <= 1.9)
	{// -0.10204	1.9	0.00907
		
		a = -0.10204;
		b = 1.9;
		c = 0.00907;
	}
	else if(InterpolatedPerInitial <= 1.95)
	{// -0.07936	1.95	0.005102
		
		a = -0.07936;
		b = 1.95;
		c = 0.005102;
	}
	else if(InterpolatedPerInitial <= 2.0)
	{// -0.05668	2	0.002268
		
		a = -0.05668;
		b = 2.0;
		c = 0.002268;
	}
	else if(InterpolatedPerInitial <= 2.05)
	{// -0.03402	2.05	0.000567
		
		a = -0.03402;
		b = 2.05;
		c = 0.000567;
	}
	else if(InterpolatedPerInitial <= 2.1)
	{// -0.01134	2.1	0
		
		a = -0.01134;
		b = 2.1;
		c = 0.0;
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
	Scalar Relativez = Position->z - BoundingBox->m_BufferedZone.m_MinBound.z /*+ EPS*/;

	RelativeX /= BoundingBox->m_CellSize;
	RelativeY /= BoundingBox->m_CellSize;
	Relativez /= BoundingBox->m_CellSize;

	i = (Integer)RelativeX;
	j = (Integer)RelativeY;
	k = (Integer)Relativez;
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
	//(Section.y * Section.z) * i + Section.z * j + k;
	const Integer CellID = (BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z) * i + BoundingBox->m_GridSize.z * j + k;
	return CellID;	
}
static __inline__ __device__ __host__ Integer GetCellID(const CGridBox* const BoundingBox, Integer i, Integer j, Integer k) 
{	
	const Integer CellID = (BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z) * i + BoundingBox->m_GridSize.z * j + k;
	return CellID;	
}

static __inline__ __device__ __host__ Integer GetGhostID(const CGridBox* const BoundingBox, const Integer Tid)
{
	Integer FirstHalfGhostCells = (BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z) * GHOST_OFFSET + BoundingBox->m_GridSize.z * GHOST_OFFSET + GHOST_OFFSET;
	if(Tid < FirstHalfGhostCells * CELLMAXPARTICLENUM)
	{
		return Tid;
	}
	Integer ExcludeSecondHalfCells = (BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z) * ( BoundingBox->m_GridSize.x - GHOST_OFFSET ) + BoundingBox->m_GridSize.z * ( BoundingBox->m_GridSize.y - GHOST_OFFSET ) + ( BoundingBox->m_GridSize.z - GHOST_OFFSET );
	return (Tid - ExcludeSecondHalfCells * CELLMAXPARTICLENUM);
}
static __inline__ 
#ifdef __CUDACC__ 
__device__ 
#else
__host__
#endif  
void GetCell(const CGridBox* const BoundingBox, const Integer Tid, Integer& Cid, Integer& Cidx, Integer& Cidy, Integer& Cidz, Integer& Pid) 
{
	Cid = Tid / CELLMAXPARTICLENUM ;
	Pid = Tid % CELLMAXPARTICLENUM;
	Integer YZ = BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z;
	Cidx = Cid / YZ ;
	Integer Rem = Cid % YZ;
	Cidy = Rem / BoundingBox->m_GridSize.z;	
	Cidz = Rem % BoundingBox->m_GridSize.z;
}

static __inline__ 
#ifdef __CUDACC__ 
__device__ 
#else
__host__
#endif  
void GetComputeCell(const CGridBox* const BoundingBox, const Integer Tid, Integer& Cid, Integer& Cidx, Integer& Cidy, Integer& Cidz, Integer& Pid) 
{
	Integer OfsetCell = (BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z) * OFFSET + BoundingBox->m_GridSize.z * OFFSET + OFFSET;
	Cid = Tid / CELLMAXPARTICLENUM + OfsetCell;
	Pid = Tid % CELLMAXPARTICLENUM;
	Integer YZ = BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z;
	Cidx = Cid / YZ ;
	Integer Rem = Cid % YZ;
	Cidy = Rem / BoundingBox->m_GridSize.z;	
	Cidz = Rem % BoundingBox->m_GridSize.z;
}
static __inline__ 
#ifdef __CUDACC__ 
__device__ 
#else
__host__
#endif  
void GetPressureCell(const CGridBox* const BoundingBox, const Integer Tid, Integer& Cid, Integer& Cidx, Integer& Cidy, Integer& Cidz, Integer& Pid) 
{
	Integer OfsetCell = (BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z) * PRESSURE_OFFSET + BoundingBox->m_GridSize.z * PRESSURE_OFFSET + PRESSURE_OFFSET;
	Cid = Tid / CELLMAXPARTICLENUM + OfsetCell;
	Pid = Tid % CELLMAXPARTICLENUM;
	Integer YZ = BoundingBox->m_GridSize.y * BoundingBox->m_GridSize.z;
	Cidx = Cid / YZ ;
	Integer Rem = Cid % YZ;
	Cidy = Rem / BoundingBox->m_GridSize.z;	
	Cidz = Rem % BoundingBox->m_GridSize.z;
}

//static __inline__ __device__ __host__ Scalar3 AddVector(const Scalar3* const P1, const Scalar3* const P2)
//{
//	Scalar3 tmp;
//	tmp.x = P1->x + P2->x;
//	tmp.y = P1->y + P2->y;
//	tmp.z = P1->z + P2->z;
//	return tmp;
//}
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
#endif
static __inline__ __host__ __device__ Scalar3 GetCellPosition(const CGridBox* const BoundingBox, Integer i, Integer j, Integer k)
{	
	Scalar3 V = make_Scalar3(0.0,0.0,0.0);
	V.x = i * BoundingBox->m_CellSize + BoundingBox->m_BufferedZone.m_MinBound.x;
	V.y = j * BoundingBox->m_CellSize + BoundingBox->m_BufferedZone.m_MinBound.y;
	V.z = k * BoundingBox->m_CellSize + BoundingBox->m_BufferedZone.m_MinBound.z;
	return V;
}
static __inline__ __host__ __device__ bool IsInclude(const CBox* const BoundingBox, const Scalar3 * const Position,const Scalar Tolerance)
{	
	if(Position->x < BoundingBox->m_MinBound.x - Tolerance)
	{
		return false;
	}
	if(Position->y < BoundingBox->m_MinBound.y - Tolerance)
	{
		return false;
	}
	if(Position->z < BoundingBox->m_MinBound.z - Tolerance)
	{
		return false;
	}

	if(Position->x >= BoundingBox->m_MaxBound.x + Tolerance) 
	{
		return false;
	}
	if(Position->y >= BoundingBox->m_MaxBound.y + Tolerance)
	{
		return false;
	}
	if(Position->z >= BoundingBox->m_MaxBound.z + Tolerance)
	{
		return false;
	}
	return true;
}

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
	const Scalar Constant = Parameter->InitialParticleDistance * Parameter->CourantNumber / Parameter->Dt;
	Scalar Ratio = static_cast<Scalar> (1.0 / 3.0);
	Scalar NumericConst = sqrt(Ratio);
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
	}		
}

static __inline__ CCTStatusType CudaSafeCall(cudaError_t Status)
{
	if(cudaSuccess != Status)
	{
		printf(cudaGetErrorString(Status));
		/*printf("\t in Function ");
		printf(Function);*/
		return CCT_CUDAERR;	
	}
	return CCT_NOERR;
}


