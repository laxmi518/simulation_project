#include "DataType.h"
#include "Constant.h"
#include "Utility_Inline.h"
static __inline__ __device__ Integer CudaGetTargetID()
{
	return blockDim.x * blockIdx.x + threadIdx.x;
}
#define OUTPUT_MACROINDEX()  int originalIndex = c_daParticleID[ID];

__global__ void CheckParticleOutsideComputeZone_Kernel(Integer * ParticleNum)
{
	int index  = CudaGetTargetID();
	if (index >= *ParticleNum)
	{
		return;
	}
	Integer PartNum = *ParticleNum;
	Scalar3 TargetPosition = c_daOutputParticlePosition[index];
	if(!IsInclude(&CONSTANT_BOUNDINGBOX,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
		c_daOutputParticlePosition[index]			= c_daOutputParticlePosition[PartNum];
		c_daOutputParticleVelocity[index]			= c_daOutputParticleVelocity[PartNum];

		c_daOutputParticlePressure[index]			= c_daOutputParticlePressure[PartNum];
		c_daOutputParticleDensity[index]			= c_daOutputParticleDensity[PartNum];
		c_daOutputParticleTemperature[index]		= c_daOutputParticleTemperature[PartNum];
		c_daOutputParticleKineticViscosity[index]	= c_daOutputParticleKineticViscosity[PartNum];
		c_daOutputParticleSolidPhaseRate[index]		= c_daOutputParticleSolidPhaseRate[PartNum];
		c_daOutputParticleType[index]				= c_daOutputParticleType[PartNum];
		c_daOutputParticleID[index]					= c_daOutputParticleID[PartNum];
		*ParticleNum--;
	}
	else
	{
		return;
	}
}
//////////////////////////////Radix sort implementation///////////////////////////////////
__global__ void calcHashD(Integer ParticleNum,int* dGridParticleHash, int* dGridParticleIndex, Scalar3* particlePosition)
{
	int index  = CudaGetTargetID();
	if (index >= ParticleNum) 
	{
		return;
	}
	Integer hash = GetCellID(&CONSTANT_BOUNDINGBOX,&particlePosition[index]);	

	dGridParticleHash[index]  = hash;
	dGridParticleIndex[index] = index;
}

// rearrange particle data into sorted order, and find the start of each cell
// in the sorted hash array
__global__ void reorderDataAndFindCellStartD( Integer	   numParticles,
											  Integer	   CellNum,
											  Integer	 * gridParticleHash, // input: sorted grid hashes
											  Integer	 * gridParticleIndex,// input: sorted particle indices
											  Integer	 * cellStart,        // output: cell start index
											  Integer	 * cellEnd          // output: cell end index
											  )
{
	extern __shared__ int sharedHash[];    // blockSize + 1 elements
	Integer index = CudaGetTargetID();
	if (index >= numParticles) 
	{
		return;
	}
	int hash	=	0;
	int phash	=	-1;

	// handle case when no. of particles not multiple of block size
	if (index < numParticles)
	{
		hash = gridParticleHash[index];	
		phash = gridParticleHash[index - 1];

		// Load hash data into shared memory so that we can look
		// at neighboring particle's hash value without loading
		// two hash values per thread
		sharedHash[threadIdx.x+1] = hash;

		if (index > 0 && threadIdx.x == 0)
		{
			// first thread in block must load neighbor particle hash
			sharedHash[0] = gridParticleHash[index-1];
		}
	}
	__syncthreads();
	if (index < numParticles)
	{
		// If this particle has a different cell index to the previous
		// particle then it must be the first particle in the cell,
		// so store the index of this particle in the cell.
		// As it isn't the first particle, it must also be the cell end of
		// the previous particle's cell

		if (index == 0 || hash != sharedHash[threadIdx.x])
		{
			cellStart[hash] = index;

			if (index > 0)
				cellEnd[sharedHash[threadIdx.x]] = index;
		}

		if (index == numParticles - 1)
		{
			cellEnd[hash] = index + 1;
		}
		// Now use the sorted index to reorder the pos and vel data
		int sortedIndex		= gridParticleIndex[index];

		c_daParticlePosition[index]			= c_daOutputParticlePosition[sortedIndex];
		c_daParticleVelocity[index]			= c_daOutputParticleVelocity[sortedIndex];

		c_daParticlePressure[index]			= c_daOutputParticlePressure[sortedIndex];
		c_daParticleDensity[index]			= c_daOutputParticleDensity[sortedIndex];
		c_daParticleTemperature[index]		= c_daOutputParticleTemperature[sortedIndex];
		c_daParticleKineticViscosity[index] = c_daOutputParticleKineticViscosity[sortedIndex];
		c_daParticleSolidPhaseRate[index]	= c_daOutputParticleSolidPhaseRate[sortedIndex];
		c_daParticleType[index]				= c_daOutputParticleType[sortedIndex];
		c_daParticleID[index]				= sortedIndex;
	}

}
static __inline__ __host__ __device__ void CalcNearestDistance(CDistance& NearestDistance,Integer& NearestID, Integer ID)
{
	NearestDistance.Magnitude = 1e8;
	NearestDistance.Direction = make_Scalar3(0,0,0);
	if(c_daOutputParticleID[ID] < -1)
	{
		c_daSTLDistance[ID] = NearestDistance;
		c_daSTLID[ID] = NearestID;
		return;
	}
	Scalar3 TargetPosition = c_daParticlePosition[ID];

	Integer TriangleID;		
	Integer CIDX, CIDY, CIDZ;

	Integer CID = GetCellID(&CONSTANT_BOUNDINGBOX,&TargetPosition,CIDX, CIDY, CIDZ);
	if(CID >= 0 && CID < c_CellNum)
	{
		Integer Range = 1;
		for(Integer k = -Range; k <= Range; ++k)
		{
			for(Integer j = -Range; j <= Range; ++j)
			{
				for(Integer i = -Range; i <= Range; ++i)
				{
					Integer MCID = GetCellID(&CONSTANT_BOUNDINGBOX,CIDX +i, CIDY + j,CIDZ + k);
					if(MCID < 0 || MCID >= c_CellNum)
					{
						continue;
					}
					unsigned int TriangleNum = c_daCell[MCID].m_TriangleNum;
					for(unsigned int l = 0; l < TriangleNum; ++l)
					{
						TriangleID = c_daCell[MCID].m_TriangleID[l];

						if(!c_daTrianglesParameters[c_daTriangles[TriangleID].ModelIDNumber].isDrag)
						{
							continue;
						}
						if( TriangleID >= 0 && TriangleID < c_TriangleNum && TriangleID != NearestID)// No need to calculate again for the same triangle
						{
							CDistance Distance ;
							Distance.Magnitude = CalcDistance(&c_daTriangles[TriangleID], &TargetPosition, &Distance.Direction);
							if(Distance.Magnitude < NearestDistance.Magnitude)
							{
								NearestDistance = Distance;
								NearestID = TriangleID;
							}
						}
					}	
				}
			}
		}
	}

}

static __inline__ __host__ __device__ void  CalcDragEffect_ExpoVelocityBoth_Kernel(DragParameter& mDragParam,CDistance& NearestDistance,Integer& NearestID, Integer ID)
{


	c_daDragAcc[ID]	=	make_Scalar3(0.0,0.0,0.0);
	Scalar pressure	=	0.0; //For Velcoity
	Scalar Enthalpy = 0.0; //For Temperature
	Scalar ScalarVelocity	=	0.0; //For Exponential Velocity
	Scalar ScalarTemperature = 0.0; //For Exponential Temperature
	Scalar3 normal = make_Scalar3(0.0,0.0,0.0);

	c_daDragAcc[ID] = make_Scalar3(0,0,0);
	c_daDragTemperature[ID] = 0;



	c_daOutputParticleType[ID] = TYPE_LIGHT_WEIGHT;

	ScalarVelocity	= Magnitude(c_daParticleVelocity[ID]);
	pressure += mDragParam.VelocityIndex * pow(ScalarVelocity,mDragParam.VelocityExponential);

	normal.x	= c_daTriangles[NearestID].Normal.x;
	normal.y	= c_daTriangles[NearestID].Normal.y;
	normal.z	= c_daTriangles[NearestID].Normal.z;

	if(mDragParam.VelocityDragPermiabilityConstant != 0 )
	{
		c_daDragAcc[ID].x = (pressure *  (normal.x) / mDragParam.VelocityDragPermiabilityConstant);
		c_daDragAcc[ID].y = (pressure *  (normal.y) / mDragParam.VelocityDragPermiabilityConstant);
		c_daDragAcc[ID].z = (pressure *  (normal.z) / mDragParam.VelocityDragPermiabilityConstant);
	}	
}


static __inline__ __host__ __device__ void  CalcDragEffect_ExpoTemperatureBoth_Kernel(DragParameter& mDragParam,CDistance& NearestDistance,Integer& NearestID, Integer ID)
{
	c_daDragAcc[ID]	=	make_Scalar3(0.0,0.0,0.0);
	Scalar pressure	=	0.0; //For Velcoity
	Scalar Enthalpy = 0.0; //For Temperature
	Scalar ScalarVelocity	=	0.0; //For Exponential Velocity
	Scalar ScalarTemperature = 0.0; //For Exponential Temperature
	Scalar3 normal = make_Scalar3(0.0,0.0,0.0);

	c_daDragAcc[ID] = make_Scalar3(0,0,0);
	c_daDragTemperature[ID] = 0;
	c_daOutputParticleType[ID] = TYPE_LIGHT_WEIGHT;

	//Formula
	//Heat Quantity = Index * Temperature ^ Exp
	ScalarVelocity	= Magnitude(c_daParticleVelocity[ID]);
	Enthalpy = mDragParam.TempIndex * pow(ScalarVelocity,mDragParam.TempExponential);

	//Formula
	//Dt = (TempDragPermiabilityConstant * Enthalpy * dt * Thermal Calculation Step)/ (Density * Velocity * Fluid Specific Heat Capacity)
	//Enthalpy Unit = {Watt/(mm^2)} or {Watt/(m^2)}

	Scalar Denominator	= Magnitude(c_daParticleVelocity[ID]) * CONSTANT_PARAMETER.Density * CONSTANT_PARAMETER.FruidSpecificHeat;
	if(Denominator != 0 )
	{						
		c_daDragTemperature[ID] = (Enthalpy * mDragParam.TempDragPermiabilityConstant * CONSTANT_PARAMETER.Dt * CONSTANT_PARAMETER.ThermalStep) / (Denominator);						
	}
	c_daOutputParticleSolidPhaseRate[ID] = c_daDragTemperature[ID];

}

static __inline__ __host__ __device__ void  CalcDragEffect_ConstantFunctionVelocity_Kernel(DragParameter& mDragParam,CDistance& NearestDistance,Integer& NearestID, Integer ID)
{
	c_daDragAcc[ID]	=	make_Scalar3(0.0,0.0,0.0);
	Scalar pressure	=	0.0; //For Velcoity
	Scalar Enthalpy = 0.0; //For Temperature
	Scalar ScalarVelocity	=	0.0; //For Exponential Velocity
	Scalar ScalarTemperature = 0.0; //For Exponential Temperature
	Scalar3 normal = make_Scalar3(0.0,0.0,0.0);

	c_daDragAcc[ID] = make_Scalar3(0,0,0);
	c_daDragTemperature[ID] = 0;

	c_daOutputParticleType[ID] = TYPE_LIGHT_WEIGHT;

	c_daOutputParticleVelocity[ID].x = mDragParam.ConstantVx;
	c_daOutputParticleVelocity[ID].y = mDragParam.ConstantVy;
	c_daOutputParticleVelocity[ID].z = mDragParam.ConstantVz;

	c_daParticleVelocity[ID].x = mDragParam.ConstantVx;
	c_daParticleVelocity[ID].y = mDragParam.ConstantVy;
	c_daParticleVelocity[ID].z = mDragParam.ConstantVz;


}


static __inline__ __host__ __device__ void  CalcDragEffect_FunctionVelocity_Kernel(DragParameter& mDragParam,CDistance& NearestDistance,Integer& NearestID, Integer ID)
{
	c_daDragAcc[ID]	=	make_Scalar3(0.0,0.0,0.0);
	Scalar pressure	=	0.0; //For Velcoity
	Scalar Enthalpy = 0.0; //For Temperature
	Scalar ScalarVelocity	=	0.0; //For Exponential Velocity
	Scalar ScalarTemperature = 0.0; //For Exponential Temperature
	Scalar3 normal = make_Scalar3(0.0,0.0,0.0);

	c_daDragAcc[ID] = make_Scalar3(0,0,0);
	c_daDragTemperature[ID] = 0;
	c_daOutputParticleType[ID] = TYPE_LIGHT_WEIGHT;

	c_daParticleVelocity[ID].x *=  mDragParam.VelocityMagnitueFactor;
	c_daParticleVelocity[ID].y *=  mDragParam.VelocityMagnitueFactor;
	c_daParticleVelocity[ID].z *=  mDragParam.VelocityMagnitueFactor;

	c_daOutputParticleVelocity[ID].x *=  mDragParam.VelocityMagnitueFactor;
	c_daOutputParticleVelocity[ID].y *=  mDragParam.VelocityMagnitueFactor;
	c_daOutputParticleVelocity[ID].z *=  mDragParam.VelocityMagnitueFactor;

}














__global__ static void CalcDragEffect_Kernel(Integer ComputeParticleNumber)
{
	Integer TID = CudaGetTargetID();
	Integer ID  = TID;
	Integer OID = TID;
	if(ID >= ComputeParticleNumber)
	{
		return ;
	}
	OUTPUT_MACROINDEX();
	//c_daOutputParticleSolidPhaseRate[originalIndex] = 0;	//Test
	
	if(c_daMagnifierCount[ID]>0)
	{
		c_daMagnifierCount[ID]++;
	}
	if(c_daMagnifierCount[ID]==10)
	{
		c_daMagnifierCount[ID]=0;
	}
	CDistance NearestDistance;
	Integer NearestID = -1;

	CalcNearestDistance(NearestDistance,NearestID,ID);

	c_daDragAcc[ID]	= make_Scalar3(0.0,0.0,0.0);
	Scalar pressure	= 0.0; //For Velcoity
	Scalar Enthalpy = 0.0; //For Temperature
	Scalar ScalarVelocity	=	0.0; //For Exponential Velocity
	Scalar ScalarTemperature = 0.0; //For Exponential Temperature
	Scalar3 normal = make_Scalar3(0.0,0.0,0.0);

	c_daDragAcc[ID] = make_Scalar3(0,0,0);
	c_daDragTemperature[ID] = 0;
	c_daOutputParticleType[ID] = TYPE_HEAVY_WEIGHT;	
	
	if(NearestID != -1 && NearestDistance.Magnitude != 1e8)
	{
		DragParameter mDragParam = c_daSTLDragParameter[c_daTrianglesParameters[c_daTriangles[NearestID].ModelIDNumber].DragTriangleID];	

		if( NearestDistance.Magnitude <= mDragParam.CoefficientofInflueneRegion * CONSTANT_PARAMETER.InitialParticleDistance)
		{
			c_daOutputParticleType[ID] = TYPE_LIGHT_WEIGHT;
			if(mDragParam.DragToken == 'E') // E is used to indicate its Exponential Function
			{
				if(mDragParam.InputVelocityToken1	==	'V' ||	mDragParam.InputVelocityToken1	==	'B')  //V For Velocity and B for Both
				{				
						CalcDragEffect_ExpoVelocityBoth_Kernel(mDragParam,NearestDistance,NearestID,ID);
	
				}
				else if(mDragParam.InputVelocityToken1	==	'T' ||	mDragParam.InputVelocityToken1	==	'B') //T for temerature and B for Both
				{
					CalcDragEffect_ExpoTemperatureBoth_Kernel(mDragParam,NearestDistance,NearestID,ID);
				}
				else
				{
					return;
				}				
			}
			else if(mDragParam.DragToken == 'C') // C is used to Indicate its Constant Function
			{
				if(mDragParam.InputVelocityToken1	==	'V' ||	mDragParam.InputVelocityToken1	==	'v')
				{
						
					CalcDragEffect_ConstantFunctionVelocity_Kernel(mDragParam,NearestDistance,NearestID,ID);
				}
				}
			}
			else if(mDragParam.DragToken == 'F') // C is used to Indicate its Constant Function
			{
				if(mDragParam.InputVelocityToken1	==	'V' ||	mDragParam.InputVelocityToken1	==	'v')
				{
				
					if(c_daMagnifierCount[ID] < 1)
					{
					CalcDragEffect_FunctionVelocity_Kernel(mDragParam,NearestDistance,NearestID,ID);
					}
				}
			}
		}		
	}	
__global__ static void CalcExplicitly_Kernel(Integer ComputeParticleNumber)
{
	Integer TID = CudaGetTargetID();
	Integer ID  = TID;
	Integer OID = TID;
	if(ID >= ComputeParticleNumber)
	{
		return ;
	}
	OUTPUT_MACROINDEX();
	//c_daOutputParticleSolidPhaseRate[originalIndex] = 0;	//Test

	Scalar SpecificDensity = CONSTANT_PARAMETER.Density;
	Scalar3 TargetPosition = c_daParticlePosition[ID];
	
	if(!IsInclude(&CONSTANT_BOUNDINGBOX,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
		c_daOutputParticleID[originalIndex]		= -1;
		c_daOutputParticleSolidPhaseRate[originalIndex] = 9999;
		return;
	}
	Integer CIDX, CIDY, CIDZ;
	Integer CID = GetCellID(&CONSTANT_BOUNDINGBOX,&TargetPosition,CIDX, CIDY, CIDZ);
	Scalar3 Laplacian	  = make_Scalar3(0,0,0);
	Scalar3 ArtificialAcc = make_Scalar3(0,0,0);
	Scalar3 CollisionAcc  = make_Scalar3(0,0,0);
	Scalar3 V0 = c_daParticleVelocity[ID];
	
	Scalar3 TurbulentKineticEnergy = make_Scalar3(0,0,0);

	//if(CID < 0 || CID >= c_CellNum)
	//{
	//	c_daOutputParticleSolidPhaseRate[originalIndex] = 9998;
	//	return;
	//}

	
	const Scalar Re = CONSTANT_PARAMETER.LaplacianInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	const Scalar Rc = CONSTANT_PARAMETER.r_coeff * CONSTANT_PARAMETER.InitialParticleDistance;
	const Scalar ReGrad = CONSTANT_PARAMETER.GradientInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;

	Scalar Density = 0.0;
	if(CID >= 0 && CID < c_CellNum)
	{
		const Scalar MinDist = Re * 1e-3;
		Integer Range = 1;
		Integer Neighbour = 0;
		
		for(Integer k = -Range; k <= Range; ++k)
		{
			for(Integer j = -Range; j <= Range; ++j)
			{
				for(Integer i = -Range; i <= Range; ++i)
				{
					Integer MCID = GetCellID(&CONSTANT_BOUNDINGBOX,CIDX +i, CIDY + j,CIDZ + k);
					if(MCID < 0 || MCID >= c_CellNum)
					{
						continue;
					}
					int startIndex = c_cellStart[MCID];
					if (startIndex != 0xffffffff)    // cell is not empty
					{       
						// iterate over particles in this cell
						int endIndex = c_cellEnd[MCID];
						for(int x = startIndex; x < endIndex; x++)
						{
							Integer NeighborID = x;
							if(c_daOutputParticleID[NeighborID] < 0)	//Test
							{
								continue;
							}
							if(NeighborID != ID)
							{				
								Scalar Distance = CalcDistance(&TargetPosition,&c_daParticlePosition[NeighborID]);
								if(Distance < Re)
								{
									if(Distance < MinDist)
									{
										Distance = MinDist;
									}
									Scalar Weight = GetWeight(Re, Distance);
									Density += Weight;
									Scalar3 V1 = c_daParticleVelocity[NeighborID];
									Laplacian.x +=  (V1.x - V0.x) * Weight;
									Laplacian.y +=  (V1.y - V0.y) * Weight;
									Laplacian.z +=  (V1.z - V0.z) * Weight;
									//TurbulenceKineticEnergy
									if(CONSTANT_PARAMETER.bTurbulance)
									{
										if(Distance < ReGrad)
										{
											Scalar Weight1 = GetWeight(ReGrad, Distance);
											Scalar StrainCoefficient = (c_daParticleStrainTensorProduct[NeighborID] - c_daParticleStrainTensorProduct[ID]) * Weight1 / (Distance * Distance);
											TurbulentKineticEnergy = TurbulentKineticEnergy + (c_daParticlePosition[NeighborID] - TargetPosition) * StrainCoefficient;							
										}
									}
									if(CONSTANT_PARAMETER.bExplicitPressure)
									{									
										if(Distance < Rc)
										{	
											Scalar Art = CONSTANT_PARAMETER.ArtificialPressure * Weight / (Distance * Distance);		
											ArtificialAcc = ArtificialAcc + (c_daParticlePosition[NeighborID] - TargetPosition) * Art;
											Scalar Col = (Rc - Distance) / Distance;
											CollisionAcc = CollisionAcc + (c_daParticlePosition[NeighborID] - TargetPosition) * Col;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	Scalar Distance = c_daSTLDistance[ID].Magnitude;
	Integer TriangleID = c_daSTLID[ID];
	Scalar WallWeight = 0.0;
	if(TriangleID >= 0 && TriangleID < c_TriangleNum)
	{
		if(Distance >0 && Distance < Re)
		{
			//Scalar WallWeight = 0.0;
			if(2 == CONSTANT_PARAMETER.Dimension)
			{
				WallWeight = GetZValueFunctionLaplacian2D(Distance, CONSTANT_PARAMETER.InitialParticleDistance);
			}
			else
			{
				WallWeight = GetZValueFunctionLaplacian(Distance, CONSTANT_PARAMETER.InitialParticleDistance);
			}
			Density += WallWeight;
			Scalar3 V1 = c_daTrianglesParameters[c_daTriangles[TriangleID].ModelIDNumber].Velocity;
			
			Laplacian.x += CONSTANT_PARAMETER.WallFrictionCoeffieicnt * (V1.x - V0.x) * WallWeight;
			Laplacian.y += CONSTANT_PARAMETER.WallFrictionCoeffieicnt * (V1.y - V0.y) * WallWeight;
			Laplacian.z += CONSTANT_PARAMETER.WallFrictionCoeffieicnt * (V1.z - V0.z) * WallWeight;				
	
		}
	}
	Scalar Coefficient = ( 2.0 * CONSTANT_PARAMETER.Dimension ) / ( CONSTANT_PARAMETER.LambdaValue * CONSTANT_PARAMETER.ConstantParticleDensity);
	Laplacian.x *= Coefficient;
	Laplacian.y *= Coefficient;
	Laplacian.z *= Coefficient;	

	Scalar3 DragVel = make_Scalar3(0.0, 0.0, 0.0);
	
	// Velocity Update	
	Scalar IA = 0;
	Scalar3 Acc = make_Scalar3(0,0,IA);
	Acc = Acc + CONSTANT_PARAMETER.GravityConstant;
	if(c_DragTriangleNum > 0)
	{
		//Acc = Acc + c_daDragAcc[ID];
		DragVel = c_daDragAcc[ID];
	}
	if(CONSTANT_PARAMETER.bTurbulance)
	{
		Laplacian = Laplacian * (c_daParticleKineticViscosity[ID] + c_daParticleTurbulaceViscosity[ID]); //Addition of turbulence viscosity to laminar viscosity
	}
	else
	{
		Laplacian = Laplacian * c_daParticleKineticViscosity[ID] ;
	}
	
	Acc = Acc + Laplacian;
	V0 = V0 + Acc * CONSTANT_PARAMETER.Dt - DragVel;
	
	//Velocity update due to TurbulenceKineticEnergy
	if (CONSTANT_PARAMETER.bTurbulance)
	{
		TurbulentKineticEnergy = TurbulentKineticEnergy * (2/3); //This step calculates del(2k/3)/del.r 		
		V0 = V0 - TurbulentKineticEnergy * CONSTANT_PARAMETER.Dt;
	}

	if(CONSTANT_PARAMETER.bExplicitPressure)
	{
		Scalar coeff = - CONSTANT_PARAMETER.Dimension / (CONSTANT_PARAMETER.ConstantParticleDensity * SpecificDensity);
		//Scalar3 Acc = make_Scalar3(0,0,0);
		Acc = ArtificialAcc * coeff + CollisionAcc * 0.01;	
		V0 = V0 + Acc * CONSTANT_PARAMETER.Dt;
	}
	if(isfinite(V0.x) && isfinite(V0.y) && isfinite(V0.z))
	{
		c_daOutputParticlePosition[originalIndex].x = TargetPosition.x + CONSTANT_PARAMETER.Dt * V0.x;
		c_daOutputParticlePosition[originalIndex].y = TargetPosition.y + CONSTANT_PARAMETER.Dt * V0.y;
		c_daOutputParticlePosition[originalIndex].z = TargetPosition.z + CONSTANT_PARAMETER.Dt * V0.z;

		c_daOutputParticleVelocity[originalIndex].x = V0.x;
		c_daOutputParticleVelocity[originalIndex].y = V0.y;
		c_daOutputParticleVelocity[originalIndex].z = V0.z;
	}
	
}

__global__ static void  RegisterTriangleTopology_Kernel(const CTriangle* const daTriangle, Integer TriangleNum, CCell* daCell, Integer CellNum)
{
	Integer ID = CudaGetTargetID();
	if(ID >= TriangleNum)
	{
		return ;
	}
	RegisterTriangle(&CONSTANT_BOUNDINGBOX,&daTriangle[ID],ID,daCell, CellNum);
	//RegisterTriangle(&CONSTANT_BOUNDINGBOX,&CONSTANT_GRIDPARAMS,daTriangle[ID].Vertex[0],daTriangle[ID].Vertex[1],daTriangle[ID].Vertex[2],ID,daCell, CellNum);
}
__global__ static void  ResetTriangleTopology_Kernel(Integer CellNum, CCell* aCell)
{
	Integer ID = CudaGetTargetID();
	if(ID >= CellNum)
	{
		return;
	}
	c_daCell[ID].m_TriangleNum = 0;
}
__global__ static void  UpdateTrianglePosition_Kernel(const Integer TriangleNum, CTriangle* daTriangles)
{
	Integer ID = CudaGetTargetID();
	if(ID >= TriangleNum)
	{
		return ;
	}
	Scalar3 Vel = c_daTrianglesParameters[c_daTriangles[ID].ModelIDNumber].Velocity;

	Scalar dt = CONSTANT_PARAMETER.Dt;

	Vel.x = Vel.x * dt;
	Vel.y = Vel.y * dt;
	Vel.z = Vel.z * dt;

	c_daTriangles[ID].Vertex[0] = c_daTriangles[ID].Vertex[0] + Vel;
	c_daTriangles[ID].Vertex[1] = c_daTriangles[ID].Vertex[1] + Vel;
	c_daTriangles[ID].Vertex[2] = c_daTriangles[ID].Vertex[2] + Vel;
}
__global__ static void  RotateTrianglePosition_Kernel(const Integer TriangleNum,CTriangle* daTriangles, const Integer analysisStep)			
{			
	Integer ID = CudaGetTargetID();		
	if(ID >= TriangleNum)		
	{		
		return ;	
	}
	CTriangle* TargetTriangle = &c_daTriangles[ID];

	//Set Triangle Parameter @Rajan 20121004
	CTriangleParameters TriangleParameters;
	TriangleParameters = c_daTrianglesParameters[TargetTriangle->ModelIDNumber];

	if(!TriangleParameters.RotatingOrNot)
	{
		return;
	}
	//To make Slope of the Rotation Starts

	if((analysisStep * CONSTANT_PARAMETER.Dt) < TriangleParameters.RotationStartTime)
	{
		return;
	}

	Scalar angleOfRotation = 0.0;
	int requiredStep = (TriangleParameters.TimeToReachRPM - TriangleParameters.RotationStartTime) / CONSTANT_PARAMETER.Dt;
	Scalar timePassed = CONSTANT_PARAMETER.Dt * analysisStep;
	if((timePassed < TriangleParameters.RotationStartTime) || (TriangleParameters.AngleOfRotationInDegree == 0))
	{
		angleOfRotation = 0.0;
	}
	else if((timePassed >= TriangleParameters.RotationStartTime) && (timePassed < TriangleParameters.TimeToReachRPM))
	{
		angleOfRotation = (TriangleParameters.AngleOfRotationInDegree / (requiredStep + 1)) * (analysisStep - (TriangleParameters.RotationStartTime/CONSTANT_PARAMETER.Dt));
	}
	else
	{
		angleOfRotation = TriangleParameters.AngleOfRotationInDegree;
	}
	//To make Slope of the Rotation Ends
	if(angleOfRotation == 0)
		return;
	

	for(int i = 0 ; i < 3 ; ++i)		
	{
		if(angleOfRotation != 0)
		{
			Scalar3 normalVector = CalcNormal(&TriangleParameters.CenterOfRotation, &TriangleParameters.SecondPointOfRotation); //Calculate SecondPoint - FirsPoint		
			Scalar theta = static_cast<Scalar>((angleOfRotation * M_PI) / 180.0);		
			Scalar a = TriangleParameters.CenterOfRotation.x;
			Scalar b = TriangleParameters.CenterOfRotation.y;
			Scalar c = TriangleParameters.CenterOfRotation.z;
			
			Scalar u = normalVector.x;
			Scalar v = normalVector.y;
			Scalar w = normalVector.z;

			Scalar x = TargetTriangle->Vertex[i].x;
			Scalar y = TargetTriangle->Vertex[i].y;
			Scalar z = TargetTriangle->Vertex[i].z;
				
			TargetTriangle->Vertex[i].x = (a * ( v*v + w*w) - (u * (b*v + c*w - u*x - v*y - w*z))) * (1 - cos(theta)) + (x * cos(theta)) + (((-c*v) + b*w - w*y + v*z) * sin(theta));	
			TargetTriangle->Vertex[i].y = (b * ( u*u + w*w) - (v * (a*u + c*w - u*x - v*y - w*z))) * (1 - cos(theta)) + (y * cos(theta)) + ((( c*u) - a*w + w*x - u*z) * sin(theta));	
			TargetTriangle->Vertex[i].z = (c * ( u*u + v*v) - (w * (a*u + b*v - u*x - v*y - w*z))) * (1 - cos(theta)) + (z * cos(theta)) + (((-b*u) + a*v - v*x + u*y) * sin(theta));	
		}
	}		
}
__global__ static void  ResetWallPosition_Kernel(const Integer TriangleNum, const Integer AnalysisStep, const CTriangle* daOriginalTriangles)
{
	Integer ID = CudaGetTargetID();
	if(ID >= TriangleNum)
	{
		return ;
	}	
	if(c_daTrianglesParameters[c_daTriangles[ID].ModelIDNumber].LinearMovableOrNot)
	{
		if(0 ==  ((AnalysisStep+1) % c_daTrianglesParameters[c_daTriangles[ID].ModelIDNumber].LinearResetInterval))
		{
			memcpy(&c_daTriangles[ID], &daOriginalTriangles[ID], sizeof(CTriangle));
		}
	}
}
__global__ static void  CalcTemperatureFactor_Kernel(Integer ComputeParticleNumber)
{
	Integer TID = CudaGetTargetID();
	Integer ID  = TID;
	Integer OID = TID;
	if(ID >= ComputeParticleNumber)
	{
		return ;
	}
	Scalar Density = CONSTANT_PARAMETER.Density;
	const Scalar EnthalpyCoefficient = CONSTANT_PARAMETER.Dt * 2 * CONSTANT_PARAMETER.Dimension / 
		( CONSTANT_PARAMETER.ConstantParticleDensity * Density * CONSTANT_PARAMETER.FruidSpecificHeat);
	Scalar Hf = 0.0, Hw = 0.0;
	const Scalar ReTemperature = CONSTANT_PARAMETER.TemperatureInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	
	Scalar  numerator = 0, denominator = 0;
	Scalar3 TargetPosition    = c_daParticlePosition[ID];
	Scalar  TargetTemperature = c_daParticleTemperature[ID];

	OUTPUT_MACROINDEX();

	if(!IsInclude(&CONSTANT_BOUNDINGBOX,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
		c_daOutputParticleID[originalIndex]		= -1;
		return;
	}

	Scalar ThermalResistance   = CONSTANT_PARAMETER.FruidThermalResistance;
	Scalar ThermalConductivity = CONSTANT_PARAMETER.FruidThermalConductivity;	

	Integer CIDX, CIDY, CIDZ;
	Integer CID = GetCellID(&CONSTANT_BOUNDINGBOX,&TargetPosition,CIDX, CIDY, CIDZ);
	if(CID < 0 || CID >= c_CellNum)
	{
		return;
	}
	const Scalar MinDist = ReTemperature * 1e-3;
	//const Scalar Re      = CONSTANT_PARAMETER.TemperatureInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	//const Scalar MinDist = Re * 1e-3;

	Integer Range = 1;
	Integer Neighbour = 0;
	for(Integer k = -Range; k <= Range; ++k)
	{
		for(Integer j = -Range; j <= Range; ++j)
		{
			for(Integer i = -Range; i <= Range; ++i)
			{
				Integer MCID = GetCellID(&CONSTANT_BOUNDINGBOX,CIDX +i, CIDY + j,CIDZ + k);
				if(MCID < 0 || MCID >= c_CellNum)
				{
					continue;
				}
				int startIndex = c_cellStart[MCID];
				if (startIndex != 0xffffffff)    // cell is not empty
				{       
					// iterate over particles in this cell
					int endIndex = c_cellEnd[MCID];
					for(int x = startIndex; x < endIndex; x++)
					{
						Integer NeighborID = x;
						if(c_daOutputParticleID[NeighborID] < 0)
						{
							continue;
						}
						if(NeighborID != ID)
						{				
							Scalar Distance = CalcDistance(&TargetPosition,&c_daParticlePosition[NeighborID]);
							if(Distance <= ReTemperature)
							{
								if(Distance < MinDist)
								{
									Distance = MinDist;
								}
								numerator = (c_daParticleTemperature[NeighborID] - TargetTemperature) * GetWeight(ReTemperature, Distance);
								denominator = Distance * Distance * (1 / CONSTANT_PARAMETER.FruidThermalConductivity + CONSTANT_PARAMETER.FruidThermalResistance / Distance);
								Hf += numerator / denominator;

								Neighbour ++;
							}
						}
					}
				}
			}
		}
	}

	ThermalResistance = CONSTANT_PARAMETER.WallThermalResistance;
	ThermalConductivity = (2 * CONSTANT_PARAMETER.FruidThermalConductivity * CONSTANT_PARAMETER.WallThermalConductivity ) / ( CONSTANT_PARAMETER.FruidThermalConductivity + CONSTANT_PARAMETER.WallThermalConductivity);

	Scalar Distance = c_daSTLDistance[ID].Magnitude;
	
	Integer TriangleID = c_daSTLID[ID];
	denominator = 0.0;
	numerator = 0.0;
	if(TriangleID >= 0 && TriangleID < c_TriangleNum)
	{
		if(Distance > 0 && Distance < ReTemperature)
		{
			ThermalResistance	=	c_daTrianglesParameters[c_daTriangles[TriangleID].ModelIDNumber].WallThermalResistivity;
 
			ThermalConductivity	=	(2 * CONSTANT_PARAMETER.FruidThermalConductivity * c_daTrianglesParameters[c_daTriangles[TriangleID].ModelIDNumber].WallThermalConductivity ) / 
									( CONSTANT_PARAMETER.FruidThermalConductivity + c_daTrianglesParameters[c_daTriangles[TriangleID].ModelIDNumber].WallThermalConductivity);
			Scalar WallWeight= 0.0;
			if(2 == CONSTANT_PARAMETER.Dimension)
			{
				WallWeight = GetZValueFunctionGradient2D(Distance, CONSTANT_PARAMETER.InitialParticleDistance);
			}
			else
			{
				WallWeight = GetZValueFunctionGradient(Distance, CONSTANT_PARAMETER.InitialParticleDistance);
			}
			numerator = (c_daTrianglesParameters[c_daTriangles[TriangleID].ModelIDNumber].Temperature - TargetTemperature) * WallWeight;
			denominator = Distance * Distance * (1 / ThermalConductivity + ThermalResistance / Distance);
			Hw = numerator / denominator;
		}
	}
	Scalar HDrag = 0.0;
	if(c_DragTriangleNum > 0)       
	{
		HDrag = c_daDragTemperature[ID];
		TargetTemperature +=HDrag;
	}
	// Ends Temprature Effect from the Drag model STL

	Scalar HFHWComb = 0.0;
	if(abs(Hf) > SCALAR_EPSILON)
	{
		HFHWComb += Hf;
	}
	if(abs(Hw) > SCALAR_EPSILON)
	{
		HFHWComb += Hw;
	}
	TargetTemperature += (EnthalpyCoefficient * (HFHWComb));
	//Limit the Temperature
	if(TargetTemperature > CONSTANT_PARAMETER.MaxWallTemperature)
	{
		TargetTemperature = CONSTANT_PARAMETER.MaxWallTemperature;
	}
	if(TargetTemperature < CONSTANT_PARAMETER.MinWallTemperature)
	{
		TargetTemperature = CONSTANT_PARAMETER.MinWallTemperature;
	}
	
	c_daOutputParticleTemperature[originalIndex] =  TargetTemperature;
	c_daOutputParticleKineticViscosity[originalIndex] = CONSTANT_PARAMETER.ViscosityUpdateParameterA * TargetTemperature + CONSTANT_PARAMETER.ViscosityUpdateParameterB;
}
__global__ static void CalcSTLDistance_Kernel(Integer ComputeParticleNumber)
{
	//const Integer TID = CudaGetTargetID();
	const Integer ID  =CudaGetTargetID(); 
	/*if(ID >= ComputeParticleNumber)
	{
		return ;
	}*/
	CDistance NearestDistance;
	Integer NearestID = -1;
	NearestDistance.Magnitude = 1e8;
	NearestDistance.Direction.x = 0;
	NearestDistance.Direction.y = 0;
	NearestDistance.Direction.z = 0;//make_Scalar3(0,0,0);
	//if(c_daOutputParticleID[ID] < -1)
	//{
	//	c_daSTLDistance[ID] = NearestDistance;
	//	c_daSTLID[ID] = NearestID;
	//	return;
	//}

	//Scalar3 TargetPosition = c_daParticlePosition[ID];

	Integer TriangleID;		
	Integer CIDX, CIDY, CIDZ;
	Integer CID = GetCellID(&CONSTANT_BOUNDINGBOX,&c_daParticlePosition[ID],CIDX, CIDY, CIDZ);
	if(CID >=0 && CID < c_CellNum)
	{
		//Integer Range = 1;
		for(Integer k = -1; k <= 1; ++k)
		{
			for(Integer j = -1; j <= 1; ++j)
			{
				for(Integer i = -1; i <= 1; ++i)
				{
					Integer MCID = GetCellID(&CONSTANT_BOUNDINGBOX,CIDX +i, CIDY + j,CIDZ + k);
					if(MCID < 0 || MCID >= c_CellNum)
					{
						continue;
					}
					unsigned int TriangleNum = c_daCell[MCID].m_TriangleNum;
					for(unsigned int l = 0; l < TriangleNum; ++l)
					{
						TriangleID = c_daCell[MCID].m_TriangleID[l];
						/*if(c_daTrianglesParameters[c_daTriangles[TriangleID].ModelIDNumber].isDrag)
						{
							continue;
						}*/

						if( TriangleID >= 0 && TriangleID < c_TriangleNum && TriangleID != NearestID)// No need to calculate again for the same triangle
						{
						CDistance Distance ;
							Distance.Magnitude = CalcDistance(&c_daTriangles[TriangleID], &c_daParticlePosition[ID], &Distance.Direction);
							if(Distance.Magnitude < NearestDistance.Magnitude)
							{
								NearestDistance = Distance;
								NearestID = TriangleID;
							}
						}
					}	
				}
			}
		}
	}
	c_daSTLDistance[ID] = NearestDistance;
	c_daSTLID[ID] = NearestID;
}
/////////////Calculate Explicit Pressure Calcuation
__global__ static void CalcExplicitPressure_Kernel(Integer ComputeParticleNumber)
{
	Integer ID = CudaGetTargetID();
	//Integer ID  =  TID;
	//Integer OID =  TID;

	if(ID >= ComputeParticleNumber)
	{
		return ;
	}
	OUTPUT_MACROINDEX();
	Scalar3 TargetPosition = c_daParticlePosition[ID];
	if(!IsInclude(&CONSTANT_BOUNDINGBOX,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
		c_daOutputParticleID[originalIndex]		= -1;
		c_daParticlePressure[ID] = 0.0;
		return;
	}
	Scalar SpecificDensity = CONSTANT_PARAMETER.Density;
	Integer CIDX, CIDY, CIDZ;
	Integer CID = GetCellID(&CONSTANT_BOUNDINGBOX,&TargetPosition,CIDX, CIDY, CIDZ);
	if(CID < 0 || CID >= c_CellNum)
	{
		c_daParticlePressure[ID] = 0.0;
		return;
	}
	const Scalar Re = CONSTANT_PARAMETER.GradientInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	const Scalar MinDist = Re * 1e-3;

	Scalar Coeffieient = CONSTANT_PARAMETER.SonicVelocity * CONSTANT_PARAMETER.SonicVelocity * SpecificDensity / CONSTANT_PARAMETER.ConstantParticleDensity;
	Scalar Density = 0.0;
	Integer Range = 1;
	Integer Neighbour = 0;
	for(Integer k = -Range; k <= Range; ++k)
	{
		for(Integer j = -Range; j <= Range; ++j)
		{
			for(Integer i = -Range; i <= Range; ++i)
			{
				Integer MCID = GetCellID(&CONSTANT_BOUNDINGBOX,CIDX +i, CIDY + j,CIDZ + k);
				if(MCID < 0 || MCID >= c_CellNum)
				{
					continue;
				}
				int startIndex = c_cellStart[MCID];
				if (startIndex != 0xffffffff)    // cell is not empty
				{       
					// iterate over particles in this cell
					int endIndex = c_cellEnd[MCID];
					for(int x = startIndex; x < endIndex; x++)
					{
						Integer NeighborID = x;
						if(c_daOutputParticleID[NeighborID] < 0)
						{
							continue;
						}
						if(NeighborID != ID)
						{				
							Scalar Distance = CalcDistance(&TargetPosition,&c_daParticlePosition[NeighborID]);
							if(Distance < Re)
							{
								if(Distance < MinDist)
								{
									Distance = MinDist;
								}
								Scalar Weight = GetWeight(Re, Distance);
								Density += Weight;

								Neighbour++;
							}
						}
					}
				}
			}
		}
	}
	//old 
	//Scalar Distance = c_daSTLDistance[ID].Magnitude;

	//Integer TriangleID = c_daSTLID[ID];
	//if(TriangleID >= 0 && TriangleID < c_TriangleNum)
	//{
	//	if(Distance >=0 && Distance <= Re)
	//	{
	//		Scalar WallWeight = 0.0;
	//		if(2 == CONSTANT_PARAMETER.Dimension)
	//		{
	//			WallWeight = GetZValueFunctionGradient2D(Distance, CONSTANT_PARAMETER.InitialParticleDistance);
	//		}
	//		else
	//		{
	//			WallWeight = GetZValueFunctionGradient(Distance, CONSTANT_PARAMETER.InitialParticleDistance);
	//		}
	//		Density += WallWeight;			
	//	}
	//}
	// old end
	
	//***********new*************
	Scalar r_iw = c_daSTLDistance[ID].Magnitude;
	Scalar WallWeight = 0.0;
	if(c_daSTLID[ID] >= 0 && c_daSTLID[ID] < c_TriangleNum)
	{
		if(r_iw >=0 && r_iw <= Re)
		{
			if(2 == CONSTANT_PARAMETER.Dimension)
			{
				WallWeight = GetZValueFunctionGradient2D(r_iw, CONSTANT_PARAMETER.InitialParticleDistance);
			}
			else
			{
				WallWeight = GetZValueFunctionGradient(r_iw, CONSTANT_PARAMETER.InitialParticleDistance);
			}
		}
	}
	Density += WallWeight;
	if(isfinite(Density))
	{
		c_daOutputParticleDensity[originalIndex]  = Density ;
	}
	//**************new end*********************

	Scalar Pressure = 0.0;
	if(Density > CONSTANT_PARAMETER.ConstantParticleDensity * CONSTANT_PARAMETER.FreeSurfaceThreshold)
	{
		Pressure = Coeffieient * (Density - CONSTANT_PARAMETER.ConstantParticleDensity);
	}
	//c_daOutputParticleDensity[originalIndex] = Density;	
	if(isfinite(Pressure))
	{
		c_daOutputParticlePressure[originalIndex] = Pressure;
		c_daParticlePressure[ID] = Pressure; //NOTE: - For using in the next step you need to insert the data with index of "ID"
	}
}

__global__ static void CalcExplicitPressureGradient_Kernel(Integer ComputeParticleNumber)
{
	Integer TID =	CudaGetTargetID();
	Integer ID	=	 TID;
	Integer OID =	 TID;
	
	if(ID >= ComputeParticleNumber)
	{
		return ;
	}
	OUTPUT_MACROINDEX();
		
	Scalar SpecificDensity = CONSTANT_PARAMETER.Density;	
	Scalar3 TargetPosition = c_daParticlePosition[ID];

	if(!IsInclude(&CONSTANT_BOUNDINGBOX,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
		c_daOutputParticleID[originalIndex]		= -1;
		return;
	}

	Integer CIDX, CIDY, CIDZ;
	Integer CID = GetCellID(&CONSTANT_BOUNDINGBOX,&TargetPosition,CIDX, CIDY, CIDZ);
	if(CID < 0 || CID >= c_CellNum)
	{
		return;
	}
	const Scalar Re = CONSTANT_PARAMETER.GradientInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	Scalar3 GradientFluid = make_Scalar3(0.0,0.0,0.0);
	Scalar3 GradientWall = make_Scalar3(0.0,0.0,0.0);
	Scalar TargetPressure  = c_daParticlePressure[ID]; //For Test only use this for now..Comment this later.
	Integer Cntr = 0;
	//Variables for Negative Pressure
	Scalar3 GradientNegativePressure = make_Scalar3(0.0,0.0,0.0);
	Scalar3 DirectionVector = make_Scalar3(0,0,0);
	Scalar ImaginaryParticleWeight = 0;
	CDistance NearestDistance;
	NearestDistance.Magnitude = 1e8;
	
	Integer NearestID = -1;
	Integer TriangleID;
	Scalar PCoefficient = -1;
	
	Integer EnableNeighborID = -1;
	Integer NeighborPressureSum = -1;
	if(CID >= 0 && CID < c_CellNum)
	{
		const Scalar MinDist = Re * 1e-3;
		//Scalar TargetPressure = c_daParticlePressure[ID] ;
		Integer Range = 1;
		for(Integer k = -Range; k <= Range; ++k)
		{
			for(Integer j = -Range; j <= Range; ++j)
			{
				for(Integer i = -Range; i <= Range; ++i)
				{
					Integer MCID = GetCellID(&CONSTANT_BOUNDINGBOX,CIDX +i, CIDY + j,CIDZ + k);
					if(MCID < 0 || MCID >= c_CellNum)
					{
						continue;
					}
					//*********for negative pressure************
					unsigned int TriangleNum = c_daCell[MCID].m_TriangleNum;
					for(unsigned int l = 0; l < TriangleNum; ++l)
					{
						TriangleID = c_daCell[MCID].m_TriangleID[l];
						if( TriangleID >= 0 && TriangleID < c_TriangleNum && TriangleID != NearestID)// No need to calculate again for the same triangle
						{
							CDistance Distance ;
							Distance.Magnitude = CalcDistance(&c_daTriangles[TriangleID], &TargetPosition, &Distance.Direction);
							if(l==0)
							{
								NearestDistance = Distance;
								NearestID = TriangleID;
							}
							if(Distance.Magnitude < NearestDistance.Magnitude)
							{
								NearestDistance = Distance;
								NearestID = TriangleID;
							}
						}
					}
					//***********for negative pressure ends**************
					int startIndex = c_cellStart[MCID];
					if (startIndex != 0xffffffff) 
					{        // cell is not empty
						// iterate over particles in this cell
						int endIndex = c_cellEnd[MCID];
						for(int x = startIndex; x < endIndex; x++)
						{
							
							Integer NeighborID = x;
							if(c_daOutputParticleID[NeighborID] < 0)
							{
								continue;
							}
							if(NeighborID != ID)
							{				
								Scalar Distance = CalcDistance(&TargetPosition, &c_daParticlePosition[NeighborID]);
								if(Distance < Re)
								{
									if(Distance < MinDist)
									{
										Distance = MinDist;
									}
									PCoefficient = 0.5 * (c_daParticlePressure[NeighborID] + TargetPressure) * GetWeight(Re, Distance) ;
									GradientFluid.x += PCoefficient * ( c_daParticlePosition[NeighborID].x - TargetPosition.x ) / (Distance * Distance);//FOR SSF * CONSTANT_PARAMETER.RefrencePressure);
									GradientFluid.y += PCoefficient * ( c_daParticlePosition[NeighborID].y - TargetPosition.y ) / (Distance * Distance);//FOR SSF * CONSTANT_PARAMETER.RefrencePressure);
									GradientFluid.z += PCoefficient * ( c_daParticlePosition[NeighborID].z - TargetPosition.z ) / (Distance * Distance);//FOR SSF * CONSTANT_PARAMETER.RefrencePressure);

									NeighborPressureSum += c_daParticlePressure[NeighborID];
									EnableNeighborID = NeighborID;
									++Cntr;							
								}								
							}
						}						
					}
				}
			}
		}
	}
	Scalar Coefficient1 = 2 * CONSTANT_PARAMETER.Dimension / (CONSTANT_PARAMETER.ConstantParticleDensity);
	Scalar Coefficient2 = SpecificDensity / (CONSTANT_PARAMETER.Dt * CONSTANT_PARAMETER.Dt);
	Scalar Coefficient3 = - CONSTANT_PARAMETER.Dt / SpecificDensity; 
	
	//Gradients for fluid
	GradientFluid.x *= Coefficient1 ;
	GradientFluid.y *= Coefficient1 ;
	GradientFluid.z *= Coefficient1 ;
	
	//Gradient for wall
	Scalar ParticlePressure = c_daParticlePressure[ID];
	Scalar WallWeight = 0.0;
	CDistance r_iw = c_daSTLDistance[ID];
	if(c_daSTLID[ID] >= 0 && c_daSTLID[ID] < c_TriangleNum)
	{
		//// Yamada et. al. wall model 
		if(r_iw.Magnitude > 0 && r_iw.Magnitude <= Re)
		{
			Scalar WallWeight = GetZValueFunctionGradient(r_iw.Magnitude,CONSTANT_PARAMETER.InitialParticleDistance);			
			Scalar WallCoefficient = 0.5 * Coefficient1 * WallWeight * TargetPressure / (r_iw.Magnitude * r_iw.Magnitude);
			GradientWall.x = WallCoefficient * (-r_iw.Direction.x);
			GradientWall.y = WallCoefficient * (-r_iw.Direction.y);
			GradientWall.z = WallCoefficient * (-r_iw.Direction.z);			
		}
		// penalty wall function harada et al. wall model 
		if(r_iw.Magnitude > 0 && r_iw.Magnitude <= CONSTANT_PARAMETER.InitialParticleDistance * CONSTANT_PARAMETER.OffSet)
		{
			Scalar dr = (CONSTANT_PARAMETER.InitialParticleDistance * CONSTANT_PARAMETER.OffSet - r_iw.Magnitude);
			Scalar myCoefficient = Coefficient2 *  dr / r_iw.Magnitude ;
			GradientWall.x += myCoefficient * (-r_iw.Direction.x);
			GradientWall.y += myCoefficient * (-r_iw.Direction.y);
			GradientWall.z += myCoefficient * (-r_iw.Direction.z);			
		}
	}
	
	
	//*********Negative pressure calculation*************
	if (CONSTANT_PARAMETER.bNegativePressure)
	{	if(NearestID != -1 && NearestDistance.Magnitude != 1e8)
		{
			if(c_daTrianglesParameters[c_daTriangles[NearestID].ModelIDNumber].RotationInRPM != 0 || c_daTrianglesParameters[c_daTriangles[NearestID].ModelIDNumber].NegativePressureCheck==true) 
			{
				if(NearestID != -1 && NearestDistance.Magnitude != 1e8)
				{
					if( NearestDistance.Magnitude <= c_daTrianglesParameters[c_daTriangles[NearestID].ModelIDNumber].InfluenceRegion * CONSTANT_PARAMETER.InitialParticleDistance)
					{
//						c_daOutputParticleType[ID] = TYPE_NEGATIVE_PRESSURE;  //for test only
						ImaginaryParticleWeight = GetZValueFunctionGradient (CONSTANT_PARAMETER.InitialParticleDistance, CONSTANT_PARAMETER.InitialParticleDistance);
						Scalar PCoefficient1 = 0.5 * (c_daTrianglesParameters[c_daTriangles[NearestID].ModelIDNumber].PressureGas + c_daParticlePressure[ID]) * ImaginaryParticleWeight / (CONSTANT_PARAMETER.InitialParticleDistance * CONSTANT_PARAMETER.InitialParticleDistance);
						GradientNegativePressure.x = PCoefficient1 * c_daTrianglesParameters[c_daTriangles[NearestID].ModelIDNumber].DirectionVector.x; 
						GradientNegativePressure.y = PCoefficient1 * c_daTrianglesParameters[c_daTriangles[NearestID].ModelIDNumber].DirectionVector.y;
						GradientNegativePressure.z = PCoefficient1 * c_daTrianglesParameters[c_daTriangles[NearestID].ModelIDNumber].DirectionVector.z; 


					}
				}
			}
		}
	}
	//Gradients for Negative Pressure
	GradientNegativePressure.x *= Coefficient1;
	GradientNegativePressure.y *= Coefficient1;
	GradientNegativePressure.z *= Coefficient1;
	//**********Negative Pressure Calculation ends********

	Scalar3 Gradient;
	//Gradient =  GradientFluid + GradientWall;
	//*******Gradient with Negative pressure
	Gradient =  GradientFluid + GradientWall + GradientNegativePressure;

	Scalar3 ImplicitVelocity;
	ImplicitVelocity.x = Coefficient3 * Gradient.x;
	ImplicitVelocity.y = Coefficient3 * Gradient.y;
	ImplicitVelocity.z = Coefficient3 * Gradient.z;	

	//ontest
	if(isfinite(ImplicitVelocity.x) && isfinite(ImplicitVelocity.y) && isfinite(ImplicitVelocity.z))
	{
		c_daOutputParticleVelocity[originalIndex].x = c_daParticleVelocity[ID].x + ImplicitVelocity.x ;
		c_daOutputParticleVelocity[originalIndex].y = c_daParticleVelocity[ID].y + ImplicitVelocity.y ;
		c_daOutputParticleVelocity[originalIndex].z = c_daParticleVelocity[ID].z + ImplicitVelocity.z ;

		CorrectVelocity(&CONSTANT_PARAMETER,c_daOutputParticleVelocity[originalIndex]);

		//modification to update position according to corrected velocity @Prashanta
		ImplicitVelocity.x = c_daOutputParticleVelocity[originalIndex].x - c_daParticleVelocity[ID].x;
		ImplicitVelocity.y = c_daOutputParticleVelocity[originalIndex].y - c_daParticleVelocity[ID].y;
		ImplicitVelocity.z = c_daOutputParticleVelocity[originalIndex].z - c_daParticleVelocity[ID].z;

		c_daOutputParticlePosition[originalIndex].x = TargetPosition.x + ImplicitVelocity.x * CONSTANT_PARAMETER.Dt;
		c_daOutputParticlePosition[originalIndex].y = TargetPosition.y + ImplicitVelocity.y * CONSTANT_PARAMETER.Dt;
		c_daOutputParticlePosition[originalIndex].z = TargetPosition.z + ImplicitVelocity.z * CONSTANT_PARAMETER.Dt;
	}
	//ontest end
}

__global__ static void CalcTurbulenceViscosity_Kernel(int ComputeParticleNum)
{
	Integer TID = CudaGetTargetID();
	Integer ID  =  TID;
	Integer OID =  TID;

	if(ID >= ComputeParticleNum)
	{
		return ;
	}
	OUTPUT_MACROINDEX();
	

	Scalar3 TargetPosition = c_daParticlePosition[ID];
	Scalar3 V0 = c_daParticleVelocity[ID];
	if(!IsInclude(&CONSTANT_BOUNDINGBOX,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
		c_daParticleTurbulaceViscosity[ID]=0.0;
		c_daOutputParticleID[originalIndex]		= -1;		
		return;
	}
	Scalar3 TensorVelocityXelement = make_Scalar3(0,0,0); // A,.,.
	Scalar3 TensorVelocityYelement = make_Scalar3(0,0,0);  // .,B,.
	Scalar3 TensorVelocityZelement = make_Scalar3(0,0,0); // .,.,C
	Scalar3 TensorNonDiagonal = make_Scalar3(0,0,0); //D,E,F components

	Scalar StrainRateTensor = 0.0;
	Scalar Coefficient = 0.0 ;
	Scalar SpecificDensity = CONSTANT_PARAMETER.Density;

	//For Vorticity
	Integer CIDX, CIDY, CIDZ;
	Integer CID = GetCellID(&CONSTANT_BOUNDINGBOX,&TargetPosition,CIDX, CIDY, CIDZ);
	if(CID < 0 || CID >= c_CellNum)
	{
		c_daParticlePressure[ID] = 0.0;
		return;
	}
	const Scalar Re = CONSTANT_PARAMETER.GradientInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	const Scalar MinDist = Re * 1e-3;
	
	Scalar Density = 0.0;
	Scalar3 VorticityVelocity = make_Scalar3(0,0,0); 
	Scalar3 VDistance = make_Scalar3(0,0,0); 
	Scalar3 Cross = make_Scalar3(0,0,0); 
	Scalar3 Vorticity = make_Scalar3(0,0,0);
	Integer Range = 1;
	Integer Neighbour = 0;
	for(Integer k = -Range; k <= Range; ++k)
	{
		for(Integer j = -Range; j <= Range; ++j)
		{
			for(Integer i = -Range; i <= Range; ++i)
			{
				Integer MCID = GetCellID(&CONSTANT_BOUNDINGBOX,CIDX +i, CIDY + j,CIDZ + k);
				if(MCID < 0 || MCID >= c_CellNum)
				{
					continue;
				}
				int startIndex = c_cellStart[MCID];
				if (startIndex != 0xffffffff)    // cell is not empty
				{       
					// iterate over particles in this cell
					int endIndex = c_cellEnd[MCID];
					for(int x = startIndex; x < endIndex; x++)
					{
						Integer NeighborID = x;
						if(c_daOutputParticleID[NeighborID] < 0)
						{
							continue;
						}
						if(NeighborID != ID)
						{				
							Scalar Distance = CalcDistance(&TargetPosition,&c_daParticlePosition[NeighborID]);
							if(Distance < Re)
							{
								if(Distance < MinDist)
								{
									Distance = MinDist;
								}
								Scalar Weight = GetWeight(Re, Distance);
								Density += Weight;
								Scalar3 V1 = c_daParticleVelocity[NeighborID];
								
								//VDistance.x = (c_daParticlePosition[NeighborID].x - TargetPosition.x);
								//VDistance.y = (c_daParticlePosition[NeighborID].y - TargetPosition.y);
								//VDistance.z = (c_daParticlePosition[NeighborID].z - TargetPosition.z);
								//
								//TensorVelocityXelement.x +=  ((V1.x - V0.x) * VDistance.x * Weight)/ (Distance * Distance); //A
								//TensorVelocityXelement.y +=  ((V1.x - V0.x) * VDistance.y * Weight)/ (Distance * Distance);
								//TensorVelocityXelement.z +=  ((V1.x - V0.x) * VDistance.z * Weight) / (Distance * Distance);
								//
								//TensorVelocityYelement.x +=  ((V1.y - V0.y) * VDistance.x * Weight) / (Distance * Distance);
								//TensorVelocityYelement.y +=  ((V1.y - V0.y) * VDistance.y * Weight) / (Distance * Distance); //B
								//TensorVelocityYelement.z +=  ((V1.y - V0.y) * VDistance.z * Weight) / (Distance * Distance);
								//
								//TensorVelocityZelement.x +=  ((V1.z - V0.z) * VDistance.x * Weight) / (Distance * Distance);
								//TensorVelocityZelement.y +=  ((V1.z - V0.z) * VDistance.y * Weight) / (Distance * Distance);
								//TensorVelocityZelement.z +=  ((V1.z - V0.z) * VDistance.z * Weight) / (Distance * Distance);//C			
								TensorVelocityXelement.x +=  ((V1.x - V0.x) * (c_daParticlePosition[NeighborID].x - TargetPosition.x) * Weight)/ (Distance * Distance); //A
								TensorVelocityXelement.y +=  ((V1.x - V0.x) * (c_daParticlePosition[NeighborID].y - TargetPosition.y) * Weight)/ (Distance * Distance);
								TensorVelocityXelement.z +=  ((V1.x - V0.x) * (c_daParticlePosition[NeighborID].z - TargetPosition.z)* Weight) / (Distance * Distance);
								
								TensorVelocityYelement.x +=  ((V1.y - V0.y) * (c_daParticlePosition[NeighborID].x - TargetPosition.x) * Weight) / (Distance * Distance);
								TensorVelocityYelement.y +=  ((V1.y - V0.y) * (c_daParticlePosition[NeighborID].y - TargetPosition.y) * Weight) / (Distance * Distance); //B
								TensorVelocityYelement.z +=  ((V1.y - V0.y) * (c_daParticlePosition[NeighborID].z - TargetPosition.z) * Weight) / (Distance * Distance);
								
								TensorVelocityZelement.x +=  ((V1.z - V0.z) * (c_daParticlePosition[NeighborID].x - TargetPosition.x) * Weight) / (Distance * Distance);
								TensorVelocityZelement.y +=  ((V1.z - V0.z) * (c_daParticlePosition[NeighborID].y -TargetPosition.y) * Weight) / (Distance * Distance);
								TensorVelocityZelement.z +=  ((V1.z - V0.z) * (c_daParticlePosition[NeighborID].z -TargetPosition.z) * Weight) / (Distance * Distance);//C														

								//Calculation of vorticity	starts
								VorticityVelocity.x = (V1.x - V0.x) ;
								VorticityVelocity.y = (V1.y - V0.y) ;
								VorticityVelocity.z = (V1.z - V0.z) ;

								Scalar VCoefficient = Weight/ (Distance * Distance);
								Cross = CrossProduct(VDistance,VorticityVelocity);
								Vorticity = Vorticity + Cross * VCoefficient;
								//calculation vorticity ends

								Neighbour++;
							}
						}
					}
				}
			}
		}
	}
	//Scalar GradientCoefficient = CONSTANT_PARAMETER.Dimension / CONSTANT_PARAMETER.ConstantParticleDensity;

//	Vorticity = Vorticity * GradientCoefficient; //calculation of vorticity vector
	Vorticity = Vorticity * (CONSTANT_PARAMETER.Dimension/CONSTANT_PARAMETER.ConstantParticleDensity);
	//c_daOutputParticleTemperature[originalIndex] = Magnitude(&Vorticity);
	
	/*TensorVelocityXelement.x = GradientCoefficient * TensorVelocityXelement.x;
	TensorVelocityXelement.y = GradientCoefficient * TensorVelocityXelement.y;
	TensorVelocityXelement.z = GradientCoefficient * TensorVelocityXelement.z;
	TensorVelocityYelement.x = GradientCoefficient * TensorVelocityYelement.x;
	TensorVelocityYelement.y = GradientCoefficient * TensorVelocityYelement.y;
	TensorVelocityYelement.z = GradientCoefficient * TensorVelocityYelement.z;
	TensorVelocityZelement.x = GradientCoefficient * TensorVelocityZelement.x;	
	TensorVelocityZelement.y = GradientCoefficient * TensorVelocityZelement.y;	
	TensorVelocityZelement.z = GradientCoefficient * TensorVelocityZelement.z;*/
	TensorVelocityXelement.x = (CONSTANT_PARAMETER.Dimension/CONSTANT_PARAMETER.ConstantParticleDensity) * TensorVelocityXelement.x;
	TensorVelocityXelement.y = (CONSTANT_PARAMETER.Dimension/CONSTANT_PARAMETER.ConstantParticleDensity) * TensorVelocityXelement.y;
	TensorVelocityXelement.z = (CONSTANT_PARAMETER.Dimension/CONSTANT_PARAMETER.ConstantParticleDensity) * TensorVelocityXelement.z;
	TensorVelocityYelement.x = (CONSTANT_PARAMETER.Dimension/CONSTANT_PARAMETER.ConstantParticleDensity) * TensorVelocityYelement.x;
	TensorVelocityYelement.y = (CONSTANT_PARAMETER.Dimension/CONSTANT_PARAMETER.ConstantParticleDensity) * TensorVelocityYelement.y;
	TensorVelocityYelement.z = (CONSTANT_PARAMETER.Dimension/CONSTANT_PARAMETER.ConstantParticleDensity) * TensorVelocityYelement.z;
	TensorVelocityZelement.x = (CONSTANT_PARAMETER.Dimension/CONSTANT_PARAMETER.ConstantParticleDensity) * TensorVelocityZelement.x;	
	TensorVelocityZelement.y = (CONSTANT_PARAMETER.Dimension/CONSTANT_PARAMETER.ConstantParticleDensity) * TensorVelocityZelement.y;	
	TensorVelocityZelement.z = (CONSTANT_PARAMETER.Dimension/CONSTANT_PARAMETER.ConstantParticleDensity) * TensorVelocityZelement.z;	

	
	TensorNonDiagonal.x = 0.5 * (TensorVelocityXelement.y + TensorVelocityYelement.x); //D
	TensorNonDiagonal.y = 0.5 * (TensorVelocityXelement.z + TensorVelocityZelement.x); //E
	TensorNonDiagonal.z = 0.5 * (TensorVelocityYelement.z + TensorVelocityZelement.y); //F

	// S = A^2+B^2+C^2+ 2D^2 + 2E^2 + 2F^2 ;
	StrainRateTensor = TensorVelocityXelement.x * TensorVelocityXelement.x + TensorVelocityYelement.y * TensorVelocityYelement.y + TensorVelocityZelement.z * TensorVelocityZelement.z 
		               + 2 * TensorNonDiagonal.x * TensorNonDiagonal.x + 2 * TensorNonDiagonal.y * TensorNonDiagonal.y + 2 * TensorNonDiagonal.z * TensorNonDiagonal.z;
	
	Scalar r_iw = c_daSTLDistance[ID].Magnitude;
	
	// Wall damping function used = 1-exp(-(r/re))
	Scalar WallDistanceRatio = - (r_iw/Re);
	Scalar dampingFunction = 1 - expf(WallDistanceRatio);
	
	Scalar SmagorinskyDamped = CONSTANT_PARAMETER.SmagorinskyConstant * dampingFunction; //modified Prashanta
	
	//Turbulent viscosity = (walldamping*smagorinskyconstant*filterwidth)*sqrt(2*S)

	Coefficient = (SmagorinskyDamped * CONSTANT_PARAMETER.FilterWidth) * (SmagorinskyDamped * CONSTANT_PARAMETER.FilterWidth); //modified Prashanta 
		
	Scalar TurbulanceViscosity = 0;
	if(StrainRateTensor > 0)
	{
		TurbulanceViscosity = Coefficient * sqrt(2 * StrainRateTensor);
	}

	//filter Turbulance viscosity @amit2013/07/07
	if(TurbulanceViscosity > CONSTANT_PARAMETER.TurbulanViscosityRatio)
	{
		TurbulanceViscosity = CONSTANT_PARAMETER.TurbulanViscosityRatio;
	}
	
	c_daParticleTurbulaceViscosity[ID] = TurbulanceViscosity;
	//c_daOutputParticleTemperature[originalIndex] = c_daParticleTurbulaceViscosity[ID]; //For debugging
	
	// calculate Turbulance Kinetic Energy @amit2013/07/07
	//float Cneu = 0.094; //0.05
	//float Cepsilon = 1.048; //1.0	
	Scalar TKCoeff = 2 * CONSTANT_PARAMETER.Cnue * CONSTANT_PARAMETER.InitialParticleDistance * CONSTANT_PARAMETER.InitialParticleDistance / CONSTANT_PARAMETER.Cepsilon ;
	Scalar TurbulanceKineticEnergy = StrainRateTensor * TKCoeff;

	if(TurbulanceKineticEnergy > CONSTANT_PARAMETER.TurbulanceIntensity)
	{
		TurbulanceKineticEnergy = CONSTANT_PARAMETER.TurbulanceIntensity;
	}

	c_daParticleStrainTensorProduct[ID] = TurbulanceKineticEnergy;
	//c_daOutputParticleTemperature[ID] = TurbulanceKineticEnergy; //for debugging	
}


