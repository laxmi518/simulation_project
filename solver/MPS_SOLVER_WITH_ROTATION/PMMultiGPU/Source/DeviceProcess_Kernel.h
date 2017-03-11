#include "DataType.h"
#include "Constant.h"
#include "Utility_Inline.h"

////In Use @Rajan
__global__ static void  RegisterParticleTopology_Kernel()
{
	const Integer ID = CudaGetTargetID();
	if( ID >= c_ParticleNum)
	{
		return;
	}
	if(c_daParticleID[ID] < 0)
	{
		c_daOutputParticleID[ID] = -1;
		return;
	}
	Scalar3 TargetPosition;
	TargetPosition = c_daParticlePosition[ID];	

	if(!IsInclude(&CONSTANT_BOUNDINGBOX.m_ComputeZone,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
		c_daOutputParticleID[ID]  = -1;
		/*
		c_daOutputParticlePosition[ID].x = CONSTANT_BOUNDINGBOX.m_ComputeZone.m_MaxBound.x + CONSTANT_PARAMETER.Tolerance;
		c_daOutputParticlePosition[ID].y = CONSTANT_BOUNDINGBOX.m_ComputeZone.m_MaxBound.y + CONSTANT_PARAMETER.Tolerance;
		c_daOutputParticlePosition[ID].z = CONSTANT_BOUNDINGBOX.m_ComputeZone.m_MaxBound.z + CONSTANT_PARAMETER.Tolerance;
		*/
	}
	else
	{
		memcpy(&c_daOutputParticleID[ID],				&c_daParticleID[ID],sizeof(Integer)); 
		memcpy(&c_daOutputParticlePosition[ID],			&c_daParticlePosition[ID],sizeof(Scalar3)); 
		memcpy(&c_daOutputParticleVelocity[ID],			&c_daParticleVelocity[ID],sizeof(Scalar3)); 
		memcpy(&c_daOutputParticlePressure[ID],			&c_daParticlePressure[ID],sizeof(Scalar)); 
		memcpy(&c_daOutputParticleDensity[ID],			&c_daParticleDensity[ID],sizeof(Scalar)); 
		memcpy(&c_daOutputParticleTemperature[ID],		&c_daParticleTemperature[ID],sizeof(Scalar)); 
		memcpy(&c_daOutputParticleKineticViscosity[ID],	&c_daParticleKineticViscosity[ID],sizeof(Scalar));
		memcpy(&c_daOutputParticleSolidPhaseRate[ID],	&c_daParticleSolidPhaseRate[ID],sizeof(Scalar));
		memcpy(&c_daOutputParticleType[ID],				&c_daParticleType[ID],sizeof(ParticleType)); 
	}	
	Integer CID = GetCellID(&CONSTANT_BOUNDINGBOX,&TargetPosition);
	if(CID < 0 || CID >= c_CellNum)
	{
		return ;
	}	
	// atomicExch prevents the race condition 
	// Add the ID int the HashID memory and get the old ID at the same time
	Integer HashID = atomicExch(&c_daCell[CID].m_HashID, ID);
	if(HashID < 0)
	{
		c_daCell[CID].m_ParticleID = ID;	
	}
	else if(HashID < c_ParticleNum)
	{
		c_dParticleHash[HashID] = ID;
	}	
}

////In Use @Rajan
__global__ static void  RelocateParticleData_Kernel()
{
	Integer ID = CudaGetTargetID();
	if(ID >= c_ParticleNum)
	{
		return ;
	}
	if(c_daOutputParticleID[ID]  < 0)
	{
		return;
	}

	Scalar3 TargetPosition = c_daOutputParticlePosition[ID] ;

	if(!IsInclude(&CONSTANT_BOUNDINGBOX.m_ComputeZone,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
		c_daParticleID[ID] = -1;
		c_daParticlePosition[ID].x = CONSTANT_BOUNDINGBOX.m_ComputeZone.m_MaxBound.x + CONSTANT_PARAMETER.Tolerance;
		c_daParticlePosition[ID].y = CONSTANT_BOUNDINGBOX.m_ComputeZone.m_MaxBound.y + CONSTANT_PARAMETER.Tolerance;
		c_daParticlePosition[ID].z = CONSTANT_BOUNDINGBOX.m_ComputeZone.m_MaxBound.z + CONSTANT_PARAMETER.Tolerance;
		return;
	}
	memcpy(&c_daParticleID[ID],					 &c_daOutputParticleID[ID],sizeof(Integer)); 
	memcpy(&c_daParticlePosition[ID],			 &c_daOutputParticlePosition[ID],sizeof(Scalar3)); 
	memcpy(&c_daParticleVelocity[ID],			 &c_daOutputParticleVelocity[ID],sizeof(Scalar3)); 
	memcpy(&c_daParticlePressure[ID],			 &c_daOutputParticlePressure[ID],sizeof(Scalar)); 
	memcpy(&c_daParticleDensity[ID],			 &c_daOutputParticleDensity[ID],sizeof(Scalar)); 
	memcpy(&c_daParticleTemperature[ID],		 &c_daOutputParticleTemperature[ID],sizeof(Scalar)); 
	memcpy(&c_daParticleKineticViscosity[ID],	 &c_daOutputParticleKineticViscosity[ID],sizeof(Scalar)); 
	memcpy(&c_daParticleSolidPhaseRate[ID],		 &c_daOutputParticleSolidPhaseRate[ID],sizeof(Scalar));
	memcpy(&c_daParticleType[ID],				 &c_daOutputParticleType[ID],sizeof(ParticleType)); 
}
////In Use @Rajan
__global__ static void  ResetParticleTopology_Kernel(Integer CellNum, CCell* aCell)
{
	Integer ID = CudaGetTargetID();
	if(ID >= CellNum)
	{
		return;
	}
	aCell[ID].m_HashID = -1;
	aCell[ID].m_ParticleID = -1;
}

//For Drag Model
__global__ static void  CalcDragEffect_Kernel()
{
	Integer ID = CudaGetTargetID();
	if(ID >= c_ParticleNum)
	{
		return ;
	}
	if(c_daOutputParticleID[ID]  < 0)
	{
		return;
	}

	CDistance NearestDistance;
	Integer NearestID = -1;
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

	c_daDragAcc[ID]	=	make_Scalar3(0.0,0.0,0.0);
	Scalar pressure	=	0.0; //For Velcoity
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
					/*if(mDragParam.InputVelocityToken1	==	'V')
					{
						c_daOutputParticleTemperature[ID] = pressure;
						c_daParticleTemperature[ID] = pressure;
					}*/
				}
				if(mDragParam.InputVelocityToken1	==	'T' ||	mDragParam.InputVelocityToken1	==	'B') //T for temerature and B for Both
				{
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
				else
				{
					return;
				}				
			}
			else if(mDragParam.DragToken == 'C') // C is used to Indicate its Constant Function
			{
				if(mDragParam.InputVelocityToken1	==	'V' ||	mDragParam.InputVelocityToken1	==	'v')
				{
					
					c_daOutputParticleVelocity[ID].x = mDragParam.ConstantVx;
					c_daOutputParticleVelocity[ID].y = mDragParam.ConstantVy;
					c_daOutputParticleVelocity[ID].z = mDragParam.ConstantVz;
					
					c_daParticleVelocity[ID].x = mDragParam.ConstantVx;
					c_daParticleVelocity[ID].y = mDragParam.ConstantVy;
					c_daParticleVelocity[ID].z = mDragParam.ConstantVz;
				}
			}
			else if(mDragParam.DragToken == 'F') // C is used to Indicate its Constant Function
			{
				if(mDragParam.InputVelocityToken1	==	'V' ||	mDragParam.InputVelocityToken1	==	'v')
				{
					/*
					//Calcuating Unit Vector
					Scalar3 UnitVector = make_Scalar3(0,0,0);
					Scalar Magnitue = Magnitude(c_daParticleVelocity[ID]);
					if( Magnitue > 0)
					{
						UnitVector.x = c_daParticleVelocity[ID].x / Magnitue;
						UnitVector.y = c_daParticleVelocity[ID].y / Magnitue;
						UnitVector.z = c_daParticleVelocity[ID].z / Magnitue;
					}
					
					c_daOutputParticleVelocity[ID].x = UnitVector.x * mDragParam.VelocityMagnitueFactor;
					c_daOutputParticleVelocity[ID].y = UnitVector.y * mDragParam.VelocityMagnitueFactor;
					c_daOutputParticleVelocity[ID].z = UnitVector.z * mDragParam.VelocityMagnitueFactor;
					*/
					c_daParticleVelocity[ID].x *=  mDragParam.VelocityMagnitueFactor;
					c_daParticleVelocity[ID].y *=  mDragParam.VelocityMagnitueFactor;
					c_daParticleVelocity[ID].z *=  mDragParam.VelocityMagnitueFactor;

					c_daOutputParticleVelocity[ID].x *=  mDragParam.VelocityMagnitueFactor;
					c_daOutputParticleVelocity[ID].y *=  mDragParam.VelocityMagnitueFactor;
					c_daOutputParticleVelocity[ID].z *=  mDragParam.VelocityMagnitueFactor;
				}
			}
		}		
	}	
}
////In Use @Rajan
__global__ static void CalcExplicitly_Kernel(Scalar InnerAcceleration)
{
	Integer ID = CudaGetTargetID();
	if(ID >= c_ParticleNum)
	{
		return ;
	}
	if(c_daOutputParticleID[ID]  < 0)
	{
		return;
	}
	Scalar SpecificDensity = CONSTANT_PARAMETER.Density;
	Scalar3 TargetPosition = c_daParticlePosition[ID] ;

	if(!IsInclude(&CONSTANT_BOUNDINGBOX.m_ComputeZone,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
		return;
	}

	//c_daParticle[ID].NeighborhoodNum = 0;
	Integer CIDX, CIDY, CIDZ;
	Integer CID = GetCellID(&CONSTANT_BOUNDINGBOX,&TargetPosition,CIDX, CIDY, CIDZ);
	//Integer Cntr = 0;
	Scalar3 Laplacian = make_Scalar3(0,0,0);
	Scalar3 ArtificialAcc = make_Scalar3(0,0,0);
	Scalar3 CollisionAcc = make_Scalar3(0,0,0);
	Scalar3 V0 = c_daParticleVelocity[ID] ;
	const Scalar Re = CONSTANT_PARAMETER.LaplacianInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	const Scalar Rc = CONSTANT_PARAMETER.r_coeff * CONSTANT_PARAMETER.InitialParticleDistance;

	Scalar Density = 0.0;
	if(CID >= 0 && CID < c_CellNum)
	{
		const Scalar MinDist = Re * 1e-3;
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
					Integer NeighborID = c_daCell[MCID].m_ParticleID;					
					while(NeighborID >= 0 && NeighborID < c_ParticleNum)
					{					
						if(NeighborID != ID)
						{				
							Scalar Distance = CalcDistance(&TargetPosition, & c_daParticlePosition[NeighborID] );
							if(Distance < Re)
							{
								if(Distance < MinDist)
								{
									Distance = MinDist;
								}
								Scalar Weight = GetWeight(Re, Distance);
								Density += Weight;
								Scalar3 V1 =  c_daParticleVelocity[NeighborID] ;
								Laplacian.x += (V1.x - V0.x) * Weight;
								Laplacian.y += (V1.y - V0.y) * Weight;
								Laplacian.z += (V1.z - V0.z) * Weight;

								if(CONSTANT_PARAMETER.bExplicitPressure)
								{
									Scalar Art = CONSTANT_PARAMETER.ArtificialPressure * Weight / (Distance * Distance);		
									ArtificialAcc = ArtificialAcc + (c_daParticlePosition[NeighborID]  - TargetPosition) * Art;
									if(Distance < Rc)
									{	
										Scalar Col = (Rc - Distance) / Distance;
										CollisionAcc = CollisionAcc + (c_daParticlePosition[NeighborID]  - TargetPosition) * Col;
									}
								}													
							}
						}
						NeighborID = c_dParticleHash[NeighborID];		
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
			if(2 == CONSTANT_PARAMETER.Dimension)
			{
				WallWeight = GetZValueFunctionLaplacian2D(Distance, CONSTANT_PARAMETER.InitialParticleDistance);
			}
			else
			{
				WallWeight = GetZValueFunctionLaplacian(Distance, CONSTANT_PARAMETER.InitialParticleDistance);
			}
			Density += WallWeight;
			//Individual Parameter
			//Scalar3 V1 = c_daTriangles[TriangleID].Velocity;
			Scalar3 V1 = c_daTrianglesParameters[c_daTriangles[TriangleID].ModelIDNumber].Velocity;

			Laplacian.x += CONSTANT_PARAMETER.WallFrictionCoeffieicnt * (V1.x - V0.x) * WallWeight;
			Laplacian.y += CONSTANT_PARAMETER.WallFrictionCoeffieicnt * (V1.y - V0.y) * WallWeight;
			Laplacian.z += CONSTANT_PARAMETER.WallFrictionCoeffieicnt * (V1.z - V0.z) * WallWeight;	
			if(CONSTANT_PARAMETER.bExplicitPressure)
			{
				Scalar Art = CONSTANT_PARAMETER.ArtificialPressure * WallWeight / (Distance * Distance);
				ArtificialAcc = ArtificialAcc + c_daSTLDistance[ID].Direction * Art;
				if(Distance < Rc)
				{
					Scalar Col = (Rc - Distance) / Distance;
					CollisionAcc = CollisionAcc +  c_daSTLDistance[ID].Direction * Col;
				}
			}
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
		DragVel = c_daDragAcc[ID];
	}
	
	Laplacian = Laplacian * c_daParticleKineticViscosity[ID] ;
	Acc = Acc + Laplacian;
	V0 = V0 + Acc * CONSTANT_PARAMETER.Dt - DragVel;   //Drag Velocity is subtracted from Main Velocity

	if(CONSTANT_PARAMETER.bExplicitPressure)
	{
		Scalar coeff = - CONSTANT_PARAMETER.Dimension / (CONSTANT_PARAMETER.ConstantParticleDensity * SpecificDensity);
		Acc = ArtificialAcc * coeff + CollisionAcc * 0.01;
		V0 = V0 + Acc * CONSTANT_PARAMETER.Dt;
	}
	if(isfinite(V0.x) && isfinite(V0.y) && isfinite(V0.z))
	{
		c_daOutputParticlePosition[ID].x = TargetPosition.x + CONSTANT_PARAMETER.Dt * V0.x;
		c_daOutputParticlePosition[ID].y = TargetPosition.y + CONSTANT_PARAMETER.Dt * V0.y;
		c_daOutputParticlePosition[ID].z = TargetPosition.z + CONSTANT_PARAMETER.Dt * V0.z;

		c_daOutputParticleVelocity[ID].x = V0.x;
		c_daOutputParticleVelocity[ID].y = V0.y;
		c_daOutputParticleVelocity[ID].z = V0.z;

		//c_daOutputParticleTemperature[ID] = WallWeight;
		//c_daOutputParticleSolidPhaseRate[ID] = c_daOutputParticleID[ID];
	}
}
////In Use @Rajan
__global__ static void  RegisterTriangleTopology_Kernel(const CTriangle* const daTriangle, Integer TriangleNum, CCell* daCell, Integer CellNum)
{
	Integer ID = CudaGetTargetID();
	if(ID >= TriangleNum)
	{
		return ;
	}
	RegisterTriangle(&CONSTANT_BOUNDINGBOX,&daTriangle[ID],ID,daCell, CellNum);
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

////In Use @Rajan
__global__ static void   UpdateTrianglePosition_Kernel(const Integer TriangleNum, CTriangle* daTriangles)
{
	Integer ID = CudaGetTargetID();
	if(ID >= TriangleNum)
	{
		return ;
	}
	//Scalar3 Vel =c_daTriangles[ID].Velocity;
	Scalar3 Vel = c_daTrianglesParameters[c_daTriangles[ID].ModelIDNumber].Velocity;

	Scalar dt = CONSTANT_PARAMETER.Dt;

	Vel.x = Vel.x * dt;
	Vel.y = Vel.y * dt;
	Vel.z = Vel.z * dt;

	c_daTriangles[ID].Vertex[0] = c_daTriangles[ID].Vertex[0] + Vel;
	c_daTriangles[ID].Vertex[1] = c_daTriangles[ID].Vertex[1] + Vel;
	c_daTriangles[ID].Vertex[2] = c_daTriangles[ID].Vertex[2] + Vel;
}

////In Use @Rajan
__global__ static void  CalcTemperatureFactor_Kernel()
{
	Integer ID = CudaGetTargetID();
	if(ID >= c_ParticleNum)
	{
		return ;
	}
	if(c_daOutputParticleID[ID]  < 0)
	{
		return;
	}
	Scalar Density = CONSTANT_PARAMETER.Density;

	const Scalar EnthalpyCoefficient = CONSTANT_PARAMETER.Dt * 2 * CONSTANT_PARAMETER.Dimension / 
		( CONSTANT_PARAMETER.ConstantParticleDensity * Density * CONSTANT_PARAMETER.FruidSpecificHeat);
	Scalar Hf = 0.0, Hw = 0.0;
	const Scalar ReTemperature = CONSTANT_PARAMETER.TemperatureInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	Scalar numerator = 0, denominator = 0;
	Scalar3 TargetPosition = c_daParticlePosition[ID] ;

	if(!IsInclude(&CONSTANT_BOUNDINGBOX.m_ComputeZone,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
		return;
	}

	Scalar	TargetTemperature = c_daParticleTemperature[ID] ;
	Scalar	ThermalResistance = 0;
	Scalar	ThermalConductivity = 0;	

	Integer CIDX, CIDY, CIDZ;
	Integer CID = GetCellID(&CONSTANT_BOUNDINGBOX,&TargetPosition,CIDX, CIDY, CIDZ);
	if(CID < 0 || CID >= c_CellNum)
	{
		return;
	}
	//const Scalar Re = CONSTANT_PARAMETER.LaplacianInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	const Scalar MinDist = ReTemperature * 1e-3;
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
				Integer NeighborID = c_daCell[MCID].m_ParticleID;					
				while(NeighborID >= 0 && NeighborID < c_ParticleNum)
				{					
					if(NeighborID != ID)
					{				
						Scalar Distance = CalcDistance(&TargetPosition, &c_daParticlePosition[NeighborID] );
						if(Distance < ReTemperature)
						{
							if(Distance < MinDist)
							{
								Distance = MinDist;
							}
							numerator = (c_daParticleTemperature[NeighborID]  - TargetTemperature) * GetWeight(ReTemperature, Distance);
							denominator = Distance * Distance * (1 / CONSTANT_PARAMETER.FruidThermalConductivity + CONSTANT_PARAMETER.FruidThermalResistance / Distance);
							Hf += numerator / denominator;
						}
					}
					NeighborID = c_dParticleHash[NeighborID];
				}
			}
		}
	}
	//ThermalResistance = CONSTANT_PARAMETER.WallThermalResistance;
	//ThermalConductivity = (2 * CONSTANT_PARAMETER.FruidThermalConductivity * CONSTANT_PARAMETER.WallThermalConductivity ) / ( CONSTANT_PARAMETER.FruidThermalConductivity + CONSTANT_PARAMETER.WallThermalConductivity);

	//c_daTrianglesParameters[c_daTriangles[ID].ModelIDNumber]

	Scalar Distance = c_daSTLDistance[ID].Magnitude;
	Integer TriangleID = c_daSTLID[ID];
	denominator = 0.0;
	numerator = 0.0;
	if(TriangleID >= 0 && TriangleID < c_TriangleNum)
	{
		if(Distance >0 && Distance < ReTemperature)
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
			//numerator = (c_daTriangles[TriangleID].Temperature - TargetTemperature) * WallWeight;
			numerator = (c_daTrianglesParameters[c_daTriangles[TriangleID].ModelIDNumber].Temperature - TargetTemperature) * WallWeight;

			denominator = Distance * Distance * (1 / ThermalConductivity + ThermalResistance / Distance);
			
			Hw = numerator / denominator;			
		}	
	}
// Temprature Effect from the Drag model STL
	Scalar HDrag = 0.0;
	if(c_DragTriangleNum > 0)
	{
		HDrag = c_daDragTemperature[ID];

		TargetTemperature += HDrag;
	}
// Ends Temprature Effect from the Drag model STL
	
	//TargetTemperature += (EnthalpyCoefficient * ( Hf + Hw + HDrag));

	Scalar HFHWComb = 0.0;
	if(abs(Hf) > SCALAR_EPSILON)
	{
		HFHWComb += Hf;
	}
	if(abs(Hw) > SCALAR_EPSILON)
	{
		HFHWComb += Hw;
	}

	TargetTemperature += (EnthalpyCoefficient * ( HFHWComb));

	//Limit the Temperature
	if(TargetTemperature > CONSTANT_PARAMETER.MaxWallTemperature)
	{
		TargetTemperature = CONSTANT_PARAMETER.MaxWallTemperature;
	}
	if(TargetTemperature < CONSTANT_PARAMETER.MinWallTemperature)
	{
		TargetTemperature = CONSTANT_PARAMETER.MinWallTemperature;
	}
	
	
	c_daOutputParticleKineticViscosity[ID]  = CONSTANT_PARAMETER.ViscosityUpdateParameterA * TargetTemperature + CONSTANT_PARAMETER.ViscosityUpdateParameterB;
	c_daOutputParticleTemperature[ID]  = TargetTemperature;	

	c_daOutputParticleSolidPhaseRate[ID] = (EnthalpyCoefficient * ( HFHWComb));
}

////In Use @Rajan
__global__ static void  CalcSTLDistance_Kernel()
{
	const Integer ID = CudaGetTargetID();
	if(ID >= c_ParticleNum)
	{
		return ;
	}
	CDistance NearestDistance;
	Integer NearestID = -1;
	NearestDistance.Magnitude = 1e8;
	NearestDistance.Direction = make_Scalar3(0,0,0);
	if(c_daOutputParticleID[ID]  < -1)
	{
		c_daSTLDistance[ID] = NearestDistance;
		c_daSTLID[ID] = NearestID;
		return;
	}

	Scalar3 TargetPosition = c_daParticlePosition[ID] ;	
	//c_daParticle[ID].NeighborhoodNum = 0;	
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

						if(c_daTrianglesParameters[c_daTriangles[TriangleID].ModelIDNumber].isDrag)
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
	c_daSTLDistance[ID] = NearestDistance;
	c_daSTLID[ID] = NearestID;
}

////In Use @Rajan
__global__ void CalcPressureExplicit_Kernel()
{
	const Integer ID = CudaGetTargetID();
	if(ID >= c_ParticleNum)
	{
		return ;
	}
	Scalar3 TargetPosition = c_daParticlePosition[ID];

	if(!IsInclude(&CONSTANT_BOUNDINGBOX.m_ComputeZone ,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
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
	/*Scalar CoefficientA = (-1) * SpecificDensity / (CONSTANT_PARAMETER.Dt * CONSTANT_PARAMETER.Dt);
	Scalar CoefficientB = 2 * CoefficientA / CONSTANT_PARAMETER.LambdaValue;*/
	Scalar Coeffieient = CONSTANT_PARAMETER.SonicVelocity * CONSTANT_PARAMETER.SonicVelocity * SpecificDensity / CONSTANT_PARAMETER.ConstantParticleDensity;
	Integer Cntr = 0;
	Scalar Density = 0.0;
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
				Integer NeighborID = c_daCell[MCID].m_ParticleID;					
				while(NeighborID >= 0 && NeighborID < c_ParticleNum)
				{					
					if(NeighborID != ID)
					{				
						Scalar Distance = CalcDistance(&TargetPosition, &c_daParticlePosition[NeighborID] );
						if(Distance < Re)
						{
							if(Distance < MinDist)
							{
								Distance = MinDist;
							}
							Scalar Weight1 = GetWeight(Re, Distance);
							Density = Density + Weight1;
							++Cntr;							
						}
					}
					NeighborID = c_dParticleHash[NeighborID];		
				}				
			}
		}
	}
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
	c_daOutputParticleDensity[ID]  = Density ;
	Scalar Pressure = 0.0;
	if(Density > CONSTANT_PARAMETER.ConstantParticleDensity * CONSTANT_PARAMETER.FreeSurfaceThreshold)
	{
		Pressure = Coeffieient * (Density - CONSTANT_PARAMETER.ConstantParticleDensity);
	}
	/*
	if(c_DragTriangleNum > 0)
	{
		//DragVel = c_daDragAcc[ID];
		Pressure -= c_daDragAcc[ID].x;
	}
	*/
	//check for finite Values
	if(isfinite(Pressure))
	{
		c_daOutputParticlePressure[ID] = Pressure;
		c_daParticlePressure[ID] = Pressure;	
	}
}

////In Use @Rajan
__global__ static void CalcImplicitExplicitly_Kernel()
{
	Integer ID = CudaGetTargetID();
	if(ID >= c_ParticleNum)
	{
		return ;
	}
	if(c_daOutputParticleID[ID] < 0)
	{
		return;
	}
	Scalar SpecificDensity = CONSTANT_PARAMETER.Density;

	Scalar3 TargetPosition = c_daParticlePosition[ID] ;

	if(!IsInclude(&CONSTANT_BOUNDINGBOX.m_ComputeZone,&TargetPosition,CONSTANT_PARAMETER.Tolerance))
	{
		return;
	}

	Integer CIDX, CIDY, CIDZ;
	Integer CID = GetCellID(&CONSTANT_BOUNDINGBOX,&TargetPosition,CIDX, CIDY, CIDZ);
	const Scalar Re = CONSTANT_PARAMETER.GradientInfluenceRadiusCoefficient * CONSTANT_PARAMETER.InitialParticleDistance;
	Scalar3 GradientFluid = make_Scalar3(0.0,0.0,0.0);
	Scalar3 GradientWall = make_Scalar3(0.0,0.0,0.0);
	Integer Cntr = 0;
	if(CID >= 0 && CID < c_CellNum)
	{
		const Scalar MinDist = Re * 1e-3;
		Scalar TargetPressure = c_daParticlePressure[ID] ;

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
					Integer NeighborID = c_daCell[MCID].m_ParticleID;					
					while(NeighborID >= 0 && NeighborID < c_ParticleNum)
					{					
						if(NeighborID != ID)
						{				
							Scalar Distance = CalcDistance(&TargetPosition, &c_daParticlePosition[NeighborID] );
							if(Distance < Re)
							{
								if(Distance < MinDist)
								{
									Distance = MinDist;
								}
								Scalar Weight = GetWeight(Re,Distance);
								Scalar PCoefficient = 0.5 * (c_daOutputParticlePressure[NeighborID] + TargetPressure) * Weight / (Distance * Distance);
								GradientFluid = GradientFluid + (c_daParticlePosition[NeighborID] - TargetPosition) * PCoefficient;
							}
						}
						NeighborID = c_dParticleHash[NeighborID];		
					}
				}
			}
		}
	}
	Scalar Coefficient1 = 2 * CONSTANT_PARAMETER.Dimension / (CONSTANT_PARAMETER.ConstantParticleDensity);
	Scalar Coefficient2 = SpecificDensity / (CONSTANT_PARAMETER.Dt * CONSTANT_PARAMETER.Dt);
	Scalar Coefficient3 = - CONSTANT_PARAMETER.Dt / SpecificDensity; 

	GradientFluid.x *= Coefficient1 ;
	GradientFluid.y *= Coefficient1 ;
	GradientFluid.z *= Coefficient1 ;

	CDistance r_iw = c_daSTLDistance[ID];
	if(c_daSTLID[ID] >= 0 && c_daSTLID[ID] < c_TriangleNum)
	{
		if(r_iw.Magnitude <= CONSTANT_PARAMETER.InitialParticleDistance && r_iw.Magnitude > 0)
		{
			Scalar dr = (CONSTANT_PARAMETER.InitialParticleDistance - r_iw.Magnitude);
			Scalar myCoefficient = Coefficient2 *  dr / r_iw.Magnitude ;
			GradientWall.x = myCoefficient * (-r_iw.Direction.x);
			GradientWall.y = myCoefficient * (-r_iw.Direction.y);
			GradientWall.z = myCoefficient * (-r_iw.Direction.z);
		}
	}
	Scalar3 Gradient;
	Gradient =  GradientFluid + GradientWall;
	Scalar3 ImplicitVelocity;

	ImplicitVelocity.x = Coefficient3 * Gradient.x;
	ImplicitVelocity.y = Coefficient3 * Gradient.y;
	ImplicitVelocity.z = Coefficient3 * Gradient.z;

	//Check for finite Values
	if(isfinite(ImplicitVelocity.x) && isfinite(ImplicitVelocity.y) && isfinite(ImplicitVelocity.z))
	{
		c_daOutputParticlePosition[ID].x = TargetPosition.x + ImplicitVelocity.x * CONSTANT_PARAMETER.Dt;
		c_daOutputParticlePosition[ID].y = TargetPosition.y + ImplicitVelocity.y * CONSTANT_PARAMETER.Dt;
		c_daOutputParticlePosition[ID].z = TargetPosition.z + ImplicitVelocity.z * CONSTANT_PARAMETER.Dt;


		c_daOutputParticleVelocity[ID].x = c_daParticleVelocity[ID].x + ImplicitVelocity.x ;
		c_daOutputParticleVelocity[ID].y = c_daParticleVelocity[ID].y + ImplicitVelocity.y ;
		c_daOutputParticleVelocity[ID].z = c_daParticleVelocity[ID].z + ImplicitVelocity.z ;

		CorrectVelocity(&CONSTANT_PARAMETER,c_daOutputParticleVelocity[ID]);
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

////In Use @Rajan
__global__ static void   RotateTrianglePosition_Kernel(const Integer TriangleNum, const Integer analysisStep)			
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
	{
		return;
	}	

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