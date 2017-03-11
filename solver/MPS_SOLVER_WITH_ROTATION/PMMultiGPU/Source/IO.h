#pragma once

namespace IO
{
	// Load the Input Particle File
	CCTStatusType Load(const std::string& FileName, Integer& ParticleNum, CParticle*& phParticle);
	
	// Load the Input Parameter File
	CCTStatusType Load(const std::string& FileName, CParameter*& pParameter);
	
	//// Load the STL File
	//CCTStatusType Load(const std::string& FileName, Integer& TriangleNum, CTriangle*& pTriangle);
	// Load the set of STL models
	CCTStatusType Load(const std::vector<std::string> &Models, std::vector<Integer>& TrianglesPerModel,Integer& TotalTriangles ,CTriangle*& pTriangle);

	// Save the particles
	/*CCTStatusType Save(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle, const std::vector<Scalar3>& Position, const std::vector<Scalar>& RotationAngle, Scalar CurrentTime, bool isAscii = true);
	CCTStatusType SaveAscii(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle, const std::vector<Scalar3>& Position, const std::vector<Scalar>& RotationAngle, Scalar CurrentTime);	
	CCTStatusType SaveBinary(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle, const std::vector<Scalar3>& Position, Scalar CurrentTime);
	*/
	// Save the particles
	CCTStatusType Save(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle, const std::vector<Scalar3>& Position, Scalar CurrentTime, bool isAscii,
					   const Integer* particleID,	const Scalar3* phParticlePosition, const Scalar3* particleVelocity, const Scalar* pressure, const Scalar* density, const Scalar* temperature,
					   const Scalar* kineticViscosity, const Scalar* solidPhaseRate, const ParticleType* type,const std::vector<Scalar>& RotationAngle);
	CCTStatusType SaveAscii(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle, const std::vector<Scalar3>& Position,const  Scalar CurrentTime,
							const Integer* particleID,	const Scalar3* phParticlePosition, const Scalar3* particleVelocity, const Scalar* pressure, const Scalar* density,const  Scalar* temperature,
					   const Scalar* kineticViscosity, const Scalar* solidPhaseRate, const ParticleType* type,const std::vector<Scalar>& RotationAngle);	
	//CCTStatusType SaveBinary(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle, const std::vector<Scalar3>& Position, Scalar CurrentTime);
	CCTStatusType SaveBinary(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle, const std::vector<Scalar3>& Position,const Scalar CurrentTime,
							const Integer* particleID,	const Scalar3* phParticlePosition, const Scalar3* particleVelocity, const Scalar* pressure, const Scalar* density, const Scalar* temperature,
							const Scalar* kineticViscosity,const  Scalar* solidPhaseRate, const ParticleType* type,const std::vector<Scalar>& RotationAngle);

	// GraphR
	CCTStatusType Dump(const std::string& FileName, const Integer ParticleNum, const CParticle* const phParticle);

	CCTStatusType Load(const std::string& FileName, CDragTriangle*& pDragTriangles, Integer& TotalTriangleNum, std::string& InnerPressureParameter);

	CCTStatusType SavePressures(const std::string& FileName, const Scalar P1, const Scalar P2 , const Scalar P3, const Scalar H_Zoint, const Scalar H_breather);

	CCTStatusType SaveHeight(const std::string& FileName, const Scalar Time, const Scalar Height);

	void saveParameterWithName(const std::string& FileName, const CParameter * pParameter,
							   Integer ModelNumber, const std::vector<Scalar3>& CenterOfRotation, const std::vector<Scalar3>& OtherPointOfRotation, const std::vector<Integer>& RotationRPM ,
							   const std::vector<Scalar>& WallFrictionCoefficient ,const std::vector<Scalar>& WallThermalConductivity ,const std::vector<Scalar>& WallThermalResistivity ,
							   const DragParameter * DragPrameter, Integer dragTriangleNum);

};
