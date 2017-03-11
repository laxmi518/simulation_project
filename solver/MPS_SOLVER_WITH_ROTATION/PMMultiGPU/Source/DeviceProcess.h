#pragma once 

#define BLOCK_MAX_DIM 256	
#define PARAMETER "CONSTANT_PARAMETER"
#define BOUNDINGBOX "CONSTANT_BOUNDINGBOX"
#define COMPUTETYPE "COMPUTE_TYPE"
#define OUTPUTPARTICLES "c_daOutputParticle"
#define PARTICLES "c_daParticle"
#define PARTICLENUM "c_ParticleNum"
#define BUCKET "c_daBucket"
#define BUCKETNUM "c_BucketNum"
#define TRIANGLES "c_daTriangles"
#define TRIANGLENUM "c_TriangleNum"

#define DRAGTRIANGLES "c_daDragTriangles"
#define DRAGTRIANGLENUM "c_DragTriangleNum"
#define DRAGACC  "c_daDragAcc"
#define DRAGTEMPERATURE "c_daDragTemperature"

#define MAXPARTICLENUM "c_MaxParticleNum"
#define STLDISTANCE "c_daSTLDistance"
#define STLID "c_daSTLID"
#define CELLS "c_daCell"
#define CELLNUM "c_CellNum"
#define PARTICLEHASH "c_dParticleHash"

#define MATB "c_dB"
#define MATX "c_dX"
#define MATX0 "c_dXo"
#define ERPS "c_dEps"
#define MINV "c_dMInv"
#define AII "c_dAii"

extern "C"
{
	CCTStatusType RegisterParticleTopology(cudaStream_t &Stream, Integer ParticleNum);
	CCTStatusType ResetParticleTopology(cudaStream_t &Stream,Integer CellNum, CCell* aCell, Integer* ParticleHash, Integer ParticleNum);
	CCTStatusType ResetCellParticleTopology(cudaStream_t &Stream,Integer CellNum, CCell* aCell);
	CCTStatusType CalcDistance(cudaStream_t &Stream,const Integer ParticleNum, CParticle* daParticle, Integer CellNum, 
									CCell* daCell,Integer* Hash ,Scalar* daParticleDistance, 
									Integer* daParticleID,Integer TriangleNum,
									CTriangle* daTriangles, CDistance* daSTLDistance,Integer* daTriangleID);
	CCTStatusType CalcParticleDistance(cudaStream_t &Stream,const Integer ParticleNum, CParticle* daParticle, Integer CellNum,
		CCell* daCell, Integer* Hash, Scalar* daParticleDistance, Integer* daParticleID);
	CCTStatusType CalcDensity(cudaStream_t &Stream,const Integer ParticleNum, CParticle* daParticle, 
									Scalar* daParticleDistance, Integer* daParticleID,Integer TriangleNum, 
									CTriangle* daTriangles, CDistance* daSTLDistance);
	CCTStatusType CalcExplicitly(cudaStream_t &Stream,Integer ParticleNum,Scalar InnerAcceleration);
	CCTStatusType CalcImplicitExplicitly(cudaStream_t &Stream,const Integer ParticleNum);
	CCTStatusType CalcPoisson(cudaStream_t &Stream,const Integer ParticleNum, CParticle* daParticle, 
									Scalar* daParticleDistance, Integer* daParticleID,Integer TriangleNum, 
									CTriangle* daTriangles, CDistance* daSTLDistance, 
									Scalar* dA, Scalar* dB, Scalar* dX, Scalar* dX0);
	CCTStatusType RegisterTriangleTopology(cudaStream_t &Stream,CTriangle* daTriangle, Integer TriangleNum, CCell* daCell, Integer CellNum);
	CCTStatusType RegisterDragTriangleTopology(cudaStream_t &Stream,CDragTriangle* daDragTriangle, Integer DragTriangleNum, CCell* daCell, Integer CellNum);
	CCTStatusType ResetTriangleTopology(cudaStream_t &Stream,Integer CellNum, CCell* aCell);
	CCTStatusType ResetWallPosition(cudaStream_t &Stream,const Integer TriangleNum,const Integer AnalysisStep,const CTriangle* daTriangles);
	CCTStatusType UpdateTrianglePosition(cudaStream_t &Stream,const Integer TriangleNum, CTriangle* daTriangles);

	CCTStatusType CalcTemperatureFactor(cudaStream_t &Stream,Integer ParticleNum);

	CCTStatusType CalcSTLDistance(cudaStream_t &Stream,Integer ParticleNum);	

	CCTStatusType RelocateParticleData(cudaStream_t &Stream,Integer ParticleNum);
	CCTStatusType RotateTrianglePosition(cudaStream_t &Stream,const Integer TriangleNum, const Integer analysisStep)	;
	CCTStatusType ResetTrianglePosition(cudaStream_t &Stream,const Integer TriangleNum, CTriangle* deviceTriangle, const CTriangle* OriginalTriangle);

	CCTStatusType CalcPressureExplicitly(cudaStream_t &Stream,Integer ParticleNum);

	CCTStatusType CalcPressureExplicit(cudaStream_t &Stream, const Integer ParticleNum);

	CCTStatusType InitializeDeviceConstOutPutParticles(Integer * OutputParticleID, Scalar3 * OutputParticlePosition, Scalar3 * OutputParticleVelocity,
		Scalar * OutputParticlePressure, Scalar * OutputParticleDensity, Scalar * OutputParticleTemperature, Scalar * OutputParticleKineticViscosity,
		Scalar * OutputParticleSolidPhaseRate, ParticleType * OutputParticleType);

	CCTStatusType InitializeDeviceConstInputParticles(Integer * InputParticleID, Scalar3 * InputParticlePosition, Scalar3 * InputParticleVelocity,
		Scalar * InputParticlePressure, Scalar * InputParticleDensity, Scalar * InputParticleTemperature, Scalar * InputParticleKineticViscosity,
		Scalar * InputParticleSolidPhaseRate, ParticleType * InputParticleType);

	CCTStatusType InitializeDeviceMemConst(CParameter Parameter,Integer ParticleNum,CTriangle * Triangles,Integer TriangleNum,CTriangleParameters * TriangleParameters,
		Integer MaxParticleNum, CDistance * STLDistance,Integer * StlID, CCell * Cell,Integer CellNum,
		Integer BucketNum,/*Integer * GridParticleIndex,	CGridParams GridParams,*/
		CGridBox BoundingBox, Integer * ParticleHash);
	CCTStatusType ParticleNumberToConst(Integer ParticleNum);

	//Drag@Rajan
	CCTStatusType DragParametersToConst(DragParameter *InputDragParameter,Scalar3 *InputDragAcc, Scalar* InputDragTemperature , Integer DragTriangleNum);
	CCTStatusType CalcDragEffect(cudaStream_t &Stream, Integer ComputeParticleNum);
}
