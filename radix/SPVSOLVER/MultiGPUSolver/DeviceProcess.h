#pragma once 

//#define BLOCK_MAX_DIM 256	

#define PARAMETER "CONSTANT_PARAMETER"

#define PARAMETERCOEFFICIENT "CONSTANT_PARAMETER_COEFFICIENT"
#define GRIDPARAMS "CONSTANT_GRIDPARAMS"

#define BOUNDINGBOX "CONSTANT_BOUNDINGBOX"
	//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
	#define OUTPUTPARTICLESID "c_daOutputParticleID"
	#define OUTPUTPARTICLESPOSITION "c_daOutputParticlePosition"
	#define OUTPUTPARTICLESVELOCITY "c_daOutputParticleVelocity"
	#define OUTPUTPARTICLESPRESSURE "c_daOutputParticlePressure"
	#define OUTPUTPARTICLESDENSITY "c_daOutputParticleDensity"
	#define OUTPUTPARTICLESTEMPERATURE "c_daOutputParticleTemperature"
	#define OUTPUTPARTICLESKINETICVISCOSITY "c_daOutputParticleKineticViscosity"
	#define OUTPUTPARTICLESSOLIDPHASERATE "c_daOutputParticleSolidPhaseRate"
	#define OUTPUTPARTICLESTYPE "c_daOutputParticleType"
	//ENDS for individual particles----------------------E&T Nepal August 2011--------------------------
	//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
	#define PARTICLESID "c_daParticleID"
	#define PARTICLESPOSITION "c_daParticlePosition"
	#define PARTICLESVELOCITY "c_daParticleVelocity"
	#define PARTICLESPRESSURE "c_daParticlePressure"
	#define PARTICLESDENSITY "c_daParticleDensity"
	#define PARTICLESTEMPERATURE "c_daParticleTemperature"
	#define PARTICLESKINETICVISCOSITY "c_daParticleKineticViscosity"
	#define PARTICLESSOLIDPHASERATE "c_daParticleSolidPhaseRate"
	#define PARTICLESTYPE "c_daParticleType"
	//ENDS for individual particles----------------------E&T Nepal August 2011--------------------------
#define PARTICLENUM "c_ParticleNum"
#define BUCKET "c_daBucket"
	//STARTS for individual particles----------------------E&T Nepal August 2011--------------------------
	#define BUCKETID "m_daBucketID"
	#define BUCKETPOSITION "m_daBucketPosition"
	#define BUCKETVELOCITY "m_daBucketVelocity"
	#define BUCKETPRESSURE "m_daBucketPressure"
	#define BUCKETDENSITY "m_daBucketDensity"
	#define BUCKETTEMPERATURE "m_daBucketTemperature"
	#define BUCKETKINETICVISCOSITY "m_daBucketKineticViscosity"
	#define BUCKETSOLIDPHASERATE "m_daBucketSolidPhaseRate"
	#define BUCKETTYPE "m_daBucketType"
	//ENDS for individual particles----------------------E&T Nepal August 2011--------------------------
#define BUCKETNUM "c_BucketNum"
#define TRIANGLES "c_daTriangles"
#define TRIANGLENUM "c_TriangleNum"

//Set Triangle Parameter @Rajan 20121004
#define TRIANGLESPARAMETER "c_daTrianglesParameters"

#define DRAGTRIANGLES "c_daDragTriangles"
#define DRAGTRIANGLENUM "c_DragTriangleNum"
#define DRAGACC  "c_daDragAcc"
#define DRAGTEMPERATURE "c_daDragTemperature"

#define MAXPARTICLENUM "c_MaxParticleNum"
#define STLDISTANCE "c_daSTLDistance"
#define STLID "c_daSTLID"
#define CELLS "c_daCell"
#define CELLNUM "c_CellNum"

#define VECTB "c_dB"
#define VECTX "c_dX"
#define VECTS "c_dS"
#define VECTR "c_dR"
#define VECTD "c_dD"
#define VECTQ "c_dQ"

#define MINV "c_dMInv"
#define AII "c_dAii"

#define CELLSTART "c_cellStart"
#define CELLEND "c_cellEnd"
#define GRIDPARTICLEINDEX "c_gridParticleIndex"

extern "C"
{	
	CCTStatusType CalcTurbulenceViscosity(cudaStream_t &Stream,Integer ComputeParticleNum); //For Turbulance
	CCTStatusType CheckParticleOutsideComputeZone(cudaStream_t &Stream, Integer *ParticleNum);
	
	CCTStatusType CalcDragEffect(cudaStream_t &Stream,Integer ComputeParticleNum);
	CCTStatusType CalcExplicitly(cudaStream_t &Stream,Integer ParticleNum);
	CCTStatusType CalcExplicitPressure(cudaStream_t &Stream,Integer ComputeParticleNum);
	CCTStatusType CalcExplicitPressureGradient(cudaStream_t &Stream,Integer ComputeParticleNum);
	CCTStatusType CalcTemperatureFactor(cudaStream_t &Stream,Integer ComputeParticleNum);

	CCTStatusType RegisterTriangleTopology(cudaStream_t &Stream,CTriangle* daTriangle, Integer TriangleNum, CCell* daCell, Integer CellNum);
	CCTStatusType ResetTriangleTopology(cudaStream_t &Stream,Integer CellNum, CCell* aCell);
	CCTStatusType ResetWallPosition(cudaStream_t &Stream,const Integer TriangleNum,const Integer AnalysisStep, const CTriangle* daTriangles);
	CCTStatusType UpdateTrianglePosition(cudaStream_t &Stream,const Integer TriangleNum, CTriangle* daTriangles);
	CCTStatusType RotateTrianglePosition(cudaStream_t &Stream,const Integer TriangleNum, CTriangle* daTriangles, const Integer analysisStep);
	CCTStatusType CalcSTLDistance(cudaStream_t &Stream,Integer ComputeParticleNum);	
	
	CCTStatusType CaculateCellIDandInitializeHash(cudaStream_t &Stream,Integer ParticleNum,Integer CellNum,int* dGridParticleHash, int* dGridParticleIndex, Scalar3* particlePosition);  //August 2011
	CCTStatusType reorderDataAndFindCellStart(cudaStream_t &Stream,Integer numParticles,Integer numCells , int*  gridParticleHash, int*  gridParticleIndex, int*  cellStart, int*  cellEnd);
	
	///Implement Soring Alrithm
	CCTStatusType SortUsingThrust(Integer MaxParticleNum, Integer * daNumberHash, Integer* daNumberIndex);
	CCTStatusType StableSortUsingThrust(Integer MaxParticleNum, Integer * daNumberHash, Integer* daNumberIndex);

	//For Explicit Calcualtion
	CCTStatusType InitializeDeviceMemConst(CParameter Parameter,Integer ParticleNum,CTriangle * Triangles,Integer TriangleNum,CTriangleParameters * TriangleParameters,
		Integer MaxParticleNum, CDistance * STLDistance,Integer * StlID, CCell * Cell,Integer CellNum,
		Integer * CellStart,Integer * CellEnd,Integer * GridParticleIndex,CGridBox BoundingBox);

	CCTStatusType InitializeDeviceConstOutPutParticles(Integer * OutputParticleID, Scalar3 * OutputParticlePosition, Scalar3 * OutputParticleVelocity,
		Scalar * OutputParticlePressure, Scalar * OutputParticleDensity, Scalar * OutputParticleTemperature, Scalar * OutputParticleKineticViscosity,
		Scalar * OutputParticleSolidPhaseRate, ParticleType * OutputParticleType);

	CCTStatusType InitializeDeviceConstInputParticles(Integer * InputParticleID, Scalar3 * InputParticlePosition, Scalar3 * InputParticleVelocity,
		Scalar * InputParticlePressure, Scalar * InputParticleDensity, Scalar * InputParticleTemperature, Scalar * InputParticleKineticViscosity,
		Scalar * InputParticleSolidPhaseRate, ParticleType * InputParticleType, Scalar* ParticleTurbulaceViscosity, Scalar* ParticleStrainTensorProduct);
	
	CCTStatusType ParticleNumberToConst(Integer ParticleNum);

	CCTStatusType DragParametersToConst(DragParameter *InputDragParameter,Scalar3 *InputDragAcc, Scalar* InputDragTemperature , Integer DragTriangleNum, Integer * MagnifierCount,CDragTriangle * DragTriangles);

}
