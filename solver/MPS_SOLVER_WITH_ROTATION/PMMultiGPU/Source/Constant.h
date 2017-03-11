#include <cuda_runtime.h>	
#include <cuda.h>

__constant__ CParameter CONSTANT_PARAMETER;		
__constant__ CGridBox CONSTANT_BOUNDINGBOX;
__constant__ ComputeType COMPUTE_TYPE;

__constant__ CParticle* c_daParticle;
__constant__ CParticle* c_daOutputParticle;
__constant__ Integer c_ParticleNum;

__constant__ CParticle* c_daBucket;
__constant__ Integer c_BucketNum;

__constant__ CTriangle * c_daTriangles;
__constant__ Integer c_TriangleNum;

__constant__ CTriangleParameters * c_daTrianglesParameters;

__constant__ CDragTriangle * c_daDragTriangles;
__constant__ Integer c_DragTriangleNum;
__constant__ Scalar3 * c_daDragAcc;
__constant__ Scalar* c_daDragTemperature;

__constant__ Integer c_MaxParticleNum;//Maximum allocated size

__constant__ CDistance *c_daSTLDistance;// Holds the distance between particle and triangle
__constant__ Integer* c_daSTLID;// Holds the ID of the neighbor triangle

__constant__ CCell* c_daCell;
__constant__ Integer c_CellNum;

__constant__ Integer* c_dParticleHash;// Stores Next ID's of the particle in the cell

__constant__ Scalar* c_dB;
__constant__ Scalar* c_dX;
__constant__ Scalar* c_dXo;
__constant__ Scalar* c_dEps;
__constant__ Scalar* c_dMInv;
__constant__ Scalar* c_dAii;

__constant__ Scalar3* c_daSubgridScaleStress;

//@Rajan STARTS for Individual Particle
__constant__	Integer*		c_daParticleID;
__constant__	Scalar3*		c_daParticlePosition;
__constant__	Scalar3*		c_daParticleVelocity;
__constant__	Scalar*			c_daParticlePressure;
__constant__	Scalar*			c_daParticleDensity;
__constant__	Scalar*			c_daParticleTemperature;
__constant__	Scalar*			c_daParticleKineticViscosity;
__constant__	Scalar*			c_daParticleSolidPhaseRate;
__constant__	ParticleType*	c_daParticleType;
//@Rajan ENDS for	 Individual Particle

//@Rajan STARTS for Individual Particle
__constant__	Integer*		c_daOutputParticleID;
__constant__	Scalar3*		c_daOutputParticlePosition;
__constant__	Scalar3*		c_daOutputParticleVelocity;
__constant__	Scalar*			c_daOutputParticlePressure;
__constant__	Scalar*			c_daOutputParticleDensity;
__constant__	Scalar*			c_daOutputParticleTemperature;
__constant__	Scalar*			c_daOutputParticleKineticViscosity;
__constant__	Scalar*			c_daOutputParticleSolidPhaseRate;
__constant__	ParticleType*	c_daOutputParticleType;
//@Rajan ENDS for individual particles

//Drag@Rajan
__constant__	DragParameter*	c_daSTLDragParameter;	//added by Arpan 20130128

