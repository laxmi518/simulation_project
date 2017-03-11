#pragma once


#define BUFFER_OFFSET(i) ((char*)NULL + (i))

bool IsGLError();
bool ISCUDAError();
bool ISCUDAError(cudaError_t Status);
//CCTStatusType CudaSafeCall(cudaError_t Status);
void Normalize(double& A, double& B, double& C);
void Swap(float& a, float&b);

__int64 SizeOfFile64( const char * szFileName ) ;
int SizeOfFile( const char * szFileName ) ;

static inline void CheckAndGrow(CBox & BoundingBox, Scalar3 & Position)
{
	BoundingBox.m_MaxBound.x = max(Position.x, BoundingBox.m_MaxBound.x);
	BoundingBox.m_MaxBound.y = max(Position.y, BoundingBox.m_MaxBound.y);
	BoundingBox.m_MaxBound.z = max(Position.z, BoundingBox.m_MaxBound.z);

	BoundingBox.m_MinBound.x = min(Position.x, BoundingBox.m_MinBound.x);
	BoundingBox.m_MinBound.y = min(Position.y, BoundingBox.m_MinBound.y);
	BoundingBox.m_MinBound.z = min(Position.z, BoundingBox.m_MinBound.z);
}

static inline void CheckAndGrow(CBox& Base, const CBox& input)
{
	Base.m_MaxBound.x = max(Base.m_MaxBound.x, input.m_MaxBound.x);
	Base.m_MaxBound.y = max(Base.m_MaxBound.y, input.m_MaxBound.y);
	Base.m_MaxBound.z = max(Base.m_MaxBound.z, input.m_MaxBound.z);

	Base.m_MinBound.x = min(Base.m_MinBound.x, input.m_MinBound.x);
	Base.m_MinBound.y = min(Base.m_MinBound.y, input.m_MinBound.y);
	Base.m_MinBound.z = min(Base.m_MinBound.z, input.m_MinBound.z);
}


static inline  void RenderQuad(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2)
{
	glPushMatrix();
	glLineWidth( 3.0 );
	glBegin(GL_LINE_STRIP);		
		glColor3f(1.0f, 0.0f, 0.0f);
		//bottom
		glVertex3f( x1, y1, z1);					
		glVertex3f(x2, y1, z1);															
		glVertex3f( x2, y2, z1);	
		glVertex3f(x1, y2, z1);
		glVertex3f( x1, y1, z1);			

		//left
		glVertex3f( x1, y1, z1);			
		glVertex3f(x1, y2, z1);
		glVertex3f(x1, y2, z2);															
		glVertex3f( x1, y1, z2);
		glVertex3f( x1, y1, z1);	

		//front
		glVertex3f( x1, y1, z1);			
		glVertex3f(x2, y1, z1);	
		glVertex3f(x2, y1, z2);															
		glVertex3f( x1, y1, z2);
		glVertex3f( x1, y1, z1);
	glEnd();
	glBegin(GL_LINE_STRIP);		
		
		//right
		glVertex3f(x2, y2, z2);	
		glVertex3f( x2, y2, z1);	
		glVertex3f( x2, y1, z1);			
		glVertex3f(x2, y1, z2);			
		glVertex3f(x2, y2, z2);

		//Back
		glVertex3f(x2, y2, z2);
		glVertex3f( x1, y2, z2);	
		glVertex3f( x1, y1, z2);			
		glVertex3f(x2, y1, z2);													
		glVertex3f(x2, y2, z2);

		//Top
		glVertex3f(x2, y2, z2);	
		glVertex3f( x2, y2, z1);
		glVertex3f( x1, y2, z1);			
		glVertex3f(x1, y2, z2);															
		glVertex3f(x2, y2, z2);
	glEnd();
	glPopMatrix();
}

//static inline DrawCone
//              (GLUquadric* quad
//              , GLdouble base
//              , GLdouble top
//              , GLdouble height
//              , GLint slices
//              , GLint stacks
//              )
static inline void DrawCone
              (GLUquadric* quad
              , GLdouble base
              , GLdouble height			 
              )
{
	gluCylinder(quad,base,0,height,50,10);
}
