#pragma once
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#define SQRT_MAGIC_F 0x5f3759df

//static   float  sqrt_stl(const float x)
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
static float sqrt_stl ( float x){
  float xhalf = 0.5f*x;
  int i = *(int*)&x;
  i = 0x5f3759df - (i>>1);
  x = *(float*)&i;
  x = x*(1.5f - xhalf*x*x);
  return 1/x;}

namespace Serializer
{
	template<typename T>
	struct Cood3D
	{
		T m_X;
		T m_Y;
		T m_Z;

		Cood3D()
			: m_X(0)
			, m_Y(0)
			, m_Z(0)
		{};

		Cood3D(T X, T Y, T Z)
			: m_X(X)
			, m_Y(Y)
			, m_Z(Z)
		{};
	};

	template<typename T>
	struct Vector3D
	{
		T m_X;
		T m_Y;
		T m_Z;

		Vector3D()
			: m_X(0)
			, m_Y(0)
			, m_Z(0)
		{};

		Vector3D(T X, T Y, T Z)
			: m_X(X)
			, m_Y(Y)
			, m_Z(Z)
		{};

		T SquaredMagnitude() const
		{
			return m_X * m_X + m_Y * m_Y + m_Z * m_Z;
		};

		T Magnitude() const
		{
			return sqrt_stl(SquaredMagnitude());
		};

		bool Normalize()
		{
			T Mag = Magnitude();
			T Epsilon = std::numeric_limits<T>::epsilon();
			if(Epsilon > Mag)
			{
				return false;
			}
			m_X /= Mag;
			m_Y /= Mag;
			m_Z /= Mag;
			return true;
		};
	};


	template<typename T>
	struct Triangle3D
	{
		Triangle3D()
		{};

		Cood3D<T> m_Coords[3];
		Vector3D<T> m_Normal;
	};

	template<typename T>
	class CTriangleContainer : public std::vector<Triangle3D<T>*>
	{
	public:
		CTriangleContainer()
		{};

		virtual ~CTriangleContainer()
		{
			Delete();
		};

		CTriangleContainer<T>& operator=(const CTriangleContainer<T>& obj)
		{
			Delete();
			for(size_t i = 0; i < obj.size(); ++i)
			{
				Triangle3D<T>* Triangle = new Triangle3D<T>();
				*Triangle = *obj[i];
				push_back(Triangle);
			}
			return (*this);
		};

		void Add(const CTriangleContainer<T>& obj)
		{
			for(size_t i = 0; i < obj.size(); ++i)
			{
				Triangle3D<T>* Triangle = new Triangle3D<T>();
				*Triangle = *obj[i];
				push_back(Triangle);
			}
		};

		void Delete()
		{
			for(size_t i = 0; i < (*this).size(); ++i)
			{
				delete (*this)[i];
			}
			clear();
		}
	};


	// reference : http://www.hiramine.com/programming/3dmodelfileformat/stlfileformat.html
	class CSTLSerializer
	{	
	public:
		static bool Load(const std::string& Filename, CTriangleContainer<float>& Triangles, bool IsAscii = false);
		static bool Write(const std::string& Filename, CTriangleContainer<float>& Triangles, bool IsAscii = false);
	};
}