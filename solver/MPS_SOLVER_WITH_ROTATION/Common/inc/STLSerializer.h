#pragma once
#include <fstream>
#include <limits>
#include <string>
#include <vector>
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
			return sqrt(SquaredMagnitude());
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