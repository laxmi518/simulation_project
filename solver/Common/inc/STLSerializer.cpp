#include "stdafx.h"
#include "STLSerializer.h"
#include <boost/algorithm/string.hpp>
#include <algorithm>
namespace Serializer
{
	namespace
	{
		void SplitAndTrim(const std::string& Line, std::vector<std::string>& vSplitted)
		{
			using namespace boost::algorithm;
			std::string LineCopy(Line);
			trim(LineCopy);
			split(vSplitted, LineCopy, boost::is_space() /*is_any_of("\t ")*/, token_compress_on);
			for(size_t i = 0; i < vSplitted.size(); ++i)
			{
				trim(vSplitted[i]);
			}
		}

		bool ReadLine(std::istream& In, std::vector<std::string>& vLine)
		{
			if(!In)
			{
				return false;
			}

			while(!In.eof())
			{
				std::string Line;
				getline(In, Line);
				vLine.push_back(Line);
			}
			return true;
		}

		bool IsKeyward(const std::string& Str, const std::string& Keyword)
		{
			std::string Lower(Str);
			transform(Lower.begin(), Lower.end(), Lower.begin(), tolower);
			return 0 == Lower.compare(Keyword);
		}

		bool IsFacetKeyword(const std::string& Str)
		{
			return IsKeyward(Str, "facet");
		}

		bool IsOuterKeyword(const std::string& Str)
		{
			return IsKeyward(Str, "outer");
		}

		bool IsVertexKeyword(const std::string& Str)
		{
			return IsKeyward(Str, "vertex");
		}

		bool HasKeyword(const std::string& Src, const std::string& Keyword)
		{
			std::string Lower(Src);
			transform(Lower.begin(), Lower.end(), Lower.begin(), tolower);
			return std::string::npos != Lower.find(Keyword);
		}

		bool HasEndLoopKeyword(const std::string& Line)
		{
			return HasKeyword(Line, "endloop");
		}

		bool HasFacetKeyword(const std::string& Line)
		{
			return HasKeyword(Line, "facet");
		}

		bool HasEndFacetKeyword(const std::string& Line)
		{
			return HasKeyword(Line, "endfacet");
		}

		bool HasSolidKeyword(const std::string& Line)
		{
			return HasKeyword(Line, "solid");
		}

		bool HasEndSolidKeyword(const std::string& Line)
		{
			return HasKeyword(Line, "endsolid");
		}

		bool LoadFacet(std::vector<std::string>& vLine, std::vector<std::string>::iterator& it, CTriangleContainer<float>& Triangles)
		{
			using namespace boost::algorithm;
			if(it == vLine.end())
			{
				return true;
			}
			while(it != vLine.end())
			{
				if((*it).empty())
				{
					++it;
					continue;
				}
				// if(std::string::npos != it->find("endsolid"))
				if(HasEndSolidKeyword(*it))
				{
					return true;
				}

				std::vector<std::string> vSplitted;

				// load facet line
				SplitAndTrim(*it, vSplitted);
				// if(5 != vSplitted.size() || 0 != vSplitted.front().compare("facet"))
				if(5 != vSplitted.size() || !IsFacetKeyword(vSplitted.front()))
				{
					return false;
				}
				float NormalX = static_cast<float>(atof(vSplitted[2].c_str()));
				float NormalY = static_cast<float>(atof(vSplitted[3].c_str()));
				float NormalZ = static_cast<float>(atof(vSplitted[4].c_str()));
				Vector3D<float> Normal(NormalX, NormalY, NormalZ);
				Normal.Normalize();
				vSplitted.clear();
				++it;

				// load outer loop line
				SplitAndTrim(*it, vSplitted);
				// if(2 != vSplitted.size() || 0 != vSplitted.front().compare("outer"))
				if(2 != vSplitted.size() || !IsOuterKeyword(vSplitted.front()))
				{
					return false;
				}
				vSplitted.clear();
				++it;

				// load vertex line
				SplitAndTrim(*it, vSplitted);
				// if(4 != vSplitted.size() || 0 != vSplitted.front().compare("vertex"))
				if(4 != vSplitted.size() || !IsVertexKeyword(vSplitted.front()))
				{
					return false;
				}
				float Vertex0X = static_cast<float>(atof(vSplitted[1].c_str()));
				float Vertex0Y = static_cast<float>(atof(vSplitted[2].c_str()));
				float Vertex0Z = static_cast<float>(atof(vSplitted[3].c_str()));
				Cood3D<float> Coord0(Vertex0X, Vertex0Y, Vertex0Z);
				vSplitted.clear();
				++it;

				SplitAndTrim(*it, vSplitted);
				// if(4 != vSplitted.size() || 0 != vSplitted.front().compare("vertex"))
				if(4 != vSplitted.size() || !IsVertexKeyword(vSplitted.front()))
				{
					return false;
				}
				float Vertex1X = static_cast<float>(atof(vSplitted[1].c_str()));
				float Vertex1Y = static_cast<float>(atof(vSplitted[2].c_str()));
				float Vertex1Z = static_cast<float>(atof(vSplitted[3].c_str()));
				Cood3D<float> Coord1(Vertex1X, Vertex1Y, Vertex1Z);
				vSplitted.clear();
				++it;

				SplitAndTrim(*it, vSplitted);
				// if(4 != vSplitted.size() || 0 != vSplitted.front().compare("vertex"))
				if(4 != vSplitted.size() || !IsVertexKeyword(vSplitted.front()))
				{
					return false;
				}
				float Vertex2X = static_cast<float>(atof(vSplitted[1].c_str()));
				float Vertex2Y = static_cast<float>(atof(vSplitted[2].c_str()));
				float Vertex2Z = static_cast<float>(atof(vSplitted[3].c_str()));
				Cood3D<float> Coord2(Vertex2X, Vertex2Y, Vertex2Z);
				vSplitted.clear();
				++it;

				// load endloop line
				// if(std::string::npos == it->find("endloop"))
				if(!HasEndLoopKeyword(*it))
				{
					return false;
				}
				++it;

				// load endloop line
				// if(std::string::npos == it->find("endfacet"))
				if(!HasEndFacetKeyword(*it))
				{
					return false;
				}
				++it;

				Triangle3D<float>* Triangle = new Triangle3D<float>();
				Triangle->m_Normal = Normal;
				Triangle->m_Coords[0] = Coord0;
				Triangle->m_Coords[1] = Coord1;
				Triangle->m_Coords[2] = Coord2;
				Triangles.push_back(Triangle);
			}
			return true;
		}

		bool LoadSolid(std::vector<std::string>& vLine, std::vector<std::string>::iterator& it, CTriangleContainer<float>& Triangles)
		{
			if(it == vLine.end())
			{
				return true;
			}
			++it;
			while(it != vLine.end())
			{
				// if(std::string::npos != it->find("facet"))
				if(HasFacetKeyword(*it))
				{
					if(!LoadFacet(vLine, it, Triangles))
					{
						return false;
					}
				}
				++it;
			}
			return true;
		}

		bool LoadAscii(std::istream& In, CTriangleContainer<float>& Triangles)
		{
			if(!In)
			{
				return false;
			}

			std::vector<std::string> vLine;
			if(!ReadLine(In, vLine))
			{
				return false;
			}

			//load solid line
			std::vector<std::string>::iterator it = vLine.begin();
			bool HasSolid = false;
			while(it != vLine.end())
			{
				// if(std::string::npos != it->find("solid"))
				if(HasSolidKeyword(*it))
				{
					if(!LoadSolid(vLine, it, Triangles))
					{
						return false;
					}
					HasSolid = true;
				}
				if(it != vLine.end())
				{
					++it;
				}
			}
			return HasSolid ? true : false;
		}

		bool LoadBinary(std::istream& In, CTriangleContainer<float>& Triangles)
		{
			if(!In)
			{
				return false;
			}

			// load header
			char Header[80];
			In.read((char*)Header, sizeof(char) * 80);

			// load num triangles
			unsigned int NumTriangle = 0;
			In.read((char*)&NumTriangle, sizeof(unsigned int));

			// load triangles
			for(unsigned int i = 0; i < NumTriangle; ++i)
			{
				float NormalX = 0;
				In.read((char*)&NormalX, sizeof(float));
				float NormalY = 0;
				In.read((char*)&NormalY, sizeof(float));
				float NormalZ = 0;
				In.read((char*)&NormalZ, sizeof(float));
				Vector3D<float> Normal(NormalX, NormalY, NormalZ);

				float Coord0X = 0;
				In.read((char*)&Coord0X, sizeof(float));
				float Coord0Y = 0;
				In.read((char*)&Coord0Y, sizeof(float));
				float Coord0Z = 0;
				In.read((char*)&Coord0Z, sizeof(float));
				Cood3D<float> Coord0(Coord0X, Coord0Y, Coord0Z);

				float Coord1X = 0;
				In.read((char*)&Coord1X, sizeof(float));
				float Coord1Y = 0;
				In.read((char*)&Coord1Y, sizeof(float));
				float Coord1Z = 0;
				In.read((char*)&Coord1Z, sizeof(float));
				Cood3D<float> Coord1(Coord1X, Coord1Y, Coord1Z);

				float Coord2X = 0;
				In.read((char*)&Coord2X, sizeof(float));
				float Coord2Y = 0;
				In.read((char*)&Coord2Y, sizeof(float));
				float Coord2Z = 0;
				In.read((char*)&Coord2Z, sizeof(float));
				Cood3D<float> Coord2(Coord2X, Coord2Y, Coord2Z);

				char Buf[2];
				In.read((char*)Buf, sizeof(char)*2);

				Triangle3D<float>* Triangle = new Triangle3D<float>();
				Triangle->m_Normal = Normal;
				Triangle->m_Coords[0] = Coord0;
				Triangle->m_Coords[1] = Coord1;
				Triangle->m_Coords[2] = Coord2;
				Triangles.push_back(Triangle);
			}
			return true;
		}
	}


	bool
		CSTLSerializer::Load(const std::string& Filename, CTriangleContainer<float>& Triangles, bool IsAscii)
	{
		bool Status = true;
		if(IsAscii)
		{
			std::ifstream In(Filename.c_str());
			if(!In)
			{
				return false;
			}
			Status = LoadAscii(In, Triangles);
			In.close();
		}
		else
		{
			std::ifstream In(Filename.c_str(), std::ios::binary);
			if(!In)
			{
				return false;
			}
			Status = LoadBinary(In, Triangles);
			In.close();
		}
		return Status;
	};


	bool
		CSTLSerializer::Write(const std::string& Filename, CTriangleContainer<float>& Triangles, bool IsAscii)
	{
		return true;
	}
}