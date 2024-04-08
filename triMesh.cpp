#include "triMesh.h" 
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <windows.h>


namespace TRIANGLE_MESH
{
	unsigned readNextData(char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize)
	{
		unsigned nIndx = 0;

		while ((pszBuf[0] == ' ') ||
			(pszBuf[0] == '\n') ||
			(pszBuf[0] == '\t') ||
			(pszBuf[0] == '\r'))
		{
			pszBuf++;
			nCount++;
		}

		while ((pszBuf[0] != ' ') &&
			(pszBuf[0] != '\n') &&
			(pszBuf[0] != '\t') &&
			(pszBuf[0] != '\r') &&
			(pszBuf[0] != '\null') &&
			(pszBuf[0] != 0) &&
			(nIndx < nMaxSize))
		{
			validData[nIndx++] = pszBuf[0];
			pszBuf++;
			nCount++;
		}

		validData[nIndx] = 0;
		return nIndx;
	};


	// 自定义顶点哈希函数：
	template <typename TV>
	class vertHash
	{
	public:
		std::size_t operator()(const triplet<TV>& v) const
		{
			return std::hash<decltype(v.x)>()(v.x) + \
				std::hash<decltype(v.y)>()(v.y) + \
				std::hash<decltype(v.z)>()(v.z);
		}
	};


	// 自定义顶点等价比较器：
	template <typename TV>
	class vertComp
	{
	public:
		bool operator()(const triplet<TV>& v0, const triplet<TV>& v1) const
		{
			double dbThreshold = 1.e-14;
			if (std::fabs(v0.x - v1.x) > dbThreshold)
				return false;
			if (std::fabs(v0.y - v1.y) > dbThreshold)
				return false;
			if (std::fabs(v0.z - v1.z) > dbThreshold)
				return false;
			return true;
		}
	};


	template <typename TV, typename TI>
	void GetVertsAndSurfs(std::vector<triplet<TV>>& vVerts, std::vector<triplet<TI>>& vSurfs, \
		const std::vector<triplet<TV>>& versIn)
	{
		// 修改说明：原先去除重复点的方法时间复杂度过高，改用hashmap
		unsigned nOldVertCnt = versIn.size();
		std::vector<unsigned> tmpTri(3, 0);
		vSurfs.resize(nOldVertCnt / 3);
		vVerts.reserve(nOldVertCnt);
		std::unordered_map<triplet<TV>, unsigned, vertHash<TV>, vertComp<TV>> mapVerts;

		for (unsigned i = 0; i < nOldVertCnt / 3; i++)
		{
			for (unsigned k = 0; k < 3; k++)
			{
				unsigned nOldIdx = i * 3 + k;
				const triplet<TV>& v = versIn[nOldIdx];
				if (0 == mapVerts.count(v))
				{
					mapVerts.insert(std::make_pair(v, vVerts.size()));
					vVerts.push_back(v);
				}
				auto iter = mapVerts.find(v);
				tmpTri[k] = iter->second;
			}

			vSurfs[i].x = static_cast<TI>(tmpTri[0]);
			vSurfs[i].y = static_cast<TI>(tmpTri[1]);
			vSurfs[i].z = static_cast<TI>(tmpTri[2]);
		}
	}


	// 三角片数据转换为有向边数据（无重复边）
	template <typename TI>
	void tris2edges(std::vector<doublet<TI>>& edges, const std::vector<triplet<TI>>& tris)
	{
		auto encodeEdge = [](const TI vaIdx, const TI vbIdx)->std::int64_t
		{
			std::int64_t a = static_cast<std::int64_t>(vaIdx);
			std::int64_t b = static_cast<std::int64_t>(vbIdx);
			std::int64_t code = 0;
			code |= (a << 32);
			code |= b;
			return code;
		};

		auto decodeEdge = [](const std::int64_t code) ->doublet<TI>
		{
			TI a = static_cast<TI>(code >> 32);
			TI b = static_cast<TI>(code - (static_cast<std::int64_t>(a) << 32));
			return doublet<TI>(a, b);
		};

		if (tris.empty())
			return;

		std::unordered_set<std::int64_t> set;
		for (const auto& t : tris)
		{
			set.insert(encodeEdge(t.x, t.y));
			set.insert(encodeEdge(t.y, t.z));
			set.insert(encodeEdge(t.z, t.x));
		}

		edges.clear();
		edges.reserve(set.size());
		for (const auto& code : set)
			edges.push_back(decodeEdge(code));
	}


	template	<typename TV>
	bool objReadVertices(std::vector<TRIANGLE_MESH::triplet<TV>>& vers, const char* fileName)
	{
		char* pTmp = NULL;
		std::ifstream ifs(fileName);			// cube bunny Eight
		if (false == ifs.is_open())
			return false;

		std::streampos   pos = ifs.tellg();			//  save   current   position   
		ifs.seekg(0, std::ios::end);
		unsigned fileLen = (unsigned)ifs.tellg();
		if (0 == fileLen)
			return false;

		ifs.seekg(pos);				  //   restore   saved   position   
		char* pFileBuf = new char[fileLen + 1];
		std::memset(pFileBuf, 0, fileLen + 1);
		ifs.read(pFileBuf, fileLen);
		char tmpBuffer[1024];
		unsigned nMaxSize = 1024;
		pTmp = pFileBuf;
		unsigned nReadLen = 0;
		unsigned nRet = 0;

		vers.clear();
		TRIANGLE_MESH::triplet<TV> ver;
		while (nReadLen < fileLen)
		{
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;

			// 顶点信息		
			if (std::strcmp(tmpBuffer, "v") == 0)
			{
				nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
				if (0 == nRet)
					break;

				ver.x = static_cast<TV>(atof(tmpBuffer));
				nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
				if (0 == nRet)
					break;
				ver.y = static_cast<TV>(atof(tmpBuffer));
				nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
				if (0 == nRet)
					break;
				ver.z = static_cast<TV>(atof(tmpBuffer));

				vers.push_back(ver);
			}
			else
				break;
		}
		delete[] pFileBuf;

		return true;
	};


	template	<typename TV>
	bool objWriteVertices(const char* fileName, const std::vector<TRIANGLE_MESH::triplet<TV>>& vers)
	{
		if (vers.empty())
			return false;

		std::ofstream dstFile(fileName);
		for (const auto& ver : vers)
		{
			char szBuf[1024] = { 0 };
			double x0 = static_cast<double>(ver.x);
			double y0 = static_cast<double>(ver.y);
			double z0 = static_cast<double>(ver.z);
			snprintf(szBuf, 1024, "v %.17g %.17g %.17g", x0, y0, z0);
			dstFile << szBuf << "\n";
		}
		dstFile.close();

		return true;
	}


	template <typename TV, typename TI>
	bool objWriteEdges(const char* fileName, const std::vector<TRIANGLE_MESH::triplet<TV>>& vers, \
		const std::vector<TRIANGLE_MESH::doublet<TI>>& edges)
	{
		if (edges.empty() || vers.empty())
			return false;

		std::ofstream dstFile(fileName);
		for (size_t i = 0; i < vers.size(); ++i)
			dstFile << "v " << vers[i].x << " " << vers[i].y << " " << vers[i].z << std::endl;

		for (size_t i = 0; i < edges.size(); ++i)
			dstFile << "l " << static_cast<int>(edges[i].x) + 1 << " "\
			<< static_cast<int>(edges[i].y) + 1 << std::endl;

		dstFile.close();

		return true;
	}


	// OBJ文件读取网格
	template	<typename TV, typename TI>
	bool objReadMesh(TRIANGLE_MESH::triMesh<TV, TI>& mesh, const char* fileName)
	{
		char* pTmp = NULL;
		std::ifstream ifs(fileName);		//cube bunny Eight
		if (false == ifs.is_open())
			return false;

		std::streampos   pos = ifs.tellg();     //   save   current   position   
		ifs.seekg(0, std::ios::end);
		unsigned fileLen = (unsigned)ifs.tellg();
		if (0 == fileLen)
			return false;

		ifs.seekg(pos);     //   restore   saved   position   
		char* pFileBuf = new char[fileLen + 1];
		std::memset(pFileBuf, 0, fileLen + 1);
		ifs.read(pFileBuf, fileLen);
		char tmpBuffer[1024];
		unsigned nMaxSize = 1024;
		pTmp = pFileBuf;
		unsigned nReadLen = 0;
		unsigned nRet = 0;

		mesh.vertices.clear();
		mesh.triangles.clear();
		TRIANGLE_MESH::triplet<TV> ver;
		TRIANGLE_MESH::triplet<TI> tri;
		while (nReadLen < fileLen)
		{
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;

			if (std::strcmp(tmpBuffer, "v") == 0)
			{
				nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
				if (0 == nRet)
					break;

				ver.x = static_cast<TV>(atof(tmpBuffer));
				nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
				if (0 == nRet)
					break;
				ver.y = static_cast<TV>(atof(tmpBuffer));
				nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
				if (0 == nRet)
					break;
				ver.z = static_cast<TV>(atof(tmpBuffer));
				mesh.vertices.push_back(ver);
			}
			else if (std::strcmp(tmpBuffer, "f") == 0)
			{
				nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
				if (0 == nRet)
					break;

				tri.x = static_cast<TI>(atoi(tmpBuffer) - 1);
				nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
				if (0 == nRet)
					break;
				tri.y = static_cast<TI>(atoi(tmpBuffer) - 1);

				nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
				if (0 == nRet)
					break;
				tri.z = static_cast<TI>(atoi(tmpBuffer) - 1);
				mesh.triangles.push_back(tri);
			}
		}
		delete[] pFileBuf;

		return true;
	}


	template	<typename TV, typename TI>
	bool objWriteMesh(const char* fileName, const TRIANGLE_MESH::triMesh<TV, TI>& mesh)
	{
		if (mesh.vertices.empty())
			return false;
		std::ofstream dstFile(fileName);
		if (false == dstFile.is_open())
			return false;

		for (const auto& ver : mesh.vertices)
		{
			char szBuf[1024] = { 0 };
			snprintf(szBuf, 1024, "v %.17g %.17g %.17g", ver.x, ver.y, ver.z);
			dstFile << szBuf << "\n";
		}

		for (const auto& tri : mesh.triangles)
		{
			char szBuf[256] = { 0 };
			snprintf(szBuf, 256, "f %d %d %d", tri.x + 1, tri.y + 1, tri.z + 1);
			dstFile << szBuf << "\n";
		}

		return true;
	}


	template	<typename TV, typename TI>
	bool stlReadMesh(triMesh<TV, TI>& mesh, const char* fileName, const bool blIsAscii = false)
	{
		using verV = triplet<TV>;
		std::vector<verV> tmpVers;

		if (blIsAscii)
		{
			std::ifstream fin(fileName, std::ios::in);

			fin.seekg(0, std::ios::end);	//seek to the end
			unsigned fileLen = (unsigned)fin.tellg();
			if (0 == fileLen)					// file is empty 
				return false;
			fin.seekg(0, std::ios::beg);	//seek to the beg

			char* pFileBuf = new char[fileLen + 1];
			std::memset(pFileBuf, 0, fileLen + 1);
			fin.read(pFileBuf, fileLen);

			char* pTemp = pFileBuf;
			char tempBuffer[1024];
			unsigned nMaxSize = 1024;
			unsigned nReadLen = 0;
			unsigned nRet = 0;

			while (nReadLen < fileLen)
			{
				nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
				if (0 == nRet)
					break;
				if (std::strcmp(tempBuffer, "vertex") == 0)    //顶点信息
				{
					verV vert;
					nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
					if (0 == nRet)
						break;
					vert.x = static_cast<TV>(atof(tempBuffer));
					nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
					if (0 == nRet)
						break;
					vert.y = static_cast<TV>(atof(tempBuffer));
					nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
					if (0 == nRet)
						break;
					vert.z = static_cast<TV>(atof(tempBuffer));
					tmpVers.push_back(vert);
				}
			}
			delete(pFileBuf);

			GetVertsAndSurfs(mesh.vertices, mesh.triangles, tmpVers);

			return true;
		}
		else
		{
			std::ifstream fin(fileName, std::ios::in | std::ios::binary);

			fin.seekg(0, std::ios::end);   //seek to the end
			unsigned fileLen = (unsigned)fin.tellg();
			if (0 == fileLen)					// file is empty 
				return false;

			fin.seekg(0, std::ios::beg);
			unsigned len = fin.tellg();
			char* buffer = new char[fileLen + 1];
			std::memset(buffer, 0, fileLen + 1);
			fin.read(buffer, fileLen);

			unsigned offset = 80;			// 跳过最开始的文件头（存贮文件名，80个字节）；
			unsigned trisCount = *(std::uint32_t*)(buffer + offset);		// 三角片数量；
			offset += 4;							// 二进制stl文件中，坐标都是REAL32、索引都是UINT32, 都是4字节；

			//从二进制文件读取顶点信息
			verV pt{ 0, 0, 0 };
			tmpVers.resize(trisCount * 3);

			for (unsigned i = 0; i < trisCount; i++)
			{
				offset += 4 * 3;					//normal
				for (unsigned k = 0; k < 3; k++)
				{
					pt.x = static_cast<TV>(*(float*)(buffer + offset));
					offset += 4;
					pt.y = static_cast<TV>(*(float*)(buffer + offset));
					offset += 4;
					pt.z = static_cast<TV>(*(float*)(buffer + offset));
					offset += 4;
					tmpVers[3 * i + k] = pt;
				}
				offset += 2;
			}
			delete(buffer);

			GetVertsAndSurfs(mesh.vertices, mesh.triangles, tmpVers);

			return true;
		}

		return true;
	}


	// 网格写入到二进制STL文件中：
	template <typename TV, typename TI>
	bool stlWriteMesh(const char* fileName, \
		const triMesh<TV, TI>& mesh, const bool blCalcNorms = false)
	{
		// lambda——计算网格所有三角片的法向：
		auto getTriNorms = [&mesh](std::vector<verF>& triNorms)->bool
		{
			triNorms.clear();
			triNorms.resize(mesh.triangles.size());
			for (int i = 0; i < mesh.triangles.size(); i++)
			{
				verF v1 = mesh.vertices[mesh.triangles[i].y] - mesh.vertices[mesh.triangles[i].x];
				verF v2 = mesh.vertices[mesh.triangles[i].z] - mesh.vertices[mesh.triangles[i].x];
				triNorms[i] = v1.cross(v2);
				float r = triNorms[i].length();
				if (r == 0)
				{
					std::cout << "Error!!! degenerate triangle detected, unable to calculate its normal dir." << std::endl;
					return false;
				}
				else
					triNorms[i] /= r;
			}
			return true;
		};

		std::ofstream fout(fileName, std::ios::binary);
		char title[80] = { 0 };
		std::int32_t trisCount = static_cast<std::int32_t>(mesh.triangles.size());
		std::vector<verF> triNorms, versF;
		versF.resize(mesh.vertices.size());
		for (int i = 0; i < mesh.vertices.size(); ++i)
			versF[i] = mesh.vertices[i].cast<float>();
		if (blCalcNorms)
			getTriNorms(triNorms);
		else
			triNorms.resize(trisCount, verF{ 0, 0, 0 });

		// 写数据到二进制stl文件：
		fout.write(title, 80);
		fout.write((char*)&trisCount, sizeof(std::int32_t));
		for (int i = 0; i < trisCount; i++)
		{
			TI vaIdx = mesh.triangles[i].x;
			TI vbIdx = mesh.triangles[i].y;
			TI vcIdx = mesh.triangles[i].z;
			fout.write((char*)&(triNorms[i].x), 4);
			fout.write((char*)&(triNorms[i].y), 4);
			fout.write((char*)&(triNorms[i].z), 4);
			fout.write((char*)&(versF[vaIdx].x), 4);
			fout.write((char*)&(versF[vaIdx].y), 4);
			fout.write((char*)&(versF[vaIdx].z), 4);
			fout.write((char*)&(versF[vbIdx].x), 4);
			fout.write((char*)&(versF[vbIdx].y), 4);
			fout.write((char*)&(versF[vbIdx].z), 4);
			fout.write((char*)&(versF[vcIdx].x), 4);
			fout.write((char*)&(versF[vcIdx].y), 4);
			fout.write((char*)&(versF[vcIdx].z), 4);

			char triAttr[2] = { 0 };
			fout.write(triAttr, 2);
		}

		return true;
	}
}
using namespace TRIANGLE_MESH;


////////////////////////////////////////////////////////////////////////////////////////////// IO methods:
bool readOBJ(std::vector<verF>& vers, const char* fileName)
{
	return TRIANGLE_MESH::objReadVertices(vers, fileName);
}


bool readOBJ(std::vector<verD>& vers, const char* fileName)
{
	return TRIANGLE_MESH::objReadVertices(vers, fileName);
}


bool readOBJ(triMeshF& mesh, const char* fileName)
{
	return TRIANGLE_MESH::objReadMesh(mesh, fileName);
}


bool readOBJ(triMeshD& mesh, const char* fileName)
{
	return TRIANGLE_MESH::objReadMesh(mesh, fileName);
}


bool writeOBJ(const char* fileName, const std::vector<verF>& vers)
{
	return TRIANGLE_MESH::objWriteVertices(fileName, vers);
}


bool writeOBJ(const char* fileName, const std::vector<verD>& vers)
{
	return TRIANGLE_MESH::objWriteVertices(fileName, vers);
}


bool writeOBJ(const char* fileName, const triMeshF& mesh)
{
	return TRIANGLE_MESH::objWriteMesh(fileName, mesh);
}


bool writeOBJ(const char* fileName, const triMeshD& mesh)
{
	return TRIANGLE_MESH::objWriteMesh(fileName, mesh);
}


bool writeOBJ(const char* fileName, const std::vector<verF>& vers, \
	const std::vector<edgeI>& edges)
{
	return TRIANGLE_MESH::objWriteEdges(fileName, vers, edges);
}


bool writeOBJ(const char* fileName, const std::vector<verD>& vers, \
	const std::vector<edgeI>& edges)
{
	return TRIANGLE_MESH::objWriteEdges(fileName, vers, edges);
}


bool readSTL(triMeshF& mesh, const char* fileName, const bool blIsAscii)
{
	return stlReadMesh(mesh, fileName, blIsAscii);
}


bool readSTL(triMeshD& mesh, const char* fileName, const bool blIsAscii)
{
	return stlReadMesh(mesh, fileName, blIsAscii);
}


bool writeSTL(const char* fileName, const triMeshF& mesh, const bool blCalcNorms)
{
	return stlWriteMesh(fileName, mesh, blCalcNorms);
}


bool writeSTL(const char* fileName, const triMeshD& mesh, const bool blCalcNorms)
{
	return stlWriteMesh(fileName, mesh, blCalcNorms);
}




////////////////////////////////////////////////////////////////////////////////////////////// other methods:
void tris2edges(std::vector<edgeI>& edges, const std::vector<triangleI>& tris)
{
	TRIANGLE_MESH::tris2edges(edges, tris);
}