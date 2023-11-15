#include "triMesh.h" 
#include <iostream>
#include <fstream>
#include <string>
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


	template <typename TV, typename TI>
	void GetVertsAndSurfs(std::vector<triplet<TV>>& vVerts, std::vector<triplet<TI>>& vSurfs, \
		const std::vector<triplet<TV>>& versIn)
	{
		auto vertHash = [](const triplet<TV>& v)->std::size_t
		{
			return std::hash<decltype(v.x)>()(v.x)	+ \
				std::hash<decltype(v.y)>()(v.y) + \
				std::hash<decltype(v.z)>()(v.z);
		};

		auto vertComp = [](const triplet<TV>& v0, const triplet<TV>& v1)->bool
		{
			double dbThreshold = 1.e-14;
			if (std::fabs(v0.x - v1.x) > dbThreshold)
				return false;
			if (std::fabs(v0.y - v1.y) > dbThreshold)
				return false;
			if (std::fabs(v0.z - v1.z) > dbThreshold)
				return false;
			return true;
		};

		// �޸�˵����ԭ��ȥ���ظ���ķ���ʱ�临�Ӷȹ��ߣ�����hashmap
		unsigned nOldVertCnt = versIn.size();
		std::vector<unsigned> tmpTri(3, 0);
		vSurfs.resize(nOldVertCnt / 3);
		vVerts.reserve(nOldVertCnt);
		std::unordered_map<triplet<TV>, unsigned, decltype(vertHash), decltype(vertComp)> mapVerts;

		for (unsigned i = 0; i < nOldVertCnt / 3; i++)
		{
			unsigned nVCnt = 0;

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

			// ������Ϣ		
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
			sprintf_s(szBuf, 1024, "v %.17g %.17g %.17g", x0, y0, z0);
			dstFile << szBuf << "\n";
		}
		dstFile.close();

		return true;
	}


	// OBJ�ļ���ȡ����
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
			sprintf_s(szBuf, 1024, "v %.17g %.17g %.17g", ver.x, ver.y, ver.z);
			dstFile << szBuf << "\n";
		}

		for (const auto& tri : mesh.triangles)
		{
			char szBuf[256] = { 0 };
			sprintf_s(szBuf, 256, "f %d %d %d", tri.x + 1, tri.y + 1, tri.z + 1);
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
				if (std::strcmp(tempBuffer, "vertex") == 0)    //������Ϣ
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

			unsigned offset = 80;			// �����ʼ���ļ�ͷ�������ļ�����80���ֽڣ���
			unsigned nVertDataCount = *(std::uint32_t*)(buffer + offset);		// ����Ƭ������
			offset += 4;							// ������stl�ļ��У����궼��REAL32����������UINT32, ����4�ֽڣ�

			//�Ӷ������ļ���ȡ������Ϣ
			verV pt{ 0, 0, 0 };
			tmpVers.resize(nVertDataCount * 3);

			for (unsigned k = 0; k < nVertDataCount; k++)
			{
				offset += 4 * 3;					//normal

				for (unsigned i = 0; i < 3; i++)
				{
					pt.x = static_cast<TV>(*(float*)(buffer + offset));
					offset += 4;
					pt.y = static_cast<TV>(*(float*)(buffer + offset));
					offset += 4;
					pt.z = static_cast<TV>(*(float*)(buffer + offset));
					offset += 4;

					tmpVers[3 * k + i] = pt;
				}
				offset += 2;
			}
			delete(buffer);

			GetVertsAndSurfs(mesh.vertices, mesh.triangles, tmpVers);

			return true;
		}

		return true;
	}
}
using namespace TRIANGLE_MESH;



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


bool readSTL(triMeshF& mesh, const char* fileName, const bool blIsAscii)
{
	return stlReadMesh(mesh, fileName, blIsAscii);
}


bool readSTL(triMeshD& mesh, const char* fileName, const bool blIsAscii)
{
	return stlReadMesh(mesh, fileName, blIsAscii);
}
