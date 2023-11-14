#include "triMesh.h" 
#include <iostream>
#include <fstream>
#include <string>
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
			sprintf_s(szBuf, 1024, "v %.17g %.17g %.17g", x0, y0, z0);
			dstFile << szBuf << "\n";
		}
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

}


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
