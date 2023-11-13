#pragma once 
#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <fstream>
#include <queue>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
#include <initializer_list>
#include <memory>
#include <thread>
#include <limits>
#include <windows.h>

#include "Eigen/Dense"
#include "Eigen/Sparse"



///////////////////////////////////////////////////////////////////////////////////////////////////// auxiliary interface:
unsigned readNextData(char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize);


///////////////////////////////////////////////////////////////////////////////////////////////////// disp interface
template <typename Derived>
void dispMat(const Eigen::PlainObjectBase<Derived>& mat);


template <typename Derived>
void dispMatBlock(const Eigen::PlainObjectBase<Derived>& mat, const int rowStart, const int rowEnd, \
	const int colStart, const int colEnd);


// 打印稀疏矩阵中的所有非零元素：
template<typename spMat>				// 如果模板参数为T, 函数参数为Eigen::SparseMatrix<T>，编译会报错，不知道什么缘故； 
void dispSpMat(const spMat& sm, const unsigned startCol, const unsigned endCol, const unsigned showElemsCount = 0);

template<typename spMat>
void dispSpMat(const spMat& sm);

template<typename T, int N>
void dispVec(const Eigen::Matrix<T, N, 1>& vec);

template<typename T, int N>
void dispVec(const Eigen::Matrix<T, 1, N>& vec);

template<typename T, int N>
void dispVecSeg(const Eigen::Matrix<T, N, 1>& vec, const int start, const int end);

template<typename T, int N>
void dispVecSeg(const Eigen::Matrix<T, 1, N>& vec, const int start, const int end);


///////////////////////////////////////////////////////////////////////////////////////////////////// IO interface
template <typename T>
void vecWriteToFile(const char* fileName, const std::vector<T>& vec);
template<typename T>
void vecReadFromFile(std::vector<T>& vec, const char* fileName, const unsigned elemCount);

template<typename T>
bool matWriteToFile(const char* fileName, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat);
template<typename T>
bool matReadFromFile(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName);


// OBJ文件读取网格
template	<typename Scalar, typename Index>
void objReadMeshMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic>& tris, const char* fileName);


#if 0
template <typename Scalar>
bool stlReadMeshMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::MatrixXi& tris, const char* fileName, bool blAscii = false)
{
	if (blAscii)
	{
#if 0
		std::ifstream fin(m_strFileName, std::ios::in);

		fin.seekg(0, std::ios::end);			//seek to the end
		unsigned fileLen = (unsigned)fin.tellg();
		if (0 == fileLen)							 // file is empty 
			return false; 
		fin.seekg(0, std::ios::beg);			//seek to the beg

		char* pFileBuf = new char[fileLen + 1];
		std::memset(pFileBuf, 0, fileLen + 1);
		fin.read(pFileBuf, fileLen);

		char* pTemp = NULL;
		pTemp = pFileBuf;
		char tempBuffer[1024];
		unsigned nMaxSize = 1024;
		unsigned nReadLen = 0;
		unsigned nRet = 0;
		while (nReadLen < fileLen)
		{
			nRet = ReadNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
			if (0 == nRet) 
				break; 
			if (std::strcmp(tempBuffer, "vertex") == 0)			//顶点信息
			{
				VFVECTOR3 vert;
				nRet = ReadNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
				if (0 == nRet) 
					break; 
				vert.x = (float)atof(tempBuffer);
				nRet = ReadNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
				if (0 == nRet) 
					break;
				}
				vert.y = (float)atof(tempBuffer);
				nRet = ReadNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
				if (0 == nRet) 
					break; 
				vert.z = (float)atof(tempBuffer);
				mVerts.push_back(vert);
			}
		}
		delete(pFileBuf);
		GetVertsAndSurfs(vVerts, vSurfs);
#endif

		return true;
	}
	else
	{
		std::ifstream fin(m_strFileName, std::ios::in | std::ios::binary);

		fin.seekg(0, std::ios::end);   //seek to the end
		unsigned fileLen = (unsigned)fin.tellg();
		if (0 == fileLen)			 // file is empty 
			return false; 

		fin.seekg(0, std::ios::beg);
		unsigned len = fin.tellg();
		char* buffer = new char[fileLen + 1];
		std::memset(buffer, 0, fileLen + 1);
		fin.read(buffer, fileLen);

		unsigned offset = 80;
		unsigned nVertDataCount = *(unsigned*)(buffer + offset);		   //获取nVertDataCount
		offset += sizeof(int32_t);

		//从二进制文件读取顶点信息
		VFVECTOR3 pt = VFVECTOR3::ZERO;
		mVerts.resize(nVertDataCount * 3);

		for (unsigned k = 0; k < nVertDataCount; k++)
		{
			offset += 4 * 3; //normal

			for (unsigned i = 0; i < 3; i++)
			{
				pt.x = *(float*)(buffer + offset);
				offset += 4;
				pt.y = *(float*)(buffer + offset);
				offset += 4;
				pt.z = *(float*)(buffer + offset);
				offset += 4;

				mVerts[3 * k + i] = pt;
			}

			offset += 2;
		}
		delete(buffer);

		GetVertsAndSurfs(vVerts, vSurfs);

		return true;
	}
}
#endif


// 网格写入到OBJ文件
template <typename DerivedV>
void objWriteMeshMat(const char* fileName, const Eigen::MatrixBase<DerivedV>& vers, \
	const Eigen::MatrixXi& tris)
{
	if (vers.cols() != 3 || tris.cols() != 3)
		return;
	std::ofstream dstFile(fileName);
	if (false == dstFile.is_open())
		return;
	 
	for (int j = 0; j < vers.rows(); j++)
	{
		char szBuf[1024] = { 0 };
		sprintf_s(szBuf, 1024, "v %.17g %.17g %.17g", vers(j, 0), vers(j, 1), vers(j, 2));
		dstFile << szBuf << "\n";
	}

	for (unsigned j = 0; j < tris.rows(); ++j)
	{
		char szBuf[256] = { 0 };
		sprintf_s(szBuf, 256, "f %d %d %d", tris(j, 0) + 1, tris(j, 1) + 1, tris(j, 2) + 1);
		dstFile << szBuf << "\n";
	}
}

// OBJ文件读取顶点数据
template	<typename Scalar>
void objReadVerticesMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const char* fileName);


// 顶点写入OBJ文件
template	<typename DerivedV>
void objWriteVerticesMat(const char* fileName, const Eigen::PlainObjectBase<DerivedV>& vers);


// XOY平面上的二维顶点写入OBJ文件；vers为两列的矩阵； 
template	<typename DerivedV>
void objWriteVerticesMat2D(const char* fileName, const Eigen::PlainObjectBase<DerivedV>& vers);


template	<typename DerivedV, typename DerivedI>
void objWriteEdgesMat(const char* pathName, const Eigen::PlainObjectBase<DerivedI>& edges, \
	const Eigen::PlainObjectBase<DerivedV>& vers);


template <typename DerivedV, typename	 DerivedI>
void objWritePath(const char* pathName, const std::vector<DerivedI>& path, const Eigen::PlainObjectBase<DerivedV>& vers);


template <typename DerivedV>
void objWriteTreePath(const char* pathName, const Eigen::VectorXi& treeVec, const Eigen::PlainObjectBase<DerivedV>& vers);


void objWriteDirection(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& dir);


void objWriteCoorSys(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& xdir, \
	const Eigen::RowVector3f& ydir, const Eigen::RowVector3f& zdir);


#include "temp.tpp"
