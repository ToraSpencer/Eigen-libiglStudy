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
 
#include "myEigenBasicMath/myEigenBasicMath.h"
#pragma comment(lib, "myEigenBasicMath.lib")

///////////////////////////////////////////////////////////////////////////////////////////////////// auxiliary types: 
template <typename T>
struct DOUBLET
{
	T x;
	T y;
};

template <typename T>
struct TRIPLET
{
	T x;
	T y;
	T z;
}; 


///////////////////////////////////////////////////////////////////////////////////////////////////// auxiliary interface:
unsigned readNextData(char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize);


// 3D顶点或三角片矩阵转换为TRIPLET向量的形式；
template <typename Derived>
std::vector<TRIPLET<typename Derived::Scalar>> mat2triplets(const Eigen::PlainObjectBase<Derived>& mat)
{
	using Scalar = typename Derived::Scalar;
	using ts = TRIPLET<Scalar>;

	std::vector<ts> vec;
	if (0 == mat.rows() || 3 != mat.cols())
		return vec;

	vec.resize(mat.rows());
	for (unsigned i = 0; i < mat.rows(); ++i)
	{
		vec[i].x = mat(i, 0);
		vec[i].y = mat(i, 1);
		vec[i].z = mat(i, 2);
	}

	return vec;
}


template <typename Derived>
TRIPLET<typename Derived::Scalar> vec2triplet(const Eigen::PlainObjectBase<Derived>& vec)
{
	assert((3 == vec.rows() * vec.cols()) && "assert!!! input vec should be in 3D space.");
	using Scalar = typename Derived::Scalar;
	using ts = TRIPLET<Scalar>;
	return ts{ vec(0), vec(1), vec(2) };
}


template <typename Derived>
std::vector<DOUBLET<typename Derived::Scalar>> mat2doublets(const Eigen::PlainObjectBase<Derived>& mat)
{
	using Scalar = typename Derived::Scalar;
	using ds = DOUBLET<Scalar>;

	std::vector<ds> vec;
	if (0 == mat.rows() || 2 != mat.cols())
		return vec;

	vec.resize(mat.rows());
	for (unsigned i = 0; i < mat.rows(); ++i)
	{
		vec[i].x = mat(i, 0);
		vec[i].y = mat(i, 1);
	}

	return vec;
}


template <typename Derived>
DOUBLET<typename Derived::Scalar> vec2doublet(const Eigen::PlainObjectBase<Derived>& vec)
{
	assert((2 == vec.rows() * vec.cols()) && "assert!!! input vec should be in 2D space.");
	using Scalar = typename Derived::Scalar;
	using ds = DOUBLET<Scalar>;
	return DOUBLET{ vec(0), vec(1) };
}



///////////////////////////////////////////////////////////////////////////////////////////////////// disp interface
template <typename Derived>
void dispMat(const Eigen::PlainObjectBase<Derived>& mat);


template <typename Derived>
void dispMatBlock(const Eigen::PlainObjectBase<Derived>& mat,\
		const int rowStart, const int rowEnd, \
		const int colStart, const int colEnd);


// 打印稀疏矩阵中的所有非零元素：
template<typename spMat>				// 如果模板参数为T, 函数参数为Eigen::SparseMatrix<T>，编译会报错，不知道什么缘故； 
void dispSpMat(const spMat& sm, const unsigned startCol, \
		const unsigned endCol, const unsigned showElemsCount = 0);

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


// OBJ文件读取顶点数据
template	<typename Scalar>
void objReadVerticesMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
		const char* fileName);


// STL文件中读取网格数据：
#if 0
template <typename DerivedV, typename DerivedI>
void stlReadMeshMat(Eigen::PlainObjectBase<DerivedV>& vers, \
	Eigen::PlainObjectBase<DerivedI>& tris, const char* fileName, const bool blIsAscii = false)
{
	using ScalarV = typename DerivedV::Scalar;
	using ScalarI = typename DerivedI::Scalar;
	using RowVector3V = Eigen::MatrixX<ScalarV, Eigen::Dynamic, Eigen::Dynamic>;
	vers.resize(0, 0);
	tris.resize(0, 0);
	if (blIsAscii)
	{
		std::ifstream fin(fileName, std::ios::in);

		fin.seekg(0, std::ios::end);			 //seek to the end
		unsigned fileLen = (unsigned)fin.tellg();
		if (0 == fileLen)							 // file is empty
			return false;

		fin.seekg(0, std::ios::beg);		//seek to the beg

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
			nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
			if (0 == nRet) 
				break; 
			if (std::strcmp(tempBuffer, "vertex") == 0)    //顶点信息
			{
				RowVector3V vert;
				nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
				if (0 == nRet) 
					break; 
				vert(0) = static_cast<ScalarV>(atof(tempBuffer));
				nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
				if (0 == nRet) 
					break; 
				vert(1) = static_cast<ScalarV>(atof(tempBuffer));
				nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
				if (0 == nRet) 
					break; 
				vert(2) = static_cast<ScalarV>(atof(tempBuffer));
				matInsertRows(vers, vert); 
			}
		}
		delete(pFileBuf);
		 
		return true;
	}
	else
	{
		std::ifstream fin(fileName, std::ios::in | std::ios::binary);

		fin.seekg(0, std::ios::end);   //seek to the end
		unsigned fileLen = (unsigned)fin.tellg();
		if (0 == fileLen)		 // file is empty 
			return false; 

		fin.seekg(0, std::ios::beg);
		unsigned len = fin.tellg();
		char* buffer = new char[fileLen + 1];
		std::memset(buffer, 0, fileLen + 1);
		fin.read(buffer, fileLen);

		unsigned offset = 80;
		unsigned nVertDataCount = *(unsigned*)(buffer + offset);   //获取nVertDataCount
		offset += sizeof(int32_t);

		//从二进制文件读取顶点信息
		VFVECTOR3 pt = VFVECTOR3::ZERO;

		mVerts.resize(nVertDataCount * 3);

		for (unsigned k = 0; k < nVertDataCount; k++)
		{
			offset += 4 * 3;										//normal

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


// 顶点写入OBJ文件
template	<typename DerivedV>
void objWriteVerticesMat(const char* fileName, \
		const Eigen::MatrixBase<DerivedV>& vers)
{
	assert(3 == vers.cols() && "Assert!!! vers mat column size should be 3.");
	std::ofstream dstFile(fileName);
	for (int j = 0; j < vers.rows(); j++)
	{
		char szBuf[1024] = { 0 };
		sprintf_s(szBuf, 1024, "v %.17g %.17g %.17g", vers(j, 0), vers(j, 1), vers(j, 2));
		dstFile << szBuf << "\n";
	}
	dstFile.close();
}


// 网格写入到OBJ文件
template <typename DerivedV>
void objWriteMeshMat(const char* fileName, \
		const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
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
 

// XOY平面上的二维顶点写入OBJ文件；vers为两列的矩阵； 
template	<typename DerivedV>
void objWriteVerticesMat2D(const char* fileName, \
		const Eigen::MatrixBase<DerivedV>& vers)
{
	assert(2 == vers.cols(), "error!!! 2D vers mat column size should be 2.");
	std::ofstream dstFile(fileName);
	for (int j = 0; j < vers.rows(); j++)
	{
		char szBuf[1024] = { 0 };
		sprintf_s(szBuf, 1024, "v %.17g %.17g %.17g", vers(j, 0), vers(j, 1), 0);
		dstFile << szBuf << "\n";
	}
	dstFile.close();
}


// 边数据写入到OBJ文件中：
template	<typename DerivedV, typename DerivedI>
void objWriteEdgesMat(const char* pathName, \
		const Eigen::MatrixBase<DerivedI>& edges, \
		const Eigen::MatrixBase<DerivedV>& vers)
{
	if (0 == edges.rows() || 0 == vers.rows())
		return;
	std::ofstream dstFile(pathName);
	for (int i = 0; i < vers.rows(); ++i)
		dstFile << "v " << vers(i, 0) << " " << vers(i, 1) << " " << vers(i, 2) << std::endl;

	for (int i = 0; i < edges.rows(); ++i)
		dstFile << "l " << edges(i, 0) + static_cast<typename DerivedI::Scalar>(1) << " " \
		<< edges(i, 1) + static_cast<typename DerivedI::Scalar>(1) << std::endl;

	dstFile.close();
}


// objWritePath() 路径数据写入到OBJ文件中：
template <typename DerivedV, typename	 IndexType>
void objWritePath(const char* pathName, const std::vector<IndexType>& path, \
	const Eigen::MatrixBase<DerivedV>& vers)
{
	if (path.size() <= 1)
		return;

	unsigned edgesCount = path.size() - 1;
	Eigen::MatrixXi pathEdges(edgesCount, 2);

	for (unsigned i = 0; i < edgesCount; ++i)
		pathEdges(i, 0) = static_cast<int>(path[i]);
	for (unsigned i = 0; i < edgesCount; ++i)
		pathEdges(i, 1) = static_cast<int>(path[i + 1]);

	objWriteEdgesMat(pathName, pathEdges, vers);
}


// objWirteTreePath()——输入树向量或路径向量，写入到OBJ文件中：
template <typename DerivedV>
void objWriteTreePath(const char* pathName, const Eigen::VectorXi& treeVec, \
		const Eigen::MatrixBase<DerivedV>& vers)
{
	// 路径被视为树的特例；

	// 树中索引为i的顶点的前驱顶点索引为treeVec(i), 若其没有前驱顶点（父节点），则treeVec(i) == -1;
	if (treeVec.size() <= 1)
		return;

	unsigned edgesCount = 0;
	for (int i = 0; i < treeVec.rows(); ++i)
		if (treeVec(i) >= 0)
			edgesCount++;

	Eigen::MatrixXi edges(edgesCount, 2);
	int rowIdx = 0;
	for (int i = 0; i < treeVec.rows(); ++i)
	{
		if (treeVec(i) >= 0)
		{
			edges(rowIdx, 0) = treeVec(i);
			edges(rowIdx, 1) = i;
			rowIdx++;
		}
	}

	objWriteEdgesMat(pathName, edges, vers);
}


template <typename DerivedVo, typename DerivedVd>
void objWriteDirection(const char* pathName, const Eigen::MatrixBase<DerivedVo>& origin, \
		const Eigen::MatrixBase<DerivedVd>& dir, const float SR = 0.5)
{
	assert((3 == origin.size() && 3 == dir.size()));

	Eigen::MatrixXf line;
	Eigen::RowVector3f dirVec{static_cast<float>(dir(0)), static_cast<float>(dir(1)), static_cast<float>(dir(2))};
	dirVec.normalize();

	const float length = 10;
	int versCount = std::round(length / SR);				// 空间采样率SR——相邻两个采样点的距离（单位mm）
	line.resize(versCount + 1, 3);
	line(0, 0) = static_cast<float>(origin(0));
	line(0, 1) = static_cast<float>(origin(1));
	line(0, 2) = static_cast<float>(origin(2));
	for (int i = 1; i <= versCount; i++)
		line.row(i) = line.row(0) + SR * dirVec * i;

	objWriteVerticesMat(pathName, line);
};


template <typename DerivedVo, typename DerivedVd>
void objWriteCoorSys(const char* pathName, const Eigen::MatrixBase<DerivedVo>& origin, \
		const Eigen::MatrixBase<DerivedVd>& xdir, const Eigen::MatrixBase<DerivedVd>& ydir, \
		const Eigen::MatrixBase<DerivedVd>& zdir)
{
	const float SR = 0.5;			// 空间采样率SR——相邻两个采样点的距离（单位mm）
	const float length = 10;
	int versCount = std::round(length / SR);
	Eigen::RowVector3f originVec{ static_cast<float>(origin(0)), static_cast<float>(origin(1)), static_cast<float>(origin(2)) };
	Eigen::RowVector3f xdirVec{ static_cast<float>(xdir(0)), static_cast<float>(xdir(1)), static_cast<float>(xdir(2)) };
	Eigen::RowVector3f ydirVec{ static_cast<float>(ydir(0)), static_cast<float>(ydir(1)), static_cast<float>(ydir(2)) };
	Eigen::RowVector3f zdirVec{ static_cast<float>(zdir(0)), static_cast<float>(zdir(1)), static_cast<float>(zdir(2)) };
	Eigen::MatrixXf line1(versCount, 3), line2(versCount, 3), line3(versCount, 3);
	for (int i = 0; i < versCount; i++)
		line1.row(i) = originVec + SR * xdirVec * (i + 1);
	for (int i = 0; i < versCount; i++)
		line2.row(i) = originVec + SR * ydirVec * (i + 1);
	for (int i = 0; i < versCount; i++)
		line3.row(i) = originVec + SR * zdirVec * (i + 1);

	Eigen::MatrixXf line = originVec;
	int rows = line.rows();
	line.conservativeResize(line.rows() + line1.rows() + line2.rows() + line3.rows(), 3);
	for (int i = 0; i < line1.rows(); ++i)
		line.row(rows + i) = line1.row(i);
	rows = line.rows();
	for (int i = 0; i < line2.rows(); ++i)
		line.row(rows + i) = line2.row(i);
	rows = line.rows();
	for (int i = 0; i < line3.rows(); ++i)
		line.row(rows + i) = line3.row(i);

	objWriteVerticesMat(pathName, line);
}



#if 0			// 当前有问题
template <typename Scalar>
bool stlReadMeshMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::MatrixXi& tris, const char* fileName, bool blAscii = false)
{
	if (blAscii)
	{
#if 0
		std::ifstream fin(fileName, std::ios::in);

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
	std::ifstream fin(fileName, std::ios::in | std::ios::binary);

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

#include "temp.tpp"
