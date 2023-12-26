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


///////////////////////////////////////////////////////////////////////////////////////////////////// auxiliary interface:
unsigned readNextData(char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize);


// 3D���������Ƭ����ת��ΪTRIPLET��������ʽ��
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


// ��ӡϡ������е����з���Ԫ�أ�
template<typename spMat>				// ���ģ�����ΪT, ��������ΪEigen::SparseMatrix<T>������ᱨ����֪��ʲôԵ�ʣ� 
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


// OBJ�ļ���ȡ����
template	<typename Scalar, typename Index>
void objReadMeshMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
		Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic>& tris, const char* fileName);


// OBJ�ļ���ȡ��������
template	<typename Scalar>
void objReadVerticesMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
		const char* fileName);
 

template	<typename DerivedVO, typename DerivedVI>
void GetVertsAndSurfs(Eigen::PlainObjectBase<DerivedVO>& versOut, \
	Eigen::MatrixXi& trisOut, const Eigen::MatrixBase<DerivedVI>& versIn)
{
	using ScalarO = typename DerivedVO::Scalar;

	auto vertHash = [](const Eigen::RowVector3f& v)->std::size_t
	{
		return std::hash<float>()(v(0)) + \
			std::hash<float>()(v(1)) + \
			std::hash<float>()(v(2));
	};

	auto vertComp = [](const Eigen::RowVector3f& v0, const Eigen::RowVector3f& v1)->bool
	{
		double dbThreshold = 1.e-14;
		if (std::fabs(v0(0) - v1(0)) > dbThreshold)
			return false;
		if (std::fabs(v0(1) - v1(1)) > dbThreshold)
			return false;
		if (std::fabs(v0(2) - v1(2)) > dbThreshold)
			return false;
		return true;
	};

	versOut.resize(0, 0);
	trisOut.resize(0, 0);

	// �޸�˵����ԭ��ȥ���ظ���ķ���ʱ�临�Ӷȹ��ߣ�����hashmap
	unsigned versCountOld = versIn.rows();
	unsigned trisCount = versCountOld / 3;
	std::vector<unsigned> tmpTri(3, 0);
	trisOut.resize(trisCount, 3);
	versOut.resize(versCountOld, 3);
	std::unordered_map<Eigen::RowVector3f, unsigned, decltype(vertHash), decltype(vertComp)> mapVerts;

	unsigned vIdx = 0;
	unsigned nOldIdx = 0;
	for (unsigned i = 0; i < trisCount; i++)
	{ 
		for (unsigned k = 0; k < 3; k++)
		{
			nOldIdx = i * 3 + k;
			Eigen::RowVector3f v = versIn.row(nOldIdx).cast<float>();
			if (0 == mapVerts.count(v))
			{
				mapVerts.insert(std::make_pair(v, vIdx));
				versOut.row(vIdx) = v.cast<ScalarO>();
				vIdx++;
			}
			auto iter = mapVerts.find(v);
			tmpTri[k] = iter->second;
		}
 
		trisOut(i, 0) = static_cast<int>(tmpTri[0]);
		trisOut(i, 1) = static_cast<int>(tmpTri[1]);
		trisOut(i, 2) = static_cast<int>(tmpTri[2]);
	}

	// shrink to fit��
	versOut.conservativeResize(vIdx, 3);
}


// STL�ļ��ж�ȡ�������ݣ�
template	<typename DerivedV>
bool stlReadMeshMat(Eigen::PlainObjectBase<DerivedV>& meshVers, Eigen::MatrixXi& meshTris,\
	const char* fileName, const bool blIsAscii = false)
{
	using ScalarV = typename DerivedV::Scalar;
	int versCount = 0;
	int index = 0;
	meshVers.resize(0, 0);
	meshTris.resize(0, 0);

	if (blIsAscii)
	{
		Eigen::MatrixXd tmpVers;
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
			index = versCount;
			versCount++;
			tmpVers.conservativeResize(versCount, 3);
			nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
			if (0 == nRet)
				break;
			if (std::strcmp(tempBuffer, "vertex") == 0)    //������Ϣ
			{
				nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
				if (0 == nRet)
					break;
				tmpVers(index, 0) = atof(tempBuffer);
				nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
				if (0 == nRet)
					break;
				tmpVers(index, 1) = atof(tempBuffer);
				nRet = readNextData(pTemp, nReadLen, tempBuffer, nMaxSize);
				if (0 == nRet)
					break;
				tmpVers(index, 2) = atof(tempBuffer);
			}
		}
		delete(pFileBuf);

		GetVertsAndSurfs(meshVers, meshTris, tmpVers);

		return true;
	}
	else     // �����Ƹ�ʽ��
	{
		Eigen::MatrixXf tmpVers;
		std::ifstream fin(fileName, std::ios::in | std::ios::binary);

		fin.seekg(0, std::ios::end);	 // seek to the end
		unsigned fileLen = (unsigned)fin.tellg();
		if (0 == fileLen)					// file is empty 
			return false;

		fin.seekg(0, std::ios::beg);
		unsigned len = fin.tellg();
		char* buffer = new char[fileLen + 1];
		std::memset(buffer, 0, fileLen + 1);
		fin.read(buffer, fileLen);

		unsigned offset = 80;			// �����ʼ���ļ�ͷ�������ļ�����80���ֽڣ���
		unsigned trisCount = *(std::uint32_t*)(buffer + offset);		// ����Ƭ������
		offset += 4;							// ������stl�ļ��У����궼��REAL32����������UINT32, ����4�ֽڣ�

		//�Ӷ������ļ���ȡ������Ϣ 
		tmpVers.resize(trisCount * 3, 3);
		for (unsigned i = 0; i < trisCount; i++)
		{
			offset += 4 * 3;								// ����������Ϣ
			for (unsigned k = 0; k < 3; k++)				// ���δ洢ÿ������Ƭ�������������ꣻ
			{
				index = 3 * i + k;
				tmpVers(index, 0) = *(float*)(buffer + offset);
				offset += 4;
				tmpVers(index, 1) = *(float*)(buffer + offset);
				offset += 4;
				tmpVers(index, 2) = *(float*)(buffer + offset);
				offset += 4;
			}
			offset += 2;
		}
		delete(buffer);

		GetVertsAndSurfs(meshVers, meshTris, tmpVers);

		return true;
	}

	return true;
}


// ����д��OBJ�ļ�
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


// ����д�뵽OBJ�ļ�
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
 

// XOYƽ���ϵĶ�ά����д��OBJ�ļ���versΪ���еľ��� 
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


// ������д�뵽OBJ�ļ��У�
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


// objWritePath() ·������д�뵽OBJ�ļ��У�
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


// objWirteTreePath()����������������·��������д�뵽OBJ�ļ��У�
template <typename DerivedV>
void objWriteTreePath(const char* pathName, const Eigen::VectorXi& treeVec, \
		const Eigen::MatrixBase<DerivedV>& vers)
{
	// ·������Ϊ����������

	// ��������Ϊi�Ķ����ǰ����������ΪtreeVec(i), ����û��ǰ�����㣨���ڵ㣩����treeVec(i) == -1;
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
	int versCount = std::round(length / SR);				// �ռ������SR������������������ľ��루��λmm��
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
	const float SR = 0.5;			// �ռ������SR������������������ľ��루��λmm��
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


 

#include "temp.tpp"
