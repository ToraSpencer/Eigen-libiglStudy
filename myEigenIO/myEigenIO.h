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
template	<typename Scalar, typename Index>

// OBJ文件读取网格
void objReadMeshMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::Matrix<Index, Eigen::Dynamic, \
	Eigen::Dynamic>& tris, const char* fileName);

// 网格写入到OBJ文件
template	<typename T>
void objWriteMeshMat(const char* fileName, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris);

// OBJ文件读取顶点数据
template	<typename Scalar>
void objReadVerticesMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, const char* fileName);

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
