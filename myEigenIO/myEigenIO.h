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

unsigned readNextData(char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize);

template <typename T>
void vecWriteToFile(const char* fileName, const std::vector<T>& vec);
template<typename T>
void vecReadFromFile(std::vector<T>& vec, const char* fileName, const unsigned elemCount);

template<typename T>
bool matWriteToFile(const char* fileName, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat);
template<typename T>
bool matReadFromFile(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName);
template	<typename Scalar, typename Index>
void objReadMeshMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::Matrix<Index, Eigen::Dynamic, \
	Eigen::Dynamic>& tris, const char* fileName);
template	<typename T>
void objWriteMeshMat(const char* fileName, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris);
template	<typename Scalar>
void objReadVerticesMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, const char* fileName);
template	<typename DerivedV>
void objWriteVerticesMat(const char* fileName, const Eigen::PlainObjectBase<DerivedV>& vers);
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
