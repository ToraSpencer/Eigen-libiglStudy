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
#define WIN32_LEAN_AND_MEAN             // 从 Windows 头文件中排除极少使用的内容
 

///////////////////////////////////////////////////////////////////////////////////////////////////////// IO 

///////////////////////////////////////////////////////////////////////////////////////////////////////// 表象变换

// 笛卡尔坐标→齐次坐标系
template <typename DerivedVo, typename DerivedVi>
bool vers2HomoVers(Eigen::PlainObjectBase<DerivedVo>& versOut, const Eigen::MatrixBase<DerivedVi>& versIn)
{
	using ScalarI = typename DerivedVi::Scalar;
	using ScalarO = typename DerivedVo::Scalar;

	const int dim = versIn.cols();
	const int versCount = versIn.rows();
	assert(2 == dim || 3 == dim, "assert!!! input vertices dimension should be 2 or 3.");

	versOut.resize(0, 0);
	if (0 == versCount)
		return true;

	versOut.resize(dim + 1, versCount);
	versOut.setOnes();
	versOut.topRows(dim) = versIn.transpose().array().cast<ScalarO>();

	return true;
}


// 齐次坐标系→笛卡尔坐标系
template <typename DerivedVo, typename DerivedVi>
bool homoVers2Vers(Eigen::PlainObjectBase<DerivedVo>& versOut, const Eigen::MatrixBase<DerivedVi>& versIn)
{
	using ScalarI = typename DerivedVi::Scalar;
	using ScalarO = typename DerivedVo::Scalar;

	const int dim = versIn.rows() - 1;
	const int versCount = versIn.cols();
	assert(2 == dim || 3 == dim, "assert!!! input vertices dimension should be 2 or 3.");

	versOut.resize(0, 0);
	if (0 == versCount)
		return true;

	versOut = versIn.transpose().leftCols(dim).array().cast<ScalarO>();

	return true;
}


template <typename DerivedV>
Eigen::MatrixXf vers2HomoVersF(const Eigen::MatrixBase<DerivedV>& versIn)
{
	Eigen::MatrixXf versOut;
	vers2HomoVers(versOut, versIn);
	return versOut;
}


template <typename DerivedV>
Eigen::MatrixXd vers2HomoVersD(const Eigen::MatrixBase<DerivedV>& versIn)
{
	Eigen::MatrixXd versOut;
	vers2HomoVers(versOut, versIn);
	return versOut;
}


template <typename DerivedV>
Eigen::MatrixXf homoVers2VersF(const Eigen::MatrixBase<DerivedV>& versIn)
{
	Eigen::MatrixXf versOut;
	homoVers2Vers(versOut, versIn);

	return versOut;
}


template <typename DerivedV>
Eigen::MatrixXf homoVers2VersD(const Eigen::MatrixBase<DerivedV>& versIn)
{
	Eigen::MatrixXd versOut;
	homoVers2Vers(versOut, versIn);

	return versOut;
}
 

 
 
