#ifndef _SMARTEE_DLL_H
#define _SMARTEE_DLL_H

#include <SDKDDKVer.h>	
#include <windows.h>					// Windows 头文件: 

#include <iostream>


// 使用eigen库中的向量、矩阵表示顶点、点云、面片等信息；
#include "Eigen/Dense"


#ifdef SMARTEE_DLL_EXPORTS
#define SMARTEE_DLL_API __declspec(dllexport) 
#else
#define SMARTEE_DLL_API __declspec(dllimport)
#endif


// 导出函数；注！！！——函数实现开头也要加上宏SMARTEE_DLL_API：

// 从OBJ文件中读点云
template	<typename Scalar>
SMARTEE_DLL_API void readOBJ(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, const char* fileName);

// 从OBJ文件中读三角网格
template	<typename Scalar>
SMARTEE_DLL_API void readOBJ(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::MatrixXi& tris, const char* fileName);

// 点云写入到OBJ文件：
template<typename DerivedV>
SMARTEE_DLL_API void writeOBJ(const char* fileName, const Eigen::PlainObjectBase<DerivedV>& vers);

// 三角网格写入到OBJ文件：
template<typename DerivedV>
SMARTEE_DLL_API void writeOBJ(const char* fileName, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::MatrixXi& tris);

// 老版的B样条插值
SMARTEE_DLL_API bool interpltnCubicBSpline(Eigen::MatrixXf& versOut, const Eigen::MatrixXf& versIn, \
	unsigned numCurvePts, bool closeFlag);

// 新版三次B样条插值拟合曲线：
SMARTEE_DLL_API bool BsplineInterp(Eigen::MatrixXf& versOut, const Eigen::MatrixXf& versIn, \
	const Eigen::RowVector3f& planeVer, const Eigen::RowVector3f& planeNorm);

// 射线求交：
struct meshRayOut
{
	std::vector<std::vector<float>> rayLen;
	std::vector<std::vector<float>> rayOppLen;
	std::vector<std::vector<unsigned>> isctTris;
	std::vector<std::vector<unsigned>> isctOppTris;
};
template <typename DerivedV1, typename DerivedV2, typename Scalar>
SMARTEE_DLL_API bool meshRayIntersect(meshRayOut& result, \
	const Eigen::PlainObjectBase<DerivedV1>& srcVers, const Eigen::Matrix<Scalar, 1, 3>& dir, \
	const Eigen::PlainObjectBase<DerivedV2>& meshVers, const Eigen::MatrixXi& meshTris);


#endif

