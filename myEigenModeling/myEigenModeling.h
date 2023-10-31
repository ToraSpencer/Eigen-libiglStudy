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

#include "myEigenHomoCoor/myEigenHomoCoor.h"
#pragma comment(lib, "myEigenHomoCoor.lib")

#include "myEigenIO/myEigenIO.h"
#pragma comment(lib, "myEigenIO.lib")

#include "myEigenPMP/myEigenPMP.h"
#pragma comment(lib, "myEigenPMP.lib")


//  64位环境下，triangulate.cpp里" x = vertexloop[0] = pointlist[coordindex++];"这一句会抛出异常，需要研究
#ifndef _WIN64
#define USE_TRIANGLE_H
#endif

#ifdef  USE_TRIANGLE_H
#include "triangulate.h"
#endif
 

///////////////////////////////////////////////////////////////////////////////////////////////////// modeling接口：

// 生成中心在原点，边长为1，三角片数为12的正方体网格；
template	<typename T>
void genCubeMesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);


// 生成轴向包围盒的三角网格；
template <typename T>
void genAABBmesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<T, 3>& aabb);


// 生成栅格采样点：
template<typename Tg, typename To>
bool genGrids(Eigen::Matrix<Tg, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, \
	const Eigen::Matrix<To, 1, 3>& origin, \
	const float step, const std::vector<unsigned>& gridCounts);


#ifdef USE_TRIANGLE_H

// 重载1：2D点云三角剖分得到面网格——可以带洞也可以不带洞
/*

	注：
		输入点云必须在XOY平面内，可以是2D的也可以是3D的；
		默认三角剖分模式为"pY"，表示三角剖分成不插点的平面直线图；


		switch string:
				-p			三角剖分生成一个平面直线图
				-r			(refine a previously generated mesh)对一个已有网格进行进一步的三角剖分；
				-q			(Quality mesh generation)后面跟一个数值，如"-q30"表示三角剖分结果中不可以存在小于30°的角；
				-a			后面跟一个数值，如"-a5"指定三角剖分结果中三角片面积不大于5mm^2;
				-Y			(Prohibits the insertion of Steiner points on the mesh boundary.)
							禁止在边缘边上插入新点；
				-YY		prohibits the insertion of Steiner points on any segment, including internal segments.
							禁止在任何原有边上插入新点；

*/
template <typename DerivedVo, typename DerivedI, typename DerivedVi, typename DerivedVh>
bool triangulateVers2Mesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const std::vector<Eigen::VectorXi>& bdryLoops, \
	const Eigen::PlainObjectBase<DerivedVh>& holeCenters, \
	const char* strSwitcher = "pY");


// 重载2：2D点云三角剖分得到面网格——不带洞
template <typename DerivedVo, typename DerivedI, typename DerivedVi>
bool triangulateVers2Mesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const std::vector<Eigen::VectorXi>& bdryLoops, \
	const char* strSwitcher = "pY");


// 三角剖分提升网格质量：
template <typename DerivedVo, typename DerivedIt, typename DerivedVi, typename DerivedIe>
bool triangulateRefineMesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedIt>& trisOut, const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const Eigen::PlainObjectBase<DerivedIt>& trisIn, const Eigen::PlainObjectBase<DerivedIe>& edges, \
	const char* strSwitcher);



//	重载1：封闭边界线点集得到面网格，可以是平面也可以是曲面，不在网格内部插点，三角片尺寸不可控。
template <typename DerivedO, typename DerivedC>
bool circuit2mesh(Eigen::PlainObjectBase<DerivedO>& vers, Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedC>& circVers);


// 重载2：会在回路内部插点，三角片面积上限由参数maxTriArea决定；
template <typename DerivedV, typename DerivedL, typename DerivedN>
bool circuit2mesh(Eigen::PlainObjectBase<DerivedV>& vers, Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedL>& versLoop, \
	const Eigen::PlainObjectBase<DerivedN>& normDir, const float maxTriArea);


// genCylinder()重载1——生成（类）柱体：
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers, \
	const bool isCovered);


// genCylinder()重载2——生成圆柱体网格
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers,\
	const float radius, const double deltaTheta, const bool isCovered);

// genCylinder()重载3——生成方柱体网格
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers,\
	const std::pair<float, float> sizePair, const float SR, const bool isCovered);


// genCylinder()重载4——输入上底面和下底面的2D边界环路 ，生成柱体：
/*
	bool genCylinder(
			MatrixXT& vers,								输出网格顶点
			Eigen::MatrixXi& tris,						输出网格三角片
			const Eigen::PlainObjectBase<DerivedVa>& axisVers,			柱体轴线
			const Eigen::PlainObjectBase<DerivedVt>& topLoop,			上底面边界环路
			const Eigen::PlainObjectBase<DerivedVb>& btmLoop,			下底面边界环路
			const bool isCovered																是否封底
			)
	注：上下底面边界环路都是在XOY平面中的点集，
					且两者顶点数需要相同；
					且从Z轴正向看环路顶点索引应该沿着顺时针方向增大；

*/
template <typename T, typename DerivedVa, typename DerivedVt, typename DerivedVb>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedVa>& axisVers, const Eigen::PlainObjectBase<DerivedVt>& topLoop, \
	const Eigen::PlainObjectBase<DerivedVb>& btmLoop, const bool isCovered);


// genCylinder()重载5——输入上底面和下底面的3D边界环路 ，生成柱体：
/*
	bool genCylinder(
			Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers,		输出网格点云
			Eigen::MatrixXi& tris,																	输出网格三角片
			const Eigen::PlainObjectBase<DerivedVt>& topLoop,				上底面边界回路
			const Eigen::PlainObjectBase<DerivedVb>& btmLoop,				下表面边界回路
			const Eigen::Matrix<ScalarN, 1, 3>& topNorm							上底面法向
			const Eigen::Matrix<ScalarN, 1, 3>& btmNorm							下底面法向
			const int layersCount,																	层数
			const float maxArea																	上下底面三角剖分的最大三角片面积
			const bool isCovered																	是否封底
			)
	注：要求上下底面边界环路顶点数需要相同；
				要求上底面边界回路生长方向的右手螺旋方向平行于上底面法向；下底面相反；
				即当柱体为竖直的圆柱时，上下底面回路右手螺旋方向都为竖直向上；
			layersCount决定了柱体被回路顶点分成了多少层
				最小为1，即在上下底面之间不插点
				上下底面之间插的回路顶点圈数为(layersCount - 1);
*/
template <typename DerivedVo, typename DerivedVt, typename DerivedVb, typename ScalarN>
bool genCylinder(Eigen::PlainObjectBase<DerivedVo>& vers, Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedVt>& topLoop, \
	const Eigen::PlainObjectBase<DerivedVb>& btmLoop, \
	const Eigen::Matrix<ScalarN, 1, 3>& topNorm, \
	const Eigen::Matrix<ScalarN, 1, 3>& btmNorm, \
	const int layersCount, const float maxArea, const bool isCovered);
 

//		生成方柱，旋转分两次，以确保侧面和XOY平面平行或垂直；
template <typename T>
bool genAlignedCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const std::pair<float, float> sizePair, const float SR, \
	const bool isCovered);

#endif



// 暂时未整理的实现：
#include "temp.tpp"