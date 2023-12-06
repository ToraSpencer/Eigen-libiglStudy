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

#include "MC_tables.h"

#include "myEigenBasicMath/myEigenBasicMath.h"
#pragma comment(lib, "myEigenBasicMath.lib")

#include "myEigenHomoCoor/myEigenHomoCoor.h"
#pragma comment(lib, "myEigenHomoCoor.lib")

#include "myEigenIO/myEigenIO.h"
#pragma comment(lib, "myEigenIO.lib")

#include "myEigenPMP/myEigenPMP.h"
#pragma comment(lib, "myEigenPMP.lib")


#define USE_TRIANGLE_H
 

#ifdef  _WIN64
#define TRI_REAL double
#else
#define TRI_REAL float
#endif  

#ifdef  USE_TRIANGLE_H
#include "triangle/triangle.h"
#pragma comment(lib, "triangle.lib")
#endif
  

///////////////////////////////////////////////////////////////////////////////////////////////////// modeling接口：

// 生成中心在原点，边长为1，三角片数为12的正方体网格；
template	<typename T>
void genCubeMesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);


// 生成轴向包围盒的三角网格；
template <typename T>
void genAABBmesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<T, 3>& aabb);


// 生成栅格采样点云
template<typename DerivedV, typename DerivedVO>
bool genGrids(Eigen::PlainObjectBase<DerivedV>& gridCenters, \
	const Eigen::MatrixBase<DerivedVO>& origin, \
	const float step, const std::vector<unsigned>& gridCounts)
{
	/*
		bool genGrids(
			gridCenters,											 输出的栅格点云
			origin,													 栅格原点，即三轴坐标最小的那个栅格点；
			step,														 采样步长
			gridCounts											 三元数组，分别是XYZ三轴上的步数；
			)
	*/
	using ScalarT = typename DerivedV::Scalar;
	using VectorXT = Eigen::Matrix<ScalarT, Eigen::Dynamic, 1>;
	using RowVector3T = Eigen::Matrix<ScalarT, 1, 3>;
	using MatrixXT = Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>;
	assert((3 == origin.size()) && "Assert!!! origin should be a 3D vertex.");
	gridCenters.resize(0, 0);

	// 整个距离场的包围盒：
	RowVector3T minp{ static_cast<ScalarT>(origin(0)), static_cast<ScalarT>(origin(1)), static_cast<ScalarT>(origin(2)) };
	RowVector3T maxp = minp + static_cast<ScalarT>(step) * RowVector3T(gridCounts[0] - 1, gridCounts[1] - 1, gridCounts[2] - 1);

	// 生成栅格：
	/*
		按索引增大排列的栅格中心点为：
		gc(000), gc(100), gc(200), gc(300),...... gc(010), gc(110), gc(210), gc(310),...... gc(001), gc(101), gc(201).....

		x坐标：
		x0, x1, x2, x3......x0, x1, x2, x3......x0, x1, x2, x3......
		周期为xCount;
		重复次数为(yCount * zCount)

		y坐标：
		y0, y0, y0...y1, y1, y1...y2, y2, y2.........y0, y0, y0...y1, y1, y1...
		周期为(xCount * yCount);
		重复次数为zCount;
		单个元素重复次数为xCount

		z坐标：
		z0, z0, z0......z1, z1, z1......z2, z2, z2......
		单个元素重复次数为(xCount * yCount)
	*/
	VectorXT xPeriod = VectorXT::LinSpaced(gridCounts[0], minp(0), maxp(0));
	VectorXT yPeriod = VectorXT::LinSpaced(gridCounts[1], minp(1), maxp(1));
	VectorXT zPeriod = VectorXT::LinSpaced(gridCounts[2], minp(2), maxp(2));

	MatrixXT tmpVec0, tmpVec1, tmpVec2, tmpVec11;
	kron(tmpVec0, VectorXT::Ones(gridCounts[1] * gridCounts[2]), xPeriod);
	kron(tmpVec11, yPeriod, VectorXT::Ones(gridCounts[0]));
	kron(tmpVec1, VectorXT::Ones(gridCounts[2]), tmpVec11);
	kron(tmpVec2, zPeriod, VectorXT::Ones(gridCounts[0] * gridCounts[1]));

	gridCenters.resize(gridCounts[0] * gridCounts[1] * gridCounts[2], 3);
	gridCenters.col(0) = tmpVec0;
	gridCenters.col(1) = tmpVec1;
	gridCenters.col(2) = tmpVec2;


	return true;
}


#ifdef USE_TRIANGLE_H

// 生成圆形面网格，用于展示三维空间中的一个平面：
template <typename DerivedVO, typename DerivedVC, typename DerivedVN>
bool genRoundSurfMesh(Eigen::PlainObjectBase<DerivedVO>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::MatrixBase<DerivedVC>& planeCenter, \
	const Eigen::MatrixBase<DerivedVN>& planeNorm, \
	const double radius = 10, const int versCount = 30);


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
	const Eigen::MatrixBase<DerivedC>& circVers)
{
	using ScalarO = typename DerivedO::Scalar;
	using RowVector3O = Eigen::Matrix<ScalarO, 1, 3>;
	using Matrix3O = Eigen::Matrix<ScalarO, 3, 3>;
	using MatrixXO = Eigen::Matrix<ScalarO, Eigen::Dynamic, Eigen::Dynamic>;
	using MatrixXR = Eigen::Matrix<TRI_REAL, Eigen::Dynamic, Eigen::Dynamic>;

	unsigned circCount = circVers.rows();
	unsigned versCount = circVers.rows();

	// 0. 边缘环路顶点坐标数据拷入到输出网格中。
	vers = circVers.array().cast<ScalarO>();

	// 输入triangulate()的点集是投影到XOY平面的二维点，原点集应该旋转到合适的角度再投影。

	// 1. 取第一个点、1/3处的点、2/3处的点的所在平面的法向量作为原点集的法向量
	RowVector3O vers1 = circVers.row(0).array().cast<ScalarO>();
	RowVector3O vers2 = circVers.row(versCount / 3).array().cast<ScalarO>();
	RowVector3O vers3 = circVers.row(2 * versCount / 3).array().cast<ScalarO>();
	RowVector3O norm = (vers1 - vers2).cross(vers3 - vers2);
	norm.normalize();

	//  2. 旋转点集使得norm平行于z轴
	Matrix3O rotation;
	getRotationMat(rotation, norm, RowVector3O{ 0, 0, 1 });
	vers = (vers * rotation.transpose()).eval();

	// 3. 旋转后的点数据写入到triangulate()接口的输入结构体中。
	MatrixXR vers2D;
	Eigen::MatrixXi edges2D;
	vers2D = vers.transpose().topRows(2).cast<TRI_REAL>();
	edges2D.resize(2, versCount);
	for (unsigned i = 0; i < versCount; i++)
	{
		edges2D(0, i) = i + 1;
		edges2D(1, i) = (i + 1) % versCount + 1;
	}

	TRIANGLE_LIB::triangulateio triIn;
	TRIANGLE_LIB::triangulateio triOut;
	triIn.numberofpoints = versCount;
	triIn.pointlist = vers2D.data();
	triIn.numberofpointattributes = 0;
	triIn.pointattributelist = NULL;
	triIn.pointmarkerlist = NULL;

	triIn.numberofsegments = versCount;
	triIn.segmentlist = (int*)edges2D.data();
	triIn.segmentmarkerlist = NULL;

	triIn.numberoftriangles = 0;
	triIn.numberofcorners = 0;
	triIn.numberoftriangleattributes = 0;

	triIn.numberofholes = 0;
	triIn.holelist = NULL;

	triIn.numberofregions = 0;
	triIn.regionlist = NULL;
	memset(&triOut, 0, sizeof(TRIANGLE_LIB::triangulateio));

	// 4. 执行二维三角剖分，得到输出网格三角片数据，不取顶点坐标数据，顶点坐标数据使用旋转操作前的。
	char triStr[256] = "pYQ";
	TRIANGLE_LIB::triangulate(triStr, &triIn, &triOut, NULL);
	tris.resize(3, triOut.numberoftriangles);
	std::memcpy(tris.data(), triOut.trianglelist, sizeof(int) * 3 * triOut.numberoftriangles);
	tris.transposeInPlace();
	tris.array() -= 1;

	return true;
}


// 重载2：会在回路内部插点，三角片面积上限由参数maxTriArea决定；
template <typename DerivedV, typename DerivedL, typename DerivedN>
bool circuit2mesh(Eigen::PlainObjectBase<DerivedV>& vers, Eigen::MatrixXi& tris, \
	const Eigen::MatrixBase<DerivedL>& versLoop, \
	const Eigen::MatrixBase<DerivedN>& normDir, \
	const float maxTriArea, const bool blRemesh = true)
{
	using ScalarV = typename DerivedV::Scalar;
	using ScalarN = typename DerivedN::Scalar;
	const int versLoopCount = versLoop.rows();

	// 1. 求一个仿射变换：
	Eigen::RowVector3d loopCenter = versLoop.colwise().mean().array().cast<double>();
	Eigen::Matrix4d rotationHomo{ Eigen::Matrix4d::Identity() };
	Eigen::Matrix4d translationHomo{ Eigen::Matrix4d::Identity() };
	Eigen::Matrix4d affineHomo, affineInvHomo;
	Eigen::Matrix3d rotation;
	getRotationMat(rotation, Eigen::RowVector3d{ 0, 0, 1 }, Eigen::RowVector3d{ normDir.array().cast<double>() });
	rotationHomo.topLeftCorner(3, 3) = rotation;
	translationHomo(0, 3) = loopCenter(0);
	translationHomo(1, 3) = loopCenter(1);
	translationHomo(2, 3) = loopCenter(2);
	affineHomo = translationHomo * rotationHomo;
	affineInvHomo = affineHomo.inverse();

	// 2. 输入回路变换到本征坐标系，做插点的三角剖分，再变换回原坐标系；
	Eigen::MatrixXd versLoopHomo0, versLoopHomo, versLoop0, vers2D, versTmp, versTmpHomo;
	Eigen::MatrixXi edges2D(2, versLoopCount);
	Eigen::VectorXi loopVec{ Eigen::VectorXi::LinSpaced(versLoopCount, 0, versLoopCount - 1) };
	char triStr[256];
	versLoop0 = homoVers2VersD(affineInvHomo * vers2HomoVersD(versLoop));
	sprintf_s(triStr, "pq30.0a%fYYQ", maxTriArea);		// "p" - 折线段; "a" - 最大三角片面积; "q" - Quality; "Y" - 边界上不插入点; "Q" - quiet;
	triangulateVers2Mesh(versTmp, tris, versLoop0, std::vector<Eigen::VectorXi>{loopVec}, triStr); 
	versTmp = homoVers2VersD(affineHomo * vers2HomoVersD(versTmp));
	versTmp.topRows(versLoopCount) = versLoop.array().cast<double>();			// 边缘点使用原始坐标；

	// 4. 此时网格内部是个平面，边界和输入边界线相同，需要在保持边界线的情形下光顺网格： 
	if (blRemesh)
	{
		std::vector<int> fixedVerIdxes(versLoopCount);
		for (int i = 0; i < versLoopCount; ++i)
			fixedVerIdxes[i] = i;
		if (!laplaceFaring(vers, versTmp, tris, 0.1, 50, fixedVerIdxes))
			return false;
	}
	else
		vers = versTmp.array().cast<ScalarV>();

	return true;
}


// genCylinder()重载1——生成（类）柱体： 
template <typename DerivedVo, typename DerivedVa, typename DerivedVb>
bool genCylinder(Eigen::PlainObjectBase<DerivedVo>& vers, Eigen::MatrixXi& tris, \
	const Eigen::MatrixBase<DerivedVa>& axisVers, const Eigen::MatrixBase<DerivedVb>& btmVers, \
	const bool isCovered = true)
{
	/*
	bool genCylinder(
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers,
		Eigen::MatrixXi& tris,
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers,			轴线；
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers,			横截面回路顶点，必须要在XOY平面内；
		const bool isCovered																						是否封底
		)


	*/

	// lambda——柱体侧面的三角片生长，会循环调用，调用一次生长一层；
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// 会重复调用，tris容器不需要为空。
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// 待生成表面的圆柱体底圈的第一个顶点的索引。
		sideTris.conservativeResize(currentTrisCount + 2 * circVersCount, 3);

		int triIdx = currentTrisCount;
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx + circVersCount - 1, startIdx, startIdx + 2 * circVersCount - 1 };
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx, startIdx + circVersCount, startIdx + 2 * circVersCount - 1 };
		for (int i = startIdx + 1; i < startIdx + circVersCount; ++i)
		{
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i - 1, i, i + circVersCount - 1 };
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i, i + circVersCount, i + circVersCount - 1 };
		}

		return true;
	};

	using T = typename DerivedVo::Scalar;
	using ScalarB = typename DerivedVb::Scalar;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	using MatrixXB = Eigen::Matrix<ScalarB, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3B = Eigen::Matrix<ScalarB, 3, 3>;
	using RowVector3B = Eigen::Matrix<ScalarB, 1, 3>;

	vers.resize(0, 0);
	tris.resize(0, 0);

	unsigned circVersCount = btmVers.rows();					// 横截面一圈的顶点数；
	unsigned circCount = axisVers.rows();							// 圈数；
	unsigned versCount = circVersCount * circCount;
	std::vector<MatrixXT> circuitsVec(circCount);									// 每个横截面的一圈顶点；
	std::vector<RowVector3T> sectionNorms(circCount);					// 每个横截面的法向；

	// 1. 计算柱体circCount个横截面的法向、每个横截面的顶点； 
	for (unsigned i = 0; i < circCount - 1; ++i)
	{
		sectionNorms[i] = (axisVers.row(i + 1) - axisVers.row(i)).array().cast<T>();
		sectionNorms[i].normalize();
		Matrix3B rotation; 
		getRotationMat(rotation, RowVector3B{ 0, 0, 1 }, \
			RowVector3B{ sectionNorms[i].array().cast<ScalarB>() });
		circuitsVec[i] = (btmVers * rotation.transpose()).array().cast<T>(); 
		 circuitsVec[i].rowwise() += RowVector3T{ axisVers.row(i).array().cast<T>() };		// mat type+=后面只能跟mat type不能跟array type
	}

	// 2. 计算最后一圈顶点：
	RowVector3T deltaNormAve{ RowVector3T::Zero() };
	for (unsigned i = 0; i < circCount - 2; ++i)
		deltaNormAve += (sectionNorms[i + 1] - sectionNorms[i]);
	deltaNormAve.array() /= (circCount - 2);
	sectionNorms[circCount - 1] = sectionNorms[circCount - 2] + deltaNormAve;
	Matrix3T rotation;
	getRotationMat(rotation, RowVector3T{ 0, 0, 1 }, sectionNorms[circCount - 1]);
	circuitsVec[circCount - 1] = btmVers * rotation.transpose();
	circuitsVec[circCount - 1].rowwise() += axisVers.row(circCount - 1);

	// 3. 生成柱体顶点：
	vers.resize(versCount, 3);
	for (unsigned i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.生成侧面三角片：
	for (unsigned i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);

	// 5. 加盖：
	if (isCovered)
	{
		MatrixXT capVers;
		Eigen::MatrixXi capTris1, capTris2;
		circuit2mesh(capVers, capTris1, btmVers);
		capTris2 = capTris1;
		for (int i = 0; i < capTris2.rows(); ++i)
		{
			int tmp = capTris2(i, 2);
			capTris2(i, 2) = capTris2(i, 1);
			capTris2(i, 1) = tmp;
		}
		capTris2.array() += versCount - circVersCount;
		matInsertRows(tris, capTris1);
		matInsertRows(tris, capTris2);
	}


	return true;
}


// genCylinder()重载2——生成圆柱体网格
template <typename DerivedVo, typename DerivedVi>
bool genCylinder(Eigen::PlainObjectBase<DerivedVo>& vers, Eigen::MatrixXi& tris, \
	const Eigen::MatrixBase<DerivedVi>& axisVers, const float radius, \
	const double deltaTheta = 2 * pi / 30, const bool isCovered = true)
{
	// 生成XOY平面上采样角度步长为的圆圈顶点：
	using T = typename DerivedVo::Scalar;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	vers.resize(0, 0);
	tris.resize(0, 0);

	MatrixXT circuit(30, 3);
	circuit.setZero();
	for (unsigned i = 0; i < 30; ++i)
	{
		double theta = deltaTheta * i;
		circuit(i, 0) = static_cast<T>(radius * cos(theta));
		circuit(i, 1) = static_cast<T>(radius * sin(theta));
	}

	genCylinder(vers, tris, axisVers, circuit);

	return true;
}


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


///////////////////////////////////////////////////////////////////////////////////////////////////// SDF相关：
		
// marching cubes中单个立方体内的操作：
template <typename DerivedGV, typename Scalar, typename Index, typename ScalarV, typename IndexF>
void handleCube(const DerivedGV& gridCenters, const Eigen::Matrix<Scalar, 8, 1>& cornerSDF, \
		const Eigen::Matrix<Index, 8, 1>& cornerIdx, const Scalar& isovalue, \
		Eigen::Matrix<ScalarV, Eigen::Dynamic, Eigen::Dynamic>& versResult, Index& curVersCount, \
		Eigen::Matrix<IndexF, Eigen::Dynamic, Eigen::Dynamic>& trisResult, Index& curTrisCount, \
		std::unordered_map<int64_t, int>& edgeIsctMap)
{
	/*
		const DerivedGV& gridCenters,															栅格数据
		const Eigen::Matrix<Scalar, 8, 1>& cornerSDF,									当前立方体八个顶点的SDF值
		const Eigen::Matrix<Index, 8, 1>& cornerIdx,									当前立方体八个顶点在栅格中的索引；
		const Scalar& isovalue,																		需要提取的等值面的SDF值
		Eigen::PlainObjectBase<DerivedV>& versResult,								输出网格的顶点
		Index& curVersCount,																			当前累计生成的输出网格顶点数
		Eigen::PlainObjectBase<DerivedF>& trisResult,									输出网格的三角片
		Index& curTrisCount,																			当前累计生成的输出网格三角片数
		std::unordered_map<int64_t, int>& edgeIsctMap								边编码-边交点索引的哈希表；

	*/

	Eigen::Matrix<Index, 12, 1> isctVerIdxes;		// 立方体边上的交点的绝对索引——是在最终输出网格中的索引；
	int cornerState = 0;											// 立方体顶点状态编码；256种情形；

// 生成无向边编码
	const auto genMCedgeCode = [](int32_t vaIdx, int32_t vbIdx)
	{
		if (vaIdx > vbIdx)
			std::swap(vaIdx, vbIdx);
		std::int64_t edgeCode = 0;
		edgeCode |= vaIdx;
		edgeCode |= static_cast<std::int64_t>(vbIdx) << 32;
		return edgeCode;
	};

	// 1. 计算当前立方体的顶点状态编码，即8个顶点在等值面的内外状态1；
	for (int i = 0; i < 8; i++)
		if (cornerSDF(i) > isovalue)
			cornerState |= 1 << i;

	// 2. 确定当前立方体中和等值面相交的边；
	int edgeState = MARCHING_CUBES_TABLES::edgeStateCodes[cornerState];		// 立方体顶点状态编码映射为相交边编码；
	if (edgeState == 0)
		return;															// 表示当前立方体整体都在等值面外部或内部，没有交点；

	// 3. 确定等值面和当前立方体的边的交点； Find the point of intersection of the surface with each edge. Then find the normal to the surface at those points
	for (int i = 0; i < 12; i++)						// 对立方体所有边的遍历
	{
		if (edgeState & (1 << i))					// 若等值面和当前边相交：
		{
			int vaIdxRela = MARCHING_CUBES_TABLES::cubeEdges[i][0];			// 当前边两端点的相对索引；
			int vbIdxRela = MARCHING_CUBES_TABLES::cubeEdges[i][1];

			// 生成边上的顶点：
			int vaIdx = cornerIdx(vaIdxRela);				// 当前边两端点的绝对索引，是栅格中的顶点索引；
			int vbIdx = cornerIdx(vbIdxRela);
			std::int64_t edgeCode = genMCedgeCode(vaIdx, vbIdx);
			const auto iter = edgeIsctMap.find(edgeCode);

			if (iter == edgeIsctMap.end())								// 若当前边交点未生成；
			{
				if (curVersCount == versResult.rows())
					versResult.conservativeResize(versResult.rows() * 2 + 1, versResult.cols());

				// 插值生成新的顶点：find crossing point assuming linear interpolation along edges
				const Scalar& SDFa = cornerSDF(vaIdxRela);			// 当前边两端点的SDF值；
				const Scalar& SDFb = cornerSDF(vbIdxRela);
				const Scalar delta = SDFb - SDFa;
				Scalar t = (isovalue - SDFa) / delta;
				versResult.row(curVersCount) = (gridCenters.row(vaIdx) + t * (gridCenters.row(vbIdx) - gridCenters.row(vaIdx))).array().cast<ScalarV>();

				isctVerIdxes[i] = curVersCount;
				edgeIsctMap[edgeCode] = isctVerIdxes[i];
				curVersCount++;
			}
			else                                                                             // 若当前边交点已生成；
				isctVerIdxes[i] = iter->second;

			assert(isctVerIdxes[i] >= 0);
			assert(isctVerIdxes[i] < curVersCount);
		}
	}

	// 4. 生成当前立方体中的三角片，一个立方体中最多生成5个三角片；
	for (int i = 0; i < 5; i++)
	{
		if (MARCHING_CUBES_TABLES::cubeTriangles[cornerState][3 * i] < 0)
			break;

		if (curTrisCount == trisResult.rows())
			trisResult.conservativeResize(trisResult.rows() * 2 + 1, trisResult.cols());

		// 新增三角片数据中的顶点索引，是相对索引；
		int vaIdxRela = MARCHING_CUBES_TABLES::cubeTriangles[cornerState][3 * i + 0];
		int vbIdxRela = MARCHING_CUBES_TABLES::cubeTriangles[cornerState][3 * i + 1];
		int vcIdxRela = MARCHING_CUBES_TABLES::cubeTriangles[cornerState][3 * i + 2];

		assert(isctVerIdxes[vaIdxRela] >= 0);
		assert(isctVerIdxes[vbIdxRela] >= 0);
		assert(isctVerIdxes[vcIdxRela] >= 0);

		// 相对索引转换为绝对索引，插入新增三角片
		trisResult.row(curTrisCount) << isctVerIdxes[vaIdxRela], isctVerIdxes[vbIdxRela], isctVerIdxes[vcIdxRela];
		curTrisCount++;
	}

}


template <typename DerivedV, typename DerivedS, typename DerivedGV>
bool marchingCubes(Eigen::PlainObjectBase<DerivedV>& versResult, \
	Eigen::MatrixXi& trisResult, const Eigen::MatrixBase<DerivedS>& SDFvec, \
	const Eigen::MatrixBase<DerivedGV>& gridCenters, \
	const unsigned nx, const unsigned ny, const unsigned nz,\
	const typename DerivedS::Scalar isovalue, const bool blSClargest)
{
	/*
		const Eigen::MatrixBase<DerivedS>& SDFvec,							符号距离场数据
		const Eigen::MatrixBase<DerivedGV>& gridCenters,					栅格数据
		const unsigned nx,																			x方向上栅格个数
		const unsigned ny,
		const unsigned nz,
		const typename DerivedS::Scalar isovalue,										需要提取的水平集的SDF值；
		Eigen::PlainObjectBase<DerivedV>& versResult,							输出网格顶点
		Eigen::PlainObjectBase<DerivedF>& trisResult								输出网格三角片
	*/
	using ScalarV = typename DerivedV::Scalar;
	using ScalarS = typename DerivedS::Scalar; 

	// lambda——栅格的三维索引映射到一维索引：
	const auto getGridIdx = [&nx, &ny, &nz](const int& x, const int& y, const int& z)->unsigned
	{
		return x + nx * (y + ny * (z));
	};

	const unsigned cornerIdxOffset[8] = { 0, 1, 1 + nx, nx, nx * ny, 1 + nx * ny, 1 + nx + nx * ny, nx + nx * ny };	// 立方体八个顶点的索引偏移量；
	std::unordered_map<int64_t, int> edgeIsctMap;							// 边编码-边交点索引的哈希表；

	unsigned curVersCount = 0;
	unsigned curTrisCount = 0;

	// 1. march over all cubes (loop order chosen to match memory)
	/*
		 Should be possible to parallelize safely if threads are "well separated".
		 Like red-black Gauss Seidel. Probably each thread need's their own edgeIsctMap, versResult, trisResult,
			   and then merge at the end.
		 Annoying part are the edges lying on the  interface between chunks.
	*/
	Eigen::Matrix<ScalarV, Eigen::Dynamic, Eigen::Dynamic> versTmp;
	Eigen::MatrixXi trisTmp;
	versTmp.resize(std::pow(nx * ny * nz, 2. / 3.), 3);
	trisTmp.resize(std::pow(nx * ny * nz, 2. / 3.), 3);
	for (int z = 0; z < nz - 1; z++)
	{
		for (int y = 0; y < ny - 1; y++)
		{
			for (int x = 0; x < nx - 1; x++)
			{
				// 1.1 计算当前栅格的索引：
				const unsigned gridIdx = getGridIdx(x, y, z);

				// 1.2 计算当前栅格对应的立方体的八个顶点的数据；
				static Eigen::Matrix<ScalarS, 8, 1> cornerSDF;						// 方块的八个顶点的SDF值
				static Eigen::Matrix<unsigned, 8, 1> cornerIdx;               // 方块的八个顶点在栅格中的索引
				for (int i = 0; i < 8; i++)
				{
					const unsigned originIdx = gridIdx + cornerIdxOffset[i];
					cornerIdx(i) = originIdx;
					cornerSDF(i) = SDFvec(originIdx);
				}

				// 1.3 生成当前立方体内的三角片
				handleCube(gridCenters, cornerSDF, cornerIdx, isovalue, versTmp, curVersCount, trisTmp, curTrisCount, edgeIsctMap);
			}
		}
	}

	// 2. shrink_to_fit();
	versTmp.conservativeResize(curVersCount, 3);
	trisTmp.conservativeResize(curTrisCount, 3);

	// 3. 提取最大联通区域（颌网格距离场生成的MC结果可能包含微小的孤立网格）
	if (blSClargest)
	{
		if (!simplyConnectedLargest(versResult, trisResult, versTmp, trisTmp))
			return false;
		else
			return true;
	}
	else
	{
		versResult = versTmp;
		trisResult = trisTmp;
	}

	return true;
}


bool SDFvec2mat(std::vector<Eigen::MatrixXf>& matLayers, \
		const std::vector<float>& SDFvec, const std::vector<unsigned>& stepsCount);

bool SDFmat2vec(std::vector<float>& SDFvec, \
	const std::vector<Eigen::MatrixXf>& matLayers, const std::vector<unsigned>& stepsCount);


// 暂时未整理的实现：
#include "temp.tpp"