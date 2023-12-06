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
  

///////////////////////////////////////////////////////////////////////////////////////////////////// modeling�ӿڣ�

// ����������ԭ�㣬�߳�Ϊ1������Ƭ��Ϊ12������������
template	<typename T>
void genCubeMesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);


// ���������Χ�е���������
template <typename T>
void genAABBmesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<T, 3>& aabb);


// ����դ���������
template<typename DerivedV, typename DerivedVO>
bool genGrids(Eigen::PlainObjectBase<DerivedV>& gridCenters, \
	const Eigen::MatrixBase<DerivedVO>& origin, \
	const float step, const std::vector<unsigned>& gridCounts)
{
	/*
		bool genGrids(
			gridCenters,											 �����դ�����
			origin,													 դ��ԭ�㣬������������С���Ǹ�դ��㣻
			step,														 ��������
			gridCounts											 ��Ԫ���飬�ֱ���XYZ�����ϵĲ�����
			)
	*/
	using ScalarT = typename DerivedV::Scalar;
	using VectorXT = Eigen::Matrix<ScalarT, Eigen::Dynamic, 1>;
	using RowVector3T = Eigen::Matrix<ScalarT, 1, 3>;
	using MatrixXT = Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>;
	assert((3 == origin.size()) && "Assert!!! origin should be a 3D vertex.");
	gridCenters.resize(0, 0);

	// �������볡�İ�Χ�У�
	RowVector3T minp{ static_cast<ScalarT>(origin(0)), static_cast<ScalarT>(origin(1)), static_cast<ScalarT>(origin(2)) };
	RowVector3T maxp = minp + static_cast<ScalarT>(step) * RowVector3T(gridCounts[0] - 1, gridCounts[1] - 1, gridCounts[2] - 1);

	// ����դ��
	/*
		�������������е�դ�����ĵ�Ϊ��
		gc(000), gc(100), gc(200), gc(300),...... gc(010), gc(110), gc(210), gc(310),...... gc(001), gc(101), gc(201).....

		x���꣺
		x0, x1, x2, x3......x0, x1, x2, x3......x0, x1, x2, x3......
		����ΪxCount;
		�ظ�����Ϊ(yCount * zCount)

		y���꣺
		y0, y0, y0...y1, y1, y1...y2, y2, y2.........y0, y0, y0...y1, y1, y1...
		����Ϊ(xCount * yCount);
		�ظ�����ΪzCount;
		����Ԫ���ظ�����ΪxCount

		z���꣺
		z0, z0, z0......z1, z1, z1......z2, z2, z2......
		����Ԫ���ظ�����Ϊ(xCount * yCount)
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

// ����Բ������������չʾ��ά�ռ��е�һ��ƽ�棺
template <typename DerivedVO, typename DerivedVC, typename DerivedVN>
bool genRoundSurfMesh(Eigen::PlainObjectBase<DerivedVO>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::MatrixBase<DerivedVC>& planeCenter, \
	const Eigen::MatrixBase<DerivedVN>& planeNorm, \
	const double radius = 10, const int versCount = 30);


// ����1��2D���������ʷֵõ������񡪡����Դ���Ҳ���Բ�����
/*

	ע��
		������Ʊ�����XOYƽ���ڣ�������2D��Ҳ������3D�ģ�
		Ĭ�������ʷ�ģʽΪ"pY"����ʾ�����ʷֳɲ�����ƽ��ֱ��ͼ��


		switch string:
				-p			�����ʷ�����һ��ƽ��ֱ��ͼ
				-r			(refine a previously generated mesh)��һ������������н�һ���������ʷ֣�
				-q			(Quality mesh generation)�����һ����ֵ����"-q30"��ʾ�����ʷֽ���в����Դ���С��30��Ľǣ�
				-a			�����һ����ֵ����"-a5"ָ�������ʷֽ��������Ƭ���������5mm^2;
				-Y			(Prohibits the insertion of Steiner points on the mesh boundary.)
							��ֹ�ڱ�Ե���ϲ����µ㣻
				-YY		prohibits the insertion of Steiner points on any segment, including internal segments.
							��ֹ���κ�ԭ�б��ϲ����µ㣻

*/
template <typename DerivedVo, typename DerivedI, typename DerivedVi, typename DerivedVh>
bool triangulateVers2Mesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const std::vector<Eigen::VectorXi>& bdryLoops, \
	const Eigen::PlainObjectBase<DerivedVh>& holeCenters, \
	const char* strSwitcher = "pY");


// ����2��2D���������ʷֵõ������񡪡�������
template <typename DerivedVo, typename DerivedI, typename DerivedVi>
bool triangulateVers2Mesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const std::vector<Eigen::VectorXi>& bdryLoops, \
	const char* strSwitcher = "pY");

// �����ʷ���������������
template <typename DerivedVo, typename DerivedIt, typename DerivedVi, typename DerivedIe>
bool triangulateRefineMesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedIt>& trisOut, const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const Eigen::PlainObjectBase<DerivedIt>& trisIn, const Eigen::PlainObjectBase<DerivedIe>& edges, \
	const char* strSwitcher);


//	����1����ձ߽��ߵ㼯�õ������񣬿�����ƽ��Ҳ���������棬���������ڲ���㣬����Ƭ�ߴ粻�ɿء�
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

	// 0. ��Ե��·�����������ݿ��뵽��������С�
	vers = circVers.array().cast<ScalarO>();

	// ����triangulate()�ĵ㼯��ͶӰ��XOYƽ��Ķ�ά�㣬ԭ�㼯Ӧ����ת�����ʵĽǶ���ͶӰ��

	// 1. ȡ��һ���㡢1/3���ĵ㡢2/3���ĵ������ƽ��ķ�������Ϊԭ�㼯�ķ�����
	RowVector3O vers1 = circVers.row(0).array().cast<ScalarO>();
	RowVector3O vers2 = circVers.row(versCount / 3).array().cast<ScalarO>();
	RowVector3O vers3 = circVers.row(2 * versCount / 3).array().cast<ScalarO>();
	RowVector3O norm = (vers1 - vers2).cross(vers3 - vers2);
	norm.normalize();

	//  2. ��ת�㼯ʹ��normƽ����z��
	Matrix3O rotation;
	getRotationMat(rotation, norm, RowVector3O{ 0, 0, 1 });
	vers = (vers * rotation.transpose()).eval();

	// 3. ��ת��ĵ�����д�뵽triangulate()�ӿڵ�����ṹ���С�
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

	// 4. ִ�ж�ά�����ʷ֣��õ������������Ƭ���ݣ���ȡ�����������ݣ�������������ʹ����ת����ǰ�ġ�
	char triStr[256] = "pYQ";
	TRIANGLE_LIB::triangulate(triStr, &triIn, &triOut, NULL);
	tris.resize(3, triOut.numberoftriangles);
	std::memcpy(tris.data(), triOut.trianglelist, sizeof(int) * 3 * triOut.numberoftriangles);
	tris.transposeInPlace();
	tris.array() -= 1;

	return true;
}


// ����2�����ڻ�·�ڲ���㣬����Ƭ��������ɲ���maxTriArea������
template <typename DerivedV, typename DerivedL, typename DerivedN>
bool circuit2mesh(Eigen::PlainObjectBase<DerivedV>& vers, Eigen::MatrixXi& tris, \
	const Eigen::MatrixBase<DerivedL>& versLoop, \
	const Eigen::MatrixBase<DerivedN>& normDir, \
	const float maxTriArea, const bool blRemesh = true)
{
	using ScalarV = typename DerivedV::Scalar;
	using ScalarN = typename DerivedN::Scalar;
	const int versLoopCount = versLoop.rows();

	// 1. ��һ������任��
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

	// 2. �����·�任����������ϵ�������������ʷ֣��ٱ任��ԭ����ϵ��
	Eigen::MatrixXd versLoopHomo0, versLoopHomo, versLoop0, vers2D, versTmp, versTmpHomo;
	Eigen::MatrixXi edges2D(2, versLoopCount);
	Eigen::VectorXi loopVec{ Eigen::VectorXi::LinSpaced(versLoopCount, 0, versLoopCount - 1) };
	char triStr[256];
	versLoop0 = homoVers2VersD(affineInvHomo * vers2HomoVersD(versLoop));
	sprintf_s(triStr, "pq30.0a%fYYQ", maxTriArea);		// "p" - ���߶�; "a" - �������Ƭ���; "q" - Quality; "Y" - �߽��ϲ������; "Q" - quiet;
	triangulateVers2Mesh(versTmp, tris, versLoop0, std::vector<Eigen::VectorXi>{loopVec}, triStr); 
	versTmp = homoVers2VersD(affineHomo * vers2HomoVersD(versTmp));
	versTmp.topRows(versLoopCount) = versLoop.array().cast<double>();			// ��Ե��ʹ��ԭʼ���ꣻ

	// 4. ��ʱ�����ڲ��Ǹ�ƽ�棬�߽������߽�����ͬ����Ҫ�ڱ��ֱ߽��ߵ������¹�˳���� 
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


// genCylinder()����1�������ɣ��ࣩ���壺 
template <typename DerivedVo, typename DerivedVa, typename DerivedVb>
bool genCylinder(Eigen::PlainObjectBase<DerivedVo>& vers, Eigen::MatrixXi& tris, \
	const Eigen::MatrixBase<DerivedVa>& axisVers, const Eigen::MatrixBase<DerivedVb>& btmVers, \
	const bool isCovered = true)
{
	/*
	bool genCylinder(
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers,
		Eigen::MatrixXi& tris,
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers,			���ߣ�
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers,			������·���㣬����Ҫ��XOYƽ���ڣ�
		const bool isCovered																						�Ƿ���
		)


	*/

	// lambda����������������Ƭ��������ѭ�����ã�����һ������һ�㣻
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// ���ظ����ã�tris��������ҪΪ�ա�
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// �����ɱ����Բ�����Ȧ�ĵ�һ�������������
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

	unsigned circVersCount = btmVers.rows();					// �����һȦ�Ķ�������
	unsigned circCount = axisVers.rows();							// Ȧ����
	unsigned versCount = circVersCount * circCount;
	std::vector<MatrixXT> circuitsVec(circCount);									// ÿ��������һȦ���㣻
	std::vector<RowVector3T> sectionNorms(circCount);					// ÿ�������ķ���

	// 1. ��������circCount�������ķ���ÿ�������Ķ��㣻 
	for (unsigned i = 0; i < circCount - 1; ++i)
	{
		sectionNorms[i] = (axisVers.row(i + 1) - axisVers.row(i)).array().cast<T>();
		sectionNorms[i].normalize();
		Matrix3B rotation; 
		getRotationMat(rotation, RowVector3B{ 0, 0, 1 }, \
			RowVector3B{ sectionNorms[i].array().cast<ScalarB>() });
		circuitsVec[i] = (btmVers * rotation.transpose()).array().cast<T>(); 
		 circuitsVec[i].rowwise() += RowVector3T{ axisVers.row(i).array().cast<T>() };		// mat type+=����ֻ�ܸ�mat type���ܸ�array type
	}

	// 2. �������һȦ���㣺
	RowVector3T deltaNormAve{ RowVector3T::Zero() };
	for (unsigned i = 0; i < circCount - 2; ++i)
		deltaNormAve += (sectionNorms[i + 1] - sectionNorms[i]);
	deltaNormAve.array() /= (circCount - 2);
	sectionNorms[circCount - 1] = sectionNorms[circCount - 2] + deltaNormAve;
	Matrix3T rotation;
	getRotationMat(rotation, RowVector3T{ 0, 0, 1 }, sectionNorms[circCount - 1]);
	circuitsVec[circCount - 1] = btmVers * rotation.transpose();
	circuitsVec[circCount - 1].rowwise() += axisVers.row(circCount - 1);

	// 3. �������嶥�㣺
	vers.resize(versCount, 3);
	for (unsigned i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.���ɲ�������Ƭ��
	for (unsigned i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);

	// 5. �Ӹǣ�
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


// genCylinder()����2��������Բ��������
template <typename DerivedVo, typename DerivedVi>
bool genCylinder(Eigen::PlainObjectBase<DerivedVo>& vers, Eigen::MatrixXi& tris, \
	const Eigen::MatrixBase<DerivedVi>& axisVers, const float radius, \
	const double deltaTheta = 2 * pi / 30, const bool isCovered = true)
{
	// ����XOYƽ���ϲ����ǶȲ���Ϊ��ԲȦ���㣺
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


// genCylinder()����3�������ɷ���������
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers,\
	const std::pair<float, float> sizePair, const float SR, const bool isCovered);


// genCylinder()����4���������ϵ�����µ����2D�߽绷· ���������壺
/*
	bool genCylinder(
			MatrixXT& vers,								������񶥵�
			Eigen::MatrixXi& tris,						�����������Ƭ
			const Eigen::PlainObjectBase<DerivedVa>& axisVers,			��������
			const Eigen::PlainObjectBase<DerivedVt>& topLoop,			�ϵ���߽绷·
			const Eigen::PlainObjectBase<DerivedVb>& btmLoop,			�µ���߽绷·
			const bool isCovered																�Ƿ���
			)
	ע�����µ���߽绷·������XOYƽ���еĵ㼯��
					�����߶�������Ҫ��ͬ��
					�Ҵ�Z�����򿴻�·��������Ӧ������˳ʱ�뷽������

*/
template <typename T, typename DerivedVa, typename DerivedVt, typename DerivedVb>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedVa>& axisVers, const Eigen::PlainObjectBase<DerivedVt>& topLoop, \
	const Eigen::PlainObjectBase<DerivedVb>& btmLoop, const bool isCovered);


// genCylinder()����5���������ϵ�����µ����3D�߽绷· ���������壺
/*
	bool genCylinder(
			Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers,		����������
			Eigen::MatrixXi& tris,																	�����������Ƭ
			const Eigen::PlainObjectBase<DerivedVt>& topLoop,				�ϵ���߽��·
			const Eigen::PlainObjectBase<DerivedVb>& btmLoop,				�±���߽��·
			const Eigen::Matrix<ScalarN, 1, 3>& topNorm							�ϵ��淨��
			const Eigen::Matrix<ScalarN, 1, 3>& btmNorm							�µ��淨��
			const int layersCount,																	����
			const float maxArea																	���µ��������ʷֵ��������Ƭ���
			const bool isCovered																	�Ƿ���
			)
	ע��Ҫ�����µ���߽绷·��������Ҫ��ͬ��
				Ҫ���ϵ���߽��·���������������������ƽ�����ϵ��淨���µ����෴��
				��������Ϊ��ֱ��Բ��ʱ�����µ����·������������Ϊ��ֱ���ϣ�
			layersCount���������屻��·����ֳ��˶��ٲ�
				��СΪ1���������µ���֮�䲻���
				���µ���֮���Ļ�·����Ȧ��Ϊ(layersCount - 1);
*/
template <typename DerivedVo, typename DerivedVt, typename DerivedVb, typename ScalarN>
bool genCylinder(Eigen::PlainObjectBase<DerivedVo>& vers, Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedVt>& topLoop, \
	const Eigen::PlainObjectBase<DerivedVb>& btmLoop, \
	const Eigen::Matrix<ScalarN, 1, 3>& topNorm, \
	const Eigen::Matrix<ScalarN, 1, 3>& btmNorm, \
	const int layersCount, const float maxArea, const bool isCovered);
 

//		���ɷ�������ת�����Σ���ȷ�������XOYƽ��ƽ�л�ֱ��
template <typename T>
bool genAlignedCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const std::pair<float, float> sizePair, const float SR, \
	const bool isCovered);

#endif


///////////////////////////////////////////////////////////////////////////////////////////////////// SDF��أ�
		
// marching cubes�е����������ڵĲ�����
template <typename DerivedGV, typename Scalar, typename Index, typename ScalarV, typename IndexF>
void handleCube(const DerivedGV& gridCenters, const Eigen::Matrix<Scalar, 8, 1>& cornerSDF, \
		const Eigen::Matrix<Index, 8, 1>& cornerIdx, const Scalar& isovalue, \
		Eigen::Matrix<ScalarV, Eigen::Dynamic, Eigen::Dynamic>& versResult, Index& curVersCount, \
		Eigen::Matrix<IndexF, Eigen::Dynamic, Eigen::Dynamic>& trisResult, Index& curTrisCount, \
		std::unordered_map<int64_t, int>& edgeIsctMap)
{
	/*
		const DerivedGV& gridCenters,															դ������
		const Eigen::Matrix<Scalar, 8, 1>& cornerSDF,									��ǰ������˸������SDFֵ
		const Eigen::Matrix<Index, 8, 1>& cornerIdx,									��ǰ������˸�������դ���е�������
		const Scalar& isovalue,																		��Ҫ��ȡ�ĵ�ֵ���SDFֵ
		Eigen::PlainObjectBase<DerivedV>& versResult,								�������Ķ���
		Index& curVersCount,																			��ǰ�ۼ����ɵ�������񶥵���
		Eigen::PlainObjectBase<DerivedF>& trisResult,									������������Ƭ
		Index& curTrisCount,																			��ǰ�ۼ����ɵ������������Ƭ��
		std::unordered_map<int64_t, int>& edgeIsctMap								�߱���-�߽��������Ĺ�ϣ��

	*/

	Eigen::Matrix<Index, 12, 1> isctVerIdxes;		// ��������ϵĽ���ľ���������������������������е�������
	int cornerState = 0;											// �����嶥��״̬���룻256�����Σ�

// ��������߱���
	const auto genMCedgeCode = [](int32_t vaIdx, int32_t vbIdx)
	{
		if (vaIdx > vbIdx)
			std::swap(vaIdx, vbIdx);
		std::int64_t edgeCode = 0;
		edgeCode |= vaIdx;
		edgeCode |= static_cast<std::int64_t>(vbIdx) << 32;
		return edgeCode;
	};

	// 1. ���㵱ǰ������Ķ���״̬���룬��8�������ڵ�ֵ�������״̬1��
	for (int i = 0; i < 8; i++)
		if (cornerSDF(i) > isovalue)
			cornerState |= 1 << i;

	// 2. ȷ����ǰ�������к͵�ֵ���ཻ�ıߣ�
	int edgeState = MARCHING_CUBES_TABLES::edgeStateCodes[cornerState];		// �����嶥��״̬����ӳ��Ϊ�ཻ�߱��룻
	if (edgeState == 0)
		return;															// ��ʾ��ǰ���������嶼�ڵ�ֵ���ⲿ���ڲ���û�н��㣻

	// 3. ȷ����ֵ��͵�ǰ������ıߵĽ��㣻 Find the point of intersection of the surface with each edge. Then find the normal to the surface at those points
	for (int i = 0; i < 12; i++)						// �����������бߵı���
	{
		if (edgeState & (1 << i))					// ����ֵ��͵�ǰ���ཻ��
		{
			int vaIdxRela = MARCHING_CUBES_TABLES::cubeEdges[i][0];			// ��ǰ�����˵�����������
			int vbIdxRela = MARCHING_CUBES_TABLES::cubeEdges[i][1];

			// ���ɱ��ϵĶ��㣺
			int vaIdx = cornerIdx(vaIdxRela);				// ��ǰ�����˵�ľ�����������դ���еĶ���������
			int vbIdx = cornerIdx(vbIdxRela);
			std::int64_t edgeCode = genMCedgeCode(vaIdx, vbIdx);
			const auto iter = edgeIsctMap.find(edgeCode);

			if (iter == edgeIsctMap.end())								// ����ǰ�߽���δ���ɣ�
			{
				if (curVersCount == versResult.rows())
					versResult.conservativeResize(versResult.rows() * 2 + 1, versResult.cols());

				// ��ֵ�����µĶ��㣺find crossing point assuming linear interpolation along edges
				const Scalar& SDFa = cornerSDF(vaIdxRela);			// ��ǰ�����˵��SDFֵ��
				const Scalar& SDFb = cornerSDF(vbIdxRela);
				const Scalar delta = SDFb - SDFa;
				Scalar t = (isovalue - SDFa) / delta;
				versResult.row(curVersCount) = (gridCenters.row(vaIdx) + t * (gridCenters.row(vbIdx) - gridCenters.row(vaIdx))).array().cast<ScalarV>();

				isctVerIdxes[i] = curVersCount;
				edgeIsctMap[edgeCode] = isctVerIdxes[i];
				curVersCount++;
			}
			else                                                                             // ����ǰ�߽��������ɣ�
				isctVerIdxes[i] = iter->second;

			assert(isctVerIdxes[i] >= 0);
			assert(isctVerIdxes[i] < curVersCount);
		}
	}

	// 4. ���ɵ�ǰ�������е�����Ƭ��һ�����������������5������Ƭ��
	for (int i = 0; i < 5; i++)
	{
		if (MARCHING_CUBES_TABLES::cubeTriangles[cornerState][3 * i] < 0)
			break;

		if (curTrisCount == trisResult.rows())
			trisResult.conservativeResize(trisResult.rows() * 2 + 1, trisResult.cols());

		// ��������Ƭ�����еĶ��������������������
		int vaIdxRela = MARCHING_CUBES_TABLES::cubeTriangles[cornerState][3 * i + 0];
		int vbIdxRela = MARCHING_CUBES_TABLES::cubeTriangles[cornerState][3 * i + 1];
		int vcIdxRela = MARCHING_CUBES_TABLES::cubeTriangles[cornerState][3 * i + 2];

		assert(isctVerIdxes[vaIdxRela] >= 0);
		assert(isctVerIdxes[vbIdxRela] >= 0);
		assert(isctVerIdxes[vcIdxRela] >= 0);

		// �������ת��Ϊ����������������������Ƭ
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
		const Eigen::MatrixBase<DerivedS>& SDFvec,							���ž��볡����
		const Eigen::MatrixBase<DerivedGV>& gridCenters,					դ������
		const unsigned nx,																			x������դ�����
		const unsigned ny,
		const unsigned nz,
		const typename DerivedS::Scalar isovalue,										��Ҫ��ȡ��ˮƽ����SDFֵ��
		Eigen::PlainObjectBase<DerivedV>& versResult,							������񶥵�
		Eigen::PlainObjectBase<DerivedF>& trisResult								�����������Ƭ
	*/
	using ScalarV = typename DerivedV::Scalar;
	using ScalarS = typename DerivedS::Scalar; 

	// lambda����դ�����ά����ӳ�䵽һά������
	const auto getGridIdx = [&nx, &ny, &nz](const int& x, const int& y, const int& z)->unsigned
	{
		return x + nx * (y + ny * (z));
	};

	const unsigned cornerIdxOffset[8] = { 0, 1, 1 + nx, nx, nx * ny, 1 + nx * ny, 1 + nx + nx * ny, nx + nx * ny };	// ������˸����������ƫ������
	std::unordered_map<int64_t, int> edgeIsctMap;							// �߱���-�߽��������Ĺ�ϣ��

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
				// 1.1 ���㵱ǰդ���������
				const unsigned gridIdx = getGridIdx(x, y, z);

				// 1.2 ���㵱ǰդ���Ӧ��������İ˸���������ݣ�
				static Eigen::Matrix<ScalarS, 8, 1> cornerSDF;						// ����İ˸������SDFֵ
				static Eigen::Matrix<unsigned, 8, 1> cornerIdx;               // ����İ˸�������դ���е�����
				for (int i = 0; i < 8; i++)
				{
					const unsigned originIdx = gridIdx + cornerIdxOffset[i];
					cornerIdx(i) = originIdx;
					cornerSDF(i) = SDFvec(originIdx);
				}

				// 1.3 ���ɵ�ǰ�������ڵ�����Ƭ
				handleCube(gridCenters, cornerSDF, cornerIdx, isovalue, versTmp, curVersCount, trisTmp, curTrisCount, edgeIsctMap);
			}
		}
	}

	// 2. shrink_to_fit();
	versTmp.conservativeResize(curVersCount, 3);
	trisTmp.conservativeResize(curTrisCount, 3);

	// 3. ��ȡ�����ͨ�����������볡���ɵ�MC������ܰ���΢С�Ĺ�������
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


// ��ʱδ�����ʵ�֣�
#include "temp.tpp"