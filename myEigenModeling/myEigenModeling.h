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


//  64λ�����£�triangulate.cpp��" x = vertexloop[0] = pointlist[coordindex++];"��һ����׳��쳣����Ҫ�о�
#ifndef _WIN64
#define USE_TRIANGLE_H
#endif

#ifdef  USE_TRIANGLE_H
#include "triangulate.h"
#endif
 

///////////////////////////////////////////////////////////////////////////////////////////////////// modeling�ӿڣ�

// ����������ԭ�㣬�߳�Ϊ1������Ƭ��Ϊ12������������
template	<typename T>
void genCubeMesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);


// ���������Χ�е���������
template <typename T>
void genAABBmesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<T, 3>& aabb);


// ����դ������㣺
template<typename Tg, typename To>
bool genGrids(Eigen::Matrix<Tg, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, \
	const Eigen::Matrix<To, 1, 3>& origin, \
	const float step, const std::vector<unsigned>& gridCounts);


#ifdef USE_TRIANGLE_H

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
	const Eigen::PlainObjectBase<DerivedC>& circVers);


// ����2�����ڻ�·�ڲ���㣬����Ƭ��������ɲ���maxTriArea������
template <typename DerivedV, typename DerivedL, typename DerivedN>
bool circuit2mesh(Eigen::PlainObjectBase<DerivedV>& vers, Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedL>& versLoop, \
	const Eigen::PlainObjectBase<DerivedN>& normDir, const float maxTriArea);


// genCylinder()����1�������ɣ��ࣩ���壺
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers, \
	const bool isCovered);


// genCylinder()����2��������Բ��������
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers,\
	const float radius, const double deltaTheta, const bool isCovered);

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



// ��ʱδ�����ʵ�֣�
#include "temp.tpp"