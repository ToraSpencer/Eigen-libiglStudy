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


//  64λ�����£�triangulate.cpp��" x = vertexloop[0] = pointlist[coordindex++];"��һ����׳��쳣����Ҫ�о�
#ifndef _WIN64
#define USE_TRIANGLE_H
#endif

#ifdef  USE_TRIANGLE_H
#include "triangulate.h"
#endif


///////////////////////////////////////////////////////////////////////////////////////////////////// modeling�ӿڣ�

// �����ֵ�õ�ֱ�߶ε���
template <typename T>
bool interpolateToLine(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::Matrix<T, 1, 3>& start, \
	const Eigen::Matrix<T, 1, 3>& end, const float SR, const bool SE = true);

// ����XOYƽ���ڵ�ԲȦ�㼯��
template <typename T>
bool getCircleVers(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, \
	const unsigned versCount = 20);

// �Ի�·�㼯���в���������ʷ֡����ֶ���Ϸ�����
template <typename IndexType>
bool circuitGetTris(Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic>& tris, \
	const std::vector<int>& indexes, const bool regularTri);

// ����������ԭ�㣬�߳�Ϊ1������Ƭ��Ϊ12������������
template	<typename T>
void genCubeMesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);

// ���������Χ�е���������
template <typename T>
void genAABBmesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<T, 3>& aabb);

// ����դ������㣺
template<typename Tg, typename To>
bool genGrids(Eigen::Matrix<Tg, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, const Eigen::Matrix<To, 1, 3>& origin, \
	const float step, const std::vector<unsigned>& gridCounts);


#ifdef USE_TRIANGLE_H
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers, \
	const bool isCovered);

template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const float radius, const double deltaTheta, const bool isCovered);

template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const std::pair<float, float> sizePair, const float SR, \
	const bool isCovered);

template <typename T>
bool genAlignedCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const std::pair<float, float> sizePair, const float SR, \
	const bool isCovered);


template <typename T>
bool circuit2mesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& circVers);
#endif



// ��ʱδ�����ʵ�֣�
#include "temp.tpp"