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


/////////////////////////////////////////////////////////////////////////////////////////////////// ����ת���ӿڣ�

// std::vector<std::pair<int, int>>��ʾ�ı�����ת��Ϊ�����ʾ��
template <typename DerivedI>
void edges2mat(Eigen::PlainObjectBase<DerivedI>& mat, const std::vector<std::pair<int, int>>& edges);

// ���������
template <typename Index>
std::int64_t encodeEdge(const Index vaIdx, const Index vbIdx);

// ���������
template <typename Index>
std::int64_t encodeUedge(const Index vaIdx, const Index vbIdx);

// ��������Ƭ
template <typename Index>
std::uint64_t encodeTriangle(const Index vaIdx, const Index vbIdx, const Index vcIdx);

// ��������߱���
std::pair<int, int> decodeEdge(const std::int64_t code);

// ��������Ƭ����
std::vector<int> decodeTrianagle(const std::uint64_t code);



/////////////////////////////////////////////////////////////////////////////////////////////////// �������ԣ�



/////////////////////////////////////////////////////////////////////////////////////////////////// ����༭��

template<typename T>
bool smoothCircuit2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& circuit, const float param);


template <typename T>
void concatMeshMat(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers1, const Eigen::MatrixXi& tris1);


template <typename IndexT>
bool removeTris(Eigen::MatrixXi& trisOut, const Eigen::MatrixXi& tris, const std::vector<IndexT>& sickTriIdxes);






// ��ʱδ�����ʵ�֣�
#include "temp.tpp"