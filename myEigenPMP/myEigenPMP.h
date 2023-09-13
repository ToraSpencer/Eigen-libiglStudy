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


/////////////////////////////////////////////////////////////////////////////////////////////////// 表象转换接口：

// std::vector<std::pair<int, int>>表示的边数据转换为矩阵表示：
template <typename DerivedI>
void edges2mat(Eigen::PlainObjectBase<DerivedI>& mat, const std::vector<std::pair<int, int>>& edges);

// 编码有向边
template <typename Index>
std::int64_t encodeEdge(const Index vaIdx, const Index vbIdx);

// 编码无向边
template <typename Index>
std::int64_t encodeUedge(const Index vaIdx, const Index vbIdx);

// 编码三角片
template <typename Index>
std::uint64_t encodeTriangle(const Index vaIdx, const Index vbIdx, const Index vcIdx);

// 解码有向边编码
std::pair<int, int> decodeEdge(const std::int64_t code);

// 解码三角片编码
std::vector<int> decodeTrianagle(const std::uint64_t code);



/////////////////////////////////////////////////////////////////////////////////////////////////// 网格属性：



/////////////////////////////////////////////////////////////////////////////////////////////////// 网格编辑：

template<typename T>
bool smoothCircuit2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& circuit, const float param);


template <typename T>
void concatMeshMat(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers1, const Eigen::MatrixXi& tris1);


template <typename IndexT>
bool removeTris(Eigen::MatrixXi& trisOut, const Eigen::MatrixXi& tris, const std::vector<IndexT>& sickTriIdxes);






// 暂时未整理的实现：
#include "temp.tpp"