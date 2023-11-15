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


// 目录：
/*
	1. 表象转换
	2. 图形属性
	3. 点云编辑
	4. 三角网格编辑
	5. 区域生长算法
	6. 网格缺陷检查和修复

*/
 

/////////////////////////////////////////////////////////////////////////////////////////////////// 表象转换接口：

// std::vector<std::pair<int, int>>表示的边数据转换为矩阵表示：
template <typename DerivedI>
void edges2mat(Eigen::PlainObjectBase<DerivedI>& mat, \
	const std::vector<std::pair<int, int>>& edges);

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


/////////////////////////////////////////////////////////////////////////////////////////////////// 生成基础图元：

// 顶点插值得到直线段点云
template <typename DerivedVo, typename ScalarVi>
bool interpolateToLine(Eigen::PlainObjectBase<DerivedVo>& vers, \
	const Eigen::Matrix<ScalarVi, 1, 3>& start, const Eigen::Matrix<ScalarVi, 1, 3>& end, \
	const float SR, const bool SE = true);


// 生成XOY平面内的圆圈点集：
template <typename T>
bool getCircleVers(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, \
	const unsigned versCount = 20);


// 对环路点集进行不插点三角剖分——手动缝合方法； 
bool circuitGetTris(Eigen::MatrixXi& tris, const std::vector<int>& indexes, const bool regularTri);



/////////////////////////////////////////////////////////////////////////////////////////////////// 图形属性：

// 输入网格三角片数据，得到有向边数据：
template <typename DerivedI>
bool getEdges(Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedI>& tris);
 

// 生成有序环路的边数据；
template <typename DerivedI, typename IndexType>
bool getSortedLoopEdges(Eigen::PlainObjectBase<DerivedI>& edges, const IndexType versCount);


// 生成边-边编码的映射表：
template <typename DerivedI>
bool getEdgeIdxMap(std::unordered_multimap<std::int64_t, int>& map, \
	const Eigen::PlainObjectBase<DerivedI>& edges);


template <typename DerivedV, typename DerivedF, typename DerivedN>
bool getTriNorms(Eigen::PlainObjectBase<DerivedN>& triNorms,\
	const Eigen::MatrixBase<DerivedV>& vers, \
	const Eigen::MatrixBase<DerivedF>& tris);


template <typename DerivedV>
void getEdgeArrows(Eigen::PlainObjectBase<DerivedV>& arrows, \
	const Eigen::MatrixXi& edges, \
	const Eigen::PlainObjectBase<DerivedV>& vers);


template <typename DerivedV, typename DerivedF, typename DerivedL>
void squared_edge_lengths(Eigen::PlainObjectBase<DerivedL>& lenMat2,\
	const Eigen::PlainObjectBase<DerivedV>& vers,
	const Eigen::MatrixBase<DerivedF>& tris);

// 计算网格中三角片三个角的余切值；
template <typename DerivedC, typename DerivedV>
void trisCotValues(Eigen::PlainObjectBase<DerivedC>& cotValues, \
	const Eigen::PlainObjectBase<DerivedV>& vers,
	const Eigen::MatrixXi& tris);

template<typename T, typename DerivedV, typename DerivedI>
bool trianglesBarycenter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& barys, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);

// 计算网格所有三角片的归一化法向量
template<typename DerivedN, typename DerivedV, typename DerivedI>
bool trianglesNorm(Eigen::PlainObjectBase<DerivedN>& triNorms, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);

template<typename DerivedV, typename DerivedI>
bool trianglesPlane(Eigen::MatrixXd& planeCoeff, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);

// 计算三角网格每个三角片的面积：
template<typename DerivedA, typename DerivedV, typename DerivedI>
bool trisArea(Eigen::PlainObjectBase<DerivedA>& trisAreaVec, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::MatrixBase<DerivedI>& tris);

// 计算水密的三角网格的体积；
template<typename DerivedV, typename DerivedI>
double meshVolume(const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);


template<typename Index>
bool adjMatrix(Eigen::SparseMatrix<Index>& adjSM_eCount, \
	Eigen::SparseMatrix<Index>& adjSM_eIdx, const Eigen::MatrixXi& tris, \
	const Eigen::MatrixXi& edges0);


template<typename DerivedI>
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, \
	const Eigen::PlainObjectBase<DerivedI>& tris);

template<typename DerivedI>
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, \
	std::vector<int>& bdryTriIdxes, const Eigen::PlainObjectBase<DerivedI>& tris);

using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool buildAdjacency(Eigen::MatrixXi& ttAdj_mnEdge, \
	std::vector<tVec>& ttAdj_nmnEdge, std::vector<tVec>& ttAdj_nmnOppEdge, \
	const Eigen::PlainObjectBase<DerivedI>& tris);

template <typename DerivedI>
bool buildAdjacency_new(Eigen::MatrixXi& ttAdj_mnEdge, \
	std::vector<tVec>& ttAdj_nmnEdge, std::vector<tVec>& ttAdj_nmnOppEdge, \
	const Eigen::PlainObjectBase<DerivedI>& tris);

// 计算余切权的laplace矩阵；
template<typename Tl, typename DerivedV>
bool cotLaplacian(Eigen::SparseMatrix<Tl>& L, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::MatrixXi& tris);


/////////////////////////////////////////////////////////////////////////////////////////////////// 点云编辑：

// 基于laplacian的回路光顺
template<typename T>
bool smoothCircuit2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& circuit, \
	const float param);


// XOY平面内的回路排序：排序成从Z轴正向看顶点随回路逆时针或顺时针增大；
/*
	bool sortLoop2D(
			Eigen::PlainObjectBase<DerivedV>& loop,			待排序的回路
			const bool blCounterClockWise,								是否排序成逆时针
			const int startIdx														指定输入回路中索引为startIdx的顶点作为整理后回路的起点；
			)

*/
template<typename DerivedV>
bool sortLoop2D(Eigen::PlainObjectBase<DerivedV>& loop, \
	const bool blCounterClockWise = true, 	const float thresholdDis = 10);


/////////////////////////////////////////////////////////////////////////////////////////////////// 三角网格编辑：

template <typename T>
void concatMeshMat(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers,\
	Eigen::MatrixXi& tris, 	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers1, \
	const Eigen::MatrixXi& tris1);


template <typename IndexT>
bool removeTris(Eigen::MatrixXi& trisOut, const Eigen::MatrixXi& tris, \
	const std::vector<IndexT>& sickTriIdxes);


// laplace光顺 
/*
	bool laplaceFaring(
		Eigen::PlainObjectBase<DerivedVo>& versOut,					输出网格顶点
		const Eigen::PlainObjectBase<DerivedVi>& vers,				输入网格顶点
		const Eigen::MatrixXi& tris,												输入网格面片
		const float deltaLB,															扩散速度参数
		const unsigned loopCount,												光顺次数
		const std::vector<int>& fixedVerIdxes							限定保持不变的顶点的索引
		)
*/
template <typename DerivedVo, typename DerivedVi>
bool laplaceFaring(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	const Eigen::PlainObjectBase<DerivedVi>& vers, \
	const Eigen::MatrixXi& tris, const float deltaLB, const unsigned loopCount, \
	const std::vector<int>& fixedVerIdxes = std::vector<int>{});



/////////////////////////////////////////////////////////////////////////////////////////////////// 区域生长算法：

template <typename DerivedV, typename DerivedI>
bool triangleGrow(Eigen::PlainObjectBase<DerivedV>& versOut,\
	Eigen::PlainObjectBase<DerivedI>& trisOut, 	\
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx);


template <typename DerivedV>
bool triangleGrowSplitMesh(std::vector<DerivedV>& meshesVersOut, \
	std::vector<Eigen::MatrixXi>& meshesTrisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris);

template <typename T>
bool triangleGrowOuterSurf(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, \
	Eigen::MatrixXi& trisOut, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers,\
	const Eigen::MatrixXi& tris, const bool blRemIso, const int startIdx);

template <typename T>
bool robustTriangleGrowOuterSurf(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, const bool blRemIso);

template <typename DerivedI>
bool uEdgeGrow(std::vector<Eigen::MatrixXi>& circs, std::vector<Eigen::MatrixXi >& segs, \
	const Eigen::PlainObjectBase<DerivedI>& uEdges);

template <typename Index>
int simplyConnectedRegion(const Eigen::SparseMatrix<Index>& adjSM,
	Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount);

template <typename T>
int simplyConnectedRegion(std::vector<int>& connectedLabels, std::vector<int>& connectedCount, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris);

template <typename DerivedI>
int simplyTrisConnectedRegion(Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount, \
	const Eigen::PlainObjectBase<DerivedI>& tris);

template <typename DerivedV, typename DerivedI>
bool simplyConnectedLargest(Eigen::PlainObjectBase<DerivedV>& versOut, \
	Eigen::PlainObjectBase<DerivedI>& trisOut, 	const Eigen::MatrixBase<DerivedV>& vers, \
	const Eigen::MatrixBase<DerivedI>& tris);

template <typename T>
bool simplyConnectedSplitMesh(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& compVers, \
	std::vector<Eigen::MatrixXi>& compTris, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const Eigen::MatrixXi& tris);



/////////////////////////////////////////////////////////////////////////////////////////////////// 网格缺陷检查和修复
template <typename DerivedI>
void findRepTris(std::vector<int>& repIdxes, const Eigen::PlainObjectBase<DerivedI>& tris);

template <typename DerivedV, typename DerivedI>
int findOverLapTris(std::vector<std::pair<int, int>>& opTrisPairs, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris, \
	const int triIdx, const double thetaThreshold);

template<typename DerivedI>
int nonManifoldEdges(Eigen::MatrixXi& nmnEdges, const Eigen::PlainObjectBase<DerivedI>& tris);


template<typename DerivedI>
int nonManifoldEdges(Eigen::MatrixXi& nmnEdges, std::vector<std::pair<std::pair<int, int>, int>>& nmnInfos, \
	const Eigen::PlainObjectBase<DerivedI>& tris);


template <typename T>
int nonManifoldUEs(Eigen::MatrixXi& nmnUedges, const Eigen::SparseMatrix<T>& adjSM_eCount, \
	const Eigen::SparseMatrix<T>& adjSM_ueCount);


template <typename DerivedI>
int nonManifoldUEs(Eigen::MatrixXi& nmnUedges, const Eigen::PlainObjectBase<DerivedI>& tris);


template<typename DerivedI>
int nonManifoldUEs(Eigen::MatrixXi& nmnUedges, std::vector<std::pair<int, std::pair<int, int>>>& nmnEdgeInfos, \
	const Eigen::PlainObjectBase<DerivedI>& tris);


template<typename DerivedV>
int nonManifoldVers(std::vector<int>& nmnVerIdxes, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris);

template <typename IndexType>
bool fillSmallHoles(Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic>& newTris, \
	const std::vector<Eigen::Matrix<IndexType, Eigen::Dynamic, 1>>& holes);

template <typename DerivedV>
std::vector<unsigned> checkIsoVers(const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris);

template <typename DerivedV>
bool removeIsoVers(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const std::vector<unsigned>& isoVerIdxes);

template <typename DerivedV>
Eigen::VectorXi checkDegEdges(const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const double threshold);

template <typename DerivedV>
int mergeDegEdges(Eigen::PlainObjectBase<DerivedV>& newVers, Eigen::MatrixXi& newTris, \
	const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const Eigen::VectorXi& degEdgeFlags);

template <typename DerivedV, typename DerivedI>
int correctTriDirs(Eigen::PlainObjectBase<DerivedI>& trisOut, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris, const double thetaThreshold);

template<typename T>
int removeSickDupTris(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);

template <typename DerivedV, typename DerivedI>
int findHoles(std::vector<std::vector<int>>& holes, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);

template <typename DerivedI>
bool checkSickTris(std::vector<int>& sickIdxes, const Eigen::PlainObjectBase<DerivedI>& tris);





// 暂时未整理的实现：
#include "temp.tpp"