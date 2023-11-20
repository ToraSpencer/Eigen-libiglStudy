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
bool getEdges(Eigen::MatrixXi& edges, const Eigen::MatrixBase<DerivedI>& tris);
 

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


// 输入有向边数据，得到顶点邻接矩阵：
/*
	bool edges2adjMat(
			Eigen::SparseMatrix<int>& adjSM_eCount,			权重为该有向边重复次数的邻接矩阵；
			Eigen::SparseMatrix<int>& adjSM_eIdx,				权重为该有向边索引的邻接矩阵；
			const Eigen::MatrixBase<DerivedI>& edges			有向边矩阵；
			)

*/
template <typename DerivedI>
bool edges2adjMat(Eigen::SparseMatrix<int>& adjSM_eCount,\
	Eigen::SparseMatrix<int>& adjSM_eIdx, const Eigen::MatrixBase<DerivedI>& edges)
{
	assert((edges.rows() > 0) && "Assert!!! edges should'nt be empty.");
	assert((2 == edges.cols()) && "Assert!!! edges should be a 2 column matrix.");

	const int versCount = edges.maxCoeff() + 1;
	const int edgesCount = edges.rows();
	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;
	smElems.reserve(edgesCount);
	smElems_weighted.reserve(edgesCount);
	for (int i = 0; i < edgesCount; ++i)
	{
		smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});
		smElems_weighted.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), i});
	}

	adjSM_eCount.resize(versCount, versCount);
	adjSM_eIdx.resize(versCount, versCount);
	adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());									// 权重为该有向边重复的次数；
	adjSM_eIdx.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());		// 权重为该有向边的索引；


	return true;
}


// 求三角网格的不同权重的邻接矩阵 
/*
	template<typename int>
	bool tris2adjMat(
		Eigen::SparseMatrix<int>& adjSM_eCount,			权重为边数目的邻接矩阵
		Eigen::SparseMatrix<int>& adjSM_eIdx,				权重为边索引的邻接矩阵；
		const Eigen::MatrixXi& tris,											输入网格的三角片		
		const Eigen::MatrixXi& edges0 = Eigen::MatrixXi{}			输入网格的边信息；默认情形下此项为空，在此函数中计算边数据；
		);

*/ 
template<typename DerivedI>
bool tris2adjMat(Eigen::SparseMatrix<int>& adjSM_eCount, \
	Eigen::SparseMatrix<int>& adjSM_eIdx, const Eigen::MatrixBase<DerivedI>& tris)
{
	Eigen::MatrixXi edges;
	getEdges(edges, tris);
	edges2adjMat(adjSM_eCount, adjSM_eIdx, edges);

	return true;
}


// 边缘有向边，重载1：
template<typename DerivedI>
bool bdryEdges(Eigen::MatrixXi& bdrys, const Eigen::MatrixBase<DerivedI>& tris)
{
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	tris2adjMat(adjSM_eCount, adjSM_eIdx, tris);
	Eigen::SparseMatrix<int> adjSM_tmp, adjSM_ueCount;

	// 1. 确定边缘无向边：
	spMatTranspose(adjSM_tmp, adjSM_eCount);
	adjSM_ueCount = adjSM_eCount + adjSM_tmp;				// 无向边邻接矩阵，权重是该边的重复次数；

	std::vector<std::pair<int, int>> bdryUeVec;						// 所有边缘无向边；约定无向边表示为前大后小的索引对；
	traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
		{
			if (1 == iter.value() && iter.row() > iter.col())
				bdryUeVec.push_back({ iter.row(), iter.col() });
		});

	// 2. 由边缘无向边确定边缘有向边：
	std::vector<std::pair<int, int>> bdryVec;
	bdryVec.reserve(bdryUeVec.size());
	for (const auto& pair : bdryUeVec)
	{
		if (adjSM_eCount.coeff(pair.first, pair.second) > 0)
			bdryVec.push_back({ pair.first, pair.second });
		else
			bdryVec.push_back({ pair.second, pair.first });
	}
	edges2mat(bdrys, bdryVec);

	return true;
}


// 边缘有向边， 重载2：
template<typename DerivedI>
bool bdryEdges(Eigen::MatrixXi& bdrys, std::vector<int>& bdryTriIdxes, \
	const Eigen::MatrixBase<DerivedI>& tris)
{
	/*
		bool bdryEdges(
				Eigen::PlainObjectBase<DerivedI>& bdrys,					bdrysCount * 2，边缘边数据，是有向边；
				std::vector<int>& bdryTriIdxes,										包含边缘边的三角片索引；
				const Eigen::PlainObjectBase<DerivedI>& tris				trisCount * 3, 网格三角片数据；
				)

	*/

	// 1. 查找网格的边缘有向边；
	if (!bdryEdges(bdrys, tris))
		return false;

	if (0 == bdrys.rows())
		return true;

	const unsigned trisCount = tris.rows();

	// 2. 确认边缘非流形边所在的三角片索引：
	Eigen::MatrixXi edgeAs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi edgeBs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi edgeCs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi vaIdxes = tris.col(0);
	Eigen::MatrixXi vbIdxes = tris.col(1);
	Eigen::MatrixXi vcIdxes = tris.col(2);

	// 3. 生成有向边编码-三角片字典：
	std::unordered_multimap<std::int64_t, unsigned> edgesMap;			// 一条有向边正常情况下只关联一个三角片，若关联多个，则是非流形边；
	for (unsigned i = 0; i < trisCount; ++i)
	{
		std::int64_t codeA = encodeEdge(vbIdxes(i), vcIdxes(i));
		std::int64_t codeB = encodeEdge(vcIdxes(i), vaIdxes(i));
		std::int64_t codeC = encodeEdge(vaIdxes(i), vbIdxes(i));
		edgesMap.insert({ codeA, i });
		edgesMap.insert({ codeB, i });
		edgesMap.insert({ codeC, i });
	}

	// 4. 在字典中查找所有边缘边所在的三角片索引：
	for (unsigned i = 0; i < bdrys.rows(); ++i)
	{
		std::int64_t code = encodeEdge(bdrys(i, 0), bdrys(i, 1));
		auto iter = edgesMap.find(code);
		unsigned codeCounts = edgesMap.count(code);
		for (unsigned k = 0; k < codeCounts; ++k)
			bdryTriIdxes.push_back((iter++)->second);
	}

	return true;
}



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

// 计算网格中所有顶点的1ring信息——1ring的顶点索引、三角片索引；
template<typename DerivedV, typename DerivedI>
bool get1ring(std::vector<std::unordered_set<int>>& vIdx1ring, std::vector<std::unordered_set<int>>& tIdx1ring,\
	const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixBase<DerivedI>& tris)
{
	const int versCount = vers.rows();
	const int trisCount = tris.rows();
	vIdx1ring.resize(versCount);
	tIdx1ring.resize(versCount);
	for (int i = 0; i < trisCount; ++i)
	{
		int vaIdx = static_cast<int>(tris(i, 0));
		int vbIdx = static_cast<int>(tris(i, 1));
		int vcIdx = static_cast<int>(tris(i, 2));
		vIdx1ring[vaIdx].insert(vbIdx);
		vIdx1ring[vaIdx].insert(vcIdx);
		vIdx1ring[vbIdx].insert(vaIdx);
		vIdx1ring[vbIdx].insert(vcIdx);
		vIdx1ring[vcIdx].insert(vaIdx);
		vIdx1ring[vcIdx].insert(vbIdx);
		tIdx1ring[vaIdx].insert(i);
		tIdx1ring[vbIdx].insert(i);
		tIdx1ring[vcIdx].insert(i);
	}
	 
	return true;
}

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
	const Eigen::MatrixBase<DerivedV>& vers, \
	const Eigen::MatrixBase<DerivedI>& tris, const int triIdx);


template <typename DerivedV>
bool triangleGrowSplitMesh(std::vector<DerivedV>& meshesVersOut, \
	std::vector<Eigen::MatrixXi>& meshesTrisOut, \
	const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris);


// 三角片区域生长得到最外层网格；
/*
	template <typename DerivedVo, typename DerivedVi>
	bool triangleGrowOuterSurf(
		Eigen::PlainObjectBase<DerivedVo>& versOut, \		输出网格
		Eigen::MatrixXi& trisOut, 
		const Eigen::MatrixBase<DerivedVi>& vers, \			输入网格
		const Eigen::MatrixXi& tris,
		const bool blRemIso,													是否去除孤立顶点
		const int startIdx															区域生长起始的三角片索引；
		)

*/
template <typename DerivedVo, typename DerivedVi>
bool triangleGrowOuterSurf(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::MatrixXi& trisOut, const Eigen::MatrixBase<DerivedVi>& vers,\
	const Eigen::MatrixXi& tris, const bool blRemIso = true, const int startIdx = -1);


// robustTriangleGrowOuterSurf()——提升了鲁棒性的三角片生长，避免了某些反向三角片的干扰； 
/*
	template <typename DerivedVo, typename DerivedVi>
	bool robustTriangleGrowOuterSurf(
		Eigen::PlainObjectBase<DerivedVo>& versOut, \		输出网格
		Eigen::MatrixXi& trisOut,
		const Eigen::MatrixBase<DerivedVi>& vers, \			输入网格
		const Eigen::MatrixXi& tris,
		const bool blRemIso,													是否去除孤立顶点
		)

*/
template <typename DerivedVo, typename DerivedVi>
bool robustTriangleGrowOuterSurf(Eigen::PlainObjectBase<DerivedVo>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::MatrixBase<DerivedVi>& vers, const Eigen::MatrixXi& tris, const bool blRemIso = true);

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
std::vector<unsigned> checkIsoVers(const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris);


// 去除网格中的孤立顶点：
template <typename DerivedVo, typename DerivedVi>
bool removeIsoVers(Eigen::PlainObjectBase<DerivedVo>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::MatrixBase<DerivedVi>& vers, const Eigen::MatrixXi& tris, \
	const std::vector<unsigned>& isoVerIdxes);


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