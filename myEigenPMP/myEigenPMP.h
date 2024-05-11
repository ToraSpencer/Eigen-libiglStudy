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
	2. 生成基础图元
	3. 图形属性
	4. 点云编辑
	5. 三角网格编辑
	6. 区域生长算法
	7. 网格缺陷检查和修复

*/
 

/////////////////////////////////////////////////////////////////////////////////////////////////// 前置声明：
template<typename DerivedV>
int removeSickDupTris(const Eigen::MatrixBase<DerivedV>& vers, Eigen::MatrixXi& tris);

template <typename DerivedVo, typename DerivedVi>
bool removeIsoVers(Eigen::PlainObjectBase<DerivedVo>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::MatrixBase<DerivedVi>& vers, const Eigen::MatrixXi& tris, \
	const std::vector<unsigned>& isoVerIdxes);

template <typename DerivedV>
std::vector<unsigned> checkIsoVers(const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris);

using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool buildAdjacency(Eigen::MatrixXi& ttAdj_mnEdge, std::vector<tVec>& ttAdj_nmnEdge, \
	std::vector<tVec>& ttAdj_nmnOppEdge, const Eigen::MatrixBase<DerivedI>& tris);


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
template <typename DerivedVo, typename DerivedVs, typename DerivedVe>
bool interpolateToLine(Eigen::PlainObjectBase<DerivedVo>& vers, \
	const Eigen::MatrixBase<DerivedVs>& start, const Eigen::MatrixBase<DerivedVe>& end, \
	const float SR, const bool SE = true)
{
	vers.resize(0, 0);
	if (vers.rows() > 0)
		return false;

	Eigen::RowVector3d startD = start.array().cast<double>();
	Eigen::RowVector3d endD = end.array().cast<double>();
	Eigen::RowVector3d dir = endD - startD;
	float length = dir.norm();
	dir.normalize();
	if (length <= SR)
		return true;

	if (SE)
		matInsertRows(vers, startD);

	float lenth0 = 0;
	for (unsigned i = 1; (length - lenth0) > 0.8 * SR; i++)			// 确保最后一个点距离终点不能太近。
	{
		Eigen::RowVector3d temp = startD + SR * i * dir;
		matInsertRows(vers, temp);
		lenth0 = SR * (i + 1);		// 下一个temp的长度。
	}

	if (SE)
		matInsertRows(vers, endD);

	return true;
};


// 生成XOY平面内的圆圈点集：
template <typename T>
bool getCircleVers(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, \
	const unsigned versCount = 20);


// 对环路点集进行不插点三角剖分——手动缝合方法； 
bool circuitGetTris(Eigen::MatrixXi& tris, const std::vector<int>& indexes, const bool regularTri);



/////////////////////////////////////////////////////////////////////////////////////////////////// 图形属性：

// 得到三角网格的有向边数据
template <typename DerivedI>
bool getEdges(Eigen::MatrixXi& edges, const Eigen::MatrixBase<DerivedI>& tris)
{
	/*
		bool getEdges(
				Eigen::MatrixXi& edges,
				const Eigen::PlainObjectBase<DerivedI>& tris
				)

		三条边在edges中的排列顺序为[bc; ca; ab]；其所对的顶点标记corner分别为0, 1, 2，即a,b,c;
		边索引到三角片索引的映射——边eIdx0所在的三角片triIdx0 == eIdx0 % trisCount;
		三角片triIdx0和corner0到边索引的映射——eIdx0 = corner0 * trisCount + triIdx0;
		边索引到corner的映射——corner0 = eIdx0 / trisCount;
	*/
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	edges = Eigen::MatrixXi::Zero(edgesCount, 2);
	Eigen::MatrixXi vaIdxes = tris.col(0).array().cast<int>();
	Eigen::MatrixXi vbIdxes = tris.col(1).array().cast<int>();
	Eigen::MatrixXi vcIdxes = tris.col(2).array().cast<int>();
	edges.block(0, 0, trisCount, 1) = vbIdxes;
	edges.block(trisCount, 0, trisCount, 1) = vcIdxes;
	edges.block(trisCount * 2, 0, trisCount, 1) = vaIdxes;
	edges.block(0, 1, trisCount, 1) = vcIdxes;
	edges.block(trisCount, 1, trisCount, 1) = vaIdxes;
	edges.block(trisCount * 2, 1, trisCount, 1) = vbIdxes;

	return true;
}
 

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
bool trianglesBarycenter(\
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& barys, \
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
	adjSM_eCount.resize(0, 0);
	adjSM_eIdx.resize(0, 0);
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


// buildAdjacency()——计算网格的三角片邻接信息：
using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool buildAdjacency(Eigen::MatrixXi& ttAdj_mnEdge, std::vector<tVec>& ttAdj_nmnEdge, \
	std::vector<tVec>& ttAdj_nmnOppEdge, const Eigen::MatrixBase<DerivedI>& tris)
{
	/*
		bool buildAdjacency(
					const Eigen::MatrixXi& tris,									输入的三角片数据

					Eigen::MatrixXi& ttAdj_mnEdge,							三角片的非边缘流形有向边邻接的三角片索引；
																									trisCount * 3
																									(i, 0)为索引为i的三角片的非边缘流形有向边(vb, vc)邻接的三角片的索引；

					std::vector<ttTuple>& ttAdj_nmnEdge,				三角片的非流形有向边所在的三角片索引；
																									size == trisCount;
																									每个ttTuple元素包含三个int向量；
																									索引为i的元素的第0个向量，为索引为i的三角片的非流形有向边(vb, vc)所在的所有三角片索引；


					std::vector<ttTuple>& ttAdj_nmnOppEdge		三角片的非流形有向边邻接的三角片索引（即对边所在的三角片的索引）；
																									size == trisCount;
																									索引为i的元素的第0个向量，为索引为i的三角片的非流形有向边(vb, vc)邻接的所有三角片索引；
					)

	*/

	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	Eigen::MatrixXi edges;							// 有向边数据；
	Eigen::MatrixXi nmnEdges;					// 非流形有向边

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_eIdx的TRIPLET数据；
	std::unordered_multimap<std::int64_t, int> edgeMap;					// 边编码——64位整型数数表示的边数据（两个端点索引）；


	// 1. 求基本的边信息、邻接矩阵：
	Eigen::SparseMatrix<int> adjSM;								// 邻接矩阵；
	Eigen::SparseMatrix<int> adjSM_eCount;					// 邻接矩阵，索引为该有向边重复的次数；
	Eigen::SparseMatrix<int> adjSM_eIdx;						// 邻接矩阵，元素若对应流形边则为该边索引，若对应非流形边则为该边所有索引的和；
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;				// adjSM_eIdx的转置，表示对边的信息；
	Eigen::SparseMatrix<int> adjSM_ueCount;				// 权重为无向边ij关联的三角片数量；
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;
	{

		// edges == [ea; eb; ec] == [vbIdxes, vcIdxes; vcIdxes, vaIdxes; vaIdxes, vbIdxes];
		/*
			若网格是流形网格，则edges列表里每条边都是unique的；
			若存在非流形边，则非流形边在edges里会重复存储，实际边数量也会少于edges行数；
			可以说使用这种representation的话，非流形边有不止一个索引，取决包含该非流形边的三角片数量；
		*/

		// 三角片三条边的索引：teIdx == [eaIdx, ebIdx, ecIdx] == [(0: trisCount-1)', (trisCount: 2*trisCount-1)', (2*trisCount, 3*trisCount-1)'];
		/*
			teIdx(i, j) = trisCount *j + i;
		*/

		// 1.1 生成有向边数据
		getEdges(edges, tris);

		// 1.2 生成三种有向边邻接矩阵：adjSM, adjSM_eCount, adjSM_eCount;
		edges2adjMat(adjSM_eCount, adjSM_eIdx, edges);
		adjSM = adjSM_eCount;																		// 有向边邻接矩阵；
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		若adjSM_eIdx(i, j)对应的是流形有向边，该权重值为该边的索引；若是非流形有向边，则该边对应两个边索引，权重为两个索引之和；
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

		//	1.3 生成无向边邻接矩阵：
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// 无向边邻接矩阵；
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.4 生成非边缘流形无向边邻接矩阵：
		smElems.clear();
		smElems.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				if (2 == iter.value())
					smElems.push_back(Eigen::Triplet<int>(iter.row(), iter.col(), 1));
			});
		adjSM_ueMN_NB.resize(versCount, versCount);
		adjSM_ueMN_NB.setFromTriplets(smElems.begin(), smElems.end());
		smElems.clear();

		//	1.5 生成非边缘流形有向边、及其对边的邻接矩阵
		adjSM_MN_NB = adjSM + adjSM_ueMN_NB;					// adjSM & adjSM_ueMN_NB
		adjSM_MN_NB.prune([](const Eigen::Index& row, const Eigen::Index& col, const int& value)->bool
			{
				if (2 == value)			// 删除非流形边
					return true;
				else
					return false;
			});
		adjSM_MN_NB /= 2;
		spMatTranspose(adjSM_MN_NB_opp, adjSM_MN_NB);
	}


	// 2. 确定非边缘流形有向边、及其对边的索引；
	std::vector<int> edgesIdx_MN_NB;							// 非边缘流形有向边的索引；
	std::vector<int> edgesIdx_MN_NB_opp;					// 非边缘流形有向边的对边的索引；
	{
		unsigned edgesCount_MN_NB = adjSM_MN_NB.sum();
		edgesIdx_MN_NB.reserve(edgesCount_MN_NB);					// 非边缘流形有向边索引；
		edgesIdx_MN_NB_opp.reserve(edgesCount_MN_NB);			// 非边缘流形有向边的对边的索引；

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB.push_back(adjSM_eIdx.coeffRef(iter.row(), iter.col()));
			});

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB_opp.push_back(adjSM_eIdx_opp.coeffRef(iter.row(), iter.col()));
			});
	}


	// 3. 流形边-三角片邻接关系，三角片邻接关系：
	std::vector<int> etInfo;						// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
	Eigen::VectorXi etAdj_mnEdge;							// 边所在三角片的索引，非流形边或边缘边写-1；
	{
		// 3.1 生成边索引 - 三角片索引映射表etInfo;
		etInfo.resize(edgesCount);							// etInfo(i)是索引为i的边所在的三角片的索引；
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

		// 3.2 求边所在三角片的索引；
		etAdj_mnEdge = -Eigen::VectorXi::Ones(edgesCount);					// 边所在三角片的索引，非流形边或边缘边写-1；
		for (unsigned i = 0; i < edgesIdx_MN_NB.size(); ++i)
		{
			const int& edgeIdx = edgesIdx_MN_NB[i];
			const int& edgesIdxOpp = edgesIdx_MN_NB_opp[i];

			if (edgesIdxOpp < 0 || edgesIdxOpp >= edgesCount)
				std::cout << "pause" << std::endl;

			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//	3.3 求三角片邻接矩阵；ttAdj_mnEdge(i, :)是索引为i的三角片三条边邻接的三个三角片，若其中有非流形边或边缘边则写为-1；
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. 计算非流形边（关联两个三角片的有向边）信息
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_map;				// 非流形有向边的索引 - 该边对应的所有索引；
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_opp_map;		// 非流形有向边的索引 - 该边的对边对应的所有索引；
	{
		//	4.1 遍历邻接矩阵找出所有非流形边（关联两个三角片的有向边）；
		std::vector<int> edgeIdx_nmn;						// 非流形有向边索引；
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)) > 1)
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();

		int neCount = edgeIdx_nmn.size();						// 非流形有向边个数；
		subFromIdxVec(nmnEdges, edges, edgeIdx_nmn);

		// 4.2 建立边编码-边索引的哈希表；
		for (int i = 0; i < edges.rows(); ++i)
		{
			std::int64_t eCode = encodeEdge(edges(i, 0), edges(i, 1));
			edgeMap.insert({ eCode, i });
		}

		// 4.3 在哈希表中搜索所有非流形边，建立非流形边及其对边的（边——边索引）映射关系；
		for (int i = 0; i < neCount; ++i)
		{
			// 注意！！存在两种以外情况： a. 当前边是非流形边，对边是流形边；	b. 当前边是非流形边，对边不存在；
			int neIdx = edgeIdx_nmn[i];		// 当前非流形边索引；
			int vaIdx = nmnEdges(i, 0);
			int vbIdx = nmnEdges(i, 1);
			std::int64_t eCode = encodeEdge(vaIdx, vbIdx);
			std::int64_t oppEcode = encodeEdge(vbIdx, vaIdx);
			auto iter = edgeMap.find(eCode);
			auto oppIter = edgeMap.find(oppEcode);

			int otherIdx = (iter->second == neIdx) ? ((++iter)->second) : (iter->second);
			edgeIdx_nmn_map.insert({ neIdx, std::vector<int>{neIdx, otherIdx} });

			if (edgeMap.end() == oppIter)			// 当前边是非流形边，对边不存在；
			{
				edgeIdx_nmn_opp_map.insert({ neIdx, std::vector<int>{} });
				continue;
			}
			int oppIdx1 = oppIter->second;
			int oppIdx2 = (++oppIter)->second;
			edgeIdx_nmn_opp_map.insert({ neIdx, std::vector<int>{oppIdx1, oppIdx2} });
		}
	}


	// 5. 求含有非流形边的三角片的邻接关系：
	{
		ttAdj_nmnEdge.resize(trisCount);
		for (auto& vec : ttAdj_nmnEdge)
			vec.resize(3);

		for (const auto& pair : edgeIdx_nmn_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;

			std::vector<int> trisIdx_nmn = pair.second;
			for (auto& index : trisIdx_nmn)
				index = etInfo[index];

			ttAdj_nmnEdge[row][col] = trisIdx_nmn;
		}

		ttAdj_nmnOppEdge.resize(trisCount);
		for (auto& vec : ttAdj_nmnOppEdge)
			vec.resize(3);
		for (const auto& pair : edgeIdx_nmn_opp_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;

			std::vector<int> trisIdx_nmn_opp = pair.second;
			for (auto& index : trisIdx_nmn_opp)
				index = etInfo[index];

			ttAdj_nmnOppEdge[row][col] = trisIdx_nmn_opp;
		}
	}

	return true;
}


// getTrisAdjacency()——计算网格的三角片邻接信息：
using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool getTrisAdjacency(Eigen::MatrixXi& ttAdj_mnEdge, std::vector<tVec>& ttAdj_nmnEdge, \
	std::vector<tVec>& ttAdj_nmnOppEdge, const Eigen::MatrixBase<DerivedI>& edges, \
	const Eigen::SparseMatrix<int>& adjSM_eCount, 	const Eigen::SparseMatrix<int>& adjSM_eIdx)
{ 
	const int edgesCount = edges.rows();
	const int trisCount = edgesCount / 3;

	// 1. 求基本的边信息、邻接矩阵：  
	Eigen::SparseMatrix<int> adjSM_eCount_opp; 
	Eigen::SparseMatrix<int> adjSM_ueCount;				// 权重为无向边ij关联的三角片数量； 
	Eigen::SparseMatrix<int> adjSM_MN;						// 流形有向边（非边缘）邻接矩阵，权重为1； 
	std::unordered_set<int> bdryIdxes;									// 边缘有向边的索引；
	{ 
		spMatTranspose(adjSM_eCount_opp, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + adjSM_eCount_opp;						// 无向边邻接矩阵，权重为无向边ij关联的三角片数量； 						 

		//	1.3 生成流形有向边（非边缘）的邻接矩阵
		adjSM_MN = adjSM_ueCount;
		traverseSparseMatrix(adjSM_MN, [&adjSM_MN](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});
		adjSM_MN.prune([&adjSM_eCount](const Eigen::Index& row, const Eigen::Index& col, const int& value)->bool
			{
				if (1 == adjSM_eCount.coeff(row, col) && 1 == adjSM_eCount.coeff(col, row))
					return true;								// 返回true的元素被保留；
				else
					return false;
			});															// 是对称矩阵；

		// 1.4 确定边缘有向边的索引： 
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				if (1 == iter.value())
				{
					int repCount = adjSM_eCount.coeff(iter.row(), iter.col());
					int repOppCount = adjSM_eCount.coeff(iter.col(), iter.row());
					int bdryEdgeIdx = (repCount > 0) ? adjSM_eIdx.coeff(iter.row(), iter.col()) : adjSM_eIdx.coeff(iter.col(), iter.row());
					bdryIdxes.insert(bdryEdgeIdx);
				}
			});
	}

	// 2. 确定流形有向边（非边缘）、及其对边的索引；
	std::vector<int> edgesIdx_MN;							// 流形有向边（非边缘）的索引；
	std::vector<int> edgesIdx_MN_opp;					// 流形有向边（非边缘）的对边的索引；
	unsigned edgesCount_MN = 0;
	{
		edgesCount_MN = adjSM_MN.sum();
		edgesIdx_MN.reserve(edgesCount_MN);					// 流形有向边（非边缘）索引；
		edgesIdx_MN_opp.reserve(edgesCount_MN);			// 流形有向边（非边缘）的对边的索引；

		traverseSparseMatrix(adjSM_MN, [&](auto& iter)
			{
				edgesIdx_MN.push_back(adjSM_eIdx.coeff(iter.row(), iter.col()));
				edgesIdx_MN_opp.push_back(adjSM_eIdx.coeff(iter.col(), iter.row()));
			});
	}

	// 3. 流形边-三角片邻接关系，三角片邻接关系：
	std::vector<int> etInfo;									// 有向边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
	Eigen::VectorXi etAdj_mnEdge;						// etAdj_mnEdge(i)是索引为i的有向边所对三角片的索引，没有的话写-1；
	{
		// 3.1 生成边索引 - 三角片索引映射表etInfo;
		etInfo.resize(edgesCount);							// etInfo(i)是索引为i的边所在的三角片的索引；
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

		// 3.2 求有向边所在三角片的索引
		etAdj_mnEdge = -Eigen::VectorXi::Ones(edgesCount);
		for (unsigned i = 0; i < edgesCount_MN; ++i)
		{
			const int& edgeIdx = edgesIdx_MN[i];
			const int& edgesIdxOpp = edgesIdx_MN_opp[i];
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		// 3.3 边缘边赋值为-2；
		for (const auto& index : bdryIdxes)
			etAdj_mnEdge(index) = -2;

		//	3.3 求三角片邻接矩阵；
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}

	// 4. 计算非流形边（关联两个三角片的有向边）信息
	Eigen::MatrixXi nmnEdges;								// 非流形有向边
	std::unordered_map<std::int64_t, std::unordered_set<int>> edgeMap;							// 边编码——64位整型数数表示的边数据（两个端点索引）；
	std::unordered_map<int, std::unordered_set<int>> edgeIdx_nmn_map;						// 非流形有向边的索引 - 该边对应的所有索引；
	std::unordered_map<int, std::unordered_set<int>> edgeIdx_nmn_opp_map;				// 非流形有向边的索引 - 该边的对边对应的所有索引；
	{
		//	4.1 遍历邻接矩阵找出所有非流形有向边——若某条有向边或其对边重复次数大于1，则这两条边都被视为非流形有向边；
		std::vector<int> edgeIdx_nmn;							// 非流形有向边索引；
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (adjSM_eCount.coeff(edges(i, 0), edges(i, 1)) > 1 || adjSM_eCount.coeff(edges(i, 1), edges(i, 0)) > 1)
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();

		int neCount = edgeIdx_nmn.size();						// 非流形有向边个数；
		subFromIdxVec(nmnEdges, edges, edgeIdx_nmn);

		// 4.2 建立边编码-边索引的哈希表；
		for (int i = 0; i < edges.rows(); ++i)
		{
			std::int64_t eCode = encodeEdge(edges(i, 0), edges(i, 1));
			auto iter = edgeMap.find(eCode);
			if (edgeMap.end() == iter)
				edgeMap.insert({ eCode, std::unordered_set<int>{} });
			edgeMap[eCode].insert(i);
		}

		// 4.3 在哈希表中搜索所有非流形边，建立非流形边及其对边的（边——边索引）映射关系；
		for (int i = 0; i < neCount; ++i)
		{
			int neIdx = edgeIdx_nmn[i];		// 当前非流形边索引；
			int vaIdx = nmnEdges(i, 0);
			int vbIdx = nmnEdges(i, 1);
			std::int64_t eCode = encodeEdge(vaIdx, vbIdx);
			std::int64_t oppEcode = encodeEdge(vbIdx, vaIdx);
			auto oppIter = edgeMap.find(oppEcode);

			// 搜索当前非流形有向边的所有索引： 
			edgeIdx_nmn_map.insert({ neIdx, edgeMap[eCode] });

			// 当前边的对边有可能不存在；
			if (edgeMap.end() != oppIter)
				edgeIdx_nmn_opp_map.insert({ neIdx, oppIter->second });
		}
	}

	// 5. 求含有非流形边的三角片的邻接关系：
	{
		ttAdj_nmnEdge.resize(trisCount);
		for (auto& vec : ttAdj_nmnEdge)
			vec.resize(3);

		for (const auto& pair : edgeIdx_nmn_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;
			for (auto& index : pair.second)
				ttAdj_nmnEdge[row][col].push_back(etInfo[index]);
		}

		ttAdj_nmnOppEdge.resize(trisCount);
		for (auto& vec : ttAdj_nmnOppEdge)
			vec.resize(3);
		for (const auto& pair : edgeIdx_nmn_opp_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;
			for (auto& index : pair.second)
				ttAdj_nmnOppEdge[row][col].push_back(etInfo[index]);
		}
	}

	return true;
}


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


// 光顺非回路曲线：
template <typename DerivedV>
bool smoothCurve(Eigen::PlainObjectBase<DerivedV>& versOut, \
	const Eigen::MatrixBase<DerivedV>& versIn, const int times, \
	const double dbWeight = 0.5)
{
	/*
		除首尾两个顶点外，对于点集中的每一个点p，设前一个点为q1，后一个点为q2。
		q1q2连线的加权中点为q——q = q1+dbWeight*(q2-q1);
		平滑操作生成的新顶点p_new为p和q的中点。

 
		const VFVECTOR3& p = versOut[k];
		const VFVECTOR3& q1 = versOut[k - 1];
		const VFVECTOR3& q2 = versOut[k + 1];
		VFVECTOR3 q = q1 + dbWeight * (q2 - q1);
		versOut[k] = 0.5 * (p + q);
 
	*/ 

	using Scalar = typename DerivedV::Scalar;

	versOut.resize(0, 0);
	const int versCount = versIn.rows();
	if (versCount < 3)
		return false;

	versOut = versIn;
	for (unsigned i = 0; i < times; ++i) 
		for (unsigned k = 1; k < versCount - 1; ++k) 
			versOut.row(k) = 0.5 * (versOut.row(k) +\
				(1 - dbWeight) * versOut.row(k - 1) + dbWeight * versOut.row(k +1)); 
	return true;
} 


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


// 环路点云均匀化：
template <typename DerivedVO, typename DerivedVI>
bool arrangeLoop(Eigen::PlainObjectBase<DerivedVO>& loopOut, \
	Eigen::MatrixBase<DerivedVI>& loopIn, const int tarVersCount)
{
	const int dim = loopIn.cols();
	using ScalarO = typename DerivedVO::Scalar;
	using MatrixXO = Eigen::Matrix<ScalarO, Eigen::Dynamic, Eigen::Dynamic>;
	using RowVectorXO = Eigen::Matrix<ScalarO, 1, Eigen::Dynamic>;

	// 1.计算输出环路的期望边长；
	const int versCount = loopIn.rows();
	float length = 0;
	float eLen = 0;
	float tarEdgeLen = 0;								// 期望的精简后的环路平均边长；
	RowVectorXO arrow;
	std::vector<ScalarO> edgeLens(versCount);		
	for (int i = 0; i < versCount - 1; ++i)
	{
		arrow = (loopIn.row(i + 1) - loopIn.row(i)).cast<ScalarO>();
		eLen = arrow.norm();
		edgeLens[i] = eLen;
		length += eLen;
	}
	arrow = (loopIn.row(0) - loopIn.row(versCount - 1)).cast<ScalarO>();
	eLen = arrow.norm();
	edgeLens[versCount - 1] = eLen;
	length += eLen;
	tarEdgeLen = length / tarVersCount;

	// 2. 对长度大于0.5 * tarEdgeLen的边进行插值： 
	int index = 0;
	int stepCount = 0;
	int versCountNew = 0;
	float threshold = 0.5 * tarEdgeLen;
	float step = 0;
	MatrixXO loopNew;											// 插值均匀化以后的回路点云；
	RowVectorXO newVer;
	loopNew.resize(10 * std::max(versCount, tarVersCount), dim);
	loopNew.row(index++) = loopIn.row(0).cast<ScalarO>(); 
	for (int i = 0; i < versCount - 1; ++i)				// 输出回路生长的循环：
	{
		eLen = edgeLens[i];
		if (eLen > threshold)				// 若原边长大于阈值，则插点：
		{
			arrow = (loopIn.row(i + 1) - loopIn.row(i)).cast<ScalarO>();
			arrow.normalize();
			stepCount = std::ceil(eLen/threshold);
			step = eLen / stepCount;
			for (int k = 0; k < stepCount; ++k)
			{
				newVer = loopIn.row(i).cast<ScalarO>() + k * step * arrow;
				loopNew.row(index++) = newVer;
			}
		}
		loopNew.row(index++) = loopIn.row(i + 1).cast<ScalarO>();
	}
	eLen = edgeLens[versCount - 1];
	if (eLen > threshold)						// 若原边长大于阈值，则插点：
	{
		arrow = (loopIn.row(0) - loopIn.row(versCount - 1)).cast<ScalarO>();
		arrow.normalize();
		stepCount = std::ceil(eLen / threshold);
		step = eLen / stepCount;
		for (int k = 0; k < stepCount; ++k)
		{
			newVer = loopIn.row(versCount - 1).cast<ScalarO>() + k * step * arrow;
			loopNew.row(index++) = newVer;
		}
	}
	versCountNew = index;
	loopNew.conservativeResize(versCountNew, dim);

	// 3. 插值之后的回路精简： 
	float accumLen = 0;
	index = 0;
	loopOut.resize(versCountNew, dim);
	loopOut.row(index++) = loopNew.row(0);
	for (int i = 1; i < versCountNew; ++i)
	{
		eLen = (loopNew.row(i) - loopNew.row(i - 1)).norm();
		accumLen += eLen;
		if (accumLen >= tarEdgeLen * index)
		{
			loopOut.row(index++) = loopNew.row(i);
			if (tarVersCount == index)
				break;
		}
	}
	loopOut.conservativeResize(index, dim);

	return true;
}



/////////////////////////////////////////////////////////////////////////////////////////////////// 三角网格编辑：

// 网格串联——合并两个孤立的网格到一个网格里 
template <typename DerivedV1, typename DerivedV2>
void concatMeshMat(Eigen::PlainObjectBase<DerivedV1>& vers, \
	Eigen::MatrixXi& tris, const Eigen::MatrixBase<DerivedV2>& vers2, \
	const Eigen::MatrixXi& tris2)
{
	using TV = typename DerivedV1::Scalar;
	const int versCount = vers.rows();
	matInsertRows(vers, vers2);

	Eigen::MatrixXi trisCopy2 = tris2;
	int* intPtr = trisCopy2.data();
	for (int i = 0; i < trisCopy2.size(); ++i)
		*(intPtr++) = versCount + *intPtr;

	matInsertRows(tris, trisCopy2);
}
 

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


// 确定三角网格中的所有单连通区域（顶点连通）；重载1——使用稀疏矩阵；
template <typename DerivedI1, typename DerivedI2>
int simplyConnectedRegion(Eigen::PlainObjectBase<DerivedI1>& connectedLabels, \
	Eigen::PlainObjectBase<DerivedI2>& connectedCount, const Eigen::SparseMatrix<int>& adjSM)
{
	/*
		int simplyConnectedRegion(												返回单连通区域的数量，出错时返回-1
					Eigen::VectorXi& connectedLabels,						同一连通区域内的顶点标签，依次为0, 1, 2,...；
																											索引为i的顶点的标签为connectedLabels(i);
					Eigen::VectorXi& connectedCount							每个标签对应的单连通区域内包含的顶点数
																											标签为i的单连通区域包含的顶点数为connectedCount(i)
					const Eigen::SparseMatrix<Index>& adjSM,			三角网格的邻接矩阵
					)
	*/
	if (adjSM.rows() == 0)
		return -1;

	if (adjSM.rows() != adjSM.cols())
		return -1;

	const int versCount = adjSM.rows();

	// 1. 初始化：
	connectedLabels.setConstant(versCount, versCount);		// 最终结果中最大可能的标签为versCount-1;
	connectedCount.setZero(versCount);

	// 2. 遍历所有顶点，搜索与其联通的顶点：
	int currentLabel = 0;
	for (int i = 0; i < versCount; ++i)
	{
		// 2.1 若当前顶点已被访问过(label < versCount)，则continue:
		if (connectedLabels(i) < versCount)
			continue;

		// 2.2 若当前顶点未被访问过，执行BFS收集其所有联通的顶点：
		std::queue<int> workingQueue;
		workingQueue.push(i);
		while (!workingQueue.empty())
		{
			const int curVerIdx = workingQueue.front();
			workingQueue.pop();

			if (connectedLabels(curVerIdx) < versCount)
				continue;

			// 2.2.1 队首顶点label赋值，label计数+1
			connectedLabels(curVerIdx) = currentLabel;
			connectedCount(currentLabel)++;

			// 2.2.2 在邻接矩阵中搜索当前顶点邻接的顶点：
			for (typename Eigen::SparseMatrix<int>::InnerIterator iter(adjSM, curVerIdx); iter; ++iter)
			{
				const int connectVerIdx = iter.row();					// 默认稀疏矩阵列优先存储；当前迭代器在第curVerIdx列；

				if (connectedLabels(connectVerIdx) < versCount)
					continue;

				workingQueue.push(connectVerIdx);
			}
		}

		// 2.3 上一个标签的顶点收集完毕，下一个循环收集下一个标签的顶点：
		currentLabel++;
	}

	// 3. shrink_to_fit()
	connectedCount.conservativeResize(currentLabel, 1);

	return currentLabel;
}


// 确定三角网格中的所有单连通区域（顶点连通）；重载2——使用std::unordered_set，不使用稀疏矩阵；
template <typename DerivedV, typename	 DerivedI>
int simplyConnectedRegion(std::vector<int>& connectedLabels, std::vector<int>& connectedCount, \
	const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixBase<DerivedI>& tris)
{
	auto collectVers = [](std::vector<int>& connectedLabels, std::vector<bool>& triCollected, \
		std::unordered_set<int>& versSet, const Eigen::MatrixBase<DerivedI>& tris, const int label)->int
	{
		const int trisCount = tris.rows();
		int collectedCount = 0;
		bool blCollectedNew = false;
		for (const auto& flag : triCollected)
		{
			if (flag)
				collectedCount++;
		}

		do
		{
			blCollectedNew = false;
			for (int i = 0; i < trisCount; ++i)
			{
				if (triCollected[i])
					continue;

				int a = static_cast<int>(tris(i, 0));
				int b = static_cast<int>(tris(i, 1));
				int c = static_cast<int>(tris(i, 2));
				auto iter1 = versSet.find(a);
				auto iter2 = versSet.find(b);
				auto iter3 = versSet.find(c);
				if (versSet.end() != iter1 || versSet.end() != iter2 || versSet.end() != iter3)
				{
					versSet.insert(a);
					versSet.insert(b);
					versSet.insert(c);
					connectedLabels[a] = label;
					connectedLabels[b] = label;
					connectedLabels[c] = label;
					triCollected[i] = true;
					collectedCount++;
					blCollectedNew = true;
				}
			}
		} while (blCollectedNew);						// 本轮有新元素被收集，则继续执行一轮

		return versSet.size();										// 返回收集的顶点个数，即被标记为label的顶点个数；
	};

	if (vers.rows() == 0 || tris.rows() == 0)
		return -1;
	connectedLabels.clear();
	connectedCount.clear();

	const int versCount = vers.rows();
	const int trisCount = tris.rows();

	// 1. 区域生长的循环：
	connectedLabels.resize(versCount, versCount);		// 最终结果中最大可能的标签为versCount-1;
	std::unordered_set<int> versSet;
	std::vector<bool> triCollected(trisCount, false);
	int currentLabel = 0;
	int collectedVersCount = 0;
	while (collectedVersCount < versCount)
	{
		// w1. 选取第一个未被标记的三角片，将其三个顶点作为区域生长的种子顶点；
		for (int i = 0; i < trisCount; ++i)
		{
			if (!triCollected[i])
			{
				int a = static_cast<int>(tris(i, 0));
				int b = static_cast<int>(tris(i, 1));
				int c = static_cast<int>(tris(i, 2));
				triCollected[i] = true;
				versSet.insert(a);
				versSet.insert(b);
				versSet.insert(c);
				connectedLabels[a] = currentLabel;
				connectedLabels[b] = currentLabel;
				connectedLabels[c] = currentLabel;
				break;
			}
		}

		// w2. 区域生长的循环：
		int currentLabelVersCount = collectVers(connectedLabels, triCollected, versSet, tris, currentLabel);

		// w3. post procedure:
		versSet.clear();
		collectedVersCount += currentLabelVersCount;
		currentLabel++;
	}

	// 2. 统计：
	int scCount = currentLabel;
	connectedCount.resize(scCount, 0);
	for (int i = 0; i < versCount; ++i)
	{
		int label = connectedLabels[i];
		connectedCount[label]++;
	}

	return scCount;
}


// 确定三角网格中的所有单连通区域（三角片连通）——比triangleGrow更强，输入网格可以带有非流形元素；
template <typename DerivedI>
int simplyTrisConnectedRegion(Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount, \
	const Eigen::MatrixBase<DerivedI>& tris)
{
	/*
	int simplyTrisConnectedRegion(											返回三角片单连通区域的数量，出错时返回-1
				Eigen::VectorXi& connectedLabels,						同一连通区域内的三角片标签，依次为0, 1, 2,...；
																										索引为i的三角片的标签为connectedLabels(i);
				Eigen::VectorXi& connectedCount							每个标签对应的单连通区域内包含的三角片数
																										标签为i的单连通区域包含的三角片数为connectedCount(i)
				const Eigen::PlainObjectBase<DerivedI>& tris
				)
*/
	const unsigned trisCount = tris.rows();
	connectedLabels = Eigen::VectorXi::Zero(trisCount);
	unsigned visitedCount = 0;

	// 1. 计算网格三角片邻接关系：
	Eigen::MatrixXi ttAdj_mnEdge, edges;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	getEdges(edges, tris);
	edges2adjMat(adjSM_eCount, adjSM_eIdx, edges);
	if (!getTrisAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, edges, adjSM_eCount, adjSM_eIdx))
		return -1;

	// 2. 邻接三角片生长的循环：
	int currentLabel = 1;
	while (visitedCount < trisCount)
	{
		// w1. 选取种子：
		int seedIdx = 0;
		for (unsigned i = 0; i < trisCount; ++i)
		{
			if (0 == connectedLabels(i))
			{
				seedIdx = i;
				break;
			}
		}

		// w2. 开始生长：
		std::unordered_set<int> workingQueue;
		workingQueue.insert(seedIdx);
		while (!workingQueue.empty())
		{
			// ww1. 队首
			int cIdx = *workingQueue.begin();
			visitedCount++;
			workingQueue.erase(cIdx);
			connectedLabels(cIdx) = currentLabel;

			// ww2. 搜索当前三角片的流形有向边相对的三角片
			for (int i = 0; i < 3; ++i)
				if (ttAdj_mnEdge(cIdx, i) >= 0)
					if (0 == connectedLabels(ttAdj_mnEdge(cIdx, i)))			// 标签为0表示该三角片未被访问过；
						workingQueue.insert(ttAdj_mnEdge(cIdx, i));					// 只插入未访问的三角片 

			// ww3. 搜索当前三角片的非流形有向边相对的三角片：
			for (int i = 0; i < 3; ++i)
				for (const auto& index : ttAdj_nmnOppEdge[cIdx][i])
					if (0 == index)
						workingQueue.insert(index);
		}

		// w3. 当前联通区域生长完成
		currentLabel++;
	}

	// 3. 
	unsigned regionCount = currentLabel - 1;
	connectedLabels.array() -= 1;
	connectedCount = Eigen::VectorXi::Zero(regionCount);
	for (unsigned i = 0; i < trisCount; ++i)
		connectedCount[connectedLabels[i]]++;


	return static_cast<int>(regionCount);
}


// 提取三角网格中最大单连通区域（顶点连通）
template <typename DerivedV, typename DerivedI>
bool simplyConnectedLargest(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixBase<DerivedI>& tris)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();

	// 1. 生成邻接矩阵：
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx, adjSM;
	tris2adjMat(adjSM_eCount, adjSM_eIdx, tris);
	adjSM = adjSM_eCount;
	traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
		{
			if (iter.value() > 0)
				iter.valueRef() = 1;
		});

	// 2. 确定网格中所有单连通区域：
	Eigen::VectorXi connectedLabels, connectedCount;
	int conCount = simplyConnectedRegion(connectedLabels, connectedCount, adjSM);
	if (conCount < 0)
		return false;

	// 3. 确定最大单连通区域（顶点数最多）；
	int mainLabel = 0;								// 最大单连通区域的标签；
	int mainLabelCount = 0;
	for (int i = 0; i < conCount; ++i)
	{
		if (connectedCount(i) > mainLabelCount)
		{
			mainLabel = i;
			mainLabelCount = connectedCount(i);
		}
	}

	// 4. 提取最大单连通区域的顶点：
	std::vector<int> mainLabelIdxes;
	mainLabelIdxes.reserve(versCount);
	for (unsigned i = 0; i < versCount; ++i)
		if (mainLabel == connectedLabels(i))
			mainLabelIdxes.push_back(i);
	mainLabelIdxes.shrink_to_fit();
	subFromIdxVec(versOut, vers, mainLabelIdxes);


	// 5. 提取最大单连通区域的三角片：

	//		5.1 生成老-新索引映射表
	std::vector<int> oldNewIdxInfo(versCount, -1);
	for (int i = 0; i < mainLabelIdxes.size(); ++i)
	{
		int oldIdx = mainLabelIdxes[i];
		oldNewIdxInfo[oldIdx] = i;
	}

	//		5.2 三角片数据中的顶点索引映射成新索引；
	DerivedI trisCopy = tris;
	int* intPtr = trisCopy.data();
	for (int i = 0; i < trisCopy.size(); ++i)
	{
		int oldIdx = *intPtr;
		*intPtr = oldNewIdxInfo[oldIdx];
		intPtr++;
	}

	//		5.3 提取最大联通区域内的三角片：
	trisOut.setZero(trisCount, 3);
	int trisCountNew = 0;
	for (int i = 0; i < trisCount; ++i)
	{
		if (trisCopy(i, 0) >= 0 && trisCopy(i, 1) >= 0 && trisCopy(i, 2) >= 0)
		{
			trisOut.row(trisCountNew) = trisCopy.row(i);
			trisCountNew++;
		}
	}

	//		5.4 shrink_to_fit();
	trisOut.conservativeResize(trisCountNew, 3);


	return true;
}


// 按顶点单连通区域分解网格：
template <typename T, typename DerivedVi>
bool simplyConnectedSplitMesh(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& compVers,\
	std::vector<Eigen::MatrixXi>& compTris, const Eigen::MatrixBase<DerivedVi>& vers, \
	const Eigen::MatrixXi& tris)
{ 
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	if (0 == vers.rows() || 0 == tris.rows())
		return false;

	compVers.clear();
	compTris.clear();

	// 0. 预处理——去除孤立顶点，去除非法、重复三角片：
	MatrixXT vers0;
	Eigen::MatrixXi tris0, trisTmp;
	trisTmp = tris;
	removeSickDupTris(vers, trisTmp);
	std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, tris);
	if (!isoVerIdxes.empty())
		removeIsoVers(vers0, tris0, vers, trisTmp, isoVerIdxes);
	else
	{
		vers0 = vers.array().cast<T>();
		tris0 = trisTmp;
	}

	// 1. 提取顶点单连通点云：
	const int versCount = vers0.rows();
	int scCount = 0;							// 顶点单连通区域个数；
	int sctCount = 0;							// 三角片单连通区域个数；
	Eigen::VectorXi connectedLabels, connectedCount;

	Eigen::SparseMatrix<int> adjSM, adjSM_eCount, adjSM_eIdx;
	tris2adjMat(adjSM_eCount, adjSM_eIdx, tris0);
	adjSM = adjSM_eCount;
	traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
		{
			if (iter.value() > 0)
				iter.valueRef() = 1;
		});
	scCount = simplyConnectedRegion(connectedLabels, connectedCount, adjSM);

	// 2. 生成顶点单连通网格：
	compVers.resize(scCount);
	compTris.resize(scCount);
	for (int i = 0; i < scCount; ++i)
	{
		Eigen::VectorXi flag = (connectedLabels.array() == i).select(Eigen::VectorXi::Ones(versCount), Eigen::VectorXi::Zero(versCount));
		std::vector<int> oldNewIdxInfo, newOldIdxInfo;
		subFromFlagVec(compVers[i], oldNewIdxInfo, newOldIdxInfo, vers0, eigenVec2Vec(flag));
		compTris[i] = tris0;
		int* ptrData = compTris[i].data();
		for (int k = 0; k < compTris[i].size(); ++k)
		{
			*ptrData = oldNewIdxInfo[*ptrData];					// 老索引映射成新索引；
			ptrData++;
		}
		removeSickDupTris(compVers[i], compTris[i]);		// 去除非法三角片；
	}

	return true;
}


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


// 检测网格中的孤立顶点
template <typename DerivedV>
std::vector<unsigned> checkIsoVers(const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();

	Eigen::VectorXi verIdxVec = Eigen::VectorXi::LinSpaced(versCount, 0, versCount - 1);

	std::vector<unsigned> isoVerIdxes;
	const int* dataPtr = tris.data();
	for (unsigned i = 0; i < tris.size(); ++i)
	{
		verIdxVec(*dataPtr) = -1;
		dataPtr++;
	}
	for (unsigned i = 0; i < versCount; ++i)
		if (verIdxVec(i) >= 0)
			isoVerIdxes.push_back(i);

	return isoVerIdxes;
}


// 去除网格中的孤立顶点：
template <typename DerivedVo, typename DerivedVi>
bool removeIsoVers(Eigen::PlainObjectBase<DerivedVo>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::MatrixBase<DerivedVi>& vers, const Eigen::MatrixXi& tris, \
	const std::vector<unsigned>& isoVerIdxes)
{
	using ScalarO = typename DerivedVo::Scalar;

	versOut.resize(0, 0);
	trisOut.resize(0, 0);
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	unsigned isoVersCount = isoVerIdxes.size();
	unsigned versCountNew = versCount - isoVersCount;

	std::vector<int> oldNewIdxInfo, newOldIdxInfo, tmpIdxVec;
	tmpIdxVec.resize(versCount);
	for (int i = 0; i < versCount; ++i)
		tmpIdxVec[i] = i;
	for (const auto& index : isoVerIdxes)
		tmpIdxVec[index] = -1;

	newOldIdxInfo.reserve(versCountNew);
	for (const auto& index : tmpIdxVec)
		if (index >= 0)
			newOldIdxInfo.push_back(index);

	int tmpIdx = 0;
	oldNewIdxInfo = tmpIdxVec;
	for (auto& index : oldNewIdxInfo)
		if (index >= 0)
			index = tmpIdx++;

	// 输出点云：
	subFromIdxVec(versOut, vers, newOldIdxInfo);

	// 输出三角片：
	trisOut = tris;
	int* dataPtr = trisOut.data();
	for (int i = 0; i < trisCount * 3; ++i)
		*dataPtr++ = oldNewIdxInfo[*dataPtr];

	return true;
}


template <typename DerivedV>
Eigen::VectorXi checkDegEdges(const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const double threshold);


template <typename DerivedV>
int mergeDegEdges(Eigen::PlainObjectBase<DerivedV>& newVers, Eigen::MatrixXi& newTris, \
	const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const Eigen::VectorXi& degEdgeFlags);


// correctTriDirs()——矫正网格的三角片朝向，基于三角片区域生长的方法；使用前需要先处理网格中的重叠三角片；
template <typename DerivedV, typename DerivedI>
int correctTriDirs(Eigen::MatrixXi& trisOut, const Eigen::MatrixBase<DerivedV>& vers, \
	const Eigen::MatrixBase<DerivedI>& tris, const double thetaThreshold = 0.99 * pi)
{
	int corrCount = 0;
	int visitedCount = 0;
	const double dotThreshold = cos(thetaThreshold);
	const int versCount = vers.rows();
	const int trisCount = tris.rows();
	const int edgesCount = 3 * trisCount;
	trisOut = tris.array().cast<int>(); 

	// 0. 求初始状态下网格所有三角片的法向：
	Eigen::MatrixXd triNorms;
	if (!getTriNorms(triNorms, vers, tris))				// 若含有退化三角片，会导致计算面片法向出错；
		return false;

	// 1. 求三角片邻接关系：
	Eigen::MatrixXi ttAdj_mnEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return false;

	// 2. 循环——队列中第一个三角片t出队 - t写入输出集合 - 求出t所有相邻三角片的索引然后入队
	int triIdx = 0;
	std::unordered_set<int> workingSet;								// 用于寻找相邻三角片的队列；
	std::vector<int> triTags(trisCount, 0);								// 0 : 未访问，1: 已访问，-1: 已访问，被标记为朝向错误三角片；
	while (visitedCount < trisCount)
	{
		for (int i = 0; i < trisCount; ++i)
		{
			if (0 == triTags[i])
			{
				triIdx = i;
				break;
			}
		}

		workingSet.insert(static_cast<int>(triIdx));						// 队列插入第一个三角片； 
		while (workingSet.size() > 0)
		{
			// w1. 首个元素出队，插入到联通三角片集合中；
			int cTriIdx = *workingSet.begin();									// 当前处理的三角片索引
			workingSet.erase(workingSet.begin());
			if (triTags[cTriIdx] != 0)
				continue;
			triTags[cTriIdx] = 1;													// 当前三角片标记为已访问；
			visitedCount++;

			std::vector<int> adjTriIdxes;									// 当前三角片的未被访问的相邻三角片的索引
			Eigen::RowVector3d cTriNorm = triNorms.row(cTriIdx);
			Eigen::MatrixXd adjTriNorms;
			unsigned adjTrisCount = 0;

			// w2. 确定当前三角片的所有相邻三角片
			for (int i = 0; i < 3; ++i)
			{
				int nbrTriIdx = ttAdj_mnEdge(cTriIdx, i);
				if (nbrTriIdx >= 0)
					if (0 == triTags[nbrTriIdx])								// a. 当前边为流形边
						adjTriIdxes.push_back(nbrTriIdx);
					else
					{
						// b. 当前边为非流形边；
						auto& vec1 = ttAdj_nmnEdge[cTriIdx][i];						// 当前非流形边所在的所有三角片；
						auto& vec2 = ttAdj_nmnOppEdge[cTriIdx][i];				// 当前非流形边的对边所在的所有三角片；
						for (const auto& index : vec1)
						{
							if (index == cTriIdx)
								continue;
							else
								if (0 == triTags[index])
									adjTriIdxes.push_back(index);
						}

						for (const auto& index : vec2)
							if (0 == triTags[index])
								adjTriIdxes.push_back(index);
					}
			}

			// 4.3 检验-矫正当前三角片的相邻三角片：
			adjTrisCount = adjTriIdxes.size();
			if (0 == adjTrisCount)
				continue;
			subFromIdxVec(adjTriNorms, triNorms, adjTriIdxes);
			for (unsigned k = 0; k < adjTrisCount; ++k)
			{
				int adjIdx = adjTriIdxes[k];
				if (cTriNorm.dot(adjTriNorms.row(k)) < dotThreshold)
				{
					// 三角片标记为-1：
					triTags[adjIdx] = -1;
					visitedCount++;
					corrCount++;												// 计数增加； 
				}
			}

			// 4.3 工作队列中插入当前三角片未被访问的相邻三角片；
			for (const auto& index : adjTriIdxes)
				if (0 == index)
					workingSet.insert(index);
		}
	} 

	return corrCount;
}


// 去除非法三角片，重复三角片：
template<typename DerivedV>
int removeSickDupTris(const Eigen::MatrixBase<DerivedV>& vers, Eigen::MatrixXi& tris)
{
	using T = typename DerivedV::Scalar;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	int versCount = vers.rows();
	int trisCount = tris.rows();

	Eigen::VectorXi flags{ Eigen::VectorXi::Ones(trisCount) };			// flag向量，需要去除的三角片标0；
	std::unordered_set<std::uint64_t> triCodeSet;

	// 1. 遍历三角片，找出非法、重复三角片索引，在flag向量中标记为0；
	for (int i = 0; i < trisCount; ++i)
	{
		// 三角片中的索引值不在0~versCount-1之间则为非法三角片；貌似有些STL文件中会出现这种情形；
		if (tris(i, 0) < 0 || tris(i, 0) > versCount - 1 || tris(i, 1) < 0 || tris(i, 1) > versCount - 1 || tris(i, 2) < 0 || tris(i, 2) > versCount - 1)
		{
			flags(i) = 0;
			continue;
		}

		std::uint64_t code = encodeTriangle(tris(i, 0), tris(i, 1), tris(i, 2));
		if (0 == code)
		{
			flags(i) = 0;				// 三角片的三个索引中有相同的，非法；
			continue;
		}

		auto retPair = triCodeSet.insert(code);
		if (!retPair.second)
			flags(i) = 0;				// 重复三角片非法；
	}

	// 2. 根据flag向量提取新的三角片数据替代原有的;
	Eigen::MatrixXi tmpMat;
	subFromFlagVec(tmpMat, tris, eigenVec2Vec(flags));
	tris = tmpMat;

	// 3. 计算去除的三角片数；
	int removeCount = trisCount - tris.rows();
	return removeCount;
}


template <typename DerivedV, typename DerivedI>
int findHoles(std::vector<std::vector<int>>& holes, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);


template <typename DerivedI>
bool checkSickTris(std::vector<int>& sickIdxes, const Eigen::PlainObjectBase<DerivedI>& tris);



// 暂时未整理的实现：
#include "temp.tpp"