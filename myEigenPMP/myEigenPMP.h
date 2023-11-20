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


// Ŀ¼��
/*
	1. ����ת��
	2. ͼ������
	3. ���Ʊ༭
	4. ��������༭
	5. ���������㷨
	6. ����ȱ�ݼ����޸�

*/
 

/////////////////////////////////////////////////////////////////////////////////////////////////// ����ת���ӿڣ�

// std::vector<std::pair<int, int>>��ʾ�ı�����ת��Ϊ�����ʾ��
template <typename DerivedI>
void edges2mat(Eigen::PlainObjectBase<DerivedI>& mat, \
	const std::vector<std::pair<int, int>>& edges);

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


/////////////////////////////////////////////////////////////////////////////////////////////////// ���ɻ���ͼԪ��

// �����ֵ�õ�ֱ�߶ε���
template <typename DerivedVo, typename ScalarVi>
bool interpolateToLine(Eigen::PlainObjectBase<DerivedVo>& vers, \
	const Eigen::Matrix<ScalarVi, 1, 3>& start, const Eigen::Matrix<ScalarVi, 1, 3>& end, \
	const float SR, const bool SE = true);


// ����XOYƽ���ڵ�ԲȦ�㼯��
template <typename T>
bool getCircleVers(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, \
	const unsigned versCount = 20);


// �Ի�·�㼯���в���������ʷ֡����ֶ���Ϸ����� 
bool circuitGetTris(Eigen::MatrixXi& tris, const std::vector<int>& indexes, const bool regularTri);



/////////////////////////////////////////////////////////////////////////////////////////////////// ͼ�����ԣ�

// ������������Ƭ���ݣ��õ���������ݣ�
template <typename DerivedI>
bool getEdges(Eigen::MatrixXi& edges, const Eigen::MatrixBase<DerivedI>& tris);
 

// ��������·�ı����ݣ�
template <typename DerivedI, typename IndexType>
bool getSortedLoopEdges(Eigen::PlainObjectBase<DerivedI>& edges, const IndexType versCount);


// ���ɱ�-�߱����ӳ���
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

// ��������������Ƭ�����ǵ�����ֵ��
template <typename DerivedC, typename DerivedV>
void trisCotValues(Eigen::PlainObjectBase<DerivedC>& cotValues, \
	const Eigen::PlainObjectBase<DerivedV>& vers,
	const Eigen::MatrixXi& tris);

template<typename T, typename DerivedV, typename DerivedI>
bool trianglesBarycenter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& barys, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);


// ����������������Ƭ�Ĺ�һ��������
template<typename DerivedN, typename DerivedV, typename DerivedI>
bool trianglesNorm(Eigen::PlainObjectBase<DerivedN>& triNorms, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);


template<typename DerivedV, typename DerivedI>
bool trianglesPlane(Eigen::MatrixXd& planeCoeff, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);


// ������������ÿ������Ƭ�������
template<typename DerivedA, typename DerivedV, typename DerivedI>
bool trisArea(Eigen::PlainObjectBase<DerivedA>& trisAreaVec, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::MatrixBase<DerivedI>& tris);


// ����ˮ�ܵ���������������
template<typename DerivedV, typename DerivedI>
double meshVolume(const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);


// ������������ݣ��õ������ڽӾ���
/*
	bool edges2adjMat(
			Eigen::SparseMatrix<int>& adjSM_eCount,			Ȩ��Ϊ��������ظ��������ڽӾ���
			Eigen::SparseMatrix<int>& adjSM_eIdx,				Ȩ��Ϊ��������������ڽӾ���
			const Eigen::MatrixBase<DerivedI>& edges			����߾���
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
	adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());									// Ȩ��Ϊ��������ظ��Ĵ�����
	adjSM_eIdx.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());		// Ȩ��Ϊ������ߵ�������


	return true;
}


// ����������Ĳ�ͬȨ�ص��ڽӾ��� 
/*
	template<typename int>
	bool tris2adjMat(
		Eigen::SparseMatrix<int>& adjSM_eCount,			Ȩ��Ϊ����Ŀ���ڽӾ���
		Eigen::SparseMatrix<int>& adjSM_eIdx,				Ȩ��Ϊ���������ڽӾ���
		const Eigen::MatrixXi& tris,											�������������Ƭ		
		const Eigen::MatrixXi& edges0 = Eigen::MatrixXi{}			��������ı���Ϣ��Ĭ�������´���Ϊ�գ��ڴ˺����м�������ݣ�
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


// ��Ե����ߣ�����1��
template<typename DerivedI>
bool bdryEdges(Eigen::MatrixXi& bdrys, const Eigen::MatrixBase<DerivedI>& tris)
{
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	tris2adjMat(adjSM_eCount, adjSM_eIdx, tris);
	Eigen::SparseMatrix<int> adjSM_tmp, adjSM_ueCount;

	// 1. ȷ����Ե����ߣ�
	spMatTranspose(adjSM_tmp, adjSM_eCount);
	adjSM_ueCount = adjSM_eCount + adjSM_tmp;				// ������ڽӾ���Ȩ���Ǹñߵ��ظ�������

	std::vector<std::pair<int, int>> bdryUeVec;						// ���б�Ե����ߣ�Լ������߱�ʾΪǰ���С�������ԣ�
	traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
		{
			if (1 == iter.value() && iter.row() > iter.col())
				bdryUeVec.push_back({ iter.row(), iter.col() });
		});

	// 2. �ɱ�Ե�����ȷ����Ե����ߣ�
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


// ��Ե����ߣ� ����2��
template<typename DerivedI>
bool bdryEdges(Eigen::MatrixXi& bdrys, std::vector<int>& bdryTriIdxes, \
	const Eigen::MatrixBase<DerivedI>& tris)
{
	/*
		bool bdryEdges(
				Eigen::PlainObjectBase<DerivedI>& bdrys,					bdrysCount * 2����Ե�����ݣ�������ߣ�
				std::vector<int>& bdryTriIdxes,										������Ե�ߵ�����Ƭ������
				const Eigen::PlainObjectBase<DerivedI>& tris				trisCount * 3, ��������Ƭ���ݣ�
				)

	*/

	// 1. ��������ı�Ե����ߣ�
	if (!bdryEdges(bdrys, tris))
		return false;

	if (0 == bdrys.rows())
		return true;

	const unsigned trisCount = tris.rows();

	// 2. ȷ�ϱ�Ե�����α����ڵ�����Ƭ������
	Eigen::MatrixXi edgeAs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi edgeBs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi edgeCs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi vaIdxes = tris.col(0);
	Eigen::MatrixXi vbIdxes = tris.col(1);
	Eigen::MatrixXi vcIdxes = tris.col(2);

	// 3. ��������߱���-����Ƭ�ֵ䣺
	std::unordered_multimap<std::int64_t, unsigned> edgesMap;			// һ����������������ֻ����һ������Ƭ����������������Ƿ����αߣ�
	for (unsigned i = 0; i < trisCount; ++i)
	{
		std::int64_t codeA = encodeEdge(vbIdxes(i), vcIdxes(i));
		std::int64_t codeB = encodeEdge(vcIdxes(i), vaIdxes(i));
		std::int64_t codeC = encodeEdge(vaIdxes(i), vbIdxes(i));
		edgesMap.insert({ codeA, i });
		edgesMap.insert({ codeB, i });
		edgesMap.insert({ codeC, i });
	}

	// 4. ���ֵ��в������б�Ե�����ڵ�����Ƭ������
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

// ��������Ȩ��laplace����
template<typename Tl, typename DerivedV>
bool cotLaplacian(Eigen::SparseMatrix<Tl>& L, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::MatrixXi& tris);

// �������������ж����1ring��Ϣ����1ring�Ķ�������������Ƭ������
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

/////////////////////////////////////////////////////////////////////////////////////////////////// ���Ʊ༭��

// ����laplacian�Ļ�·��˳
template<typename T>
bool smoothCircuit2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& circuit, \
	const float param);


// XOYƽ���ڵĻ�·��������ɴ�Z�����򿴶������·��ʱ���˳ʱ������
/*
	bool sortLoop2D(
			Eigen::PlainObjectBase<DerivedV>& loop,			������Ļ�·
			const bool blCounterClockWise,								�Ƿ��������ʱ��
			const int startIdx														ָ�������·������ΪstartIdx�Ķ�����Ϊ������·����㣻
			)

*/
template<typename DerivedV>
bool sortLoop2D(Eigen::PlainObjectBase<DerivedV>& loop, \
	const bool blCounterClockWise = true, 	const float thresholdDis = 10);


/////////////////////////////////////////////////////////////////////////////////////////////////// ��������༭��

template <typename T>
void concatMeshMat(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers,\
	Eigen::MatrixXi& tris, 	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers1, \
	const Eigen::MatrixXi& tris1);


template <typename IndexT>
bool removeTris(Eigen::MatrixXi& trisOut, const Eigen::MatrixXi& tris, \
	const std::vector<IndexT>& sickTriIdxes);


// laplace��˳ 
/*
	bool laplaceFaring(
		Eigen::PlainObjectBase<DerivedVo>& versOut,					������񶥵�
		const Eigen::PlainObjectBase<DerivedVi>& vers,				�������񶥵�
		const Eigen::MatrixXi& tris,												����������Ƭ
		const float deltaLB,															��ɢ�ٶȲ���
		const unsigned loopCount,												��˳����
		const std::vector<int>& fixedVerIdxes							�޶����ֲ���Ķ��������
		)
*/
template <typename DerivedVo, typename DerivedVi>
bool laplaceFaring(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	const Eigen::PlainObjectBase<DerivedVi>& vers, \
	const Eigen::MatrixXi& tris, const float deltaLB, const unsigned loopCount, \
	const std::vector<int>& fixedVerIdxes = std::vector<int>{});



/////////////////////////////////////////////////////////////////////////////////////////////////// ���������㷨��

template <typename DerivedV, typename DerivedI>
bool triangleGrow(Eigen::PlainObjectBase<DerivedV>& versOut,\
	Eigen::PlainObjectBase<DerivedI>& trisOut, 	\
	const Eigen::MatrixBase<DerivedV>& vers, \
	const Eigen::MatrixBase<DerivedI>& tris, const int triIdx);


template <typename DerivedV>
bool triangleGrowSplitMesh(std::vector<DerivedV>& meshesVersOut, \
	std::vector<Eigen::MatrixXi>& meshesTrisOut, \
	const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris);


// ����Ƭ���������õ����������
/*
	template <typename DerivedVo, typename DerivedVi>
	bool triangleGrowOuterSurf(
		Eigen::PlainObjectBase<DerivedVo>& versOut, \		�������
		Eigen::MatrixXi& trisOut, 
		const Eigen::MatrixBase<DerivedVi>& vers, \			��������
		const Eigen::MatrixXi& tris,
		const bool blRemIso,													�Ƿ�ȥ����������
		const int startIdx															����������ʼ������Ƭ������
		)

*/
template <typename DerivedVo, typename DerivedVi>
bool triangleGrowOuterSurf(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::MatrixXi& trisOut, const Eigen::MatrixBase<DerivedVi>& vers,\
	const Eigen::MatrixXi& tris, const bool blRemIso = true, const int startIdx = -1);


// robustTriangleGrowOuterSurf()����������³���Ե�����Ƭ������������ĳЩ��������Ƭ�ĸ��ţ� 
/*
	template <typename DerivedVo, typename DerivedVi>
	bool robustTriangleGrowOuterSurf(
		Eigen::PlainObjectBase<DerivedVo>& versOut, \		�������
		Eigen::MatrixXi& trisOut,
		const Eigen::MatrixBase<DerivedVi>& vers, \			��������
		const Eigen::MatrixXi& tris,
		const bool blRemIso,													�Ƿ�ȥ����������
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



/////////////////////////////////////////////////////////////////////////////////////////////////// ����ȱ�ݼ����޸�
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


// ȥ�������еĹ������㣺
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





// ��ʱδ�����ʵ�֣�
#include "temp.tpp"