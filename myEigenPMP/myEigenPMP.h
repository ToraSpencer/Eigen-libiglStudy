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
	2. ���ɻ���ͼԪ
	3. ͼ������
	4. ���Ʊ༭
	5. ��������༭
	6. ���������㷨
	7. ����ȱ�ݼ����޸�

*/
 

/////////////////////////////////////////////////////////////////////////////////////////////////// ǰ��������
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
	for (unsigned i = 1; (length - lenth0) > 0.8 * SR; i++)			// ȷ�����һ��������յ㲻��̫����
	{
		Eigen::RowVector3d temp = startD + SR * i * dir;
		matInsertRows(vers, temp);
		lenth0 = SR * (i + 1);		// ��һ��temp�ĳ��ȡ�
	}

	if (SE)
		matInsertRows(vers, endD);

	return true;
};


// ����XOYƽ���ڵ�ԲȦ�㼯��
template <typename T>
bool getCircleVers(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, \
	const unsigned versCount = 20);


// �Ի�·�㼯���в���������ʷ֡����ֶ���Ϸ����� 
bool circuitGetTris(Eigen::MatrixXi& tris, const std::vector<int>& indexes, const bool regularTri);



/////////////////////////////////////////////////////////////////////////////////////////////////// ͼ�����ԣ�

// �õ�������������������
template <typename DerivedI>
bool getEdges(Eigen::MatrixXi& edges, const Eigen::MatrixBase<DerivedI>& tris)
{
	/*
		bool getEdges(
				Eigen::MatrixXi& edges,
				const Eigen::PlainObjectBase<DerivedI>& tris
				)

		��������edges�е�����˳��Ϊ[bc; ca; ab]�������ԵĶ�����corner�ֱ�Ϊ0, 1, 2����a,b,c;
		������������Ƭ������ӳ�䡪����eIdx0���ڵ�����ƬtriIdx0 == eIdx0 % trisCount;
		����ƬtriIdx0��corner0����������ӳ�䡪��eIdx0 = corner0 * trisCount + triIdx0;
		��������corner��ӳ�䡪��corner0 = eIdx0 / trisCount;
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
bool trianglesBarycenter(\
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& barys, \
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


// buildAdjacency()�����������������Ƭ�ڽ���Ϣ��
using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool buildAdjacency(Eigen::MatrixXi& ttAdj_mnEdge, std::vector<tVec>& ttAdj_nmnEdge, \
	std::vector<tVec>& ttAdj_nmnOppEdge, const Eigen::MatrixBase<DerivedI>& tris)
{
	/*
		bool buildAdjacency(
					const Eigen::MatrixXi& tris,									���������Ƭ����

					Eigen::MatrixXi& ttAdj_mnEdge,							����Ƭ�ķǱ�Ե����������ڽӵ�����Ƭ������
																									trisCount * 3
																									(i, 0)Ϊ����Ϊi������Ƭ�ķǱ�Ե���������(vb, vc)�ڽӵ�����Ƭ��������

					std::vector<ttTuple>& ttAdj_nmnEdge,				����Ƭ�ķ�������������ڵ�����Ƭ������
																									size == trisCount;
																									ÿ��ttTupleԪ�ذ�������int������
																									����Ϊi��Ԫ�صĵ�0��������Ϊ����Ϊi������Ƭ�ķ����������(vb, vc)���ڵ���������Ƭ������


					std::vector<ttTuple>& ttAdj_nmnOppEdge		����Ƭ�ķ�����������ڽӵ�����Ƭ���������Ա����ڵ�����Ƭ����������
																									size == trisCount;
																									����Ϊi��Ԫ�صĵ�0��������Ϊ����Ϊi������Ƭ�ķ����������(vb, vc)�ڽӵ���������Ƭ������
					)

	*/

	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	Eigen::MatrixXi edges;							// ��������ݣ�
	Eigen::MatrixXi nmnEdges;					// �����������

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// ��������ϡ�����adjSM_eCount, adjSM_eIdx��TRIPLET���ݣ�
	std::unordered_multimap<std::int64_t, int> edgeMap;					// �߱��롪��64λ����������ʾ�ı����ݣ������˵���������


	// 1. ������ı���Ϣ���ڽӾ���
	Eigen::SparseMatrix<int> adjSM;								// �ڽӾ���
	Eigen::SparseMatrix<int> adjSM_eCount;					// �ڽӾ�������Ϊ��������ظ��Ĵ�����
	Eigen::SparseMatrix<int> adjSM_eIdx;						// �ڽӾ���Ԫ������Ӧ���α���Ϊ�ñ�����������Ӧ�����α���Ϊ�ñ����������ĺͣ�
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;				// adjSM_eIdx��ת�ã���ʾ�Աߵ���Ϣ��
	Eigen::SparseMatrix<int> adjSM_ueCount;				// Ȩ��Ϊ�����ij����������Ƭ������
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;
	{

		// edges == [ea; eb; ec] == [vbIdxes, vcIdxes; vcIdxes, vaIdxes; vaIdxes, vbIdxes];
		/*
			������������������edges�б���ÿ���߶���unique�ģ�
			�����ڷ����αߣ�������α���edges����ظ��洢��ʵ�ʱ�����Ҳ������edges������
			����˵ʹ������representation�Ļ��������α��в�ֹһ��������ȡ�������÷����αߵ�����Ƭ������
		*/

		// ����Ƭ�����ߵ�������teIdx == [eaIdx, ebIdx, ecIdx] == [(0: trisCount-1)', (trisCount: 2*trisCount-1)', (2*trisCount, 3*trisCount-1)'];
		/*
			teIdx(i, j) = trisCount *j + i;
		*/

		// 1.1 �������������
		getEdges(edges, tris);

		// 1.2 ��������������ڽӾ���adjSM, adjSM_eCount, adjSM_eCount;
		edges2adjMat(adjSM_eCount, adjSM_eIdx, edges);
		adjSM = adjSM_eCount;																		// ������ڽӾ���
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		��adjSM_eIdx(i, j)��Ӧ������������ߣ���Ȩ��ֵΪ�ñߵ����������Ƿ���������ߣ���ñ߶�Ӧ������������Ȩ��Ϊ��������֮�ͣ�
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

		//	1.3 ����������ڽӾ���
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// ������ڽӾ���
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.4 ���ɷǱ�Ե����������ڽӾ���
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

		//	1.5 ���ɷǱ�Ե��������ߡ�����Աߵ��ڽӾ���
		adjSM_MN_NB = adjSM + adjSM_ueMN_NB;					// adjSM & adjSM_ueMN_NB
		adjSM_MN_NB.prune([](const Eigen::Index& row, const Eigen::Index& col, const int& value)->bool
			{
				if (2 == value)			// ɾ�������α�
					return true;
				else
					return false;
			});
		adjSM_MN_NB /= 2;
		spMatTranspose(adjSM_MN_NB_opp, adjSM_MN_NB);
	}


	// 2. ȷ���Ǳ�Ե��������ߡ�����Աߵ�������
	std::vector<int> edgesIdx_MN_NB;							// �Ǳ�Ե��������ߵ�������
	std::vector<int> edgesIdx_MN_NB_opp;					// �Ǳ�Ե��������ߵĶԱߵ�������
	{
		unsigned edgesCount_MN_NB = adjSM_MN_NB.sum();
		edgesIdx_MN_NB.reserve(edgesCount_MN_NB);					// �Ǳ�Ե���������������
		edgesIdx_MN_NB_opp.reserve(edgesCount_MN_NB);			// �Ǳ�Ե��������ߵĶԱߵ�������

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB.push_back(adjSM_eIdx.coeffRef(iter.row(), iter.col()));
			});

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB_opp.push_back(adjSM_eIdx_opp.coeffRef(iter.row(), iter.col()));
			});
	}


	// 3. ���α�-����Ƭ�ڽӹ�ϵ������Ƭ�ڽӹ�ϵ��
	std::vector<int> etInfo;						// ������ - ����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
	Eigen::VectorXi etAdj_mnEdge;							// ����������Ƭ�������������α߻��Ե��д-1��
	{
		// 3.1 ���ɱ����� - ����Ƭ����ӳ���etInfo;
		etInfo.resize(edgesCount);							// etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

		// 3.2 �����������Ƭ��������
		etAdj_mnEdge = -Eigen::VectorXi::Ones(edgesCount);					// ����������Ƭ�������������α߻��Ե��д-1��
		for (unsigned i = 0; i < edgesIdx_MN_NB.size(); ++i)
		{
			const int& edgeIdx = edgesIdx_MN_NB[i];
			const int& edgesIdxOpp = edgesIdx_MN_NB_opp[i];

			if (edgesIdxOpp < 0 || edgesIdxOpp >= edgesCount)
				std::cout << "pause" << std::endl;

			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//	3.3 ������Ƭ�ڽӾ���ttAdj_mnEdge(i, :)������Ϊi������Ƭ�������ڽӵ���������Ƭ���������з����α߻��Ե����дΪ-1��
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. ��������αߣ�������������Ƭ������ߣ���Ϣ
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_map;				// ����������ߵ����� - �ñ߶�Ӧ������������
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_opp_map;		// ����������ߵ����� - �ñߵĶԱ߶�Ӧ������������
	{
		//	4.1 �����ڽӾ����ҳ����з����αߣ�������������Ƭ������ߣ���
		std::vector<int> edgeIdx_nmn;						// �����������������
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)) > 1)
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();

		int neCount = edgeIdx_nmn.size();						// ����������߸�����
		subFromIdxVec(nmnEdges, edges, edgeIdx_nmn);

		// 4.2 �����߱���-�������Ĺ�ϣ��
		for (int i = 0; i < edges.rows(); ++i)
		{
			std::int64_t eCode = encodeEdge(edges(i, 0), edges(i, 1));
			edgeMap.insert({ eCode, i });
		}

		// 4.3 �ڹ�ϣ�����������з����αߣ����������α߼���Աߵģ��ߡ�����������ӳ���ϵ��
		for (int i = 0; i < neCount; ++i)
		{
			// ע�⣡������������������� a. ��ǰ���Ƿ����αߣ��Ա������αߣ�	b. ��ǰ���Ƿ����αߣ��Ա߲����ڣ�
			int neIdx = edgeIdx_nmn[i];		// ��ǰ�����α�������
			int vaIdx = nmnEdges(i, 0);
			int vbIdx = nmnEdges(i, 1);
			std::int64_t eCode = encodeEdge(vaIdx, vbIdx);
			std::int64_t oppEcode = encodeEdge(vbIdx, vaIdx);
			auto iter = edgeMap.find(eCode);
			auto oppIter = edgeMap.find(oppEcode);

			int otherIdx = (iter->second == neIdx) ? ((++iter)->second) : (iter->second);
			edgeIdx_nmn_map.insert({ neIdx, std::vector<int>{neIdx, otherIdx} });

			if (edgeMap.end() == oppIter)			// ��ǰ���Ƿ����αߣ��Ա߲����ڣ�
			{
				edgeIdx_nmn_opp_map.insert({ neIdx, std::vector<int>{} });
				continue;
			}
			int oppIdx1 = oppIter->second;
			int oppIdx2 = (++oppIter)->second;
			edgeIdx_nmn_opp_map.insert({ neIdx, std::vector<int>{oppIdx1, oppIdx2} });
		}
	}


	// 5. ���з����αߵ�����Ƭ���ڽӹ�ϵ��
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


// getTrisAdjacency()�����������������Ƭ�ڽ���Ϣ��
using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool getTrisAdjacency(Eigen::MatrixXi& ttAdj_mnEdge, std::vector<tVec>& ttAdj_nmnEdge, \
	std::vector<tVec>& ttAdj_nmnOppEdge, const Eigen::MatrixBase<DerivedI>& edges, \
	const Eigen::SparseMatrix<int>& adjSM_eCount, 	const Eigen::SparseMatrix<int>& adjSM_eIdx)
{ 
	const int edgesCount = edges.rows();
	const int trisCount = edgesCount / 3;

	// 1. ������ı���Ϣ���ڽӾ���  
	Eigen::SparseMatrix<int> adjSM_eCount_opp; 
	Eigen::SparseMatrix<int> adjSM_ueCount;				// Ȩ��Ϊ�����ij����������Ƭ������ 
	Eigen::SparseMatrix<int> adjSM_MN;						// ��������ߣ��Ǳ�Ե���ڽӾ���Ȩ��Ϊ1�� 
	std::unordered_set<int> bdryIdxes;									// ��Ե����ߵ�������
	{ 
		spMatTranspose(adjSM_eCount_opp, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + adjSM_eCount_opp;						// ������ڽӾ���Ȩ��Ϊ�����ij����������Ƭ������ 						 

		//	1.3 ������������ߣ��Ǳ�Ե�����ڽӾ���
		adjSM_MN = adjSM_ueCount;
		traverseSparseMatrix(adjSM_MN, [&adjSM_MN](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});
		adjSM_MN.prune([&adjSM_eCount](const Eigen::Index& row, const Eigen::Index& col, const int& value)->bool
			{
				if (1 == adjSM_eCount.coeff(row, col) && 1 == adjSM_eCount.coeff(col, row))
					return true;								// ����true��Ԫ�ر�������
				else
					return false;
			});															// �ǶԳƾ���

		// 1.4 ȷ����Ե����ߵ������� 
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

	// 2. ȷ����������ߣ��Ǳ�Ե��������Աߵ�������
	std::vector<int> edgesIdx_MN;							// ��������ߣ��Ǳ�Ե����������
	std::vector<int> edgesIdx_MN_opp;					// ��������ߣ��Ǳ�Ե���ĶԱߵ�������
	unsigned edgesCount_MN = 0;
	{
		edgesCount_MN = adjSM_MN.sum();
		edgesIdx_MN.reserve(edgesCount_MN);					// ��������ߣ��Ǳ�Ե��������
		edgesIdx_MN_opp.reserve(edgesCount_MN);			// ��������ߣ��Ǳ�Ե���ĶԱߵ�������

		traverseSparseMatrix(adjSM_MN, [&](auto& iter)
			{
				edgesIdx_MN.push_back(adjSM_eIdx.coeff(iter.row(), iter.col()));
				edgesIdx_MN_opp.push_back(adjSM_eIdx.coeff(iter.col(), iter.row()));
			});
	}

	// 3. ���α�-����Ƭ�ڽӹ�ϵ������Ƭ�ڽӹ�ϵ��
	std::vector<int> etInfo;									// ��������� - ����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
	Eigen::VectorXi etAdj_mnEdge;						// etAdj_mnEdge(i)������Ϊi���������������Ƭ��������û�еĻ�д-1��
	{
		// 3.1 ���ɱ����� - ����Ƭ����ӳ���etInfo;
		etInfo.resize(edgesCount);							// etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

		// 3.2 ���������������Ƭ������
		etAdj_mnEdge = -Eigen::VectorXi::Ones(edgesCount);
		for (unsigned i = 0; i < edgesCount_MN; ++i)
		{
			const int& edgeIdx = edgesIdx_MN[i];
			const int& edgesIdxOpp = edgesIdx_MN_opp[i];
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		// 3.3 ��Ե�߸�ֵΪ-2��
		for (const auto& index : bdryIdxes)
			etAdj_mnEdge(index) = -2;

		//	3.3 ������Ƭ�ڽӾ���
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}

	// 4. ��������αߣ�������������Ƭ������ߣ���Ϣ
	Eigen::MatrixXi nmnEdges;								// �����������
	std::unordered_map<std::int64_t, std::unordered_set<int>> edgeMap;							// �߱��롪��64λ����������ʾ�ı����ݣ������˵���������
	std::unordered_map<int, std::unordered_set<int>> edgeIdx_nmn_map;						// ����������ߵ����� - �ñ߶�Ӧ������������
	std::unordered_map<int, std::unordered_set<int>> edgeIdx_nmn_opp_map;				// ����������ߵ����� - �ñߵĶԱ߶�Ӧ������������
	{
		//	4.1 �����ڽӾ����ҳ����з���������ߡ�����ĳ������߻���Ա��ظ���������1�����������߶�����Ϊ����������ߣ�
		std::vector<int> edgeIdx_nmn;							// �����������������
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (adjSM_eCount.coeff(edges(i, 0), edges(i, 1)) > 1 || adjSM_eCount.coeff(edges(i, 1), edges(i, 0)) > 1)
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();

		int neCount = edgeIdx_nmn.size();						// ����������߸�����
		subFromIdxVec(nmnEdges, edges, edgeIdx_nmn);

		// 4.2 �����߱���-�������Ĺ�ϣ��
		for (int i = 0; i < edges.rows(); ++i)
		{
			std::int64_t eCode = encodeEdge(edges(i, 0), edges(i, 1));
			auto iter = edgeMap.find(eCode);
			if (edgeMap.end() == iter)
				edgeMap.insert({ eCode, std::unordered_set<int>{} });
			edgeMap[eCode].insert(i);
		}

		// 4.3 �ڹ�ϣ�����������з����αߣ����������α߼���Աߵģ��ߡ�����������ӳ���ϵ��
		for (int i = 0; i < neCount; ++i)
		{
			int neIdx = edgeIdx_nmn[i];		// ��ǰ�����α�������
			int vaIdx = nmnEdges(i, 0);
			int vbIdx = nmnEdges(i, 1);
			std::int64_t eCode = encodeEdge(vaIdx, vbIdx);
			std::int64_t oppEcode = encodeEdge(vbIdx, vaIdx);
			auto oppIter = edgeMap.find(oppEcode);

			// ������ǰ����������ߵ����������� 
			edgeIdx_nmn_map.insert({ neIdx, edgeMap[eCode] });

			// ��ǰ�ߵĶԱ��п��ܲ����ڣ�
			if (edgeMap.end() != oppIter)
				edgeIdx_nmn_opp_map.insert({ neIdx, oppIter->second });
		}
	}

	// 5. ���з����αߵ�����Ƭ���ڽӹ�ϵ��
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


// ��˳�ǻ�·���ߣ�
template <typename DerivedV>
bool smoothCurve(Eigen::PlainObjectBase<DerivedV>& versOut, \
	const Eigen::MatrixBase<DerivedV>& versIn, const int times, \
	const double dbWeight = 0.5)
{
	/*
		����β���������⣬���ڵ㼯�е�ÿһ����p����ǰһ����Ϊq1����һ����Ϊq2��
		q1q2���ߵļ�Ȩ�е�Ϊq����q = q1+dbWeight*(q2-q1);
		ƽ���������ɵ��¶���p_newΪp��q���е㡣

 
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


// ��·���ƾ��Ȼ���
template <typename DerivedVO, typename DerivedVI>
bool arrangeLoop(Eigen::PlainObjectBase<DerivedVO>& loopOut, \
	Eigen::MatrixBase<DerivedVI>& loopIn, const int tarVersCount)
{
	const int dim = loopIn.cols();
	using ScalarO = typename DerivedVO::Scalar;
	using MatrixXO = Eigen::Matrix<ScalarO, Eigen::Dynamic, Eigen::Dynamic>;
	using RowVectorXO = Eigen::Matrix<ScalarO, 1, Eigen::Dynamic>;

	// 1.���������·�������߳���
	const int versCount = loopIn.rows();
	float length = 0;
	float eLen = 0;
	float tarEdgeLen = 0;								// �����ľ����Ļ�·ƽ���߳���
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

	// 2. �Գ��ȴ���0.5 * tarEdgeLen�ı߽��в�ֵ�� 
	int index = 0;
	int stepCount = 0;
	int versCountNew = 0;
	float threshold = 0.5 * tarEdgeLen;
	float step = 0;
	MatrixXO loopNew;											// ��ֵ���Ȼ��Ժ�Ļ�·���ƣ�
	RowVectorXO newVer;
	loopNew.resize(10 * std::max(versCount, tarVersCount), dim);
	loopNew.row(index++) = loopIn.row(0).cast<ScalarO>(); 
	for (int i = 0; i < versCount - 1; ++i)				// �����·������ѭ����
	{
		eLen = edgeLens[i];
		if (eLen > threshold)				// ��ԭ�߳�������ֵ�����㣺
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
	if (eLen > threshold)						// ��ԭ�߳�������ֵ�����㣺
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

	// 3. ��ֵ֮��Ļ�·���� 
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



/////////////////////////////////////////////////////////////////////////////////////////////////// ��������༭��

// �����������ϲ���������������һ�������� 
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


// ȷ�����������е����е���ͨ���򣨶�����ͨ��������1����ʹ��ϡ�����
template <typename DerivedI1, typename DerivedI2>
int simplyConnectedRegion(Eigen::PlainObjectBase<DerivedI1>& connectedLabels, \
	Eigen::PlainObjectBase<DerivedI2>& connectedCount, const Eigen::SparseMatrix<int>& adjSM)
{
	/*
		int simplyConnectedRegion(												���ص���ͨ���������������ʱ����-1
					Eigen::VectorXi& connectedLabels,						ͬһ��ͨ�����ڵĶ����ǩ������Ϊ0, 1, 2,...��
																											����Ϊi�Ķ���ı�ǩΪconnectedLabels(i);
					Eigen::VectorXi& connectedCount							ÿ����ǩ��Ӧ�ĵ���ͨ�����ڰ����Ķ�����
																											��ǩΪi�ĵ���ͨ��������Ķ�����ΪconnectedCount(i)
					const Eigen::SparseMatrix<Index>& adjSM,			����������ڽӾ���
					)
	*/
	if (adjSM.rows() == 0)
		return -1;

	if (adjSM.rows() != adjSM.cols())
		return -1;

	const int versCount = adjSM.rows();

	// 1. ��ʼ����
	connectedLabels.setConstant(versCount, versCount);		// ���ս���������ܵı�ǩΪversCount-1;
	connectedCount.setZero(versCount);

	// 2. �������ж��㣬����������ͨ�Ķ��㣺
	int currentLabel = 0;
	for (int i = 0; i < versCount; ++i)
	{
		// 2.1 ����ǰ�����ѱ����ʹ�(label < versCount)����continue:
		if (connectedLabels(i) < versCount)
			continue;

		// 2.2 ����ǰ����δ�����ʹ���ִ��BFS�ռ���������ͨ�Ķ��㣺
		std::queue<int> workingQueue;
		workingQueue.push(i);
		while (!workingQueue.empty())
		{
			const int curVerIdx = workingQueue.front();
			workingQueue.pop();

			if (connectedLabels(curVerIdx) < versCount)
				continue;

			// 2.2.1 ���׶���label��ֵ��label����+1
			connectedLabels(curVerIdx) = currentLabel;
			connectedCount(currentLabel)++;

			// 2.2.2 ���ڽӾ�����������ǰ�����ڽӵĶ��㣺
			for (typename Eigen::SparseMatrix<int>::InnerIterator iter(adjSM, curVerIdx); iter; ++iter)
			{
				const int connectVerIdx = iter.row();					// Ĭ��ϡ����������ȴ洢����ǰ�������ڵ�curVerIdx�У�

				if (connectedLabels(connectVerIdx) < versCount)
					continue;

				workingQueue.push(connectVerIdx);
			}
		}

		// 2.3 ��һ����ǩ�Ķ����ռ���ϣ���һ��ѭ���ռ���һ����ǩ�Ķ��㣺
		currentLabel++;
	}

	// 3. shrink_to_fit()
	connectedCount.conservativeResize(currentLabel, 1);

	return currentLabel;
}


// ȷ�����������е����е���ͨ���򣨶�����ͨ��������2����ʹ��std::unordered_set����ʹ��ϡ�����
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
		} while (blCollectedNew);						// ��������Ԫ�ر��ռ��������ִ��һ��

		return versSet.size();										// �����ռ��Ķ���������������Ϊlabel�Ķ��������
	};

	if (vers.rows() == 0 || tris.rows() == 0)
		return -1;
	connectedLabels.clear();
	connectedCount.clear();

	const int versCount = vers.rows();
	const int trisCount = tris.rows();

	// 1. ����������ѭ����
	connectedLabels.resize(versCount, versCount);		// ���ս���������ܵı�ǩΪversCount-1;
	std::unordered_set<int> versSet;
	std::vector<bool> triCollected(trisCount, false);
	int currentLabel = 0;
	int collectedVersCount = 0;
	while (collectedVersCount < versCount)
	{
		// w1. ѡȡ��һ��δ����ǵ�����Ƭ����������������Ϊ�������������Ӷ��㣻
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

		// w2. ����������ѭ����
		int currentLabelVersCount = collectVers(connectedLabels, triCollected, versSet, tris, currentLabel);

		// w3. post procedure:
		versSet.clear();
		collectedVersCount += currentLabelVersCount;
		currentLabel++;
	}

	// 2. ͳ�ƣ�
	int scCount = currentLabel;
	connectedCount.resize(scCount, 0);
	for (int i = 0; i < versCount; ++i)
	{
		int label = connectedLabels[i];
		connectedCount[label]++;
	}

	return scCount;
}


// ȷ�����������е����е���ͨ��������Ƭ��ͨ��������triangleGrow��ǿ������������Դ��з�����Ԫ�أ�
template <typename DerivedI>
int simplyTrisConnectedRegion(Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount, \
	const Eigen::MatrixBase<DerivedI>& tris)
{
	/*
	int simplyTrisConnectedRegion(											��������Ƭ����ͨ���������������ʱ����-1
				Eigen::VectorXi& connectedLabels,						ͬһ��ͨ�����ڵ�����Ƭ��ǩ������Ϊ0, 1, 2,...��
																										����Ϊi������Ƭ�ı�ǩΪconnectedLabels(i);
				Eigen::VectorXi& connectedCount							ÿ����ǩ��Ӧ�ĵ���ͨ�����ڰ���������Ƭ��
																										��ǩΪi�ĵ���ͨ�������������Ƭ��ΪconnectedCount(i)
				const Eigen::PlainObjectBase<DerivedI>& tris
				)
*/
	const unsigned trisCount = tris.rows();
	connectedLabels = Eigen::VectorXi::Zero(trisCount);
	unsigned visitedCount = 0;

	// 1. ������������Ƭ�ڽӹ�ϵ��
	Eigen::MatrixXi ttAdj_mnEdge, edges;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	getEdges(edges, tris);
	edges2adjMat(adjSM_eCount, adjSM_eIdx, edges);
	if (!getTrisAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, edges, adjSM_eCount, adjSM_eIdx))
		return -1;

	// 2. �ڽ�����Ƭ������ѭ����
	int currentLabel = 1;
	while (visitedCount < trisCount)
	{
		// w1. ѡȡ���ӣ�
		int seedIdx = 0;
		for (unsigned i = 0; i < trisCount; ++i)
		{
			if (0 == connectedLabels(i))
			{
				seedIdx = i;
				break;
			}
		}

		// w2. ��ʼ������
		std::unordered_set<int> workingQueue;
		workingQueue.insert(seedIdx);
		while (!workingQueue.empty())
		{
			// ww1. ����
			int cIdx = *workingQueue.begin();
			visitedCount++;
			workingQueue.erase(cIdx);
			connectedLabels(cIdx) = currentLabel;

			// ww2. ������ǰ����Ƭ�������������Ե�����Ƭ
			for (int i = 0; i < 3; ++i)
				if (ttAdj_mnEdge(cIdx, i) >= 0)
					if (0 == connectedLabels(ttAdj_mnEdge(cIdx, i)))			// ��ǩΪ0��ʾ������Ƭδ�����ʹ���
						workingQueue.insert(ttAdj_mnEdge(cIdx, i));					// ֻ����δ���ʵ�����Ƭ 

			// ww3. ������ǰ����Ƭ�ķ������������Ե�����Ƭ��
			for (int i = 0; i < 3; ++i)
				for (const auto& index : ttAdj_nmnOppEdge[cIdx][i])
					if (0 == index)
						workingQueue.insert(index);
		}

		// w3. ��ǰ��ͨ�����������
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


// ��ȡ���������������ͨ���򣨶�����ͨ��
template <typename DerivedV, typename DerivedI>
bool simplyConnectedLargest(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixBase<DerivedI>& tris)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();

	// 1. �����ڽӾ���
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx, adjSM;
	tris2adjMat(adjSM_eCount, adjSM_eIdx, tris);
	adjSM = adjSM_eCount;
	traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
		{
			if (iter.value() > 0)
				iter.valueRef() = 1;
		});

	// 2. ȷ�����������е���ͨ����
	Eigen::VectorXi connectedLabels, connectedCount;
	int conCount = simplyConnectedRegion(connectedLabels, connectedCount, adjSM);
	if (conCount < 0)
		return false;

	// 3. ȷ�������ͨ���򣨶�������ࣩ��
	int mainLabel = 0;								// �����ͨ����ı�ǩ��
	int mainLabelCount = 0;
	for (int i = 0; i < conCount; ++i)
	{
		if (connectedCount(i) > mainLabelCount)
		{
			mainLabel = i;
			mainLabelCount = connectedCount(i);
		}
	}

	// 4. ��ȡ�����ͨ����Ķ��㣺
	std::vector<int> mainLabelIdxes;
	mainLabelIdxes.reserve(versCount);
	for (unsigned i = 0; i < versCount; ++i)
		if (mainLabel == connectedLabels(i))
			mainLabelIdxes.push_back(i);
	mainLabelIdxes.shrink_to_fit();
	subFromIdxVec(versOut, vers, mainLabelIdxes);


	// 5. ��ȡ�����ͨ���������Ƭ��

	//		5.1 ������-������ӳ���
	std::vector<int> oldNewIdxInfo(versCount, -1);
	for (int i = 0; i < mainLabelIdxes.size(); ++i)
	{
		int oldIdx = mainLabelIdxes[i];
		oldNewIdxInfo[oldIdx] = i;
	}

	//		5.2 ����Ƭ�����еĶ�������ӳ�����������
	DerivedI trisCopy = tris;
	int* intPtr = trisCopy.data();
	for (int i = 0; i < trisCopy.size(); ++i)
	{
		int oldIdx = *intPtr;
		*intPtr = oldNewIdxInfo[oldIdx];
		intPtr++;
	}

	//		5.3 ��ȡ�����ͨ�����ڵ�����Ƭ��
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


// �����㵥��ͨ����ֽ�����
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

	// 0. Ԥ������ȥ���������㣬ȥ���Ƿ����ظ�����Ƭ��
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

	// 1. ��ȡ���㵥��ͨ���ƣ�
	const int versCount = vers0.rows();
	int scCount = 0;							// ���㵥��ͨ���������
	int sctCount = 0;							// ����Ƭ����ͨ���������
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

	// 2. ���ɶ��㵥��ͨ����
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
			*ptrData = oldNewIdxInfo[*ptrData];					// ������ӳ�����������
			ptrData++;
		}
		removeSickDupTris(compVers[i], compTris[i]);		// ȥ���Ƿ�����Ƭ��
	}

	return true;
}


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


// ��������еĹ�������
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


// ȥ�������еĹ������㣺
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

	// ������ƣ�
	subFromIdxVec(versOut, vers, newOldIdxInfo);

	// �������Ƭ��
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


// correctTriDirs()�����������������Ƭ���򣬻�������Ƭ���������ķ�����ʹ��ǰ��Ҫ�ȴ��������е��ص�����Ƭ��
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

	// 0. ���ʼ״̬��������������Ƭ�ķ���
	Eigen::MatrixXd triNorms;
	if (!getTriNorms(triNorms, vers, tris))				// �������˻�����Ƭ���ᵼ�¼�����Ƭ�������
		return false;

	// 1. ������Ƭ�ڽӹ�ϵ��
	Eigen::MatrixXi ttAdj_mnEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return false;

	// 2. ѭ�����������е�һ������Ƭt���� - tд��������� - ���t������������Ƭ������Ȼ�����
	int triIdx = 0;
	std::unordered_set<int> workingSet;								// ����Ѱ����������Ƭ�Ķ��У�
	std::vector<int> triTags(trisCount, 0);								// 0 : δ���ʣ�1: �ѷ��ʣ�-1: �ѷ��ʣ������Ϊ�����������Ƭ��
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

		workingSet.insert(static_cast<int>(triIdx));						// ���в����һ������Ƭ�� 
		while (workingSet.size() > 0)
		{
			// w1. �׸�Ԫ�س��ӣ����뵽��ͨ����Ƭ�����У�
			int cTriIdx = *workingSet.begin();									// ��ǰ���������Ƭ����
			workingSet.erase(workingSet.begin());
			if (triTags[cTriIdx] != 0)
				continue;
			triTags[cTriIdx] = 1;													// ��ǰ����Ƭ���Ϊ�ѷ��ʣ�
			visitedCount++;

			std::vector<int> adjTriIdxes;									// ��ǰ����Ƭ��δ�����ʵ���������Ƭ������
			Eigen::RowVector3d cTriNorm = triNorms.row(cTriIdx);
			Eigen::MatrixXd adjTriNorms;
			unsigned adjTrisCount = 0;

			// w2. ȷ����ǰ����Ƭ��������������Ƭ
			for (int i = 0; i < 3; ++i)
			{
				int nbrTriIdx = ttAdj_mnEdge(cTriIdx, i);
				if (nbrTriIdx >= 0)
					if (0 == triTags[nbrTriIdx])								// a. ��ǰ��Ϊ���α�
						adjTriIdxes.push_back(nbrTriIdx);
					else
					{
						// b. ��ǰ��Ϊ�����αߣ�
						auto& vec1 = ttAdj_nmnEdge[cTriIdx][i];						// ��ǰ�����α����ڵ���������Ƭ��
						auto& vec2 = ttAdj_nmnOppEdge[cTriIdx][i];				// ��ǰ�����αߵĶԱ����ڵ���������Ƭ��
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

			// 4.3 ����-������ǰ����Ƭ����������Ƭ��
			adjTrisCount = adjTriIdxes.size();
			if (0 == adjTrisCount)
				continue;
			subFromIdxVec(adjTriNorms, triNorms, adjTriIdxes);
			for (unsigned k = 0; k < adjTrisCount; ++k)
			{
				int adjIdx = adjTriIdxes[k];
				if (cTriNorm.dot(adjTriNorms.row(k)) < dotThreshold)
				{
					// ����Ƭ���Ϊ-1��
					triTags[adjIdx] = -1;
					visitedCount++;
					corrCount++;												// �������ӣ� 
				}
			}

			// 4.3 ���������в��뵱ǰ����Ƭδ�����ʵ���������Ƭ��
			for (const auto& index : adjTriIdxes)
				if (0 == index)
					workingSet.insert(index);
		}
	} 

	return corrCount;
}


// ȥ���Ƿ�����Ƭ���ظ�����Ƭ��
template<typename DerivedV>
int removeSickDupTris(const Eigen::MatrixBase<DerivedV>& vers, Eigen::MatrixXi& tris)
{
	using T = typename DerivedV::Scalar;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	int versCount = vers.rows();
	int trisCount = tris.rows();

	Eigen::VectorXi flags{ Eigen::VectorXi::Ones(trisCount) };			// flag��������Ҫȥ��������Ƭ��0��
	std::unordered_set<std::uint64_t> triCodeSet;

	// 1. ��������Ƭ���ҳ��Ƿ����ظ�����Ƭ��������flag�����б��Ϊ0��
	for (int i = 0; i < trisCount; ++i)
	{
		// ����Ƭ�е�����ֵ����0~versCount-1֮����Ϊ�Ƿ�����Ƭ��ò����ЩSTL�ļ��л�����������Σ�
		if (tris(i, 0) < 0 || tris(i, 0) > versCount - 1 || tris(i, 1) < 0 || tris(i, 1) > versCount - 1 || tris(i, 2) < 0 || tris(i, 2) > versCount - 1)
		{
			flags(i) = 0;
			continue;
		}

		std::uint64_t code = encodeTriangle(tris(i, 0), tris(i, 1), tris(i, 2));
		if (0 == code)
		{
			flags(i) = 0;				// ����Ƭ����������������ͬ�ģ��Ƿ���
			continue;
		}

		auto retPair = triCodeSet.insert(code);
		if (!retPair.second)
			flags(i) = 0;				// �ظ�����Ƭ�Ƿ���
	}

	// 2. ����flag������ȡ�µ�����Ƭ�������ԭ�е�;
	Eigen::MatrixXi tmpMat;
	subFromFlagVec(tmpMat, tris, eigenVec2Vec(flags));
	tris = tmpMat;

	// 3. ����ȥ��������Ƭ����
	int removeCount = trisCount - tris.rows();
	return removeCount;
}


template <typename DerivedV, typename DerivedI>
int findHoles(std::vector<std::vector<int>>& holes, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);


template <typename DerivedI>
bool checkSickTris(std::vector<int>& sickIdxes, const Eigen::PlainObjectBase<DerivedI>& tris);



// ��ʱδ�����ʵ�֣�
#include "temp.tpp"