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
bool getEdges(Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedI>& tris);
 

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

// ��������Ȩ��laplace����
template<typename Tl, typename DerivedV>
bool cotLaplacian(Eigen::SparseMatrix<Tl>& L, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::MatrixXi& tris);


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





// ��ʱδ�����ʵ�֣�
#include "temp.tpp"