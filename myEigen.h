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
 
#define USE_TRIANGLE_H

#ifdef USE_TRIANGLE_H
// ��algorithm����һ����ʹ�õ����ȣ�libigl���з�װ�������ʷ�ʹ�õ���˫���ȣ�
#define ANSI_DECLARATORS
#define REAL DOUBLE
#define VOID int
#include "triangulate.h"
#endif

#include "Eigen/Dense"
#include "Eigen/Sparse"


const double pi = 3.14159265359;


/////////////////////////////////////////////////////////////////////////////////////////////////// �����ṹ��
template <typename T>
struct triplet 
{
	T x;
	T y;
	T z;
};

template <typename T>
struct doublet
{
	T x;
	T y;
};


/////////////////////////////////////////////////////////////////////////////////////////////////// ǰ������

template <typename T>
std::pair<double, double> cart2polar(const T x, const T y);
template <typename T>
std::pair<double, double> polar2cart(const T radius, const T theta);
template<typename Func>
void PARALLEL_FOR(unsigned int  beg, unsigned int  end, const Func& func, const unsigned int serial_if_less_than);
template<typename Func, typename paramTuple>
void PARALLEL_FOR(unsigned int  beg, unsigned int  end, const Func& func, const paramTuple& pt, const unsigned int serial_if_less_than);
template<typename T, typename F>
void traverseSTL(T& con, F f);
template<typename T, typename F>
void revTraverseSTL(T& con, F f);
template <typename Derived, typename F>
void traverseMatrix(Eigen::PlainObjectBase<Derived>& m, F f);
template<typename spMat, typename F>
void traverseSparseMatrix(spMat& sm, F f);

template <typename T1, typename T2>
void dispPair(const std::pair<T1, T2>& pair);
template<typename T>
void dispQuat(const Eigen::Quaternion<T>& q);
template <typename Derived>
void dispMat(const Eigen::PlainObjectBase<Derived>& mat);
template <typename Derived>
void dispMatBlock(const Eigen::PlainObjectBase<Derived>& mat, const int rowStart, const int rowEnd, const int colStart, const int colEnd);

template<typename spMat>				// ���ģ�����ΪT, ��������ΪEigen::SparseMatrix<T>������ᱨ����֪��ʲôԵ�ʣ�
void dispSpMat(const spMat& sm, const unsigned startCol, const unsigned endCol);
template<typename spMat>
void dispSpMat(const spMat& sm, const unsigned startCol, const unsigned endCol, const unsigned showElemsCount);
template<typename spMat>
void dispSpMat(const spMat& sm);

template<typename T, int N>
void dispVec(const Eigen::Matrix<T, N, 1>& vec);
template<typename T, int N>
void dispVec(const Eigen::Matrix<T, 1, N>& vec);
template<typename T, int N>
void dispVecSeg(const Eigen::Matrix<T, N, 1>& vec, const int start, const int end);
template<typename T, int N>
void dispVecSeg(const Eigen::Matrix<T, 1, N>& vec, const int start, const int end);
template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& axisArrow, const float theta);
template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& originArrow, const Eigen::Matrix<T, 1, 3>& targetArrow);

template <typename T>
bool interpolateToLine(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::Matrix<T, 1, 3>& start, \
	const Eigen::Matrix<T, 1, 3>& end, const float SR, const bool SE);
template <typename T>
bool getCircleVers(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, const unsigned versCount);
template	<typename DerivedV, typename DerivedI>
void genCubeMesh(Eigen::Matrix<DerivedV, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::Matrix<DerivedI, Eigen::Dynamic, Eigen::Dynamic>& tris);
template <typename _Scalar, int _AmbientDim>
void genAABBmesh(const Eigen::AlignedBox<_Scalar, _AmbientDim>& aabb, Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::MatrixXi& tris);
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

template <typename Derived>
bool subFromIdxVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec);
template <typename Derived, typename Index>
bool subFromIdxVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const std::vector<Index>& vec);
template <typename Derived, typename IndexContainer>
bool subFromIdxCon(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const IndexContainer& con);
template <typename Derived>
bool subFromFlagVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec);
template<typename T, int N>
std::vector<T>  vec2Vec(const Eigen::Matrix<T, N, 1>& vIn);
template<typename T, int N>
Eigen::Matrix<T, N, 1> vec2Vec(const std::vector<T>& vIn);
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> vec2Vec(const std::vector<T>& vIn);
template <typename DerivedI>
void edges2mat(Eigen::PlainObjectBase<DerivedI>& mat, const std::vector<std::pair<int, int>>& edges);
template <typename Index>
std::int64_t encodeEdge(const Index vaIdx, const Index vbIdx);
template <typename Index>
std::int64_t encodeUedge(const Index vaIdx, const Index vbIdx);
std::pair<int, int> decodeEdge(const std::int64_t code);
template <typename Index>
std::uint64_t encodeTriangle(const Index vaIdx, const Index vbIdx, const Index vcIdx);
std::vector<int> decodeTrianagle(const std::uint64_t code);

template <typename DerivedI>
void findRepTris(std::vector<int>& repIdxes, const Eigen::PlainObjectBase<DerivedI>& tris);

template<typename T>
bool vecInsertNum(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const T num);
template<typename T>
bool vecInsertVec(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec2);
template<typename T>
bool matInsertRows(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat1);
template<typename T, int N>					// NΪ������
bool matInsertRows(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const Eigen::Matrix<T, 1, N>& rowVec);
template <typename Derived1, typename Derived2>
Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Derived1>& mat, const Eigen::PlainObjectBase<Derived2>& rowVec);
template <typename DerivedV, typename DerivedI>
void concatMeshMat(Eigen::PlainObjectBase<DerivedV>& vers, Eigen::PlainObjectBase<DerivedI>& tris, \
	const Eigen::PlainObjectBase<DerivedV>& vers1, const Eigen::PlainObjectBase<DerivedI>& tris1);

unsigned readNextData(char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize);
template <typename T>
void vecWriteToFile(const char* fileName, const std::vector<T>& vec);
template<typename T>
void vecReadFromFile(std::vector<T>& vec, const char* fileName, const unsigned elemCount);
template<typename T>
bool matWriteToFile(const char* fileName, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat);
template<typename T>
bool matReadFromFile(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName);
template	<typename Scalar, typename Index>
void objReadMeshMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::Matrix<Index, Eigen::Dynamic, \
	Eigen::Dynamic>& tris, const char* fileName);
template	<typename T>
void objWriteMeshMat(const char* fileName, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris);
template	<typename Scalar>
void objReadVerticesMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, const char* fileName);
template	<typename DerivedV>
void objWriteVerticesMat(const char* fileName, const Eigen::PlainObjectBase<DerivedV>& vers);
void printDirEigen(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& dir);
void printCoordinateEigen(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& xdir, \
	const Eigen::RowVector3f& ydir, const Eigen::RowVector3f& zdir);
template	<typename DerivedV, typename DerivedI>
void objWriteEdgesMat(const char* pathName, const Eigen::PlainObjectBase<DerivedI>& edges, \
	const Eigen::PlainObjectBase<DerivedV>& vers);
template <typename DerivedV, typename	 DerivedI>
void objWritePath(const char* pathName, const std::vector<DerivedI>& path, const Eigen::Matrix<DerivedV, Eigen::Dynamic, Eigen::Dynamic>& vers);
template <typename DerivedV>
void objWriteTreePath(const char* pathName, const Eigen::VectorXi& treeVec, const Eigen::Matrix<DerivedV, Eigen::Dynamic, Eigen::Dynamic>& vers);
void objReadVerticesHomoMat(Eigen::MatrixXf& vers, const char* fileName);

void objWriteVerticesHomoMat(const char* fileName, const Eigen::MatrixXf& vers);

void vers2homoVers(Eigen::MatrixXf& homoVers, const Eigen::MatrixXf& vers);

Eigen::MatrixXf vers2homoVers(const Eigen::MatrixXf& vers);

void homoVers2vers(Eigen::MatrixXf& vers, const Eigen::MatrixXf& homoVers);

Eigen::MatrixXf homoVers2vers(const Eigen::MatrixXf& homoVers);
template<typename T>
bool spMatTranspose(Eigen::SparseMatrix<T>& smOut, const Eigen::SparseMatrix<T>& smIn);
template<typename T, int N>
bool solveLinearEquation(Eigen::Matrix<T, N, 1>& x, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, const Eigen::Matrix<T, N, 1>& b);
template <typename T>
bool solveLinearEquations(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& X, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& B);
template<typename T, int N>
double hornersPoly(const Eigen::Matrix<T, N, 1>& coeffs, const double x);
template<typename T, int N>
float polyDiff(const Eigen::Matrix<T, N, 1>& coeffs, const float x);
template <typename T, typename Derived1, typename Derived2>
void kron(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& result, \
	const Eigen::MatrixBase<Derived1>& mm1, \
	const Eigen::MatrixBase<Derived2>& mm2);
template <typename Derived1, typename Derived2>
Eigen::MatrixXd kron(const Eigen::MatrixBase<Derived1>& mm1, \
	const Eigen::MatrixBase<Derived2>& mm2);
void polyInterpolation();
void gaussInterpolation();
void leastSquarePolyFitting();
template <typename T>
void ridgeRegressionPolyFitting(Eigen::VectorXd& theta, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, unsigned);
template<typename T>
Eigen::VectorXd fittingStandardEllipse(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& sampleVers);

template <typename DerivedI>
bool getEdges(Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedI>& tris);

template <typename DerivedV, typename DerivedF, typename DerivedN>
bool getTriNorms(Eigen::PlainObjectBase<DerivedN>& triNorms, const Eigen::MatrixBase<DerivedV>& vers, \
	const Eigen::MatrixBase<DerivedF>& tris);
template <typename DerivedV>
void getEdgeArrows(Eigen::PlainObjectBase<DerivedV>& arrows, const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers);
template<typename T, typename DerivedV, typename DerivedI>
bool trianglesBarycenter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& barys, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);
template<typename T, typename DerivedV, typename DerivedI>
bool trianglesNorm(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& triNorms, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);
template<typename DerivedV, typename DerivedI>
bool trianglesPlane(Eigen::MatrixXd& planeCoeff, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris);
template<typename Scalar, typename Ta, typename DerivedI>
bool trisArea(Eigen::Matrix<Ta, Eigen::Dynamic, 1>& trisAreaVec, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const Eigen::MatrixBase<DerivedI>& tris);
template<typename DerivedV, typename DerivedI>
double meshVolume(const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris);
template<typename Index>
bool adjMatrix(Eigen::SparseMatrix<Index>& adjSM_eCount, \
	Eigen::SparseMatrix<Index>& adjSM_eIdx, const Eigen::MatrixXi& tris, const Eigen::MatrixXi& edges0);
template<typename DerivedI>
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, const Eigen::PlainObjectBase<DerivedI>& tris);
template<typename DerivedI>
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, std::vector<int>& bdryTriIdxes, const Eigen::PlainObjectBase<DerivedI>& tris);
template<typename DerivedI>
bool nonManifoldEdges(Eigen::MatrixXi& nmnEdges, const Eigen::PlainObjectBase<DerivedI>& tris);
template<typename DerivedI>
int nonManifoldEdges(Eigen::MatrixXi& nmnEdges, std::vector<std::pair<int, std::pair<int, int>>>& nmnEdgeInfos, \
	const Eigen::PlainObjectBase<DerivedI>& tris);
using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool buildAdjacency(Eigen::MatrixXi& ttAdj_nmEdge, std::vector<tVec>& ttAdj_nmnEdge, std::vector<tVec>& ttAdj_nmnOppEdge, \
	const Eigen::PlainObjectBase<DerivedI>& tris);
template <typename DerivedV, typename DerivedI>
bool triangleGrow(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx);
template <typename DerivedV>
bool triangleGrowSplitMesh(std::vector<DerivedV>& meshesVersOut, std::vector<Eigen::MatrixXi>& meshesTrisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris);
template <typename T>
bool triangleGrowOuterSurf(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, const bool blRemIso);
template <typename DerivedV, typename DerivedI>
int correctTriDirs(Eigen::PlainObjectBase<DerivedI>& trisOut, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx, const double thetaThreshold);
template <typename DerivedV, typename DerivedI>
int findOverLapTris(std::vector<std::pair<int, int>>& opTrisPairs, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx, const double thetaThreshold);
template <typename DerivedI>
bool uEdgeGrow(std::vector<Eigen::MatrixXi>& circs, std::vector<Eigen::MatrixXi >& segs, \
	const Eigen::PlainObjectBase<DerivedI>& uEdges);
template <typename DerivedV, typename DerivedI>
int findHolesBdrySegs(std::vector<std::vector<int>>& holes, std::vector<std::vector<int>>& bdrySegs, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris);
template <typename DerivedI>
bool checkSickTris(std::vector<unsigned>& sickIdxes, const Eigen::PlainObjectBase<DerivedI>& tris);
template <typename Index>
int simplyConnectedRegion(const Eigen::SparseMatrix<Index>& adjSM,
	Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount);
template <typename DerivedI>
int simplyTrisConnectedRegion(Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount, \
	const Eigen::PlainObjectBase<DerivedI>& tris);
template <typename DerivedV, typename DerivedI>
bool simplyConnectedLargest(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris);
template <typename IndexT>
bool removeTris(Eigen::MatrixXi& trisOut, const Eigen::MatrixXi& tris, const std::vector<IndexT>& sickTriIdxes);
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
template<typename T>
bool genGrids(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, const Eigen::Matrix<T, 1, 3>& origin, \
	const float step, const std::vector<unsigned>& gridCounts);
template<typename DerivedV, typename DerivedI>
int removeSickDupTris(const Eigen::PlainObjectBase<DerivedV>& vers, Eigen::PlainObjectBase<DerivedI>& tris);
template <typename Tv, typename DerivedF, typename DerivedL>
void squared_edge_lengths(const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers,
	const Eigen::MatrixBase<DerivedF>& tris, Eigen::PlainObjectBase<DerivedL>& lenMat2);
template <typename Tv, typename DerivedF, typename Tl>
void trisAreaDB(const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixBase<DerivedF>& tris,
	Eigen::Matrix<Tl, Eigen::Dynamic, 1>& dbArea);
template <typename Tl>
void trisAreaDB(const Eigen::Matrix<Tl, Eigen::Dynamic, 3>& lenMat,
	const double nan_replacement, Eigen::Matrix<Tl, Eigen::Dynamic, 1>& dbArea);
template <typename Tv, typename Tl>
void trisCotValues(const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers,
	const Eigen::MatrixXi& tris, Eigen::Matrix<Tl, Eigen::Dynamic, Eigen::Dynamic>& cotValues);
template<typename Tv, typename Tl>
bool cotLaplacian(Eigen::SparseMatrix<Tl>& L, const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris);
template<typename T>
bool linearSpatialFilter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matOut, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matIn, \
	const Eigen::MatrixXd& mask);


/////////////////////////////////////////////////////////////////////////////////////////////////// debug�ӿڣ�

// ���������Ƭ����ת��Ϊtriplet��������ʽ��
template <typename T>
std::vector<triplet<T>> mat2triplets(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat)
{
	std::vector<triplet<T>> vec;
	if (0 == mat.rows() || 3 != mat.cols())
		return vec;

	vec.resize(mat.rows());
	for (unsigned i = 0; i < mat.rows(); ++i) 
	{
		vec[i].x = mat(i, 0);
		vec[i].y = mat(i, 1);
		vec[i].z = mat(i, 2);
	}

	return vec;
}

template <typename T>
std::vector<doublet<T>> mat2doublets(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat)
{
	std::vector<doublet<T>> vec;
	if (0 == mat.rows() || 2 != mat.cols())
		return vec;

	vec.resize(mat.rows());
	for (unsigned i = 0; i < mat.rows(); ++i)
	{
		vec[i].x = mat(i, 0);
		vec[i].y = mat(i, 1);
	}

	return vec;
}


// ��������triplet������������Ϊindex��tripletԪ��ʹ��ν��f����ֵΪtrue���򷵻�index; ���Ҳ���������򷵻�-1��
template <typename T, typename F>
int findTriplet(const std::vector<triplet<T>>& trips, F f)		
{
	// ν��F����ʽΪbool foo(const triplet<T>& t);
	if (trips.empty())
		return -1;
	for (unsigned i = 0; i < trips.size(); ++i)
		if (f(trips[i]))
			return static_cast<int>(i);

	return -1;
}

template <typename T, typename F>
int findTriplet(const std::vector<doublet<T>>& doubs, F f)
{
	// ν��F����ʽΪbool foo(const doublet<T>& d);
	if (doubs.empty())
		return -1;
	for (unsigned i = 0; i < doubs.size(); ++i)
		if (f(doubs[i]))
			return static_cast<int>(i);

	return -1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////// basic math interface

//		��ά�ѿ�������ϵת��Ϊ������ϵ
template <typename T>
std::pair<double, double> cart2polar(const T x, const T y)
{
	// theta���귶ΧΪ[0, 2*pi]��
	double x0 = static_cast<double>(x);
	double y0 = static_cast<double>(y);
	const double pi = 3.14159;
	double eps = 10e-9;

	double radius = std::sqrt(x0 * x0 + y0 * y0);

	if (std::abs(x0) < eps)
	{
		if (std::abs(y0) < eps)
			return { radius, NAN };
		else if (y > 0)
			return { radius, pi / 2 };
		else
			return { radius, 3 * pi / 2 };
	}

	if (std::abs(y0) < eps)
	{
		if (x > 0)
			return { radius, 0 };
		else
			return { radius, -pi };
	}

	// ��radius�ӽ�0�� �����Ӧ��theta��0��
	if (radius < eps)
		return { 0, 0 };

	double theta = std::acos(x / radius);
	if (y < 0)
		theta = 2 * pi - theta;

	return { radius, theta };
}


template <typename T>
std::pair<double, double> polar2cart(const T radius, const T theta)
{
	return { radius * cos(theta), radius * sin(theta) };
}



/////////////////////////////////////////////////////////////////////////////////////////////////// ����̨��ӡ�ӿ�

// ����forѭ��
template<typename Func>
void PARALLEL_FOR(unsigned int  beg, unsigned int  end, const Func& func, const unsigned int serial_if_less_than = 12)
{
	/*
		PARALLEL_FOR(
			unsigned int   beg,                                     ��ʼԪ������
			unsigned int  end,                                      βԪ������
			const unsigned int  serial_if_less_than,      �����Ҫ�����Ԫ�ز�С�ڸ�ֵ����ʹ�ò��л���
			const Func & func										����Ԫ�صĺ����ӣ�
			)
	*/
	unsigned int elemCount = end - beg + 1;

	if (elemCount < serial_if_less_than)
		for (unsigned int i = beg; i < end; ++i)
			func(i);
	else
	{
		// ȷ������߳�����
		const static unsigned n_threads_hint = std::thread::hardware_concurrency();
		const static unsigned n_threads = (n_threads_hint == 0u) ? 8u : n_threads_hint;

		// forѭ���ķ�Χ�ֽ�ɼ��Σ�
		unsigned int slice = (unsigned int)std::round(elemCount / static_cast<double>(n_threads));
		slice = std::max(slice, 1u);

		// �̺߳�����
		auto subTraverse = [&func](unsigned int head, unsigned int tail)
		{
			for (unsigned int k = head; k < tail; ++k)
				func(k);
		};

		// �����̳߳أ�ִ�в���forѭ����
		std::vector<std::thread> pool;              // �̳߳أ�
		pool.reserve(n_threads);
		unsigned int head = beg;
		unsigned int tail = std::min(beg + slice, end);
		for (unsigned int i = 0; i + 1 < n_threads && head < end; ++i)
		{
			pool.emplace_back(subTraverse, head, tail);
			head = tail;
			tail = std::min(tail + slice, end);
		}
		if (head < end)
			pool.emplace_back(subTraverse, head, end);

		// �߳�ͬ����
		for (std::thread& t : pool)
		{
			if (t.joinable())
				t.join();
		}
	}
}


// ��β���forѭ��������������Ĳ���ʹ��std::tuple���룻
template<typename Func, typename paramTuple>
void PARALLEL_FOR(unsigned int  beg, unsigned int  end, const Func& func, const paramTuple& pt, const unsigned int serial_if_less_than = 12)
{
	/*
		PARALLEL_FOR(
			unsigned int   beg,                                     ��ʼԪ������
			unsigned int  end,                                      βԪ������
			const unsigned int  serial_if_less_than,      �����Ҫ�����Ԫ�ز�С�ڸ�ֵ����ʹ�ò��л���
			const paramTuple& pt								�������������������
			const Func & func										����Ԫ�صĺ����ӣ�
			)
	*/
	unsigned int elemCount = end - beg + 1;

	if (elemCount < serial_if_less_than)
		for (unsigned int i = beg; i < end; ++i)
			func(i, pt);
	else
	{
		// ȷ������߳�����
		const static unsigned n_threads_hint = std::thread::hardware_concurrency();
		const static unsigned n_threads = (n_threads_hint == 0u) ? 8u : n_threads_hint;

		// forѭ���ķ�Χ�ֽ�ɼ��Σ�
		unsigned int slice = (unsigned int)std::round(elemCount / static_cast<double>(n_threads));
		slice = std::max(slice, 1u);

		// �̺߳�����
		auto subTraverse = [&func, &pt](unsigned int head, unsigned int tail)
		{
			for (unsigned int k = head; k < tail; ++k)
				func(k, pt);
		};

		// �����̳߳أ�ִ�в���forѭ����
		std::vector<std::thread> pool;              // �̳߳أ�
		pool.reserve(n_threads);
		unsigned int head = beg;
		unsigned int tail = std::min(beg + slice, end);
		for (unsigned int i = 0; i + 1 < n_threads && head < end; ++i)
		{
			pool.emplace_back(subTraverse, head, tail);
			head = tail;
			tail = std::min(tail + slice, end);
		}
		if (head < end)
			pool.emplace_back(subTraverse, head, end);

		// �߳�ͬ����
		for (std::thread& t : pool)
		{
			if (t.joinable())
				t.join();
		}
	}
}


// ���뺯���ӻ���ָ�����stl����
template<typename T, typename F>
void traverseSTL(T& con, F f)
{
	std::for_each(con.begin(), con.end(), f);
	std::cout << std::endl;
}


// �������
template<typename T, typename F>
void revTraverseSTL(T& con, F f)
{
	std::for_each(con.rbegin(), con.rend(), f);
	std::cout << std::endl;
}


template <typename Derived, typename F>
void traverseMatrix(Eigen::PlainObjectBase<Derived>& m, F f) 
{
	auto dataPtr = m.data();
	int elemCount = m.rows() * m.cols();
	for (int i = 0; i < elemCount; ++i)
		f(*dataPtr++);
}


// ���뺯���ӱ���ϡ������еķ���Ԫ�أ������ӽ��ܵĲ�����Eigen::SparseMatrix<T>::InnerIterator&
template<typename spMat, typename F>
void traverseSparseMatrix(spMat& sm, F f)
{
	for (unsigned i = 0; i < sm.outerSize(); ++i)
		for (auto iter = spMat::InnerIterator(sm, i); iter; ++iter)
			f(iter);
}


// lambda������ӡstd::cout֧�ֵ����ͱ�����
template <typename T>
auto disp = [](const T& arg)
{
	std::cout << arg << ", ";
};


template <typename T1, typename T2>
void dispPair(const std::pair<T1, T2>& pair) 
{
	std::cout << "(" << pair.first << ", " << pair.second << ")" << std::endl;
}


template<typename T>
void dispQuat(const Eigen::Quaternion<T>& q)
{
	std::cout << q.w() << std::endl;
	std::cout << q.vec() << std::endl;
}


template <typename Derived>
void dispMat(const Eigen::PlainObjectBase<Derived>& mat)
{
	std::cout << ": rows == " << mat.rows() << ",  cols == " << mat.cols() << std::endl;
	for (int i = 0; i < mat.rows(); ++i)
	{
		std::cout << i << " ---\t";
		for (int j = 0; j < mat.cols(); ++j)
			std::cout << mat.coeff(i, j) << ", ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

 
template <typename Derived>
void dispMatBlock(const Eigen::PlainObjectBase<Derived>& mat, const int rowStart, const int rowEnd, const int colStart, const int colEnd)
{
	if (rowEnd > mat.rows() - 1 || rowStart < 0 || rowStart >= rowEnd)
		return;
	if (colEnd > mat.cols() - 1 || colStart < 0 || colStart >= colEnd)
		return;

	std::cout << ": rows == " << mat.rows() << ",  cols == " << mat.cols() << std::endl;
	for (int i = rowStart; i <= rowEnd; ++i)
	{
		std::cout << i << " ---\t";
		for (int j = colStart; j <= colEnd; ++j)
			std::cout << mat(i, j) << ", ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


// ��ӡϡ������е����з���Ԫ�أ�
template<typename spMat>				// ���ģ�����ΪT, ��������ΪEigen::SparseMatrix<T>������ᱨ����֪��ʲôԵ�ʣ�
void dispSpMat(const spMat& sm, const unsigned startCol, const unsigned endCol)
{
	std::cout << "rows == " << sm.rows() << ", cols == " << sm.cols() << std::endl;
	for (unsigned i = startCol; i <= endCol; ++i)
		for (auto iter = spMat::InnerIterator(sm, i); iter; ++iter)
			std::cout << "(" << iter.row() << ", " << iter.col() << ") " << iter.value() << std::endl;
	std::cout << std::endl;
}


template<typename spMat>
void dispSpMat(const spMat& sm, const unsigned startCol, const unsigned endCol, const unsigned showElemsCount)
{
	unsigned count = 0;
	std::cout << "rows == " << sm.rows() << ", cols == " << sm.cols() << std::endl;
	for (unsigned i = startCol; i <= endCol; ++i)
		for (auto iter = spMat::InnerIterator(sm, i); iter; ++iter)
		{
			std::cout << "(" << iter.row() << ", " << iter.col() << ") " << iter.value() << std::endl;
			count++;
			if (count >= showElemsCount)
				break;
		}
	std::cout << std::endl;
}


template<typename spMat>
void dispSpMat(const spMat& sm)
{
	dispSpMat(sm, 0, sm.rows() - 1);
}
 

template<typename T, int N>
void dispVec(const Eigen::Matrix<T, N, 1>& vec)
{
	std::cout << ": rows == " << vec.rows() << std::endl;
	for (int i = 0; i < vec.rows(); ++i)
		std::cout << i << "---\t" << vec(i) << std::endl;
	std::cout << std::endl;
}


template<typename T, int N>
void dispVec(const Eigen::Matrix<T, 1, N>& vec)
{
	std::cout << ": cols == " << vec.cols() << std::endl;
	for (int i = 0; i < vec.cols(); ++i)
		std::cout << i << "---\t" << vec(i) << std::endl;
	std::cout << std::endl;
}


template<typename T, int N>
void dispVecSeg(const Eigen::Matrix<T, N, 1>& vec, const int start, const int end)
{
	if (start < 0 || start >= end)
	{
		std::cout << "input error." << std::endl;
		return;
	}

	int endIdx = (vec.rows() - 1 < end) ? (vec.rows() - 1) : end;
	std::cout << ": rows == " << vec.rows() << std::endl;
	for (int i = start; i <= endIdx; ++i)
		std::cout << i << "---\t" << vec(i) << std::endl;
	std::cout << std::endl;
}


template<typename T, int N>
void dispVecSeg(const Eigen::Matrix<T, 1, N>& vec, const int start, const int end)
{
	if (start < 0  || start >= end)
	{
		std::cout << "input error." << std::endl;
		return;
	}
	int endIdx = (vec.cols() - 1 < end) ? (vec.cols() -  1) : end;
	std::cout << ": cols == " << vec.cols() << std::endl;
	for (int i = start; i <= endIdx; ++i)
		std::cout << i << "---\t" << vec(i) << std::endl;
	std::cout << std::endl;
}



///////////////////////////////////////////////////////////////////////////////////////////////////// �ռ�任�ӿڣ�

// ע���������������ת����ȫ��Ĭ��Ϊ����������������v,uΪ��������r,sΪ��������R * v == u; �ȼ��� r * R.transpose() == s;

// getRotationMat() ����1����������ת����������ת�Ƕȣ�������ת����
template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& axisArrow, const float theta) 
{
	Eigen::Matrix<T, 3, 1> axis = axisArrow.transpose().normalized();
	return Eigen::AngleAxis<T>(theta, axis).toRotationMatrix();
}


// getRotationMat() ����2�����õ���originArrow��ת��targetArrow����ת����
template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& originArrow, const Eigen::Matrix<T, 1, 3>& targetArrow)
{
	Eigen::Matrix<T, 3, 3> rotation = Eigen::Matrix<T, 3, 3>::Zero();
	if (0 == originArrow.norm() || 0 == targetArrow.norm())
		return rotation;

	Eigen::Matrix<T, 1, 3> axisArrow = originArrow.cross(targetArrow);		// ��ת�᣻

	if (0 == axisArrow.norm())
		return Eigen::Matrix<T, 3, 3>::Identity();

	axisArrow.normalize();
	T x0 = axisArrow(0);
	T y0 = axisArrow(1);
	T z0 = axisArrow(2);
	T cosTheta = originArrow.dot(targetArrow) / (originArrow.norm() * targetArrow.norm());
	T sinTheta = std::sqrt(1 - cosTheta * cosTheta);

	// �ȼ���Eigen::AngleAxis<T>(theta, axis).toRotationMatrix()��������������������תtheta�Ƕȣ�
	rotation << cosTheta + (1 - cosTheta) * x0 * x0, (1 - cosTheta)* x0* y0 - sinTheta * z0, (1 - cosTheta)* x0* z0 + sinTheta * y0, \
		(1 - cosTheta)* y0* x0 + sinTheta * z0, cosTheta + (1 - cosTheta) * y0 * y0, (1 - cosTheta)* y0* z0 - sinTheta * x0, \
		(1 - cosTheta)* z0* x0 - sinTheta * y0, (1 - cosTheta)* z0* y0 + sinTheta * x0, cosTheta + (1 - cosTheta) * z0 * z0;
	return rotation;
}
 


///////////////////////////////////////////////////////////////////////////////////////////////////// ͼ�����ɽӿڣ�
 
// ������㡢�յ㡢�ռ�����ʣ���ֵ����һ��ֱ�ߵ��ƣ�
template <typename T>
bool interpolateToLine(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::Matrix<T, 1, 3>& start, \
			const Eigen::Matrix<T, 1, 3>& end, const float SR, const bool SE = true)
{
	if (vers.rows() > 0)
		return false;

	Eigen::Matrix<T, 1, 3> dir = end - start;
	float length = dir.norm();
	dir.normalize();

	if (length <= SR)
		return true;

	if (SE)
		matInsertRows<T, 3>(vers, start);

	float lenth0 = 0;
	for (unsigned i = 1; (length - lenth0) > 0.8 * SR; i++)			// ȷ�����һ��������յ㲻��̫����
	{
		Eigen::Matrix<T, 1, 3> temp = start + SR * i * dir;
		matInsertRows<T, 3>(vers, temp);
		lenth0 = SR * (i + 1);		// ��һ��temp�ĳ��ȡ�
	}

	if (SE)
		matInsertRows<T, 3>(vers, end);

	return true;
};


// ����XOYƽ���ڵ�ԲȦ�㼯��
template <typename T>
bool getCircleVers(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, const unsigned versCount = 20)
{
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	double deltaTheta = 2 * pi / versCount;
	vers.resize(versCount, 3);
	vers.setZero();
	for (unsigned i = 0; i < versCount; ++i)
	{
		double theta = deltaTheta * i;
		vers(i, 0) = radius * cos(theta);
		vers(i, 1) = radius * sin(theta);
	}

	return true;
}


// ����������ԭ�㣬�߳�Ϊ1������Ƭ��Ϊ12������������
template	<typename DerivedV, typename DerivedI>
void genCubeMesh(Eigen::Matrix<DerivedV, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::Matrix<DerivedI, Eigen::Dynamic, Eigen::Dynamic>& tris)
{
	vers.resize(8, 3);
	vers << -0.5000000, -0.5000000, -0.5000000, -0.5000000, 0.5000000, -0.5000000, \
		0.5000000, -0.5000000, -0.5000000, 0.5000000, 0.5000000, -0.5000000, \
		0.5000000, -0.5000000, 0.5000000, 0.5000000, 0.5000000, 0.5000000, \
		- 0.5000000, -0.5000000, 0.5000000, -0.5000000, 0.5000000, 0.5000000;

	tris.resize(12, 3);
	tris << 1, 2, 0, 1, 3, 2, 3, 4, 2, 3, 5, 4, 0, 4, 6, 0, 2, 4, 7, 3, 1, 7, 5, 3, 7, 0, 6, 7, 1, 0, 5, 6, 4, 5, 7, 6;
}


// ���������Χ�е���������
template <typename _Scalar, int _AmbientDim>
void genAABBmesh(const Eigen::AlignedBox<_Scalar, _AmbientDim>& aabb, Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::MatrixXi& tris)
{
	Eigen::Matrix<_Scalar, _AmbientDim, 1> minp = aabb.min();
	Eigen::Matrix<_Scalar, _AmbientDim, 1> maxp = aabb.max();
	Eigen::Matrix<_Scalar, _AmbientDim, 1> newOri = (minp + maxp) / 2.0;
	Eigen::Matrix<_Scalar, _AmbientDim, 1> sizeVec = maxp - minp;

	genCubeMesh(vers, tris);
	vers.col(0) *= sizeVec(0);
	vers.col(1) *= sizeVec(1);
	vers.col(2) *= sizeVec(2);
	vers.rowwise() += newOri.transpose();
}


#ifdef USE_TRIANGLE_H
// genCylinder()����1�������ɣ��ࣩ���壺
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers, \
		const bool isCovered = true)
{
	/*
	bool genCylinder(
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers,
		Eigen::MatrixXi& tris, 
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers,			���ߣ�
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers,			������·���㣬����Ҫ��XOYƽ���ڣ�
		const bool isCovered																						�Ƿ���
		)
	
	
	*/

	// lambda����������������Ƭ��������ѭ�����ã�����һ������һ�㣻
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// ���ظ����ã�tris��������ҪΪ�ա�
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// �����ɱ����Բ�����Ȧ�ĵ�һ�������������
		sideTris.conservativeResize(currentTrisCount + 2 * circVersCount, 3);

		int triIdx = currentTrisCount;
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx + circVersCount - 1, startIdx, startIdx + 2 * circVersCount - 1 };
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx, startIdx + circVersCount, startIdx + 2 * circVersCount - 1 };
		for (int i = startIdx + 1; i < startIdx + circVersCount; ++i)
		{
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i - 1, i, i + circVersCount - 1 };
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i, i + circVersCount, i + circVersCount - 1 };
		}

		return true;
	};

	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	unsigned circVersCount = btmVers.rows();					// �����һȦ�Ķ�������
	unsigned circCount = axisVers.rows();							// Ȧ����
	unsigned versCount = circVersCount * circCount;
	std::vector<MatrixXT> circuitsVec(circCount);									// ÿ��������һȦ���㣻
	std::vector<RowVector3T> sectionNorms(circCount);					// ÿ�������ķ���

	// 1. ��������circCount�������ķ���ÿ�������Ķ��㣻
	for (unsigned i = 0; i < circCount - 1; ++i)
	{
		sectionNorms[i] = axisVers.row(i + 1) - axisVers.row(i);
		sectionNorms[i].normalize();
		Matrix3T rotation = getRotationMat(RowVector3T{0, 0, 1}, sectionNorms[i]);
		circuitsVec[i] = btmVers * rotation.transpose();
		circuitsVec[i].rowwise() += axisVers.row(i);
	}

	// 2. �������һȦ���㣺
	RowVector3T deltaNormAve{RowVector3T::Zero()};
	for (unsigned i = 0; i < circCount - 2; ++i)
		deltaNormAve += (sectionNorms[i + 1] - sectionNorms[i]);
	deltaNormAve.array() /= (circCount - 2);
	sectionNorms[circCount - 1] = sectionNorms[circCount - 2] + deltaNormAve;
	Matrix3T rotation = getRotationMat(RowVector3T{ 0, 0, 1 }, sectionNorms[circCount - 1]);
	circuitsVec[circCount - 1] = btmVers * rotation.transpose();
	circuitsVec[circCount - 1].rowwise() += axisVers.row(circCount - 1);

	// 3. �������嶥�㣺
	vers.resize(versCount, 3);
	for (unsigned i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.���ɲ�������Ƭ��
	for (unsigned i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);

	
	// 5. �Ӹǣ�
	if (isCovered)
	{
		MatrixXT capVers;
		Eigen::MatrixXi capTris1, capTris2;
		circuit2mesh(capVers, capTris1, btmVers);
		capTris2 = capTris1;
		for (int i = 0; i < capTris1.rows(); ++i)
		{
			int tmp = capTris1(i, 2);
			capTris1(i, 2) = capTris1(i, 1);
			capTris1(i, 1) = tmp;
		}
		capTris2.array() += versCount - circVersCount;
		matInsertRows(tris, capTris1);
		matInsertRows(tris, capTris2);
	}

	return true;
}


// genCylinder()����2��������Բ��������
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const float radius, const double deltaTheta = 2 * pi / 30, const bool isCovered = true)
{
	// ����XOYƽ���ϲ����ǶȲ���Ϊ��ԲȦ���㣺
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	MatrixXT circuit(30, 3);
	circuit.setZero();
	for (unsigned i = 0; i < 30; ++i)
	{
		double theta = deltaTheta * i;
		circuit(i, 0) = radius * cos(theta);
		circuit(i, 1) = radius * sin(theta);
	}

	genCylinder(vers, tris, axisVers, circuit);

	return true;
}


// genCylinder()����3�������ɷ���������
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const std::pair<float, float> sizePair, const float SR = 0.5, \
	const bool isCovered = true)
{
	// ����XOYƽ���ϲ����ǶȲ���Ϊ��ԲȦ���㣺
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// ����XOYƽ���ڵķ��򶥵㣺
	float length = sizePair.first;
	float width = sizePair.second;
	std::vector<RowVector3T> corners(4);
	corners[0] = RowVector3T{ length / 2, width / 2, 0 };
	corners[1] = RowVector3T{ -length / 2, width / 2, 0 };
	corners[2] = RowVector3T{ -length / 2, -width / 2, 0 };
	corners[3] = RowVector3T{ length / 2, -width / 2, 0 };

	MatrixXT circuit, tmpVers1, tmpVers2, tmpVers3, tmpVers4;
	interpolateToLine(tmpVers1, corners[0], corners[1], SR, true);
	interpolateToLine(tmpVers2, corners[1], corners[2], SR, false);
	interpolateToLine(tmpVers3, corners[2], corners[3], SR, true);
	interpolateToLine(tmpVers4, corners[3], corners[0], SR, false);
	matInsertRows(circuit, tmpVers1);
	matInsertRows(circuit, tmpVers2);
	matInsertRows(circuit, tmpVers3);
	matInsertRows(circuit, tmpVers4);
 
	genCylinder(vers, tris, axisVers, circuit);

	return true;
}


//		���ɷ�������ת�����Σ���ȷ�������XOYƽ��ƽ�л�ֱ��
template <typename T>
bool genAlignedCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const std::pair<float, float> sizePair, const float SR, \
	const bool isCovered = true)
{
	// ����XOYƽ���ϲ����ǶȲ���Ϊ��ԲȦ���㣺
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// ����XOYƽ���ڵķ��򶥵㣺
	float length = sizePair.first;
	float width = sizePair.second;
	std::vector<RowVector3T> corners(4);
	corners[0] = RowVector3T{ length / 2, width / 2, 0 };
	corners[1] = RowVector3T{ -length / 2, width / 2, 0 };
	corners[2] = RowVector3T{ -length / 2, -width / 2, 0 };
	corners[3] = RowVector3T{ length / 2, -width / 2, 0 };

	MatrixXT circuit, tmpVers1, tmpVers2, tmpVers3, tmpVers4;
	interpolateToLine(tmpVers1, corners[0], corners[1], SR, true);
	interpolateToLine(tmpVers2, corners[1], corners[2], SR, false);
	interpolateToLine(tmpVers3, corners[2], corners[3], SR, true);
	interpolateToLine(tmpVers4, corners[3], corners[0], SR, false);
	matInsertRows(circuit, tmpVers1);
	matInsertRows(circuit, tmpVers2);
	matInsertRows(circuit, tmpVers3);
	matInsertRows(circuit, tmpVers4);

	// lambda����������������Ƭ��������ѭ�����ã�����һ������һ�㣻
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// ���ظ����ã�tris��������ҪΪ�ա�
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// �����ɱ����Բ�����Ȧ�ĵ�һ�������������
		sideTris.conservativeResize(currentTrisCount + 2 * circVersCount, 3);

		int triIdx = currentTrisCount;
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx + circVersCount - 1, startIdx, startIdx + 2 * circVersCount - 1 };
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx, startIdx + circVersCount, startIdx + 2 * circVersCount - 1 };
		for (int i = startIdx + 1; i < startIdx + circVersCount; ++i)
		{
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i - 1, i, i + circVersCount - 1 };
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i, i + circVersCount, i + circVersCount - 1 };
		}

		return true;
	};

	unsigned circVersCount = circuit.rows();					// �����һȦ�Ķ�������
	unsigned circCount = axisVers.rows();							// Ȧ����
	unsigned versCount = circVersCount * circCount;
	std::vector<MatrixXT> circuitsVec(circCount);									// ÿ��������һȦ���㣻
	std::vector<RowVector3T> sectionNorms(circCount);					// ÿ�������ķ���

	// 1. ��������circCount�������ķ���ÿ�������Ķ��㣻
	for (unsigned i = 0; i < circCount - 1; ++i)
	{
		sectionNorms[i] = axisVers.row(i + 1) - axisVers.row(i);
		sectionNorms[i].normalize();
		Matrix3T rotation1 = getRotationMat(RowVector3T{ 0, 0, 1 }, RowVector3T{ 0, 1, 0 });
		Matrix3T rotation2 = getRotationMat(RowVector3T{ 0, 1, 0 }, sectionNorms[i]);
		Matrix3T rotation = rotation2 * rotation1;
		circuitsVec[i] = circuit * rotation.transpose();
		circuitsVec[i].rowwise() += axisVers.row(i);
	}

	// 2. �������һȦ���㣺
	RowVector3T deltaNormAve{ RowVector3T::Zero() };
	for (unsigned i = 0; i < circCount - 2; ++i)
		deltaNormAve += (sectionNorms[i + 1] - sectionNorms[i]);
	deltaNormAve.array() /= (circCount - 2);
	sectionNorms[circCount - 1] = sectionNorms[circCount - 2] + deltaNormAve;
	Matrix3T rotation1 = getRotationMat(RowVector3T{ 0, 0, 1 }, RowVector3T{ 0, 1, 0 });
	Matrix3T rotation2 = getRotationMat(RowVector3T{ 0, 1, 0 }, sectionNorms[circCount - 1]);
	Matrix3T rotation = rotation2 * rotation1;
	circuitsVec[circCount - 1] = circuit * rotation.transpose();
	circuitsVec[circCount - 1].rowwise() += axisVers.row(circCount - 1);

	// 3. �������嶥�㣺
	vers.resize(versCount, 3);
	for (unsigned i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.���ɲ�������Ƭ��
	for (unsigned i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);


	// 5. �Ӹǣ�
	if (isCovered)
	{
		MatrixXT capVers;
		Eigen::MatrixXi capTris1, capTris2;
		circuit2mesh(capVers, capTris1, circuit);
		capTris2 = capTris1;
		for (int i = 0; i < capTris1.rows(); ++i)
		{
			int tmp = capTris1(i, 2);
			capTris1(i, 2) = capTris1(i, 1);
			capTris1(i, 1) = tmp;
		}
		capTris2.array() += versCount - circVersCount;
		matInsertRows(tris, capTris1);
		matInsertRows(tris, capTris2);
	}

	return true;
}



//			circuitToMesh()����1��triangle�������ʷ֡�����ձ߽��ߵ㼯�õ������񣬿�����ƽ��Ҳ���������棬����Ƭ�ߴ粻�ɿأ������������ڲ���㡣
template <typename T>
bool circuit2mesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& circVers)
{
	// ������ò�Ƶ�ǰ�����ʷ������⣡��
	using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;

	unsigned circCount = circVers.rows();
	unsigned versCount = circVers.rows();

	// 0. ��Ե��·�����������ݿ��뵽��������С�
	vers = circVers;

	// ����triangulate()�ĵ㼯��ͶӰ��XOYƽ��Ķ�ά�㣬ԭ�㼯Ӧ����ת�����ʵĽǶ���ͶӰ��

	// 1. ȡ��һ���㡢1/3���ĵ㡢2/3���ĵ������ƽ��ķ�������Ϊԭ�㼯�ķ�����
	RowVector3T vers1 = circVers.row(0);
	RowVector3T vers2 = circVers.row(versCount / 3);
	RowVector3T vers3 = circVers.row(2 * versCount / 3);
	RowVector3T norm = (vers1 - vers2).cross(vers3 - vers2);
	norm.normalize();

	//  2. ��ת�㼯ʹ��normƽ����z��
	Matrix3T rotation = getRotationMat(norm, RowVector3T{0, 0, 1});
	vers = (vers * rotation.transpose()).eval();
 
	// 3. ��ת��ĵ�����д�뵽triangulate()�ӿڵ�����ṹ���С�
	Eigen::MatrixXd vers2D;
	Eigen::MatrixXi edges2D;
	vers2D = vers.transpose().topRows(2).cast<double>();
	edges2D.resize(2, versCount);
	for (unsigned i = 0; i < versCount; i++)
	{
		edges2D(0, i) = i + 1;
		edges2D(1, i) = (i + 1) % versCount + 1;
	}

	triangulateio triIn, triOut;
	triIn.numberofpoints = versCount;
	triIn.pointlist = (double*)vers2D.data();
	triIn.numberofpointattributes = 0;
	triIn.pointattributelist = NULL;
	triIn.pointmarkerlist = NULL;

	triIn.numberofsegments = versCount;
	triIn.segmentlist = (int*)edges2D.data();
	triIn.segmentmarkerlist = NULL;

	triIn.numberoftriangles = 0;
	triIn.numberofcorners = 0;
	triIn.numberoftriangleattributes = 0;

	triIn.numberofholes = 0;
	triIn.holelist = NULL;

	triIn.numberofregions = 0;
	triIn.regionlist = NULL;
	memset(&triOut, 0, sizeof(triangulateio));

	// 4. ִ�ж�ά�����ʷ֣��õ������������Ƭ���ݣ���ȡ�����������ݣ�������������ʹ����ת����ǰ�ġ�
	char triStr[256] = "pY";
	triangulate(triStr, &triIn, &triOut, NULL);
	tris.resize(3, triOut.numberoftriangles);
	MatrixXT norms(versCount, 3);
	std::memcpy(tris.data(), triOut.trianglelist, sizeof(int) * 3 * triOut.numberoftriangles);
	tris.transposeInPlace();
	tris.array() -= 1;
 
	return true;
}
#endif


///////////////////////////////////////////////////////////////////////////////////////////////////// ��ͬ�������͵ı任

// ��������������Դ��������ȡԪ�������������
template <typename Derived>
bool subFromIdxVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec)
{
	Derived& matOut = matBaseOut.derived();
	matOut.resize(vec.size(), matBaseIn.cols());
	for (unsigned i = 0; i < vec.rows(); ++i)
	{
		const int& index = vec(i);
		matOut.row(i) = matBaseIn.row(index);
	}

	return true;
}
 

template <typename Derived, typename Index>
bool subFromIdxVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const std::vector<Index>& vec)
{
	Derived& matOut = matBaseOut.derived();
	matOut.resize(vec.size(), matBaseIn.cols());
	for (unsigned i = 0; i < vec.size(); ++i)
	{
		const Eigen::Index& index = vec[i];
		matOut.row(i) = matBaseIn.row(index);
	}

	return true;
}


// ���������е�������Դ��������ȡԪ�������������
template <typename Derived, typename IndexContainer>
bool subFromIdxCon(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const IndexContainer& con)
{
	Derived& matOut = matBaseOut.derived();
	matOut.resize(con.size(), matBaseIn.cols());

	auto iter = con.begin();
	for (unsigned i = 0; iter != con.end(); ++i)
	{
		const Eigen::Index& index = *iter++;
		matOut.row(i) = matBaseIn.row(index);
	}

	return true;
}
 

// ����flag������Դ��������ȡԪ�������������
template <typename Derived>
bool subFromFlagVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec)
{
	Derived& matOut = matBaseOut.derived();
	matOut.resize(vec.sum(), matBaseIn.cols());

	int count = 0;
	for (unsigned i = 0; i < vec.rows(); ++i)
		if (vec(i) > 0)
			matOut.row(count++) = matBaseIn.row(i);

	return true;
}


// eigen��������std::vector<T>�໥ת��
template<typename T, int N>
std::vector<T>  vec2Vec(const Eigen::Matrix<T, N, 1>& vIn)
{
	unsigned elemCount = vIn.rows();
	std::vector<T> vOut(elemCount, 1);
	std::memcpy(&vOut[0], vIn.data(), sizeof(T) * elemCount);

	return vOut;
}

 
template<typename T, int N>
Eigen::Matrix<T, N, 1> vec2Vec(const std::vector<T>& vIn)
{
	unsigned elemCount = vIn.size();
	Eigen::Matrix<T, N, 1> vOut;
	vOut.resize(elemCount, 1);
	std::memcpy(vOut.data(), &vIn[0], sizeof(T) * elemCount);

	return vOut;
}
 

template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> vec2Vec(const std::vector<T>& vIn)
{
	unsigned elemCount = vIn.size();
	Eigen::Matrix<T, Eigen::Dynamic, 1> vOut;
	vOut.resize(elemCount, 1);
	std::memcpy(vOut.data(), &vIn[0], sizeof(T) * elemCount);

	return vOut;
}


template <typename DerivedI>
void edges2mat(Eigen::PlainObjectBase<DerivedI>& mat, const std::vector<std::pair<int, int>>& edges)
{
	mat.resize(edges.size(), 2);
	for (unsigned i = 0; i < edges.size(); ++i)
	{
		mat(i, 0) = edges[i].first;
		mat(i, 1) = edges[i].second;
	}
}


// ���ɱ߱��롪���������˵��������ӳ���һ��64λ������
template <typename Index>
std::int64_t encodeEdge(const Index vaIdx, const Index vbIdx)
{
	std::int64_t a = static_cast<std::int64_t>(vaIdx);
	std::int64_t b = static_cast<std::int64_t>(vbIdx);
	std::int64_t code = 0;
	code |= (a << 32);
	code |= b;
	return code;
}


// ��������߱��롪���������˵�������������ǰС����˳��ӳ���һ��64λ������
template <typename Index>
std::int64_t encodeUedge(const Index vaIdx, const Index vbIdx)
{
	Index vaIdx0 = (vaIdx < vbIdx ? vaIdx : vbIdx);
	Index vbIdx0 = (vaIdx < vbIdx ? vbIdx : vaIdx);
	std::int64_t a = static_cast<std::int64_t>(vaIdx0);
	std::int64_t b = static_cast<std::int64_t>(vbIdx0);
	std::int64_t code = 0;
	code |= (a << 32);
	code |= b;
	return code;
}


// ����߱��룻
std::pair<int, int> decodeEdge(const std::int64_t code);


// ��������Ƭ���롪�������������������ӳ���64λ�޷���������
template <typename Index>
std::uint64_t encodeTriangle(const Index vaIdx, const Index vbIdx, const Index vcIdx)
{
	unsigned long triIdxLimit = 0x1FFFFF;								// ���Ϊ21λȫ1, 0x1FFFFF ==  2097181, ���ٶ�������Ƭ��
	if (vaIdx > triIdxLimit || vbIdx > triIdxLimit || vcIdx > triIdxLimit)
		return 0;			// ����������Χ
	if (vaIdx == vbIdx || vaIdx == vcIdx || vbIdx == vcIdx || vaIdx < 0 || vbIdx < 0 || vcIdx < 0)
		return 0;			// �Ƿ�����Ƭ

	// ���������桪����(a, b, c)��(a, c, b)��ӳ��ɲ�ͬ�ı��룻����(a,b,c)(b,c,a)(c,a,b)ӳ�����ͬ�ı��룻
	Index A, B, C;
	if (vaIdx < vbIdx && vaIdx < vcIdx)
	{
		A = vaIdx; B = vbIdx; C = vcIdx;			// abc
	}						
	else if (vbIdx < vaIdx && vbIdx < vcIdx)
	{
		A = vbIdx; B = vcIdx; C = vaIdx;						// bca
	}
	else
	{
		A = vcIdx; B = vaIdx; C = vbIdx;						// cab;
	}
	std::uint64_t code = 0;
	code |= (static_cast<std::uint64_t>(A) << 42);
	code |= (static_cast<std::uint64_t>(B) << 21);
	code |= static_cast<std::uint64_t>(C);
	return code;
	 
}


// ��������Ƭ���룺
std::vector<int> decodeTrianagle(const std::uint64_t code);


template <typename DerivedI>
void findRepTris(std::vector<int>& repIdxes, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	std::unordered_set<std::int64_t> codeSet;
	unsigned trisCount = tris.rows();
	repIdxes.reserve(trisCount);

	for (unsigned i = 0; i<trisCount; ++i) 
	{
		std::int64_t code = encodeTriangle(tris(i, 0), tris(i, 1), tris(i, 2));
		auto retPair = codeSet.insert(code);
		if (!retPair.second)
			repIdxes.push_back(static_cast<int>(i));
	}
	repIdxes.shrink_to_fit();
}




///////////////////////////////////////////////////////////////////////////////////////////////////////// �������ɾ���

// ������������
template<typename T>
bool vecInsertNum(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const T num)
{
	vec.conservativeResize(vec.rows() + 1, 1);
	vec(vec.rows() - 1) = num;
	return true;
}


// ������������
template<typename T>
bool vecInsertVec(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec2)
{
	unsigned currentRows = vec1.rows();
	unsigned addRows = vec2.rows();
	unsigned finalRows = currentRows + addRows;
	vec1.conservativeResize(finalRows);

	// �������ݣ�
	T* dataPtr = vec1.data();
	dataPtr += currentRows;
	std::memcpy(dataPtr, vec2.data(), sizeof(T) * addRows);

	return true;
}


//	����ĩ�˲������
template<typename T>
bool matInsertRows(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat1)
{
	unsigned cols = mat1.cols();
	unsigned currentRows = mat.rows();
	unsigned addRows = mat1.rows();
	mat.conservativeResize(currentRows + addRows, cols);
	for (unsigned i = 0; i < addRows; ++i)
		mat.row(currentRows + i) = mat1.row(i);

	return true;
}


//	����ĩ�˲���������
template<typename T, int N>					// NΪ������
bool matInsertRows(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const Eigen::Matrix<T, 1, N>& rowVec)
{
	unsigned cols = N;
	unsigned currentRows = mat.rows();

	if (0 == rowVec.cols())
		return false;

	if (N != mat.cols() && mat.rows() != 0)
		return false;

	mat.conservativeResize(currentRows + 1, cols);
	mat.row(currentRows) = rowVec;

	return true;
}


// ����һ��flag������retVec����mat�ĵ�i�к�������vec��ȣ���retVec(i)==1���������0�������������retVec����Ԫ��Ϊ-1
template <typename Derived1, typename Derived2>
Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Derived1>& mat, const Eigen::PlainObjectBase<Derived2>& rowVec)
{
	int rows = mat.rows();
	int cols = mat.cols();

	Eigen::VectorXi retVec(rows);
	if (rowVec.cols() != cols)
	{
		retVec = -Eigen::VectorXi::Ones(rows);
		return retVec;
	}

	// ���бȽϣ�
	Eigen::MatrixXi tempMat(rows, cols);
	for (int i = 0; i < cols; ++i)
		tempMat.col(i) = (mat.col(i).array() == rowVec(i)).select(Eigen::VectorXi::Ones(rows), Eigen::VectorXi::Zero(rows));

	retVec = tempMat.col(0);

	if (cols > 1)
		for (int i = 1; i < cols; ++i)
			retVec = retVec.array() * tempMat.col(i).array();			// ������ˣ�

	return retVec;
}


// �����������ϲ���������������һ��������
template <typename DerivedV, typename DerivedI>
void concatMeshMat(Eigen::PlainObjectBase<DerivedV>& vers, Eigen::PlainObjectBase<DerivedI>& tris, \
	const Eigen::PlainObjectBase<DerivedV>& vers1, const Eigen::PlainObjectBase<DerivedI>& tris1)
{
	int versCount = vers.rows();
	matInsertRows(vers, vers1);

	DerivedI trisCopy1 = tris1;
	int* intPtr = trisCopy1.data();
	for (int i = 0; i < trisCopy1.size(); ++i)
		*(intPtr++) = versCount + *intPtr;

	matInsertRows(tris, trisCopy1);
};

// MATLAB����repmat();
template <typename DerivedA, typename DerivedB>
void repmat(Eigen::PlainObjectBase<DerivedB>& B, const Eigen::PlainObjectBase<DerivedA>& A, const int repRows, const int repCols)
{
	assert(repRows > 0);
	assert(repCols > 0);
	B.resize(repRows * A.rows(), repCols * A.cols());

	// copy tiled blocks
	for (int i = 0; i < repRows; i++)
		for (int j = 0; j < repCols; j++)
			B.block(i * A.rows(), j * A.cols(), A.rows(), A.cols()) = A;
}

// MATLAB����sparse()
template <class IndexVectorI, class IndexVectorJ,	class ValueVector,	typename T>
void sparse(Eigen::SparseMatrix<T>& SM, const IndexVectorI& I, const IndexVectorJ& J,
		const ValueVector& values, const size_t m, const size_t n	)
{
	assert((int)I.maxCoeff() < (int)m);
	assert((int)I.minCoeff() >= 0);
	assert((int)J.maxCoeff() < (int)n);
	assert((int)J.minCoeff() >= 0);
	assert(I.size() == J.size());
	assert(J.size() == values.size());

	// Really we just need .size() to be the same, but this is safer
	assert(I.rows() == J.rows());
	assert(J.rows() == values.rows());
	assert(I.cols() == J.cols());
	assert(J.cols() == values.cols());

	std::vector<Eigen::Triplet<T>> IJV;
	IJV.reserve(I.size());
	for (int x = 0; x < I.size(); x++)
		IJV.push_back(Eigen::Triplet<T>(I(x), J(x), values(x)));
	SM.resize(m, n);
	SM.setFromTriplets(IJV.begin(), IJV.end());
}



//////////////////////////////////////////////////////////////////////////////////////////// IO�ӿ�

unsigned readNextData(char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize);


// ��std::vector�е�����д�뵽�ļ��У���������ʽ��vector�е�Ԫ�ز�����pointer-like����
template <typename T>
void vecWriteToFile(const char* fileName, const std::vector<T>& vec)
{
	std::ofstream file(fileName, std::ios_base::out | std::ios_base::binary);
	file.write((char*)(&vec[0]), sizeof(T) * vec.size());
	file.close();
}


// to be optimized����Ԫ����Ӧ�ò���Ҫ���뺯����Ӧ�ÿ������ļ�����С�ƶϳ�����
template<typename T>
void vecReadFromFile(std::vector<T>& vec, const char* fileName, const unsigned elemCount)
{
	std::ifstream file(fileName, std::ios_base::in | std::ios_base::binary);
	vec.resize(elemCount);
	file.read((char*)(&vec[0]), sizeof(T) * vec.size());
	file.close();
}


template<typename T>
bool matWriteToFile(const char* fileName, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat)
{
	std::ofstream file(fileName);
	file << "row " << mat.rows() << std::endl;
	file << "col " << mat.cols() << std::endl;
	const T* ptr = mat.data();
	unsigned elemCount = mat.rows() * mat.cols();
	for (unsigned i = 0; i < elemCount; ++i)
	{
		file << *ptr++ << std::endl;
	}
	file.close();

	return true;
};


template<typename T>
bool matReadFromFile(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName)
{
	std::ifstream file(fileName);
	const unsigned LINE_LENGTH = 100;
	char cstr[100];
	unsigned lineOrder = 0;
	unsigned row = 0;
	unsigned col = 0;
	std::string str1, str2;

	// ������ߴ���Ϣ
	file.getline(cstr, LINE_LENGTH);
	str1 = cstr;
	str2.insert(str2.end(), str1.begin(), str1.begin() + 3);
	if (std::strcmp(str2.c_str(), "row") == 0)
	{
		str2.clear();
		str2.insert(str2.end(), str1.begin() + 3, str1.end());
		row = std::stoi(str2);
	}
	else
		return false;

	file.getline(cstr, LINE_LENGTH);
	str1 = cstr;
	str2.clear();
	str2.insert(str2.end(), str1.begin(), str1.begin() + 3);
	if (std::strcmp(str2.c_str(), "col") == 0)
	{
		str2.clear();
		str2.insert(str2.end(), str1.begin() + 3, str1.end());
		col = std::stoi(str2);
	}
	else
		return false;

	// ������Ԫ��
	mat.resize(row, col);
	T* ptr = mat.data();
	while (!file.eof())
	{
		file.getline(cstr, LINE_LENGTH);
		str1 = cstr;
		if (str1.size() == 0)
			break;
		std::string::iterator iter = str1.begin();
		for (unsigned j = 0; j < 3; ++j)
		{
			if (iter == str1.end())
				break;

			if (*iter == ' ')						// ���ź�������пո���Ҫȥ��
				iter = str1.erase(iter);
			iter++;
		}
		*ptr++ = std::stoi(str1);
	}
	file.close();
};


template	<typename Scalar, typename Index>
void objReadMeshMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic>& tris, const char* fileName)
{
	char* pTmp = NULL;
	std::ifstream ifs(fileName);		//cube bunny Eight
	if (false == ifs.is_open())
		return;

	std::streampos   pos = ifs.tellg();     //   save   current   position   
	ifs.seekg(0, std::ios::end);
	unsigned fileLen = (unsigned)ifs.tellg();
	if (0 == fileLen)
		return;

	ifs.seekg(pos);     //   restore   saved   position   
	char* pFileBuf = new char[fileLen + 1];
	std::memset(pFileBuf, 0, fileLen + 1);
	ifs.read(pFileBuf, fileLen);
	char tmpBuffer[1024];
	unsigned nMaxSize = 1024;
	pTmp = pFileBuf;
	unsigned nReadLen = 0;
	unsigned nRet = 0;

	vers.resize(0, 0);
	tris.resize(0, 0);
	int versCols = vers.cols();
	int trisCols = tris.cols();

	while (nReadLen < fileLen)
	{
		nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
		if (0 == nRet)
			break;

		if (std::strcmp(tmpBuffer, "v") == 0)
		{
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;

			Eigen::Matrix<Scalar, 3, 1> ver;
			ver(0) = (Scalar)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(1) = (Scalar)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(2) = (Scalar)atof(tmpBuffer);

			vers.conservativeResize(3, ++versCols);
			vers.col(versCols - 1) = ver;
		}
		else if (std::strcmp(tmpBuffer, "f") == 0)
		{
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;

			Eigen::Matrix<Index, 3, 1> tri;
			// Vector3i tri;
			tri(0) = atoi(tmpBuffer) - 1;
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			tri(1) = atoi(tmpBuffer) - 1;

			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			tri(2) = atoi(tmpBuffer) - 1;
			tris.conservativeResize(3, ++trisCols);
			tris.col(trisCols - 1) = tri;
		}
	}
	vers.transposeInPlace();
	tris.transposeInPlace();
	delete[] pFileBuf;
};


template	<typename T>
void objWriteMeshMat(const char* fileName, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
{
	std::ofstream dstFile(fileName);
	if (vers.cols() != 3 || tris.cols() != 3)
		return;

	for (int j = 0; j < vers.rows(); j++)
	{
		char szBuf[1024] = { 0 };
		sprintf_s(szBuf, 1024, "v %.17g %.17g %.17g", vers(j, 0), vers(j, 1), vers(j, 2));
		dstFile << szBuf << "\n";
	}

	for (unsigned j = 0; j < tris.rows(); ++j)
	{
		char szBuf[256] = { 0 };
		sprintf_s(szBuf, 256, "f %d %d %d", tris(j, 0) + 1, tris(j, 1) + 1, tris(j, 2) + 1);
		dstFile << szBuf << "\n";
	}
};
 

template	<typename Scalar>
void objReadVerticesMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, const char* fileName)
{
	char* pTmp = NULL;
	std::ifstream ifs(fileName);			// cube bunny Eight
	if (false == ifs.is_open())
		return;
	
	std::streampos   pos = ifs.tellg();			//  save   current   position   
	ifs.seekg(0, std::ios::end);
	unsigned fileLen = (unsigned)ifs.tellg();
	if (0 == fileLen)
		return;

	ifs.seekg(pos);				  //   restore   saved   position   
	char* pFileBuf = new char[fileLen + 1];
	std::memset(pFileBuf, 0, fileLen + 1);
	ifs.read(pFileBuf, fileLen);
	char tmpBuffer[1024];
	unsigned nMaxSize = 1024;
	pTmp = pFileBuf;
	unsigned nReadLen = 0;
	unsigned nRet = 0;

	vers.resize(0, 0);
	int cols = vers.cols();
	while (nReadLen < fileLen)
	{
		nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
		if (0 == nRet)
			break;

		// ������Ϣ		
		if (std::strcmp(tmpBuffer, "v") == 0)
		{
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;

			Eigen::Matrix<Scalar, 3, 1> ver(Eigen::Matrix<Scalar, 3, 1>::Zero());
 
			ver(0) = (Scalar)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(1) = (Scalar)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(2) = (Scalar)atof(tmpBuffer);

			vers.conservativeResize(3, ++cols);
			vers.col(cols - 1) = ver;
		}
		else
			break;
	}
	vers.transposeInPlace();
	delete[] pFileBuf;
};


template	<typename DerivedV>
void objWriteVerticesMat(const char* fileName, const Eigen::PlainObjectBase<DerivedV>& vers)
{
	std::ofstream dstFile(fileName);
	for (int i = 0; i < vers.rows(); i++)
		dstFile << "v " << vers.coeffRef(i, 0) << " " << vers.coeffRef(i, 1) << " " << vers.coeffRef(i, 2) << std::endl;
	dstFile.close();
};


void printDirEigen(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& dir);


void printCoordinateEigen(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& xdir, \
	const Eigen::RowVector3f& ydir, const Eigen::RowVector3f& zdir);


// ������д�뵽OBJ�ļ��У�
template	<typename DerivedV, typename DerivedI>
void objWriteEdgesMat(const char* pathName, const Eigen::PlainObjectBase<DerivedI>& edges, \
	const Eigen::PlainObjectBase<DerivedV>& vers)
{
	if (0 == edges.rows() || 0 == vers.rows())
		return;
	std::ofstream dstFile(pathName);
	for (int i = 0; i < vers.rows(); ++i)
		dstFile << "v " << vers.coeffRef(i, 0) << " " << vers.coeffRef(i, 1) << " " << vers.coeffRef(i, 2) << std::endl;

	for (int i = 0; i < edges.rows(); ++i)
		dstFile << "l " << edges.coeffRef(i, 0) + 1 << " " << edges.coeffRef(i, 1) + 1 << std::endl;

	dstFile.close();
}



// objWritePath() ·������д�뵽OBJ�ļ��У�
template <typename DerivedV, typename	 DerivedI>
void objWritePath(const char* pathName, const std::vector<DerivedI>& path, const Eigen::Matrix<DerivedV, Eigen::Dynamic, Eigen::Dynamic>& vers)
{
	if (path.size() <= 1)
		return;

	unsigned edgesCount = path.size() - 1;
	Eigen::Matrix<DerivedI, Eigen::Dynamic, Eigen::Dynamic> pathEdges(edgesCount, 2);

	for (unsigned i = 0; i < edgesCount; ++i)
		pathEdges(i, 0) = path[i];
	for (unsigned i = 0; i < edgesCount; ++i)
		pathEdges(i, 1) = path[i + 1];
	objWriteEdgesMat(pathName, pathEdges, vers);
}
 

// objWirteTreePath()����������������·��������д�뵽OBJ�ļ��У�
template <typename DerivedV>
void objWriteTreePath(const char* pathName, const Eigen::VectorXi& treeVec, const Eigen::Matrix<DerivedV, Eigen::Dynamic, Eigen::Dynamic>& vers)
{
	// ·������Ϊ����������
	
	// ��������Ϊi�Ķ����ǰ����������ΪtreeVec(i), ����û��ǰ�����㣨���ڵ㣩����treeVec(i) == -1;
	if (treeVec.size() <= 1)
		return;

	unsigned edgesCount = 0;
	for (int i = 0; i < treeVec.rows(); ++i)
		if (treeVec(i) >= 0)
			edgesCount++;
	 
	Eigen::MatrixXi edges(edgesCount, 2);
	int rowIdx = 0;
	for (int i = 0; i < treeVec.rows(); ++i)
	{
		if (treeVec(i) >= 0)
		{
			edges(rowIdx, 0) = treeVec(i);
			edges(rowIdx, 1) = i;
			rowIdx++;
		}
	}
 
	objWriteEdgesMat(pathName, edges, vers);
}

 


///////////////////////////////////////////////////////////////////////////////////////////////////////// �������ϵ��ؽӿ�
void objReadVerticesHomoMat(Eigen::MatrixXf& vers, const char* fileName);

void objWriteVerticesHomoMat(const char* fileName, const Eigen::MatrixXf& vers);

void vers2homoVers(Eigen::MatrixXf& homoVers, const Eigen::MatrixXf& vers);

Eigen::MatrixXf vers2homoVers(const Eigen::MatrixXf& vers);

void homoVers2vers(Eigen::MatrixXf& vers, const Eigen::MatrixXf& homoVers);

Eigen::MatrixXf homoVers2vers(const Eigen::MatrixXf& homoVers);





//////////////////////////////////////////////////////////////////////////////////////////// ��ѧ�������

// ϡ�����ת�ã�Eigen::SparseMatrix�Դ���transpose()����̫������
template<typename T>
bool spMatTranspose(Eigen::SparseMatrix<T>& smOut, const Eigen::SparseMatrix<T>& smIn) 
{
	smOut.resize(smIn.cols(), smIn.rows());
	std::vector<Eigen::Triplet<T>> trips;
	trips.reserve(smIn.nonZeros());
	traverseSparseMatrix(smIn, [&smIn, &trips](auto& iter) 
		{
			trips.push_back(Eigen::Triplet<T>{static_cast<int>(iter.col()), static_cast<int>(iter.row()), iter.value()});
		});
	smOut.setFromTriplets(trips.begin(), trips.end());

	return true;
}


// ��ǡ���ĳ������Է�����Ax == b;
template<typename T, int N>
bool solveLinearEquation(Eigen::Matrix<T, N, 1>& x, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, const Eigen::Matrix<T, N, 1>& b)
{
	// �����Է�����Ax == b;
	Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	x = svd.solve(b);

	return true;
}


// ��һϵ��ǡ���ĳ������Է�����AX == B;
template <typename T>
bool solveLinearEquations(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& X, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, \
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& B)
{
	if (A.rows() != B.rows())
		return false;

	Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	X.resize(A.cols(), B.cols());
	for (int i = 0; i < B.cols(); ++i)
	{
		Eigen::Matrix < T, Eigen::Dynamic, 1> x = svd.solve(B.col(i));
		X.col(i) = x;
	}

	return true;
}


// ���ɷ������ؾ����㷨�������ʽ��ֵ
template<typename T, int N>
double hornersPoly(const Eigen::Matrix<T, N, 1>& coeffs, const double x)
{
	// coeffs��{a0, a1, a2, ..., an}��ɵ�(n+1)������������ʽΪp == a0 + a1*x + a2*x^2 + ... + an* x^n; 
	int n = coeffs.rows() - 1;
	if (n < 0)
		return NAN;

	double result = coeffs(n);
	for (int i = n; i - 1 >= 0; i--)
	{
		result *= x;
		result += coeffs(i - 1);
	}

	return result;
}


// �����ʽ��һ��΢��
template<typename T, int N>
float polyDiff(const Eigen::Matrix<T, N, 1>& coeffs, const float x)
{
	// ����ʽp == a0 + a1*x + a2*x^2 + ... + an* x^n һ��΢��Ϊ��p' == a1 + 2*a2*x + 3*a3*x^2 ... n*an*x^(n-1)
	int coeffsCount = coeffs.rows() * coeffs.cols();
	if (coeffsCount <= 0)
		return NAN;

	if (coeffsCount == 1)
		return 1.0f;

	Eigen::VectorXf diffCoeffs(coeffsCount - 1);
	for (int i = 0; i < diffCoeffs.rows(); ++i)
		diffCoeffs(i) = (i + 1) * coeffs(i + 1);

	float result = hornersPoly(diffCoeffs, x);

	return result;
}


// ������ܾ���Ŀ����ڿ˻���Kronecker product��
template <typename T, typename Derived1, typename Derived2>
void kron(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& result, \
	const Eigen::MatrixBase<Derived1>& mm1, \
	const Eigen::MatrixBase<Derived2>& mm2)
{
	// A(m��n) tensor B(p��q) == C(m*p �� n*q);
	/*
		����C.block(i * m, j* n, p, q) == Aij * B;
		�磺
		1 2 3
		4 5 6

		tensor

		100 10
		1		0.1

		==
		100,	10,	200,		20,	 300,		 30,
		1,		0.1,		2,		0.2,		 3,		0.3,
		400,	 40,	500,		 50,	 600, 	60,
		4,		 0.4,		5,		0.5,		 6,		0.6,

	*/
	const Derived1& m1 = mm1.derived();
	const Derived2& m2 = mm2.derived();

	unsigned m = m1.rows();
	unsigned n = m1.cols();
	unsigned p = m2.rows();
	unsigned q = m2.cols();

	result.resize(0, 0);						// ���·����ڴ棻
	result.resize(m * p, n * q);
	Eigen::VectorXi rowIdxVec = Eigen::VectorXi::LinSpaced(0, m - 1, m - 1);
	auto calcKronCol = [&](const unsigned colIdx, const std::tuple<unsigned>& pt)
	{
		unsigned rowIdx0 = std::get<0>(pt);
		result.block(rowIdx0 * p, colIdx * q, p, q) = static_cast<T>(m1(rowIdx0, colIdx)) * m2.array().cast<T>();		// ת������������е��������ͣ�
	};

	auto calcKronRow = [&](const unsigned rowIdx)
	{
		std::tuple<unsigned> pt0 = std::make_tuple(rowIdx);
		PARALLEL_FOR(0, n, calcKronCol, pt0);
	};

	PARALLEL_FOR(0, m, calcKronRow);
}


template <typename Derived1, typename Derived2>
Eigen::MatrixXd kron(const Eigen::MatrixBase<Derived1>& mm1, \
	const Eigen::MatrixBase<Derived2>& mm2)
{
	Eigen::MatrixXd result;
	kron(result, mm1, mm2);
	return result;
}


// ����ʽ��ֵ
void polyInterpolation();


// ��˹��ֵ
void gaussInterpolation();


// ��С���˶���ʽ������ߣ�
void leastSquarePolyFitting();


// ��ع����ʽ�������
template <typename T>
void ridgeRegressionPolyFitting(Eigen::VectorXd& theta, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, unsigned m = 4)
{
	/*
		void ridgeRegressionPolyFitting(
					Eigen::VectorXd & theta,							��ϵĶ���ʽ����
					const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers			��ɢ������
					const unsigned m						��ϵĶ���ʽ������
					)
	*/

	const double lambda = 0.1;

	int versCount = vers.rows();
	if (versCount == 0)
		return;
	if (m >= versCount)
		m = versCount - 1;

	Eigen::MatrixXd X(versCount, m);
	for (int i = 0; i < versCount; ++i)
		for (int j = 0; j < m; ++j)
			X(i, j) = std::powf(static_cast<double>(vers(i, 0)), j);

	Eigen::VectorXd Y = vers.col(1).cast<double>();
	Eigen::MatrixXd Id(m, m);
	Id.setIdentity();
	theta = (X.transpose() * X + Id * lambda).inverse() * X.transpose() * Y;
}



// ��С���˷���ϣ��ƽ�����׼��Բ���������xy��������룩
template<typename T>
Eigen::VectorXd fittingStandardEllipse(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& sampleVers)
{
	/*
		Eigen::VectorXd fittingStandardEllipse(								// ����������(a,c,d,e,f)��Ϊ��׼��Բ���̵�ϵ����
					const Eigen::MatrixXf& sampleVers					// ����������㣬��������XOYƽ���ϵĵ㣻
		)
	*/
	const double epsilon = 1e-8;							// ����������ֵС�ڴ�ֵʱ��ΪΪ0��
	using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

	// ��׼��Բ���̣�a*x^2 + c*y^2 + d*x + e*y + f  = 0������a*c > 0
	Eigen::VectorXd x(Eigen::VectorXd::Zero(5));

	unsigned m = sampleVers.rows();				// sample count;
	Eigen::VectorXd x0(Eigen::VectorXd::Zero(m));
	Eigen::VectorXd y0(Eigen::VectorXd::Zero(m));
	for (unsigned i = 0; i < m; ++i)
	{
		x0(i) = static_cast<double>(sampleVers(i, 0));
		y0(i) = static_cast<double>(sampleVers(i, 1));
	}

	// alpha = [x^2, y^2, x, y, 1]; ������Ϣ����A = [alpha1; alpha2; .... alpham]; ��Բ����дΪ��A*x = 0;
	Eigen::MatrixXd A = Eigen::MatrixXd::Ones(m, 5);
	A.col(0) = x0.array() * x0.array();
	A.col(1) = y0.array() * y0.array();
	A.col(2) = x0;
	A.col(3) = y0;

	Eigen::MatrixXd ATA = A.transpose().eval() * A;
	Eigen::MatrixXd B(Eigen::MatrixXd::Zero(5, 5));
	B(0, 1) = 1;
	B(1, 0) = 1;
	Eigen::MatrixXd S = ATA.inverse() * B.transpose();

	// ��S������ֵ������������
	Eigen::EigenSolver<Eigen::MatrixXd> es(S);
	Eigen::MatrixXd D = es.pseudoEigenvalueMatrix();			// �Խ���Ԫ��������ֵ
	Eigen::MatrixXd V = es.pseudoEigenvectors();					// ÿһ����������������������

	// Ѱ������ֵ��Ϊ0��������Լ������������������
	for (unsigned i = 0; i < V.cols(); ++i)
	{
		double eigenValue = D(i, i);
		if (std::abs(eigenValue) < epsilon)
			continue;
		x = V.col(i);
		double a = x(0);
		double c = x(1);
		if (a * c > 0)
			break;
	}

	return x;
}





//////////////////////////////////////////////////////////////////////////////////////////////// ����������

// �õ�������������������
template <typename DerivedI>
bool getEdges(Eigen::MatrixXi& edges, 	const Eigen::PlainObjectBase<DerivedI>& tris)
{
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


template <typename DerivedI>
bool getUedges(Eigen::MatrixXi& uEdges, const Eigen::PlainObjectBase<DerivedI>& edges)
{
	unsigned edgesCount = edges.rows();
	uEdges.resize(0, 0);
	if (0 == edges.size() || 2 != edges.rows())
		return false;

	std::unordered_set<std::int64_t> uEdgeSet;
	for (unsigned i = 0; i < edgesCount; ++i)
		uEdgeSet.insert(encodeUedge(edges(i, 0), edges(i, 1)));
	unsigned uEdgesCount = uEdgeSet.size();
	uEdges.resize(uEdgesCount, 2);

	int index = 0;
	for (auto& code : uEdgeSet)
	{
		auto pair = decodeEdge(code);
		uEdges(index, 0) = pair.first;
		uEdges(index, 1) = pair.second;
		index++;
	}

	return true;
}


// ����������������Ƭ����
template <typename DerivedV, typename DerivedF, typename DerivedN>
bool getTriNorms(Eigen::PlainObjectBase<DerivedN>& triNorms, const Eigen::MatrixBase<DerivedV>& vers, \
	const Eigen::MatrixBase<DerivedF>& tris)
{
	triNorms.resize(tris.rows(), 3);
	int Frows = tris.rows();
	for (int i = 0; i < Frows; i++)
	{
		const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v1 = vers.row(tris(i, 1)) - vers.row(tris(i, 0));
		const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v2 = vers.row(tris(i, 2)) - vers.row(tris(i, 0));
		triNorms.row(i) = v1.cross(v2);        
		typename DerivedV::Scalar r = triNorms.row(i).norm();
		if (r == 0)
			return false;
		else
			triNorms.row(i) /= r;
	}

	return true;
}
 

// �õ�����������ʾ������ߣ�
template <typename DerivedV>
void getEdgeArrows(Eigen::PlainObjectBase<DerivedV>& arrows, const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers)
{
	DerivedV vas, vbs;
	Eigen::VectorXi vaIdxes = edges.col(0);
	Eigen::VectorXi vbIdxes = edges.col(1);
	subFromIdxVec(vas, vers, vaIdxes);
	subFromIdxVec(vbs, vers, vbIdxes);
	arrows = vbs - vas;
}


// ����������������Ƭ������
template<typename T, typename DerivedV, typename DerivedI>
bool trianglesBarycenter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& barys, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	const double eps = 1e-6;
	int versCount = vers.rows();
	int trisCount = tris.rows();
	if (versCount == 0 || trisCount == 0)
		return false;

	barys.resize(trisCount, 3);
	for (int i = 0; i < trisCount; ++i)
	{
		int vaIdx, vbIdx, vcIdx;
		vaIdx = tris(i, 0);
		vbIdx = tris(i, 1);
		vcIdx = tris(i, 2);
		Eigen::Matrix<T, 3, 3> tmpMat;
		tmpMat.row(0) = vers.row(vaIdx).array().cast<T>();
		tmpMat.row(1) = vers.row(vbIdx).array().cast<T>();
		tmpMat.row(2) = vers.row(vcIdx).array().cast<T>();
		barys.row(i) = tmpMat.colwise().mean();
	}

	return true;
}


// ����������������Ƭ�Ĺ�һ��������
template<typename T, typename DerivedV, typename DerivedI>
bool trianglesNorm(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& triNorms, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	const double eps = 1e-12;
	int versCount = vers.rows();
	int trisCount = tris.rows();
	if (versCount == 0 || trisCount == 0)
		return false;

	triNorms.resize(trisCount, 3);
	for (int i = 0; i < trisCount; ++i)
	{
		int vaIdx, vbIdx, vcIdx;
		Eigen::RowVector3d va, vb, vc, arrow1, arrow2, norm;
		vaIdx = tris(i, 0);
		vbIdx = tris(i, 1);
		vcIdx = tris(i, 2);
		va = vers.row(vaIdx).array().cast<double>();
		vb = vers.row(vbIdx).array().cast<double>();
		vc = vers.row(vcIdx).array().cast<double>();
		arrow1 = vb - va;
		arrow2 = vc - va;
		triNorms.row(i) = (arrow1.cross(arrow2)).array().cast<T>();
	}

	// ��������һ�����������˻�����Ƭ����дΪinf
	for (int i = 0; i < trisCount; ++i)
	{
		double length = triNorms.row(i).norm();
		if (abs(length) < eps)
		{
			triNorms.row(i) = Eigen::Matrix<T, 1, 3>{ INFINITY, INFINITY, INFINITY };
			continue;
		}
		Eigen::Matrix<T, 1, 3> tmpVec = triNorms.row(i) / length;
		triNorms.row(i) = tmpVec;
	}

	return true;
}


// ����������������Ƭ����ƽ��
template<typename DerivedV, typename DerivedI>
bool trianglesPlane(Eigen::MatrixXd& planeCoeff, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	int versCount = vers.rows();
	int trisCount = tris.rows();
	if (versCount == 0 || trisCount == 0)
		return false;

	Eigen::MatrixXd triNorms;
	if (!trianglesNorm(triNorms, vers, tris))
		return false;

	planeCoeff.resize(trisCount, 4);
	planeCoeff.leftCols(3) = triNorms;
	for (int i = 0; i < trisCount; ++i)
	{
		Eigen::RowVector3d norm = triNorms.row(i);

		// �������˻�����Ƭ����ƽ��ϵ����дΪNAN:
		if (std::isinf(norm(0)))
			planeCoeff(i, 3) = INFINITY;

		Eigen::RowVector3d va = vers.row(tris(i, 0)).array().cast<double>();

		// va����ƽ���� �� va�㵽ƽ��ľ���Ϊ0 �� norm.dot(va) +d == 0; �� d = -norm.dot(va);
		planeCoeff(i, 3) = -norm.dot(va);
	}

	return true;
}


// ������������ÿ������Ƭ�������
template<typename Scalar, typename Ta, typename DerivedI>
bool trisArea(Eigen::Matrix<Ta, Eigen::Dynamic, 1>& trisAreaVec, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const Eigen::MatrixBase<DerivedI>& tris)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> vas, vbs, vcs, arrows1, arrows2, arrows3;
	Eigen::VectorXd lens1, lens2, lens3, s;
	Eigen::VectorXi vaIdxes = tris.col(0);
	Eigen::VectorXi vbIdxes = tris.col(1);
	Eigen::VectorXi vcIdxes = tris.col(2);

	subFromIdxVec(vas, vers, vaIdxes);
	subFromIdxVec(vbs, vers, vbIdxes);
	subFromIdxVec(vcs, vers, vcIdxes);
	arrows1 = vbs - vas;
	arrows2 = vcs - vas;
	arrows3 = vbs - vcs;
	lens1 = arrows1.rowwise().norm().cast<double>();
	lens2 = arrows2.rowwise().norm().cast<double>();
	lens3 = arrows3.rowwise().norm().cast<double>();
	s = (lens1 + lens2 + lens3) / 2.0;				// ����Ƭ�ܳ���һ�룻

	// ʹ��Heron��ʽ��������Ƭ���:
	Eigen::VectorXd tmpVec = s.array() * (s - lens1).array() * (s - lens2).array() * (s - lens3).array();
	tmpVec = tmpVec.cwiseSqrt();
	trisAreaVec = tmpVec.array().cast<Ta>();

	return true;
}


template <typename Tv, typename DerivedF, typename DerivedL>
void edge_lengths(Eigen::PlainObjectBase<DerivedL>& lenMat, const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers,
	const Eigen::MatrixBase<DerivedF>& tris){
	const int trisCount = tris.rows();
	lenMat.resize(trisCount, 3);
	PARALLEL_FOR(0, trisCount, [&vers, &tris, &lenMat](const int i)
		{
			lenMat(i, 0) = (vers.row(tris(i, 1)) - vers.row(tris(i, 2))).norm();
			lenMat(i, 1) = (vers.row(tris(i, 2)) - vers.row(tris(i, 0))).norm();
			lenMat(i, 2) = (vers.row(tris(i, 0)) - vers.row(tris(i, 1))).norm();
		});
}


template <typename Tv, typename DerivedF, typename DerivedL>
void squared_edge_lengths(const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers,
	const Eigen::MatrixBase<DerivedF>& tris, Eigen::PlainObjectBase<DerivedL>& lenMat2)
{
	const int trisCount = tris.rows();

	lenMat2.resize(trisCount, 3);
	PARALLEL_FOR(0, trisCount, [&vers, &tris, &lenMat2](const int i)
		{
			lenMat2(i, 0) = (vers.row(tris(i, 1)) - vers.row(tris(i, 2))).squaredNorm();
			lenMat2(i, 1) = (vers.row(tris(i, 2)) - vers.row(tris(i, 0))).squaredNorm();
			lenMat2(i, 2) = (vers.row(tris(i, 0)) - vers.row(tris(i, 1))).squaredNorm();
		});
}
 

// ������������������
template<typename DerivedV, typename DerivedI>
double meshVolume(const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	double volume = -1.0;
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	if (vers.cols() != 3 || tris.cols() != 3)
		return volume;

	Eigen::VectorXi vaIdxes = tris.col(0);
	Eigen::VectorXi vbIdxes = tris.col(1);
	Eigen::VectorXi vcIdxes = tris.col(2);
	DerivedV versA, versB, versC;
	subFromIdxVec(versA, vers, vaIdxes);
	subFromIdxVec(versB, vers, vbIdxes);
	subFromIdxVec(versC, vers, vcIdxes);

	// һ������Ƭ��Ӧ��������ķ��������signed volume�� V(t0) == (-x3y2z1 + x2y3z1 + x3y1z2 - x1y3z2 - x2y1z3 + x1y2z3) / 6.0;
	auto volumeVec = \
		(-versC.col(0).array() * versB.col(1).array() * versA.col(2).array() + versB.col(0).array() * versC.col(1).array() * versA.col(2).array()\
			+ versC.col(0).array() * versA.col(1).array() * versB.col(2).array() - versA.col(0).array() * versC.col(1).array() * versB.col(2).array()\
			- versB.col(0).array() * versA.col(1).array() * versC.col(2).array() + versA.col(0).array() * versB.col(1).array() * versC.col(2).array()) / 6.0;

	volume = volumeVec.sum();

	return volume;
}


// ����������Ĳ�ͬȨ�ص��ڽӾ��� 
template<typename Index>
bool adjMatrix(Eigen::SparseMatrix<Index>& adjSM_eCount, \
	Eigen::SparseMatrix<Index>& adjSM_eIdx, const Eigen::MatrixXi& tris, const Eigen::MatrixXi& edges0 = Eigen::MatrixXi{})
{
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	std::vector<Eigen::Triplet<Index>> smElems, smElems_weighted;
	smElems.reserve(edgesCount);
	smElems_weighted.reserve(edgesCount);

	// 1. �󶥵��ڽӹ�ϵ��

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
	if (0 == edges0.size())
	{
		Eigen::MatrixXi edges = Eigen::MatrixXi::Zero(edgesCount, 2);
		Eigen::MatrixXi vaIdxes = tris.col(0);
		Eigen::MatrixXi vbIdxes = tris.col(1);
		Eigen::MatrixXi vcIdxes = tris.col(2);
		edges.block(0, 0, trisCount, 1) = vbIdxes;
		edges.block(trisCount, 0, trisCount, 1) = vcIdxes;
		edges.block(trisCount * 2, 0, trisCount, 1) = vaIdxes;
		edges.block(0, 1, trisCount, 1) = vcIdxes;
		edges.block(trisCount, 1, trisCount, 1) = vaIdxes;
		edges.block(trisCount * 2, 1, trisCount, 1) = vbIdxes;
		for (int i = 0; i < edgesCount; ++i)
		{
			smElems.push_back(Eigen::Triplet<Index>{edges(i, 0), edges(i, 1), 1});
			smElems_weighted.push_back(Eigen::Triplet<Index>{edges(i, 0), edges(i, 1), i});
		}
	}
	else
	{
		for (int i = 0; i < edgesCount; ++i)
		{
			smElems.push_back(Eigen::Triplet<Index>{edges0(i, 0), edges0(i, 1), 1});
			smElems_weighted.push_back(Eigen::Triplet<Index>{edges0(i, 0), edges0(i, 1), i});
		}
	}
 
	adjSM_eCount.resize(versCount, versCount);
	adjSM_eIdx.resize(versCount, versCount);
	adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());										// Ȩ��Ϊ��������ظ��Ĵ�����
	adjSM_eIdx.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());		// Ȩ��Ϊ������ߵ�������

	return true;
}



// ��Ե����ߣ�����1��
template<typename DerivedI>
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
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
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, std::vector<int>& bdryTriIdxes, const Eigen::PlainObjectBase<DerivedI>& tris)
{
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
		edgesMap.insert({ codeA, i});
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


// �����������еķ����ΰ�ߣ�����1
template<typename DerivedI>
bool nonManifoldEdges(Eigen::MatrixXi& nmnEdges, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
	
	unsigned nmnEdgesCount = 0;
	traverseSparseMatrix(adjSM_eCount, [&nmnEdgesCount](auto& iter)
		{
			if (iter.value() > 1)
				nmnEdgesCount++;
		});

	nmnEdges.resize(nmnEdgesCount, 2);
	int index = 0;
	traverseSparseMatrix(adjSM_eCount, [&index, &nmnEdges](auto& iter) 
		{
			if (iter.value() > 1)
			{
				nmnEdges(index, 0) = iter.row();
				nmnEdges(index, 1) = iter.col();
				index++;
			}
		});

	return true;
}


 
// �����������еķ����ΰ�ߣ�����2������������αߵ���ϸ��Ϣ��
template<typename DerivedI>
int nonManifoldEdges(Eigen::MatrixXi& nmnEdges, std::vector<std::pair<int, std::pair<int, int>>>& nmnEdgeInfos, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	adjMatrix(adjSM_eCount, adjSM_eIdx, tris);

	unsigned nmnEdgesCount = 0;
	traverseSparseMatrix(adjSM_eCount, [&nmnEdgesCount](auto& iter)
		{
			if (iter.value() > 1)
				nmnEdgesCount++;
		});

	nmnEdges.resize(nmnEdgesCount, 2);
	nmnEdgeInfos.reserve(nmnEdgesCount);
	int index = 0;
	traverseSparseMatrix(adjSM_eCount, [&index, &nmnEdges, &nmnEdgeInfos](auto& iter)
		{
			if (iter.value() > 1)
			{
				nmnEdges(index, 0) = iter.row();
				nmnEdges(index, 1) = iter.col();
				std::pair<int, int> edgePair = std::make_pair(iter.row(), iter.col());
				std::pair<int, std::pair<int, int>> infoPair = std::make_pair(static_cast<int>(iter.value()), edgePair);
				nmnEdgeInfos.push_back(infoPair);

				index++;
			}
		});

	return nmnEdgesCount;
}
 

// buildAdjacency()�����������������Ƭ�ڽ���Ϣ��
using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool buildAdjacency(Eigen::MatrixXi& ttAdj_nmEdge, std::vector<tVec>& ttAdj_nmnEdge, std::vector<tVec>& ttAdj_nmnOppEdge, \
		const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		bool buildAdjacency(
					const Eigen::MatrixXi& tris,									���������Ƭ����

					Eigen::MatrixXi& ttAdj_nmEdge,							����Ƭ�ķǱ�Ե����������ڽӵ�����Ƭ������
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


	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// ��������ϡ�����adjSM_eCount, adjSM_eIdx��triplet���ݣ�

	std::unordered_multimap<std::int64_t, int> edgeMap;			// �߱��롪��64λ����������ʾ�ı����ݣ������˵���������


	// 1. ������ı���Ϣ���ڽӾ���
	Eigen::SparseMatrix<int> adjSM;								// �ڽӾ���
	Eigen::SparseMatrix<int> adjSM_eCount;					// �ڽӾ�������Ϊ��������ظ��Ĵ�����
	Eigen::SparseMatrix<int> adjSM_eIdx;				// �ڽӾ���Ԫ������Ӧ���α���Ϊ�ñ�����������Ӧ�����α���Ϊ�ñ����������ĺͣ�
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;		// adjSM_eIdx��ת�ã���ʾ�Աߵ���Ϣ��
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
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris, edges);
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
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//	3.3 ������Ƭ�ڽӾ���ttAdj_mnEdge(i, :)������Ϊi������Ƭ�������ڽӵ���������Ƭ���������з����α߻��Ե����дΪ-1��
		ttAdj_nmEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
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


// triangleGrow()������ָ������Ƭ��ʼ�����������������񲻿����з����α�
template <typename DerivedV, typename DerivedI>
bool triangleGrow(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
		const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx)
{
	/*
		
			bool  triangleGrow(
						versOut,							// �������ĵ���
						trisOut,							// ������������Ƭ
						vers,								// ��������ĵ���
						tris,								// �������������Ƭ
						triIdx							// ����Ƭ�������������Ƭ������
						)

			�������������ڷ����αߣ����return false
	*/
	int versCount = vers.rows();
	int trisCount = tris.rows();
	int edgesCount = 3 * trisCount;

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// ��������ϡ�����adjSM_eCount, adjSM_eIdx��triplet���ݣ�	

	// 1. ������ı���Ϣ���ڽӾ���
	Eigen::SparseMatrix<int> adjSM;								// �ڽӾ���
	Eigen::SparseMatrix<int> adjSM_eCount;					// �ڽӾ�������Ϊ��������ظ��Ĵ�����
	Eigen::SparseMatrix<int> adjSM_eIdx;				// �ڽӾ���Ԫ������Ӧ���α���Ϊ�ñ�����������Ӧ�����α���Ϊ�ñ����������ĺͣ�
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;		// adjSM_eIdx��ת�ã���ʾ�Աߵ���Ϣ��
	Eigen::SparseMatrix<int> adjSM_ueCount;				// Ȩ��Ϊ�����ij����������Ƭ������
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;				// �Ǳ�Ե����������ڽӾ���
	Eigen::SparseMatrix<int> adjSM_MN_NB;					// �Ǳ�Ե����������ڽӾ���
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;			// �Ǳ�Ե��������߶Աߵ��ڽӾ���
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

		// 1.1 ��������������ڽӾ���adjSM, adjSM_eCount, smElems_weighted
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
		adjSM = adjSM_eCount;																		// ������ڽӾ���
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		��adjSM_eIdx(i, j)��Ӧ������������ߣ���Ȩ��ֵΪ�ñߵ����������Ƿ���������ߣ���ñ߶�Ӧ������������Ȩ��Ϊ��������֮�ͣ�
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

		//	1.2 ����������ڽӾ���
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// ������ڽӾ���
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.3 ���ɷǱ�Ե����������ڽӾ���
		smElems.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				if (2 == iter.value())
					smElems.push_back(Eigen::Triplet<int>(iter.row(), iter.col(), 1));
			});
		adjSM_ueMN_NB.resize(versCount, versCount);
		adjSM_ueMN_NB.setFromTriplets(smElems.begin(), smElems.end());
		smElems.clear();

		//	1.4 ���ɷǱ�Ե��������ߡ�����Աߵ��ڽӾ���
		adjSM_MN_NB = adjSM + adjSM_ueMN_NB;					// adjSM & adjSM_ueMN_NB
		adjSM_MN_NB.prune([](const Eigen::Index& row, const Eigen::Index& col, const int& value)->bool
			{
				if (2 == value)			// ɾ�������α�
					return true;
				else
					return false;
			});
		adjSM_MN_NB /= 2;
		adjSM_MN_NB_opp = adjSM_MN_NB.transpose();
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
	std::vector<int> etInfo;												// ������ - ����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
	Eigen::VectorXi etAdj_mnEdge;							// ����������Ƭ�������������α߻��Ե��д-1��
	Eigen::MatrixXi ttAdj_nmEdge;					// ����Ƭ�ķǱ�Ե����������ڽӵ�����Ƭ������	trisCount * 3(i, 0)Ϊ����Ϊi������Ƭ�ķǱ�Ե���������(vb, vc)�ڽӵ�����Ƭ��������
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
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//	3.3 ������Ƭ�ڽӾ���ttAdj_mnEdge(i, :)������Ϊi������Ƭ�������ڽӵ���������Ƭ���������з����α߻��Ե����дΪ-1��
		ttAdj_nmEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. ѭ�����������е�һ������Ƭt���� - tд��������� - ���t������������Ƭ������Ȼ�����
	std::unordered_set<int> finalTriIdx;						// ��ͨ����Ƭ���ϣ�
	std::unordered_set<int> workingSet;						// ����Ѱ����������Ƭ�Ķ��У�
	workingSet.insert(static_cast<int>(triIdx));							// ���в����һ������Ƭ��
	while (workingSet.size() > 0)
	{
		// 4.1 �׸�Ԫ�س��ӣ����뵽��ͨ����Ƭ�����У�
		int cTriIdx = *workingSet.begin();									// ��ǰ���������Ƭ����
		std::vector<int> adjTriIdx;												// ��ǰ����Ƭ��������Ƭ������
		workingSet.erase(workingSet.begin());
		finalTriIdx.insert(cTriIdx);

		// 4.2 ȷ����ǰ����Ƭ����������Ƭ
		for (int i = 0; i < 3; ++i)
			if (ttAdj_nmEdge(cTriIdx, i) >= 0)
				adjTriIdx.push_back(ttAdj_nmEdge(cTriIdx, i));

		// 4.3 ����ͨ����Ƭ������û�е�ǰ������Ƭ�����������Ƭ��������Ҳ���������Ƭ��
		for (const auto& index : adjTriIdx)
		{
			auto retPair = finalTriIdx.insert(index);
			if (retPair.second)
				workingSet.insert(index);
		}
	}


	// 5. ��ѡ��������Ƭȡ���㣺
	std::set<int> finalVerIdx;
	Eigen::VectorXi finalVerIdxVec;													// ����������������ԭ�����е�������
	Eigen::VectorXi oldNewIdxInfo = -Eigen::VectorXi::Ones(vers.rows());		// �������������±�Ϊ��������Ԫ��ֵΪ�������������������еĶ���дΪ-1
	for (const auto& index : finalTriIdx)
	{
		finalVerIdx.insert(tris(static_cast<int>(index), 0));
		finalVerIdx.insert(tris(static_cast<int>(index), 1));
		finalVerIdx.insert(tris(static_cast<int>(index), 2));
	}

	auto iter = finalVerIdx.begin();
	finalVerIdxVec.resize(finalVerIdx.size());
	for (int i = 0; i < finalVerIdx.size(); ++i)
		finalVerIdxVec(i) = static_cast<int>(*iter++);

	int setOrder = 0;
	for (const auto& index : finalVerIdx)
		oldNewIdxInfo[static_cast<int>(index)] = setOrder++;
	subFromIdxVec(versOut, vers, finalVerIdxVec);


	// 6. ԭ������������Ƭ�еĶ����������ϵĻ�Ϊ�µġ�
	Eigen::MatrixXi tempMat = tris;
	int* intPtr = tempMat.data();
	for (int i = 0; i < tempMat.size(); ++i)
	{
		int oldIdx = *intPtr;
		*intPtr = oldNewIdxInfo[oldIdx];
		intPtr++;
	}

	// 7. ȥ���Ƿ�����Ƭ��������-1Ԫ�ص�����Ƭ��
	trisOut.resize(tris.rows(), 3);
	int count = 0;
	for (int i = 0; i < tris.rows(); ++i)
	{
		Eigen::RowVector3i currentTri = tempMat.row(i);
		if ((currentTri.array() >= 0).all())
			trisOut.row(count++) = currentTri;
	}
	trisOut.conservativeResize(count, 3);

	return true;
}
 

// triangleGrowSplitMesh()�����������������������Ϊ�������ͨ�����������񲻿����з����α�
template <typename DerivedV>
bool triangleGrowSplitMesh(std::vector<DerivedV>& meshesVersOut, std::vector<Eigen::MatrixXi>& meshesTrisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	int versCount = vers.rows();
	int trisCount = tris.rows();
	int edgesCount = 3 * trisCount;

	// ��������е�triIdxΪ������ߵ������Ƭ������
	std::unordered_set<int> finalTriIdx;						// ��ͨ����Ƭ���ϣ�
	std::unordered_set<int> workingSet;						// ����Ѱ����������Ƭ�Ķ��У�

	Eigen::MatrixXi ttAdj_nmEdge;		// ����Ƭ�ķǱ�Ե����������ڽӵ�����Ƭ������	trisCount * 3(i, 0)Ϊ����Ϊi������Ƭ�ķǱ�Ե���������(vb, vc)�ڽӵ�����Ƭ��������

	Eigen::SparseMatrix<int> adjSM;								// �ڽӾ���
	Eigen::SparseMatrix<int> adjSM_eCount;					// �ڽӾ�������Ϊ��������ظ��Ĵ�����
	Eigen::SparseMatrix<int> adjSM_eIdx;				// �ڽӾ���Ԫ������Ӧ���α���Ϊ�ñ�����������Ӧ�����α���Ϊ�ñ����������ĺͣ�
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;		// adjSM_eIdx��ת�ã���ʾ�Աߵ���Ϣ��
	Eigen::SparseMatrix<int> adjSM_ueCount;				// Ȩ��Ϊ�����ij����������Ƭ������
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;				// �Ǳ�Ե����������ڽӾ���
	Eigen::SparseMatrix<int> adjSM_MN_NB;					// �Ǳ�Ե����������ڽӾ���
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;			// �Ǳ�Ե��������߶Աߵ��ڽӾ���

	Eigen::VectorXi etAdj_mnEdge;							// ����������Ƭ�������������α߻��Ե��д-1��

	std::vector<int> etInfo;												// ������ - ����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
	std::vector<int> edgesIdx_MN_NB;							// �Ǳ�Ե��������ߵ�������
	std::vector<int> edgesIdx_MN_NB_opp;					// �Ǳ�Ե��������ߵĶԱߵ�������
	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// ��������ϡ�����adjSM_eCount, adjSM_eIdx��triplet���ݣ�	

	std::vector<bool> trisFlag(trisCount, true);		// true��ʾ������Ƭ���ռ���false��ʾ���ռ���


	// 1. ������ı���Ϣ���ڽӾ���
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

		// 1.1 ��������������ڽӾ���adjSM, adjSM_eCount, smElems_weighted
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
		adjSM = adjSM_eCount;																		// ������ڽӾ���
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		��adjSM_eIdx(i, j)��Ӧ������������ߣ���Ȩ��ֵΪ�ñߵ����������Ƿ���������ߣ���ñ߶�Ӧ������������Ȩ��Ϊ��������֮�ͣ�
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

		//	1.2 ����������ڽӾ���
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// ������ڽӾ���
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.3 ���ɷǱ�Ե����������ڽӾ���
		smElems.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				if (2 == iter.value())
					smElems.push_back(Eigen::Triplet<int>(iter.row(), iter.col(), 1));
			});
		adjSM_ueMN_NB.resize(versCount, versCount);
		adjSM_ueMN_NB.setFromTriplets(smElems.begin(), smElems.end());
		smElems.clear();

		//	1.4 ���ɷǱ�Ե��������ߡ�����Աߵ��ڽӾ���
		adjSM_MN_NB = adjSM + adjSM_ueMN_NB;					// adjSM & adjSM_ueMN_NB
		adjSM_MN_NB.prune([](const Eigen::Index& row, const Eigen::Index& col, const int& value)->bool
			{
				if (2 == value)			// ɾ�������α�
					return true;
				else
					return false;
			});
		adjSM_MN_NB /= 2;
		adjSM_MN_NB_opp = adjSM_MN_NB.transpose();
	}


	// 2. ȷ���Ǳ�Ե��������ߡ�����Աߵ�������
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
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//	3.3 ������Ƭ�ڽӾ���ttAdj_mnEdge(i, :)������Ϊi������Ƭ�������ڽӵ���������Ƭ���������з����α߻��Ե����дΪ-1��
		ttAdj_nmEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. ѭ�����������е�һ������Ƭt���� - tд��������� - ���t������������Ƭ������Ȼ�����
	std::vector<std::unordered_set<int>> connectedTriIdxes;
	int collectedTrisCount = 0;
	while (collectedTrisCount < trisCount) 
	{
		finalTriIdx.clear();

		// ȡ��ǰ�׸�δ�ռ�������Ƭ��Ϊ��������Ƭ��
		int seedTriIdx = 0;
		for (int i = 0; i < trisCount; ++i)
			if (trisFlag[i])
				seedTriIdx = i;

		workingSet.insert(static_cast<int>(seedTriIdx));					// ���в�����������Ƭ��
		while (workingSet.size() > 0)
		{
			// 4.1 �׸�Ԫ�س��ӣ����뵽��ͨ����Ƭ�����У�
			int cTriIdx = *workingSet.begin();									// ��ǰ���������Ƭ����
			trisFlag[cTriIdx] = false;

			std::vector<int> adjTriIdx;												// ��ǰ����Ƭ��������Ƭ������
			workingSet.erase(workingSet.begin());
			finalTriIdx.insert(cTriIdx);

			// 4.2 ȷ����ǰ����Ƭ����������Ƭ
			for (int i = 0; i < 3; ++i)
				if (ttAdj_nmEdge(cTriIdx, i) >= 0)
					adjTriIdx.push_back(ttAdj_nmEdge(cTriIdx, i));

			// 4.3 ����ͨ����Ƭ������û�е�ǰ������Ƭ�����������Ƭ��������Ҳ���������Ƭ��
			for (const auto& index : adjTriIdx)
			{
				auto retPair = finalTriIdx.insert(index);
				if (retPair.second && trisFlag[index])
				{
					workingSet.insert(index);
					trisFlag[index] = false;
				}
			}
		}

		connectedTriIdxes.push_back(finalTriIdx);
		collectedTrisCount += finalTriIdx.size();
	}


	// 5. ������е���ͨ����
	int meshesCount = connectedTriIdxes.size();
	meshesVersOut.resize(meshesCount);
	meshesTrisOut.resize(meshesCount);
	for (int i = 0; i < meshesCount; ++i) 
	{
		DerivedV& cMeshVers = meshesVersOut[i];
		Eigen::MatrixXi& cMeshTris = meshesTrisOut[i];

		// 5.1 ��connectedTriIdxes����ͨ������Ƭ����ȡ���㣺
		std::set<int> cMeshVerIdx;																	// ��ǰ�ĵ���ͨ�����а����Ķ����������
		for (const auto& index : connectedTriIdxes[i])
		{
			cMeshVerIdx.insert(tris(static_cast<int>(index), 0));
			cMeshVerIdx.insert(tris(static_cast<int>(index), 1));
			cMeshVerIdx.insert(tris(static_cast<int>(index), 2));
		}

		Eigen::VectorXi cMeshVerIdxVec(cMeshVerIdx.size());									// ����������������ԭ�����е�������
		Eigen::VectorXi oldNewIdxInfo = -Eigen::VectorXi::Ones(vers.rows());			// �������������±�Ϊ��������Ԫ��ֵΪ�������������������еĶ���дΪ-1
		auto iter = cMeshVerIdx.begin();
		for (int i = 0; i < cMeshVerIdx.size(); ++i)
			cMeshVerIdxVec(i) = static_cast<int>(*iter++);

		int setOrder = 0;
		for (const auto& index : cMeshVerIdx)
			oldNewIdxInfo[static_cast<int>(index)] = setOrder++;
		subFromIdxVec(cMeshVers, vers, cMeshVerIdxVec);

		// 5.2 ȷ����ǰ�ĵ���ͨ���������Ƭ��
		std::vector<int> cTriIdxes;
		cTriIdxes.insert(cTriIdxes.end(), connectedTriIdxes[i].begin(), connectedTriIdxes[i].end());
		subFromIdxVec(cMeshTris, tris, cTriIdxes);
		int* verIdxPtr = cMeshTris.data();
		for (int k = 0; k < 3 * cMeshTris.rows(); ++k)				// ����Ƭ�е��ϵĶ���������Ϊ�µģ� 
		{
			int oldIdx = *verIdxPtr;
			*verIdxPtr = oldNewIdxInfo[oldIdx];
			verIdxPtr++;
		}
	}

	return true;
}


// triangleGrowOuterSurf()������ָ������Ƭ���������ⲿ������Ƭ����ʼ������������ȡ�����浥��ͨ��������
template <typename T>
bool triangleGrowOuterSurf(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, const bool blRemIso = true)
{
	// ������Ŀǰ�ٶ������������������Ƭ������ȷ����û���˻�����Ƭ��

	int versCount = vers.rows();
	int trisCount = tris.rows();
	int edgesCount = 3 * trisCount;
	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// ��������ϡ�����adjSM_eCount, adjSM_eIdx��triplet���ݣ�	

	// 0. ������������Ƭ����
	Eigen::MatrixXd triNorms;
	if (!getTriNorms(triNorms, vers, tris))				// �������˻�����Ƭ���ᵼ�¼�����Ƭ�������
		return false;

	// 1. ������Ƭ�ڽӹ�ϵ��
	Eigen::MatrixXi ttAdj_nmEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency(ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return false;

	// 2. ѡȡһ�������������αߵ�����Ƭ��Ϊ��ɢ����������Ƭ��
	int seedIdx = -1;
	for (int i = 0; i<trisCount; ++i) 
	{
		if (ttAdj_nmEdge(i, 0) >= 0 && ttAdj_nmEdge(i, 1) >= 0 && ttAdj_nmEdge(i, 2) >= 0)
		{
			seedIdx = i;
			break;
		}
	}
	if (seedIdx < 0)
		return false;

	// 3. ѭ�����������е�һ������Ƭt���� - tд��������� - ���t������������Ƭ������Ȼ�����
#ifdef LOCAL_DEBUG
	std::unordered_set<int> deleTriIdxes, seleTriIdxes;
#endif
	std::unordered_set<int> finalTriIdxes;							 
	std::unordered_set<int> workingSet;								// ����Ѱ����������Ƭ�Ķ��У�
	std::vector<int> triTags(trisCount, 0);								// 0 : δ���ʣ�1: ����¼�� -1: ��ɾ����

	// ������ע��tagд1�Ժ����п��ܱ���дΪ-1�ģ�Ŀǰ���ǵ�һ�η��ʵ�ʱ�ı�ǩ�������Ƿ񱣴�����Ƭ��
	workingSet.insert(static_cast<int>(seedIdx));						// ���в����һ������Ƭ��
	triTags[seedIdx] = 1;
	while (workingSet.size() > 0)
	{
		// w.1 �׸�Ԫ�س��ӣ����뵽��ͨ����Ƭ�����У�
		int cTriIdx = *workingSet.begin();									// ��ǰ���������Ƭ����
		workingSet.erase(workingSet.begin());		
		if (triTags[cTriIdx] < 0)											// �п������ʱ��δ����״̬�����֮�󱻱��Ϊɾ������Ƭ��
			continue;
		finalTriIdxes.insert(cTriIdx);

		// w.2 ȷ����ǰ����Ƭ���ڵ���������Ƭ���������������ɾ��
		std::vector<int> adjTriIdx;												// ��ǰ����Ƭ��������Ƭ������
		for (int i = 0; i < 3; ++i)
		{
			int nmAdjTriIdx = ttAdj_nmEdge(cTriIdx, i);
			if (nmAdjTriIdx >= 0)
					adjTriIdx.push_back(nmAdjTriIdx);						// wf.a. ���α߹�������������Ƭ��������
			else
			{
				// wf.b.1. ɾ����ǰ�����α߹�������������Ƭ��
				auto& vec1 = ttAdj_nmnEdge[cTriIdx][i];				// ��ǰ�����α����ڵ�����Ƭ��������ǰ����Ƭ����
				for (const auto& index : vec1)
				{
					if (index == cTriIdx)
						continue;
					else
					{
						triTags[index] = -1;
#ifdef LOCAL_DEBUG
						deleTriIdxes.insert(index);
#endif
					}
				}

				// wf.b.2. ��ǰ�����αߵĶԱ����ڵ�����Ƭ�У�ѡȡ����н����(��dotֵ��С)���Ǹ���
				auto& vec2 = ttAdj_nmnOppEdge[cTriIdx][i];			// ��ǰ�����αߵĶԱ����ڵ���������Ƭ
				if (vec2.empty())
					continue;
				else if (1 == vec2.size())
						adjTriIdx.push_back(vec2[0]);
				else
				{
					unsigned nopTrisCount = vec2.size();
					Eigen::RowVector3d cTriNorm = triNorms.row(cTriIdx);
					Eigen::MatrixXd nopTriNorms;
					Eigen::VectorXd dotValues(nopTrisCount);
					subFromIdxVec(nopTriNorms, triNorms, vec2);
					for (unsigned k = 0; k < nopTrisCount; ++k)
						dotValues(k) = cTriNorm.dot(nopTriNorms.row(k));
					Eigen::Index minIdx;
					dotValues.minCoeff(&minIdx);
					int selTriIdx = vec2[minIdx];
					adjTriIdx.push_back(selTriIdx);

					for (const auto& index : vec2)				// δ��ѡȡ�Ķ�����������Ƭȫ��ɾ�� ��
					{
						if (selTriIdx != index)
						{
							triTags[index] = -1;
#ifdef LOCAL_DEBUG
							deleTriIdxes.insert(index);
#endif
						}
					}
				}

			}
		}

		// w.3 finalTriIdxes�в��뱣����������������Ƭ����������Ƭ֮ǰδ���ʣ��������У�
		for (const auto& index : adjTriIdx)
		{
			if (0 == triTags[index])
			{
				auto retPair = finalTriIdxes.insert(index);
				if (retPair.second)
				{
					workingSet.insert(index);
					triTags[index] = 1;
				}
			}
		}
	}

	// for debug:
#ifdef LOCAL_DEBUG
	Eigen::MatrixXi deleTris, seleTris;
	Eigen::MatrixXd versCopy = vers.array().cast<double>();
	subFromIdxCon(deleTris, tris, deleTriIdxes);
	objWriteMeshMat("E:/deleTris.obj", versCopy, deleTris);
	subFromIdxCon(seleTris, tris, seleTriIdxes);
	objWriteMeshMat("E:/seleTris.obj", versCopy, seleTris);
#endif

	// 4. ��ȡѡ�е�����Ƭ
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vers1 = vers;
	Eigen::MatrixXi tris1;
	subFromIdxCon(tris1, tris, finalTriIdxes); 

	// 5. ȥ���������㣻
	if (blRemIso)
	{
		std::vector<unsigned> isoVerIdxes = checkIsoVers(vers1, tris1);
		if (isoVerIdxes.empty())
		{
			versOut = vers1;
			trisOut = tris1;
		}
		else
			removeIsoVers(versOut, trisOut, vers1, tris1, isoVerIdxes);
	}
	else
	{
		versOut = vers1;
		trisOut = tris1;
	} 

	return true;
}


// correctTriDirs()�����������������Ƭ���򣬻�������Ƭ���������ķ�����ʹ��ǰ��Ҫ�ȴ��������е��ص�����Ƭ��
template <typename DerivedV, typename DerivedI>
int correctTriDirs(Eigen::PlainObjectBase<DerivedI>& trisOut, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx, const double thetaThreshold = 0.99 * pi)
{
	int corrCount = 0;
	const double dotThreshold = cos(thetaThreshold);
	const int versCount = vers.rows();
	const int trisCount = tris.rows();
	const int edgesCount = 3 * trisCount;
	trisOut = tris;

#ifdef  LOCAL_DEBUG
	std::vector<int> corrTriIdxes;
#endif 
 
	// 0. ���ʼ״̬��������������Ƭ�ķ���
	Eigen::MatrixXd triNorms;
	if (!getTriNorms(triNorms, vers, tris))				// �������˻�����Ƭ���ᵼ�¼�����Ƭ�������
		return false;

	// 1. ������Ƭ�ڽӹ�ϵ��
	Eigen::MatrixXi ttAdj_nmEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency(ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return false;

	// 2. ѭ�����������е�һ������Ƭt���� - tд��������� - ���t������������Ƭ������Ȼ�����
	std::unordered_set<int> workingSet;								// ����Ѱ����������Ƭ�Ķ��У�
	std::vector<int> triTags(trisCount, 0);								// 0 : δ���ʣ�1: �ѷ���
	workingSet.insert(static_cast<int>(triIdx));						// ���в����һ������Ƭ�� 
	while (workingSet.size() > 0)
	{
		// 4.1 �׸�Ԫ�س��ӣ����뵽��ͨ����Ƭ�����У�
		int cTriIdx = *workingSet.begin();									// ��ǰ���������Ƭ����
		std::vector<int> adjTriIdxes;												// ��ǰ����Ƭ��δ�����ʵ���������Ƭ������
		Eigen::RowVector3d cTriNorm = triNorms.row(cTriIdx);
		Eigen::MatrixXd adjTriNorms;
		unsigned adjTrisCount = 0;
		workingSet.erase(workingSet.begin());
		triTags[cTriIdx] = 1;											// ��ǰ����Ƭ���Ϊ�ѷ��ʣ�

		// 4.2 ȷ����ǰ����Ƭ��������������Ƭ
		for (int i = 0; i < 3; ++i)
		{
			int nbrTriIdx = ttAdj_nmEdge(cTriIdx, i);
			if (nbrTriIdx >= 0)
			{
				// a. ��ǰ��Ϊ���α�
				if (0 == triTags[nbrTriIdx])
					adjTriIdxes.push_back(nbrTriIdx);
			}
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
				//		��������Ƭ����
				int tmp = trisOut(adjIdx, 1);
				trisOut(adjIdx, 1) = trisOut(adjIdx, 2);
				trisOut(adjIdx, 2) = tmp;
				triNorms.row(adjIdx).array() *= -1;				//	���·���
				corrCount++;											// �������ӣ�
#ifdef LOCAL_DEBUG
				corrTriIdxes.push_back(adjIdx);
#endif
			}			
		}

		// 4.3 ���������в��뵱ǰ����Ƭδ�����ʵ���������Ƭ��
		for (const auto& index : adjTriIdxes)
				workingSet.insert(index);
	}

#ifdef LOCAL_DEBUG
	Eigen::MatrixXi sickTris;
	Eigen::MatrixXd versCopy = vers.array().cast<double>();
	subFromIdxVec(sickTris, tris, corrTriIdxes);
	objWriteMeshMat("E:/sickTris.obj", versCopy, sickTris);
#endif

	return corrCount;
}


// findOverlapTris()����Ѱ�������е��ص�����Ƭ����������Ƭ���������ķ�����
template <typename DerivedV, typename DerivedI>
int findOverLapTris(std::vector<std::pair<int, int>>& opTrisPairs, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx = 0, const double thetaThreshold = 0.99 * pi)
{
	int olCount = 0;
	const double dotThreshold1 = cos(thetaThreshold);					
	const double dotThreshold2 = cos(pi - thetaThreshold);
	const int versCount = vers.rows();
	const int trisCount = tris.rows();
	const int edgesCount = 3 * trisCount;
	std::unordered_set<std::int64_t> opTriPairCodes;			// һ���ص�������Ƭ������ֵ��ǰС������У������std::int64_t

#ifdef  LOCAL_DEBUG
	std::set<int> olTriIdxes;
#endif 

	// 0. ���ʼ״̬��������������Ƭ�ķ���
	Eigen::MatrixXd triNorms;
	if (!getTriNorms(triNorms, vers, tris))				// �������˻�����Ƭ���ᵼ�¼�����Ƭ�������
		return false;

	// 1. ������Ƭ�ڽӹ�ϵ��
	Eigen::MatrixXi ttAdj_nmEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency(ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return false;

	// 2. ����Ƭ��������-Ѱ�����ڵ��ص�����Ƭ�Ե�ѭ����
	std::unordered_set<int> workingSet;								// ����Ѱ����������Ƭ�Ķ��У�
	std::vector<int> triTags(trisCount, 0);								// 0 : δ���ʣ�1: �ѷ���
	workingSet.insert(static_cast<int>(triIdx));						// ���в����һ������Ƭ�� 
	while (workingSet.size() > 0)
	{
		// 2.1 �׸�Ԫ�س��ӣ����뵽��ͨ����Ƭ�����У�
		int cTriIdx = *workingSet.begin();									// ��ǰ���������Ƭ����
		std::vector<int> adjTriIdxes;												// ��ǰ����Ƭ��δ�����ʵ���������Ƭ������
		Eigen::RowVector3d cTriNorm = triNorms.row(cTriIdx);
		unsigned adjTrisCount = 0;
		workingSet.erase(workingSet.begin());
		triTags[cTriIdx] = 1;											// ��ǰ����Ƭ���Ϊ�ѷ��ʣ�

		// 2.2 ȷ����ǰ����Ƭ��������������Ƭ
		for (int i = 0; i < 3; ++i)
		{
			int nbrTriIdx = ttAdj_nmEdge(cTriIdx, i);
			if (nbrTriIdx >= 0)
			{
				// a. ��ǰ��Ϊ���α�
				if (0 == triTags[nbrTriIdx])
					adjTriIdxes.push_back(nbrTriIdx);
			}
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

		// 2.3 ����-�ж���ǰ����Ƭ���ܱ�����Ƭ�Ƿ��ص���
		adjTrisCount = adjTriIdxes.size();
		if (0 == adjTrisCount)
			continue;
		for (unsigned k = 0; k < adjTrisCount; ++k)
		{
			int adjIdx = adjTriIdxes[k];
			double dotValue = cTriNorm.dot(triNorms.row(adjIdx));
			if (dotValue > dotThreshold2 || dotValue < dotThreshold1)			
			{
				// ����ǰ����Ƭ����������Ƭ����ӽ���ͬ���෴�������������ƬͶӰ���䷨��ƽ�����Ƿ����ص����� 
				int A, B, C1, C2;				// AB�������������Ƭ����Ķ��㣬C1, C2�����Ǹ�������Ķ��㣻
				std::vector<int> ctv{ tris(cTriIdx, 0), tris(cTriIdx, 1), tris(cTriIdx, 2) };
				std::vector<int> adjtv{ tris(adjIdx, 0), tris(adjIdx, 1), tris(adjIdx, 2) };
				auto iterC = ctv.begin();
				for (; iterC != ctv.end(); iterC++) 
				{
					int cIdx = *iterC;
					if (cIdx != adjtv[0] && cIdx != adjtv[1] && cIdx != adjtv[2])
					{
						C1 = cIdx;
						break;
					}
				}
				for (const auto& aIdx : adjtv)
				{
					if (aIdx != ctv[0] && aIdx != ctv[1] && aIdx != ctv[2])
					{
						C2 = aIdx;
						break;
					}
				}
				ctv.erase(iterC);
				A = ctv[0];					
				B = ctv[1];

				// �ж�C1, C2�Ƿ���ͬ�ࣺ
				Eigen::RowVector3d AB = vers.row(B).array().cast<double>() - vers.row(A).array().cast<double>();
				Eigen::RowVector3d AC1 = vers.row(C1).array().cast<double>() - vers.row(A).array().cast<double>();
				Eigen::RowVector3d AC2 = vers.row(C2).array().cast<double>() - vers.row(A).array().cast<double>();
				Eigen::RowVector3d crossArrow1 = AB.cross(AC1);
				Eigen::RowVector3d crossArrow2 = AB.cross(AC2);
				if (crossArrow1.dot(cTriNorm) * crossArrow2.dot(cTriNorm) > 0)				// C1, C2�ڹ�����AB��ͬһ�ࣻ
				{
					auto retPair = opTriPairCodes.insert(encodeUedge(cTriIdx, adjIdx));
					if (retPair.second)
					{
						opTrisPairs.push_back(std::make_pair(cTriIdx, adjIdx));
						olCount++;								// �������ӣ�
#ifdef LOCAL_DEBUG
						olTriIdxes.insert(cTriIdx);
						olTriIdxes.insert(adjIdx);
#endif
					}
				}							
			}
		}

		// 2.3 ���������в��뵱ǰ����Ƭδ�����ʵ���������Ƭ��
		for (const auto& index : adjTriIdxes)
			workingSet.insert(index);
	}
 
#ifdef LOCAL_DEBUG
	if (olTriIdxes.size() > 0) 
	{
		Eigen::MatrixXd versCopy = vers.array().cast<double>();
		Eigen::MatrixXi olTris;
		subFromIdxCon(olTris, tris, olTriIdxes);
		objWriteMeshMat("E:/overlapTris.obj", versCopy, olTris);
	}
#endif 

	return olCount;
}


// ����ߵ�������������ȡ��·�ߣ�����ı߾�������������Ҳ����������ߣ�
template <typename DerivedI>
bool uEdgeGrow(std::vector<Eigen::MatrixXi>& circs, std::vector<Eigen::MatrixXi >& segs, \
	const Eigen::PlainObjectBase<DerivedI>& uEdges)
{
	std::list<std::list<std::pair<int, int>>> circList, segList;
	Eigen::SparseMatrix<int> adjSM;										// �ڽӾ���
	std::vector<Eigen::Triplet<int>> smElems;

	// 1. ����������ڽӾ���
	int maxIdx = 0;
	const int* idxPtr = uEdges.data();
	for (int i = 0; i < uEdges.size(); i++, idxPtr++)
		if (*idxPtr > maxIdx)
			maxIdx = *idxPtr;
	smElems.reserve(2 * uEdges.rows());
	for (unsigned i = 0; i < uEdges.rows(); ++i)
	{
		smElems.push_back(Eigen::Triplet<int>{uEdges(i, 0), uEdges(i, 1), 1});
		smElems.push_back(Eigen::Triplet<int>{uEdges(i, 1), uEdges(i, 0), 1});
	}
	adjSM.resize(maxIdx + 1, maxIdx + 1);
	adjSM.setFromTriplets(smElems.begin(), smElems.end());
	smElems.clear();
	traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
		{
			iter.valueRef() = 1;
		});

	// 2. ͳ�Ʋ��ظ��������������
	std::unordered_set<std::int64_t> tmpSet;
	for (unsigned i = 0; i < uEdges.rows(); ++i)
		tmpSet.insert(encodeUedge(uEdges(i, 0), uEdges(i, 1)));
	int remainCount = tmpSet.size();			// (3,4) (4,3)��ʾͬһ������ߣ�
	
	// 3. �����������ѭ����
	while (remainCount > 0)
	{
		std::list<std::pair<int, int>> currentSeg;
	
		// w1. ��ȡ�ڽӱ��е�һ������ߣ�
		bool breakFlag = false;
		for (unsigned i = 0; !breakFlag & i < adjSM.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, i); !breakFlag & iter; ++iter)
			{
				if (iter.valueRef() > 0)
				{
					currentSeg.push_back({ iter.row(), iter.col() });
					iter.valueRef() = -1;												// ��ʾɾ����Ԫ�أ�
					adjSM.coeffRef(iter.col(), iter.row()) = -1;
					remainCount--;
					breakFlag = true;
				}
			}
		}

		// w2. ������һ������߹����ıߣ�ֱ�����������������γɻ�·Ϊֹ��
		int head = currentSeg.front().first;
		int tail = currentSeg.front().second;
		while (1)
		{
			int otherEnd = -1;
			breakFlag = false;
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, head); !breakFlag & iter; ++iter)	// �Ե�head�еı�����
			{
				if (iter.valueRef() > 0)
				{
					otherEnd = iter.row();
					currentSeg.push_back({ head, otherEnd });
					iter.valueRef() = -1;												// ��ʾɾ����Ԫ�أ�
					adjSM.coeffRef(iter.col(), iter.row()) = -1;
					remainCount--;
					head = otherEnd;
					breakFlag = true;
				}
			}

			if (otherEnd == tail)
			{
				circList.push_back(currentSeg);
				break;
			}

			if (otherEnd < 0)
			{
				segList.push_back(currentSeg);
				break;
			}
		}
	}

	// 4. �����
	circs.reserve(circList.size());
	segs.reserve(segList.size());
	for (const auto& list : circList)
	{
		unsigned pairCount = list.size();
		Eigen::MatrixXi circEdges(pairCount, 2);
		auto iter = list.begin();
		for (unsigned i = 0; i < pairCount; ++i, iter++)
		{
			circEdges(i, 0) = iter->first;
			circEdges(i, 1) = iter->second;
		}
		circs.push_back(circEdges);
	}

	for (const auto& list : segList)
	{
		unsigned pairCount = list.size();
		Eigen::MatrixXi segEdges(pairCount, 2);
		auto iter = list.begin();
		for (unsigned i = 0; i < pairCount; ++i, iter++)
		{
			segEdges(i, 0) = iter->first;
			segEdges(i, 1) = iter->second;
		}
		segs.push_back(segEdges);
	}

	return true;
}


template <typename DerivedV, typename DerivedI>
int findHolesBdrySegs(std::vector<std::vector<int>>& holes, std::vector<std::vector<int>>& bdrySegs, \
		const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		int findHolesBdrySegs(															�ɹ�����(����+��Ե������)��ʧ�ܷ���-1;
			std::vector<std::vector<int>>& holes,								�������պ�Ϊ��·�ı�Ե���㣬ÿȦ���㰴˳ʱ�����У�
			std::vector<std::vector<int>>& bdrySegs,							��Ե���ߣ������պϵı�Ե����
			const Eigen::PlainObjectBase<DerivedV>& vers, 
			const Eigen::PlainObjectBase<DerivedI>& tris
			)
			
	*/

	// ȷ�����ıߡ���ֻ����һ������Ƭ�ıߣ�
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	std::vector<Eigen::MatrixXi> circs, segs;
	Eigen::MatrixXi bdrys;
	bdryEdges(bdrys, tris);
	uEdgeGrow(circs, segs, bdrys);

#ifdef LOCAL_DEBUG
	Eigen::MatrixXd versCopy = vers.array().cast<double>();
	if (circs.size() > 0)
	{
		Eigen::MatrixXi tmpMat;
		for (const auto& mat : circs)
			matInsertRows(tmpMat, mat);
		objWriteEdgesMat("E:/circEdges.obj", tmpMat, versCopy);
	}
	if (segs.size() > 0)
	{
		Eigen::MatrixXi tmpMat;
		for (const auto& mat : segs)
			matInsertRows(tmpMat, mat);
		objWriteEdgesMat("E:/segEdges.obj", tmpMat, versCopy);
	}
#endif


	return holes.size() + bdrySegs.size();
}



// ����������Ƿ��зǷ�����Ƭ������������������������ͬ�ģ�
template <typename DerivedI>
bool checkSickTris(std::vector<unsigned>& sickIdxes, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	bool retFlag = false;
	for (unsigned i = 0; i < tris.rows(); ++i)
	{
		if (tris(i, 0) == tris(i, 1) || tris(i, 0) == tris(i, 2) || tris(i, 1) == tris(i, 2))
		{
			sickIdxes.push_back(i);
			retFlag = true;
		}
	}
	return retFlag;
}


// ȷ�����������е����е���ͨ���򣨶�����ͨ��
template <typename Index>
int simplyConnectedRegion(const Eigen::SparseMatrix<Index>& adjSM,
	Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount)
{
	/*
		int simplyConnectedRegion(												���ص���ͨ���������������ʱ����-1
					const Eigen::SparseMatrix<Index>& adjSM,			����������ڽӾ���
					Eigen::VectorXi& connectedLabels,						ͬһ��ͨ�����ڵĶ����ǩ������Ϊ0, 1, 2,...��
																											����Ϊi�Ķ���ı�ǩΪconnectedLabels(i);
					Eigen::VectorXi& connectedCount							ÿ����ǩ��Ӧ�ĵ���ͨ�����ڰ����Ķ�����
																											��ǩΪi�ĵ���ͨ��������Ķ�����ΪconnectedCount(i)
					)
	*/
	if (adjSM.rows() == 0)
		return -1;

	if (adjSM.rows() != adjSM.cols())
		return -1;

	const unsigned versCount = adjSM.rows();

	// 1. ��ʼ����
	connectedLabels.setConstant(versCount, versCount);		// ���ս���������ܵı�ǩΪversCount-1;
	connectedCount.setZero(versCount);

	// 2. �������ж��㣬����������ͨ�Ķ��㣺
	int currentLabel = 0;
	for (unsigned i = 0; i < versCount; ++i) 
	{
		// 2.1 ����ǰ�����ѱ����ʹ�(label < versCount)����continue:
		if (connectedLabels(i) < versCount)
			continue;

		// 2.2 ����ǰ����δ�����ʹ���ִ��BFS�ռ���������ͨ�Ķ��㣺
		std::queue<Index> workingQueue;
		workingQueue.push(i);
		while (!workingQueue.empty())
		{
			const Index curVerIdx = workingQueue.front();
			workingQueue.pop();
			
			if (connectedLabels(curVerIdx) < versCount)
				continue;

			// 2.2.1 ���׶���label��ֵ��label����+1
			connectedLabels(curVerIdx) = currentLabel;
			connectedCount(currentLabel)++;

			// 2.2.2 ���ڽӾ�����������ǰ�����ڽӵĶ��㣺
			for (typename Eigen::SparseMatrix<Index>::InnerIterator iter(adjSM, curVerIdx); iter; ++iter)
			{
				const Index connectVerIdx = iter.row();					// Ĭ��ϡ����������ȴ洢����ǰ�������ڵ�curVerIdx�У�

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


// ȷ�����������е����е���ͨ��������Ƭ��ͨ��������triangleGrow��ǿ������������Դ��з�����Ԫ�أ�
template <typename DerivedI>
int simplyTrisConnectedRegion(	Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	// 
	const unsigned trisCount = tris.rows();
	connectedLabels = Eigen::VectorXi::Zero(trisCount);
	unsigned visitedCount = 0;

	Eigen::MatrixXi ttAdj_nmEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency(ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return -1;

	// �ڽ�����Ƭ������ѭ����
	int currentLabel = 1;
	while (visitedCount < trisCount)
	{
		int seedIdx = 0;
		for (unsigned i = 0; i < trisCount; ++i)
		{
			if (0 == connectedLabels(i))
			{
				seedIdx = i;
				break;
			}
		}

		std::unordered_set<int> workingQueue;
		workingQueue.insert(seedIdx);

		while (!workingQueue.empty())
		{
			int cIdx = *workingQueue.begin();
			visitedCount++;
			workingQueue.erase(cIdx);
			connectedLabels(cIdx) = currentLabel;

			for (unsigned i = 0; i < 3; ++i)
				if (ttAdj_nmEdge(cIdx, i) >= 0)
					if (0 == connectedLabels(ttAdj_nmEdge(cIdx, i)))
						workingQueue.insert(ttAdj_nmEdge(cIdx, i));					// ֻ����δ���ʵ�����Ƭ

			for (const auto& vec : ttAdj_nmnEdge[cIdx])
				for (const auto& index : vec)
					if (0 == connectedLabels(index))
						workingQueue.insert(index);
 
		}

		currentLabel++;
	}

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
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();

	// 1. �����ڽӾ���
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx, adjSM;
	adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
	adjSM = adjSM_eCount;
	traverseSparseMatrix(adjSM, [&adjSM](auto& iter) 
		{
			if (iter.value() > 0)
				iter.valueRef() = 1;
		});

	// 2. ȷ�����������е���ͨ����
	Eigen::VectorXi connectedLabels, connectedCount;
	int conCount = simplyConnectedRegion(adjSM, connectedLabels, connectedCount);
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


// ȥ������Ƭ��
template <typename IndexT>
bool removeTris(Eigen::MatrixXi& trisOut, const Eigen::MatrixXi& tris, const std::vector<IndexT>& sickTriIdxes)
{
	const unsigned trisCount = tris.rows();
	if (sickTriIdxes.size() == 0)
		return false;

	Eigen::VectorXi tmpVec = Eigen::VectorXi::Ones(trisCount);
	for (const auto& index : sickTriIdxes)
	{
		if (index < IndexT(0) || index >= trisCount)
			return false;
		tmpVec(index) = -1;
	}

	std::vector<unsigned> selectedTriIdxes;
	selectedTriIdxes.reserve(trisCount);
	for (unsigned i = 0; i < trisCount; ++i)
		if (tmpVec(i) > 0)
			selectedTriIdxes.push_back(i);
	selectedTriIdxes.shrink_to_fit();

	trisOut.resize(0, 0);
	subFromIdxVec(trisOut, tris, selectedTriIdxes);

	return true;
}


// ��������еĹ�������
template <typename DerivedV>
std::vector<unsigned> checkIsoVers(const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
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

	// for debug:
	if (0) 
	{
		Eigen::MatrixXd isoVers;
		if (!isoVerIdxes.empty())
		{
			std::cout << "isoVerIdxes.size() == " << isoVerIdxes.size() << std::endl;
			subFromIdxVec(isoVers, vers, isoVerIdxes);
			objWriteVerticesMat("E:/isoVers.obj", isoVers);
		}
	}
 
	return isoVerIdxes;
}


// ȥ�������еĹ������㣺
template <typename DerivedV>
bool removeIsoVers(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::MatrixXi& trisOut, \
		const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const std::vector<unsigned>& isoVerIdxes)
{
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


// ����˻��ߣ��߳����̵ıߣ�
template <typename DerivedV>
Eigen::VectorXi checkDegEdges(const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
		const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const double threshold = 1e-3)
{
	/*
		Eigen::VectorXi checkDegEdges(														// ���ر���˻��ߵ�flag����
				const Eigen::MatrixXi& edges,												// ���������
				const Eigen::PlainObjectBase<DerivedV>& edgeArrows,		// ����ߵı���������
				const Eigen::PlainObjectBase<DerivedV>& vers, 
				const Eigen::MatrixXi& tris, 
				const double threshold = 1e-3												// �ж��˻��ߵı߳���ֵ��
				)
	*/
	Eigen::VectorXd edgesLen = edgeArrows.rowwise().norm();
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	unsigned edgesCount = edges.rows();

	Eigen::VectorXi flags = (edgesLen.array() <= threshold).select(Eigen::VectorXi::Ones(edgesCount), Eigen::VectorXi::Zero(edgesCount));
	
	return flags;
}


// ȥ�������˻��ߣ��߳����̵ıߣ������������ںϣ�;
template <typename DerivedV>
int mergeDegEdges(Eigen::PlainObjectBase<DerivedV>& newVers, Eigen::MatrixXi& newTris, \
	const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const Eigen::VectorXi& degEdgeFlags)
{
	/*
		int mergeDegEdges(																				����ȥ�����Ķ���������ʧ���򷵻�-1��
				Eigen::PlainObjectBase<DerivedV>& newVers, 
				Eigen::MatrixXi& newTris, 
				const Eigen::MatrixXi& edges, 
				const Eigen::PlainObjectBase<DerivedV>& edgeArrows, 
				const Eigen::PlainObjectBase<DerivedV>& vers, 
				const Eigen::MatrixXi& tris,	
				const Eigen::VectorXi& degEdgeFlags								����˻��ߵ�flag����
				)
	
	*/
	int repVersCount = -1;
	int* dataPtr = nullptr;
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	unsigned edgesCount = edges.rows();

	// 1. ��ȡ�˻������
	int degCount = degEdgeFlags.sum();				// �˻�����ߵ�����
	Eigen::MatrixXi degEdges;
	subFromFlagVec(degEdges, edges, degEdgeFlags);

	// for debug:
	if (1)
	{
		Eigen::MatrixXd degEdgeVers;
		std::unordered_set<int> tmpSet;
		std::vector<int> degVerIdxes;

		objWriteEdgesMat("E:/degEdges.obj", degEdges, vers);
		for (unsigned i = 0; i < degCount; ++i)
		{
			tmpSet.insert(degEdges(i, 0));
			tmpSet.insert(degEdges(i, 1));
		}
		degVerIdxes.insert(degVerIdxes.end(), tmpSet.begin(), tmpSet.end());
		subFromIdxVec(degEdgeVers, vers, degVerIdxes);
		objWriteVerticesMat("E:/degEdgeVers.obj", degEdgeVers);
	}

	// �����˻����У�����������С�Ķ��㣬������Ķ����дΪ����С�Ķ��㣺
	Eigen::MatrixXi revDegEdges = degEdges;
	for (unsigned i = 0; i < degCount; ++i)
	{
		int vaIdx = degEdges(i, 0);
		int vbIdx = degEdges(i, 1);
		if (degEdges(i, 0) < degEdges(i, 1))
			revDegEdges.row(i) = Eigen::RowVector2i(vbIdx, vaIdx);
		else
			degEdges.row(i) = Eigen::RowVector2i(vbIdx, vaIdx);
	}

	std::list<std::set<int>> clusters;
	for (unsigned i = 0; i < degCount; ++i)
	{
		int vaIdx = degEdges(i, 0);
		int vbIdx = degEdges(i, 1);
		for (auto& clusterSet : clusters)
		{
			auto retIter1 = clusterSet.find(vaIdx);
			auto retIter2 = clusterSet.find(vbIdx);
			if (clusterSet.end() == retIter1 && clusterSet.end() == retIter2)
				continue;
			else       // ��ǰ�ߵĶ�����֮ǰ�Ĵ��еĶ�����ͬ�����뵽�ô��У�
			{
				clusterSet.insert(vaIdx);
				clusterSet.insert(vbIdx);
				break;
			}
		}

		// �����µĴأ�
		std::set<int> tmpCluster;
		tmpCluster.insert(vaIdx);
		tmpCluster.insert(vbIdx);
		clusters.push_back(tmpCluster);
	}

	// ���ɶ�������ӳ���ֵ䡪��һ�����������������㶼ӳ��Ϊ������С���Ǹ����㣺
	std::map<int, int> clusterMap;
	for (const auto clusterSet: clusters) 
	{
		auto iter = clusterSet.begin();
		int headIdx = *iter;
		iter++;
		for (; iter != clusterSet.end(); ++iter)
			clusterMap.insert({*iter, headIdx});
	}

	std::vector<int> repIdxes;							// ��Ҫȥ���Ķ����������
	repIdxes.reserve(degCount);
	for (const auto& pair : clusterMap)
		repIdxes.push_back(pair.second);
	repIdxes.shrink_to_fit();
	repVersCount = repIdxes.size();
 	
	std::map<int, int> fullMap = clusterMap;
	for (int i = 0; i < versCount; ++i)
		fullMap.insert({i, i});

	// ����ӳ�䣺
	Eigen::MatrixXi trisCopy = tris;
	dataPtr = trisCopy.data();
	for (unsigned i = 0;  i < 3 * trisCount; ++i ) 
	{
		int index = *dataPtr;
		*dataPtr = fullMap[index];
		dataPtr++;
	}

	// ɾ���Ƿ�����Ƭ��
	std::vector<unsigned> sickIdxes;
	checkSickTris(sickIdxes, trisCopy);
	removeTris(newTris, trisCopy, sickIdxes);
 
	// ɾ���������㣺
	trisCopy = newTris;
	std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, trisCopy);
	newVers.resize(0, 0);
	newTris.resize(0, 0);
	removeIsoVers(newVers, newTris, vers, trisCopy, isoVerIdxes);

	return degCount;
}


// ����դ�����
template<typename T>
bool genGrids(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, const Eigen::Matrix<T, 1, 3>& origin, \
		const float step, const std::vector<unsigned>& gridCounts)
{
	/*
		bool genGrids(
			Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& gridCenters,			// �����դ�����
			const Eigen::Matrix<T, 1, 3>& origin,														// դ��ԭ�㣬������������С���Ǹ�դ��㣻
			const float step,																						// ��������
			const std::vector<unsigned>& gridCounts											// ��Ԫ���飬�ֱ���XYZ�����ϵĲ�����
			)	
	*/

	using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

	// �������볡�İ�Χ�У�
	RowVector3T minp = origin;
	RowVector3T maxp = origin + step * RowVector3T(gridCounts[0] - 1, gridCounts[1] - 1, gridCounts[2] - 1);

	// ����դ��
	/*
		�������������е�դ�����ĵ�Ϊ��
		gc(000), gc(100), gc(200), gc(300),...... gc(010), gc(110), gc(210), gc(310),...... gc(001), gc(101), gc(201).....

		x���꣺
		x0, x1, x2, x3......x0, x1, x2, x3......x0, x1, x2, x3......
		����ΪxCount;
		�ظ�����Ϊ(yCount * zCount)

		y���꣺
		y0, y0, y0...y1, y1, y1...y2, y2, y2.........y0, y0, y0...y1, y1, y1...
		����Ϊ(xCount * yCount);
		�ظ�����ΪzCount;
		����Ԫ���ظ�����ΪxCount

		z���꣺
		z0, z0, z0......z1, z1, z1......z2, z2, z2......
		����Ԫ���ظ�����Ϊ(xCount * yCount)
	*/
	VectorXT xPeriod = VectorXT::LinSpaced(gridCounts[0], minp(0), maxp(0));
	VectorXT yPeriod = VectorXT::LinSpaced(gridCounts[1], minp(1), maxp(1));
	VectorXT zPeriod = VectorXT::LinSpaced(gridCounts[2], minp(2), maxp(2));

	MatrixXT tmpVec0, tmpVec1, tmpVec2, tmpVec11;
	kron(tmpVec0, VectorXT::Ones(gridCounts[1] * gridCounts[2]), xPeriod);
	kron(tmpVec11, yPeriod, VectorXT::Ones(gridCounts[0]));
	kron(tmpVec1, VectorXT::Ones(gridCounts[2]), tmpVec11);
	kron(tmpVec2, zPeriod, VectorXT::Ones(gridCounts[0] * gridCounts[1]));
	gridCenters.resize(gridCounts[0] * gridCounts[1] * gridCounts[2], 3);
	gridCenters.col(0) = tmpVec0;
	gridCenters.col(1) = tmpVec1;
	gridCenters.col(2) = tmpVec2;

	return true;
}


// ȥ���Ƿ�����Ƭ���ظ�����Ƭ��
template<typename DerivedV, typename DerivedI>
int removeSickDupTris(const Eigen::PlainObjectBase<DerivedV>& vers, Eigen::PlainObjectBase<DerivedI>& tris)
{
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();

	Eigen::VectorXi flags{Eigen::VectorXi::Ones(trisCount)};			// flag��������Ҫȥ��������Ƭ��0��
	std::unordered_set<std::uint64_t> triCodeSet;

	// ��������Ƭ
	for (unsigned i = 0; i < trisCount; ++i) 
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
		if(!retPair.second)
			flags(i) = 0;				// �ظ�����Ƭ�Ƿ���
	}

	DerivedI tmpMat;
	subFromFlagVec(tmpMat, tris, flags);
	tris = tmpMat;

	int removeCount = trisCount - tris.rows();
	return removeCount;
}


// ��������������Ƭ�����ǵ�����ֵ��
template <typename Tv, typename Tl>
void trisCotValues(const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers,
	const Eigen::MatrixXi& tris, Eigen::Matrix<Tl, Eigen::Dynamic, Eigen::Dynamic>& cotValues)
{
	// 1/2*cotangents corresponding angles. for triangles, columns correspond to edges bc, ca, ab;
	const int trisCount = tris.rows();

	// ����ÿ������Ƭ�������������
	Eigen::Matrix<Tl, Eigen::Dynamic, 3> lenMat2;
	Eigen::Matrix<Tl, Eigen::Dynamic, 3> lenMat;
	Eigen::Matrix<Tl, Eigen::Dynamic, 1> dblA;
	squared_edge_lengths(vers, tris, lenMat2);
	lenMat = lenMat2.array().sqrt();
	trisArea(dblA, vers, tris);
	dblA.array() *=  2;
	cotValues.resize(trisCount, 3);
	for (int i = 0; i < trisCount; i++)
	{
		cotValues(i, 0) = (lenMat2(i, 1) + lenMat2(i, 2) - lenMat2(i, 0)) / dblA(i) / 4.0;
		cotValues(i, 1) = (lenMat2(i, 2) + lenMat2(i, 0) - lenMat2(i, 1)) / dblA(i) / 4.0;
		cotValues(i, 2) = (lenMat2(i, 0) + lenMat2(i, 1) - lenMat2(i, 2)) / dblA(i) / 4.0;
	} 
}


// ��������Ȩ��laplacian (ͬʱҲ�ǸնȾ���(stiffness matrix))�����������ģ�Polygon Laplacian Made Simple [Bunge et al. 2020]
template<typename Tv, typename Tl>
bool cotLaplacian(Eigen::SparseMatrix<Tl>& L, const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();
	L.resize(versCount, versCount);
	L.reserve(10 * versCount);

	Eigen::MatrixXi edges(3, 2);
	edges << 1, 2, 2, 0, 0, 1;
 
	// Gather cotangents
	Eigen::Matrix<Tl, Eigen::Dynamic, Eigen::Dynamic> C;
	trisCotValues(vers, tris, C);

	std::vector<Eigen::Triplet<Tl> > IJV;
	IJV.reserve(tris.rows() * edges.rows() * 4);

	// Loop over triangles
	for (int i = 0; i < tris.rows(); i++)
	{
		// loop over edges of element
		for (int e = 0; e < edges.rows(); e++)
		{
			int source = tris(i, edges(e, 0));
			int dest = tris(i, edges(e, 1));
			IJV.push_back(Eigen::Triplet<Tl>(source, dest, C(i, e)));
			IJV.push_back(Eigen::Triplet<Tl>(dest, source, C(i, e)));
			IJV.push_back(Eigen::Triplet<Tl>(source, source, -C(i, e)));
			IJV.push_back(Eigen::Triplet<Tl>(dest, dest, -C(i, e)));
		}
	}
	L.setFromTriplets(IJV.begin(), IJV.end());

	return true;
}


// ������������
template <typename Tm, typename Tv>
bool massMatrix_baryCentric(Eigen::SparseMatrix<Tm>& MM, const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const Eigen::MatrixXi& tris)
{
	const int versCount = vers.rows();
	const int trisCount = tris.rows();
	Eigen::MatrixXd lenMat;
	edge_lengths(lenMat, vers, tris);
	
	Eigen::VectorXd dblA;
	trisArea(dblA, vers, tris);
	dblA.array() *= 2;
	Eigen::MatrixXi MI;
	Eigen::MatrixXi MJ;
	Eigen::MatrixXd MV;

	// diagonal entries for each face corner
	MI.resize(trisCount * 3, 1); 
	MJ.resize(trisCount * 3, 1); 
	MV.resize(trisCount * 3, 1);
	MI.block(0 * trisCount, 0, trisCount, 1) = tris.col(0);
	MI.block(1 * trisCount, 0, trisCount, 1) = tris.col(1);
	MI.block(2 * trisCount, 0, trisCount, 1) = tris.col(2);
	MJ = MI;
	repmat(MV, dblA, 3, 1);
	MV.array() /= 6.0;
	sparse(MM, MI, MJ, MV, versCount, versCount);

	return true;
}


// laplace��˳ 
template <typename T>
bool laplaceFaring(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, const float deltaLB, const unsigned loopCount)
{
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// 1. ����laplace����
	Eigen::SparseMatrix<T> Lmat;
	cotLaplacian(Lmat, vers, tris);

	// 2. ��˳��ѭ����
	MatrixXT versCopy = vers;
	for (unsigned i = 0; i < loopCount; ++i)
	{
		// f1. ���㵱ǰ��������
		Eigen::SparseMatrix<T> mass;
		massMatrix_baryCentric(mass, vers, tris);

		// f2. �����Է����� (mass - delta*L) * newVers = mass * newVers
		const auto& S = (mass - deltaLB * Lmat);
		Eigen::SimplicialLLT<Eigen::SparseMatrix<T>> solver(S);
		assert(solver.info() == Eigen::Success);
		versOut = solver.solve(mass * versCopy).eval();
		versCopy = versOut;
	}

	return true;
}

 
//////////////////////////////////////////////////////////////////////////////////////////////// ͼ������أ�

// ����Ŀ��������˲�����������ע���˲�mask�ߴ����Ϊ����������
template<typename T>
bool linearSpatialFilter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matOut, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matIn, \
	const Eigen::MatrixXd& mask)
{
	unsigned rows = matIn.rows();
	unsigned cols = matIn.cols();
	matOut.resize(rows, cols);
	matOut.setZero(rows, cols);

	unsigned sizeMask = mask.rows();
	unsigned im = (sizeMask + 1) / 2;					// ��Ĥ���ĵ��±ꣻ
	unsigned km = im;
	unsigned offset = (sizeMask - 1) / 2;

	// 1. �����������б�Ե���أ�������չ����
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matExt = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows + 2*offset, cols + 2 * offset);
	matExt.block(offset, offset, rows, cols) = matIn;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rowHead = matIn.block(0, 0, offset, cols);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rowTail = matIn.block(rows - offset, 0, offset, cols);

	matExt.block(0, offset, offset, cols) = rowHead;
	matExt.block(rows + offset, offset, offset, cols) = rowTail;

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> colHead = matExt.block(0, offset, rows + 2*offset, offset);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> colTail = matExt.block(0, cols, rows + 2 * offset, offset);
	matExt.block(0, 0, rows + 2*offset, offset) = colHead;
	matExt.block(0, cols + offset, rows + 2 * offset, offset) = colTail;

	// 2. �������������
	PARALLEL_FOR(0, rows, [&](int i)
		{
			for (unsigned k = 0; k < cols; ++k)
			{
				unsigned ie = i + offset;					// ��ʱmask����Ԫ����matExt�е�λ���±ꣻ
				unsigned ke = k + offset;
				Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> coveredElems = matExt.block(ie - offset, ke - offset, sizeMask, sizeMask);
				Eigen::MatrixXd tmpMat = coveredElems.array().cast<double>().array() * mask.array();
				matOut(i, k) = static_cast<T>(tmpMat.sum());
			}
		});
 

	return true;
}
 
 

///////////////////////////////////////////////////////////////////////////////////////// debug �ӿڣ�

// �Զ����ʱ����ʹ��WINDOWS��ʱAPI
class tiktok
{
private:
	tiktok() = default;
	tiktok(const tiktok&) {}
	~tiktok() = default;

public:
	DWORD startTik;
	DWORD endTik;

	static tiktok& getInstance()
	{
		static tiktok tt_instance;
		return tt_instance;
	}

	void start()
	{
		this->startTik = GetTickCount();
	}

	void endCout(const char* str)
	{
		this->endTik = GetTickCount();
		std::cout << str << endTik - startTik << std::endl;
	}


	bool endWrite(const char* fileName, const char* str)
	{
		this->endTik = GetTickCount();
		std::ofstream file(fileName, std::ios_base::out | std::ios_base::app);
		if (!file)
		{
			return false;
		}

		file << str << endTik - startTik << std::endl;
		file.close();
		return true;
	}
};

 

namespace TEST_MYEIGEN
{
	const double pi = 3.14159;
	void test0();
	void test1();
	void test11();
	void test111();
	void test1111();
	void test2();
	void test3();
	void test4();
	void test5();
	void test6();
	void test7();
	void test8();
	void test9();
	void test10();

	void test12();
	void test13();
	void test14();
	void test15();
	void test16();
}


// ����ͼ����
namespace TEST_DIP 
{
	void test0();

}




