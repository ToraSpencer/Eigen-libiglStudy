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
 
//#define USE_TRIANGLE_H

#ifdef USE_TRIANGLE_H
// 和algorithm工程一样，使用单精度；libigl库中封装的三角剖分使用的是双精度；
#define ANSI_DECLARATORS
#define REAL DOUBLE
#define VOID int
#include "triangulate.h"
#endif

#include "Eigen/Dense"
#include "Eigen/Sparse"


/////////////////////////////////////////////////////////////////////////////////////////////////// 自己写的各个静态库：
#include "myEigenIO/myEigenIO.h"
#pragma comment(lib,"myEigenIO.lib")	

#include "myEigenBasicMath/myEigenBasicMath.h"
#pragma comment(lib, "myEigenBasicMath.lib")

#include "myEigenModeling/myEigenModeling.h"
#pragma comment(lib, "myEigenModeling.lib")

#include "myEigenPMP/myEigenPMP.h"
#pragma comment(lib, "myEigenPMP.lib")


/////////////////////////////////////////////////////////////////////////////////////////////////// debug全局变量
static Eigen::MatrixXd g_debugVers;


/////////////////////////////////////////////////////////////////////////////////////////////////// general tools: 
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


/////////////////////////////////////////////////////////////////////////////////////////////////// 辅助结构：
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


/////////////////////////////////////////////////////////////////////////////////////////////////// debug接口：
namespace MY_DEBUG
{
	// 顶点或三角片矩阵转换为triplet向量的形式；
	template <typename T>
	std::vector<triplet<T>> mat2triplets(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat);

	template <typename T>
	std::vector<doublet<T>> mat2doublets(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat);

	template <typename T>
	triplet<T> vec2triplet(const Eigen::Matrix<T, 1, 3>& vec);

	template <typename T>
	doublet<T> vec2doublet(const Eigen::Matrix<T, 1, 2>& vec);


	// 遍历搜索triplet向量，若索引为index的triplet元素使得谓词f返回值为true，则返回index; 若找不到或出错则返回-1；
	template <typename T, typename F>
	int findTriplet(const std::vector<triplet<T>>& trips, F f);

	template <typename T, typename F>
	int findTriplet(const std::vector<doublet<T>>& doubs, F f);

	template <typename T1, typename T2>
	void dispPair(const std::pair<T1, T2>& pair);
	template<typename T>
	void dispQuat(const Eigen::Quaternion<T>& q);
	template <typename Derived>
	void dispMat(const Eigen::PlainObjectBase<Derived>& mat);
	template <typename Derived>
	void dispMatBlock(const Eigen::PlainObjectBase<Derived>& mat, const int rowStart, const int rowEnd, const int colStart, const int colEnd);

	template<typename spMat>				// 如果模板参数为T, 函数参数为Eigen::SparseMatrix<T>，编译会报错，不知道什么缘故；
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

	// 自定义计时器，使用WINDOWS计时API
	class tiktok
	{
	private:
		tiktok() = default;
		tiktok(const tiktok&) {}
		~tiktok() = default;

	public:
		DWORD startTik;
		DWORD endTik;
		unsigned recordCount;
		std::vector<DWORD> records;

		static tiktok& getInstance()
		{
			static tiktok tt_instance;
			return tt_instance;
		}

		void start()
		{
			this->startTik = GetTickCount();
			this->recordCount = 0;
			this->records.clear();
		}

		void endCout(const char* str)
		{
			this->endTik = GetTickCount();
			std::cout << str << endTik - startTik << std::endl;
		}

		DWORD endGetCount()
		{
			this->endTik = GetTickCount();
			return endTik - startTik;
		}

		bool endWrite(const char* fileName, const char* str)
		{
			this->endTik = GetTickCount();
			std::ofstream file(fileName, std::ios_base::out | std::ios_base::app);
			if (!file)
				return false;

			file << str << endTik - startTik << std::endl;
			file.close();
			return true;
		}

		void takeArecord()
		{
			this->records.push_back(GetTickCount());
			recordCount++;
		}
	};
}
using namespace MY_DEBUG;



/////////////////////////////////////////////////////////////////////////////////////////////////// 暂时不知道如何分类：
template<typename Tv, typename Tl>
bool cotLaplacian(Eigen::SparseMatrix<Tl>& L, const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris);


// 生成质量矩阵：
template <typename Tm, typename Tv>
bool massMatrix_baryCentric(Eigen::SparseMatrix<Tm>& MM, const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const Eigen::MatrixXi& tris);


// laplace光顺 
template <typename T>
bool laplaceFaring(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, \
	const float deltaLB, const unsigned loopCount);


template<typename T>
bool linearSpatialFilter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matIn, const Eigen::MatrixXd& mask);

 

/////////////////////////////////////////////////////////////////////////////////////////////////// 网格缺陷检查和修复
template <typename DerivedI>
void findRepTris(std::vector<int>& repIdxes, const Eigen::PlainObjectBase<DerivedI>& tris);

template <typename DerivedV, typename DerivedI>
int findOverLapTris(std::vector<std::pair<int, int>>& opTrisPairs, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx, const double thetaThreshold);

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
int nonManifoldVers(std::vector<int>& nmnVerIdxes, const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris);


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






/////////////////////////////////////////////////////////////////////////////////////////////////// 区域生长算法：

template <typename DerivedV, typename DerivedI>
bool triangleGrow(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx);


template <typename DerivedV>
bool triangleGrowSplitMesh(std::vector<DerivedV>& meshesVersOut, std::vector<Eigen::MatrixXi>& meshesTrisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris);


template <typename T>
bool triangleGrowOuterSurf(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, const bool blRemIso, const int startIdx);


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
bool simplyConnectedLargest(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris);


template <typename T>
bool simplyConnectedSplitMesh(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& compVers,\
	std::vector<Eigen::MatrixXi>& compTris, 	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const Eigen::MatrixXi& tris);



/////////////////////////////////////////////////////////////////////////////////////////////////// 网格属性：

// 输入网格三角片数据，得到有向边数据：
template <typename DerivedI>
bool getEdges(Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedI>& tris);


template <typename DerivedI>
bool getEdgeIdxMap(std::unordered_multimap<std::int64_t, int>& map, const Eigen::PlainObjectBase<DerivedI>& edges);


template <typename DerivedV, typename DerivedF, typename DerivedN>
bool getTriNorms(Eigen::PlainObjectBase<DerivedN>& triNorms, const Eigen::MatrixBase<DerivedV>& vers, \
	const Eigen::MatrixBase<DerivedF>& tris);


template <typename DerivedV>
void getEdgeArrows(Eigen::PlainObjectBase<DerivedV>& arrows, const Eigen::MatrixXi& edges, \
	const Eigen::PlainObjectBase<DerivedV>& vers);


template <typename Tv, typename DerivedF, typename DerivedL>
void squared_edge_lengths(const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers,
	const Eigen::MatrixBase<DerivedF>& tris, Eigen::PlainObjectBase<DerivedL>& lenMat2);


template <typename Tv, typename Tl>
void trisCotValues(const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers,
	const Eigen::MatrixXi& tris, Eigen::Matrix<Tl, Eigen::Dynamic, Eigen::Dynamic>& cotValues);


template<typename T, typename DerivedV, typename DerivedI>
bool trianglesBarycenter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& barys, \
	const Eigen::PlainObjectBase<DerivedV>& vers, 	const Eigen::PlainObjectBase<DerivedI>& tris);


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



using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool buildAdjacency(Eigen::MatrixXi& ttAdj_mnEdge, std::vector<tVec>& ttAdj_nmnEdge, std::vector<tVec>& ttAdj_nmnOppEdge, \
	const Eigen::PlainObjectBase<DerivedI>& tris);


template <typename DerivedI>
bool buildAdjacency_new(Eigen::MatrixXi& ttAdj_mnEdge, std::vector<tVec>& ttAdj_nmnEdge, std::vector<tVec>& ttAdj_nmnOppEdge, \
	const Eigen::PlainObjectBase<DerivedI>& tris);




// class——方向包围盒类OBB（本质上是中点在原点的AABB加上一个仿射变换）；
template <typename _Scalar>
class OBB : public Eigen::AlignedBox<_Scalar, 3>
{
	// 成员数据：注！！！继承自Eigen::AlignedBox中的m_min, m_max是列向量，下面两个是行向量；
public:
	Eigen::Matrix<_Scalar, 1, 3> m_dir;
	Eigen::Matrix<_Scalar, 1, 3> m_center;
	Eigen::Matrix3d m_rotation;

	// constructor:
public:
	inline OBB() :m_dir(0, 0, 1), m_center(Eigen::Matrix<_Scalar, 1, 3>::Zero(3)), m_rotation(Eigen::Matrix3d::Identity()) {}


	inline OBB(const Eigen::AlignedBox<_Scalar, 3>& aabb)
	{
		auto minp = aabb.m_min;
		auto maxp = aabb.m_max;
		this->center = (minp + maxp) / 2.0;
		this->m_dir = Eigen::Matrix<_Scalar, 1, 3>(0, 0, 1);
		this->m_rotation = Eigen::Matrix3d::Identity();

		this->m_min = minp - this->center.transpose();
		this->m_max = maxp - this->center.transpose();
	}


	inline OBB(const Eigen::AlignedBox<_Scalar, 3>& aabb, const Eigen::Matrix<_Scalar, 1, 3>& dir, \
		const Eigen::Matrix<_Scalar, 1, 3>& center) : m_dir(dir), m_center(center)
	{
		// 注：方向包围盒的方向取dir和-dir没有区别，所以控制dir和z轴正向夹角不大于pi/2;
		const double epsilon = 1e-6;
		assert((dir.norm() - 1) < epsilon && "input dir is invalid.");

		double dotValue = Eigen::Matrix<_Scalar, 1, 3>(0, 0, 1).dot(dir);
		if (dotValue < 0)
			this->m_dir = -dir;
		Eigen::Matrix<_Scalar, 1, 3> crossValue = Eigen::Matrix<_Scalar, 1, 3>(0, 0, 1).cross(this->m_dir);		// 旋转轴；

		auto minp = aabb.min();
		auto maxp = aabb.max();
		Eigen::Matrix<_Scalar, 1, 3> centerOri = (minp + maxp) / 2.0;
		this->m_min = minp - centerOri.transpose();
		this->m_max = maxp - centerOri.transpose();

		// 若指定方向和z轴平行：
		if (crossValue.norm() < epsilon)
		{
			this->m_rotation = Eigen::Matrix3d::Identity();
			return;
		}

		Eigen::Matrix<_Scalar, 1, 3> zNew = this->m_dir;
		Eigen::Matrix<_Scalar, 1, 3> xNew = -crossValue.normalized();
		Eigen::Matrix<_Scalar, 1, 3> yNew = zNew.cross(xNew);

		for (int i = 0; i < 3; ++i)
		{
			this->m_rotation(i, 0) = static_cast<double>(xNew(i));
			this->m_rotation(i, 1) = static_cast<double>(yNew(i));
			this->m_rotation(i, 2) = static_cast<double>(zNew(i));
		}
	}


	inline OBB(const Eigen::AlignedBox<_Scalar, 3>& aabb, const Eigen::Matrix<_Scalar, 1, 3>& center, \
		const Eigen::Matrix3d& rotation) : m_center(center), m_rotation(rotation)
	{
		const double epsilon = 1e-6;
		auto minp = aabb.min();
		auto maxp = aabb.max();
		Eigen::Matrix<_Scalar, 1, 3> centerOri = (minp + maxp) / 2.0;
		this->m_min = minp - centerOri.transpose();
		this->m_max = maxp - centerOri.transpose();

		this->m_dir = (rotation * Eigen::Matrix<_Scalar, 3, 1>(0, 0, 1)).transpose();
	}


	template<typename Derived>
	inline bool contains(const Eigen::MatrixBase<Derived>& ver) const
	{
		Eigen::Vector3d minp = this->m_min.cast<double>();
		Eigen::Vector3d maxp = this->m_max.cast<double>();

		const Derived& ver0 = ver.derived();
		assert((ver0.rows() == 1 && ver0.cols() == 3) || (ver0.cols() == 1 && ver0.rows() == 3) && "Input vertice must be a 3D vector!");

		// 顶点逆仿射变换到AABB的空间：
		auto dataPtr = ver0.data();
		Eigen::Vector3d p00(static_cast<double>(dataPtr[0]), static_cast<double>(dataPtr[1]), static_cast<double>(dataPtr[2]));
		Eigen::Vector3d p0 = this->m_rotation.inverse() * (p00 - this->m_center.transpose());
		return (minp.array() <= p0.array()).all() && (p0.array() <= maxp.array()).all();
	}
};


template <typename _Scalar>
void genOBBmesh(const OBB<_Scalar>& obb, Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::MatrixXi& tris)
{
	Eigen::Matrix<_Scalar, 3, 1> minp = obb.min();
	Eigen::Matrix<_Scalar, 3, 1> maxp = obb.max();
	Eigen::Matrix<_Scalar, 3, 1> sizeVec = maxp - minp;

	genCubeMesh(vers, tris);
	vers.col(0) *= sizeVec(0);
	vers.col(1) *= sizeVec(1);
	vers.col(2) *= sizeVec(2);

	vers = ((obb.m_rotation * vers.transpose()).eval()).transpose().eval();
	vers.rowwise() += obb.m_center;
}



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
	void test99();
	void test10();

	void test12();
	void test13();
	void test14();
	void test15();
	void test16();
}


// 数字图像处理
namespace TEST_DIP 
{
	void test0();

}


#include "myEigen.tpp"




