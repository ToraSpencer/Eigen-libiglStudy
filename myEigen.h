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


// 和algorithm工程一样，使用单精度；libigl库中封装的三角剖分使用的是双精度；
#define ANSI_DECLARATORS
#define REAL float
#define VOID int
#include "triangulate.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "tmesh.h"								// 拓扑网格类TMESH
#include "detectIntersections.h"			// TMESHS的自相交检测功能；

const double pi = 3.14159265359;

/////////////////////////////////////////////////////////////////////////////////////////////////// 前置声明
template<typename T>
bool matInsertRows(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat1);
template<typename T, int N>
bool matInsertRows(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const Eigen::Matrix<T, 1, N>& rowVec);
template <typename T>
void vecWriteToFile(const char* fileName, const std::vector<T>& vec);
template<typename T>
void vecReadFromFile(std::vector<T>& vec, const char* fileName, const unsigned elemCount);
template<typename T>
bool matWriteToFile(const char* fileName, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat);
template<typename T>
bool matReadFromFile(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName);
template	<typename Scalar, typename Index>
void objReadMeshMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic>& tris, const char* fileName);
template	<typename T>
void objWriteMeshMat(const char* fileName, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris);
template	<typename Scalar>
void objReadVerticesMat(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, const char* fileName);
template	<typename DerivedV>
void objWriteVerticesMat(const char* fileName, const Eigen::PlainObjectBase<DerivedV>& vers);
template	<typename DerivedV, typename DerivedI>
void objWriteEdgesMat(const char* pathName, const Eigen::PlainObjectBase<DerivedI>& edges, \
	const Eigen::PlainObjectBase<DerivedV>& vers);
template <typename DerivedV, typename	 DerivedI>
void objWritePath(const char* pathName, const std::vector<DerivedI>& path, \
	const Eigen::Matrix<DerivedV, Eigen::Dynamic, Eigen::Dynamic>& vers);
template <typename DerivedV>
void objWriteTreePath(const char* pathName, const Eigen::VectorXi& treeVec, \
	const Eigen::Matrix<DerivedV, Eigen::Dynamic, Eigen::Dynamic>& vers);
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
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const std::pair<float, float> sizePair, const float SR, const bool isCovered);
template <typename T>
bool circuit2mesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& circVers);


/////////////////////////////////////////////////////////////////////////////////////////////////// basic math interface

//		二维笛卡尔坐标系转换为极坐标系
template <typename T>
std::pair<double, double> cart2polar(const T x, const T y)
{
	// theta坐标范围为[0, 2*pi]；
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

	// 若radius接近0， 将其对应的theta置0；
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



/////////////////////////////////////////////////////////////////////////////////////////////////// 控制台打印接口

// 并行for循环
template<typename Func>
void PARALLEL_FOR(unsigned int  beg, unsigned int  end, const Func& func, const unsigned int serial_if_less_than = 12)
{
	/*
		PARALLEL_FOR(
			unsigned int   beg,                                     起始元素索引
			unsigned int  end,                                      尾元素索引
			const unsigned int  serial_if_less_than,      如果需要处理的元素不小于该值，则使用并行化；
			const Func & func										操作元素的函数子；
			)
	*/
	unsigned int elemCount = end - beg + 1;

	if (elemCount < serial_if_less_than)
		for (unsigned int i = beg; i < end; ++i)
			func(i);
	else
	{
		// 确定起的线程数；
		const static unsigned n_threads_hint = std::thread::hardware_concurrency();
		const static unsigned n_threads = (n_threads_hint == 0u) ? 8u : n_threads_hint;

		// for循环的范围分解成几段；
		unsigned int slice = (unsigned int)std::round(elemCount / static_cast<double>(n_threads));
		slice = std::max(slice, 1u);

		// 线程函数：
		auto subTraverse = [&func](unsigned int head, unsigned int tail)
		{
			for (unsigned int k = head; k < tail; ++k)
				func(k);
		};

		// 生成线程池，执行并发for循环；
		std::vector<std::thread> pool;              // 线程池；
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

		// 线程同步；
		for (std::thread& t : pool)
		{
			if (t.joinable())
				t.join();
		}
	}
}


// 变参并行for循环——索引以外的参数使用std::tuple传入；
template<typename Func, typename paramTuple>
void PARALLEL_FOR(unsigned int  beg, unsigned int  end, const Func& func, const paramTuple& pt, const unsigned int serial_if_less_than = 12)
{
	/*
		PARALLEL_FOR(
			unsigned int   beg,                                     起始元素索引
			unsigned int  end,                                      尾元素索引
			const unsigned int  serial_if_less_than,      如果需要处理的元素不小于该值，则使用并行化；
			const paramTuple& pt								索引以外的其他参数；
			const Func & func										操作元素的函数子；
			)
	*/
	unsigned int elemCount = end - beg + 1;

	if (elemCount < serial_if_less_than)
		for (unsigned int i = beg; i < end; ++i)
			func(i, pt);
	else
	{
		// 确定起的线程数；
		const static unsigned n_threads_hint = std::thread::hardware_concurrency();
		const static unsigned n_threads = (n_threads_hint == 0u) ? 8u : n_threads_hint;

		// for循环的范围分解成几段；
		unsigned int slice = (unsigned int)std::round(elemCount / static_cast<double>(n_threads));
		slice = std::max(slice, 1u);

		// 线程函数：
		auto subTraverse = [&func, &pt](unsigned int head, unsigned int tail)
		{
			for (unsigned int k = head; k < tail; ++k)
				func(k, pt);
		};

		// 生成线程池，执行并发for循环；
		std::vector<std::thread> pool;              // 线程池；
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

		// 线程同步；
		for (std::thread& t : pool)
		{
			if (t.joinable())
				t.join();
		}
	}
}


// 传入函数子或函数指针遍历stl容器
template<typename T, typename F>
void traverseSTL(T& con, F f)
{
	std::for_each(con.begin(), con.end(), f);
	std::cout << std::endl;
}


// 反向遍历
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


// 传入函数子遍历稀疏矩阵中的非零元素，函数子接受的参数是Eigen::SparseMatrix<T>::InnerIterator&
template<typename spMat, typename F>
void traverseSparseMatrix(spMat& sm, F f)
{
	for (unsigned i = 0; i < sm.outerSize(); ++i)
		for (auto iter = spMat::InnerIterator(sm, i); iter; ++iter)
			f(iter);
}


// lambda——打印std::cout支持的类型变量。
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


// 打印稀疏矩阵中的所有非零元素：
template<typename spMat>				// 如果模板参数为T, 函数参数为Eigen::SparseMatrix<T>，编译会报错，不知道什么缘故；
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



///////////////////////////////////////////////////////////////////////////////////////////////////// 空间变换接口：

// 注！！！计算出的旋转矩阵全部默认为作用于列向量；若v,u为列向量，r,s为行向量：R * v == u; 等价于 r * R.transpose() == s;

// getRotationMat() 重载1——输入旋转轴向量，旋转角度，返回旋转矩阵：
template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& axisArrow, const float theta) 
{
	Eigen::Matrix<T, 3, 1> axis = axisArrow.transpose().normalized();
	return Eigen::AngleAxis<T>(theta, axis).toRotationMatrix();
}


// getRotationMat() 重载2——得到将originArrow旋转到targetArrow的旋转矩阵
template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& originArrow, const Eigen::Matrix<T, 1, 3>& targetArrow)
{
	Eigen::Matrix<T, 3, 3> rotation = Eigen::Matrix<T, 3, 3>::Zero();
	if (0 == originArrow.norm() || 0 == targetArrow.norm())
		return rotation;

	Eigen::Matrix<T, 1, 3> axisArrow = originArrow.cross(targetArrow);		// 旋转轴；

	if (0 == axisArrow.norm())
		return Eigen::Matrix<T, 3, 3>::Identity();

	axisArrow.normalize();
	T x0 = axisArrow(0);
	T y0 = axisArrow(1);
	T z0 = axisArrow(2);
	T cosTheta = originArrow.dot(targetArrow) / (originArrow.norm() * targetArrow.norm());
	T sinTheta = std::sqrt(1 - cosTheta * cosTheta);

	// 等价于Eigen::AngleAxis<T>(theta, axis).toRotationMatrix()，计算绕任意轴向量旋转theta角度；
	rotation << cosTheta + (1 - cosTheta) * x0 * x0, (1 - cosTheta)* x0* y0 - sinTheta * z0, (1 - cosTheta)* x0* z0 + sinTheta * y0, \
		(1 - cosTheta)* y0* x0 + sinTheta * z0, cosTheta + (1 - cosTheta) * y0 * y0, (1 - cosTheta)* y0* z0 - sinTheta * x0, \
		(1 - cosTheta)* z0* x0 - sinTheta * y0, (1 - cosTheta)* z0* y0 + sinTheta * x0, cosTheta + (1 - cosTheta) * z0 * z0;
	return rotation;
}
 


///////////////////////////////////////////////////////////////////////////////////////////////////// 图形生成接口：
 
// 输入起点、终点、空间采样率，插值生成一条直线点云；
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
	for (unsigned i = 1; (length - lenth0) > 0.8 * SR; i++)			// 确保最后一个点距离终点不能太近。
	{
		Eigen::Matrix<T, 1, 3> temp = start + SR * i * dir;
		matInsertRows<T, 3>(vers, temp);
		lenth0 = SR * (i + 1);		// 下一个temp的长度。
	}

	if (SE)
		matInsertRows<T, 3>(vers, end);

	return true;
};


// 生成XOY平面内的圆圈点集：
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


// 生成中心在原点，边长为1，三角片数为12的正方体网格；
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


// 生成轴向包围盒的三角网格；
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


// genCylinder()重载1——生成（类）柱体：
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers, \
		const bool isCovered = true)
{
	/*
	bool genCylinder(
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers,
		Eigen::MatrixXi& tris, 
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers,			轴线；
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers,			横截面回路顶点，必须要在XOY平面内；
		const bool isCovered																						是否封底
		)
	
	
	*/

	// lambda——柱体侧面的三角片生长，会循环调用，调用一次生长一层；
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// 会重复调用，tris容器不需要为空。
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// 待生成表面的圆柱体底圈的第一个顶点的索引。
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
	unsigned circVersCount = btmVers.rows();					// 横截面一圈的顶点数；
	unsigned circCount = axisVers.rows();							// 圈数；
	unsigned versCount = circVersCount * circCount;
	std::vector<MatrixXT> circuitsVec(circCount);									// 每个横截面的一圈顶点；
	std::vector<RowVector3T> sectionNorms(circCount);					// 每个横截面的法向；

	// 1. 计算柱体circCount个横截面的法向、每个横截面的顶点；
	for (unsigned i = 0; i < circCount - 1; ++i)
	{
		sectionNorms[i] = axisVers.row(i + 1) - axisVers.row(i);
		sectionNorms[i].normalize();
		Matrix3T rotation = getRotationMat(RowVector3T{0, 0, 1}, sectionNorms[i]);
		circuitsVec[i] = btmVers * rotation.transpose();
		circuitsVec[i].rowwise() += axisVers.row(i);
	}

	// 2. 计算最后一圈顶点：
	RowVector3T deltaNormAve{RowVector3T::Zero()};
	for (unsigned i = 0; i < circCount - 2; ++i)
		deltaNormAve += (sectionNorms[i + 1] - sectionNorms[i]);
	deltaNormAve.array() /= (circCount - 2);
	sectionNorms[circCount - 1] = sectionNorms[circCount - 2] + deltaNormAve;
	Matrix3T rotation = getRotationMat(RowVector3T{ 0, 0, 1 }, sectionNorms[circCount - 1]);
	circuitsVec[circCount - 1] = btmVers * rotation.transpose();
	circuitsVec[circCount - 1].rowwise() += axisVers.row(circCount - 1);

	// 3. 生成柱体顶点：
	vers.resize(versCount, 3);
	for (unsigned i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.生成侧面三角片：
	for (unsigned i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);

	
	// 5. 加盖：
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


// genCylinder()重载2——生成圆柱体网格
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const float radius, const double deltaTheta = 2 * pi / 30, const bool isCovered = true)
{
	// 生成XOY平面上采样角度步长为的圆圈顶点：
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


// genCylinder()重载3——生成方柱体网格
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const std::pair<float, float> sizePair, const float SR = 0.5, \
	const bool isCovered = true)
{
	// 生成XOY平面上采样角度步长为的圆圈顶点：
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// 生成XOY平面内的方框顶点：
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


//			circuitToMesh()重载1：triangle库三角剖分——封闭边界线点集得到面网格，可以是平面也可以是曲面，三角片尺寸不可控，不会在网格内部插点。
template <typename T>
bool circuit2mesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& circVers)
{
	using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;

	unsigned circCount = circVers.rows();
	unsigned versCount = circVers.rows();

	// 0. 边缘环路顶点坐标数据拷入到输出网格中。
	vers = circVers;

	// 输入triangulate()的点集是投影到XOY平面的二维点，原点集应该旋转到合适的角度再投影。

	// 1. 取第一个点、1/3处的点、2/3处的点的所在平面的法向量作为原点集的法向量
	RowVector3T vers1 = circVers.row(0);
	RowVector3T vers2 = circVers.row(versCount / 3);
	RowVector3T vers3 = circVers.row(2 * versCount / 3);
	RowVector3T norm = (vers1 - vers2).cross(vers3 - vers2);
	norm.normalize();

	//  2. 旋转点集使得norm平行于z轴
	Matrix3T rotation = getRotationMat(norm, RowVector3T{0, 0, 1});
	vers = (vers * rotation.transpose()).eval();

	// for debug;
	objWriteVerticesMat("E:/versafterrotation.obj", vers);

	// 3. 旋转后的点数据写入到triangulate()接口的输入结构体中。
	Eigen::MatrixXf vers2D;
	Eigen::MatrixXi edges2D;
	vers2D = vers.transpose().topRows(2).cast<float>();
	edges2D.resize(2, versCount);
	for (unsigned i = 0; i < versCount; i++)
	{
		edges2D(0, i) = i + 1;
		edges2D(1, i) = (i + 1) % versCount + 1;
	}

	triangulateio triIn, triOut;
	triIn.numberofpoints = versCount;
	triIn.pointlist = (float*)vers2D.data();
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

	// 4. 执行二维三角剖分，得到输出网格三角片数据，不取顶点坐标数据，顶点坐标数据使用旋转操作前的。
	char triStr[256] = "pY";
	triangulate(triStr, &triIn, &triOut, NULL);
	tris.resize(3, triOut.numberoftriangles);
	MatrixXT norms(versCount, 3);
	std::memcpy(tris.data(), triOut.trianglelist, sizeof(int) * 3 * triOut.numberoftriangles);
	tris.transposeInPlace();
	tris.array() -= 1;
 
	return true;
}
 


///////////////////////////////////////////////////////////////////////////////////////////////////// 不同数据类型的变换

// 根据索引向量从源矩阵中提取元素生成输出矩阵。
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


// 根据flag向量从源矩阵中提取元素生成输出矩阵。
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


// eigen的向量和std::vector<T>相互转换
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


// 生成边编码——边两个端点的索引对映射成一个64位整型数
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


// 生成无向边编码——边两个端点的索引对整理成前小后大的顺序，映射成一个64位整型数
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


// 解码边编码；
std::pair<int, int> decodeEdge(const std::int64_t code);


// 生成三角片编码——三个顶点索引排序后映射成64位无符号整型数
template <typename Index>
std::uint64_t encodeTriangle(const Index vaIdx, const Index vbIdx, const Index vcIdx)
{
	unsigned long triIdxLimit = 0x1FFFFF;								// 最多为21位全1, 0x1FFFFF ==  2097181, 两百多万三角片；
	if (vaIdx > triIdxLimit || vbIdx > triIdxLimit || vcIdx > triIdxLimit)
		return 0;			// 索引超出范围
	if (vaIdx == vbIdx || vaIdx == vcIdx || vbIdx == vcIdx || vaIdx < 0 || vbIdx < 0 || vcIdx < 0)
		return 0;			// 非法三角片

	// 区分正反面——即(a, b, c)和(a, c, b)会映射成不同的编码；但是(a,b,c)(b,c,a)(c,a,b)映射成相同的编码；
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


// 解码三角片编码：
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




///////////////////////////////////////////////////////////////////////////////////////////////////////// 矩阵的增删查改

// 向量插入数据
template<typename T>
bool vecInsertNum(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const T num)
{
	vec.conservativeResize(vec.rows() + 1, 1);
	vec(vec.rows() - 1) = num;
	return true;
}


// 向量插入向量
template<typename T>
bool vecInsertVec(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec2)
{
	unsigned currentRows = vec1.rows();
	unsigned addRows = vec2.rows();
	unsigned finalRows = currentRows + addRows;
	vec1.conservativeResize(finalRows);

	// 拷贝数据：
	T* dataPtr = vec1.data();
	dataPtr += currentRows;
	std::memcpy(dataPtr, vec2.data(), sizeof(T) * addRows);

	return true;
}


//	矩阵末端插入矩阵
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


//	矩阵末端插入行向量
template<typename T, int N>					// N为列数；
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


// 返回一个flag列向量retVec，若mat的第i行和行向量vec相等，则retVec(i)==1，否则等于0；若程序出错则retVec所有元素为-1
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

	// 逐列比较：
	Eigen::MatrixXi tempMat(rows, cols);
	for (int i = 0; i < cols; ++i)
		tempMat.col(i) = (mat.col(i).array() == rowVec(i)).select(Eigen::VectorXi::Ones(rows), Eigen::VectorXi::Zero(rows));

	retVec = tempMat.col(0);

	if (cols > 1)
		for (int i = 1; i < cols; ++i)
			retVec = retVec.array() * tempMat.col(i).array();			// 逐列相乘：

	return retVec;
}


// 网格串联——合并两个孤立的网格到一个网格里
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





//////////////////////////////////////////////////////////////////////////////////////////// IO接口

unsigned readNextData(char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize);


// 将std::vector中的数据写入到文件中，二进制形式，vector中的元素不能是pointer-like类型
template <typename T>
void vecWriteToFile(const char* fileName, const std::vector<T>& vec)
{
	std::ofstream file(fileName, std::ios_base::out | std::ios_base::binary);
	file.write((char*)(&vec[0]), sizeof(T) * vec.size());
	file.close();
}


// to be optimized——元素数应该不需要传入函数，应该可以由文件流大小推断出来。
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

	// 读矩阵尺寸信息
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

	// 读矩阵元素
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

			if (*iter == ' ')						// 负号后面可能有空格，需要去除
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
		char szBuf[256] = { 0 };
		sprintf_s(szBuf, 256, "v %f %f %f", static_cast<float>(vers(j, 0)), static_cast<float>(vers(j, 1)), static_cast<float>(vers(j, 2)));
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

		// 顶点信息		
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


// 边数据写入到OBJ文件中：
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



// objWritePath() 路径数据写入到OBJ文件中：
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
 

// objWirteTreePath()——输入树向量或路径向量，写入到OBJ文件中：
template <typename DerivedV>
void objWriteTreePath(const char* pathName, const Eigen::VectorXi& treeVec, const Eigen::Matrix<DerivedV, Eigen::Dynamic, Eigen::Dynamic>& vers)
{
	// 路径被视为树的特例；
	
	// 树中索引为i的顶点的前驱顶点索引为treeVec(i), 若其没有前驱顶点（父节点），则treeVec(i) == -1;
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

 


///////////////////////////////////////////////////////////////////////////////////////////////////////// 齐次坐标系相关接口
void objReadVerticesHomoMat(Eigen::MatrixXf& vers, const char* fileName);

void objWriteVerticesHomoMat(const char* fileName, const Eigen::MatrixXf& vers);

void vers2homoVers(Eigen::MatrixXf& homoVers, const Eigen::MatrixXf& vers);

Eigen::MatrixXf vers2homoVers(const Eigen::MatrixXf& vers);

void homoVers2vers(Eigen::MatrixXf& vers, const Eigen::MatrixXf& homoVers);

Eigen::MatrixXf homoVers2vers(const Eigen::MatrixXf& homoVers);





//////////////////////////////////////////////////////////////////////////////////////////// 科学计算相关

// 稀疏矩阵转置；Eigen::SparseMatrix自带的transpose()方法太垃圾了
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


// 解恰定的稠密线性方程组Ax == b;
template<typename T, int N>
bool solveLinearEquation(Eigen::Matrix<T, N, 1>& x, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, const Eigen::Matrix<T, N, 1>& b)
{
	// 解线性方程组Ax == b;
	Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	x = svd.solve(b);

	return true;
}


// 解一系列恰定的稠密线性方程组AX == B;
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


// 霍纳方法（秦九昭算法）求多项式的值
template<typename T, int N>
double hornersPoly(const Eigen::Matrix<T, N, 1>& coeffs, const double x)
{
	// coeffs是{a0, a1, a2, ..., an}组成的(n+1)列向量，多项式为p == a0 + a1*x + a2*x^2 + ... + an* x^n; 
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


// 求多项式的一阶微分
template<typename T, int N>
float polyDiff(const Eigen::Matrix<T, N, 1>& coeffs, const float x)
{
	// 多项式p == a0 + a1*x + a2*x^2 + ... + an* x^n 一阶微分为：p' == a1 + 2*a2*x + 3*a3*x^2 ... n*an*x^(n-1)
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


// 计算稠密矩阵的克罗内克积（Kronecker product）
template <typename T, typename Derived1, typename Derived2>
void kron(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& result, \
	const Eigen::MatrixBase<Derived1>& mm1, \
	const Eigen::MatrixBase<Derived2>& mm2)
{
	// A(m×n) tensor B(p×q) == C(m*p × n*q);
	/*
		其中C.block(i * m, j* n, p, q) == Aij * B;
		如：
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

	result.resize(0, 0);						// 重新分配内存；
	result.resize(m * p, n * q);
	Eigen::VectorXi rowIdxVec = Eigen::VectorXi::LinSpaced(0, m - 1, m - 1);
	auto calcKronCol = [&](const unsigned colIdx, const std::tuple<unsigned>& pt)
	{
		unsigned rowIdx0 = std::get<0>(pt);
		result.block(rowIdx0 * p, colIdx * q, p, q) = static_cast<T>(m1(rowIdx0, colIdx)) * m2.array().cast<T>();		// 转换成输出矩阵中的数据类型；
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


// 多项式插值
void polyInterpolation();


// 高斯插值
void gaussInterpolation();


// 最小二乘多项式拟合曲线：
void leastSquarePolyFitting();


// 岭回归多项式拟合曲线
template <typename T>
void ridgeRegressionPolyFitting(Eigen::VectorXd& theta, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, unsigned m = 4)
{
	/*
		void ridgeRegressionPolyFitting(
					Eigen::VectorXd & theta,							拟合的多项式函数
					const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers			离散样本点
					const unsigned m						拟合的多项式项数；
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



// 最小二乘法拟合（逼近）标准椭圆（长短轴和xy坐标轴对齐）
template<typename T>
Eigen::VectorXd fittingStandardEllipse(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& sampleVers)
{
	/*
		Eigen::VectorXd fittingStandardEllipse(								// 返回列向量(a,c,d,e,f)，为标准椭圆方程的系数；
					const Eigen::MatrixXf& sampleVers					// 输入的样本点，必须是在XOY平面上的点；
		)
	*/
	const double epsilon = 1e-8;							// 浮点数绝对值小于此值时认为为0；
	using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

	// 标准椭圆方程：a*x^2 + c*y^2 + d*x + e*y + f  = 0，其中a*c > 0
	Eigen::VectorXd x(Eigen::VectorXd::Zero(5));

	unsigned m = sampleVers.rows();				// sample count;
	Eigen::VectorXd x0(Eigen::VectorXd::Zero(m));
	Eigen::VectorXd y0(Eigen::VectorXd::Zero(m));
	for (unsigned i = 0; i < m; ++i)
	{
		x0(i) = static_cast<double>(sampleVers(i, 0));
		y0(i) = static_cast<double>(sampleVers(i, 1));
	}

	// alpha = [x^2, y^2, x, y, 1]; 样本信息矩阵：A = [alpha1; alpha2; .... alpham]; 椭圆方程写为：A*x = 0;
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

	// 求S的特征值，特征向量：
	Eigen::EigenSolver<Eigen::MatrixXd> es(S);
	Eigen::MatrixXd D = es.pseudoEigenvalueMatrix();			// 对角线元素是特征值
	Eigen::MatrixXd V = es.pseudoEigenvectors();					// 每一个列向量都是特征向量。

	// 寻找特征值不为0，且满足约束条件的特征向量：
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





//////////////////////////////////////////////////////////////////////////////////////////////// 三角网格处理：

// 得到三角网格的有向边数据
template <typename DerivedV, typename DerivedI>
void getEdges(Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
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
}
 

// 得到坐标向量表示的有向边：
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


// 计算网格所有三角片的重心
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


// 计算网格所有三角片的归一化法向量
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

	// 法向量归一化，若存在退化三角片，则写为inf
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


// 计算网格所有三角片所在平面
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

		// 若存在退化三角片，则平面系数都写为NAN:
		if (std::isinf(norm(0)))
			planeCoeff(i, 3) = INFINITY;

		Eigen::RowVector3d va = vers.row(tris(i, 0)).array().cast<double>();

		// va点在平面上 → va点到平面的距离为0 → norm.dot(va) +d == 0; → d = -norm.dot(va);
		planeCoeff(i, 3) = -norm.dot(va);
	}

	return true;
}


// 计算三角网格每个三角片的面积：
template<typename Scalar, typename DerivedI>
bool trisArea(Eigen::VectorXd& trisAreaVec, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
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
	s = (lens1 + lens2 + lens3) / 2.0;				// 三角片周长的一半；

	// 使用Heron公式计算三角片面积:
	trisAreaVec = s.array() * (s - lens1).array() * (s - lens2).array() * (s - lens3).array();
	trisAreaVec = trisAreaVec.cwiseSqrt();
	
	return true;
}


// 计算三角网格的体积：
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

	// 一个三角片对应的四面体的符号体积（signed volume） V(t0) == (-x3y2z1 + x2y3z1 + x3y1z2 - x1y3z2 - x2y1z3 + x1y2z3) / 6.0;
	auto volumeVec = \
		(-versC.col(0).array() * versB.col(1).array() * versA.col(2).array() + versB.col(0).array() * versC.col(1).array() * versA.col(2).array()\
			+ versC.col(0).array() * versA.col(1).array() * versB.col(2).array() - versA.col(0).array() * versC.col(1).array() * versB.col(2).array()\
			- versB.col(0).array() * versA.col(1).array() * versC.col(2).array() + versA.col(0).array() * versB.col(1).array() * versC.col(2).array()) / 6.0;

	volume = volumeVec.sum();

	return volume;
}


// 求三角网格的不同权重的邻接矩阵 
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

	// 1. 求顶点邻接关系：

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
	adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());										// 权重为该有向边重复的次数；
	adjSM_eIdx.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());		// 权重为该有向边的索引；

	return true;
}



// 边缘有向边，重载1：
template<typename DerivedI>
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
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
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, std::vector<int>& bdryTriIdxes, const Eigen::PlainObjectBase<DerivedI>& tris)
{
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
		edgesMap.insert({ codeA, i});
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


// 求三角网格中的非流形半边：
template<typename DerivedI>
bool nonManifoldEdges(const Eigen::PlainObjectBase<DerivedI>& tris, Eigen::MatrixXi& nmnEdges)
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
 

// buildAdjacency()——计算网格的三角片邻接信息：
using ttTuple = std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>;
bool buildAdjacency(const Eigen::MatrixXi& tris, Eigen::MatrixXi& ttAdj_nmEdge, \
	std::vector<ttTuple>& ttAdj_nmnEdge, std::vector<ttTuple>& ttAdj_nmnOppEdge);


// triangleGrow()——从指定三角片开始区域生长；输入网格不可以有非流形边
template <typename DerivedV, typename DerivedI>
bool triangleGrow(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
		const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx)
{
	/*
		
			bool  triangleGrow(
						versOut,							// 输出网格的点云
						trisOut,							// 输出网格的三角片
						vers,								// 输入网格的点云
						tris,								// 输入网格的三角片
						triIdx							// 三角片生长的起点三角片索引。
						)

			输入网格若存在非流形边，则会return false
	*/
	int versCount = vers.rows();
	int trisCount = tris.rows();
	int edgesCount = 3 * trisCount;

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_eIdx的triplet数据；	

	// 1. 求基本的边信息、邻接矩阵：
	Eigen::SparseMatrix<int> adjSM;								// 邻接矩阵；
	Eigen::SparseMatrix<int> adjSM_eCount;					// 邻接矩阵，索引为该有向边重复的次数；
	Eigen::SparseMatrix<int> adjSM_eIdx;				// 邻接矩阵，元素若对应流形边则为该边索引，若对应非流形边则为该边所有索引的和；
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;		// adjSM_eIdx的转置，表示对边的信息；
	Eigen::SparseMatrix<int> adjSM_ueCount;				// 权重为无向边ij关联的三角片数量；
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;				// 非边缘流形无向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB;					// 非边缘流形有向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;			// 非边缘流形有向边对边的邻接矩阵
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

		// 1.1 生成三种有向边邻接矩阵：adjSM, adjSM_eCount, smElems_weighted
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
		adjSM = adjSM_eCount;																		// 有向边邻接矩阵；
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		若adjSM_eIdx(i, j)对应的是流形有向边，该权重值为该边的索引；若是非流形有向边，则该边对应两个边索引，权重为两个索引之和；
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

		//	1.2 生成无向边邻接矩阵：
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// 无向边邻接矩阵；
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.3 生成非边缘流形无向边邻接矩阵：
		smElems.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				if (2 == iter.value())
					smElems.push_back(Eigen::Triplet<int>(iter.row(), iter.col(), 1));
			});
		adjSM_ueMN_NB.resize(versCount, versCount);
		adjSM_ueMN_NB.setFromTriplets(smElems.begin(), smElems.end());
		smElems.clear();

		//	1.4 生成非边缘流形有向边、及其对边的邻接矩阵
		adjSM_MN_NB = adjSM + adjSM_ueMN_NB;					// adjSM & adjSM_ueMN_NB
		adjSM_MN_NB.prune([](const Eigen::Index& row, const Eigen::Index& col, const int& value)->bool
			{
				if (2 == value)			// 删除非流形边
					return true;
				else
					return false;
			});
		adjSM_MN_NB /= 2;
		adjSM_MN_NB_opp = adjSM_MN_NB.transpose();
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
	std::vector<int> etInfo;												// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
	Eigen::VectorXi etAdj_mnEdge;							// 边所在三角片的索引，非流形边或边缘边写-1；
	Eigen::MatrixXi ttAdj_nmEdge;					// 三角片的非边缘流形有向边邻接的三角片索引；	trisCount * 3(i, 0)为索引为i的三角片的非边缘流形有向边(vb, vc)邻接的三角片的索引；
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
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//	3.3 求三角片邻接矩阵；ttAdj_mnEdge(i, :)是索引为i的三角片三条边邻接的三个三角片，若其中有非流形边或边缘边则写为-1；
		ttAdj_nmEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. 循环——队列中第一个三角片t出队 - t写入输出集合 - 求出t所有相邻三角片的索引然后入队
	std::unordered_set<int> finalTriIdx;						// 联通三角片集合；
	std::unordered_set<int> workingSet;						// 用于寻找相邻三角片的队列；
	workingSet.insert(static_cast<int>(triIdx));							// 队列插入第一个三角片；
	while (workingSet.size() > 0)
	{
		// 4.1 首个元素出队，插入到联通三角片集合中；
		int cTriIdx = *workingSet.begin();									// 当前处理的三角片索引
		std::vector<int> adjTriIdx;												// 当前三角片相邻三角片的索引
		workingSet.erase(workingSet.begin());
		finalTriIdx.insert(cTriIdx);

		// 4.2 确定当前三角片的相邻三角片
		for (int i = 0; i < 3; ++i)
			if (ttAdj_nmEdge(cTriIdx, i) >= 0)
				adjTriIdx.push_back(ttAdj_nmEdge(cTriIdx, i));

		// 4.3 若联通三角片集合中没有当前的三角片，插入该三角片，队列中也插入该三角片；
		for (const auto& index : adjTriIdx)
		{
			auto retPair = finalTriIdx.insert(index);
			if (retPair.second)
				workingSet.insert(index);
		}
	}


	// 5. 由选定的三角片取顶点：
	std::set<int> finalVerIdx;
	Eigen::VectorXi finalVerIdxVec;													// 输出顶点的索引，是原网格中的索引；
	Eigen::VectorXi oldNewIdxInfo = -Eigen::VectorXi::Ones(vers.rows());		// 新老索引表——下标为老索引，元素值为新索引，不在新网格中的顶点写为-1
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


	// 6. 原网格所有三角片中的顶点索引由老的换为新的。
	Eigen::MatrixXi tempMat = tris;
	int* intPtr = tempMat.data();
	for (int i = 0; i < tempMat.size(); ++i)
	{
		int oldIdx = *intPtr;
		*intPtr = oldNewIdxInfo[oldIdx];
		intPtr++;
	}

	// 7. 去除非法三角片，即含有-1元素的三角片；
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
 

// triangleGrowSplitMesh()——区域生长将输入网格分为多个单连通网格；输入网格不可以有非流形边
template <typename DerivedV>
bool triangleGrowSplitMesh(std::vector<DerivedV>& meshesVersOut, std::vector<Eigen::MatrixXi>& meshesTrisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	int versCount = vers.rows();
	int trisCount = tris.rows();
	int edgesCount = 3 * trisCount;

	// 输入参数中的triIdx为包含最高点的三角片索引；
	std::unordered_set<int> finalTriIdx;						// 联通三角片集合；
	std::unordered_set<int> workingSet;						// 用于寻找相邻三角片的队列；

	Eigen::MatrixXi ttAdj_nmEdge;		// 三角片的非边缘流形有向边邻接的三角片索引；	trisCount * 3(i, 0)为索引为i的三角片的非边缘流形有向边(vb, vc)邻接的三角片的索引；

	Eigen::SparseMatrix<int> adjSM;								// 邻接矩阵；
	Eigen::SparseMatrix<int> adjSM_eCount;					// 邻接矩阵，索引为该有向边重复的次数；
	Eigen::SparseMatrix<int> adjSM_eIdx;				// 邻接矩阵，元素若对应流形边则为该边索引，若对应非流形边则为该边所有索引的和；
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;		// adjSM_eIdx的转置，表示对边的信息；
	Eigen::SparseMatrix<int> adjSM_ueCount;				// 权重为无向边ij关联的三角片数量；
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;				// 非边缘流形无向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB;					// 非边缘流形有向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;			// 非边缘流形有向边对边的邻接矩阵

	Eigen::VectorXi etAdj_mnEdge;							// 边所在三角片的索引，非流形边或边缘边写-1；

	std::vector<int> etInfo;												// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
	std::vector<int> edgesIdx_MN_NB;							// 非边缘流形有向边的索引；
	std::vector<int> edgesIdx_MN_NB_opp;					// 非边缘流形有向边的对边的索引；
	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_eIdx的triplet数据；	

	std::vector<bool> trisFlag(trisCount, true);		// true表示该三角片可收集，false表示已收集；


	// 1. 求基本的边信息、邻接矩阵：
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

		// 1.1 生成三种有向边邻接矩阵：adjSM, adjSM_eCount, smElems_weighted
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
		adjSM = adjSM_eCount;																		// 有向边邻接矩阵；
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		若adjSM_eIdx(i, j)对应的是流形有向边，该权重值为该边的索引；若是非流形有向边，则该边对应两个边索引，权重为两个索引之和；
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

		//	1.2 生成无向边邻接矩阵：
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// 无向边邻接矩阵；
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.3 生成非边缘流形无向边邻接矩阵：
		smElems.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				if (2 == iter.value())
					smElems.push_back(Eigen::Triplet<int>(iter.row(), iter.col(), 1));
			});
		adjSM_ueMN_NB.resize(versCount, versCount);
		adjSM_ueMN_NB.setFromTriplets(smElems.begin(), smElems.end());
		smElems.clear();

		//	1.4 生成非边缘流形有向边、及其对边的邻接矩阵
		adjSM_MN_NB = adjSM + adjSM_ueMN_NB;					// adjSM & adjSM_ueMN_NB
		adjSM_MN_NB.prune([](const Eigen::Index& row, const Eigen::Index& col, const int& value)->bool
			{
				if (2 == value)			// 删除非流形边
					return true;
				else
					return false;
			});
		adjSM_MN_NB /= 2;
		adjSM_MN_NB_opp = adjSM_MN_NB.transpose();
	}


	// 2. 确定非边缘流形有向边、及其对边的索引；
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
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//	3.3 求三角片邻接矩阵；ttAdj_mnEdge(i, :)是索引为i的三角片三条边邻接的三个三角片，若其中有非流形边或边缘边则写为-1；
		ttAdj_nmEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. 循环——队列中第一个三角片t出队 - t写入输出集合 - 求出t所有相邻三角片的索引然后入队
	std::vector<std::unordered_set<int>> connectedTriIdxes;
	int collectedTrisCount = 0;
	while (collectedTrisCount < trisCount) 
	{
		finalTriIdx.clear();

		// 取当前首个未收集的三角片作为种子三角片：
		int seedTriIdx = 0;
		for (int i = 0; i < trisCount; ++i)
			if (trisFlag[i])
				seedTriIdx = i;

		workingSet.insert(static_cast<int>(seedTriIdx));					// 队列插入种子三角片；
		while (workingSet.size() > 0)
		{
			// 4.1 首个元素出队，插入到联通三角片集合中；
			int cTriIdx = *workingSet.begin();									// 当前处理的三角片索引
			trisFlag[cTriIdx] = false;

			std::vector<int> adjTriIdx;												// 当前三角片相邻三角片的索引
			workingSet.erase(workingSet.begin());
			finalTriIdx.insert(cTriIdx);

			// 4.2 确定当前三角片的相邻三角片
			for (int i = 0; i < 3; ++i)
				if (ttAdj_nmEdge(cTriIdx, i) >= 0)
					adjTriIdx.push_back(ttAdj_nmEdge(cTriIdx, i));

			// 4.3 若联通三角片集合中没有当前的三角片，插入该三角片，队列中也插入该三角片；
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


	// 5. 输出所有单连通网格
	int meshesCount = connectedTriIdxes.size();
	meshesVersOut.resize(meshesCount);
	meshesTrisOut.resize(meshesCount);
	for (int i = 0; i < meshesCount; ++i) 
	{
		DerivedV& cMeshVers = meshesVersOut[i];
		Eigen::MatrixXi& cMeshTris = meshesTrisOut[i];

		// 5.1 由connectedTriIdxes中连通的三角片索引取顶点：
		std::set<int> cMeshVerIdx;																	// 当前的单连通网格中包含的顶点的索引；
		for (const auto& index : connectedTriIdxes[i])
		{
			cMeshVerIdx.insert(tris(static_cast<int>(index), 0));
			cMeshVerIdx.insert(tris(static_cast<int>(index), 1));
			cMeshVerIdx.insert(tris(static_cast<int>(index), 2));
		}

		Eigen::VectorXi cMeshVerIdxVec(cMeshVerIdx.size());									// 输出顶点的索引，是原网格中的索引；
		Eigen::VectorXi oldNewIdxInfo = -Eigen::VectorXi::Ones(vers.rows());			// 新老索引表——下标为老索引，元素值为新索引，不在新网格中的顶点写为-1
		auto iter = cMeshVerIdx.begin();
		for (int i = 0; i < cMeshVerIdx.size(); ++i)
			cMeshVerIdxVec(i) = static_cast<int>(*iter++);

		int setOrder = 0;
		for (const auto& index : cMeshVerIdx)
			oldNewIdxInfo[static_cast<int>(index)] = setOrder++;
		subFromIdxVec(cMeshVers, vers, cMeshVerIdxVec);

		// 5.2 确定当前的单连通网格的三角片：
		std::vector<int> cTriIdxes;
		cTriIdxes.insert(cTriIdxes.end(), connectedTriIdxes[i].begin(), connectedTriIdxes[i].end());
		subFromIdxVec(cMeshTris, tris, cTriIdxes);
		int* verIdxPtr = cMeshTris.data();
		for (int k = 0; k < 3 * cMeshTris.rows(); ++k)				// 三角片中的老的顶点索引改为新的： 
		{
			int oldIdx = *verIdxPtr;
			*verIdxPtr = oldNewIdxInfo[oldIdx];
			verIdxPtr++;
		}
	}

	return true;
}


// triangleGrowOuterSurf()——从指定三角片（必须是外部的三角片）开始区域生长，提取外层表面单连通流形网格
template <typename DerivedV, typename DerivedI>
bool triangleGrowOuterSurf(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx)
{
	// ！！！目前假定输入网格的所有三角片朝向正确

	int versCount = vers.rows();
	int trisCount = tris.rows();
	int edgesCount = 3 * trisCount;
	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_eIdx的triplet数据；	

	// 1. 求三角片邻接关系：
	Eigen::MatrixXi ttAdj_nmEdge;
	std::vector<ttTuple> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	bool retFlag = buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);

	// 2. 循环——队列中第一个三角片t出队 - t写入输出集合 - 求出t所有相邻三角片的索引然后入队
	std::unordered_set<int> finalTriIdx;								// 最终保留的三角片索引集合；
	std::unordered_set<int> workingSet;								// 用于寻找相邻三角片的队列；
	std::unordered_set<int> deleTriIdx;							// 删除的三角片索引；
	workingSet.insert(static_cast<int>(triIdx));							// 队列插入第一个三角片；
	while (workingSet.size() > 0)
	{
		// 4.1 首个元素出队，插入到联通三角片集合中；
		int cTriIdx = *workingSet.begin();									// 当前处理的三角片索引
		std::vector<int> adjTriIdx;												// 当前三角片相邻三角片的索引
		workingSet.erase(workingSet.begin());
		finalTriIdx.insert(cTriIdx);

		// 4.2 确定当前三角片的流形边关联的相邻三角片
		for (int i = 0; i < 3; ++i)
		{
			if (ttAdj_nmEdge(cTriIdx, i) >= 0)
				adjTriIdx.push_back(ttAdj_nmEdge(cTriIdx, i));
			else
			{
				// 若当前边不是流形边——选取同侧的，

			}
		}

		// 4.3 若联通三角片集合中没有当前的三角片，插入该三角片，队列中也插入该三角片；
		for (const auto& index : adjTriIdx)
		{
			auto retPair = finalTriIdx.insert(index);
			if (retPair.second)
				workingSet.insert(index);
		}
	}


	// 3. 由选定的三角片取顶点：
	std::set<int> finalVerIdx;
	Eigen::VectorXi finalVerIdxVec;													// 输出顶点的索引，是原网格中的索引；
	Eigen::VectorXi oldNewIdxInfo = -Eigen::VectorXi::Ones(vers.rows());		// 新老索引表——下标为老索引，元素值为新索引，不在新网格中的顶点写为-1
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


	// 6. 原网格所有三角片中的顶点索引由老的换为新的。
	Eigen::MatrixXi tempMat = tris;
	int* intPtr = tempMat.data();
	for (int i = 0; i < tempMat.size(); ++i)
	{
		int oldIdx = *intPtr;
		*intPtr = oldNewIdxInfo[oldIdx];
		intPtr++;
	}

	// 7. 去除非法三角片，即含有-1元素的三角片；
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


// 无向边的区域生长，提取环路边；输入的边矩阵可以是有向边也可以是无向边；
template <typename DerivedI>
bool uEdgeGrow(std::vector<Eigen::MatrixXi>& circs, std::vector<Eigen::MatrixXi >& segs, \
	const Eigen::PlainObjectBase<DerivedI>& uEdges)
{
	std::list<std::list<std::pair<int, int>>> circList, segList;
	Eigen::SparseMatrix<int> adjSM;										// 邻接矩阵；
	std::vector<Eigen::Triplet<int>> smElems;

	// 1. 建立无向边邻接矩阵：
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

	// 2. 统计不重复的无向边数量；
	std::unordered_set<std::int64_t> tmpSet;
	for (unsigned i = 0; i < uEdges.rows(); ++i)
		tmpSet.insert(encodeUedge(uEdges(i, 0), uEdges(i, 1)));
	int remainCount = tmpSet.size();			// (3,4) (4,3)表示同一条无向边；
	

	// 3. 无向边生长的循环：
	while (remainCount > 0)
	{
		std::list<std::pair<int, int>> currentSeg;
	
		// w1. 提取邻接表中第一个无向边：
		bool breakFlag = false;
		for (unsigned i = 0; !breakFlag & i < adjSM.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, i); !breakFlag & iter; ++iter)
			{
				if (iter.valueRef() > 0)
				{
					currentSeg.push_back({ iter.row(), iter.col() });
					iter.valueRef() = -1;												// 表示删除该元素；
					adjSM.coeffRef(iter.col(), iter.row()) = -1;
					remainCount--;
					breakFlag = true;
				}
			}
		}

		// w2. 搜索第一条无向边关联的边，直到搜索不到，或者形成环路为止。
		int head = currentSeg.front().first;
		int tail = currentSeg.front().second;
		while (1)
		{
			int otherEnd = -1;
			breakFlag = false;
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, head); !breakFlag & iter; ++iter)	// 对第head列的遍历：
			{
				if (iter.valueRef() > 0)
				{
					otherEnd = iter.row();
					currentSeg.push_back({ head, otherEnd });
					iter.valueRef() = -1;												// 表示删除该元素；
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

	
	// 4. 输出：
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


// 检测网格中是否有非法三角片（三个顶点索引中有两个相同的）
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


// 确定三角网格中的所有单连通区域
template <typename Index>
int simplyConnectedRegion(const Eigen::SparseMatrix<Index>& adjSM,
	Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount)
{
	/*
		int simplyConnectedRegion(									返回单连通区域的数量，出错时返回-1
					const Eigen::SparseMatrix<Index>& adjSM,			三角网格的邻接矩阵
					Eigen::VectorXi& connectedLabels,						同一连通区域内的顶点标签，依次为0, 1, 2,...；
																						索引为i的顶点的标签为connectedLabels(i);
					Eigen::VectorXi& connectedCount							每个标签对应的单连通区域内包含的顶点数
																						标签为i的单连通区域包含的顶点数为connectedCount(i)
					)
	*/
	if (adjSM.rows() == 0)
		return -1;

	if (adjSM.rows() != adjSM.cols())
		return -1;

	const unsigned versCount = adjSM.rows();

	// 1. 初始化：
	connectedLabels.setConstant(versCount, versCount);		// 最终结果中最大可能的标签为versCount-1;
	connectedCount.setZero(versCount);

	// 2. 遍历所有顶点，搜索与其联通的顶点：
	int currentLabel = 0;
	for (unsigned i = 0; i < versCount; ++i) 
	{
		// 2.1 若当前顶点已被访问过(label < versCount)，则continue:
		if (connectedLabels(i) < versCount)
			continue;

		// 2.2 若当前顶点未被访问过，执行BFS收集其所有联通的顶点：
		std::queue<Index> workingQueue;
		workingQueue.push(i);
		while (!workingQueue.empty())
		{
			const Index curVerIdx = workingQueue.front();
			workingQueue.pop();
			
			if (connectedLabels(curVerIdx) < versCount)
				continue;

			// 2.2.1 队首顶点label赋值，label计数+1
			connectedLabels(curVerIdx) = currentLabel;
			connectedCount(currentLabel)++;

			// 2.2.2 在邻接矩阵中搜索当前顶点邻接的顶点：
			for (typename Eigen::SparseMatrix<Index>::InnerIterator iter(adjSM, curVerIdx); iter; ++iter)
			{
				const Index connectVerIdx = iter.row();					// 默认稀疏矩阵列优先存储；当前迭代器在第curVerIdx列；

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


// 提取三角网格中最大单连通区域
template <typename DerivedV, typename DerivedI>
bool simplyConnectedLargest(const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris, \
	Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();

	// 1. 生成邻接矩阵：
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx, adjSM;
	adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
	adjSM = adjSM_eCount;
	traverseSparseMatrix(adjSM, [&adjSM](auto& iter) 
		{
			if (iter.value() > 0)
				iter.valueRef() = 1;
		});

	// 2. 确定网格中所有单连通区域：
	Eigen::VectorXi connectedLabels, connectedCount;
	int conCount = simplyConnectedRegion(adjSM, connectedLabels, connectedCount);
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


// 去除三角片：
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


// 检测网格中的孤立顶点
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


// 去除网格中的孤立顶点：
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

	// 输出点云：
	subFromIdxVec(versOut, vers, newOldIdxInfo);

	// 输出三角片：
	trisOut = tris;
	int* dataPtr = trisOut.data();
	for (int i = 0; i < trisCount * 3; ++i)
		*dataPtr++ = oldNewIdxInfo[*dataPtr];

	return true;
}


// 检测退化边（边长过短的边）
template <typename DerivedV>
Eigen::VectorXi checkDegEdges(const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
		const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const double threshold = 1e-3)
{
	/*
		Eigen::VectorXi checkDegEdges(														// 返回标记退化边的flag向量
				const Eigen::MatrixXi& edges,												// 网格有向边
				const Eigen::PlainObjectBase<DerivedV>& edgeArrows,		// 有向边的边向量数据
				const Eigen::PlainObjectBase<DerivedV>& vers, 
				const Eigen::MatrixXi& tris, 
				const double threshold = 1e-3												// 判定退化边的边长阈值；
				)
	*/
	Eigen::VectorXd edgesLen = edgeArrows.rowwise().norm();
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	unsigned edgesCount = edges.rows();

	Eigen::VectorXi flags = (edgesLen.array() <= threshold).select(Eigen::VectorXi::Ones(edgesCount), Eigen::VectorXi::Zero(edgesCount));
	
	return flags;
}


// 去除网格退化边（边长过短的边，将两个顶点融合）;
template <typename DerivedV>
int mergeDegEdges(Eigen::PlainObjectBase<DerivedV>& newVers, Eigen::MatrixXi& newTris, \
	const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const Eigen::VectorXi& degEdgeFlags)
{
	/*
		int mergeDegEdges(																				返回去除掉的顶点数，若失败则返回-1；
				Eigen::PlainObjectBase<DerivedV>& newVers, 
				Eigen::MatrixXi& newTris, 
				const Eigen::MatrixXi& edges, 
				const Eigen::PlainObjectBase<DerivedV>& edgeArrows, 
				const Eigen::PlainObjectBase<DerivedV>& vers, 
				const Eigen::MatrixXi& tris,	
				const Eigen::VectorXi& degEdgeFlags								标记退化边的flag向量
				)
	
	*/
	int repVersCount = -1;
	int* dataPtr = nullptr;
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	unsigned edgesCount = edges.rows();

	// 1. 提取退化有向边
	int degCount = degEdgeFlags.sum();				// 退化有向边的数量
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

	// 所有退化边中，保留索引较小的顶点，索引大的顶点改写为索引小的顶点：
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
			else       // 当前边的顶点与之前的簇中的顶点相同，插入到该簇中；
			{
				clusterSet.insert(vaIdx);
				clusterSet.insert(vbIdx);
				break;
			}
		}

		// 生成新的簇；
		std::set<int> tmpCluster;
		tmpCluster.insert(vaIdx);
		tmpCluster.insert(vbIdx);
		clusters.push_back(tmpCluster);
	}

	// 生成顶点聚类的映射字典——一个簇中所有其他顶点都映射为索引最小的那个顶点：
	std::map<int, int> clusterMap;
	for (const auto clusterSet: clusters) 
	{
		auto iter = clusterSet.begin();
		int headIdx = *iter;
		iter++;
		for (; iter != clusterSet.end(); ++iter)
			clusterMap.insert({*iter, headIdx});
	}

	std::vector<int> repIdxes;							// 需要去除的顶点的索引；
	repIdxes.reserve(degCount);
	for (const auto& pair : clusterMap)
		repIdxes.push_back(pair.second);
	repIdxes.shrink_to_fit();
	repVersCount = repIdxes.size();
 	
	std::map<int, int> fullMap = clusterMap;
	for (int i = 0; i < versCount; ++i)
		fullMap.insert({i, i});

	// 聚类映射：
	Eigen::MatrixXi trisCopy = tris;
	dataPtr = trisCopy.data();
	for (unsigned i = 0;  i < 3 * trisCount; ++i ) 
	{
		int index = *dataPtr;
		*dataPtr = fullMap[index];
		dataPtr++;
	}

	// 删除非法三角片：
	std::vector<unsigned> sickIdxes;
	checkSickTris(sickIdxes, trisCopy);
	removeTris(newTris, trisCopy, sickIdxes);
 
	// 删除孤立顶点：
	trisCopy = newTris;
	std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, trisCopy);
	newVers.resize(0, 0);
	newTris.resize(0, 0);
	removeIsoVers(newVers, newTris, vers, trisCopy, isoVerIdxes);

	return degCount;
}


// 生成栅格点云
template<typename T>
bool genGrids(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, const Eigen::Matrix<T, 1, 3>& origin, \
		const float step, const std::vector<unsigned>& gridCounts)
{
	/*
		bool genGrids(
			Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& gridCenters,			// 输出的栅格点云
			const Eigen::Matrix<T, 1, 3>& origin,														// 栅格原点，即三轴坐标最小的那个栅格点；
			const float step,																						// 采样步长
			const std::vector<unsigned>& gridCounts											// 三元数组，分别是XYZ三轴上的步数；
			)	
	*/

	using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

	// 整个距离场的包围盒：
	RowVector3T minp = origin;
	RowVector3T maxp = origin + step * RowVector3T(gridCounts[0] - 1, gridCounts[1] - 1, gridCounts[2] - 1);

	// 生成栅格：
	/*
		按索引增大排列的栅格中心点为：
		gc(000), gc(100), gc(200), gc(300),...... gc(010), gc(110), gc(210), gc(310),...... gc(001), gc(101), gc(201).....

		x坐标：
		x0, x1, x2, x3......x0, x1, x2, x3......x0, x1, x2, x3......
		周期为xCount;
		重复次数为(yCount * zCount)

		y坐标：
		y0, y0, y0...y1, y1, y1...y2, y2, y2.........y0, y0, y0...y1, y1, y1...
		周期为(xCount * yCount);
		重复次数为zCount;
		单个元素重复次数为xCount

		z坐标：
		z0, z0, z0......z1, z1, z1......z2, z2, z2......
		单个元素重复次数为(xCount * yCount)
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


// 去除非法三角片，重复三角片：
template<typename DerivedV, typename DerivedI>
int removeSickDupTris(const Eigen::PlainObjectBase<DerivedV>& vers, Eigen::PlainObjectBase<DerivedI>& tris)
{
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();

	Eigen::VectorXi flags{Eigen::VectorXi::Ones(trisCount)};			// flag向量，需要去除的三角片标0；
	std::unordered_set<std::uint64_t> triCodeSet;

	// 遍历三角片
	for (unsigned i = 0; i < trisCount; ++i) 
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
		if(!retPair.second)
			flags(i) = 0;				// 重复三角片非法；
	}

	DerivedI tmpMat;
	subFromFlagVec(tmpMat, tris, flags);
	tris = tmpMat;

	int removeCount = trisCount - tris.rows();
	return removeCount;
}



//////////////////////////////////////////////////////////////////////////////////////////////// 图像处理相关：

// 矩阵的空域线性滤波——！！！注：滤波mask尺寸必须为奇数！！！
template<typename T>
bool linearSpatialFilter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matOut, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matIn, \
	const Eigen::MatrixXd& mask)
{
	unsigned rows = matIn.rows();
	unsigned cols = matIn.cols();
	matOut.resize(rows, cols);
	matOut.setZero(rows, cols);

	unsigned sizeMask = mask.rows();
	unsigned im = (sizeMask + 1) / 2;					// 掩膜中心的下标；
	unsigned km = im;
	unsigned offset = (sizeMask - 1) / 2;

	// 1. 生成拓展矩阵：
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matExt = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows + 2*offset, cols + 2 * offset);
	matExt.block(offset, offset, rows, cols) = matIn;

	// 2. 滑动领域操作：
	PARALLEL_FOR(0, rows, [&](int i)
		{
			for (unsigned k = 0; k < cols; ++k)
			{
				unsigned ie = i + offset;					// 此时mask中心元素在matExt中的位置下标；
				unsigned ke = k + offset;
				Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> coveredElems = matExt.block(ie - offset, ke - offset, sizeMask, sizeMask);
				Eigen::MatrixXd tmpMat = coveredElems.array().cast<double>().array() * mask.array();
				matOut(i, k) = static_cast<T>(tmpMat.sum());
			}
		});
 

	return true;
}
 



//////////////////////////////////////////////////////////////////////////////////////////////////TMESH相关接口：
 
// 对TMESH顶点的只读遍历；
template <typename F>
void traverseVersList(const T_MESH::List& list, F f) 
{
	const T_MESH::Vertex* verPtr = nullptr;					// 底层const指针
	const T_MESH::Node* nodePtr = nullptr;
	const T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), verPtr = nodePtr ? ((const T_MESH::Vertex*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), verPtr = (nodePtr) ? ((const T_MESH::Vertex*)nodePtr->data) : nullptr)
	{
		f(verPtr);
	}
}


// 对TMESH顶点的非只读遍历；
template <typename F>
void traverseVersList(T_MESH::List& list, F f)
{
	T_MESH::Vertex* verPtr = nullptr;					 
	T_MESH::Node* nodePtr = nullptr;
	T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), verPtr = nodePtr ? ((T_MESH::Vertex*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), verPtr = (nodePtr) ? ((T_MESH::Vertex*)nodePtr->data) : nullptr)
	{
		f(verPtr);
	}
}



template <typename F>
void traverseEdgesList(const T_MESH::List& list, F f)
{
	const T_MESH::Edge* ePtr = nullptr;
	const T_MESH::Node* nodePtr = nullptr;
	const T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), ePtr = nodePtr ? ((const T_MESH::Edge*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), ePtr = (nodePtr) ? ((const T_MESH::Edge*)nodePtr->data) : nullptr)
	{
		f(ePtr);
	}
}

template <typename F>
void traverseEdgesList(T_MESH::List& list, F f)
{
	T_MESH::Edge* ePtr = nullptr;
	T_MESH::Node* nodePtr = nullptr;
	T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), ePtr = nodePtr ? ((T_MESH::Edge*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), ePtr = (nodePtr) ? ((T_MESH::Edge*)nodePtr->data) : nullptr)
	{
		f(ePtr);
	}
}


template <typename F>
void traverseTrisList(const T_MESH::List& list, F f)
{
	const T_MESH::Triangle* triPtr = nullptr;
	const T_MESH::Node* nodePtr = nullptr;
	const T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), triPtr = nodePtr ? ((const T_MESH::Triangle*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), triPtr = (nodePtr) ? ((const T_MESH::Triangle*)nodePtr->data) : nullptr)
	{
		f(triPtr);
	}
}


template <typename F>
void traverseTrisList(T_MESH::List& list, F f)
{
	T_MESH::Triangle* triPtr = nullptr;
	T_MESH::Node* nodePtr = nullptr;
	T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), triPtr = nodePtr ? ((T_MESH::Triangle*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), triPtr = (nodePtr) ? ((T_MESH::Triangle*)nodePtr->data) : nullptr)
	{
		f(triPtr);
	}
}


// TMESH网格转换为矩阵： 
template <typename T>
void TMesh2MeshMat(T_MESH::Basic_TMesh& mesh, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris)
{
	const int versCount = mesh.V.numels();
	const int trisCount = mesh.T.numels();
	int index = 0;
	vers.resize(versCount, 3);
	tris.resize(trisCount, 3);

	// 1. 生成顶点矩阵：
	traverseVersList(mesh.V, [&](T_MESH::Vertex* verPtr)
		{
			vers(index, 0) = verPtr->x;
			vers(index, 1) = verPtr->y;
			vers(index, 2) = verPtr->z;
			index++;
		});

	// 2. 顶点节点的x坐标暂存入xValues中，然后将其改写为顶点索引；
	std::vector<double> xValues;
	xValues.reserve(versCount);

	index = 0;
	traverseVersList(mesh.V, [&](T_MESH::Vertex* verPtr)
		{
			xValues.push_back(verPtr->x);
			verPtr->x = static_cast<double>(index++);
		});
 

	// 3. 生成三角片数据，存入矩阵：
	index = 0;
	traverseTrisList(mesh.T, [&](T_MESH::Triangle* triPtr)
		{
			tris(index, 0) = static_cast<int>(triPtr->v1()->x);
			tris(index, 1) = static_cast<int>(triPtr->v2()->x);
			tris(index, 2) = static_cast<int>(triPtr->v3()->x);
			index++;
		});

	// 4. 顶点节点中的x数据还原：
	index = 0;
	traverseVersList(mesh.V, [&](T_MESH::Vertex* verPtr)
		{
			verPtr->x = xValues[index++];
		});
}


// 计算四个3D点围成的四面体的符号体积×6，或以v2,v3,v4围成的三角形为底面的棱柱的符号体积；返回正值：v1在三角形正侧；负值：v1在三角形负侧；0：四点共面；
/*
	Returns a positive value if the point d lies above the plane passing through a, b, and c, meaning that a, b,
			and c appear in counterclockwise order when viewed from d.

	Returns a negative value if d lies below the plane.

	Returns zero if the points are coplanar.

	The result is also an approximation of six times the signed volume of the tetrahedron defined by the four points.

*/
template <typename T>
double orient3D(const Eigen::Matrix<T, 1, 3>& v1, const Eigen::Matrix<T, 1, 3>& v2, const Eigen::Matrix<T, 1, 3>& v3, const Eigen::Matrix<T, 1, 3>& v4)
{
	std::vector<double> p1{ v1(0), v1(1), v1(2) };
	std::vector<double> p2{ v2(0), v2(1), v2(2) };
	std::vector<double> p3{ v3(0), v3(1), v3(2) };
	std::vector<double> p4{ v4(0), v4(1), v4(2) };
	return orient3d(&p1[0], &p2[0], &p3[0], &p4[0]);
}


// 计算三个2D点围成三角形的符号面积×2，可三个点的位置关系——正数: 逆时针； 负数：顺时针； 0：共线； 
/*
	Returns a positive value if the points a, b, and c occur in counterclockwise order
			(c lies to the left of the directed line defined by points a and b).

	Returns a negative value if they occur in clockwise order (c lies to the right of the directed line ab).

	Returns zero if they are collinear.

	The result is also an approximation of twice the signed area of the triangle defined by the three points.
*/
template <typename T>
double orient2D(const Eigen::Matrix<T, 1, 2>& v1, const Eigen::Matrix<T, 1, 2>& v2, const Eigen::Matrix<T, 1, 2>& v3)
{
	std::vector<double> p1{ v1(0), v1(1) };
	std::vector<double> p2{ v2(0), v2(1) };
	std::vector<double> p3{ v3(0), v3(1) };
	return orient2d(&p1[0], &p2[0], &p3[0]);
}


// 生成cell对应的轴向包围盒网格；
void genAABBmesh(const T_MESH::di_cell& cell, Eigen::MatrixXd& vers, Eigen::MatrixXi& tris);



///////////////////////////////////////////////////////////////////////////////////////// debug 接口：

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
	void test2();
	void test3();
	void test4();
	void test5();
	void test6();
	void test7();
	void test8();
	void test9();
	void test10();
	void test11();
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


// 测试拓扑网格类tmesh
namespace TEST_TMESH
{
	void test0();
	void test1();
	void test2();
	void test3();
	void test4();
	void test5();
	void test6();
	void test7();
	void test8();
	void test9();
	void test10();
	void test11();
	void test12();
	void test13();
	void test14();
	void test15();
	void test16();
}


