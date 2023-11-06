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

const double pi = 3.14159265359;

// Index:
/*
	auxiliary interfaces
	basic math tools
	矩阵的增删查改和运算
	拟合、插值

*/


/////////////////////////////////////////////////////////////////////////////////////////////////// auxiliary interfaces:

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
void PARALLEL_FOR(unsigned int  beg, unsigned int  end, const Func& func, \
	const paramTuple& pt, const unsigned int serial_if_less_than = 12)
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




/////////////////////////////////////////////////////////////////////////////////////////////////// basic math tools:

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


// 二维笛卡尔坐标转换为极坐标；
template <typename T>
std::pair<double, double> cart2polar(const T x, const T y);


// 极坐标转换为二维笛卡尔坐标：
template <typename T>
std::pair<double, double> polar2cart(const T radius, const T theta);


// 霍纳方法（秦九昭算法）求多项式的值
template<typename T, int N>
double hornersPoly(const Eigen::Matrix<T, N, 1>& coeffs, const double x);


// 求多项式的一阶微分
template<typename T, int N>
float polyDiff(const Eigen::Matrix<T, N, 1>& coeffs, const float x);


// eigen向量转换为std向量, 重载1
template<typename T, typename Derived>
void eigenVec2Vec(std::vector<T>& vOut, const Eigen::PlainObjectBase<Derived>& vIn);


// eigen向量转换为std向量，重载2
template<typename Derived>
std::vector<typename Derived::Scalar> eigenVec2Vec(const Eigen::PlainObjectBase<Derived>& vIn);


// std::向量转换为eigen列向量， 重载1
template<typename Derived, typename T>
void vec2EigenVec(Eigen::PlainObjectBase<Derived>& vOut, const std::vector<T>& vIn);


// std::向量转换为eigen列向量， 重载2；
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<T>& vIn);


// 输入旋转信息得到旋转矩阵，重载1——输入旋转轴向量，旋转角度，返回旋转矩阵：
/*
	template <typename ScalarVO, typename ScalarVI>
	void getRotationMat(
			Eigen::Matrix<ScalarVO, 3, 3>& rotation,							输出的旋转矩阵
			const Eigen::Matrix<ScalarVI, 1, 3>& axisArrow,				旋转轴向量
			const double theta																逆时针旋转角度
			)


		注！！！计算出的旋转矩阵全部默认为作用于列向量；
				若v,u为列向量，r,s为行向量：R * v == u; 等价于 r * R.transpose() == s;
*/
template <typename ScalarVO, typename DerivedVA>
void getRotationMat(Eigen::Matrix<ScalarVO, 3, 3>& rotation, \
	const Eigen::PlainObjectBase<DerivedVA>& axisArrow, const double theta);


// 输入旋转信息得到旋转矩阵，重载2——得到将originArrow旋转到targetArrow的旋转矩阵
/*
	template <typename ScalarVO, typename ScalarV1, typename ScalarV2>
	void getRotationMat(
			Eigen::Matrix<ScalarVO, 3, 3>& rotation, \
			const Eigen::Matrix<ScalarV1, 1, 3>& originArrow, \
			const Eigen::Matrix<ScalarV2, 1, 3>& targetArrow
			)


		注！！！计算出的旋转矩阵全部默认为作用于列向量；
				若v,u为列向量，r,s为行向量：R * v == u; 等价于 r * R.transpose() == s;
*/
template <typename ScalarVO, typename DerivedV1, typename DerivedV2>
void getRotationMat(Eigen::Matrix<ScalarVO, 3, 3>& rotation, \
	const Eigen::PlainObjectBase<DerivedV1>& originArrow, \
	const Eigen::PlainObjectBase<DerivedV2>& targetArrow);


//
template <typename Derived>
bool subFromIdxVec(Eigen::MatrixBase<Derived>& matBaseOut, \
	const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec);


//
template <typename Derived, typename Index>
bool subFromIdxVec(Eigen::MatrixBase<Derived>& matBaseOut,\
	const Eigen::MatrixBase<Derived>& matBaseIn, const std::vector<Index>& vec);


//
template <typename Derived, typename IndexContainer>
bool subFromIdxCon(Eigen::MatrixBase<Derived>& matBaseOut, \
	const Eigen::MatrixBase<Derived>& matBaseIn, const IndexContainer& con);


//
template <typename Derived>
bool subFromFlagVec(Eigen::MatrixBase<Derived>& matBaseOut, \
	const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec);


//
template <typename Derived>
bool subFromFlagVec(Eigen::MatrixBase<Derived>& matBaseOut, \
	std::vector<int>& oldNewIdxInfo, std::vector<int>& newOldIdxInfo, \
	const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec);
 

//
Eigen::VectorXi flagVec2IdxVec(const Eigen::VectorXi& flag);


//
Eigen::VectorXi IdxVec2FlagVec(const Eigen::VectorXi& idxVec, const unsigned size);



///////////////////////////////////////////////////////////////////////////////////////////////////////// 矩阵的增删查改和运算

// 列向量尾后插入元素
template<typename T>
bool vecInsertNum(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const T num);

// 列向量尾后插入向量
template<typename T>
bool vecInsertVec(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec2);
 

// 矩阵末尾插入矩阵/行向量：
template <typename Derived1, typename Derived2>
bool matInsertRows(Eigen::PlainObjectBase<Derived1>& mat, \
	const Eigen::PlainObjectBase<Derived2>& mat1)
{
	assert((0 == mat.cols()) || (mat.cols() == mat1.cols()), "Error!!! Matrix size not match.");
	unsigned cols = mat1.cols();
	unsigned currentRows = mat.rows();
	unsigned addRows = mat1.rows();
	mat.conservativeResize(currentRows + addRows, cols);
	for (unsigned i = 0; i < addRows; ++i)
		mat.row(currentRows + i) = mat1.row(i);

	return true;
}


template <typename Derived, typename Scalar, int N>
bool matInsertRows(Eigen::PlainObjectBase<Derived>& mat, \
	const Eigen::Matrix<Scalar, 1, N>& rVec)
{
	using ScalarO = typename Derived::Scalar;
	assert((0 == mat.cols()) || (mat.cols() == rVec.cols()), "Error!!! Matrix size not match.");
	unsigned cols = rVec.cols();
	unsigned currentRows = mat.rows(); 
	mat.conservativeResize(currentRows + 1, cols);
	mat.row(currentRows) = rVec.array().cast<ScalarO>();

	return true;
}
 
// 
template <typename Derived1, typename Derived2>
Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Derived1>& mat, const Eigen::PlainObjectBase<Derived2>& rowVec);

// 生成稀疏矩阵，类似于MATLAB中的sparse()
template <class ValueVector, typename T>
void sparse(Eigen::SparseMatrix<T>& SM, const Eigen::VectorXi& I, const Eigen::VectorXi& J,
	const ValueVector& values, const size_t m, const size_t n);

// 稀疏矩阵转置；Eigen::SparseMatrix自带的transpose()方法太垃圾了
template<typename T>
bool spMatTranspose(Eigen::SparseMatrix<T>& smOut, const Eigen::SparseMatrix<T>& smIn);

// 计算稠密矩阵的克罗内克积（Kronecker product）：重载1——可能的参数组合太多了，懒得写模板特化了 
template <typename Derived0, typename Derived1, typename Derived2>
bool kron(Eigen::PlainObjectBase<Derived0>& result, const Eigen::MatrixBase<Derived1>& m11, \
	const Eigen::MatrixBase<Derived2>& m22)
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
	using Scalar = typename Derived0::Scalar;							// 使用从属名称，需要加上"typename"前缀；
	Derived1 m1 = m11.derived();
	Derived2 m2 = m22.derived();;
	unsigned m = m1.rows();
	unsigned n = m1.cols();
	unsigned p = m2.rows();
	unsigned q = m2.cols();
	unsigned rows = m * p;
	unsigned cols = n * q;

	result.resize(0, 0);
	result.resize(rows, cols);
	Eigen::VectorXi rowIdxVec = Eigen::VectorXi::LinSpaced(0, m - 1, m - 1);
	auto calcKronCol = [&](const unsigned colIdx, const std::tuple<unsigned>& pt)
	{
		unsigned rowIdx0 = std::get<0>(pt);
		result.block(rowIdx0 * p, colIdx * q, p, q) = \
			static_cast<Scalar>(m1(rowIdx0, colIdx)) * m2.array().cast<Scalar>();		// 转换成输出矩阵中的数据类型；
	};

	auto calcKronRow = [&](const unsigned rowIdx)
	{
		std::tuple<unsigned> pt0 = std::make_tuple(rowIdx);
		PARALLEL_FOR(0, n, calcKronCol, pt0);
	};

	PARALLEL_FOR(0, m, calcKronRow);

	return true;
}


// 计算稠密矩阵的克罗内克积（Kronecker product）：重载2
template <typename Derived1, typename Derived2>
Eigen::MatrixXd kron(const Eigen::MatrixBase<Derived1>& m11, \
	const Eigen::MatrixBase<Derived2>& m22)
{ 
	Eigen::MatrixXd result;
	kron(result, m11, m22);

	return result;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////// 线性方程组：

// 解线性方程组
template<typename T, int N>
bool solveLinearEquation(Eigen::Matrix<T, N, 1>& x, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<T, N, 1>& b);

// 解多组线性方程组
template <typename T>
bool solveLinearEquations(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& X, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& B);

// 通过SVD求稠密矩阵的条件数：
template<typename T>
double calcCondNum(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m);



///////////////////////////////////////////////////////////////////////////////////////////////////////// 插值、拟合：

//
void polyInterpolation();

//
void gaussInterpolation();

// 
void leastSquarePolyFitting();

// 岭回归多项式拟合曲线
template <typename T>
void ridgeRegressionPolyFitting(Eigen::VectorXd& theta, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	unsigned m = 4);

// 最小二乘法拟合（逼近）标准椭圆（长短轴和xy坐标轴对齐）
template<typename T>
Eigen::VectorXd fittingStandardEllipse(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& sampleVers);
 

#include "tmp.tpp"