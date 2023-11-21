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

// 传入函数子遍历稠密矩阵中的元素，函数子接受的参数是类型为typename Derived::Scalar的矩阵元素对象
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
			const Eigen::MatrixBase<DerivedVA>& axisArrow,				旋转轴向量
			const double theta																逆时针旋转角度
			)


		注！！！计算出的旋转矩阵全部默认为作用于列向量；
				若v,u为列向量，r,s为行向量：R * v == u; 等价于 r * R.transpose() == s;
*/
template <typename ScalarVO, typename DerivedVA>
void getRotationMat(Eigen::Matrix<ScalarVO, 3, 3>& rotation, \
	const Eigen::MatrixBase<DerivedVA>& axisArrow, const double theta)
{
	assert((3 == axisArrow.rows() * axisArrow.cols()) && "assert!!! the input axisArrow should be in 3D space.");
	Eigen::Matrix<ScalarVO, 3, 1> axis = axisArrow.transpose().normalized().array().cast<ScalarVO>();
	rotation = Eigen::AngleAxis<ScalarVO>(theta, axis).toRotationMatrix();
}


// 输入旋转信息得到旋转矩阵，重载2——得到将originArrow旋转到targetArrow的旋转矩阵，仅限3D情形；
/*
	template <typename ScalarVO, typename ScalarV1, typename ScalarV2>
	void getRotationMat(
			Eigen::Matrix<ScalarVO, 3, 3>& rotation, \						 输出的旋转矩阵，(3×3)
			const Eigen::MatrixBase<DerivedV1>& originArrow,			 旋转前的列向量
			const Eigen::MatrixBase<DerivedV2>& targetArrow			 旋转后的列向量；
			)
		注！！！计算出的旋转矩阵全部默认为作用于列向量；
				若v,u为列向量，r,s为行向量：R * v == u; 等价于 r * R.transpose() == s;
*/
template <typename ScalarVO, typename DerivedV1, typename DerivedV2>
void getRotationMat(Eigen::Matrix<ScalarVO, 3, 3>& rotation, \
	const Eigen::MatrixBase<DerivedV1>& originArrow, \
	const Eigen::MatrixBase<DerivedV2>& targetArrow)
{
	using RowVector3O = Eigen::Matrix<ScalarVO, 1, 3>;
	const double eps = 1e-5;
	assert((3 == originArrow.rows() * originArrow.cols()) && "assert!!! the input originArrow should be in 3D space.");
	assert((3 == targetArrow.rows() * targetArrow.cols()) && "assert!!! the input targetArrow should be in 3D space.");
	rotation = Eigen::Matrix<ScalarVO, 3, 3>::Identity();
	if (0 == originArrow.norm() || 0 == targetArrow.norm())
		return;

	// 1. 若原方向与目标方向平行：
	RowVector3O originArrowCast{ static_cast<ScalarVO>(originArrow(0)), \
			static_cast<ScalarVO>(originArrow(1)), static_cast<ScalarVO>(originArrow(2))};
	RowVector3O targetArrowCast{ static_cast<ScalarVO>(targetArrow(0)), \
			static_cast<ScalarVO>(targetArrow(1)), static_cast<ScalarVO>(targetArrow(2)) };
	RowVector3O axisArrowCast = originArrowCast.cross(targetArrowCast);			// 旋转轴；
	if (std::abs(axisArrowCast.norm()) < eps)
	{
		// 1.1 若原方向与目标方向相同：
		if (originArrowCast.dot(targetArrowCast) > 0)
			return;
		else         // 1.2 若原方向与目标方向相反：
		{
			// 求以originArrow为法向的平面内的一个向量：
			axisArrowCast = RowVector3O{ 1, 1, 1 };
			ScalarVO x0 = originArrowCast(0);
			ScalarVO y0 = originArrowCast(1);
			ScalarVO z0 = originArrowCast(2);

			// 若向量o最大分量为x分量，则设定a向量的另外两个分量为1, x分量待定，求解o.dot(a)解出x分量，再归一化得到旋转轴向量a
			if (std::abs(x0) >= std::abs(y0) && std::abs(x0) >= std::abs(z0))
				axisArrowCast(0) = -(y0 + z0) / x0;
			else if (std::abs(y0) >= std::abs(x0) && std::abs(y0) >= std::abs(z0))
				axisArrowCast(1) = -(x0 + z0) / y0;
			else
				axisArrowCast(2) = -(x0 + y0) / z0;
			axisArrowCast.normalize();
			ScalarVO x1 = axisArrowCast(0);
			ScalarVO y1 = axisArrowCast(1);
			ScalarVO z1 = axisArrowCast(2);

			rotation << 2 * x1 * x1 - 1, 2 * x1 * y1, 2 * x1 * z1, \
				2 * x1 * y1, 2 * y1 * y1 - 1, 2 * y1 * z1, \
				2 * x1 * z1, 2 * y1 * z1, 2 * z1 * z1 - 1;

			return;
		}
	}

	// 2. 原方向与目标方向不平行的一般情形：
	axisArrowCast.normalize();
	ScalarVO x0 = axisArrowCast(0);
	ScalarVO y0 = axisArrowCast(1);
	ScalarVO z0 = axisArrowCast(2);
	ScalarVO cosTheta = originArrowCast.dot(targetArrowCast) / (originArrowCast.norm() * targetArrowCast.norm());
	ScalarVO sinTheta = std::sqrt(1 - cosTheta * cosTheta);

	//		note: 等价于Eigen::AngleAxis<T>(theta, axis).toRotationMatrix()，计算绕任意轴向量旋转theta角度；
	rotation << 	cosTheta + (1 - cosTheta) * x0 * x0, (1 - cosTheta)* x0* y0 - sinTheta * z0, (1 - cosTheta)* x0* z0 + sinTheta * y0, \
		(1 - cosTheta)* y0* x0 + sinTheta * z0, cosTheta + (1 - cosTheta) * y0 * y0, (1 - cosTheta)* y0* z0 - sinTheta * x0, \
		(1 - cosTheta)* z0* x0 - sinTheta * y0, (1 - cosTheta)* z0* y0 + sinTheta * x0, cosTheta + (1 - cosTheta) * z0 * z0;
}


// 按照索引向量从输入矩阵中提取元素，重载1——索引向量为Eigen向量
template <typename DerivedO, typename DerivedI, typename DerivedVI>
bool subFromIdxVec(Eigen::PlainObjectBase<DerivedO>& matOut, \
	const Eigen::MatrixBase<DerivedI>& matIn, const Eigen::MatrixBase<DerivedVI>& vec)
{ 
	using ScalarO = typename DerivedO::Scalar;
	int maxIdx = static_cast<int>(vec.maxCoeff());
	assert((maxIdx < matIn.rows()) && "Assesrt!!! index out of range in matIn.");
	assert((1 == vec.rows() || 1 == vec.cols()) && "Assesrt!!! input vec is not a vector.");

	matOut.resize(0, 0);
	matOut.resize(vec.size(), matIn.cols());
	for (unsigned i = 0; i < vec.rows(); ++i)
	{
		const int& index = vec(i);
		matOut.row(i) = matIn.row(index).array().cast<ScalarO>();
	}

	return true;
}


// 按照索引向量从输入矩阵中提取元素，重载2——索引向量为std::vector<>
template <typename DerivedO, typename DerivedI, typename Index>
bool subFromIdxVec(Eigen::PlainObjectBase<DerivedO>& matOut,\
	const Eigen::MatrixBase<DerivedI>& matIn, const std::vector<Index>& vec)
{ 
	using ScalarO = typename DerivedO::Scalar;
	assert((!vec.empty()) && "Assert!!! input vec should not be empty.");
	int maxIdx = static_cast<int>(*std::max_element(vec.begin(), vec.end()));
	assert((maxIdx < matIn.rows()) && "Assesrt!!! index out of range in matIn.");

	matOut.resize(vec.size(), matIn.cols());
	for (unsigned i = 0; i < vec.size(); ++i)
	{
		const Eigen::Index& index = vec[i];
		matOut.row(i) = matIn.row(index).array().cast<ScalarO>();
	}

	return true;
}


// 按照索引容器（只能是STL容器）从输入矩阵中提取元素
template <typename DerivedO, typename DerivedI, typename IndexContainer>
bool subFromIdxCon(Eigen::PlainObjectBase<DerivedO>& matOut, \
	const Eigen::MatrixBase<DerivedI>& matIn, const IndexContainer& con)
{
	using ScalarO = typename DerivedO::Scalar;
	assert((!con.empty()) && "Assert!!! input vec should not be empty.");
	int maxIdx = static_cast<int>(*std::max_element(con.begin(), con.end()));
	assert((maxIdx < matIn.rows()) && "Assesrt!!! index out of range in matIn.");

	matOut.resize(con.size(), matIn.cols());

	auto iter = con.begin();
	for (unsigned i = 0; iter != con.end(); ++i)
	{
		const Eigen::Index& index = *iter++;
		matOut.row(i) = matIn.row(index).array().cast<ScalarO>();
	}

	return true;
}


// 根据flag向量从源矩阵中提取元素生成输出矩阵，重载1：
template <typename DerivedO, typename DerivedI, typename DerivedVI>
bool subFromFlagVec(Eigen::PlainObjectBase<DerivedO>& matOut, \
	const Eigen::MatrixBase<DerivedI>& matIn, const Eigen::MatrixBase<DerivedVI>& vec)
{ 
	using ScalarO = typename DerivedO::Scalar;
	assert((1 == vec.rows() || 1 == vec.cols()) && "Assesrt!!! input vec is not a vector.");
	assert((vec.size() == matIn.rows()) && "Assert!!! input vec and matIn do not match.");

	matOut.resize(0, 0);
	matOut.resize(vec.sum(), matIn.cols());

	int count = 0;
	for (unsigned i = 0; i < vec.rows(); ++i)
		if (vec(i) > 0)
			matOut.row(count++) = matIn.row(i).array().cast<ScalarO>();

	return true;
}


// 根据flag向量从源矩阵中提取元素生成输出矩阵，重载2：输出
template <typename DerivedO, typename DerivedI, typename DerivedVI>
bool subFromFlagVec(Eigen::PlainObjectBase<DerivedO>& matOut, \
	std::vector<int>& oldNewIdxInfo, std::vector<int>& newOldIdxInfo, \
	const Eigen::MatrixBase<DerivedI>& matIn, const Eigen::MatrixBase<DerivedVI>& vec)
{
	using ScalarO = typename DerivedO::Scalar;
	assert((1 == vec.rows() || 1 == vec.cols()) && "Assesrt!!! input vec is not a vector.");
	assert((vec.size() == matIn.rows()) && "Assert!!! input vec and matIn do not match.");
	matOut.resize(0, 0);

	// 1. 抽取vec标记为1的层数：
	const int N = vec.rows();
	const int M = vec.sum();
	int count = 0; 
	oldNewIdxInfo.clear();
	newOldIdxInfo.clear();
	matOut.resize(M, matIn.cols());
	for (unsigned i = 0; i < vec.rows(); ++i)
		if (vec(i) > 0)
			matOut.row(count++) = matIn.row(i).array().cast<ScalarO>();

	// 2. 生成新老索引映射表：
	oldNewIdxInfo.resize(N, -1);
	newOldIdxInfo.reserve(M);
	int index = 0;
	for (int k = 0; k < N; ++k)
	{
		int oldIdx = k;
		if (vec(oldIdx) > 0)
		{
			int newIdx = index++;
			oldNewIdxInfo[oldIdx] = newIdx;
			newOldIdxInfo.push_back(oldIdx);
		}
	}

	return true;
}
 

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
bool vecInsertVec(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec1, \
		const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec2);
 

// 矩阵末尾插入矩阵/行向量：
template <typename Derived1, typename Derived2>
bool matInsertRows(Eigen::PlainObjectBase<Derived1>& mat, \
	const Eigen::MatrixBase<Derived2>& mat1)
{
	using Scalar = typename Derived1::Scalar;
	assert((0 == mat.cols()) || (mat.cols() == mat1.cols()) && "Assert!!! Matrix size not match.");
	const int cols = mat1.cols();
	const int currentRows = mat.rows();
	const int addRows = mat1.rows();
	mat.conservativeResize(currentRows + addRows, cols);
	mat.bottomRows(addRows) = mat1.array().cast<Scalar>();

	return true;
}


// 矩阵右边插入矩阵/列向量：
template <typename Derived1, typename Derived2>
bool matInsertCols(Eigen::PlainObjectBase<Derived1>& mat, \
	const Eigen::MatrixBase<Derived2>& mat1)
{
	using Scalar = typename Derived1::Scalar;
	assert((0 == mat.rows()) || (mat.rows() == mat1.rows()), "Assert!!! Matrix size not match.");
	const int rows = mat1.rows();
	const int currentCols = mat.cols();
	const int addCols = mat1.cols();
	mat.conservativeResize(rows, currentCols + addCols);
	mat.rightCols(addCols) = mat1.array().cast<Scalar>();

	return true;
}


template <typename Derived, typename Scalar, int N>
bool matInsertCols(Eigen::PlainObjectBase<Derived>& mat, \
	const Eigen::Matrix<Scalar, N, 1>& vec)
{
	using ScalarO = typename Derived::Scalar;
	assert((0 == mat.rows()) || (mat.rows() == vec.rows()), "Error!!! Matrix size not match.");
	unsigned rows = vec.rows();
	unsigned currentCols = mat.cols();
	mat.conservativeResize(rows, currentCols + 1);
	mat.col(currentCols) = vec.array().cast<ScalarO>();

	return true;
}


// 输入flag向量，得到新老索引的映射关系； 
/*
	void flagVec2oldNewIdxInfo(
			std::vector<int>& oldNewIdxInfo,
			std::vector<int>& newOldIdxInfo,
			const Eigen::MatrixBase<DerivedVI>& flagVec
			)

	flagVec中，若1 == flagVec(i)表示索引为i的点被选中，若0 == flagVec(j)表示索引为j的点没有被选中；
	-1 == oldNewIdxInfo(i)表示原点云中索引为i的点没有被选中，在新点云中没有对应；
	index0 == newOldIdxInfo(i)表示新点云中索引为i的点，在原点云中的索引为index0;
*/
template <typename DerivedVI>
void flagVec2oldNewIdxInfo(std::vector<int>& oldNewIdxInfo, \
		std::vector<int>& newOldIdxInfo, const Eigen::MatrixBase<DerivedVI>& flagVec) 
{
	const int versCount = flagVec.size();
	assert((1 == flagVec.rows() || 1 == flagVec.cols()) && "assert!!! input flagVec is not a Eigen vector.");
	const int versCountNew = static_cast<int>(flagVec.sum());
	assert((versCount > 0) && "assert!!! input flagVec should not be empty.");
	oldNewIdxInfo.clear();
	newOldIdxInfo.clear();

	oldNewIdxInfo.resize(versCount, -1);
	newOldIdxInfo.reserve(versCount);
	int index = 0;
	for (int i = 0; i < versCount; ++i)
	{
		if (static_cast<int>(flagVec(i)) > 0)
		{
			oldNewIdxInfo[i] = index++;
			newOldIdxInfo.push_back(i);
		}
	}
	newOldIdxInfo.shrink_to_fit();
}


// 返回一个flag列向量retVec，若mat的第i行和行向量vec相等，则retVec(i)==1，否则等于0；
template <typename Derived1, typename Derived2>
Eigen::VectorXi rowInMat(const Eigen::MatrixBase<Derived1>& mat, \
	const Eigen::MatrixBase<Derived2>& rowVec)
{
	using T = typename Derived1::Scalar;
	const int rows = mat.rows();
	const int cols = mat.cols();
	Eigen::VectorXi retVec(rows);
	assert(rowVec.cols() == cols, "Error!!! Tow mats size do not match.");

	// 逐列比较：
	Eigen::MatrixXi tempMat(rows, cols);
	for (int i = 0; i < cols; ++i)
	{
		T value = static_cast<T>(rowVec(i));
		tempMat.col(i) = (mat.col(i).array() == value).select(\
				Eigen::VectorXi::Ones(rows), Eigen::VectorXi::Zero(rows));
	}

	retVec = tempMat.col(0);
	if (cols > 1)
		for (int i = 1; i < cols; ++i)
			retVec = retVec.array() * tempMat.col(i).array();			// 逐列相乘：

	return retVec;
}


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
 

///////////////////////////////////////////////////////////////////////////////////////////////////////// 滤波：

// 矩阵的空域线性滤波——！！！注：滤波mask尺寸必须为奇数！！！
template<typename T>
bool linearSpatialFilter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matIn, \
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

	// 1. 对输入矩阵进行边缘延拓，生成拓展矩阵：
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matExt = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows + 2 * offset, cols + 2 * offset);
	matExt.block(offset, offset, rows, cols) = matIn;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rowHead = matIn.block(0, 0, offset, cols);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rowTail = matIn.block(rows - offset, 0, offset, cols);

	matExt.block(0, offset, offset, cols) = rowHead;
	matExt.block(rows + offset, offset, offset, cols) = rowTail;

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> colHead = matExt.block(0, offset, rows + 2 * offset, offset);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> colTail = matExt.block(0, cols, rows + 2 * offset, offset);
	matExt.block(0, 0, rows + 2 * offset, offset) = colHead;
	matExt.block(0, cols + offset, rows + 2 * offset, offset) = colTail;

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

#include "tmp.tpp"