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
#include <numeric>
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
	线性方程组
	拟合、插值
	欧氏几何
	滤波
*/


///////////////////////////////////////////////////////////////////////////////////////////////////// auxiliary types and interfaces
template <typename T>
struct DOUBLET
{
	T x;
	T y;
};

template <typename T>
struct TRIPLET
{
	T x;
	T y;
	T z;
};


// 3D顶点或三角片矩阵转换为TRIPLET向量的形式；
template <typename Derived>
std::vector<TRIPLET<typename Derived::Scalar>> mat2triplets(\
	const Eigen::PlainObjectBase<Derived>& mat)
{
	using Scalar = typename Derived::Scalar;
	using ts = TRIPLET<Scalar>;

	std::vector<ts> vec;
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


template <typename Derived>
TRIPLET<typename Derived::Scalar> vec2triplet(\
	const Eigen::PlainObjectBase<Derived>& vec)
{
	assert((3 == vec.rows() * vec.cols()) && "assert!!! input vec should be in 3D space.");
	using Scalar = typename Derived::Scalar;
	using ts = TRIPLET<Scalar>;
	return ts{ vec(0), vec(1), vec(2) };
}


template <typename Derived>
std::vector<DOUBLET<typename Derived::Scalar>> mat2doublets(\
	const Eigen::PlainObjectBase<Derived>& mat)
{
	using Scalar = typename Derived::Scalar;
	using ds = DOUBLET<Scalar>;

	std::vector<ds> vec;
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


template <typename Derived>
DOUBLET<typename Derived::Scalar> vec2doublet(\
	const Eigen::PlainObjectBase<Derived>& vec)
{
	assert((2 == vec.rows() * vec.cols()) && "assert!!! input vec should be in 2D space.");
	using Scalar = typename Derived::Scalar;
	using ds = DOUBLET<Scalar>;
	return DOUBLET{ vec(0), vec(1) };
}


// 并行for循环
template<typename Func>
void PARALLEL_FOR(unsigned int  beg, unsigned int  end, \
	const Func& func, const unsigned int serial_if_less_than = 12)
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


// Eigen向量转换为std::vector，重载1: 
template<typename T, typename Derived>
void eigenVec2Vec(std::vector<T>& vOut, const Eigen::MatrixBase<Derived>& vIn)
{
	assert(1 == vIn.rows() || 1 == vIn.cols(), "assert!!! Input arg vIn should be a vector.");
	unsigned elemCount = vIn.rows() * vIn.cols();
	vOut.clear();
	vOut.resize(elemCount);
	for (unsigned i = 0; i < elemCount; ++i)
		vOut[i] = static_cast<T>(vIn(i));
}


// Eigen向量转换为std::vector，重载2: 
template<typename Derived>
std::vector<typename Derived::Scalar>  eigenVec2Vec(const Eigen::MatrixBase<Derived>& vIn)
{
	using Scalar = typename Derived::Scalar;
	assert(1 == vIn.rows() || 1 == vIn.cols(), "assert!!! Input arg vIn should be a vector.");

	unsigned elemCount = vIn.rows() * vIn.cols();
	std::vector<Scalar> vOut(elemCount, 1);
	std::memcpy(&vOut[0], &vIn(0), sizeof(Scalar) * elemCount);

	return vOut;
}


// std::向量转换为eigen向量， 重载1；
template<typename Derived, typename T>
void vec2EigenVec(Eigen::PlainObjectBase<Derived>& vOut, const std::vector<T>& vIn)
{
	using Scalar = typename Derived::Scalar;
	unsigned elemCount = vIn.size();
	vOut.resize(0);
	vOut.resize(elemCount);
	for (unsigned i = 0; i < elemCount; ++i)
		vOut(i) = static_cast<Scalar>(vIn[i]);
}


// std::向量转换为eigen向量， 重载2——返回列向量；
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<T>& vIn)
{
	unsigned elemCount = vIn.size();
	Eigen::Matrix<T, Eigen::Dynamic, 1> vOut;
	vOut.resize(elemCount, 1);
	std::memcpy(vOut.data(), &vIn[0], sizeof(T) * elemCount);

	return vOut;
}


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


// 得到齐次坐标系下绕坐标轴逆时针旋转theta弧度后的旋转矩阵；
template <typename T>
Eigen::Matrix<T, 4, 4> getRotationX(const T theta) 
{
	const T cosTheta = cos(theta);
	const T sinTheta = sin(theta);
	Eigen::Matrix<T, 4, 4> rotationX;
	rotationX << 1, 0, 0, 0,  \
						  0, cosTheta, -sinTheta, 0,  \
						  0, sinTheta, cosTheta, 0, \
						  0, 0, 0, 1;
	return rotationX;
}

template <typename T>
Eigen::Matrix<T, 4, 4> getRotationY(const T theta)
{
	const T cosTheta = cos(theta);
	const T sinTheta = sin(theta);
	Eigen::Matrix<T, 4, 4> rotationY;
	rotationY << cosTheta, 0, sinTheta, 0, \
		0, 1, 0, 0, \
		-sinTheta, 0, cosTheta, 0, \
		0, 0, 0, 1;
	return rotationY;
}

template <typename T>
Eigen::Matrix<T, 4, 4> getRotationZ(const T theta)
{
	const T cosTheta = cos(theta);
	const T sinTheta = sin(theta);
	Eigen::Matrix<T, 4, 4> rotationZ;
	rotationZ << cosTheta, -sinTheta, 0, 0, \
		sinTheta, cosTheta, 0, 0, \
		0, 0, 1, 0, \
		0, 0, 0, 1;
	return rotationZ;
}


// 计算局部坐标系下的点云变换到世界坐标系的仿射变换矩阵；
template <typename T, typename DerivedV>
bool getAffineL2G(Eigen::Matrix<T, 4, 4>& affine, \
	const Eigen::MatrixBase<DerivedV>& versLocal, \
	const Eigen::MatrixBase<DerivedV>& versGlobal, \
	const int sampleCount = 10)
{
	/*
	bool getAffineL2G(
		Eigen::Matrix<T, 4, 4>& affine,										输出的齐次坐标系下仿射变换矩阵						
		const Eigen::MatrixBase<DerivedV>& versLocal,			局部坐标系中的点云
		const Eigen::MatrixBase<DerivedV>& versGlobal,			世界坐标系中的点云
		const int sampleCount = 10											计算使用的样本点数目
		)
		样本点数目不可以太多也不可以太少；
				太少的话可能会导致线性方程组的A矩阵不满秩；
				太多的话可能导致A矩阵元素的值上溢出；
	*/
	using Scalar = DerivedV::Scalar;
	using Matrix3S = Eigen::Matrix<Scalar, 3, 3>;
	using Matrix4S = Eigen::Matrix<Scalar, 4, 4>;
	using MatrixXS = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
	using RowVector3S = Eigen::Matrix<Scalar, 1, 3>;
	using Vector3S = Eigen::Matrix<Scalar, 3, 1>;	
	using VectorXS = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
	const double eps = 1e-6;
	affine.setIdentity();

	// 1. 问题建模：
	MatrixXS V, VT, VTV, U, A, versRotated;
	Matrix3S RT;
	RowVector3S center;
	VectorXS u0, u1, u2;
	Vector3S	rt0, rt1, rt2;
	center = versGlobal.colwise().mean();
	versRotated = versGlobal.rowwise() - center;
	V = versLocal.topRows(sampleCount);
	U = versRotated.topRows(sampleCount);
	VT = V.transpose();
	VTV = VT * V;
	u0 = U.col(0);
	u1 = U.col(1);
	u2 = U.col(2);

	// 2. 解线性方程组：(VT * V) * RT = (VT * U);
	Eigen::FullPivLU<MatrixXS> solverLU;			
	if (std::abs(VTV.determinant()) < eps)				// 样本点不够或仿射变换不唯一；
		return false;
	solverLU.compute(VTV);
	rt0 = solverLU.solve(VT * u0);
	rt1 = solverLU.solve(VT * u1);
	rt2 = solverLU.solve(VT * u2);

	// 3. 输出：
	RT.col(0) = rt0;
	RT.col(1) = rt1;
	RT.col(2) = rt2; 
	affine.block(0, 0, 3, 3) = RT.transpose();
	affine(0, 3) = center(0);
	affine(1, 3) = center(1);
	affine(2, 3) = center(2); 

	return true;
}


// 按照索引向量从输入矩阵中提取元素，重载1——索引向量为Eigen向量
template <typename DerivedO, typename DerivedI, typename DerivedVI>
bool subFromIdxVec(Eigen::PlainObjectBase<DerivedO>& matOut, \
	const Eigen::MatrixBase<DerivedI>& matIn, const Eigen::MatrixBase<DerivedVI>& vec)
{
	using ScalarO = typename DerivedO::Scalar;
	matOut.resize(0, 0);
	if (0 == vec.rows())
		return true;

	int maxIdx = static_cast<int>(vec.maxCoeff());
	assert((maxIdx < matIn.rows()) && "Assesrt!!! index out of range in matIn.");
	assert((1 == vec.rows() || 1 == vec.cols()) && "Assesrt!!! input vec is not a vector.");

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
bool subFromIdxVec(Eigen::PlainObjectBase<DerivedO>& matOut, \
	const Eigen::MatrixBase<DerivedI>& matIn, const std::vector<Index>& vec)
{
	using ScalarO = typename DerivedO::Scalar;
	matOut.resize(0, 0);
	if (vec.empty())
		return true;

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
	matOut.resize(0, 0);
	if (con.empty())
		return true;

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
template <typename DerivedO, typename DerivedI, typename TI>
bool subFromFlagVec(Eigen::PlainObjectBase<DerivedO>& matOut, \
	const Eigen::MatrixBase<DerivedI>& matIn, const std::vector<TI>& vec)
{ 
	using ScalarO = typename DerivedO::Scalar;
	assert((vec.size() == matIn.rows()) && "Assert!!! input vec and matIn do not match.");
	matOut.resize(0, 0);

	int rowsOut = static_cast<int>(std::accumulate(vec.begin(), vec.end(), 0));
	int colsOut = matIn.cols(); 
	matOut.resize(rowsOut, colsOut);

	int count = 0;
	for (unsigned i = 0; i < vec.size(); ++i)
		if (vec[i] > 0)
			matOut.row(count++) = matIn.row(i).array().cast<ScalarO>();

	return true;
}


// 根据flag向量从源矩阵中提取元素生成输出矩阵，重载2：输出
template <typename DerivedO, typename DerivedI, typename TI>
bool subFromFlagVec(Eigen::PlainObjectBase<DerivedO>& matOut, \
	std::vector<int>& oldNewIdxInfo, std::vector<int>& newOldIdxInfo, \
	const Eigen::MatrixBase<DerivedI>& matIn, const std::vector<TI>& vec)
{
	using ScalarO = typename DerivedO::Scalar; 
	assert((vec.size() == matIn.rows()) && "Assert!!! input vec and matIn do not match.");
	matOut.resize(0, 0);

	// 1. 抽取vec标记为1的层数：
	const int N = vec.size();
	const int M = static_cast<int>(std::accumulate(vec.begin(), vec.end(), 0));
	int count = 0; 
	oldNewIdxInfo.clear();
	newOldIdxInfo.clear();
	matOut.resize(M, matIn.cols());
	for (unsigned i = 0; i < vec.size(); ++i)
		if (vec[i] > 0)
			matOut.row(count++) = matIn.row(i).array().cast<ScalarO>();

	// 2. 生成新老索引映射表：
	oldNewIdxInfo.resize(N, -1);
	newOldIdxInfo.reserve(M);
	int index = 0;
	for (int k = 0; k < N; ++k)
	{
		int oldIdx = k;
		if (vec[oldIdx] > 0)
		{
			int newIdx = index++;
			oldNewIdxInfo[oldIdx] = newIdx;
			newOldIdxInfo.push_back(oldIdx);
		}
	}

	return true;
}
 

// 根据flag向量从源向量中提取元素生成输出向量
template <typename DerivedO, typename DerivedI, typename TI>
bool subVecFromFlagVec(Eigen::PlainObjectBase<DerivedO>& vecOut, \
	const Eigen::MatrixBase<DerivedI>& vecIn, const std::vector<TI>& vec)
{
	using ScalarO = typename DerivedO::Scalar; 
	assert((vec.size() == vecIn.rows()) && "Assert!!! input vec and vecIn do not match.");
	vecOut.resize(0);

	int rowsOut = static_cast<int>(std::accumulate(vec.begin(), vec.end(), 0));
	vecOut.resize(rowsOut);

	int count = 0;
	for (int i = 0; i < vec.size(); ++i)
		if (vec[i] > 0)
			vecOut.row(count++) = vecIn.row(i).array().cast<ScalarO>();

	return true;
}


//
Eigen::VectorXi flagVec2IdxVec(const Eigen::VectorXi& flag);


//
Eigen::VectorXi IdxVec2FlagVec(const Eigen::VectorXi& idxVec, const unsigned size);


// 三维情形下的主成分分析：
template <typename DerivedV>
bool PCA3d(std::vector<Eigen::RowVector3d>& resultVecs,\
	std::vector<double>& resultVals, const Eigen::MatrixBase<DerivedV>& samples)
{
	assert(3 == samples.cols() && "Assert!!! Samples should be in 3-dimension space.");

	// 1. 所有顶点作为三维样本数据，存入矩阵A中，一行一个样本，
	const int  N = samples.rows();											// 样本容量； 
	Eigen::MatrixXd A = samples.array().cast<double>();
	Eigen::RowVector3d miu = A.colwise().mean();				// 样本均值——miu向量，也是所有顶点的重心；

	// 2. 均值点平移到原点的样本矩阵
	Eigen::MatrixXd B = A - Eigen::MatrixXd::Ones(N, 1) * miu;

	// 3. 求样本的协方差矩阵：
	Eigen::Matrix3d C = B.transpose() * B / (N - 1);

	// 4. 求协方差矩阵的特征值和特征向量
	Eigen::EigenSolver<Eigen::Matrix3d> es(C);
	Eigen::Matrix3d D = es.pseudoEigenvalueMatrix();			// 对角线元素是特征值
	Eigen::Matrix3d V = es.pseudoEigenvectors();					// 每一个列向量都是特征向量。
 
	// 5. 输出：
	std::map<double, Eigen::RowVector3d> tmpMap;
	tmpMap.insert(std::make_pair(D(0, 0), Eigen::RowVector3d{V.col(0).transpose()}));
	tmpMap.insert(std::make_pair(D(1, 1), Eigen::RowVector3d{ V.col(1).transpose() }));
	tmpMap.insert(std::make_pair(D(2, 2), Eigen::RowVector3d{ V.col(2).transpose() }));
	auto iter = tmpMap.begin();
	resultVecs.reserve(3);
	resultVals.reserve(3); 
	for (int i = 0; i < 3; ++i)
	{
		resultVals.push_back(iter->first);
		resultVecs.push_back(iter->second);
		iter++;
	}

	return true;
}




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
	if (0 == mat1.rows())
		return true;

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
	if (0 == mat1.rows())
		return true;

	assert((0 == mat.rows()) || (mat.rows() == mat1.rows()), "Assert!!! Matrix size not match.");
	const int rows = mat1.rows();
	const int currentCols = mat.cols();
	const int addCols = mat1.cols();
	mat.conservativeResize(rows, currentCols + addCols);
	mat.rightCols(addCols) = mat1.array().cast<Scalar>();

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
template <typename TI>
void flagVec2oldNewIdxInfo(std::vector<int>& oldNewIdxInfo, \
		std::vector<int>& newOldIdxInfo, const std::vector<TI>& flagVec) 
{
	const int versCount = flagVec.size(); 
	const int versCountNew = static_cast<int>(std::accumulate(flagVec.begin(), flagVec.end(), 0));
	assert((versCount > 0) && "assert!!! input flagVec should not be empty.");
	oldNewIdxInfo.clear();
	newOldIdxInfo.clear();

	oldNewIdxInfo.resize(versCount, -1);
	newOldIdxInfo.reserve(versCount);
	int index = 0;
	for (int i = 0; i < versCount; ++i)
	{
		if (static_cast<int>(flagVec[i]) > 0)
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

// 解恰定的稠密线性方程组Ax == b;
template<typename DerivedX, typename DerivedA, typename DerivedB>
bool solveLinearEquation(Eigen::PlainObjectBase<DerivedX>& x, const Eigen::MatrixBase<DerivedA>& A, \
	const Eigen::MatrixBase<DerivedB>& b)
{
	using ScalarX = typename DerivedX::Scalar;
	using MatrixXx = Eigen::Matrix<ScalarX, Eigen::Dynamic, Eigen::Dynamic>;
	using VectorXx = Eigen::Matrix<ScalarX, Eigen::Dynamic, 1>;

	// 解线性方程组Ax == b;
	VectorXx b0 = b.array().cast<ScalarX>();
	Eigen::JacobiSVD<MatrixXx> solverSVD(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	x = solverSVD.solve(b0);

	return true;
}


// 解一系列恰定的稠密线性方程组AX == B;
template <typename T>
bool solveLinearEquations(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& X, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& B)
{
	if (A.rows() != B.rows())
		return false;

	Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> solverSVD(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	X.resize(A.cols(), B.cols());
	for (int i = 0; i < B.cols(); ++i)
	{
		Eigen::Matrix < T, Eigen::Dynamic, 1> x = solverSVD.solve(B.col(i));
		X.col(i) = x;
	}

	return true;
}


// 通过SVD求稠密矩阵的条件数：
template<typename T>
double calcCondNum(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m);




///////////////////////////////////////////////////////////////////////////////////////////////////////// 插值、拟合：

// 输入三维点云，输出拟合的平面方程系数(a,b,c,d)，其中(a,b,c)是平面的单位法向量；
template <typename T, typename DerivedV>
bool versFitPlane(std::vector<T>& coeffs, const Eigen::MatrixBase<DerivedV>& versSample) 
{
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Vector4T = Eigen::Matrix<T, 4, 1>;
	assert(versSample.rows() > 0 && "Assert!!! input vertices mat is empty.");
	assert(3 == versSample.cols() && "Assert!!! input vertices should be in 3D space.");
	const int versCount = versSample.rows();
	MatrixXT A(versCount, 4);
	A.leftCols(3) = versSample.cast<T>();
	A.col(3).setOnes();
	
	// A矩阵奇异值分解：
	Eigen::JacobiSVD<MatrixXT> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	MatrixXT V = svd.matrixV();				// 右奇异矩阵；
	Vector4T x = V.col(3);
	double length = std::sqrt(x(0) * x(0) + x(1) * x(1) + x(2)* x(2));
	x.array() /= length;
	eigenVec2Vec(coeffs, x);

	return true;
}


//
void polyInterpolation();

//
void gaussInterpolation();

// 
void leastSquarePolyFitting();


// 岭回归多项式拟合曲线
template <typename T>
void ridgeRegressionPolyFitting(Eigen::VectorXd& theta, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	unsigned m = 4);


// 最小二乘法拟合（逼近）标准椭圆（长短轴和xy坐标轴对齐）
template<typename DerivedV>
Eigen::VectorXd fittingStandardEllipse(const Eigen::MatrixBase<DerivedV>& sampleVers)
{
	/*
		Eigen::VectorXd fittingStandardEllipse(								// 返回列向量(a,c,d,e,f)，为标准椭圆方程的系数；
					const Eigen::MatrixXf& sampleVers					// 输入的样本点，必须是在XOY平面上的点；
		)
	*/
	const double epsilon = 1e-8;							// 浮点数绝对值小于此值时认为为0；
	using T = typename DerivedV::Scalar;
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
 




///////////////////////////////////////////////////////////////////////////////////////////////////////// 欧氏几何

// 输入平面上一点和平面法向，得到平面方程系数；
/*
	平面方程: ax +by+cz+d=0;
			其中(a, b, c)是平面正向法向量，是单位向量；

*/
template <typename T, typename DerivedV>
bool getPlaneCoeff(std::vector<T>& coeffs, const Eigen::MatrixBase<DerivedV>& planeVer, \
	const Eigen::MatrixBase<DerivedV>& planeNorm)
{	
	coeffs.clear(); 
	T a = static_cast<T>(planeNorm(0));
	T b = static_cast<T>(planeNorm(1));
	T c = static_cast<T>(planeNorm(2));
	T d = static_cast<T>(-a * planeVer(0) - b * planeVer(1) - c*planeVer(2));
	coeffs = std::vector<T>{a, b, c, d};
	return true;
}


// 平面到顶点的符号距离：
template<typename T, typename DerivedV>
double disPlane2Vert(const std::vector<T>& coeffs, const Eigen::MatrixBase<DerivedV>& ver)
{
	return static_cast<double>(coeffs[0] * ver(0) + coeffs[1] * ver(1) + coeffs[2] * ver(2) + coeffs[3]);
}

 


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