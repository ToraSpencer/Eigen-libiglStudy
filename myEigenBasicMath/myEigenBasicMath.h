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
	�������ɾ��ĺ�����
	���Է�����
	��ϡ���ֵ
	ŷ�ϼ���
	�˲�
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


// 3D���������Ƭ����ת��ΪTRIPLET��������ʽ��
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


// ����forѭ��
template<typename Func>
void PARALLEL_FOR(unsigned int  beg, unsigned int  end, \
	const Func& func, const unsigned int serial_if_less_than = 12)
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
void PARALLEL_FOR(unsigned int  beg, unsigned int  end, const Func& func, \
	const paramTuple& pt, const unsigned int serial_if_less_than = 12)
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




/////////////////////////////////////////////////////////////////////////////////////////////////// basic math tools:

// ���뺯���ӱ������ܾ����е�Ԫ�أ������ӽ��ܵĲ���������Ϊtypename Derived::Scalar�ľ���Ԫ�ض���
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


// ��ά�ѿ�������ת��Ϊ�����ꣻ
template <typename T>
std::pair<double, double> cart2polar(const T x, const T y);


// ������ת��Ϊ��ά�ѿ������꣺
template <typename T>
std::pair<double, double> polar2cart(const T radius, const T theta);


// ���ɷ������ؾ����㷨�������ʽ��ֵ
template<typename T, int N>
double hornersPoly(const Eigen::Matrix<T, N, 1>& coeffs, const double x);


// �����ʽ��һ��΢��
template<typename T, int N>
float polyDiff(const Eigen::Matrix<T, N, 1>& coeffs, const float x);


// Eigen����ת��Ϊstd::vector������1: 
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


// Eigen����ת��Ϊstd::vector������2: 
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


// std::����ת��Ϊeigen������ ����1��
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


// std::����ת��Ϊeigen������ ����2����������������
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<T>& vIn)
{
	unsigned elemCount = vIn.size();
	Eigen::Matrix<T, Eigen::Dynamic, 1> vOut;
	vOut.resize(elemCount, 1);
	std::memcpy(vOut.data(), &vIn[0], sizeof(T) * elemCount);

	return vOut;
}


// ������ת��Ϣ�õ���ת��������1����������ת����������ת�Ƕȣ�������ת����
/*
	template <typename ScalarVO, typename ScalarVI>
	void getRotationMat(
			Eigen::Matrix<ScalarVO, 3, 3>& rotation,							�������ת����
			const Eigen::MatrixBase<DerivedVA>& axisArrow,				��ת������
			const double theta																��ʱ����ת�Ƕ�
			)


		ע���������������ת����ȫ��Ĭ��Ϊ��������������
				��v,uΪ��������r,sΪ��������R * v == u; �ȼ��� r * R.transpose() == s;
*/
template <typename ScalarVO, typename DerivedVA>
void getRotationMat(Eigen::Matrix<ScalarVO, 3, 3>& rotation, \
	const Eigen::MatrixBase<DerivedVA>& axisArrow, const double theta)
{
	assert((3 == axisArrow.rows() * axisArrow.cols()) && "assert!!! the input axisArrow should be in 3D space.");
	Eigen::Matrix<ScalarVO, 3, 1> axis = axisArrow.transpose().normalized().array().cast<ScalarVO>();
	rotation = Eigen::AngleAxis<ScalarVO>(theta, axis).toRotationMatrix();
}


// ������ת��Ϣ�õ���ת��������2�����õ���originArrow��ת��targetArrow����ת���󣬽���3D���Σ�
/*
	template <typename ScalarVO, typename ScalarV1, typename ScalarV2>
	void getRotationMat(
			Eigen::Matrix<ScalarVO, 3, 3>& rotation, \						 �������ת����(3��3)
			const Eigen::MatrixBase<DerivedV1>& originArrow,			 ��תǰ��������
			const Eigen::MatrixBase<DerivedV2>& targetArrow			 ��ת�����������
			)
		ע���������������ת����ȫ��Ĭ��Ϊ��������������
				��v,uΪ��������r,sΪ��������R * v == u; �ȼ��� r * R.transpose() == s;
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

	// 1. ��ԭ������Ŀ�귽��ƽ�У�
	RowVector3O originArrowCast{ static_cast<ScalarVO>(originArrow(0)), \
			static_cast<ScalarVO>(originArrow(1)), static_cast<ScalarVO>(originArrow(2))};
	RowVector3O targetArrowCast{ static_cast<ScalarVO>(targetArrow(0)), \
			static_cast<ScalarVO>(targetArrow(1)), static_cast<ScalarVO>(targetArrow(2)) };
	RowVector3O axisArrowCast = originArrowCast.cross(targetArrowCast);			// ��ת�᣻
	if (std::abs(axisArrowCast.norm()) < eps)
	{
		// 1.1 ��ԭ������Ŀ�귽����ͬ��
		if (originArrowCast.dot(targetArrowCast) > 0)
			return;
		else         // 1.2 ��ԭ������Ŀ�귽���෴��
		{
			// ����originArrowΪ�����ƽ���ڵ�һ��������
			axisArrowCast = RowVector3O{ 1, 1, 1 };
			ScalarVO x0 = originArrowCast(0);
			ScalarVO y0 = originArrowCast(1);
			ScalarVO z0 = originArrowCast(2);

			// ������o������Ϊx���������趨a������������������Ϊ1, x�������������o.dot(a)���x�������ٹ�һ���õ���ת������a
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

	// 2. ԭ������Ŀ�귽��ƽ�е�һ�����Σ�
	axisArrowCast.normalize();
	ScalarVO x0 = axisArrowCast(0);
	ScalarVO y0 = axisArrowCast(1);
	ScalarVO z0 = axisArrowCast(2);
	ScalarVO cosTheta = originArrowCast.dot(targetArrowCast) / (originArrowCast.norm() * targetArrowCast.norm());
	ScalarVO sinTheta = std::sqrt(1 - cosTheta * cosTheta);

	//		note: �ȼ���Eigen::AngleAxis<T>(theta, axis).toRotationMatrix()��������������������תtheta�Ƕȣ�
	rotation << 	cosTheta + (1 - cosTheta) * x0 * x0, (1 - cosTheta)* x0* y0 - sinTheta * z0, (1 - cosTheta)* x0* z0 + sinTheta * y0, \
		(1 - cosTheta)* y0* x0 + sinTheta * z0, cosTheta + (1 - cosTheta) * y0 * y0, (1 - cosTheta)* y0* z0 - sinTheta * x0, \
		(1 - cosTheta)* z0* x0 - sinTheta * y0, (1 - cosTheta)* z0* y0 + sinTheta * x0, cosTheta + (1 - cosTheta) * z0 * z0;
}


// �õ��������ϵ������������ʱ����תtheta���Ⱥ����ת����
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


// ����ֲ�����ϵ�µĵ��Ʊ任����������ϵ�ķ���任����
template <typename T, typename DerivedV>
bool getAffineL2G(Eigen::Matrix<T, 4, 4>& affine, \
	const Eigen::MatrixBase<DerivedV>& versLocal, \
	const Eigen::MatrixBase<DerivedV>& versGlobal, \
	const int sampleCount = 10)
{
	/*
	bool getAffineL2G(
		Eigen::Matrix<T, 4, 4>& affine,										������������ϵ�·���任����						
		const Eigen::MatrixBase<DerivedV>& versLocal,			�ֲ�����ϵ�еĵ���
		const Eigen::MatrixBase<DerivedV>& versGlobal,			��������ϵ�еĵ���
		const int sampleCount = 10											����ʹ�õ���������Ŀ
		)
		��������Ŀ������̫��Ҳ������̫�٣�
				̫�ٵĻ����ܻᵼ�����Է������A�������ȣ�
				̫��Ļ����ܵ���A����Ԫ�ص�ֵ�������
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

	// 1. ���⽨ģ��
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

	// 2. �����Է����飺(VT * V) * RT = (VT * U);
	Eigen::FullPivLU<MatrixXS> solverLU;			
	if (std::abs(VTV.determinant()) < eps)				// �����㲻�������任��Ψһ��
		return false;
	solverLU.compute(VTV);
	rt0 = solverLU.solve(VT * u0);
	rt1 = solverLU.solve(VT * u1);
	rt2 = solverLU.solve(VT * u2);

	// 3. �����
	RT.col(0) = rt0;
	RT.col(1) = rt1;
	RT.col(2) = rt2; 
	affine.block(0, 0, 3, 3) = RT.transpose();
	affine(0, 3) = center(0);
	affine(1, 3) = center(1);
	affine(2, 3) = center(2); 

	return true;
}


// �������������������������ȡԪ�أ�����1������������ΪEigen����
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


// �������������������������ȡԪ�أ�����2������������Ϊstd::vector<>
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


// ��������������ֻ����STL�������������������ȡԪ��
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


// ����flag������Դ��������ȡԪ�����������������1��
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


// ����flag������Դ��������ȡԪ�����������������2�����
template <typename DerivedO, typename DerivedI, typename TI>
bool subFromFlagVec(Eigen::PlainObjectBase<DerivedO>& matOut, \
	std::vector<int>& oldNewIdxInfo, std::vector<int>& newOldIdxInfo, \
	const Eigen::MatrixBase<DerivedI>& matIn, const std::vector<TI>& vec)
{
	using ScalarO = typename DerivedO::Scalar; 
	assert((vec.size() == matIn.rows()) && "Assert!!! input vec and matIn do not match.");
	matOut.resize(0, 0);

	// 1. ��ȡvec���Ϊ1�Ĳ�����
	const int N = vec.size();
	const int M = static_cast<int>(std::accumulate(vec.begin(), vec.end(), 0));
	int count = 0; 
	oldNewIdxInfo.clear();
	newOldIdxInfo.clear();
	matOut.resize(M, matIn.cols());
	for (unsigned i = 0; i < vec.size(); ++i)
		if (vec[i] > 0)
			matOut.row(count++) = matIn.row(i).array().cast<ScalarO>();

	// 2. ������������ӳ���
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
 

// ����flag������Դ��������ȡԪ�������������
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


// ��ά�����µ����ɷַ�����
template <typename DerivedV>
bool PCA3d(std::vector<Eigen::RowVector3d>& resultVecs,\
	std::vector<double>& resultVals, const Eigen::MatrixBase<DerivedV>& samples)
{
	assert(3 == samples.cols() && "Assert!!! Samples should be in 3-dimension space.");

	// 1. ���ж�����Ϊ��ά�������ݣ��������A�У�һ��һ��������
	const int  N = samples.rows();											// ���������� 
	Eigen::MatrixXd A = samples.array().cast<double>();
	Eigen::RowVector3d miu = A.colwise().mean();				// ������ֵ����miu������Ҳ�����ж�������ģ�

	// 2. ��ֵ��ƽ�Ƶ�ԭ�����������
	Eigen::MatrixXd B = A - Eigen::MatrixXd::Ones(N, 1) * miu;

	// 3. ��������Э�������
	Eigen::Matrix3d C = B.transpose() * B / (N - 1);

	// 4. ��Э������������ֵ����������
	Eigen::EigenSolver<Eigen::Matrix3d> es(C);
	Eigen::Matrix3d D = es.pseudoEigenvalueMatrix();			// �Խ���Ԫ��������ֵ
	Eigen::Matrix3d V = es.pseudoEigenvectors();					// ÿһ����������������������
 
	// 5. �����
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




///////////////////////////////////////////////////////////////////////////////////////////////////////// �������ɾ��ĺ�����

// ������β�����Ԫ��
template<typename T>
bool vecInsertNum(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const T num);

// ������β���������
template<typename T>
bool vecInsertVec(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec1, \
		const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec2);
 

// ����ĩβ�������/��������
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


// �����ұ߲������/��������
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


// ����flag�������õ�����������ӳ���ϵ�� 
/*
	void flagVec2oldNewIdxInfo(
			std::vector<int>& oldNewIdxInfo,
			std::vector<int>& newOldIdxInfo,
			const Eigen::MatrixBase<DerivedVI>& flagVec
			)

	flagVec�У���1 == flagVec(i)��ʾ����Ϊi�ĵ㱻ѡ�У���0 == flagVec(j)��ʾ����Ϊj�ĵ�û�б�ѡ�У�
	-1 == oldNewIdxInfo(i)��ʾԭ����������Ϊi�ĵ�û�б�ѡ�У����µ�����û�ж�Ӧ��
	index0 == newOldIdxInfo(i)��ʾ�µ���������Ϊi�ĵ㣬��ԭ�����е�����Ϊindex0;
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


// ����һ��flag������retVec����mat�ĵ�i�к�������vec��ȣ���retVec(i)==1���������0��
template <typename Derived1, typename Derived2>
Eigen::VectorXi rowInMat(const Eigen::MatrixBase<Derived1>& mat, \
	const Eigen::MatrixBase<Derived2>& rowVec)
{
	using T = typename Derived1::Scalar;
	const int rows = mat.rows();
	const int cols = mat.cols();
	Eigen::VectorXi retVec(rows);
	assert(rowVec.cols() == cols, "Error!!! Tow mats size do not match.");

	// ���бȽϣ�
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
			retVec = retVec.array() * tempMat.col(i).array();			// ������ˣ�

	return retVec;
}


// ����ϡ�����������MATLAB�е�sparse()
template <class ValueVector, typename T>
void sparse(Eigen::SparseMatrix<T>& SM, const Eigen::VectorXi& I, const Eigen::VectorXi& J,
	const ValueVector& values, const size_t m, const size_t n);


// ϡ�����ת�ã�Eigen::SparseMatrix�Դ���transpose()����̫������
template<typename T>
bool spMatTranspose(Eigen::SparseMatrix<T>& smOut, const Eigen::SparseMatrix<T>& smIn);


// ������ܾ���Ŀ����ڿ˻���Kronecker product��������1�������ܵĲ������̫���ˣ�����дģ���ػ��� 
template <typename Derived0, typename Derived1, typename Derived2>
bool kron(Eigen::PlainObjectBase<Derived0>& result, const Eigen::MatrixBase<Derived1>& m11, \
	const Eigen::MatrixBase<Derived2>& m22)
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
	using Scalar = typename Derived0::Scalar;							// ʹ�ô������ƣ���Ҫ����"typename"ǰ׺��
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
			static_cast<Scalar>(m1(rowIdx0, colIdx)) * m2.array().cast<Scalar>();		// ת������������е��������ͣ�
	};

	auto calcKronRow = [&](const unsigned rowIdx)
	{
		std::tuple<unsigned> pt0 = std::make_tuple(rowIdx);
		PARALLEL_FOR(0, n, calcKronCol, pt0);
	};

	PARALLEL_FOR(0, m, calcKronRow);

	return true;
}


// ������ܾ���Ŀ����ڿ˻���Kronecker product��������2
template <typename Derived1, typename Derived2>
Eigen::MatrixXd kron(const Eigen::MatrixBase<Derived1>& m11, \
	const Eigen::MatrixBase<Derived2>& m22)
{ 
	Eigen::MatrixXd result;
	kron(result, m11, m22);

	return result;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////// ���Է����飺

// ��ǡ���ĳ������Է�����Ax == b;
template<typename DerivedX, typename DerivedA, typename DerivedB>
bool solveLinearEquation(Eigen::PlainObjectBase<DerivedX>& x, const Eigen::MatrixBase<DerivedA>& A, \
	const Eigen::MatrixBase<DerivedB>& b)
{
	using ScalarX = typename DerivedX::Scalar;
	using MatrixXx = Eigen::Matrix<ScalarX, Eigen::Dynamic, Eigen::Dynamic>;
	using VectorXx = Eigen::Matrix<ScalarX, Eigen::Dynamic, 1>;

	// �����Է�����Ax == b;
	VectorXx b0 = b.array().cast<ScalarX>();
	Eigen::JacobiSVD<MatrixXx> solverSVD(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	x = solverSVD.solve(b0);

	return true;
}


// ��һϵ��ǡ���ĳ������Է�����AX == B;
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


// ͨ��SVD����ܾ������������
template<typename T>
double calcCondNum(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m);




///////////////////////////////////////////////////////////////////////////////////////////////////////// ��ֵ����ϣ�

// ������ά���ƣ������ϵ�ƽ�淽��ϵ��(a,b,c,d)������(a,b,c)��ƽ��ĵ�λ��������
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
	
	// A��������ֵ�ֽ⣺
	Eigen::JacobiSVD<MatrixXT> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	MatrixXT V = svd.matrixV();				// ���������
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


// ��ع����ʽ�������
template <typename T>
void ridgeRegressionPolyFitting(Eigen::VectorXd& theta, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	unsigned m = 4);


// ��С���˷���ϣ��ƽ�����׼��Բ���������xy��������룩
template<typename DerivedV>
Eigen::VectorXd fittingStandardEllipse(const Eigen::MatrixBase<DerivedV>& sampleVers)
{
	/*
		Eigen::VectorXd fittingStandardEllipse(								// ����������(a,c,d,e,f)��Ϊ��׼��Բ���̵�ϵ����
					const Eigen::MatrixXf& sampleVers					// ����������㣬��������XOYƽ���ϵĵ㣻
		)
	*/
	const double epsilon = 1e-8;							// ����������ֵС�ڴ�ֵʱ��ΪΪ0��
	using T = typename DerivedV::Scalar;
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
 




///////////////////////////////////////////////////////////////////////////////////////////////////////// ŷ�ϼ���

// ����ƽ����һ���ƽ�淨�򣬵õ�ƽ�淽��ϵ����
/*
	ƽ�淽��: ax +by+cz+d=0;
			����(a, b, c)��ƽ�������������ǵ�λ������

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


// ƽ�浽����ķ��ž��룺
template<typename T, typename DerivedV>
double disPlane2Vert(const std::vector<T>& coeffs, const Eigen::MatrixBase<DerivedV>& ver)
{
	return static_cast<double>(coeffs[0] * ver(0) + coeffs[1] * ver(1) + coeffs[2] * ver(2) + coeffs[3]);
}

 


///////////////////////////////////////////////////////////////////////////////////////////////////////// �˲���

// ����Ŀ��������˲�����������ע���˲�mask�ߴ����Ϊ����������
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
	unsigned im = (sizeMask + 1) / 2;					// ��Ĥ���ĵ��±ꣻ
	unsigned km = im;
	unsigned offset = (sizeMask - 1) / 2;

	// 1. �����������б�Ե���أ�������չ����
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

#include "tmp.tpp"