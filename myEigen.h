#pragma once
#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <fstream>
#include <set>
#include <unordered_set>
#include <set>
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


using namespace std;
using namespace Eigen;


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
	cout << endl;
}


// 反向遍历
template<typename T, typename F>
void revTraverseSTL(T& con, F f)
{
	std::for_each(con.rbegin(), con.rend(), f);
	cout << endl;
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


// lambda——打印cout支持的类型变量。
template <typename T>
auto disp = [](const T& arg)
{
	cout << arg << ", ";
};


template<typename T>
void dispQuat(const Quaternion<T>& q)
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
void dispSpMat(const spMat& sm, const unsigned startRow, const unsigned endRow)
{
	std::cout << "rows == " << sm.rows() << ", cols == " << sm.cols() << std::endl;
	for (unsigned i = startRow; i <= endRow; ++i)
		for (auto iter = spMat::InnerIterator(sm, i); iter; ++iter)
			std::cout << "(" << iter.row() << ", " << iter.col() << ") " << iter.value() << std::endl;
	std::cout << std::endl;
}


template<typename spMat>
void dispSpMat(const spMat& sm, const unsigned startRow, const unsigned endRow, const unsigned showElemsCount)
{
	unsigned count = 0;
	std::cout << "rows == " << sm.rows() << ", cols == " << sm.cols() << std::endl;
	for (unsigned i = startRow; i <= endRow; ++i)
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
void dispVec(const Matrix<T, N, 1>& vec)
{
	std::cout << ": rows == " << vec.rows() << std::endl;
	for (int i = 0; i < vec.rows(); ++i)
	{
		std::cout << i << "---\t" << vec(i) << std::endl;
	}
	std::cout << std::endl;
}


template<typename T, int N>
void dispVec(const Matrix<T, 1, N>& vec)
{
	std::cout << ": cols == " << vec.cols() << std::endl;
	for (int i = 0; i < vec.cols(); ++i)
	{
		std::cout << i << "---\t" << vec(i) << std::endl;
	}
	std::cout << std::endl;
}


template<typename T, int N>
void dispVecSeg(const Matrix<T, N, 1>& vec, const int start, const int end)
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
void dispVecSeg(const Matrix<T, 1, N>& vec, const int start, const int end)
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




/////////////////////////////////////// 图形生成接口：
bool interpolateToLine(MatrixXf& vers, const RowVector3f& start, const RowVector3f& end, const float SR, const bool SE);


// 生成中心在原点，边长为1，三角片数为12的正方体网格；
template	<typename DerivedV, typename DerivedI>
void genCubeMesh(Eigen::Matrix<DerivedV, Dynamic, Dynamic>& vers, Eigen::Matrix<DerivedI, Dynamic, Dynamic>& tris)
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
void genAABBmesh(const Eigen::AlignedBox<_Scalar, _AmbientDim>& aabb, Eigen::Matrix<_Scalar, Dynamic, Dynamic>& vers, \
	Eigen::MatrixXi& tris)
{
	Matrix<_Scalar, _AmbientDim, 1> minp = aabb.min();
	Matrix<_Scalar, _AmbientDim, 1> maxp = aabb.max();
	Matrix<_Scalar, _AmbientDim, 1> newOri = (minp + maxp) / 2.0;
	Matrix<_Scalar, _AmbientDim, 1> sizeVec = maxp - minp;

	genCubeMesh(vers, tris);
	vers.col(0) *= sizeVec(0);
	vers.col(1) *= sizeVec(1);
	vers.col(2) *= sizeVec(2);
	vers.rowwise() += newOri.transpose();
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
 

template <typename Derived>
bool subFromIdxVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const std::vector<int>& vec)
{
	Derived& matOut = matBaseOut.derived();
	matOut.resize(vec.size(), matBaseIn.cols());
	for (unsigned i = 0; i < vec.size(); ++i)
	{
		const int& index = vec[i];
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
	{
		if (vec(i) > 0)
			matOut.row(count++) = matBaseIn.row(i);
	}

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


// 解码边编码；
std::pair<int, int> decodeEdge(const std::int64_t code);


// 生成三角片编码——三个顶点索引排序后映射成64位无符号整型数
template <typename Index>
std::uint64_t encodeTriangle(const Index vaIdx, const Index vbIdx, const Index vcIdx)
{
	// 不区分正反面，即(a, b, c)和(a, c, b)映射为同一个编码；
	unsigned long triIdxLimit = 0x1FFFFF;								// 最多为21位全1, 0x1FFFFF ==  2097181, 两百多万三角片；
	if (vaIdx > triIdxLimit || vbIdx > triIdxLimit || vcIdx > triIdxLimit)
		return 0;			// 索引超出范围
	if (vaIdx == vbIdx || vaIdx == vcIdx || vbIdx == vcIdx || vaIdx < 0 || vbIdx < 0 || vcIdx < 0)
		return 0;			// 非法三角片

	std::vector<std::uint64_t> vec{static_cast<std::uint64_t>(vaIdx), static_cast<std::uint64_t>(vbIdx), static_cast<std::uint64_t>(vcIdx)};
	std::sort(vec.begin(), vec.end());
	std::uint64_t code = 0;
	code |= (vec[0] << 42);
	code |= (vec[1] << 21);
	code |= vec[2];
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
template<typename T, int N>
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
	{
		for (int i = 1; i < cols; ++i)
			retVec = retVec.array() * tempMat.col(i).array();			// 逐列相乘：
	}

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
bool matWriteToFile(const char* fileName, const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat)
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
bool matReadFromFile(Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName)
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
	{
		return false;
	}

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
	{
		return false;
	}

	// 读矩阵元素
	mat.resize(row, col);
	T* ptr = mat.data();
	while (!file.eof())
	{
		file.getline(cstr, LINE_LENGTH);
		str1 = cstr;
		if (str1.size() == 0)
		{
			break;
		}
		std::string::iterator iter = str1.begin();
		for (unsigned j = 0; j < 3; ++j)
		{
			if (iter == str1.end())
			{
				break;
			}

			if (*iter == ' ')						// 负号后面可能有空格，需要去除
			{
				iter = str1.erase(iter);
			}
			iter++;
		}
		*ptr++ = std::stoi(str1);
	}
	file.close();

};


template	<typename Scalar, typename Index>
void objReadMeshMat(Eigen::Matrix<Scalar, Dynamic, Dynamic>& vers, Eigen::Matrix<Index, Dynamic, Dynamic>& tris, const char* fileName)
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
			// Vector3f ver;
			ver(0) = (float)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(1) = (float)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(2) = (float)atof(tmpBuffer);

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
void objReadVerticesMat(Eigen::Matrix<Scalar, Dynamic, Dynamic>& vers, const char* fileName)
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
			// Vector3f ver(Vector3f::Zero());
			ver(0) = (float)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(1) = (float)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(2) = (float)atof(tmpBuffer);

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


void printDirEigen(const char* pathName, const RowVector3f& origin, const RowVector3f& dir);


void printCoordinateEigen(const char* pathName, const RowVector3f& origin, const RowVector3f& xdir, \
	const RowVector3f& ydir, const RowVector3f& zdir);


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
void objWritePath(const char* pathName, const std::vector<DerivedI>& path, const Eigen::Matrix<DerivedV, Dynamic, Dynamic>& vers)
{
	if (path.size() <= 1)
		return;

	unsigned edgesCount = path.size() - 1;
	Eigen::Matrix<DerivedI, Dynamic, Dynamic> pathEdges(edgesCount, 2);

	for (unsigned i = 0; i < edgesCount; ++i)
		pathEdges(i, 0) = path[i];
	for (unsigned i = 0; i < edgesCount; ++i)
		pathEdges(i, 1) = path[i + 1];
	objWriteEdgesMat(pathName, pathEdges, vers);
}
 

// objWirteTreePath()——输入树向量或路径向量，写入到OBJ文件中：
template <typename DerivedV>
void objWriteTreePath(const char* pathName, const Eigen::VectorXi& treeVec, const Eigen::Matrix<DerivedV, Dynamic, Dynamic>& vers)
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
void objReadVerticesHomoMat(MatrixXf& vers, const char* fileName);

void objWriteVerticesHomoMat(const char* fileName, const MatrixXf& vers);

void vers2homoVers(MatrixXf& homoVers, const MatrixXf& vers);

MatrixXf vers2homoVers(const MatrixXf& vers);

void homoVers2vers(MatrixXf& vers, const MatrixXf& homoVers);

MatrixXf homoVers2vers(const MatrixXf& homoVers);


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
bool solveLinearEquation(Matrix<T, N, 1>& x, const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, const Matrix<T, N, 1>& b)
{
	// 解线性方程组Ax == b;
	JacobiSVD<Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(A, ComputeThinU | ComputeThinV);
	x = svd.solve(b);

	return true;
}


// 解一系列恰定的稠密线性方程组AX == B;
template <typename T>
bool solveLinearEquations(Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& X, const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& B)
{
	if (A.rows() != B.rows())
	{
		return false;
	}

	JacobiSVD<Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(A, ComputeThinU | ComputeThinV);
	X.resize(A.cols(), B.cols());
	for (int i = 0; i < B.cols(); ++i)
	{
		Matrix < T, Eigen::Dynamic, 1> x = svd.solve(B.col(i));
		X.col(i) = x;
	}

	return true;
}


// 霍纳方法（秦九昭算法）求多项式的值
template<typename T, int N>
float hornersPoly(const Eigen::Matrix<T, N, 1>& coeffs, const float x)
{
	// coeffs是{a0, a1, a2, ..., an}组成的(n+1)列向量，多项式为p == a0 + a1*x + a2*x^2 + ... + an* x^n; 
	int n = coeffs.rows() - 1;
	if (n < 0)
	{
		return NAN;
	}

	float result = coeffs(n);
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
	{
		return NAN;
	}

	if (coeffsCount == 1)
	{
		return 1.0f;
	}

	VectorXf diffCoeffs(coeffsCount - 1);
	for (int i = 0; i < diffCoeffs.rows(); ++i)
	{
		diffCoeffs(i) = (i + 1) * coeffs(i + 1);
	}

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
void ridgeRegressionPolyFitting(VectorXf& theta, const MatrixXf& vers);


Matrix3f getRotationMat(const RowVector3f& originArrow, const RowVector3f& targetArrow);



// 最小二乘法拟合（逼近）标准椭圆（长短轴和xy坐标轴对齐）
VectorXd fittingStandardEllipse(const MatrixXf& sampleVers);


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
		RowVector3d va, vb, vc, arrow1, arrow2, norm;
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
		RowVector3d norm = triNorms.row(i);

		// 若存在退化三角片，则平面系数都写为NAN:
		if (std::isinf(norm(0)))
			planeCoeff(i, 3) = INFINITY;

		RowVector3d va = vers.row(tris(i, 0)).array().cast<double>();

		// va点在平面上 → va点到平面的距离为0 → norm.dot(va) +d == 0; → d = -norm.dot(va);
		planeCoeff(i, 3) = -norm.dot(va);
	}

	return true;
}


// 计算三角网格的体积：
template<typename DerivedV, typename DerivedI>
double meshVolume(const Eigen::MatrixBase<DerivedV>& V, const Eigen::MatrixBase<DerivedI>& T)
{
	double volume = -1.0;
	const DerivedV& vers = V.derived();
	const DerivedI& tris = T.derived();
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
template<typename DerivedI>
bool adjMatrix(const Eigen::PlainObjectBase<DerivedI>& tris, Eigen::SparseMatrix<int>& adjSM_eCount, \
		Eigen::SparseMatrix<int>& adjSM_eIdx)
{
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

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
	adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());										// 权重为该有向边重复的次数；
	adjSM_eIdx.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());		// 权重为该有向边的索引；

	return true;
}


// 边缘边：
template<typename DerivedI>
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdryEdges, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	adjMatrix(tris, adjSM_eCount, adjSM_eIdx);
	Eigen::SparseMatrix<int> adjSM_tmp, adjSM_ueCount;
	spMatTranspose(adjSM_tmp, adjSM_eCount);
	adjSM_ueCount = adjSM_eCount + adjSM_tmp;
	std::vector<std::pair<int, int>> edgesVec;
	traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
		{
			if (1 == iter.value() && iter.row() > iter.col())
				edgesVec.push_back({ iter.row(), iter.col() });
		});
	edges2mat(bdryEdges, edgesVec);

	return true;
}


// 求三角网格中的非流形半边：
template<typename DerivedI>
bool nonManifoldEdges(const Eigen::PlainObjectBase<DerivedI>& tris, Eigen::MatrixXi& nmnEdges)
{
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	adjMatrix(tris, adjSM_eCount, adjSM_eIdx);
	
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
 

// buildAdjacency()——计算流形网格的三角片邻接信息：
using ttTuple = std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>;

bool buildAdjacency(const Eigen::MatrixXi& tris, Eigen::MatrixXi& ttAdj_nmEdge, \
	std::vector<ttTuple>& ttAdj_nmnEdge, std::vector<ttTuple>& ttAdj_nmnOppEdge);


// triangleGrow()——区域生长：
template <typename DerivedV, typename DerivedI>
bool triangleGrow(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
		const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx)
{
	/*
		
			bool VSRootmerging::triangleGrow(
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

	// 输入参数中的triIdx为包含最高点的三角片索引；
	std::unordered_set<int> finalTriIdx;						// 联通三角片集合；
	std::unordered_set<int> workingSet;						// 用于寻找相邻三角片的队列；

	Eigen::MatrixXi edges;
	Eigen::MatrixXi adjTrisIdxes = -Eigen::MatrixXi::Ones(trisCount, 3);			// 三角片三条边相邻的三角片索引，不存在则写-1；
	Eigen::MatrixXi ttAdj_nmEdge;		// 三角片的非边缘流形有向边邻接的三角片索引；	trisCount * 3(i, 0)为索引为i的三角片的非边缘流形有向边(vb, vc)邻接的三角片的索引；

	Eigen::SparseMatrix<int> adjSM;								// 邻接矩阵；
	Eigen::SparseMatrix<int> adjSM_eCount;					// 邻接矩阵，索引为该有向边重复的次数；
	Eigen::SparseMatrix<int> adjSM_weighted;				// 邻接矩阵，元素若对应流形边则为该边索引，若对应非流形边则为该边所有索引的和；
	Eigen::SparseMatrix<int> adjSM_weighted_opp;		// adjSM_weighted的转置，表示对边的信息；
	Eigen::SparseMatrix<int> adjSM_ueCount;				// 权重为无向边ij关联的三角片数量；
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;				// 非边缘流形无向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB;					// 非边缘流形有向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;			// 非边缘流形有向边对边的邻接矩阵

	Eigen::VectorXi etAdj_mnEdge;							// 边所在三角片的索引，非流形边或边缘边写-1；

	std::vector<int> etInfo;												// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
	std::vector<int> edgesIdx_MN_NB;							// 非边缘流形有向边的索引；
	std::vector<int> edgesIdx_MN_NB_opp;					// 非边缘流形有向边的对边的索引；
	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_weighted的triplet数据；	

	edges.resize(edgesCount, 2);
	Eigen::MatrixXi vaIdxes = tris.col(0);
	Eigen::MatrixXi vbIdxes = tris.col(1);
	Eigen::MatrixXi vcIdxes = tris.col(2);
	edges.block(0, 0, trisCount, 1) = vbIdxes;
	edges.block(trisCount, 0, trisCount, 1) = vcIdxes;
	edges.block(trisCount * 2, 0, trisCount, 1) = vaIdxes;
	edges.block(0, 1, trisCount, 1) = vcIdxes;
	edges.block(trisCount, 1, trisCount, 1) = vaIdxes;
	edges.block(trisCount * 2, 1, trisCount, 1) = vbIdxes;

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

		// 1.1 生成有向边数据
		edges = Eigen::MatrixXi::Zero(edgesCount, 2);
		Eigen::MatrixXi vaIdxes = tris.col(0);
		Eigen::MatrixXi vbIdxes = tris.col(1);
		Eigen::MatrixXi vcIdxes = tris.col(2);
		edges.block(0, 0, trisCount, 1) = vbIdxes;
		edges.block(trisCount, 0, trisCount, 1) = vcIdxes;
		edges.block(trisCount * 2, 0, trisCount, 1) = vaIdxes;
		edges.block(0, 1, trisCount, 1) = vcIdxes;
		edges.block(trisCount, 1, trisCount, 1) = vaIdxes;
		edges.block(trisCount * 2, 1, trisCount, 1) = vbIdxes;

		// 1.2 生成三种有向边邻接矩阵：adjSM, adjSM_eCount, smElems_weighted
		smElems.reserve(edgesCount);
		smElems_weighted.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
		{
			smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});
			smElems_weighted.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), i});
		}

		adjSM_eCount.resize(versCount, versCount);
		adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());	// 权重为该有向边重复的次数；
		adjSM = adjSM_eCount;																		// 有向边邻接矩阵；
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		若adjSM_weighted(i, j)对应的是流形有向边，该权重值为该边的索引；若是非流形有向边，则该边对应两个边索引，权重为两个索引之和；
		adjSM_weighted.resize(versCount, versCount);
		adjSM_weighted.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());
		adjSM_weighted_opp = adjSM_weighted.transpose();

		//	1.3 生成无向边邻接矩阵：
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// 无向边邻接矩阵；
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.4 生成非边缘流形无向边邻接矩阵：
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

		//	1.5 生成非边缘流形有向边、及其对边的邻接矩阵
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
				edgesIdx_MN_NB.push_back(adjSM_weighted.coeffRef(iter.row(), iter.col()));
			});

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB_opp.push_back(adjSM_weighted_opp.coeffRef(iter.row(), iter.col()));
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
	for (const auto& index : finalTriIdx)
	{
		finalVerIdx.insert(tris(static_cast<int>(index), 0));
		finalVerIdx.insert(tris(static_cast<int>(index), 1));
		finalVerIdx.insert(tris(static_cast<int>(index), 2));
	}

	VectorXi finalVerIdxVec(finalVerIdx.size());								// 输出顶点的索引，是原网格中的索引；
	VectorXi oldNewIdxInfo = -VectorXi::Ones(vers.rows());		// 新老索引表——下标为老索引，元素值为新索引，不在新网格中的顶点写为-1
	auto iter = finalVerIdx.begin();
	for (int i = 0; i < finalVerIdx.size(); ++i)
		finalVerIdxVec(i) = static_cast<int>(*iter++);

	int setOrder = 0;
	for (const auto& index : finalVerIdx)
		oldNewIdxInfo[static_cast<int>(index)] = setOrder++;
	subFromIdxVec(versOut, vers, finalVerIdxVec);

	// 6. 原网格所有三角片中的顶点索引由老的换为新的。
	MatrixXi tempMat = tris;
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
		RowVector3i currentTri = tempMat.row(i);
		if ((currentTri.array() >= 0).all())
			trisOut.row(count++) = currentTri;
	}
	trisOut.conservativeResize(count, 3);

	return true;
}
 

// triangleGrowSplitMesh()——区域生长将输入网格分为多个单连通网格：
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

	Eigen::MatrixXi edges;
	Eigen::MatrixXi adjTrisIdxes = -Eigen::MatrixXi::Ones(trisCount, 3);			// 三角片三条边相邻的三角片索引，不存在则写-1；
	Eigen::MatrixXi ttAdj_nmEdge;		// 三角片的非边缘流形有向边邻接的三角片索引；	trisCount * 3(i, 0)为索引为i的三角片的非边缘流形有向边(vb, vc)邻接的三角片的索引；

	Eigen::SparseMatrix<int> adjSM;								// 邻接矩阵；
	Eigen::SparseMatrix<int> adjSM_eCount;					// 邻接矩阵，索引为该有向边重复的次数；
	Eigen::SparseMatrix<int> adjSM_weighted;				// 邻接矩阵，元素若对应流形边则为该边索引，若对应非流形边则为该边所有索引的和；
	Eigen::SparseMatrix<int> adjSM_weighted_opp;		// adjSM_weighted的转置，表示对边的信息；
	Eigen::SparseMatrix<int> adjSM_ueCount;				// 权重为无向边ij关联的三角片数量；
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;				// 非边缘流形无向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB;					// 非边缘流形有向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;			// 非边缘流形有向边对边的邻接矩阵

	Eigen::VectorXi etAdj_mnEdge;							// 边所在三角片的索引，非流形边或边缘边写-1；

	std::vector<int> etInfo;												// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
	std::vector<int> edgesIdx_MN_NB;							// 非边缘流形有向边的索引；
	std::vector<int> edgesIdx_MN_NB_opp;					// 非边缘流形有向边的对边的索引；
	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_weighted的triplet数据；	

	edges.resize(edgesCount, 2);
	Eigen::MatrixXi vaIdxes = tris.col(0);
	Eigen::MatrixXi vbIdxes = tris.col(1);
	Eigen::MatrixXi vcIdxes = tris.col(2);
	edges.block(0, 0, trisCount, 1) = vbIdxes;
	edges.block(trisCount, 0, trisCount, 1) = vcIdxes;
	edges.block(trisCount * 2, 0, trisCount, 1) = vaIdxes;
	edges.block(0, 1, trisCount, 1) = vcIdxes;
	edges.block(trisCount, 1, trisCount, 1) = vaIdxes;
	edges.block(trisCount * 2, 1, trisCount, 1) = vbIdxes;

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

		// 1.1 生成有向边数据
		edges = Eigen::MatrixXi::Zero(edgesCount, 2);
		Eigen::MatrixXi vaIdxes = tris.col(0);
		Eigen::MatrixXi vbIdxes = tris.col(1);
		Eigen::MatrixXi vcIdxes = tris.col(2);
		edges.block(0, 0, trisCount, 1) = vbIdxes;
		edges.block(trisCount, 0, trisCount, 1) = vcIdxes;
		edges.block(trisCount * 2, 0, trisCount, 1) = vaIdxes;
		edges.block(0, 1, trisCount, 1) = vcIdxes;
		edges.block(trisCount, 1, trisCount, 1) = vaIdxes;
		edges.block(trisCount * 2, 1, trisCount, 1) = vbIdxes;

		// 1.2 生成三种有向边邻接矩阵：adjSM, adjSM_eCount, smElems_weighted
		smElems.reserve(edgesCount);
		smElems_weighted.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
		{
			smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});
			smElems_weighted.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), i});
		}

		adjSM_eCount.resize(versCount, versCount);
		adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());		// 权重为该有向边重复的次数；

		
		//	若输入网格中存在非流形边（多个三角片包含同一条有向边），则return false;
		bool MNcheckFlag = true;
		traverseSparseMatrix(adjSM_eCount, [&adjSM_eCount, &MNcheckFlag](auto& iter)
			{
				if (iter.value() > 1)
					MNcheckFlag = false;
			});
		if (!MNcheckFlag)
			return false;

		adjSM = adjSM_eCount;																		// 有向边邻接矩阵；

		//		若adjSM_weighted(i, j)对应的是流形有向边，该权重值为该边的索引；若是非流形有向边，则该边对应两个边索引，权重为两个索引之和；
		adjSM_weighted.resize(versCount, versCount);
		adjSM_weighted.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());
		adjSM_weighted_opp = adjSM_weighted.transpose();

		//	1.3 生成无向边邻接矩阵：
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// 无向边邻接矩阵；
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.4 生成非边缘流形无向边邻接矩阵：
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

		//	1.5 生成非边缘流形有向边、及其对边的邻接矩阵
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
				edgesIdx_MN_NB.push_back(adjSM_weighted.coeffRef(iter.row(), iter.col()));
			});

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB_opp.push_back(adjSM_weighted_opp.coeffRef(iter.row(), iter.col()));
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


