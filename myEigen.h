#pragma once
#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <fstream>
#include <set>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
#include <initializer_list>
#include <memory>
#include <windows.h>

#include "Eigen/Dense"
#include "Eigen/Sparse"


using namespace std;
using namespace Eigen;



/////////////////////////////////////////////////////////////////////////////////////////////////// 控制台打印接口

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


template<typename T>
void dispMat(const Matrix<T, Dynamic, Dynamic>& mat)
{
	std::cout << ": rows == " << mat.rows() << ",  cols == " << mat.cols() << std::endl;
	for (int i = 0; i < mat.rows(); ++i)
	{
		std::cout << i << " ---\t";
		for (int j = 0; j < mat.cols(); ++j)
		{
			std::cout << mat(i, j) << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


template<typename T>
void dispMatBlock(const Matrix<T, Dynamic, Dynamic>& mat, const int rowStart, const int rowEnd, const int colStart, const int colEnd)
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


#if 0
template<typename T>
void dispSpMat(const SparseMatrix<T>& mat, const int showElems = -1)
{
	std::cout << ": rows == " << mat.rows() << ",  cols == " << mat.cols() << std::endl;
	std::cout << "非零元素数：" << mat.nonZeros() << std::endl;
	int count = 0;
	for (int k = 0; k < mat.outerSize(); ++k)
	{
		//for(decltype(mat)::InnerIterator it(mat, k); ; ++it)
		for (SparseMatrix<T>::InnerIterator it(mat, k); it; ++it)// 64位时，貌似T类型在编译器不确定的话，此句会报错
		{
			std::cout << "(" << it.row() << ", " << it.col() << ") ---\t" << it.value() << std::endl;
			count++;
			if (count >= showElems && showElems >= 0)
				return;
		}
	}
}
#endif

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



///////////////////////////////////////////////////////////////////////////////////////////////////// 不同数据类型的变换

// 根据索引向量从源矩阵中提取元素生成输出矩阵。
template <typename T>
bool subFromIdxVec(Matrix<T, Dynamic, Dynamic>& matOut, const Matrix<T, Dynamic, Dynamic>& matIn, const VectorXi& vec)
{
	matOut.resize(vec.rows(), matIn.cols());
	for (unsigned i = 0; i < vec.rows(); ++i)
	{
		const int& index = vec(i);
		matOut.row(i) = matIn.row(index);
	}

	return true;
}


// 根据flag向量从源矩阵中提取元素生成输出矩阵。
template<typename T>
bool subFromFlagVec(Matrix<T, Dynamic, Dynamic>& matOut, const Matrix<T, Dynamic, Dynamic>& matIn, const VectorXi& vec)
{
	matOut.resize(vec.sum(), matIn.cols());
	int count = 0;
	for (unsigned i = 0; i < vec.rows(); ++i)
	{
		if (vec(i) > 0)
		{
			matOut.row(count++) = matIn.row(i);
		}
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
std::vector<T>  vec2Vec(const Eigen::Matrix<T, Dynamic, 1>& vIn)		// 注：Matrix<T, Dynamic, 1> == VectorXT;
{
	unsigned elemCount = vIn.rows();
	std::vector<T> vOut(elemCount, 1);
	std::memcpy(&vOut[0], vIn.data(), sizeof(T) * elemCount);

	return vOut;
}

template<typename T>
Eigen::Matrix<T, Dynamic, 1> vec2Vec(const std::vector<T>& vIn)
{
	unsigned elemCount = vIn.size();
	Eigen::Matrix<T, Dynamic, 1> vOut;
	vOut.resize(elemCount, 1);
	std::memcpy(vOut.data(), &vIn[0], sizeof(T) * elemCount);

	return vOut;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////// 矩阵的增删查改

// 向量插入数据
template<typename T>
bool vecInsertNum(Eigen::Matrix<T, Dynamic, 1>& vec, const T num)
{
	vec.conservativeResize(vec.rows() + 1, 1);
	vec(vec.rows() - 1) = num;
	return true;
}


// 向量插入向量
template<typename T>
bool vecInsertVec(Eigen::Matrix<T, Dynamic, 1>& vec1, const Eigen::Matrix<T, Dynamic, 1>& vec2)
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
bool matInsertRows(Matrix<T, Dynamic, Dynamic>& mat, const Matrix<T, Dynamic, Dynamic>& mat1)
{
	unsigned cols = mat1.cols();
	unsigned currentRows = mat.rows();
	unsigned addRows = mat1.rows();
	mat.conservativeResize(currentRows + addRows, cols);
	for (unsigned i = 0; i < addRows; ++i)
	{
		mat.row(currentRows + i) = mat1.row(i);
	}


	return true;
}


//	矩阵末端插入行向量
template<typename T, int N>
bool matInsertRows(Matrix<T, Dynamic, Dynamic>& mat, const Matrix<T, 1, N>& rowVec)
{
	unsigned cols = rowVec.cols();
	unsigned currentRows = mat.rows();

	if (0 == cols)
	{
		return false;
	}

	if (rowVec.cols() != mat.cols() && rowVec.cols() != 0)
	{
		return false;
	}

	mat.conservativeResize(currentRows + 1, cols);
	mat.row(currentRows) = rowVec;

	return true;
}


// 返回一个类似于matlab索引向量的int列向量retVec，若mat的第i行和行向量vec相等，则retVec(i)==1，否则等于0；若程序出错则retVec所有元素为-1
template<typename T, int N>
VectorXi vecInMat(const Matrix<T, Dynamic, Dynamic>& mat, const Matrix<T, 1, N>& vec)
{
	int rows = mat.rows();
	int cols = mat.cols();

	VectorXi retVec(rows);
	if (vec.cols() != cols)
	{
		retVec = -VectorXi::Ones(rows);
		return retVec;
	}

	// 逐列比较：
	MatrixXi tempMat(rows, cols);
	for (int i = 0; i < cols; ++i)
	{
		tempMat.col(i) = (mat.col(i).array() == vec(i)).select(VectorXi::Ones(rows), VectorXi::Zero(rows));
	}

	retVec = tempMat.col(0);

	if (cols > 1)
	{
		for (int i = 1; i < cols; ++i)
		{
			retVec = retVec.array() * tempMat.col(i).array();			// 逐列相乘：
		}
	}

	return retVec;
}

template<typename T>
VectorXi vecInMat(const Matrix<T, Dynamic, Dynamic>& mat, const Matrix<T, 1, Dynamic>& vec)
{
	int rows = mat.rows();
	int cols = mat.cols();

	VectorXi retVec(rows);
	if (vec.cols() != cols)
	{
		retVec = -VectorXi::Ones(rows);
		return retVec;
	}

	// 逐列比较：
	MatrixXi tempMat(rows, cols);
	for (int i = 0; i < cols; ++i)
		tempMat.col(i) = (mat.col(i).array() == vec(i)).select(VectorXi::Ones(rows), VectorXi::Zero(rows));

	retVec = tempMat.col(0);

	if (cols > 1)
	{
		for (int i = 1; i < cols; ++i)
			retVec = retVec.array() * tempMat.col(i).array();			// 逐列相乘：
	}

	return retVec;
}


// 网格串联——合并两个孤立的网格到一个网格里
void concatMeshMat(MatrixXf& vers, MatrixXi& tris, const MatrixXf& vers1, const MatrixXi& tris1);



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
bool matWriteToFile(const char* fileName, const Matrix<T, Dynamic, Dynamic>& mat)
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
bool matReadFromFile(Matrix<T, Dynamic, Dynamic>& mat, const char* fileName)
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


template	<typename DerivedV, typename DerivedI>
void objReadMeshMat(Eigen::Matrix<DerivedV, Dynamic, Dynamic>& vers, Eigen::Matrix<DerivedI, Dynamic, Dynamic>& tris, const char* fileName)
{
	char* pTmp = NULL;
	std::ifstream ifs(fileName);//cube bunny Eight
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

			Vector3f ver;
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

			Vector3i tri;
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


template	<typename DerivedV, typename DerivedI>
void objWriteMeshMat(const char* fileName, const Eigen::Matrix<DerivedV, Dynamic, Dynamic>& vers, \
	const Eigen::Matrix<DerivedI, Dynamic, Dynamic>& tris)
{
	std::ofstream dstFile(fileName);

	for (int j = 0; j < vers.rows(); j++)
	{
		char szBuf[256] = { 0 };
		sprintf_s(szBuf, 256, "v %f %f %f", vers(j, 0), vers(j, 1), vers(j, 2));
		dstFile << szBuf << "\n";
	}

	for (unsigned j = 0; j < tris.rows(); ++j)
	{
		char szBuf[256] = { 0 };
		sprintf_s(szBuf, 256, "f %d %d %d", tris(j, 0) + 1, tris(j, 1) + 1, tris(j, 2) + 1);
		dstFile << szBuf << "\n";
	}
};


template	<typename DerivedV>
void objReadVerticesMat(Eigen::Matrix<DerivedV, Dynamic, Dynamic>& vers, const char* fileName)
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

			Vector3f ver(Vector3f::Zero());
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
void objWriteVerticesMat(const char* fileName, const Eigen::Matrix<DerivedV, Dynamic, Dynamic>& vers)
{
	std::ofstream dstFile(fileName);
	for (int i = 0; i < vers.rows(); i++)
		dstFile << "v " << vers(i, 0) << " " << vers(i, 1) << " " << vers(i, 2) << std::endl;
	dstFile.close();
};


void printDirEigen(const char* pathName, const RowVector3f& origin, const RowVector3f& dir);


void printCoordinateEigen(const char* pathName, const RowVector3f& origin, const RowVector3f& xdir, \
	const RowVector3f& ydir, const RowVector3f& zdir);


// 边数据写入到OBJ文件中：
template	<typename DerivedV, typename DerivedI>
void objWriteEdgesMat(const char* pathName, const Eigen::Matrix<DerivedI, Dynamic, Dynamic>& edges, \
	const Eigen::Matrix<DerivedV, Dynamic, Dynamic>& vers)
{
	std::ofstream dstFile(pathName);
	for (int i = 0; i < vers.rows(); ++i)
		dstFile << "v " << vers(i, 0) << " " << vers(i, 1) << " " << vers(i, 2) << std::endl;

	for (int i = 0; i < edges.rows(); ++i)
		dstFile << "l " << edges(i, 0) + 1 << " " << edges(i, 1) + 1 << std::endl;

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

// 解恰定的稠密线性方程组Ax == b;
template<typename T, int N>
bool solveLinearEquation(Matrix<T, N, 1>& x, const Matrix<T, Dynamic, Dynamic>& A, const Matrix<T, N, 1>& b)
{
	// 解线性方程组Ax == b;
	JacobiSVD<Matrix<T, Dynamic, Dynamic>> svd(A, ComputeThinU | ComputeThinV);
	x = svd.solve(b);

	return true;
}


// 解一系列恰定的稠密线性方程组AX == B;
template <typename T>
bool solveLinearEquations(Matrix<T, Dynamic, Dynamic>& X, const Matrix<T, Dynamic, Dynamic>& A, const Matrix<T, Dynamic, Dynamic>& B)
{
	if (A.rows() != B.rows())
	{
		return false;
	}

	JacobiSVD<Matrix<T, Dynamic, Dynamic>> svd(A, ComputeThinU | ComputeThinV);
	X.resize(A.cols(), B.cols());
	for (int i = 0; i < B.cols(); ++i)
	{
		Matrix < T, Dynamic, 1> x = svd.solve(B.col(i));
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


// 多项式插值
void polyInterpolation();


// 高斯插值
void gaussInterpolation();


// 最小二乘多项式拟合曲线：
void leastSquarePolyFitting();


// 岭回归多项式拟合曲线
void ridgeRegressionPolyFitting(VectorXf& theta, const MatrixXf& vers);


Matrix3f getRotationMat(const RowVector3f& originArrow, const RowVector3f& targetArrow);


/////////////////////////////////////// 图形生成接口：
bool interpolateToLine(MatrixXf& vers, const RowVector3f& start, const RowVector3f& end, const float SR, const bool SE);



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





