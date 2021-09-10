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


#define _DEBUG
#define PI 3.14159
#define EPS	0.000001

using namespace std;
using namespace Eigen;

//声明***********************************************************************************************************************************************
 
auto readNextData = [](char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize)->unsigned
{
	unsigned nIndx = 0;

	while ((pszBuf[0] == ' ') ||
		(pszBuf[0] == '\n') ||
		(pszBuf[0] == '\t') ||
		(pszBuf[0] == '\r'))
	{
		pszBuf++;
		nCount++;
	}

	while ((pszBuf[0] != ' ') &&
		(pszBuf[0] != '\n') &&
		(pszBuf[0] != '\t') &&
		(pszBuf[0] != '\r') &&
		(pszBuf[0] != '\null') &&
		(pszBuf[0] != 0) &&
		(nIndx < nMaxSize))
	{
		validData[nIndx++] = pszBuf[0];
		pszBuf++;
		nCount++;
	}

	validData[nIndx] = 0;
	return nIndx;
};



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
	std:: cout << q.vec() << std::endl;
}

  
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

 

// 文件中读写矩阵
auto matIwriteToFile = [](const char* fileName, const MatrixXi& mat)
{
	std::ofstream file(fileName);
	file << "row " << mat.rows() << std::endl;
	file << "col " << mat.cols() << std::endl;
	const int* intPtr = mat.data();
	unsigned elemCount = mat.rows() * mat.cols();
	for (unsigned i = 0; i < elemCount; ++i)
	{
		file << *intPtr++ << std::endl;
	}
	file.close();

	return true;
};

auto matFwriteToFile = [](const char* fileName, const MatrixXf& mat)
{
	std::ofstream file(fileName);
	file << "row " << mat.rows() << std::endl;
	file << "col " << mat.cols() << std::endl;
	const float* fPtr = mat.data();
	unsigned elemCount = mat.rows() * mat.cols();
	for (unsigned i = 0; i < elemCount; ++i)
	{
		file << *fPtr++ << std::endl;
	}
	file.close();

	return true;
};

auto matIreadFromFile = [](MatrixXi& mat, const char* fileName)
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
	int* intPtr = mat.data();
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
		*intPtr++ = std::stoi(str1);
	}
	file.close();

};

auto matFreadFromFile = [](MatrixXf& mat, const char* fileName)
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
	float* fPtr = mat.data();
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
		*fPtr++ = std::stof(str1);
	}
	file.close();

};

auto vecIwriteToFile = [](const char* fileName, const VectorXi& vec)
{
	std::ofstream file(fileName);
	file << "row " << vec.rows() << std::endl;
	const int* intPtr = vec.data();
	unsigned elemCount = vec.rows();
	for (unsigned i = 0; i < elemCount; ++i)
	{
		file << *intPtr++ << std::endl;
	}
	file.close();

	return true;
};

auto vecIreadFromFile = [](VectorXi& vec, const char* fileName)
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

	// 读矩阵元素
	vec.resize(row);
	int* intPtr = vec.data();
	while (!file.eof())
	{
		file.getline(cstr, LINE_LENGTH);
		str1 = cstr;
		str2.clear();
		str2.insert(str2.end(), str1.begin(), str1.begin() + 3);
		if (std::strcmp(str2.c_str(), "col") == 0)
		{
			continue;
		}

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
		*intPtr++ = std::stoi(str1);
	}
	file.close();
};


// eigen的向量和std::vector<T>相互转换
template<typename T, size_t N>
std::vector<T>  vec2Vec(const Eigen::Matrix<T, N, 1>& vIn)
{
	unsigned elemCount = vIn.rows();
	std::vector<T> vOut(elemCount, 1);
	std::memcpy(&vOut[0], vIn.data(), sizeof(T) * elemCount);

	return vOut;
}

template<typename T, size_t N>
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
	{
		tempMat.col(i) = (mat.col(i).array() == vec(i)).select(VectorXi::Ones(rows), VectorXi::Zero(rows));
	}

	if (cols == 1)
	{
		retVec = tempMat.col(0);
	}
	else
	{
		for (int i = 0; i < cols - 1; ++i)
		{
			retVec = tempMat.col(i).array() * tempMat.col(i + 1).array();
		}
	}

	return retVec;
}


// 打印矩阵信息：
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
	{
		return;
	}
	if (colEnd > mat.cols() - 1 || colStart < 0 || colStart >= colEnd)
	{
		return;
	}

	std::cout << ": rows == " << mat.rows() << ",  cols == " << mat.cols() << std::endl;
	for (int i = rowStart; i <= rowEnd; ++i)
	{
		std::cout << i << " ---\t";
		for (int j = colStart; j <= colEnd; ++j)
		{
			std::cout << mat(i, j) << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


template<typename T>
void dispVec(const Matrix<T, Dynamic, 1>& vec)
{
	std::cout << ": rows == " << vec.rows() << std::endl;
	for (int i = 0; i < vec.rows(); ++i)
	{
		std::cout << i << "---\t" << vec(i) << std::endl;
	}
	std::cout << std::endl;
}

template<typename T>
void dispVec(const Matrix<T, 1, Dynamic>& vec)
{
	std::cout << ": cols == " << vec.cols() << std::endl;
	for (int i = 0; i < vec.cols(); ++i)
	{
		std::cout << i << "---\t" << vec(i) << std::endl;
	}
	std::cout << std::endl;
}


template<typename T>
void dispVecSeg(const Matrix<T, Dynamic, 1>& vec, const int start, const int end)
{
	if (start < 0 || end > vec.rows() - 1 || start >= end)
	{
		std::cout << "input error." << std::endl;
		return;
	}

	std::cout << ": rows == " << vec.rows() << std::endl;
	for (int i = start; i <= end; ++i)
	{
		std::cout << i << "---\t" << vec(i) << std::endl;
	}
	std::cout << std::endl;
}


template<typename T>
void dispVecSeg(const Matrix<T, 1, Dynamic>& vec, const int start, const int end)
{
	if (start < 0 || end > vec.cols() - 1 || start >= end)
	{
		std::cout << "input error." << std::endl;
		return;
	}

	std::cout << ": cols == " << vec.cols() << std::endl;
	for (int i = start; i <= end; ++i)
	{
		std::cout << i << "---\t" << vec(i) << std::endl;
	}
	std::cout << std::endl;
}

auto objReadMeshMat = [](MatrixXf& vers, MatrixXi& tris, const char* fileName)
{
	char* pTmp = NULL;
	std::ifstream ifs(fileName);//cube bunny Eight
	if (false == ifs.is_open())
	{
		return;
	}
	std::streampos   pos = ifs.tellg();     //   save   current   position   
	ifs.seekg(0, std::ios::end);
	unsigned fileLen = (unsigned)ifs.tellg();
	if (0 == fileLen)
	{
		return;
	}
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

auto objWriteMeshMat = [](const char* fileName, const MatrixXf& vers, const MatrixXi& tris)
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

auto objReadVerticesMat = [](MatrixXf& vers, const char* fileName)
{
	char* pTmp = NULL;
	std::ifstream ifs(fileName);			// cube bunny Eight
	if (false == ifs.is_open())
	{
		return;
	}

	std::streampos   pos = ifs.tellg();			//  save   current   position   
	ifs.seekg(0, std::ios::end);
	unsigned fileLen = (unsigned)ifs.tellg();
	if (0 == fileLen)
	{
		return;
	}

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

			vers.conservativeResize(3, ++cols);
			vers.col(cols - 1) = ver;
		}
		else
			break;
	}
	vers.transposeInPlace();
	delete[] pFileBuf;
};

auto objWriteVerticesMat = [](const char* fileName, const MatrixXf& vers)
{
	std::ofstream dstFile(fileName);
	for (int i = 0; i < vers.rows(); i++)
	{
		dstFile << "v " << vers(i, 0) << " " << vers(i, 1) << " " << vers(i, 2) << std::endl;
	}
	dstFile.close();
};

auto printDirEigen = [](const char* pathName, const RowVector3f& origin, const RowVector3f& dir)
{
	MatrixXf line;

	const float SR = 0.5;			// 空间采样率SR——相邻两个采样点的距离（单位mm）
	const float length = 10;
	int versCount = std::round(length / SR);
	line.resize(versCount + 1, 3);
	line.row(0) = origin;
	for (int i = 1; i <= versCount; i++)
	{
		line.row(i) = line.row(0) + SR * dir * i;
	}

	objWriteVerticesMat(pathName, line);
};


// 网格串联——合并两个孤立的网格到一个网格里
auto concatMeshMat = [](MatrixXf& vers, MatrixXi& tris, const MatrixXf& vers1, const MatrixXi& tris1)
{
	int versCount = vers.rows();
	matInsertRows<float>(vers, vers1);

	MatrixXi trisCopy1 = tris1;
	int* intPtr = trisCopy1.data();
	for (int i = 0; i<trisCopy1.size(); ++i) 
	{
		*(intPtr++) = versCount + *intPtr;
	}

	matInsertRows<int>(tris, trisCopy1);
};


// 解恰定的稠密线性方程组Ax == b;
template<typename T, size_t N>
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





