#include "myEigenIO.h"

static const double pi = 3.14159265359;


unsigned readNextData(char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize)
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
	for (int j = 0; j < vers.rows(); j++)
	{
		char szBuf[1024] = { 0 };
		sprintf_s(szBuf, 1024, "v %.17g %.17g %.17g", vers(j, 0), vers(j, 1), vers(j, 2));
		dstFile << szBuf << "\n";
	}
	dstFile.close();
}


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
		dstFile << "l " << edges.coeffRef(i, 0) + static_cast<typename DerivedI::Scalar>(1) << " " \
		<< edges.coeffRef(i, 1) + static_cast<typename DerivedI::Scalar>(1) << std::endl;

	dstFile.close();
}


// objWritePath() 路径数据写入到OBJ文件中：
template <typename DerivedV, typename	 IndexType>
void objWritePath(const char* pathName, const std::vector<IndexType>& path, \
	const Eigen::PlainObjectBase<DerivedV>& vers)
{
	if (path.size() <= 1)
		return;

	unsigned edgesCount = path.size() - 1;
	Eigen::MatrixXi pathEdges(edgesCount, 2);

	for (unsigned i = 0; i < edgesCount; ++i)
		pathEdges(i, 0) = static_cast<int>(path[i]);
	for (unsigned i = 0; i < edgesCount; ++i)
		pathEdges(i, 1) = static_cast<int>(path[i + 1]);

	objWriteEdgesMat(pathName, pathEdges, vers);
}


// objWirteTreePath()——输入树向量或路径向量，写入到OBJ文件中：
template <typename DerivedV>
void objWriteTreePath(const char* pathName, const Eigen::VectorXi& treeVec, const Eigen::PlainObjectBase<DerivedV>& vers)
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
 

void objWriteDirection(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& dir)
{
	Eigen::MatrixXf line;

	const float SR = 0.5;			// 空间采样率SR——相邻两个采样点的距离（单位mm）
	const float length = 10;
	int versCount = std::round(length / SR);
	line.resize(versCount + 1, 3);
	line.row(0) = origin;
	for (int i = 1; i <= versCount; i++)
		line.row(i) = line.row(0) + SR * dir * i;

	objWriteVerticesMat(pathName, line);
};


void objWriteCoorSys(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& xdir, \
	const Eigen::RowVector3f& ydir, const Eigen::RowVector3f& zdir)
{
	const float SR = 0.5;			// 空间采样率SR——相邻两个采样点的距离（单位mm）
	const float length = 10;
	int versCount = std::round(length / SR);
	Eigen::MatrixXf line1(versCount, 3), line2(versCount, 3), line3(versCount, 3);
	for (int i = 0; i < versCount; i++)
		line1.row(i) = origin + SR * xdir * (i + 1);
	for (int i = 0; i < versCount; i++)
		line2.row(i) = origin + SR * ydir * (i + 1);
	for (int i = 0; i < versCount; i++)
		line3.row(i) = origin + SR * zdir * (i + 1);

	Eigen::MatrixXf line = origin;
	int rows = line.rows();
	line.conservativeResize(line.rows() + line1.rows() + line2.rows() + line3.rows(), 3);
	for (int i = 0; i < line1.rows(); ++i)
		line.row(rows + i) = line1.row(i);
	rows = line.rows();
	for (int i = 0; i < line2.rows(); ++i)
		line.row(rows + i) = line2.row(i);
	rows = line.rows();
	for (int i = 0; i < line3.rows(); ++i)
		line.row(rows + i) = line3.row(i);
	 
	objWriteVerticesMat(pathName, line);
}



// 模板函数需要特化之后才能在静态库中输出： 

/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化：
template void vecWriteToFile<int>(const char* fileName, const std::vector<int>& vec);
template void vecWriteToFile<float>(const char* fileName, const std::vector<float>& vec);
template void vecWriteToFile<double>(const char* fileName, const std::vector<double>& vec);

template void vecReadFromFile<int>(std::vector<int>& vec, const char* fileName, const unsigned elemCount);
template void vecReadFromFile<float>(std::vector<float>& vec, const char* fileName, const unsigned elemCount);
template void vecReadFromFile<double>(std::vector<double>& vec, const char* fileName, const unsigned elemCount);

template bool matWriteToFile<int>(const char* fileName, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& mat);
template bool matWriteToFile<float>(const char* fileName, const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& mat);
template bool matWriteToFile<double>(const char* fileName, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mat);

template bool matReadFromFile<int>(Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName);
template bool matReadFromFile<float>(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName);
template bool matReadFromFile<double>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName);

template void objReadMeshMat<float, int>(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::Matrix<int, Eigen::Dynamic, \
	Eigen::Dynamic>& tris, const char* fileName);
template void objReadMeshMat<double, int>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::Matrix<int, Eigen::Dynamic, \
	Eigen::Dynamic>& tris, const char* fileName);

template void objWriteMeshMat<float>(const char* fileName, const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const Eigen::MatrixXi& tris);
template void objWriteMeshMat<double>(const char* fileName, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const Eigen::MatrixXi& tris);

template void objReadVerticesMat<float>(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, const char* fileName);
template void objReadVerticesMat<double>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, const char* fileName);

template	void objWriteVerticesMat<Eigen::Matrix<float, -1, -1, 0, -1, -1>>(const char* fileName, \
	const Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>& vers);
template	void objWriteVerticesMat<Eigen::Matrix<double, -1, -1, 0, -1, -1>>(const char* fileName, \
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>& vers);
template	void objWriteVerticesMat<Eigen::Matrix<float, 1, 3, 1, 1, 3>>(const char* fileName, \
	const Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3>>& vers);
template	void objWriteVerticesMat<Eigen::Matrix<double, 1, 3, 1, 1, 3>>(const char* fileName, \
	const Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3>>& vers);

template void objWriteEdgesMat(const char* pathName, const Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>>& edges, \
	const Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>& vers);
template void objWriteEdgesMat(const char* pathName, const Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>>& edges, \
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>& vers);
template void objWriteEdgesMat(const char* pathName, const Eigen::PlainObjectBase<Eigen::RowVector2i>& edges, \
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>& vers);

template void objWritePath(const char* pathName, const std::vector<int>& path, \
	const Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>& vers);
template void objWritePath(const char* pathName, const std::vector<int>& path, \
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>& vers);
template void objWritePath(const char* pathName, const std::vector<unsigned>& path, \
	const Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>& vers);
template void objWritePath(const char* pathName, const std::vector<unsigned>& path, \
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>& vers);


template void objWriteTreePath(const char* pathName, const Eigen::VectorXi& treeVec, \
	const Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>& vers);
template void objWriteTreePath(const char* pathName, const Eigen::VectorXi& treeVec, \
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>& vers);