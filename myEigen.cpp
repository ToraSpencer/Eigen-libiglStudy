#include "myEigen.h"
 
using namespace std;
using namespace Eigen;
 
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
 
// 文件中读写矩阵
void objReadMeshMat(MatrixXf& vers, MatrixXi& tris, const char* fileName)
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

void objWriteMeshMat(const char* fileName, const MatrixXf& vers, const MatrixXi& tris)
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

void objReadVerticesMat(MatrixXf& vers, const char* fileName)
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

void objWriteVerticesMat(const char* fileName, const MatrixXf& vers)
{
	std::ofstream dstFile(fileName);
	for (int i = 0; i < vers.rows(); i++)
	{
		dstFile << "v " << vers(i, 0) << " " << vers(i, 1) << " " << vers(i, 2) << std::endl;
	}
	dstFile.close();
};

// 读取点云OBJ文件中的数据，存储到齐次坐标系表示的点云矩阵中；注：顶点在矩阵中是列表示的，第四维的元素始终为1；
void objReadVerticesHomoMat(MatrixXf& vers, const char* fileName)

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
	int rows = vers.cols();
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

			Vector4f ver(Vector4f::Ones());
			ver(0) = (float)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(1) = (float)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(2) = (float)atof(tmpBuffer);

			vers.conservativeResize(++rows, 4);
			vers.row(rows - 1) = ver;
		}
		else
			break;
	}
	vers.transposeInPlace();
	delete[] pFileBuf;
};


// 齐次坐标系点云数据写入到OBJ文件中；注：顶点在矩阵中是列表示的，第四维的元素始终为1；
void objWriteVerticesHomoMat(const char* fileName, const MatrixXf& vers)
{
	std::ofstream dstFile(fileName);
	for (int i = 0; i < vers.cols(); i++)
	{
		dstFile << "v " << vers(0, i) << " " << vers(1, i) << " " << vers(2, i) << std::endl;
	}
	dstFile.close();

}


//	普通点云矩阵转换为齐次坐标系下的点云矩阵：
void vers2homoVers(MatrixXf& homoVers, const MatrixXf& vers)
{
	homoVers = MatrixXf::Ones(4, vers.rows());
	homoVers.topRows(3) = vers.transpose();
}

MatrixXf vers2homoVers(const MatrixXf& vers)
{
	MatrixXf homoVers = MatrixXf::Ones(4, vers.rows());
	homoVers.topRows(3) = vers.transpose();
	return homoVers;
}


// 齐次坐标系下的点云矩阵变换为普通点云矩阵
void homoVers2vers(MatrixXf& vers, const MatrixXf& homoVers)
{
	MatrixXf tempMat = homoVers.transpose();
	vers = tempMat.leftCols(3);
}

MatrixXf homoVers2vers(const MatrixXf& homoVers)
{
	MatrixXf tempMat = homoVers.transpose();
	MatrixXf vers = tempMat.leftCols(3);
	return vers;
}


void printDirEigen(const char* pathName, const RowVector3f& origin, const RowVector3f& dir)
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

void printCoordinateEigen(const char* pathName, const RowVector3f& origin, const RowVector3f& xdir, \
	const RowVector3f& ydir, const RowVector3f& zdir)
{
	const float SR = 0.5;			// 空间采样率SR——相邻两个采样点的距离（单位mm）
	const float length = 10;
	int versCount = std::round(length / SR);
	MatrixXf line1(versCount, 3), line2(versCount, 3), line3(versCount, 3); 
	for (int i = 0; i < versCount; i++)
	{
		line1.row(i) = origin + SR * xdir * (i + 1);
	}
	for (int i = 0; i < versCount; i++)
	{
		line2.row(i) = origin + SR * ydir * (i + 1);
	}
	for (int i = 0; i < versCount; i++)
	{
		line3.row(i) = origin + SR * zdir * (i + 1);
	}

	MatrixXf line = origin;
	matInsertRows<float>(line, line1);
	matInsertRows<float>(line, line2);
	matInsertRows<float>(line, line3);
	objWriteVerticesMat(pathName, line);
}

// 网格串联——合并两个孤立的网格到一个网格里
void concatMeshMat(MatrixXf& vers, MatrixXi& tris, const MatrixXf& vers1, const MatrixXi& tris1)
{
	int versCount = vers.rows();
	matInsertRows<float>(vers, vers1);

	MatrixXi trisCopy1 = tris1;
	int* intPtr = trisCopy1.data();
	for (int i = 0; i < trisCopy1.size(); ++i)
	{
		*(intPtr++) = versCount + *intPtr;
	}

	matInsertRows<int>(tris, trisCopy1);
};
 

// 多项式插值
void polyInterpolation() 
{

}


// 高斯插值
void gaussInterpolation() {}



// 最小二乘多项式拟合曲线：
void leastSquarePolyFitting()
{}



 // 岭回归多项式拟合曲线
void ridgeRegressionPolyFitting(VectorXf& theta, const MatrixXf& vers)
{
	/*
		void ridgeRegressionPolyFitting(
					VectorXf & theta,				拟合的多项式函数
					const MatrixXf & vers			离散样本点
					)
	*/

	const float lambda = 0.1;
	int m = 4;									// 多项式最多项数

	int n = vers.rows();
	if (n == 0)
	{
		return;
	}
	if (m >= n)
	{
		m = n - 1;
	}

	Eigen::MatrixXf X(n, m);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			X(i, j) = std::powf(vers(i, 0), j);
		}
	}

	Eigen::VectorXf Y = vers.col(1);
	Eigen::MatrixXf I(m, m);
	I.setIdentity();
	theta = (X.transpose() * X + I * lambda).inverse() * X.transpose() * Y;
}


