#include "myEigenHomoCoor.h"

 
// 读取点云OBJ文件中的数据，存储到齐次坐标系表示的点云矩阵中；注：顶点在矩阵中是列表示的，第四维的元素始终为1；
void objReadVerticesHomoMat(Eigen::MatrixXf& vers, const char* fileName)
{
	auto readNextData = [](char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize)
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

			Eigen::Vector4f ver(Eigen::Vector4f::Ones());
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
void objWriteVerticesHomoMat(const char* fileName, const Eigen::MatrixXf& vers)
{
	std::ofstream dstFile(fileName);
	for (int i = 0; i < vers.cols(); i++) 
		dstFile << "v " << vers(0, i) << " " << vers(1, i) << " " << vers(2, i) << std::endl; 

	dstFile.close();
}


// 齐次坐标系下的点云矩阵变换为普通点云矩阵
void homoVers2vers(Eigen::MatrixXf& vers, const Eigen::MatrixXf& homoVers)
{
	Eigen::MatrixXf tempMat = homoVers.transpose();
	vers = tempMat.leftCols(3);
}


Eigen::MatrixXf homoVers2vers(const Eigen::MatrixXf& homoVers)
{
	Eigen::MatrixXf tempMat = homoVers.transpose();
	Eigen::MatrixXf vers = tempMat.leftCols(3);
	return vers;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板函数实现：

template	<typename T>
void objWriteHomoMeshMat(const char* fileName, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& homoVers, const Eigen::MatrixXi& tris)
{
	std::ofstream dstFile(fileName);
	if (homoVers.rows() != 4 || tris.cols() != 3)
		return;

	for (int j = 0; j < homoVers.cols(); j++)
	{
		char szBuf[1024] = { 0 };
		sprintf_s(szBuf, 1024, "v %.17g %.17g %.17g", homoVers(0, j), homoVers(1, j), homoVers(2, j));
		dstFile << szBuf << "\n";
	}

	for (unsigned j = 0; j < tris.rows(); ++j)
	{
		char szBuf[256] = { 0 };
		sprintf_s(szBuf, 256, "f %d %d %d", tris(j, 0) + 1, tris(j, 1) + 1, tris(j, 2) + 1);
		dstFile << szBuf << "\n";
	}
};


template <typename Derived>
Eigen::MatrixXf vers2homoVersF(const Eigen::PlainObjectBase<Derived>& vers)
{
	const int versCount = vers.rows();
	Eigen::MatrixXf tmpVers = vers.transpose().cast<float>();

	Eigen::MatrixXf homoVers(4, versCount);
	homoVers.setOnes();
	homoVers.topRows(3) = tmpVers;

	return homoVers;
}


template <typename Derived>
Eigen::MatrixXd vers2homoVersD(const Eigen::PlainObjectBase<Derived>& vers)
{
	const int versCount = vers.rows();
	Eigen::MatrixXd tmpVers = vers.transpose().cast<double>();

	Eigen::MatrixXd homoVers(4, versCount);
	homoVers.setOnes();
	homoVers.topRows(3) = tmpVers;

	return homoVers;
}






#include "templateSpecialization.cpp"