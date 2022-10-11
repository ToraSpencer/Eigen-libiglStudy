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


// 布尔向量转化为索引向量；
VectorXi flagVec2IdxVec(const VectorXi& flag)
{
	VectorXi idxVec;

	for (unsigned i = 0; i < flag.rows(); ++i)
	{
		if (0 != flag(i) && 1 != flag(i))
			return idxVec;
	}

	std::vector<int> tempVec;
	tempVec.reserve(flag.rows());
	for (unsigned i = 0; i< flag.rows(); ++i) 
	{
		if (flag(i) > 0)
			tempVec.push_back(i);
	}
	idxVec = vec2Vec<int>(tempVec);
	
	return idxVec;
}


// 索引向量转化为布尔向量；
VectorXi IdxVec2FlagVec(const VectorXi& idxVec, const unsigned size)
{
	VectorXi flag;

	if (idxVec.rows() > size)
		return flag;

	for (unsigned i = 0; i < idxVec.rows(); ++i) 
	{
		if (idxVec(i) >= size)
			return flag;
	}

	flag = VectorXi::Zero(size);
	for (unsigned i = 0; i < idxVec.rows(); ++i)
	{
		flag(idxVec(i)) = 1;
	}

	return flag;
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


// 输入起点、终点、空间采样率，插值生成一条直线点云；
bool interpolateToLine(MatrixXf& vers, const RowVector3f& start, const RowVector3f& end, const float SR, const bool SE)
{
	if (vers.rows() > 0)
		return false;

	RowVector3f dir = end - start;
	float length = dir.norm();
	dir.normalize();

	if (length <= SR)
		return true;

	if (SE)
		matInsertRows<float>(vers, start);

	float lenth0 = 0;
	for (unsigned i = 1; (length - lenth0) > 0.8 * SR; i++)			// 确保最后一个点距离终点不能太近。
	{
		RowVector3f temp = start + SR * i * dir;
		matInsertRows<float>(vers, temp);
		lenth0 = SR * (i + 1);		// 下一个temp的长度。
	}

	if (SE)
		matInsertRows<float>(vers, end);

	return true;
};


// 得到将originArrow旋转到targetArrow的旋转矩阵
Matrix3f getRotationMat(const RowVector3f& originArrow, const RowVector3f& targetArrow)
{
	Matrix3f rotation = Matrix3f::Zero();
	if (0 == originArrow.norm() || 0 == targetArrow.norm())
		return rotation;

	RowVector3f axisArrow = originArrow.cross(targetArrow);		// 旋转轴；

	if (0 == axisArrow.norm())
		return Matrix3f::Identity();

	axisArrow.normalize();
	float x0 = axisArrow(0);
	float y0 = axisArrow(1);
	float z0 = axisArrow(2);
	float cosTheta = originArrow.dot(targetArrow) / (originArrow.norm() * targetArrow.norm());
	float sinTheta = std::sqrt(1 - cosTheta*cosTheta);
	rotation << cosTheta + (1 - cosTheta) * x0 * x0, (1 - cosTheta)* x0* y0 - sinTheta * z0, (1 - cosTheta)* x0* z0 + sinTheta * y0, \
		(1 - cosTheta)* y0* x0 + sinTheta * z0, cosTheta + (1 - cosTheta) * y0 * y0, (1 - cosTheta)* y0* z0 - sinTheta * x0, \
		(1 - cosTheta)* z0* x0 - sinTheta * y0, (1 - cosTheta)* z0* y0 + sinTheta * x0, cosTheta + (1 - cosTheta) * z0 * z0;
	return rotation;
}


bool buildAdjacency(const Eigen::MatrixXi& tris, Eigen::MatrixXi& ttAdj_nmEdge, \
	std::vector<ttTuple>& ttAdj_nmnEdge, std::vector<ttTuple>& ttAdj_nmnOppEdge)
{
	/*
		bool buildAdjacency(
					const Eigen::MatrixXi& tris,												输入的三角片数据
					Eigen::MatrixXi& ttAdj_nmEdge,										三角片的非边缘流形有向边邻接的三角片索引；
					std::vector<ttTuple>& ttAdj_nmnEdge,							三角片的非流形有向边所在的三角片索引；
					std::vector<ttTuple>& ttAdj_nmnOppEdge					三角片的非流形有向边邻接的三角片索引（即对边所在的三角片的索引）；
					)

	*/
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

	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_weighted;
	adjSM_eCount.resize(versCount, versCount);
	adjSM_weighted.resize(versCount, versCount);
	adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());		// 权重为该有向边重复的次数；
	adjSM_weighted.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());		// 权重为该有向边的索引；
	Eigen::SparseMatrix<int> adjSM_weighted_opp = adjSM_weighted.transpose();

	Eigen::SparseMatrix<int> adjSM = adjSM_eCount;		// 有向边邻接矩阵；
	for (unsigned i = 0; i < adjSM.outerSize(); ++i)
	{
		for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, i); iter; ++iter)
		{
			if (iter.value() > 0)
				adjSM.coeffRef(iter.row(), iter.col()) = 1;
		}
	}


	// 2. 确定所有非边缘的流形有向边：
	Eigen::SparseMatrix<int> adjSM_eCount_ND = adjSM_eCount + Eigen::SparseMatrix<int>(adjSM_eCount.transpose());		// 权重为无向边ij关联的三角片数量；

	//		非边缘流形无向边：
	Eigen::SparseMatrix<int> adjSM_MNnonBdry_ND = adjSM_eCount_ND;
	for (unsigned i = 0; i < adjSM_MNnonBdry_ND.outerSize(); ++i)
	{
		for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_MNnonBdry_ND, i); iter; ++iter)
		{
			if (iter.value() > 0)
			{
				if (iter.value() == 2)
					adjSM_MNnonBdry_ND.coeffRef(iter.row(), iter.col()) = 1;
				else
					adjSM_MNnonBdry_ND.coeffRef(iter.row(), iter.col()) = 0;
			}
		}
	}

	//		非边缘流形有向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MNnonBdry = adjSM + adjSM_MNnonBdry_ND;		// adjSM & adjSM_MNnonBdry_ND
	adjSM_MNnonBdry.prune([](const Eigen::Index& row, const Eigen::Index& col, const float& value)->bool
		{
			if (2 == value)
				return true;
			else
				return false;
		});
	adjSM_MNnonBdry /= 2;

	unsigned edgesCount_MNnonBdry = adjSM_MNnonBdry.sum();
	Eigen::MatrixXi edges_MNnonBdry(edgesCount_MNnonBdry, 2);		// 非边缘流形有向边
	std::vector<int> edgesIdx_MNnonBdry;												// 非边缘流形有向边的索引；
	edgesIdx_MNnonBdry.reserve(edgesCount_MNnonBdry);
	unsigned index = 0;
	for (unsigned i = 0; i < adjSM_MNnonBdry.outerSize(); ++i)
	{
		for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_MNnonBdry, i); iter; ++iter)
		{
			edges_MNnonBdry(index, 0) = iter.row();
			edges_MNnonBdry(index, 1) = iter.col();
			edgesIdx_MNnonBdry.push_back(adjSM_weighted.coeffRef(iter.row(), iter.col()));
			index++;
		}
	}

	Eigen::SparseMatrix<int> adjSM_MNnonBdry_opp = adjSM_MNnonBdry.transpose();
	Eigen::MatrixXi edges_MNnonBdry_opp(edgesCount_MNnonBdry, 2);		// 非边缘流形有向边的对边；
	std::vector<int> edgesIdx_MNnonBdry_opp;								// 非边缘流形有向边的对边的索引；
	edgesIdx_MNnonBdry_opp.reserve(edgesCount_MNnonBdry);
	index = 0;
	for (unsigned i = 0; i < adjSM_MNnonBdry_opp.outerSize(); ++i)
	{
		for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_MNnonBdry_opp, i); iter; ++iter)
		{
			edges_MNnonBdry_opp(index, 0) = iter.row();
			edges_MNnonBdry_opp(index, 1) = iter.col();
			edgesIdx_MNnonBdry_opp.push_back(adjSM_weighted_opp.coeffRef(iter.row(), iter.col()));
			index++;
		}
	}

	// 3. 流形边-三角片邻接关系，三角片邻接关系：
	std::vector<int> etInfo(edgesCount);				// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
	for (int i = 0; i < edgesCount; ++i)
		etInfo[i] = i % trisCount;

	Eigen::VectorXi etAdj_mnEdge(-Eigen::VectorXi::Ones(edgesCount));	// 边所在三角片的索引，非流形边或边缘边写-1；
	for (unsigned i = 0; i < edgesIdx_MNnonBdry.size(); ++i)
	{
		const int& edgeIdx = edgesIdx_MNnonBdry[i];
		const int& edgesIdxOpp = edgesIdx_MNnonBdry_opp[i];
		etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
	}

	//			三角片邻接矩阵，ttAdj_mnEdge(i, :)是索引为i的三角片三条边邻接的三个三角片，若其中有非流形边或边缘边则写为-1；
	ttAdj_nmEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);

	// 4. 非流形边（关联两个三角片的有向边）信息
	index = 0;
	std::vector<int> edgeIdx_nmn;
	edgeIdx_nmn.reserve(edgesCount);
	for (int i = 0; i < edgesCount; ++i)
		if (2 == adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)))
			edgeIdx_nmn.push_back(i);
	edgeIdx_nmn.shrink_to_fit();


	//			找出同一条非流形有向边对应的多个边索引：map<pair表示的非流形边数据， 边索引vector>
	std::unordered_map<std::pair<int, int>, std::vector<int>, edgeHash, edgeComparator> edges_nmn_map;
	for (const auto& eIdx : edgeIdx_nmn)
	{
		std::pair<int, int> edge{ edges(eIdx, 0), edges(eIdx, 1) };
		auto retPair = edges_nmn_map.insert({ edge , std::vector<int>{eIdx} });
		if (!retPair.second)		// 若插入失败，则说明已有此键；
		{
			auto iter = edges_nmn_map.find(edge);
			iter->second.push_back(eIdx);
		}
	}

	//			map<非流形边Idx, 该边的对边对应的多个边索引的vector>
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_map;
	for (const auto& eIdx : edgeIdx_nmn)
	{
		std::pair<int, int> edge{ edges(eIdx, 0), edges(eIdx, 1) };
		const std::vector<int>& commonEidxes = edges_nmn_map.find(edge)->second;
		edgeIdx_nmn_map.insert({ eIdx, commonEidxes });
	}
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_opp_map;
	for (const auto& pair : edgeIdx_nmn_map)
	{
		int eIdx = pair.first;
		int vaIdx = edges(eIdx, 0);
		int vbIdx = edges(eIdx, 1);
		int eOppIdx = -1;
		for (int i = 0; i < edgesCount; ++i)
			if (edges(i, 0) == vbIdx && edges(i, 1) == vaIdx)
			{
				eOppIdx = i;
				break;
			}
		auto iter = edgeIdx_nmn_map.find(eOppIdx);
		edgeIdx_nmn_opp_map.insert({ eIdx, iter->second });
	}


	// 5. 含有非流形边的三角片的邻接关系：
	ttAdj_nmnEdge.resize(trisCount);
	for (const auto& pair : edgeIdx_nmn_map)
	{
		int eIdx = pair.first;
		int row = eIdx % trisCount;
		int col = eIdx / trisCount;

		std::vector<int> trisIdx_nmn = pair.second;
		for (auto& index : trisIdx_nmn)
			index = etInfo[index];

		switch (col)
		{
		case 0:	std::get<0>(ttAdj_nmnEdge[row]) = trisIdx_nmn;
			break;
		case 1:	std::get<1>(ttAdj_nmnEdge[row]) = trisIdx_nmn;
			break;
		case 2:	std::get<2>(ttAdj_nmnEdge[row]) = trisIdx_nmn;
			break;
		default:
			return false;
		}
	}

	ttAdj_nmnOppEdge.resize(trisCount);
	for (const auto& pair : edgeIdx_nmn_opp_map)
	{
		int eIdx = pair.first;
		int row = eIdx % trisCount;
		int col = eIdx / trisCount;

		std::vector<int> trisIdx_nmn_opp = pair.second;
		for (auto& index : trisIdx_nmn_opp)
			index = etInfo[index];

		switch (col)
		{
		case 0:	std::get<0>(ttAdj_nmnOppEdge[row]) = trisIdx_nmn_opp;
			break;
		case 1:	std::get<1>(ttAdj_nmnOppEdge[row]) = trisIdx_nmn_opp;
			break;
		case 2:	std::get<2>(ttAdj_nmnOppEdge[row]) = trisIdx_nmn_opp;
			break;
		default:
			return false;
		}
	}


	return true;
}