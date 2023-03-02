#include "myEigen.h"
 
 
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
void objReadVerticesHomoMat(Eigen::MatrixXf& vers, const char* fileName)
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
	{
		dstFile << "v " << vers(0, i) << " " << vers(1, i) << " " << vers(2, i) << std::endl;
	}
	dstFile.close();

}


//	普通点云矩阵转换为齐次坐标系下的点云矩阵：
void vers2homoVers(Eigen::MatrixXf& homoVers, const Eigen::MatrixXf& vers)
{
	homoVers = Eigen::MatrixXf::Ones(4, vers.rows());
	homoVers.topRows(3) = vers.transpose();
}


Eigen::MatrixXf vers2homoVers(const Eigen::MatrixXf& vers)
{
	Eigen::MatrixXf homoVers = Eigen::MatrixXf::Ones(4, vers.rows());
	homoVers.topRows(3) = vers.transpose();
	return homoVers;
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


void printDirEigen(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& dir)
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


void printCoordinateEigen(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& xdir, \
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
	matInsertRows<float>(line, line1);
	matInsertRows<float>(line, line2);
	matInsertRows<float>(line, line3);
	objWriteVerticesMat(pathName, line);
}


// 布尔向量转化为索引向量；
Eigen::VectorXi flagVec2IdxVec(const Eigen::VectorXi& flag)
{
	Eigen::VectorXi idxVec;

	for (unsigned i = 0; i < flag.rows(); ++i)
		if (0 != flag(i) && 1 != flag(i))
			return idxVec;

	std::vector<int> tempVec;
	tempVec.reserve(flag.rows());
	for (unsigned i = 0; i< flag.rows(); ++i) 
		if (flag(i) > 0)
			tempVec.push_back(i);
	idxVec = vec2Vec<int>(tempVec);
	
	return idxVec;
}


// 索引向量转化为布尔向量；
Eigen::VectorXi IdxVec2FlagVec(const Eigen::VectorXi& idxVec, const unsigned size)
{
	Eigen::VectorXi flag;

	if (idxVec.rows() > size)
		return flag;

	for (unsigned i = 0; i < idxVec.rows(); ++i) 
		if (idxVec(i) >= size)
			return flag;

	flag = Eigen::VectorXi::Zero(size);
	for (unsigned i = 0; i < idxVec.rows(); ++i)
		flag(idxVec(i)) = 1;

	return flag;
}


// 多项式插值
void polyInterpolation() 
{

}


// 高斯插值
void gaussInterpolation() {}


// 最小二乘多项式拟合曲线：
void leastSquarePolyFitting()
{}








Eigen::VectorXd fittingStandardEllipse(const Eigen::MatrixXf& sampleVers)
{
	/*
		Eigen::VectorXd fittingStandardEllipse(								// 返回列向量(a,c,d,e,f)，为标准椭圆方程的系数；
					const Eigen::MatrixXf& sampleVers					// 输入的样本点，必须是在XOY平面上的点；
		)
	*/
	const double epsilon = 1e-8;				// 浮点数绝对值小于此值时认为为0；

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


// 计算网格中边的三角片的邻接关系（非流形有向边关联的三角片最多只能为两个）：
bool buildAdjacency(const Eigen::MatrixXi& tris, Eigen::MatrixXi& ttAdj_nmEdge, \
	std::vector<ttTuple>& ttAdj_nmnEdge, std::vector<ttTuple>& ttAdj_nmnOppEdge)
{
	/*
		bool buildAdjacency(
					const Eigen::MatrixXi& tris,									输入的三角片数据

					Eigen::MatrixXi& ttAdj_nmEdge,							三角片的非边缘流形有向边邻接的三角片索引；
																									trisCount * 3
																									(i, 0)为索引为i的三角片的非边缘流形有向边(vb, vc)邻接的三角片的索引；

					std::vector<ttTuple>& ttAdj_nmnEdge,				三角片的非流形有向边所在的三角片索引；
																									size == trisCount;
																									每个ttTuple元素包含三个int向量；
																									索引为i的元素的第0个向量，为索引为i的三角片的非流形有向边(vb, vc)所在的所有三角片索引；


					std::vector<ttTuple>& ttAdj_nmnOppEdge		三角片的非流形有向边邻接的三角片索引（即对边所在的三角片的索引）；
																									size == trisCount;
																									索引为i的元素的第0个向量，为索引为i的三角片的非流形有向边(vb, vc)邻接的所有三角片索引；
					)

	*/
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	Eigen::MatrixXi edges;							// 有向边数据；
	Eigen::MatrixXi nmnEdges;					// 非流形有向边
	std::vector<int> etInfo;						// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；

	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_map;				// 非流形有向边的索引 - 该边对应的所有索引；
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_opp_map;		// 非流形有向边的索引 - 该边的对边对应的所有索引；

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_weighted的triplet数据；

	Eigen::SparseMatrix<int> adjSM;								// 邻接矩阵；
	Eigen::SparseMatrix<int> adjSM_eCount;					// 邻接矩阵，索引为该有向边重复的次数；
	Eigen::SparseMatrix<int> adjSM_weighted;				// 邻接矩阵，元素若对应流形边则为该边索引，若对应非流形边则为该边所有索引的和；
	Eigen::SparseMatrix<int> adjSM_weighted_opp;		// adjSM_weighted的转置，表示对边的信息；

	Eigen::SparseMatrix<int> adjSM_ueCount;				// 权重为无向边ij关联的三角片数量；
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;

	std::vector<int> edgesIdx_MN_NB;							// 非边缘流形有向边的索引；
	std::vector<int> edgesIdx_MN_NB_opp;					// 非边缘流形有向边的对边的索引；

	std::unordered_multimap<std::int64_t, int> edgeMap;			// 边编码——64位整型数数表示的边数据（两个端点索引）；

	Eigen::VectorXi etAdj_mnEdge;							// 边所在三角片的索引，非流形边或边缘边写-1；

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


	// 4. 计算非流形边（关联两个三角片的有向边）信息
	{
		//	4.1 遍历邻接矩阵找出所有非流形边（关联两个三角片的有向边）；
		std::vector<int> edgeIdx_nmn;						// 非流形有向边索引；
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (2 == adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)))
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();

		int neCount = edgeIdx_nmn.size();
		subFromIdxVec(nmnEdges, edges, edgeIdx_nmn);

		// 4.2 建立边编码-边索引的哈希表；
		for (int i = 0; i < edges.rows(); ++i)
		{
			std::int64_t eCode = encodeEdge(edges(i, 0), edges(i, 1));
			edgeMap.insert({ eCode, i });
		}

		// 4.3 在哈希表中搜索所有非流形边，建立非流形边及其对边的（边——边索引）映射关系；
		for (int i = 0; i < neCount; ++i)
		{
			int neIdx = edgeIdx_nmn[i];		// 当前非流形边索引；
			int vaIdx = nmnEdges(i, 0);
			int vbIdx = nmnEdges(i, 1);
			std::int64_t eCode = encodeEdge(vaIdx, vbIdx);
			std::int64_t oppEcode = encodeEdge(vbIdx, vaIdx);
			auto iter = edgeMap.find(eCode);
			auto oppIter = edgeMap.find(oppEcode);
			int otherIdx = (iter->second == neIdx) ? ((++iter)->second) : (iter->second);
			int oppIdx1 = oppIter->second;
			int oppIdx2 = (++oppIter)->second;
			edgeIdx_nmn_map.insert({ neIdx, std::vector<int>{neIdx, otherIdx} });
			edgeIdx_nmn_opp_map.insert({ neIdx, std::vector<int>{oppIdx1, oppIdx2} });
		}
	}


	// 5. 求含有非流形边的三角片的邻接关系：
	{
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
	}

	return true;
}


std::pair<int, int> decodeEdge(const std::int64_t code)
{
	int a = static_cast<int>(code >> 32);
	int b = static_cast<int>(code - (static_cast<std::int64_t>(a) << 32));
	return std::make_pair(a, b);
}


std::vector<int> decodeTrianagle(const std::uint64_t code)
{
	std::uint64_t a = code >> 42;
	std::uint64_t resi = code - (a << 42);
	std::uint64_t b = resi >> 21;
	std::uint64_t c = resi - (b << 21);
	return std::vector<int>{static_cast<int>(a), static_cast<int>(b), static_cast<int>(c)};
}


// 生成T_MESH::di_cell对应的轴向包围盒网格；
void genAABBmesh(const T_MESH::di_cell& cell, Eigen::MatrixXd& vers, Eigen::MatrixXi& tris)
{
	vers.resize(0, 0);
	tris.resize(0, 0);
	T_MESH::Point arrow = cell.Mp - cell.mp;
	T_MESH::Point newOri = (cell.mp + cell.Mp) / 2.0;

	genCubeMesh(vers, tris);
	vers.col(0) *= arrow.x;
	vers.col(1) *= arrow.y;
	vers.col(2) *= arrow.z;
	vers.rowwise() += Eigen::RowVector3d{ newOri.x, newOri.y, newOri.z };
}



/////////////////////////////////////////////////////////////////////////////////////////////////DEBUG接口：
static std::string g_debugPath = "E:/";


template<typename DerivedV>
static void debugWriteVers(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteVerticesMat(path, vers);
}


static void debugWriteMesh(const char* name, T_MESH::Basic_TMesh& mesh)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	mesh.save(path);
}


template<typename T>
static void debugWriteMesh(const char* name, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteMeshMat(path, vers, tris);
}
 


// 测试myEigen中的接口
namespace TEST_MYEIGEN 
{
	// 测试编码解码：
	void test0() 
	{
		std::vector<int> retVec = decodeTrianagle(encodeTriangle(65535, 65534, 65533));
		traverseSTL(retVec, disp<int>);
		std::cout << std::hex << encodeTriangle(65535, 65534, 65533) << std::endl;
		std::cout << std::dec << encodeTriangle(1, 1, 2) << std::endl;;

		retVec.clear();
		retVec = decodeTrianagle(encodeTriangle(651135, 522534, 653533));
		traverseSTL(retVec, disp<int>);

		std::cout << "finished." << std::endl;
	}


	// 测试区域生长triangleGrow();
	void test1() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/材料/meshInnerRev.obj");
		objWriteMeshMat("E:/triangleGrowInput.obj", vers, tris);

		std::vector<Eigen::MatrixXd> meshesVersOut;
		std::vector<Eigen::MatrixXi> meshesTrisOut;
		triangleGrowSplitMesh(meshesVersOut, meshesTrisOut, vers, tris);

		for (int i = 0; i<meshesVersOut.size(); ++i) 
		{
			char str[256];
			sprintf_s(str, "E:/triangleGrowOut%d.obj", i);
			objWriteMeshMat(str, meshesVersOut[i], meshesTrisOut[i]);
		}

		//triangleGrow(versOut, trisOut, vers, tris, 0);
		//objWriteMeshMat("E:/triangleGrowOut.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// 检测重复三角片，重复顶点，
	void test2() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/材料/shrinkedMeshDirty.obj");

		Eigen::MatrixXi bdrys;
		bdryEdges(bdrys, tris);
		std::cout << "boundary edges count is " << bdrys.rows() << std::endl;
		objWriteEdgesMat("E:/bdrys.obj", bdrys, vers);

		std::vector<int> repIdxes;
		findRepTris(repIdxes, tris);


		std::cout << "finished." << std::endl;
	}


	// 测试获取三角网格基础属性的接口：
	void test3() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		bool retFlag = true;
		objReadMeshMat(vers, tris, "E:/材料/jawCore.obj");

		const unsigned versCount = vers.rows();
		const unsigned trisCount = tris.rows();

		// 计算网格每个三角片的面积：
		Eigen::VectorXd trisAreaVec;
		if (!trisArea(trisAreaVec, vers, tris))
			return;

		// 手动计算第一个三角片的面积：
		Eigen::RowVector3d arrow1, arrow2, dir1, dir2, va, vb, vc;
		va = vers.row(tris(0,0));
		vb = vers.row(tris(0, 1));
		vc = vers.row(tris(0, 2));
		arrow1 = vb - va;
		arrow2 = vc - va;
		dir1 = arrow1.normalized();
		dir2 = arrow2.normalized();
		double cosTheta = dir1.dot(dir2);
		double sinTheta = std::sqrt(1 - cosTheta*cosTheta);
		double area0 = 0.5 * arrow1.norm() * arrow2.norm() * sinTheta;
		std::cout << "area0 == " << area0 << std::endl;
		dispVecSeg(trisAreaVec, 0, 9);

		std::cout << "finished." << std::endl;
	}


	// 测试基本数学接口
	void test4()
	{
		dispPair(cart2polar(5.0, -5.0));
		dispPair(polar2cart(20.0, 3.14159/3));

		std::cout << "finished." << std::endl;
	}


	// 测试生成基础图形的接口：
	void test5() 
	{
		Eigen::MatrixXf axis(15, 3);
		axis.setZero();
		double deltaTheta = pi / 10;
		for (unsigned i = 0; i < 15; ++i)
		{
			double theta = deltaTheta * i;
			axis(i, 0) = 50 * cos(theta);
			axis(i, 1) = 50 * sin(theta);
		}
 
		debugWriteVers("axis", axis);
		std::cout << "finished." << std::endl;

		// 三角剖分：
		Eigen::MatrixXd surfVers, circleVers;
		Eigen::MatrixXi surfTris;
		objReadVerticesMat(circleVers, "E:/材料/circleVers.obj");
		circuit2mesh(surfVers, surfTris, circleVers);
		debugWriteMesh("surfMesh", surfVers, surfTris);

		// 生成柱体：
		Eigen::MatrixXf cylinderVers;
		Eigen::MatrixXi cylinderTris;
		genCylinder(cylinderVers, cylinderTris, axis, 10);				// 生成圆柱
		debugWriteMesh("cylinder", cylinderVers, cylinderTris);
		cylinderVers.resize(0, 0);
		cylinderTris.resize(0, 0);
		genCylinder(cylinderVers, cylinderTris, axis, std::make_pair(10.0, 20.0));
		debugWriteMesh("pillar", cylinderVers, cylinderTris);



		std::cout << "finished." << std::endl;
	}

	
	// 测试空间变换相关的接口：
	void test6()
	{
		Eigen::Matrix3d rotation;
		Eigen::RowVector3d axisArrow;
		Eigen::RowVector3d arrow1, arrow2, arrowRotted;
		float theta = pi;

		// 绕任意轴旋转：
		theta = pi / 3;
		axisArrow = Eigen::RowVector3d{20, 0, 0};
		rotation = getRotationMat(axisArrow, theta);
		arrowRotted = arrow1 * rotation.transpose();
		dispVec<double, 3>(arrowRotted);

		// 旋转到目标向量：
		arrow1.setRandom();
		arrow2.setRandom();
		arrow1.normalize();
		arrow2.normalize();
		dispVec<double, 3>(arrow1);
		dispVec<double, 3>(arrow2);
		rotation = getRotationMat(arrow1, arrow2);
		arrowRotted = arrow1 * rotation.transpose();
		dispVec<double, 3>(arrowRotted);

		std::cout << "finished." << std::endl;
	}

}



// 测试图像处理：
namespace TEST_DIP
{
	void test0()
	{
		Eigen::MatrixXi mat1(9, 9);
		mat1.setRandom();
		
		Eigen::MatrixXd mask(3, 3);
		mask.setOnes();

		Eigen::MatrixXi matOut;
		linearSpatialFilter(matOut, mat1, mask);

		dispMat(mat1);
		dispMat(matOut);

		std::cout << "finished." << std::endl;
	}

}

 

// 测试拓扑网格类tmesh
namespace TEST_TMESH
{
	// TMESH的IO，基本功能
	void test0()
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();												// ？？？This is mandatory
		T_MESH::Basic_TMesh tmesh1;
		T_MESH::Node* nodePtr = nullptr;
		tmesh1.load("E:/材料/tooth.obj");
		tmesh1.save("E:/meshInput.obj");

		const int versCount = tmesh1.V.numels();
		const int edgesCount = tmesh1.E.numels();
		const int trisCount = tmesh1.T.numels();

		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		TMesh2MeshMat(tmesh1, vers, tris);
		debugWriteMesh("meshInputCopy", vers, tris);

		// 顶点数据：
		nodePtr = tmesh1.V.head();
		nodePtr = nodePtr->next();					// 第二个顶点
		T_MESH::Vertex* vPtr = reinterpret_cast<T_MESH::Vertex*>(nodePtr->data);
		T_MESH::Vertex& ver0 = *vPtr;

		// 边数据；
		nodePtr =tmesh1.E.head();
		T_MESH::Edge* ePtr = reinterpret_cast<T_MESH::Edge*>(nodePtr->data);
		T_MESH::Edge& e0 = *ePtr;

		// 三角片数据：
		nodePtr = tmesh1.T.head();
		T_MESH::Triangle* tPtr = reinterpret_cast<T_MESH::Triangle*>(nodePtr->data);
		T_MESH::Triangle& t0 = *tPtr;
 
		std::cout << "finished." << std::endl;
	}


	// 测试T_MESH::di_cell类
	void test1() 
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh mesh;
		T_MESH::Node* nodePtr = nullptr;
		mesh.load("E:/材料/tooth.obj");				// 必须是流形网格；
		debugWriteMesh("meshInput", mesh);

		T_MESH::di_cell cell(&mesh);			// 输入网格，生成di_cell对象；
		Eigen::MatrixXd aabbVers;
		Eigen::MatrixXi aabbTris;

		// 打印cell对应的包围盒；
		genAABBmesh(cell, aabbVers, aabbTris);
		debugWriteMesh("aabb", aabbVers, aabbTris);

		// fork()方法：
		T_MESH::di_cell cell2 = *cell.fork();
		genAABBmesh(cell, aabbVers, aabbTris);
		debugWriteMesh("aabb11", aabbVers, aabbTris);
		genAABBmesh(cell2, aabbVers, aabbTris);
		debugWriteMesh("aabb12", aabbVers, aabbTris);

		std::cout << "finished." << std::endl;
	}


	// 测试TMESH中的几何元素标记功能：
	void test2() 
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/材料/tooth.obj");
		debugWriteMesh("meshInput", tMesh);

		// 1.遍历顶点，标记x坐标小于0的顶点；
		traverseVersList(tMesh.V, [&](T_MESH::Vertex* verPtr) 
			{
				if (verPtr->x < 0)
					MARK_VISIT(verPtr);
			});

		// 2. 选中所有包含选中顶点的三角片；
		tMesh.growSelection();

		// 3. 删除所有选中三角片，及其包含的点和边；
		tMesh.removeSelectedTriangles();
		debugWriteMesh("meshOut", tMesh);

		std::cout << "finished." << std::endl;
	}


	// 拓扑性质、拓扑修复
	void test3() 
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/材料/meshArranged.obj");
		debugWriteMesh("meshInput", tMesh);

		// 检测网格拓扑性质：
		const char* retStr = tMesh.checkConnectivity();				// 若有拓扑问题则会返回字符串；
		if (NULL != retStr)
			std::cout << retStr << std::endl;

		std::cout << "finished." << std::endl;
	}
}