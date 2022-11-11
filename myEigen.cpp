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


VectorXd fittingStandardEllipse(const MatrixXf& sampleVers)
{
	/*
		VectorXd fittingStandardEllipse(								// 返回列向量(a,c,d,e,f)，为标准椭圆方程的系数；
					const MatrixXf& sampleVers					// 输入的样本点，必须是在XOY平面上的点；
		)
	*/
	const double epsilon = 1e-8;				// 浮点数绝对值小于此值时认为为0；

	// 标准椭圆方程：a*x^2 + c*y^2 + d*x + e*y + f  = 0，其中a*c > 0
	VectorXd x(VectorXd::Zero(5));

	unsigned m = sampleVers.rows();				// sample count;
	VectorXd x0(VectorXd::Zero(m));
	VectorXd y0(VectorXd::Zero(m));
	for (unsigned i = 0; i < m; ++i)
	{
		x0(i) = static_cast<double>(sampleVers(i, 0));
		y0(i) = static_cast<double>(sampleVers(i, 1));
	}

	// alpha = [x^2, y^2, x, y, 1]; 样本信息矩阵：A = [alpha1; alpha2; .... alpham]; 椭圆方程写为：A*x = 0;
	MatrixXd A = MatrixXd::Ones(m, 5);
	A.col(0) = x0.array() * x0.array();
	A.col(1) = y0.array() * y0.array();
	A.col(2) = x0;
	A.col(3) = y0;

	MatrixXd ATA = A.transpose().eval() * A;
	MatrixXd B(MatrixXd::Zero(5, 5));
	B(0, 1) = 1;
	B(1, 0) = 1;
	MatrixXd S = ATA.inverse() * B.transpose();

	// 求S的特征值，特征向量：
	EigenSolver<MatrixXd> es(S);
	MatrixXd D = es.pseudoEigenvalueMatrix();			// 对角线元素是特征值
	MatrixXd V = es.pseudoEigenvectors();				// 每一个列向量都是特征向量。

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
					const Eigen::MatrixXi& tris,												输入的三角片数据
					Eigen::MatrixXi& ttAdj_nmEdge,										三角片的非边缘流形有向边邻接的三角片索引；
					std::vector<ttTuple>& ttAdj_nmnEdge,							三角片的非流形有向边所在的三角片索引；
					std::vector<ttTuple>& ttAdj_nmnOppEdge					三角片的非流形有向边邻接的三角片索引（即对边所在的三角片的索引）；
					)

	*/
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

#ifdef LOCAL_DEBUG
	tiktok& tt = tiktok::getInstance();
#endif

	Eigen::MatrixXi edges;							// 有向边数据；
	std::vector<int> etInfo;						// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；

	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_map;				// 非流形有向边的索引 - 该边对应的所有索引；
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_opp_map;		// 非流形有向边的索引 - 该边的对边对应的所有索引；

	Eigen::SparseMatrix<int> adjSM;								// 邻接矩阵；
	Eigen::SparseMatrix<int> adjSM_eCount;					// 邻接矩阵，索引为该有向边重复的次数；
	Eigen::SparseMatrix<int> adjSM_weighted;				// 邻接矩阵，元素若对应流形边则为该边索引，若对应非流形边则为该边所有索引的和；
	Eigen::SparseMatrix<int> adjSM_weighted_opp;		// adjSM_weighted的转置，表示对边的信息；

	Eigen::SparseMatrix<int> adjSM_eCount_ND;				// 权重为无向边ij关联的三角片数量；
	Eigen::SparseMatrix<int> adjSM_MNnonBdry_ND;
	Eigen::SparseMatrix<int> adjSM_MNnonBdry;
	Eigen::SparseMatrix<int> adjSM_MNnonBdry_opp;

	std::vector<int> edgesIdx_MNnonBdry;							// 非边缘流形有向边的索引；
	std::vector<int> edgesIdx_MNnonBdry_opp;					// 非边缘流形有向边的对边的索引；

	// 1. 求基本的边信息、邻接矩阵：
	{
#ifdef LOCAL_DEBUG
		tt.start();
#endif

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

		std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;
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
		for (unsigned i = 0; i < adjSM.outerSize(); ++i)
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, i); iter; ++iter)
				if (iter.value() > 0)
					adjSM.coeffRef(iter.row(), iter.col()) = 1;

		// 若adjSM_weighted(i, j)对应的是流形有向边，该权重值为该边的索引；若是非流形有向边，则该边对应两个边索引，权重为两个索引之和；
		adjSM_weighted.resize(versCount, versCount);
		adjSM_weighted.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());
		adjSM_weighted_opp = adjSM_weighted.transpose();

		// 若存在重复次数大于2的有向边，则输入不合法；
		for (unsigned i = 0; i < adjSM_eCount.outerSize(); ++i)
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_eCount, i); iter; ++iter)
				if (iter.value() > 2)
					return false;

#ifdef LOCAL_DEBUG
		tt.endCout("Elapsed time of calculating adjacency matrices is : ");
#endif
	}


	// 2. 确定所有非边缘的流形有向边的索引，及其对边的索引；
	{
#ifdef LOCAL_DEBUG
		tt.start();
#endif
		//		所有无向边：
		adjSM_eCount_ND = adjSM_eCount + Eigen::SparseMatrix<int>(adjSM_eCount.transpose());		// 权重为无向边ij关联的三角片数量；

		//		非边缘流形无向边：
		adjSM_MNnonBdry_ND = adjSM_eCount_ND;
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

		//		非边缘流形有向边、及其对边的邻接矩阵
		adjSM_MNnonBdry = adjSM + adjSM_MNnonBdry_ND;			// adjSM & adjSM_MNnonBdry_ND
		adjSM_MNnonBdry.prune([](const Eigen::Index& row, const Eigen::Index& col, const float& value)->bool
			{
				if (2 == value)			// 删除非流形边
					return true;
				else
					return false;
			});
		adjSM_MNnonBdry /= 2;
		adjSM_MNnonBdry_opp = adjSM_MNnonBdry.transpose();

		unsigned edgesCount_MNnonBdry = adjSM_MNnonBdry.sum();
		edgesIdx_MNnonBdry.reserve(edgesCount_MNnonBdry);					// 非边缘流形有向边索引；
		edgesIdx_MNnonBdry_opp.reserve(edgesCount_MNnonBdry);

		for (unsigned i = 0; i < adjSM_MNnonBdry.outerSize(); ++i)
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_MNnonBdry, i); iter; ++iter)
				edgesIdx_MNnonBdry.push_back(adjSM_weighted.coeffRef(iter.row(), iter.col()));

		for (unsigned i = 0; i < adjSM_MNnonBdry_opp.outerSize(); ++i)
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_MNnonBdry_opp, i); iter; ++iter)
				edgesIdx_MNnonBdry_opp.push_back(adjSM_weighted_opp.coeffRef(iter.row(), iter.col()));

#ifdef LOCAL_DEBUG
		tt.endCout("Elapsed time of step 2 is : ");
#endif
	}


	// 3. 流形边-三角片邻接关系，三角片邻接关系：
	{
#ifdef LOCAL_DEBUG
		tt.start();
#endif

		// 3.1 生成边索引 - 三角片索引映射表etInfo;
		etInfo.resize(edgesCount);					// etInfo(i)是索引为i的边所在的三角片的索引；
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

		// 3.2 求边所在三角片的索引；
		Eigen::VectorXi etAdj_mnEdge(-Eigen::VectorXi::Ones(edgesCount));	// 边所在三角片的索引，非流形边或边缘边写-1；
		for (unsigned i = 0; i < edgesIdx_MNnonBdry.size(); ++i)
		{
			const int& edgeIdx = edgesIdx_MNnonBdry[i];
			const int& edgesIdxOpp = edgesIdx_MNnonBdry_opp[i];
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//	3.3 求三角片邻接矩阵；ttAdj_mnEdge(i, :)是索引为i的三角片三条边邻接的三个三角片，若其中有非流形边或边缘边则写为-1；
		ttAdj_nmEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);

#ifdef LOCAL_DEBUG
		tt.endCout("Elapsed time of step 3 is : ");
#endif
	}


	// 4. 计算非流形边（关联两个三角片的有向边）信息
	{
#ifdef LOCAL_DEBUG
		tt.start();
#endif

		//	4.1 遍历邻接矩阵找出所有非流形边（关联两个三角片的有向边）；
		std::vector<int> edgeIdx_nmn;						// 非流形有向边索引；
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (2 == adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)))
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();
		int neCount = edgeIdx_nmn.size();
		Eigen::MatrixXi nmnEdges;
		subFromIdxVec(nmnEdges, edges, edgeIdx_nmn);

		// 4.2 建立边-边索引的哈希表；
		std::unordered_multimap<double, int> edgeMap;			// 浮点数表示的边数据——边索引；
		for (int i = 0; i < edges.rows(); ++i)
		{
			double key = edges(i, 0) + 1e-10 * edges(i, 1);
			edgeMap.insert({ key, i });
		}

		// 4.3 在哈希表中搜索所有非流形边，建立非流形边及其对边的（边——边索引）映射关系；
		for (int i = 0; i < neCount; ++i)
		{
			int neIdx = edgeIdx_nmn[i];		// 当前非流形边索引；
			int vaIdx = nmnEdges(i, 0);
			int vbIdx = nmnEdges(i, 1);
			double key = vaIdx + 1e-10 * vbIdx;
			double oppKey = vbIdx + 1e-10 * vaIdx;
			auto iter = edgeMap.find(key);
			auto oppIter = edgeMap.find(oppKey);
			int otherIdx = (iter->second == neIdx) ? ((++iter)->second) : (iter->second);
			int oppIdx1 = oppIter->second;
			int oppIdx2 = (++oppIter)->second;
			edgeIdx_nmn_map.insert({ neIdx, std::vector<int>{neIdx, otherIdx} });
			edgeIdx_nmn_opp_map.insert({ neIdx, std::vector<int>{oppIdx1, oppIdx2} });
		}

#ifdef LOCAL_DEBUG
		tt.endCout("elapsed time of step 4 is : ");
#endif
	}


	// 5. 求含有非流形边的三角片的邻接关系：
	{
#ifdef LOCAL_DEBUG
		tt.start();
#endif

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
#ifdef LOCAL_DEBUG
		tt.endCout("Elapsed time of step 5 is : ");
#endif
	}

	return true;
}