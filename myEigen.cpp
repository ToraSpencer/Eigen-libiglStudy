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


// ��ȡ����OBJ�ļ��е����ݣ��洢���������ϵ��ʾ�ĵ��ƾ����У�ע�������ھ��������б�ʾ�ģ�����ά��Ԫ��ʼ��Ϊ1��
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

		// ������Ϣ		
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


// �������ϵ��������д�뵽OBJ�ļ��У�ע�������ھ��������б�ʾ�ģ�����ά��Ԫ��ʼ��Ϊ1��
void objWriteVerticesHomoMat(const char* fileName, const MatrixXf& vers)
{
	std::ofstream dstFile(fileName);
	for (int i = 0; i < vers.cols(); i++)
	{
		dstFile << "v " << vers(0, i) << " " << vers(1, i) << " " << vers(2, i) << std::endl;
	}
	dstFile.close();

}


//	��ͨ���ƾ���ת��Ϊ�������ϵ�µĵ��ƾ���
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


// �������ϵ�µĵ��ƾ���任Ϊ��ͨ���ƾ���
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

	const float SR = 0.5;			// �ռ������SR������������������ľ��루��λmm��
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
	const float SR = 0.5;			// �ռ������SR������������������ľ��루��λmm��
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


// ��������ת��Ϊ����������
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


// ��������ת��Ϊ����������
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


// ����ʽ��ֵ
void polyInterpolation() 
{

}


// ��˹��ֵ
void gaussInterpolation() {}


// ��С���˶���ʽ������ߣ�
void leastSquarePolyFitting()
{}


 // ��ع����ʽ�������
void ridgeRegressionPolyFitting(VectorXf& theta, const MatrixXf& vers)
{
	/*
		void ridgeRegressionPolyFitting(
					VectorXf & theta,				��ϵĶ���ʽ����
					const MatrixXf & vers			��ɢ������
					)
	*/

	const float lambda = 0.1;
	int m = 4;									// ����ʽ�������

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


// ������㡢�յ㡢�ռ�����ʣ���ֵ����һ��ֱ�ߵ��ƣ�
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
	for (unsigned i = 1; (length - lenth0) > 0.8 * SR; i++)			// ȷ�����һ��������յ㲻��̫����
	{
		RowVector3f temp = start + SR * i * dir;
		matInsertRows<float>(vers, temp);
		lenth0 = SR * (i + 1);		// ��һ��temp�ĳ��ȡ�
	}

	if (SE)
		matInsertRows<float>(vers, end);

	return true;
};


// �õ���originArrow��ת��targetArrow����ת����
Matrix3f getRotationMat(const RowVector3f& originArrow, const RowVector3f& targetArrow)
{
	Matrix3f rotation = Matrix3f::Zero();
	if (0 == originArrow.norm() || 0 == targetArrow.norm())
		return rotation;

	RowVector3f axisArrow = originArrow.cross(targetArrow);		// ��ת�᣻

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
		VectorXd fittingStandardEllipse(								// ����������(a,c,d,e,f)��Ϊ��׼��Բ���̵�ϵ����
					const MatrixXf& sampleVers					// ����������㣬��������XOYƽ���ϵĵ㣻
		)
	*/
	const double epsilon = 1e-8;				// ����������ֵС�ڴ�ֵʱ��ΪΪ0��

	// ��׼��Բ���̣�a*x^2 + c*y^2 + d*x + e*y + f  = 0������a*c > 0
	VectorXd x(VectorXd::Zero(5));

	unsigned m = sampleVers.rows();				// sample count;
	VectorXd x0(VectorXd::Zero(m));
	VectorXd y0(VectorXd::Zero(m));
	for (unsigned i = 0; i < m; ++i)
	{
		x0(i) = static_cast<double>(sampleVers(i, 0));
		y0(i) = static_cast<double>(sampleVers(i, 1));
	}

	// alpha = [x^2, y^2, x, y, 1]; ������Ϣ����A = [alpha1; alpha2; .... alpham]; ��Բ����дΪ��A*x = 0;
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

	// ��S������ֵ������������
	EigenSolver<MatrixXd> es(S);
	MatrixXd D = es.pseudoEigenvalueMatrix();			// �Խ���Ԫ��������ֵ
	MatrixXd V = es.pseudoEigenvectors();				// ÿһ����������������������

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


// ���������бߵ�����Ƭ���ڽӹ�ϵ������������߹���������Ƭ���ֻ��Ϊ��������
bool buildAdjacency(const Eigen::MatrixXi& tris, Eigen::MatrixXi& ttAdj_nmEdge, \
	std::vector<ttTuple>& ttAdj_nmnEdge, std::vector<ttTuple>& ttAdj_nmnOppEdge)
{
	/*
		bool buildAdjacency(
					const Eigen::MatrixXi& tris,									���������Ƭ����

					Eigen::MatrixXi& ttAdj_nmEdge,							����Ƭ�ķǱ�Ե����������ڽӵ�����Ƭ������
																									trisCount * 3
																									(i, 0)Ϊ����Ϊi������Ƭ�ķǱ�Ե���������(vb, vc)�ڽӵ�����Ƭ��������

					std::vector<ttTuple>& ttAdj_nmnEdge,				����Ƭ�ķ�������������ڵ�����Ƭ������
																									size == trisCount;
																									ÿ��ttTupleԪ�ذ�������int������
																									����Ϊi��Ԫ�صĵ�0��������Ϊ����Ϊi������Ƭ�ķ����������(vb, vc)���ڵ���������Ƭ������


					std::vector<ttTuple>& ttAdj_nmnOppEdge		����Ƭ�ķ�����������ڽӵ�����Ƭ���������Ա����ڵ�����Ƭ����������
																									size == trisCount;
																									����Ϊi��Ԫ�صĵ�0��������Ϊ����Ϊi������Ƭ�ķ����������(vb, vc)�ڽӵ���������Ƭ������
					)

	*/
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	Eigen::MatrixXi edges;							// ��������ݣ�
	Eigen::MatrixXi nmnEdges;					// �����������
	std::vector<int> etInfo;						// ������ - ����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������

	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_map;				// ����������ߵ����� - �ñ߶�Ӧ������������
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_opp_map;		// ����������ߵ����� - �ñߵĶԱ߶�Ӧ������������

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// ��������ϡ�����adjSM_eCount, adjSM_weighted��triplet���ݣ�

	Eigen::SparseMatrix<int> adjSM;								// �ڽӾ���
	Eigen::SparseMatrix<int> adjSM_eCount;					// �ڽӾ�������Ϊ��������ظ��Ĵ�����
	Eigen::SparseMatrix<int> adjSM_weighted;				// �ڽӾ���Ԫ������Ӧ���α���Ϊ�ñ�����������Ӧ�����α���Ϊ�ñ����������ĺͣ�
	Eigen::SparseMatrix<int> adjSM_weighted_opp;		// adjSM_weighted��ת�ã���ʾ�Աߵ���Ϣ��

	Eigen::SparseMatrix<int> adjSM_ueCount;				// Ȩ��Ϊ�����ij����������Ƭ������
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;

	std::vector<int> edgesIdx_MN_NB;							// �Ǳ�Ե��������ߵ�������
	std::vector<int> edgesIdx_MN_NB_opp;					// �Ǳ�Ե��������ߵĶԱߵ�������

	std::unordered_multimap<std::int64_t, int> edgeMap;			// �߱��롪��64λ����������ʾ�ı����ݣ������˵���������

	Eigen::VectorXi etAdj_mnEdge;							// ����������Ƭ�������������α߻��Ե��д-1��

	// 1. ������ı���Ϣ���ڽӾ���
	{

		// edges == [ea; eb; ec] == [vbIdxes, vcIdxes; vcIdxes, vaIdxes; vaIdxes, vbIdxes];
		/*
			������������������edges�б���ÿ���߶���unique�ģ�
			�����ڷ����αߣ�������α���edges����ظ��洢��ʵ�ʱ�����Ҳ������edges������
			����˵ʹ������representation�Ļ��������α��в�ֹһ��������ȡ�������÷����αߵ�����Ƭ������
		*/

		// ����Ƭ�����ߵ�������teIdx == [eaIdx, ebIdx, ecIdx] == [(0: trisCount-1)', (trisCount: 2*trisCount-1)', (2*trisCount, 3*trisCount-1)'];
		/*
			teIdx(i, j) = trisCount *j + i;
		*/

		// 1.1 �������������
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

		// 1.2 ��������������ڽӾ���adjSM, adjSM_eCount, smElems_weighted
		smElems.reserve(edgesCount);
		smElems_weighted.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
		{
			smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});
			smElems_weighted.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), i});
		}

		adjSM_eCount.resize(versCount, versCount);
		adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());	// Ȩ��Ϊ��������ظ��Ĵ�����
		adjSM = adjSM_eCount;																		// ������ڽӾ���
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		��adjSM_weighted(i, j)��Ӧ������������ߣ���Ȩ��ֵΪ�ñߵ����������Ƿ���������ߣ���ñ߶�Ӧ������������Ȩ��Ϊ��������֮�ͣ�
		adjSM_weighted.resize(versCount, versCount);
		adjSM_weighted.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());
		adjSM_weighted_opp = adjSM_weighted.transpose();

		//	1.3 ����������ڽӾ���
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// ������ڽӾ���
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.4 ���ɷǱ�Ե����������ڽӾ���
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

		//	1.5 ���ɷǱ�Ե��������ߡ�����Աߵ��ڽӾ���
		adjSM_MN_NB = adjSM + adjSM_ueMN_NB;					// adjSM & adjSM_ueMN_NB
		adjSM_MN_NB.prune([](const Eigen::Index& row, const Eigen::Index& col, const int& value)->bool
			{
				if (2 == value)			// ɾ�������α�
					return true;
				else
					return false;
			});
		adjSM_MN_NB /= 2;
		adjSM_MN_NB_opp = adjSM_MN_NB.transpose();
	}

	// 2. ȷ���Ǳ�Ե��������ߡ�����Աߵ�������
	{
		unsigned edgesCount_MN_NB = adjSM_MN_NB.sum();
		edgesIdx_MN_NB.reserve(edgesCount_MN_NB);					// �Ǳ�Ե���������������
		edgesIdx_MN_NB_opp.reserve(edgesCount_MN_NB);			// �Ǳ�Ե��������ߵĶԱߵ�������

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB.push_back(adjSM_weighted.coeffRef(iter.row(), iter.col()));
			});

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB_opp.push_back(adjSM_weighted_opp.coeffRef(iter.row(), iter.col()));
			});
	}

	// 2+.ͨ���߼������жϵ�ǰ�����Ƿ�������
	{
#ifdef LOCAL_DEBUG
		tt.start();

		std::unordered_map<std::int64_t, unsigned> eCountMap;
		std::unordered_map<std::int64_t, unsigned> ueCountMap;

		// ��������߼����ڽӱ�
		traverseSparseMatrix(adjSM_eCount, [&](auto& iter)
			{
				std::int64_t edgeCode = encodeEdge(iter.row(), iter.col());
				eCountMap.insert({ edgeCode, iter.value() });
			});

		// ��������߼����ڽӱ�
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				std::int64_t edgeCode = encodeEdge(iter.row(), iter.col());
				ueCountMap.insert({ edgeCode, iter.value() });
			});

		// ���:
		std::vector<int> sickEdgeVerIdxes;
		std::vector<std::pair<unsigned, unsigned>> sickCount;
		sickEdgeVerIdxes.reserve(this->m_meshVers.size());
		sickCount.reserve(this->m_meshVers.size());
		for (const auto& pair : ueCountMap)
		{
			std::int64_t edgeCode = pair.first;
			unsigned ueCount = pair.second;
			unsigned eCount = eCountMap[edgeCode];
			bool mnFlag = (eCount == 1 && ueCount == 2);						// ���α�
			bool normalNMNflag = (eCount == 2 && ueCount == 4);		// ����ķ����αߣ�

			if (!(mnFlag || normalNMNflag))			// ��Ե�ߣ�(eCount == 1 && ueCount == 1), (eCount == 0, ueCount == 1)Ҳ��Ϊ������
			{
				std::pair<int, int> sickVers = decodeEdge(edgeCode);
				sickEdgeVerIdxes.push_back(sickVers.first);
				sickEdgeVerIdxes.push_back(sickVers.second);
				sickCount.push_back({ eCount, ueCount });
			}
		}
		sickEdgeVerIdxes.shrink_to_fit();
		sickCount.shrink_to_fit();

		if (!sickEdgeVerIdxes.empty())
		{
			// ��ӡ�쳣�ı�
			Eigen::VectorXi sickEdgeVersVec = vec2Vec<int, -1>(sickEdgeVerIdxes);
			Eigen::MatrixXi sickEdges = Eigen::Map<Eigen::MatrixXi>(sickEdgeVersVec.data(), 2, sickEdgeVersVec.size() / 2).transpose();
			objWriteEdgesMat("E:/sickEdges.obj", sickEdges, this->m_meshVers);
			std::cout << "unnormal edge detected!" << std::endl;
			tt.endCout("Elapsed time of step 2+ is : ");
			return false;
		}

		tt.endCout("Elapsed time of step 2+ is : ");
#endif
	}


	// 3. ���α�-����Ƭ�ڽӹ�ϵ������Ƭ�ڽӹ�ϵ��
	{
		// 3.1 ���ɱ����� - ����Ƭ����ӳ���etInfo;
		etInfo.resize(edgesCount);							// etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

		// 3.2 �����������Ƭ��������
		etAdj_mnEdge = -Eigen::VectorXi::Ones(edgesCount);					// ����������Ƭ�������������α߻��Ե��д-1��
		for (unsigned i = 0; i < edgesIdx_MN_NB.size(); ++i)
		{
			const int& edgeIdx = edgesIdx_MN_NB[i];
			const int& edgesIdxOpp = edgesIdx_MN_NB_opp[i];
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//	3.3 ������Ƭ�ڽӾ���ttAdj_mnEdge(i, :)������Ϊi������Ƭ�������ڽӵ���������Ƭ���������з����α߻��Ե����дΪ-1��
		ttAdj_nmEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. ��������αߣ�������������Ƭ������ߣ���Ϣ
	{
		//	4.1 �����ڽӾ����ҳ����з����αߣ�������������Ƭ������ߣ���
		std::vector<int> edgeIdx_nmn;						// �����������������
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (2 == adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)))
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();

		int neCount = edgeIdx_nmn.size();
		subFromIdxVec(nmnEdges, edges, edgeIdx_nmn);

		// 4.2 �����߱���-�������Ĺ�ϣ��
		for (int i = 0; i < edges.rows(); ++i)
		{
			std::int64_t eCode = encodeEdge(edges(i, 0), edges(i, 1));
			edgeMap.insert({ eCode, i });
		}

		// 4.3 �ڹ�ϣ�����������з����αߣ����������α߼���Աߵģ��ߡ�����������ӳ���ϵ��
		for (int i = 0; i < neCount; ++i)
		{
			int neIdx = edgeIdx_nmn[i];		// ��ǰ�����α�������
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


	// 5. ���з����αߵ�����Ƭ���ڽӹ�ϵ��
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
	std::uint64_t a = code >> 32;
	std::uint64_t resi = code - (a << 32);
	std::uint64_t b = resi >> 16;
	std::uint64_t c = resi - (b << 16);
	return std::vector<int>{static_cast<int>(a), static_cast<int>(b), static_cast<int>(c)};
}



namespace TEST_MYEIGEN 
{
	// ���Ա�����룺
	void test0() 
	{
		std::vector<int> retVec = decodeTrianagle(encodeTriangle(65535, 65534, 65533));
		traverseSTL(retVec, disp<int>);
		std::cout << std::hex << encodeTriangle(65535, 65534, 65533) << std::endl;
		std::cout << std::dec << encodeTriangle(1, 1, 2) << std::endl;;


		std::cout << "finished." << std::endl;
	}


	// ������������triangleGrow();
	void test1() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/����/meshInnerRev.obj");
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


}
