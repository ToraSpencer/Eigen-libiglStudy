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


// �����������ϲ���������������һ��������
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


bool buildAdjacency(const Eigen::MatrixXi& tris, Eigen::MatrixXi& ttAdj_nmEdge, \
	std::vector<ttTuple>& ttAdj_nmnEdge, std::vector<ttTuple>& ttAdj_nmnOppEdge)
{
	/*
		bool buildAdjacency(
					const Eigen::MatrixXi& tris,												���������Ƭ����
					Eigen::MatrixXi& ttAdj_nmEdge,										����Ƭ�ķǱ�Ե����������ڽӵ�����Ƭ������
					std::vector<ttTuple>& ttAdj_nmnEdge,							����Ƭ�ķ�������������ڵ�����Ƭ������
					std::vector<ttTuple>& ttAdj_nmnOppEdge					����Ƭ�ķ�����������ڽӵ�����Ƭ���������Ա����ڵ�����Ƭ����������
					)

	*/
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	// 1. �󶥵��ڽӹ�ϵ��

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
	adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());		// Ȩ��Ϊ��������ظ��Ĵ�����
	adjSM_weighted.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());		// Ȩ��Ϊ������ߵ�������
	Eigen::SparseMatrix<int> adjSM_weighted_opp = adjSM_weighted.transpose();

	Eigen::SparseMatrix<int> adjSM = adjSM_eCount;		// ������ڽӾ���
	for (unsigned i = 0; i < adjSM.outerSize(); ++i)
	{
		for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, i); iter; ++iter)
		{
			if (iter.value() > 0)
				adjSM.coeffRef(iter.row(), iter.col()) = 1;
		}
	}


	// 2. ȷ�����зǱ�Ե����������ߣ�
	Eigen::SparseMatrix<int> adjSM_eCount_ND = adjSM_eCount + Eigen::SparseMatrix<int>(adjSM_eCount.transpose());		// Ȩ��Ϊ�����ij����������Ƭ������

	//		�Ǳ�Ե��������ߣ�
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

	//		�Ǳ�Ե����������ڽӾ���
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
	Eigen::MatrixXi edges_MNnonBdry(edgesCount_MNnonBdry, 2);		// �Ǳ�Ե���������
	std::vector<int> edgesIdx_MNnonBdry;												// �Ǳ�Ե��������ߵ�������
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
	Eigen::MatrixXi edges_MNnonBdry_opp(edgesCount_MNnonBdry, 2);		// �Ǳ�Ե��������ߵĶԱߣ�
	std::vector<int> edgesIdx_MNnonBdry_opp;								// �Ǳ�Ե��������ߵĶԱߵ�������
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

	// 3. ���α�-����Ƭ�ڽӹ�ϵ������Ƭ�ڽӹ�ϵ��
	std::vector<int> etInfo(edgesCount);				// ������ - ����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
	for (int i = 0; i < edgesCount; ++i)
		etInfo[i] = i % trisCount;

	Eigen::VectorXi etAdj_mnEdge(-Eigen::VectorXi::Ones(edgesCount));	// ����������Ƭ�������������α߻��Ե��д-1��
	for (unsigned i = 0; i < edgesIdx_MNnonBdry.size(); ++i)
	{
		const int& edgeIdx = edgesIdx_MNnonBdry[i];
		const int& edgesIdxOpp = edgesIdx_MNnonBdry_opp[i];
		etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
	}

	//			����Ƭ�ڽӾ���ttAdj_mnEdge(i, :)������Ϊi������Ƭ�������ڽӵ���������Ƭ���������з����α߻��Ե����дΪ-1��
	ttAdj_nmEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);

	// 4. �����αߣ�������������Ƭ������ߣ���Ϣ
	index = 0;
	std::vector<int> edgeIdx_nmn;
	edgeIdx_nmn.reserve(edgesCount);
	for (int i = 0; i < edgesCount; ++i)
		if (2 == adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)))
			edgeIdx_nmn.push_back(i);
	edgeIdx_nmn.shrink_to_fit();


	//			�ҳ�ͬһ������������߶�Ӧ�Ķ����������map<pair��ʾ�ķ����α����ݣ� ������vector>
	std::unordered_map<std::pair<int, int>, std::vector<int>, edgeHash, edgeComparator> edges_nmn_map;
	for (const auto& eIdx : edgeIdx_nmn)
	{
		std::pair<int, int> edge{ edges(eIdx, 0), edges(eIdx, 1) };
		auto retPair = edges_nmn_map.insert({ edge , std::vector<int>{eIdx} });
		if (!retPair.second)		// ������ʧ�ܣ���˵�����д˼���
		{
			auto iter = edges_nmn_map.find(edge);
			iter->second.push_back(eIdx);
		}
	}

	//			map<�����α�Idx, �ñߵĶԱ߶�Ӧ�Ķ����������vector>
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


	// 5. ���з����αߵ�����Ƭ���ڽӹ�ϵ��
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