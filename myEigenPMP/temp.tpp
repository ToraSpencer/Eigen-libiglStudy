
/////////////////////////////////////////////////////////////////////////////////////////////////// ͼ�����ԣ�



// �õ�����߱���-������������ֵ䣺
template <typename DerivedI>
bool getEdgeIdxMap(std::unordered_multimap<std::int64_t, int>& map, const Eigen::PlainObjectBase<DerivedI>& edges)
{
	int edgesCount = edges.rows();
	for (int i = 0; i < edgesCount; ++i)
	{
		std::int64_t code = encodeEdge(edges(i, 0), edges(i, 1));
		map.insert({ code, i });
	}

	return true;
}


template <typename DerivedI>
bool getUedges(Eigen::MatrixXi& uEdges, const Eigen::PlainObjectBase<DerivedI>& edges)
{
	unsigned edgesCount = edges.rows();
	uEdges.resize(0, 0);
	if (0 == edges.size() || 2 != edges.cols())
		return false;

	std::unordered_set<std::int64_t> uEdgeSet;
	for (unsigned i = 0; i < edgesCount; ++i)
		uEdgeSet.insert(encodeUedge(edges(i, 0), edges(i, 1)));
	unsigned uEdgesCount = uEdgeSet.size();
	uEdges.resize(uEdgesCount, 2);

	int index = 0;
	for (auto& code : uEdgeSet)
	{
		auto pair = decodeEdge(code);
		uEdges(index, 0) = pair.first;
		uEdges(index, 1) = pair.second;
		index++;
	}

	return true;
}


template <typename DerivedI>
bool getUedges(Eigen::MatrixXi& uEdges, Eigen::VectorXi& edgeUeInfo, Eigen::MatrixXi& ueEdgeInfo, const Eigen::PlainObjectBase<DerivedI>& edges)
{
	// ���������������񲻿�ʹ�ã�
	/*
		bool getUedges(
				Eigen::MatrixXi& uEdges,													�����
				Eigen::VectorXi& edgeUeInfo,												edgeUeInfo(i)������Ϊi������߶�Ӧ������ߵ�������
				Eigen::MatrixXi& ueEdgeInfo												ueCount * 2�ľ���
																													ueEdgeInfo(i, 0)��ueEdgeInfo(i, 1)������Ϊi������߶�Ӧ����������ߵ�������
																													�����������ǰ������������ں�
																													����ǰ������Ǳ�Ե�ߣ���ֻ��Ӧһ������ߣ�ueEdgeInfo(i, 1)==-1
				const Eigen::PlainObjectBase<DerivedI>& edges
				)
	*/
	edgeUeInfo.resize(0);
	if (!getUedges(uEdges, edges))
		return false;

	const unsigned uEdgesCount = uEdges.rows();
	const unsigned edgesCount = edges.rows();

	// ȷ������������������������ӳ�䣻
	std::unordered_map<std::int64_t, int> codeMap;
	for (unsigned i = 0; i < uEdgesCount; ++i)
		codeMap.insert({ encodeUedge(uEdges(i, 0), uEdges(i, 1)), i });
	edgeUeInfo.resize(edgesCount);
	for (unsigned i = 0; i < edgesCount; ++i)
	{
		std::int64_t code = encodeUedge(edges(i, 0), edges(i, 1));
		edgeUeInfo(i) = codeMap[code];
	}

	// ȷ������������������������ӳ�䣻
	ueEdgeInfo.resize(uEdgesCount, 2);
	ueEdgeInfo.setConstant(-1);
	for (int i = 0; i < edgesCount; ++i)
	{
		int ueIdx = edgeUeInfo(i);
		if (ueEdgeInfo(ueIdx, 0) < 0)
			ueEdgeInfo(ueIdx, 0) = i;
		else
			ueEdgeInfo(ueIdx, 1) = i;
	}

	// ueEdgeInfo��ÿ��������������򡪡�������ǰ�����ں�
	for (int i = 0; i < uEdgesCount; ++i)
	{
		if (ueEdgeInfo(i, 1) < 0)
			continue;
		int eIdx0 = ueEdgeInfo(i, 0);
		if (edges(eIdx0, 0) > edges(eIdx0, 1))
		{
			int tmp = ueEdgeInfo(i, 0);
			ueEdgeInfo(i, 0) = ueEdgeInfo(i, 1);
			ueEdgeInfo(i, 1) = tmp;
		}
	}

	return true;
}


template <typename DerivedI>
bool getUeInfos(Eigen::MatrixXi& UeTrisInfo, Eigen::MatrixXi& UeCornersInfo, const Eigen::VectorXi& edgeUeInfo, \
	const Eigen::PlainObjectBase<DerivedI>& uEdges, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		bool getUeInfos(
				Eigen::MatrixXi& UeTrisInfo,						ueCount��2������߹���������Ƭ��������
																					������Ϊi�������Ϊ(A,B)����UeTrisInfo(i, 0)Ϊ�����(A,B)���ڵ�����Ƭ������
																							UeTrisInfo(i, 1)Ϊ�����(B, A)���ڵ�����Ƭ������
				Eigen::MatrixXi& UeCornersInfo,				ueCount��2����������ԵĶ����ǣ�
																					�����ǣ�0-a, 1-b, 2-c
																					������Ϊi�������Ϊ(A,B)����UeCornersInfo(i, 0)Ϊ�����(A,B)���ԵĶ���ı�ǣ�
																							��UeCornersInfo(i, 1)Ϊ�����(B,A)���ԵĶ���ı�ǡ�

				const Eigen::VectorXi& edgeUeInfo,										����������������������ӳ�䣻
				const Eigen::PlainObjectBase<DerivedI>& uEdges,				���������
				const Eigen::PlainObjectBase<DerivedI>& tris						����Ƭ����
				)

		ע�⣺�����ͳһʹ�������ʾ����vaIdx < vbIdx;
	*/
	const unsigned uEdgesCount = edgeUeInfo.maxCoeff() + 1;
	UeTrisInfo.resize(uEdgesCount, 2);
	UeCornersInfo.resize(uEdgesCount, 2);
	UeTrisInfo.setConstant(-1);
	UeCornersInfo.setConstant(-1);

	for (int i = 0; i < tris.rows(); i++)
	{
		for (int k = 0; k < 3; k++)
		{
			const int eIdx = k * tris.rows() + i;
			const int ueIdx = edgeUeInfo(eIdx);

			if (tris(i, (k + 1) % 3) == uEdges(ueIdx, 0) && tris(i, (k + 2) % 3) == uEdges(ueIdx, 1))
			{
				// ����ǰ�����ab�������AB��ͬ
				UeTrisInfo(ueIdx, 0) = i;
				UeCornersInfo(ueIdx, 0) = k;
			}
			else
			{
				// ����ǰ�����ab�������AB�෴
				UeTrisInfo(ueIdx, 1) = i;
				UeCornersInfo(ueIdx, 1) = k;
			}
		}
	}

	return true;
}


// ����������������Ƭ����
template <typename DerivedV, typename DerivedF, typename DerivedN>
bool getTriNorms(Eigen::PlainObjectBase<DerivedN>& triNorms, const Eigen::MatrixBase<DerivedV>& vers, \
	const Eigen::MatrixBase<DerivedF>& tris)
{
	triNorms.resize(tris.rows(), 3);
	int Frows = tris.rows();
	for (int i = 0; i < Frows; i++)
	{
		const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v1 = vers.row(tris(i, 1)) - vers.row(tris(i, 0));
		const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v2 = vers.row(tris(i, 2)) - vers.row(tris(i, 0));
		triNorms.row(i) = v1.cross(v2);
		typename DerivedV::Scalar r = triNorms.row(i).norm();
		if (r == 0)
		{
			std::cout << "error!!! degenerate triangle detected, unable to calculate its normal dir." << std::endl;
			return false;
		}
		else
			triNorms.row(i) /= r;
	}

	return true;
}


// �õ�����������ʾ������ߣ�
template <typename DerivedV>
void getEdgeArrows(Eigen::PlainObjectBase<DerivedV>& arrows, const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers)
{
	DerivedV vas, vbs;
	Eigen::VectorXi vaIdxes = edges.col(0);
	Eigen::VectorXi vbIdxes = edges.col(1);
	subFromIdxVec(vas, vers, vaIdxes);
	subFromIdxVec(vbs, vers, vbIdxes);
	arrows = vbs - vas;
}


// ����������������Ƭ������
template<typename T, typename DerivedV, typename DerivedI>
bool trianglesBarycenter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& barys, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	const double eps = 1e-6;
	int versCount = vers.rows();
	int trisCount = tris.rows();
	if (versCount == 0 || trisCount == 0)
		return false;

	barys.resize(trisCount, 3);
	for (int i = 0; i < trisCount; ++i)
	{
		int vaIdx, vbIdx, vcIdx;
		vaIdx = tris(i, 0);
		vbIdx = tris(i, 1);
		vcIdx = tris(i, 2);
		Eigen::Matrix<T, 3, 3> tmpMat;
		tmpMat.row(0) = vers.row(vaIdx).array().cast<T>();
		tmpMat.row(1) = vers.row(vbIdx).array().cast<T>();
		tmpMat.row(2) = vers.row(vcIdx).array().cast<T>();
		barys.row(i) = tmpMat.colwise().mean();
	}

	return true;
}


// ����������������Ƭ�Ĺ�һ��������
template<typename T, typename DerivedV, typename DerivedI>
bool trianglesNorm(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& triNorms, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	const double eps = 1e-12;
	int versCount = vers.rows();
	int trisCount = tris.rows();
	if (versCount == 0 || trisCount == 0)
		return false;

	triNorms.resize(trisCount, 3);
	for (int i = 0; i < trisCount; ++i)
	{
		int vaIdx, vbIdx, vcIdx;
		Eigen::RowVector3d va, vb, vc, arrow1, arrow2, norm;
		vaIdx = tris(i, 0);
		vbIdx = tris(i, 1);
		vcIdx = tris(i, 2);
		va = vers.row(vaIdx).array().cast<double>();
		vb = vers.row(vbIdx).array().cast<double>();
		vc = vers.row(vcIdx).array().cast<double>();
		arrow1 = vb - va;
		arrow2 = vc - va;
		triNorms.row(i) = (arrow1.cross(arrow2)).array().cast<T>();
	}

	// ��������һ�����������˻�����Ƭ����дΪinf
	for (int i = 0; i < trisCount; ++i)
	{
		double length = triNorms.row(i).norm();
		if (abs(length) < eps)
		{
			triNorms.row(i) = Eigen::Matrix<T, 1, 3>{ INFINITY, INFINITY, INFINITY };
			continue;
		}
		Eigen::Matrix<T, 1, 3> tmpVec = triNorms.row(i) / length;
		triNorms.row(i) = tmpVec;
	}

	return true;
}


// ����������������Ƭ����ƽ��
template<typename DerivedV, typename DerivedI>
bool trianglesPlane(Eigen::MatrixXd& planeCoeff, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	int versCount = vers.rows();
	int trisCount = tris.rows();
	if (versCount == 0 || trisCount == 0)
		return false;

	Eigen::MatrixXd triNorms;
	if (!trianglesNorm(triNorms, vers, tris))
		return false;

	planeCoeff.resize(trisCount, 4);
	planeCoeff.leftCols(3) = triNorms;
	for (int i = 0; i < trisCount; ++i)
	{
		Eigen::RowVector3d norm = triNorms.row(i);

		// �������˻�����Ƭ����ƽ��ϵ����дΪNAN:
		if (std::isinf(norm(0)))
			planeCoeff(i, 3) = INFINITY;

		Eigen::RowVector3d va = vers.row(tris(i, 0)).array().cast<double>();

		// va����ƽ���� �� va�㵽ƽ��ľ���Ϊ0 �� norm.dot(va) +d == 0; �� d = -norm.dot(va);
		planeCoeff(i, 3) = -norm.dot(va);
	}

	return true;
}


// ������������ÿ������Ƭ�������
template<typename Scalar, typename Ta, typename DerivedI>
bool trisArea(Eigen::Matrix<Ta, Eigen::Dynamic, 1>& trisAreaVec, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const Eigen::MatrixBase<DerivedI>& tris)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> vas, vbs, vcs, arrows1, arrows2, arrows3;
	Eigen::VectorXd lens1, lens2, lens3, s;
	Eigen::VectorXi vaIdxes = tris.col(0);
	Eigen::VectorXi vbIdxes = tris.col(1);
	Eigen::VectorXi vcIdxes = tris.col(2);

	subFromIdxVec(vas, vers, vaIdxes);
	subFromIdxVec(vbs, vers, vbIdxes);
	subFromIdxVec(vcs, vers, vcIdxes);
	arrows1 = vbs - vas;
	arrows2 = vcs - vas;
	arrows3 = vbs - vcs;
	lens1 = arrows1.rowwise().norm().cast<double>();
	lens2 = arrows2.rowwise().norm().cast<double>();
	lens3 = arrows3.rowwise().norm().cast<double>();
	s = (lens1 + lens2 + lens3) / 2.0;				// ����Ƭ�ܳ���һ�룻

	// ʹ��Heron��ʽ��������Ƭ���:
	Eigen::VectorXd tmpVec = s.array() * (s - lens1).array() * (s - lens2).array() * (s - lens3).array();
	tmpVec = tmpVec.cwiseSqrt();
	trisAreaVec = tmpVec.array().cast<Ta>();

	return true;
}


template <typename Tv, typename DerivedF, typename DerivedL>
void edge_lengths(Eigen::PlainObjectBase<DerivedL>& lenMat, const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers,
	const Eigen::MatrixBase<DerivedF>& tris) {
	const int trisCount = tris.rows();
	lenMat.resize(trisCount, 3);
	PARALLEL_FOR(0, trisCount, [&vers, &tris, &lenMat](const int i)
		{
			lenMat(i, 0) = (vers.row(tris(i, 1)) - vers.row(tris(i, 2))).norm();
			lenMat(i, 1) = (vers.row(tris(i, 2)) - vers.row(tris(i, 0))).norm();
			lenMat(i, 2) = (vers.row(tris(i, 0)) - vers.row(tris(i, 1))).norm();
		});
}


template <typename Tv, typename DerivedF, typename DerivedL>
void squared_edge_lengths(const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers,
	const Eigen::MatrixBase<DerivedF>& tris, Eigen::PlainObjectBase<DerivedL>& lenMat2)
{
	const int trisCount = tris.rows();

	lenMat2.resize(trisCount, 3);
	PARALLEL_FOR(0, trisCount, [&vers, &tris, &lenMat2](const int i)
		{
			lenMat2(i, 0) = (vers.row(tris(i, 1)) - vers.row(tris(i, 2))).squaredNorm();
			lenMat2(i, 1) = (vers.row(tris(i, 2)) - vers.row(tris(i, 0))).squaredNorm();
			lenMat2(i, 2) = (vers.row(tris(i, 0)) - vers.row(tris(i, 1))).squaredNorm();
		});
}


// ������������������
template<typename DerivedV, typename DerivedI>
double meshVolume(const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	double volume = -1.0;
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	if (vers.cols() != 3 || tris.cols() != 3)
		return volume;

	Eigen::VectorXi vaIdxes = tris.col(0);
	Eigen::VectorXi vbIdxes = tris.col(1);
	Eigen::VectorXi vcIdxes = tris.col(2);
	DerivedV versA, versB, versC;
	subFromIdxVec(versA, vers, vaIdxes);
	subFromIdxVec(versB, vers, vbIdxes);
	subFromIdxVec(versC, vers, vcIdxes);

	// һ������Ƭ��Ӧ��������ķ��������signed volume�� V(t0) == (-x3y2z1 + x2y3z1 + x3y1z2 - x1y3z2 - x2y1z3 + x1y2z3) / 6.0;
	auto volumeVec = \
		(-versC.col(0).array() * versB.col(1).array() * versA.col(2).array() + versB.col(0).array() * versC.col(1).array() * versA.col(2).array()\
			+ versC.col(0).array() * versA.col(1).array() * versB.col(2).array() - versA.col(0).array() * versC.col(1).array() * versB.col(2).array()\
			- versB.col(0).array() * versA.col(1).array() * versC.col(2).array() + versA.col(0).array() * versB.col(1).array() * versC.col(2).array()) / 6.0;

	volume = volumeVec.sum();

	return volume;
}


// ����������Ĳ�ͬȨ�ص��ڽӾ��� 
template<typename Index>
bool adjMatrix(Eigen::SparseMatrix<Index>& adjSM_eCount, \
	Eigen::SparseMatrix<Index>& adjSM_eIdx, const Eigen::MatrixXi& tris, const Eigen::MatrixXi& edges0 = Eigen::MatrixXi{})
{
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	std::vector<Eigen::Triplet<Index>> smElems, smElems_weighted;
	smElems.reserve(edgesCount);
	smElems_weighted.reserve(edgesCount);

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
	if (0 == edges0.size())
	{
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
		for (int i = 0; i < edgesCount; ++i)
		{
			smElems.push_back(Eigen::Triplet<Index>{edges(i, 0), edges(i, 1), 1});
			smElems_weighted.push_back(Eigen::Triplet<Index>{edges(i, 0), edges(i, 1), i});
		}
	}
	else
	{
		for (int i = 0; i < edgesCount; ++i)
		{
			smElems.push_back(Eigen::Triplet<Index>{edges0(i, 0), edges0(i, 1), 1});
			smElems_weighted.push_back(Eigen::Triplet<Index>{edges0(i, 0), edges0(i, 1), i});
		}
	}

	adjSM_eCount.resize(versCount, versCount);
	adjSM_eIdx.resize(versCount, versCount);
	adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());										// Ȩ��Ϊ��������ظ��Ĵ�����
	adjSM_eIdx.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());		// Ȩ��Ϊ������ߵ�������

	return true;
}


// ��Ե����ߣ�����1��
template<typename DerivedI>
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
	Eigen::SparseMatrix<int> adjSM_tmp, adjSM_ueCount;

	// 1. ȷ����Ե����ߣ�
	spMatTranspose(adjSM_tmp, adjSM_eCount);
	adjSM_ueCount = adjSM_eCount + adjSM_tmp;				// ������ڽӾ���Ȩ���Ǹñߵ��ظ�������

	std::vector<std::pair<int, int>> bdryUeVec;						// ���б�Ե����ߣ�Լ������߱�ʾΪǰ���С�������ԣ�
	traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
		{
			if (1 == iter.value() && iter.row() > iter.col())
				bdryUeVec.push_back({ iter.row(), iter.col() });
		});

	// 2. �ɱ�Ե�����ȷ����Ե����ߣ�
	std::vector<std::pair<int, int>> bdryVec;
	bdryVec.reserve(bdryUeVec.size());
	for (const auto& pair : bdryUeVec)
	{
		if (adjSM_eCount.coeff(pair.first, pair.second) > 0)
			bdryVec.push_back({ pair.first, pair.second });
		else
			bdryVec.push_back({ pair.second, pair.first });
	}
	edges2mat(bdrys, bdryVec);

	return true;
}


// ��Ե����ߣ� ����2��
template<typename DerivedI>
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, std::vector<int>& bdryTriIdxes, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		bool bdryEdges(
				Eigen::PlainObjectBase<DerivedI>& bdrys,					bdrysCount * 2����Ե�����ݣ�������ߣ�
				std::vector<int>& bdryTriIdxes,										������Ե�ߵ�����Ƭ������
				const Eigen::PlainObjectBase<DerivedI>& tris				trisCount * 3, ��������Ƭ���ݣ�
				)

	*/

	// 1. ��������ı�Ե����ߣ�
	if (!bdryEdges(bdrys, tris))
		return false;

	if (0 == bdrys.rows())
		return true;

	const unsigned trisCount = tris.rows();

	// 2. ȷ�ϱ�Ե�����α����ڵ�����Ƭ������
	Eigen::MatrixXi edgeAs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi edgeBs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi edgeCs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi vaIdxes = tris.col(0);
	Eigen::MatrixXi vbIdxes = tris.col(1);
	Eigen::MatrixXi vcIdxes = tris.col(2);

	// 3. ��������߱���-����Ƭ�ֵ䣺
	std::unordered_multimap<std::int64_t, unsigned> edgesMap;			// һ����������������ֻ����һ������Ƭ����������������Ƿ����αߣ�
	for (unsigned i = 0; i < trisCount; ++i)
	{
		std::int64_t codeA = encodeEdge(vbIdxes(i), vcIdxes(i));
		std::int64_t codeB = encodeEdge(vcIdxes(i), vaIdxes(i));
		std::int64_t codeC = encodeEdge(vaIdxes(i), vbIdxes(i));
		edgesMap.insert({ codeA, i });
		edgesMap.insert({ codeB, i });
		edgesMap.insert({ codeC, i });
	}

	// 4. ���ֵ��в������б�Ե�����ڵ�����Ƭ������
	for (unsigned i = 0; i < bdrys.rows(); ++i)
	{
		std::int64_t code = encodeEdge(bdrys(i, 0), bdrys(i, 1));
		auto iter = edgesMap.find(code);
		unsigned codeCounts = edgesMap.count(code);
		for (unsigned k = 0; k < codeCounts; ++k)
			bdryTriIdxes.push_back((iter++)->second);
	}

	return true;
}


// ��������������Ƭ�����ǵ�����ֵ��
template <typename Tv, typename Tl>
void trisCotValues(const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers,
	const Eigen::MatrixXi& tris, Eigen::Matrix<Tl, Eigen::Dynamic, Eigen::Dynamic>& cotValues)
{
	// 1/2*cotangents corresponding angles. for triangles, columns correspond to edges bc, ca, ab;
	const int trisCount = tris.rows();

	// ����ÿ������Ƭ�������������
	Eigen::Matrix<Tl, Eigen::Dynamic, 3> lenMat2;
	Eigen::Matrix<Tl, Eigen::Dynamic, 3> lenMat;
	Eigen::Matrix<Tl, Eigen::Dynamic, 1> dblA;
	squared_edge_lengths(vers, tris, lenMat2);
	lenMat = lenMat2.array().sqrt();
	trisArea(dblA, vers, tris);
	dblA.array() *= 2;
	cotValues.resize(trisCount, 3);
	for (int i = 0; i < trisCount; i++)
	{
		cotValues(i, 0) = (lenMat2(i, 1) + lenMat2(i, 2) - lenMat2(i, 0)) / dblA(i) / 4.0;
		cotValues(i, 1) = (lenMat2(i, 2) + lenMat2(i, 0) - lenMat2(i, 1)) / dblA(i) / 4.0;
		cotValues(i, 2) = (lenMat2(i, 0) + lenMat2(i, 1) - lenMat2(i, 2)) / dblA(i) / 4.0;
	}
}



// buildAdjacency()�����������������Ƭ�ڽ���Ϣ��
using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool buildAdjacency(Eigen::MatrixXi& ttAdj_mnEdge, std::vector<tVec>& ttAdj_nmnEdge, std::vector<tVec>& ttAdj_nmnOppEdge, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		bool buildAdjacency(
					const Eigen::MatrixXi& tris,									���������Ƭ����

					Eigen::MatrixXi& ttAdj_mnEdge,							����Ƭ�ķǱ�Ե����������ڽӵ�����Ƭ������
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

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// ��������ϡ�����adjSM_eCount, adjSM_eIdx��triplet���ݣ�
	std::unordered_multimap<std::int64_t, int> edgeMap;					// �߱��롪��64λ����������ʾ�ı����ݣ������˵���������


	// 1. ������ı���Ϣ���ڽӾ���
	Eigen::SparseMatrix<int> adjSM;								// �ڽӾ���
	Eigen::SparseMatrix<int> adjSM_eCount;					// �ڽӾ�������Ϊ��������ظ��Ĵ�����
	Eigen::SparseMatrix<int> adjSM_eIdx;						// �ڽӾ���Ԫ������Ӧ���α���Ϊ�ñ�����������Ӧ�����α���Ϊ�ñ����������ĺͣ�
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;				// adjSM_eIdx��ת�ã���ʾ�Աߵ���Ϣ��
	Eigen::SparseMatrix<int> adjSM_ueCount;				// Ȩ��Ϊ�����ij����������Ƭ������
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;
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
		getEdges(edges, tris);

		// 1.2 ��������������ڽӾ���adjSM, adjSM_eCount, adjSM_eCount;
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris, edges);
		adjSM = adjSM_eCount;																		// ������ڽӾ���
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		��adjSM_eIdx(i, j)��Ӧ������������ߣ���Ȩ��ֵΪ�ñߵ����������Ƿ���������ߣ���ñ߶�Ӧ������������Ȩ��Ϊ��������֮�ͣ�
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

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
		spMatTranspose(adjSM_MN_NB_opp, adjSM_MN_NB);
	}


	// 2. ȷ���Ǳ�Ե��������ߡ�����Աߵ�������
	std::vector<int> edgesIdx_MN_NB;							// �Ǳ�Ե��������ߵ�������
	std::vector<int> edgesIdx_MN_NB_opp;					// �Ǳ�Ե��������ߵĶԱߵ�������
	{
		unsigned edgesCount_MN_NB = adjSM_MN_NB.sum();
		edgesIdx_MN_NB.reserve(edgesCount_MN_NB);					// �Ǳ�Ե���������������
		edgesIdx_MN_NB_opp.reserve(edgesCount_MN_NB);			// �Ǳ�Ե��������ߵĶԱߵ�������

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB.push_back(adjSM_eIdx.coeffRef(iter.row(), iter.col()));
			});

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB_opp.push_back(adjSM_eIdx_opp.coeffRef(iter.row(), iter.col()));
			});
	}


	// 3. ���α�-����Ƭ�ڽӹ�ϵ������Ƭ�ڽӹ�ϵ��
	std::vector<int> etInfo;						// ������ - ����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
	Eigen::VectorXi etAdj_mnEdge;							// ����������Ƭ�������������α߻��Ե��д-1��
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

			if (edgesIdxOpp < 0 || edgesIdxOpp >= edgesCount)
				std::cout << "pause" << std::endl;

			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//	3.3 ������Ƭ�ڽӾ���ttAdj_mnEdge(i, :)������Ϊi������Ƭ�������ڽӵ���������Ƭ���������з����α߻��Ե����дΪ-1��
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. ��������αߣ�������������Ƭ������ߣ���Ϣ
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_map;				// ����������ߵ����� - �ñ߶�Ӧ������������
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_opp_map;		// ����������ߵ����� - �ñߵĶԱ߶�Ӧ������������
	{
		//	4.1 �����ڽӾ����ҳ����з����αߣ�������������Ƭ������ߣ���
		std::vector<int> edgeIdx_nmn;						// �����������������
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)) > 1)
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();

		int neCount = edgeIdx_nmn.size();						// ����������߸�����
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
			// ע�⣡������������������� a. ��ǰ���Ƿ����αߣ��Ա������αߣ�	b. ��ǰ���Ƿ����αߣ��Ա߲����ڣ�
			int neIdx = edgeIdx_nmn[i];		// ��ǰ�����α�������
			int vaIdx = nmnEdges(i, 0);
			int vbIdx = nmnEdges(i, 1);
			std::int64_t eCode = encodeEdge(vaIdx, vbIdx);
			std::int64_t oppEcode = encodeEdge(vbIdx, vaIdx);
			auto iter = edgeMap.find(eCode);
			auto oppIter = edgeMap.find(oppEcode);

			int otherIdx = (iter->second == neIdx) ? ((++iter)->second) : (iter->second);
			edgeIdx_nmn_map.insert({ neIdx, std::vector<int>{neIdx, otherIdx} });

			if (edgeMap.end() == oppIter)			// ��ǰ���Ƿ����αߣ��Ա߲����ڣ�
			{
				edgeIdx_nmn_opp_map.insert({ neIdx, std::vector<int>{} });
				continue;
			}
			int oppIdx1 = oppIter->second;
			int oppIdx2 = (++oppIter)->second;
			edgeIdx_nmn_opp_map.insert({ neIdx, std::vector<int>{oppIdx1, oppIdx2} });
		}
	}


	// 5. ���з����αߵ�����Ƭ���ڽӹ�ϵ��
	{
		ttAdj_nmnEdge.resize(trisCount);
		for (auto& vec : ttAdj_nmnEdge)
			vec.resize(3);

		for (const auto& pair : edgeIdx_nmn_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;

			std::vector<int> trisIdx_nmn = pair.second;
			for (auto& index : trisIdx_nmn)
				index = etInfo[index];

			ttAdj_nmnEdge[row][col] = trisIdx_nmn;
		}

		ttAdj_nmnOppEdge.resize(trisCount);
		for (auto& vec : ttAdj_nmnOppEdge)
			vec.resize(3);
		for (const auto& pair : edgeIdx_nmn_opp_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;

			std::vector<int> trisIdx_nmn_opp = pair.second;
			for (auto& index : trisIdx_nmn_opp)
				index = etInfo[index];

			ttAdj_nmnOppEdge[row][col] = trisIdx_nmn_opp;
		}
	}

	return true;
}


// buildAdjacency_new()�����������������Ƭ�ڽ���Ϣ��
using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool buildAdjacency_new(Eigen::MatrixXi& ttAdj_mnEdge, std::vector<tVec>& ttAdj_nmnEdge, std::vector<tVec>& ttAdj_nmnOppEdge, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		bool buildAdjacency_new(
					const Eigen::MatrixXi& tris,									���������Ƭ����

					Eigen::MatrixXi& ttAdj_mnEdge,							ttAdj_mnEdge(i, :)������Ϊi������Ƭ�������ڽӵ���������Ƭ
																									trisCount * 3
																									(i, 0)Ϊ����Ϊi������Ƭ�����������(vb, vc)���Ǳ�Ե���ڽӵ�����Ƭ��������
																									��(vb, vc)Ϊ����������ߣ���(i, 0)Ϊ-1��
																									��(vb, vc)Ϊ��Ե����ߣ���(i, 0)Ϊ-2��

					std::vector<tVec>& ttAdj_nmnEdge,				����Ƭ�ķ�������������ڵ�����Ƭ������
																									size == trisCount;
																									ÿ��tVecԪ�ذ�������int������
																									ttAdj_nmnEdge[i][0]��Ϊ����Ϊi������Ƭ�ķ����������(vb, vc)���ڵ���������Ƭ������


					std::vector<tVec>& ttAdj_nmnOppEdge		����Ƭ�ķ�����������ڽӵ�����Ƭ���������Ա����ڵ�����Ƭ����������
																									size == trisCount;
																									ttAdj_nmnOppEdge[i][0]��Ϊ����Ϊi������Ƭ�ķ����������(vb, vc)�ڽӵ���������Ƭ������
					)

	*/
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	Eigen::MatrixXi edges;							// ��������ݣ�
	Eigen::MatrixXi nmnEdges;					// �����������

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// ��������ϡ�����adjSM_eCount, adjSM_eIdx��triplet���ݣ�

	// 1. ������ı���Ϣ���ڽӾ��� 
	Eigen::SparseMatrix<int> adjSM_eCount;					// �ڽӾ�������Ϊ��������ظ��Ĵ�����
	Eigen::SparseMatrix<int> adjSM_eCount_opp;
	Eigen::SparseMatrix<int> adjSM_eIdx;						// �ڽӾ���Ԫ������Ӧ���α���Ϊ�ñ�����������Ӧ�����α���Ϊ�ñ����������ĺͣ�
	Eigen::SparseMatrix<int> adjSM_ueCount;				// Ȩ��Ϊ�����ij����������Ƭ������ 
	Eigen::SparseMatrix<int> adjSM_MN;						// ��������ߣ��Ǳ�Ե���ڽӾ���Ȩ��Ϊ1�� 
	std::unordered_set<int> bdryIdxes;									// ��Ե����ߵ�������
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
		getEdges(edges, tris);

		// 1.2 ��������������ڽӾ���adjSM, adjSM_eIdx, adjSM_eCount;
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris, edges);
		spMatTranspose(adjSM_eCount_opp, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + adjSM_eCount_opp;						// ������ڽӾ���Ȩ��Ϊ�����ij����������Ƭ������ 						 

		//	1.3 ������������ߣ��Ǳ�Ե�����ڽӾ���
		adjSM_MN = adjSM_ueCount;
		traverseSparseMatrix(adjSM_MN, [&adjSM_MN](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});
		adjSM_MN.prune([&adjSM_eCount](const Eigen::Index& row, const Eigen::Index& col, const int& value)->bool
			{
				if (1 == adjSM_eCount.coeffRef(row, col) && 1 == adjSM_eCount.coeffRef(col, row))
					return true;								// ����true��Ԫ�ر�������
				else
					return false;
			});															// �ǶԳƾ���

		// 1.4 ȷ����Ե����ߵ������� 
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				if (1 == iter.value())
				{
					int repCount = adjSM_eCount.coeff(iter.row(), iter.col());
					int repOppCount = adjSM_eCount.coeff(iter.col(), iter.row());
					int bdryEdgeIdx = (repCount > 0) ? adjSM_eIdx.coeff(iter.row(), iter.col()) : adjSM_eIdx.coeff(iter.col(), iter.row());
					bdryIdxes.insert(bdryEdgeIdx);
				}
			});
	}


	// 2. ȷ����������ߣ��Ǳ�Ե��������Աߵ�������
	std::vector<int> edgesIdx_MN;							// ��������ߣ��Ǳ�Ե����������
	std::vector<int> edgesIdx_MN_opp;					// ��������ߣ��Ǳ�Ե���ĶԱߵ�������
	unsigned edgesCount_MN = 0;
	{
		edgesCount_MN = adjSM_MN.sum();
		edgesIdx_MN.reserve(edgesCount_MN);					// ��������ߣ��Ǳ�Ե��������
		edgesIdx_MN_opp.reserve(edgesCount_MN);			// ��������ߣ��Ǳ�Ե���ĶԱߵ�������

		traverseSparseMatrix(adjSM_MN, [&](auto& iter)
			{
				edgesIdx_MN.push_back(adjSM_eIdx.coeffRef(iter.row(), iter.col()));
				edgesIdx_MN_opp.push_back(adjSM_eIdx.coeffRef(iter.col(), iter.row()));
			});
	}

	// 3. ���α�-����Ƭ�ڽӹ�ϵ������Ƭ�ڽӹ�ϵ��
	std::vector<int> etInfo;									// ��������� - ����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
	Eigen::VectorXi etAdj_mnEdge;						// etAdj_mnEdge(i)������Ϊi���������������Ƭ��������û�еĻ�д-1��
	{
		// 3.1 ���ɱ����� - ����Ƭ����ӳ���etInfo;
		etInfo.resize(edgesCount);							// etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

		// 3.2 ���������������Ƭ������
		etAdj_mnEdge = -Eigen::VectorXi::Ones(edgesCount);
		for (unsigned i = 0; i < edgesCount_MN; ++i)
		{
			const int& edgeIdx = edgesIdx_MN[i];
			const int& edgesIdxOpp = edgesIdx_MN_opp[i];
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		// 3.3 ��Ե�߸�ֵΪ-2��
		for (const auto& index : bdryIdxes)
			etAdj_mnEdge(index) = -2;

		//	3.3 ������Ƭ�ڽӾ���
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}

	// 4. ��������αߣ�������������Ƭ������ߣ���Ϣ
	std::unordered_map<std::int64_t, std::unordered_set<int>> edgeMap;							// �߱��롪��64λ����������ʾ�ı����ݣ������˵���������
	std::unordered_map<int, std::unordered_set<int>> edgeIdx_nmn_map;						// ����������ߵ����� - �ñ߶�Ӧ������������
	std::unordered_map<int, std::unordered_set<int>> edgeIdx_nmn_opp_map;				// ����������ߵ����� - �ñߵĶԱ߶�Ӧ������������
	{
		//	4.1 �����ڽӾ����ҳ����з���������ߡ�����ĳ������߻���Ա��ظ���������1�����������߶�����Ϊ����������ߣ�
		std::vector<int> edgeIdx_nmn;							// �����������������
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)) > 1 || adjSM_eCount.coeffRef(edges(i, 1), edges(i, 0)) > 1)
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();

		int neCount = edgeIdx_nmn.size();						// ����������߸�����
		subFromIdxVec(nmnEdges, edges, edgeIdx_nmn);

		// 4.2 �����߱���-�������Ĺ�ϣ��
		for (int i = 0; i < edges.rows(); ++i)
		{
			std::int64_t eCode = encodeEdge(edges(i, 0), edges(i, 1));
			auto iter = edgeMap.find(eCode);
			if (edgeMap.end() == iter)
				edgeMap.insert({ eCode, std::unordered_set<int>{} });
			edgeMap[eCode].insert(i);
		}

		// 4.3 �ڹ�ϣ�����������з����αߣ����������α߼���Աߵģ��ߡ�����������ӳ���ϵ��
		for (int i = 0; i < neCount; ++i)
		{
			int neIdx = edgeIdx_nmn[i];		// ��ǰ�����α�������
			int vaIdx = nmnEdges(i, 0);
			int vbIdx = nmnEdges(i, 1);
			std::int64_t eCode = encodeEdge(vaIdx, vbIdx);
			std::int64_t oppEcode = encodeEdge(vbIdx, vaIdx);
			auto oppIter = edgeMap.find(oppEcode);

			// ������ǰ����������ߵ����������� 
			edgeIdx_nmn_map.insert({ neIdx, edgeMap[eCode] });

			// ��ǰ�ߵĶԱ��п��ܲ����ڣ�
			if (edgeMap.end() != oppIter)
				edgeIdx_nmn_opp_map.insert({ neIdx, oppIter->second });
		}
	}

	// 5. ���з����αߵ�����Ƭ���ڽӹ�ϵ��
	{
		ttAdj_nmnEdge.resize(trisCount);
		for (auto& vec : ttAdj_nmnEdge)
			vec.resize(3);

		for (const auto& pair : edgeIdx_nmn_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;
			for (auto& index : pair.second)
				ttAdj_nmnEdge[row][col].push_back(etInfo[index]);
		}

		ttAdj_nmnOppEdge.resize(trisCount);
		for (auto& vec : ttAdj_nmnOppEdge)
			vec.resize(3);
		for (const auto& pair : edgeIdx_nmn_opp_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;
			for (auto& index : pair.second)
				ttAdj_nmnOppEdge[row][col].push_back(etInfo[index]);
		}
	}

	return true;
}



/////////////////////////////////////////////////////////////////////////////////////////////////// ͼ�α༭��

// ���д�Ļ�·��˳�ӿڣ�
template <typename T>
bool smoothCircuit2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& circuit, const float param = 0.1)
{
	int versCount = circuit.rows();
	Eigen::MatrixXd circDouble(circuit.rows(), circuit.cols());
	circDouble.array() = circuit.array().cast<double>();

	// 1. ��·���ĵ��ƶ���ԭ�㣬�����·̫С���ʵ��Ŵ�
	double scale = 0.0;
	Eigen::RowVector3d center = circDouble.colwise().mean();
	circDouble.rowwise() -= center;
	Eigen::VectorXd length = circDouble.rowwise().norm();
	double maxLen = length.maxCoeff();
	if (maxLen < 10.0)
	{
		scale = 10.0 / maxLen;
		circDouble.array() *= scale;
	}

	// 2. ����һ�����ݣ�		A = [A1; p * A2];		B = [B1; p * B2];
	Eigen::SparseMatrix<double> W, A, AT, A1, A2, B, B1, B2, ATA, ATB, tempMat, temp1, temp2;
	W.resize(versCount, versCount);
	std::vector<Eigen::Triplet<double>> data;
	data.reserve(2 * versCount);
	for (int i = 0; i < versCount; ++i)
	{
		if (i == 0)
		{
			data.push_back(Eigen::Triplet<double>(i, versCount - 1, 1));
			data.push_back(Eigen::Triplet<double>(i, i + 1, 1));
		}
		else if (i == versCount - 1)
		{
			data.push_back(Eigen::Triplet<double>(i, i - 1, 1));
			data.push_back(Eigen::Triplet<double>(i, 0, 1));
		}
		else
		{
			data.push_back(Eigen::Triplet<double>(i, i - 1, 1));
			data.push_back(Eigen::Triplet<double>(i, i + 1, 1));
		}
	}
	W.setFromTriplets(data.begin(), data.end());

	A1.resize(versCount, versCount);
	B1.resize(circuit.rows(), circuit.cols());
	A2.resize(versCount, versCount);
	A1.setIdentity();
	A1 -= 0.5 * W;
	B1.setZero();
	A2.setIdentity();
	B2 = circDouble.sparseView();
	tempMat.resize(versCount, 2 * versCount);

	spMatTranspose(temp1, A1);
	spMatTranspose(temp2, A2);
	tempMat.leftCols(versCount) = temp1;
	tempMat.rightCols(versCount) = param * temp2;
	spMatTranspose(A, tempMat);
	spMatTranspose(AT, A);

	tempMat.resize(3, 2 * versCount);
	spMatTranspose(temp1, B1);
	spMatTranspose(temp2, B2);
	tempMat.leftCols(versCount) = temp1;
	tempMat.rightCols(versCount) = param * temp2;
	spMatTranspose(B, tempMat);

	ATA = AT * A;
	ATB = AT * B;
	ATA.makeCompressed();
	ATB.makeCompressed();
	W.resize(0, 0);
	A.resize(0, 0);
	A1.resize(0, 0);
	A2.resize(0, 0);
	B.resize(0, 0);
	B1.resize(0, 0);
	B2.resize(0, 0);
	tempMat.resize(0, 0);
	temp1.resize(0, 0);
	temp2.resize(0, 0);

	// 3. ��С���˷��ⳬ�����Է�����A*X = B;
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(ATA);
	if (solver.info() != Eigen::Success)
	{
		std::cout << "error!!! matrix decomposition failed" << std::endl;
		return false;
	}
	circDouble = solver.solve(ATB);
	if (solver.info() != Eigen::Success)
	{
		std::cout << "error!!! solving linear equations failed" << std::endl;
		return false;
	}

	// 4. ��·����λ�Ƶ�ԭλ�������
	if (scale > 0)
		circDouble.array() /= scale;
	circDouble.rowwise() += center;
	circuit.array() = circDouble.array().cast<T>();

	return true;
}

// �����������ϲ���������������һ��������
template <typename T>
void concatMeshMat(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers1, const Eigen::MatrixXi& tris1)
{
	int versCount = vers.rows();
	matInsertRows(vers, vers1);

	Eigen::MatrixXi trisCopy1 = tris1;
	int* intPtr = trisCopy1.data();
	for (int i = 0; i < trisCopy1.size(); ++i)
		*(intPtr++) = versCount + *intPtr;

	matInsertRows(tris, trisCopy1);
};

// ȥ������Ƭ��
template <typename IndexT>
bool removeTris(Eigen::MatrixXi& trisOut, const Eigen::MatrixXi& tris, const std::vector<IndexT>& sickTriIdxes)
{
	const unsigned trisCount = tris.rows();
	if (sickTriIdxes.size() == 0)
		return false;

	Eigen::VectorXi tmpVec = Eigen::VectorXi::Ones(trisCount);
	for (const auto& index : sickTriIdxes)
	{
		if (index < IndexT(0) || index >= trisCount)
			return false;
		tmpVec(index) = -1;
	}

	std::vector<unsigned> selectedTriIdxes;
	selectedTriIdxes.reserve(trisCount);
	for (unsigned i = 0; i < trisCount; ++i)
		if (tmpVec(i) > 0)
			selectedTriIdxes.push_back(i);
	selectedTriIdxes.shrink_to_fit();

	trisOut.resize(0, 0);
	subFromIdxVec(trisOut, tris, selectedTriIdxes);

	return true;
}



/////////////////////////////////////////////////////////////////////////////////////////////////// ���������㷨��

// triangleGrow()������ָ������Ƭ��ʼ�����������������񲻿����з����α�
template <typename DerivedV, typename DerivedI>
bool triangleGrow(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx)
{
	/*

			bool  triangleGrow(
						versOut,							// �������ĵ���
						trisOut,							// ������������Ƭ
						vers,								// ��������ĵ���
						tris,								// �������������Ƭ
						triIdx							// ����Ƭ�������������Ƭ������
						)

			�������������ڷ����αߣ����return false
	*/
	int versCount = vers.rows();
	int trisCount = tris.rows();
	int edgesCount = 3 * trisCount;

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// ��������ϡ�����adjSM_eCount, adjSM_eIdx��triplet���ݣ�	

	// 1. ������ı���Ϣ���ڽӾ���
	Eigen::SparseMatrix<int> adjSM;								// �ڽӾ���
	Eigen::SparseMatrix<int> adjSM_eCount;					// �ڽӾ�������Ϊ��������ظ��Ĵ�����
	Eigen::SparseMatrix<int> adjSM_eIdx;				// �ڽӾ���Ԫ������Ӧ���α���Ϊ�ñ�����������Ӧ�����α���Ϊ�ñ����������ĺͣ�
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;		// adjSM_eIdx��ת�ã���ʾ�Աߵ���Ϣ��
	Eigen::SparseMatrix<int> adjSM_ueCount;				// Ȩ��Ϊ�����ij����������Ƭ������
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;				// �Ǳ�Ե����������ڽӾ���
	Eigen::SparseMatrix<int> adjSM_MN_NB;					// �Ǳ�Ե����������ڽӾ���
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;			// �Ǳ�Ե��������߶Աߵ��ڽӾ���
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

		// 1.1 ��������������ڽӾ���adjSM, adjSM_eCount, smElems_weighted
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
		adjSM = adjSM_eCount;																		// ������ڽӾ���
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		��adjSM_eIdx(i, j)��Ӧ������������ߣ���Ȩ��ֵΪ�ñߵ����������Ƿ���������ߣ���ñ߶�Ӧ������������Ȩ��Ϊ��������֮�ͣ�
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

		//	1.2 ����������ڽӾ���
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// ������ڽӾ���
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.3 ���ɷǱ�Ե����������ڽӾ���
		smElems.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				if (2 == iter.value())
					smElems.push_back(Eigen::Triplet<int>(iter.row(), iter.col(), 1));
			});
		adjSM_ueMN_NB.resize(versCount, versCount);
		adjSM_ueMN_NB.setFromTriplets(smElems.begin(), smElems.end());
		smElems.clear();

		//	1.4 ���ɷǱ�Ե��������ߡ�����Աߵ��ڽӾ���
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
	std::vector<int> edgesIdx_MN_NB;							// �Ǳ�Ե��������ߵ�������
	std::vector<int> edgesIdx_MN_NB_opp;					// �Ǳ�Ե��������ߵĶԱߵ�������
	{
		unsigned edgesCount_MN_NB = adjSM_MN_NB.sum();
		edgesIdx_MN_NB.reserve(edgesCount_MN_NB);					// �Ǳ�Ե���������������
		edgesIdx_MN_NB_opp.reserve(edgesCount_MN_NB);			// �Ǳ�Ե��������ߵĶԱߵ�������

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB.push_back(adjSM_eIdx.coeffRef(iter.row(), iter.col()));
			});

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB_opp.push_back(adjSM_eIdx_opp.coeffRef(iter.row(), iter.col()));
			});
	}


	// 3. ���α�-����Ƭ�ڽӹ�ϵ������Ƭ�ڽӹ�ϵ��
	std::vector<int> etInfo;												// ������ - ����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
	Eigen::VectorXi etAdj_mnEdge;							// ����������Ƭ�������������α߻��Ե��д-1��
	Eigen::MatrixXi ttAdj_mnEdge;					// ����Ƭ�ķǱ�Ե����������ڽӵ�����Ƭ������	trisCount * 3(i, 0)Ϊ����Ϊi������Ƭ�ķǱ�Ե���������(vb, vc)�ڽӵ�����Ƭ��������
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
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. ѭ�����������е�һ������Ƭt���� - tд��������� - ���t������������Ƭ������Ȼ�����
	std::unordered_set<int> finalTriIdx;						// ��ͨ����Ƭ���ϣ�
	std::unordered_set<int> workingSet;						// ����Ѱ����������Ƭ�Ķ��У�
	workingSet.insert(static_cast<int>(triIdx));							// ���в����һ������Ƭ��
	while (workingSet.size() > 0)
	{
		// 4.1 �׸�Ԫ�س��ӣ����뵽��ͨ����Ƭ�����У�
		int cTriIdx = *workingSet.begin();									// ��ǰ���������Ƭ����
		std::vector<int> adjTriIdx;												// ��ǰ����Ƭ��������Ƭ������
		workingSet.erase(workingSet.begin());
		finalTriIdx.insert(cTriIdx);

		// 4.2 ȷ����ǰ����Ƭ����������Ƭ
		for (int i = 0; i < 3; ++i)
			if (ttAdj_mnEdge(cTriIdx, i) >= 0)
				adjTriIdx.push_back(ttAdj_mnEdge(cTriIdx, i));

		// 4.3 ����ͨ����Ƭ������û�е�ǰ������Ƭ�����������Ƭ��������Ҳ���������Ƭ��
		for (const auto& index : adjTriIdx)
		{
			auto retPair = finalTriIdx.insert(index);
			if (retPair.second)
				workingSet.insert(index);
		}
	}


	// 5. ��ѡ��������Ƭȡ���㣺
	std::set<int> finalVerIdx;
	Eigen::VectorXi finalVerIdxVec;													// ����������������ԭ�����е�������
	Eigen::VectorXi oldNewIdxInfo = -Eigen::VectorXi::Ones(vers.rows());		// �������������±�Ϊ��������Ԫ��ֵΪ�������������������еĶ���дΪ-1
	for (const auto& index : finalTriIdx)
	{
		finalVerIdx.insert(tris(static_cast<int>(index), 0));
		finalVerIdx.insert(tris(static_cast<int>(index), 1));
		finalVerIdx.insert(tris(static_cast<int>(index), 2));
	}

	auto iter = finalVerIdx.begin();
	finalVerIdxVec.resize(finalVerIdx.size());
	for (int i = 0; i < finalVerIdx.size(); ++i)
		finalVerIdxVec(i) = static_cast<int>(*iter++);

	int setOrder = 0;
	for (const auto& index : finalVerIdx)
		oldNewIdxInfo[static_cast<int>(index)] = setOrder++;
	subFromIdxVec(versOut, vers, finalVerIdxVec);


	// 6. ԭ������������Ƭ�еĶ����������ϵĻ�Ϊ�µġ�
	Eigen::MatrixXi tempMat = tris;
	int* intPtr = tempMat.data();
	for (int i = 0; i < tempMat.size(); ++i)
	{
		int oldIdx = *intPtr;
		*intPtr = oldNewIdxInfo[oldIdx];
		intPtr++;
	}

	// 7. ȥ���Ƿ�����Ƭ��������-1Ԫ�ص�����Ƭ��
	trisOut.resize(tris.rows(), 3);
	int count = 0;
	for (int i = 0; i < tris.rows(); ++i)
	{
		Eigen::RowVector3i currentTri = tempMat.row(i);
		if ((currentTri.array() >= 0).all())
			trisOut.row(count++) = currentTri;
	}
	trisOut.conservativeResize(count, 3);

	return true;
}


// triangleGrowSplitMesh()�����������������������Ϊ�������ͨ�����������񲻿����з����α�
template <typename DerivedV>
bool triangleGrowSplitMesh(std::vector<DerivedV>& meshesVersOut, std::vector<Eigen::MatrixXi>& meshesTrisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	int versCount = vers.rows();
	int trisCount = tris.rows();
	int edgesCount = 3 * trisCount;

	// ��������е�triIdxΪ������ߵ������Ƭ������
	std::unordered_set<int> finalTriIdx;						// ��ͨ����Ƭ���ϣ�
	std::unordered_set<int> workingSet;						// ����Ѱ����������Ƭ�Ķ��У�

	Eigen::MatrixXi ttAdj_mnEdge;		// ����Ƭ�ķǱ�Ե����������ڽӵ�����Ƭ������	trisCount * 3(i, 0)Ϊ����Ϊi������Ƭ�ķǱ�Ե���������(vb, vc)�ڽӵ�����Ƭ��������

	Eigen::SparseMatrix<int> adjSM;								// �ڽӾ���
	Eigen::SparseMatrix<int> adjSM_eCount;					// �ڽӾ�������Ϊ��������ظ��Ĵ�����
	Eigen::SparseMatrix<int> adjSM_eIdx;				// �ڽӾ���Ԫ������Ӧ���α���Ϊ�ñ�����������Ӧ�����α���Ϊ�ñ����������ĺͣ�
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;		// adjSM_eIdx��ת�ã���ʾ�Աߵ���Ϣ��
	Eigen::SparseMatrix<int> adjSM_ueCount;				// Ȩ��Ϊ�����ij����������Ƭ������
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;				// �Ǳ�Ե����������ڽӾ���
	Eigen::SparseMatrix<int> adjSM_MN_NB;					// �Ǳ�Ե����������ڽӾ���
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;			// �Ǳ�Ե��������߶Աߵ��ڽӾ���

	Eigen::VectorXi etAdj_mnEdge;							// ����������Ƭ�������������α߻��Ե��д-1��

	std::vector<int> etInfo;												// ������ - ����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
	std::vector<int> edgesIdx_MN_NB;							// �Ǳ�Ե��������ߵ�������
	std::vector<int> edgesIdx_MN_NB_opp;					// �Ǳ�Ե��������ߵĶԱߵ�������
	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// ��������ϡ�����adjSM_eCount, adjSM_eIdx��triplet���ݣ�	

	std::vector<bool> trisFlag(trisCount, true);		// true��ʾ������Ƭ���ռ���false��ʾ���ռ���


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

		// 1.1 ��������������ڽӾ���adjSM, adjSM_eCount, smElems_weighted
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
		adjSM = adjSM_eCount;																		// ������ڽӾ���
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		��adjSM_eIdx(i, j)��Ӧ������������ߣ���Ȩ��ֵΪ�ñߵ����������Ƿ���������ߣ���ñ߶�Ӧ������������Ȩ��Ϊ��������֮�ͣ�
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

		//	1.2 ����������ڽӾ���
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// ������ڽӾ���
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.3 ���ɷǱ�Ե����������ڽӾ���
		smElems.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				if (2 == iter.value())
					smElems.push_back(Eigen::Triplet<int>(iter.row(), iter.col(), 1));
			});
		adjSM_ueMN_NB.resize(versCount, versCount);
		adjSM_ueMN_NB.setFromTriplets(smElems.begin(), smElems.end());
		smElems.clear();

		//	1.4 ���ɷǱ�Ե��������ߡ�����Աߵ��ڽӾ���
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
				edgesIdx_MN_NB.push_back(adjSM_eIdx.coeffRef(iter.row(), iter.col()));
			});

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB_opp.push_back(adjSM_eIdx_opp.coeffRef(iter.row(), iter.col()));
			});
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
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. ѭ�����������е�һ������Ƭt���� - tд��������� - ���t������������Ƭ������Ȼ�����
	std::vector<std::unordered_set<int>> connectedTriIdxes;
	int collectedTrisCount = 0;
	while (collectedTrisCount < trisCount)
	{
		finalTriIdx.clear();

		// ȡ��ǰ�׸�δ�ռ�������Ƭ��Ϊ��������Ƭ��
		int seedTriIdx = 0;
		for (int i = 0; i < trisCount; ++i)
			if (trisFlag[i])
				seedTriIdx = i;

		workingSet.insert(static_cast<int>(seedTriIdx));					// ���в�����������Ƭ��
		while (workingSet.size() > 0)
		{
			// 4.1 �׸�Ԫ�س��ӣ����뵽��ͨ����Ƭ�����У�
			int cTriIdx = *workingSet.begin();									// ��ǰ���������Ƭ����
			trisFlag[cTriIdx] = false;

			std::vector<int> adjTriIdx;												// ��ǰ����Ƭ��������Ƭ������
			workingSet.erase(workingSet.begin());
			finalTriIdx.insert(cTriIdx);

			// 4.2 ȷ����ǰ����Ƭ����������Ƭ
			for (int i = 0; i < 3; ++i)
				if (ttAdj_mnEdge(cTriIdx, i) >= 0)
					adjTriIdx.push_back(ttAdj_mnEdge(cTriIdx, i));

			// 4.3 ����ͨ����Ƭ������û�е�ǰ������Ƭ�����������Ƭ��������Ҳ���������Ƭ��
			for (const auto& index : adjTriIdx)
			{
				auto retPair = finalTriIdx.insert(index);
				if (retPair.second && trisFlag[index])
				{
					workingSet.insert(index);
					trisFlag[index] = false;
				}
			}
		}

		connectedTriIdxes.push_back(finalTriIdx);
		collectedTrisCount += finalTriIdx.size();
	}


	// 5. ������е���ͨ����
	int meshesCount = connectedTriIdxes.size();
	meshesVersOut.resize(meshesCount);
	meshesTrisOut.resize(meshesCount);
	for (int i = 0; i < meshesCount; ++i)
	{
		DerivedV& cMeshVers = meshesVersOut[i];
		Eigen::MatrixXi& cMeshTris = meshesTrisOut[i];

		// 5.1 ��connectedTriIdxes����ͨ������Ƭ����ȡ���㣺
		std::set<int> cMeshVerIdx;																	// ��ǰ�ĵ���ͨ�����а����Ķ����������
		for (const auto& index : connectedTriIdxes[i])
		{
			cMeshVerIdx.insert(tris(static_cast<int>(index), 0));
			cMeshVerIdx.insert(tris(static_cast<int>(index), 1));
			cMeshVerIdx.insert(tris(static_cast<int>(index), 2));
		}

		Eigen::VectorXi cMeshVerIdxVec(cMeshVerIdx.size());									// ����������������ԭ�����е�������
		Eigen::VectorXi oldNewIdxInfo = -Eigen::VectorXi::Ones(vers.rows());			// �������������±�Ϊ��������Ԫ��ֵΪ�������������������еĶ���дΪ-1
		auto iter = cMeshVerIdx.begin();
		for (int i = 0; i < cMeshVerIdx.size(); ++i)
			cMeshVerIdxVec(i) = static_cast<int>(*iter++);

		int setOrder = 0;
		for (const auto& index : cMeshVerIdx)
			oldNewIdxInfo[static_cast<int>(index)] = setOrder++;
		subFromIdxVec(cMeshVers, vers, cMeshVerIdxVec);

		// 5.2 ȷ����ǰ�ĵ���ͨ���������Ƭ��
		std::vector<int> cTriIdxes;
		cTriIdxes.insert(cTriIdxes.end(), connectedTriIdxes[i].begin(), connectedTriIdxes[i].end());
		subFromIdxVec(cMeshTris, tris, cTriIdxes);
		int* verIdxPtr = cMeshTris.data();
		for (int k = 0; k < 3 * cMeshTris.rows(); ++k)				// ����Ƭ�е��ϵĶ���������Ϊ�µģ� 
		{
			int oldIdx = *verIdxPtr;
			*verIdxPtr = oldNewIdxInfo[oldIdx];
			verIdxPtr++;
		}
	}

	return true;
}


// triangleGrowOuterSurf()������ָ������Ƭ���������ⲿ������Ƭ����ʼ������������ȡ�����浥��ͨ��������
template <typename T>
bool triangleGrowOuterSurf(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, const bool blRemIso = true, const int startIdx = -1)
{
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	versOut.resize(0, 0);
	trisOut.resize(0, 0);

	int versCount = vers.rows();
	int trisCount = tris.rows();
	int edgesCount = 3 * trisCount;
	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// ��������ϡ�����adjSM_eCount, adjSM_eIdx��triplet���ݣ�	

	// 0. ������������Ƭ�����������˻�����Ƭ�����Ϊ0������ȥ����
	Eigen::MatrixXi trisClean;
	Eigen::MatrixXd tmpNorms(tris.rows(), 3);
	Eigen::MatrixXd triNorms;
	Eigen::VectorXi flagValidTri(Eigen::VectorXi::Ones(trisCount));
	bool blDegTris = false;
	for (int i = 0; i < trisCount; i++)
	{
		RowVector3T arrow1 = vers.row(tris(i, 1)) - vers.row(tris(i, 0));
		RowVector3T arrow2 = vers.row(tris(i, 2)) - vers.row(tris(i, 0));
		tmpNorms.row(i) = (arrow1.cross(arrow2)).array().cast<double>();
		double areaDB = tmpNorms.row(i).norm();
		if (0 == areaDB)								// �������˻�����Ƭ���ᵼ�¼�����Ƭ�������
		{
			blDegTris = true;
			flagValidTri(i) = 0;
		}
		else
			tmpNorms.row(i) /= areaDB;
	}
	if (blDegTris)
	{
		subFromFlagVec(trisClean, tris, flagValidTri);
		subFromFlagVec(triNorms, tmpNorms, flagValidTri);
		trisCount = trisClean.rows();
		edgesCount = 3 * trisCount;
	}
	else
	{
		trisClean = tris;
		triNorms = tmpNorms;
	}

	// 0+.�����쳣�ķ����αߡ�����Щ������meshArrange�����Խ������ɹ�����������Ƭ������ߣ�����һ������߹���һ������Ա߹���������
	Eigen::MatrixXi nmnEdges, edges, trisTmp;
	std::vector<std::pair<std::pair<int, int>, int>> nmnInfos;
	std::unordered_map<std::int64_t, int> nmnMap;						// �����������code-���������ֵ䣻
	std::unordered_set<std::int64_t> sickCodes;							// ��ҪԤ����ķ���������ߵ�code��
	std::vector<std::int64_t> edgeCodes;
	std::unordered_map<std::int64_t, int> nmnTriIdxMap;
	std::unordered_map<std::int64_t, double> nmnTriAreaMap;
	std::vector<int> sickTriIdxes;
	getEdges(edges, tris);
	edgeCodes.resize(edgesCount);
	int nmnCount = nonManifoldEdges(nmnEdges, nmnInfos, trisClean);
	for (const auto& pair : nmnInfos)
	{
		std::int64_t code = encodeEdge(pair.first.first, pair.first.second);
		nmnMap.insert({ code, pair.second });
	}
	for (const auto& pair : nmnMap)
	{
		std::int64_t code = pair.first;
		auto retPair = decodeEdge(code);
		int repCount = pair.second;
		std::int64_t codeOpp = encodeEdge(retPair.second, retPair.first);
		int repCountOpp = nmnMap[codeOpp];
		if (2 == repCount && 1 == repCountOpp)
			sickCodes.insert(code);
		if (1 == repCount && 2 == repCountOpp)
			sickCodes.insert(codeOpp);
	}
	if (!sickCodes.empty())
	{
		for (int i = 0; i < edgesCount; ++i)
			edgeCodes[i] = encodeEdge(edges(i, 0), edges(i, 1));
		for (int i = 0; i < edgesCount; ++i)
		{
			std::int64_t& code = edgeCodes[i];
			auto iter = sickCodes.find(code);
			if (sickCodes.end() != iter)
			{
				int triIdx = i % trisCount;							// ��ǰ���α����ڵ�����Ƭ������
				RowVector3T arrow1 = vers.row(trisClean(triIdx, 1)) - vers.row(trisClean(triIdx, 0));
				RowVector3T arrow2 = vers.row(trisClean(triIdx, 2)) - vers.row(trisClean(triIdx, 0));
				double area = std::abs(0.5 * arrow1.dot(arrow2));
				auto areaIter = nmnTriAreaMap.find(code);
				if (nmnTriAreaMap.end() == areaIter)
				{
					nmnTriAreaMap.insert({ code, area });
					nmnTriIdxMap.insert({ code, triIdx });
				}
				else
				{
					double areaOld = areaIter->second;
					if (area < areaOld)
					{
						areaIter->second = area;
						nmnTriIdxMap[code] = triIdx;
					}
				}
			}
		}
		for (const auto& pair : nmnTriIdxMap)
			sickTriIdxes.push_back(pair.second);

		// ����������߹�������������Ƭ��ɾ�������С���Ǹ���
		if (!removeTris(trisTmp, trisClean, sickTriIdxes))
		{
			std::cout << "error!!! removeTris() failed." << std::endl;
			return false;
		}
		trisClean = trisTmp;
	}

	// 1. ������Ƭ�ڽӹ�ϵ��
	Eigen::MatrixXi ttAdj_mnEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency_new(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, trisClean))
	{
		std::cout << "error!!! buildAdjacency() failed." << std::endl;
		return false;
	}

	// 2. ��δָ����������Ƭ��������ѡȡһ�������������αߵ�����Ƭ��Ϊ��ɢ����������Ƭ��
	int seedIdx = -1;
	if (startIdx < 0)
	{
		for (int i = 0; i < trisCount; ++i)
		{
			if (ttAdj_mnEdge(i, 0) >= 0 && ttAdj_mnEdge(i, 1) >= 0 && ttAdj_mnEdge(i, 2) >= 0)
			{
				seedIdx = i;
				break;
			}
		}
		if (seedIdx < 0)
			return false;
	}
	else
		seedIdx = startIdx;

	// 3. ѭ�����������е�һ������Ƭt���� - tд��������� - ���t������������Ƭ������Ȼ����� 
	std::unordered_set<int> finalTriIdxes;
	std::unordered_set<int> workingSet;								// ����Ѱ����������Ƭ�Ķ��У�
	std::vector<int> triTags(trisCount, 0);								// 0 : δ���ʣ�1: ����¼�� -1: ��ɾ����

	// ������ע��tagд1�Ժ����п��ܱ���дΪ-1�ģ�Ŀǰ���ǵ�һ�η��ʵ�ʱ�ı�ǩ�������Ƿ񱣴�����Ƭ��
	workingSet.insert(static_cast<int>(seedIdx));						// ���в����һ������Ƭ��
	triTags[seedIdx] = 1;
	while (workingSet.size() > 0)
	{
		// w.1 �׸�Ԫ�س��ӣ����뵽��ͨ����Ƭ�����У�
		int cTriIdx = *workingSet.begin();									// ��ǰ���������Ƭ����
		workingSet.erase(workingSet.begin());
		if (triTags[cTriIdx] < 0)											// �п������ʱ��δ����״̬�����֮�󱻱��Ϊɾ������Ƭ��
			continue;
		finalTriIdxes.insert(cTriIdx);

		// w.2 ȷ����ǰ����Ƭ���ڵ���������Ƭ���������������ɾ��
		std::vector<int> adjTriIdx;												// ��ǰ����Ƭ��������Ƭ������
		for (int i = 0; i < 3; ++i)
		{
			int nmAdjTriIdx = ttAdj_mnEdge(cTriIdx, i);
			if (nmAdjTriIdx >= 0)										// ���αߣ��Ǳ�Ե����
				adjTriIdx.push_back(nmAdjTriIdx);						// wf.a. ���α߹�������������Ƭ��������
			else if (nmAdjTriIdx == -1)								// �����α�
			{
				// wf.b.1. ɾ����ǰ�����α߹�������������Ƭ��
				auto& vec1 = ttAdj_nmnEdge[cTriIdx][i];				// ��ǰ�����α����ڵ�����Ƭ��������ǰ����Ƭ����
				for (const auto& index : vec1)
				{
					if (index == cTriIdx)
						continue;
					else
						triTags[index] = -1;
				}

				// wf.b.2. ��ǰ�����αߵĶԱ����ڵ�����Ƭ�У�ѡȡ����н����(��dotֵ��С)���Ǹ���
				auto& vec2 = ttAdj_nmnOppEdge[cTriIdx][i];			// ��ǰ�����αߵĶԱ����ڵ���������Ƭ
				if (vec2.empty())
					continue;
				else if (1 == vec2.size())
					adjTriIdx.push_back(vec2[0]);
				else
				{
					unsigned nopTrisCount = vec2.size();
					Eigen::RowVector3d cTriNorm = triNorms.row(cTriIdx);
					Eigen::MatrixXd nopTriNorms;
					Eigen::VectorXd dotValues(nopTrisCount);
					subFromIdxVec(nopTriNorms, triNorms, vec2);
					for (unsigned k = 0; k < nopTrisCount; ++k)
						dotValues(k) = cTriNorm.dot(nopTriNorms.row(k));
					Eigen::Index minIdx;
					dotValues.minCoeff(&minIdx);
					int selTriIdx = vec2[minIdx];
					adjTriIdx.push_back(selTriIdx);

					for (const auto& index : vec2)				// δ��ѡȡ�Ķ�����������Ƭȫ��ɾ�� �� 
						if (selTriIdx != index)
							triTags[index] = -1;
				}

			}

			// nmAdjTriIdx == -2��Ե�ߣ�ֱ��������
		}

		// w.3 finalTriIdxes�в��뱣����������������Ƭ����������Ƭ֮ǰδ���ʣ��������У�
		for (const auto& index : adjTriIdx)
		{
			if (0 == triTags[index])
			{
				auto retPair = finalTriIdxes.insert(index);
				if (retPair.second)
				{
					workingSet.insert(index);
					triTags[index] = 1;
				}
			}
		}
	}

	// 4. ��ȡѡ�е�����Ƭ
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vers1 = vers;
	Eigen::MatrixXi tris1;
	subFromIdxCon(tris1, trisClean, finalTriIdxes);

	// 5. ȥ���������㣻
	if (blRemIso)
	{
		std::vector<unsigned> isoVerIdxes = checkIsoVers(vers1, tris1);
		if (isoVerIdxes.empty())
		{
			versOut = vers1;
			trisOut = tris1;
		}
		else
			removeIsoVers(versOut, trisOut, vers1, tris1, isoVerIdxes);
	}
	else
	{
		versOut = vers1;
		trisOut = tris1;
	}

	return true;
}


// robustTriangleGrowOuterSurf()����������³���Ե�����Ƭ������������ĳЩ��������Ƭ�ĸ��ţ�
template <typename T>
bool robustTriangleGrowOuterSurf(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, const bool blRemIso = true)
{
	if (0 == vers.rows() || 0 == tris.rows())
		return false;

	const int trisCount = tris.rows();
	const int maxIter = 5;
	int startIdx = 0;
	bool subFlag = false;
	std::default_random_engine e;
	std::uniform_int_distribution<int>	UID_i(0, trisCount - 1);

	// ��ʱ�з�������Ƭ�ĸ��ţ�������������������ֻ��ԭ�����С��һ���֣���Ҫ���Զ�α������ִ���
	for (int i = 0; i < maxIter; ++i)
	{
		versOut.resize(0, 0);
		trisOut.resize(0, 0);
		startIdx = UID_i(e);
		triangleGrowOuterSurf(versOut, trisOut, vers, tris, blRemIso, startIdx);
		int trisCountNew = trisOut.rows();
		if (trisCountNew > trisCount / 2)						// �������Ƭ������˵���д���������
		{
			subFlag = true;
			break;
		}
	}

	return subFlag;
}



// ����ߵ�������������ȡ��·�ߺͱ����ߣ�����ı߾�������������Ҳ����������ߣ�
template <typename DerivedI>
bool uEdgeGrow(std::vector<Eigen::MatrixXi>& circs, std::vector<Eigen::MatrixXi >& segs, \
	const Eigen::PlainObjectBase<DerivedI>& uEdges)
{
	/*
	bool uEdgeGrow(
			std::vector<Eigen::MatrixXi>& circs,									Ԫ�������һ����·�ı�Ե������
			std::vector<Eigen::MatrixXi >& segs,									Ԫ�������һ���Ƿ�����ߵı�Ե������
			const Eigen::PlainObjectBase<DerivedI>& uEdges			ueCount * 2, ���������
			)

	*/
	std::list<std::list<std::pair<int, int>>> circList, segList;
	Eigen::SparseMatrix<int> adjSM;										// �ڽӾ���
	std::vector<Eigen::Triplet<int>> smElems;

	// 1. ����������ڽӾ���
	int maxIdx = 0;
	const int* idxPtr = uEdges.data();
	for (int i = 0; i < uEdges.size(); i++, idxPtr++)
		if (*idxPtr > maxIdx)
			maxIdx = *idxPtr;
	smElems.reserve(2 * uEdges.rows());
	for (unsigned i = 0; i < uEdges.rows(); ++i)
	{
		smElems.push_back(Eigen::Triplet<int>{uEdges(i, 0), uEdges(i, 1), 1});
		smElems.push_back(Eigen::Triplet<int>{uEdges(i, 1), uEdges(i, 0), 1});
	}
	adjSM.resize(maxIdx + 1, maxIdx + 1);
	adjSM.setFromTriplets(smElems.begin(), smElems.end());
	smElems.clear();
	traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
		{
			iter.valueRef() = 1;
		});

	// 2. ͳ�Ʋ��ظ��������������
	std::unordered_set<std::int64_t> tmpSet;
	for (unsigned i = 0; i < uEdges.rows(); ++i)
		tmpSet.insert(encodeUedge(uEdges(i, 0), uEdges(i, 1)));
	int remainCount = tmpSet.size();			// (3,4) (4,3)��ʾͬһ������ߣ�

	// 3. �����������ѭ����
	while (remainCount > 0)
	{
		std::list<std::pair<int, int>> currentSeg;

		// w1. ��ȡ�ڽӱ��е�һ������ߣ�
		bool breakFlag = false;
		for (unsigned i = 0; !breakFlag & i < adjSM.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, i); !breakFlag & iter; ++iter)
			{
				if (iter.valueRef() > 0)
				{
					currentSeg.push_back({ iter.row(), iter.col() });
					iter.valueRef() = -1;												// ��ʾɾ����Ԫ�أ�
					adjSM.coeffRef(iter.col(), iter.row()) = -1;
					remainCount--;
					breakFlag = true;
				}
			}
		}

		// w2. ������һ������߹����ıߣ�ֱ�����������������γɻ�·Ϊֹ��
		int head = currentSeg.front().first;
		int tail = currentSeg.front().second;
		while (1)
		{
			int otherEnd = -1;
			breakFlag = false;
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, head); !breakFlag & iter; ++iter)	// �Ե�head�еı�����
			{
				if (iter.valueRef() > 0)
				{
					otherEnd = iter.row();
					currentSeg.push_back({ head, otherEnd });
					iter.valueRef() = -1;												// ��ʾɾ����Ԫ�أ�
					adjSM.coeffRef(iter.col(), iter.row()) = -1;
					remainCount--;
					head = otherEnd;
					breakFlag = true;
				}
			}

			if (otherEnd == tail)
			{
				circList.push_back(currentSeg);
				break;
			}

			if (otherEnd < 0)
			{
				segList.push_back(currentSeg);
				break;
			}
		}
	}

	// 4. �����
	circs.reserve(circList.size());
	segs.reserve(segList.size());
	for (const auto& list : circList)
	{
		unsigned pairCount = list.size();
		Eigen::MatrixXi circEdges(pairCount, 2);
		auto iter = list.begin();
		for (unsigned i = 0; i < pairCount; ++i, iter++)
		{
			circEdges(i, 0) = iter->first;
			circEdges(i, 1) = iter->second;
		}
		circs.push_back(circEdges);
	}

	for (const auto& list : segList)
	{
		unsigned pairCount = list.size();
		Eigen::MatrixXi segEdges(pairCount, 2);
		auto iter = list.begin();
		for (unsigned i = 0; i < pairCount; ++i, iter++)
		{
			segEdges(i, 0) = iter->first;
			segEdges(i, 1) = iter->second;
		}
		segs.push_back(segEdges);
	}

	return true;
}

// ȷ�����������е����е���ͨ���򣨶�����ͨ������ʹ��ϡ�����
template <typename Index>
int simplyConnectedRegion(const Eigen::SparseMatrix<Index>& adjSM,
	Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount)
{
	/*
		int simplyConnectedRegion(												���ص���ͨ���������������ʱ����-1
					const Eigen::SparseMatrix<Index>& adjSM,			����������ڽӾ���
					Eigen::VectorXi& connectedLabels,						ͬһ��ͨ�����ڵĶ����ǩ������Ϊ0, 1, 2,...��
																											����Ϊi�Ķ���ı�ǩΪconnectedLabels(i);
					Eigen::VectorXi& connectedCount							ÿ����ǩ��Ӧ�ĵ���ͨ�����ڰ����Ķ�����
																											��ǩΪi�ĵ���ͨ��������Ķ�����ΪconnectedCount(i)
					)
	*/
	if (adjSM.rows() == 0)
		return -1;

	if (adjSM.rows() != adjSM.cols())
		return -1;

	const unsigned versCount = adjSM.rows();

	// 1. ��ʼ����
	connectedLabels.setConstant(versCount, versCount);		// ���ս���������ܵı�ǩΪversCount-1;
	connectedCount.setZero(versCount);

	// 2. �������ж��㣬����������ͨ�Ķ��㣺
	int currentLabel = 0;
	for (unsigned i = 0; i < versCount; ++i)
	{
		// 2.1 ����ǰ�����ѱ����ʹ�(label < versCount)����continue:
		if (connectedLabels(i) < versCount)
			continue;

		// 2.2 ����ǰ����δ�����ʹ���ִ��BFS�ռ���������ͨ�Ķ��㣺
		std::queue<Index> workingQueue;
		workingQueue.push(i);
		while (!workingQueue.empty())
		{
			const Index curVerIdx = workingQueue.front();
			workingQueue.pop();

			if (connectedLabels(curVerIdx) < versCount)
				continue;

			// 2.2.1 ���׶���label��ֵ��label����+1
			connectedLabels(curVerIdx) = currentLabel;
			connectedCount(currentLabel)++;

			// 2.2.2 ���ڽӾ�����������ǰ�����ڽӵĶ��㣺
			for (typename Eigen::SparseMatrix<Index>::InnerIterator iter(adjSM, curVerIdx); iter; ++iter)
			{
				const Index connectVerIdx = iter.row();					// Ĭ��ϡ����������ȴ洢����ǰ�������ڵ�curVerIdx�У�

				if (connectedLabels(connectVerIdx) < versCount)
					continue;

				workingQueue.push(connectVerIdx);
			}
		}

		// 2.3 ��һ����ǩ�Ķ����ռ���ϣ���һ��ѭ���ռ���һ����ǩ�Ķ��㣺
		currentLabel++;
	}

	// 3. shrink_to_fit()
	connectedCount.conservativeResize(currentLabel, 1);

	return currentLabel;
}


// ȷ�����������е����е���ͨ���򣨶�����ͨ������ʹ��std::unordered_set����ʹ��ϡ�����
template <typename T>
int simplyConnectedRegion(std::vector<int>& connectedLabels, std::vector<int>& connectedCount, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
{
	auto collectVers = [](std::vector<int>& connectedLabels, std::vector<bool>& triCollected, \
		std::unordered_set<int>& versSet, const Eigen::MatrixXi& tris, const int label)->int
	{
		const int trisCount = tris.rows();
		int collectedCount = 0;
		bool blCollectedNew = false;
		for (const auto& flag : triCollected)
		{
			if (flag)
				collectedCount++;
		}

		do
		{
			blCollectedNew = false;
			for (int i = 0; i < trisCount; ++i)
			{
				if (triCollected[i])
					continue;

				int a = tris(i, 0);
				int b = tris(i, 1);
				int c = tris(i, 2);
				auto iter1 = versSet.find(a);
				auto iter2 = versSet.find(b);
				auto iter3 = versSet.find(c);
				if (versSet.end() != iter1 || versSet.end() != iter2 || versSet.end() != iter3)
				{
					versSet.insert(a);
					versSet.insert(b);
					versSet.insert(c);
					connectedLabels[a] = label;
					connectedLabels[b] = label;
					connectedLabels[c] = label;
					triCollected[i] = true;
					collectedCount++;
					blCollectedNew = true;
				}
			}
		} while (blCollectedNew);						// ��������Ԫ�ر��ռ��������ִ��һ��

		return versSet.size();										// �����ռ��Ķ���������������Ϊlabel�Ķ��������
	};

	if (vers.rows() == 0 || tris.rows() == 0)
		return -1;

	const int versCount = vers.rows();
	const int trisCount = tris.rows();

	// 1. ����������ѭ����
	connectedLabels.resize(versCount, versCount);		// ���ս���������ܵı�ǩΪversCount-1;
	std::unordered_set<int> versSet;
	std::vector<bool> triCollected(trisCount, false);
	int currentLabel = 0;
	int collectedVersCount = 0;
	while (collectedVersCount < versCount)
	{
		// w1. ѡȡ��һ��δ����ǵ�����Ƭ����������������Ϊ�������������Ӷ��㣻
		for (int i = 0; i < trisCount; ++i)
		{
			if (!triCollected[i])
			{
				int a = tris(i, 0);
				int b = tris(i, 1);
				int c = tris(i, 2);
				triCollected[i] = true;
				versSet.insert(a);
				versSet.insert(b);
				versSet.insert(c);
				connectedLabels[a] = currentLabel;
				connectedLabels[b] = currentLabel;
				connectedLabels[c] = currentLabel;
				break;
			}
		}

		// w2. ����������ѭ����
		int currentLabelVersCount = collectVers(connectedLabels, triCollected, versSet, tris, currentLabel);

		// w3. post procedure:
		versSet.clear();
		collectedVersCount += currentLabelVersCount;
		currentLabel++;
	}

	// 2. ͳ�ƣ�
	int scCount = currentLabel;
	connectedCount.resize(scCount, 0);
	for (int i = 0; i < versCount; ++i)
	{
		int label = connectedLabels[i];
		connectedCount[label]++;
	}

	return scCount;
}



// ȷ�����������е����е���ͨ��������Ƭ��ͨ��������triangleGrow��ǿ������������Դ��з�����Ԫ�أ�
template <typename DerivedI>
int simplyTrisConnectedRegion(Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
	int simplyTrisConnectedRegion(											��������Ƭ����ͨ���������������ʱ����-1
				Eigen::VectorXi& connectedLabels,						ͬһ��ͨ�����ڵ�����Ƭ��ǩ������Ϊ0, 1, 2,...��
																										����Ϊi������Ƭ�ı�ǩΪconnectedLabels(i);
				Eigen::VectorXi& connectedCount							ÿ����ǩ��Ӧ�ĵ���ͨ�����ڰ���������Ƭ��
																										��ǩΪi�ĵ���ͨ�������������Ƭ��ΪconnectedCount(i)
				const Eigen::PlainObjectBase<DerivedI>& tris
				)
*/
	const unsigned trisCount = tris.rows();
	connectedLabels = Eigen::VectorXi::Zero(trisCount);
	unsigned visitedCount = 0;

	// 1. ������������Ƭ�ڽӹ�ϵ��
	Eigen::MatrixXi ttAdj_mnEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency_new(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return -1;

	// 2. �ڽ�����Ƭ������ѭ����
	int currentLabel = 1;
	while (visitedCount < trisCount)
	{
		// w1. ѡȡ���ӣ�
		int seedIdx = 0;
		for (unsigned i = 0; i < trisCount; ++i)
		{
			if (0 == connectedLabels(i))
			{
				seedIdx = i;
				break;
			}
		}

		// w2. ��ʼ������
		std::unordered_set<int> workingQueue;
		workingQueue.insert(seedIdx);
		while (!workingQueue.empty())
		{
			// ww1. ����
			int cIdx = *workingQueue.begin();
			visitedCount++;
			workingQueue.erase(cIdx);
			connectedLabels(cIdx) = currentLabel;

			// ww2. ������ǰ����Ƭ�������������Ե�����Ƭ
			for (int i = 0; i < 3; ++i)
				if (ttAdj_mnEdge(cIdx, i) >= 0)
					if (0 == connectedLabels(ttAdj_mnEdge(cIdx, i)))			// ��ǩΪ0��ʾ������Ƭδ�����ʹ���
						workingQueue.insert(ttAdj_mnEdge(cIdx, i));					// ֻ����δ���ʵ�����Ƭ 

			// ww3. ������ǰ����Ƭ�ķ������������Ե�����Ƭ��
			for (int i = 0; i < 3; ++i)
				for (const auto& index : ttAdj_nmnOppEdge[cIdx][i])
					if (0 == index)
						workingQueue.insert(index);
		}

		// w3. ��ǰ��ͨ�����������
		currentLabel++;
	}

	// 3. 
	unsigned regionCount = currentLabel - 1;
	connectedLabels.array() -= 1;
	connectedCount = Eigen::VectorXi::Zero(regionCount);
	for (unsigned i = 0; i < trisCount; ++i)
		connectedCount[connectedLabels[i]]++;


	return static_cast<int>(regionCount);
}



// ��ȡ���������������ͨ���򣨶�����ͨ��
template <typename DerivedV, typename DerivedI>
bool simplyConnectedLargest(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();

	// 1. �����ڽӾ���
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx, adjSM;
	adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
	adjSM = adjSM_eCount;
	traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
		{
			if (iter.value() > 0)
				iter.valueRef() = 1;
		});

	// 2. ȷ�����������е���ͨ����
	Eigen::VectorXi connectedLabels, connectedCount;
	int conCount = simplyConnectedRegion(adjSM, connectedLabels, connectedCount);
	if (conCount < 0)
		return false;

	// 3. ȷ�������ͨ���򣨶�������ࣩ��
	int mainLabel = 0;								// �����ͨ����ı�ǩ��
	int mainLabelCount = 0;
	for (int i = 0; i < conCount; ++i)
	{
		if (connectedCount(i) > mainLabelCount)
		{
			mainLabel = i;
			mainLabelCount = connectedCount(i);
		}
	}

	// 4. ��ȡ�����ͨ����Ķ��㣺
	std::vector<int> mainLabelIdxes;
	mainLabelIdxes.reserve(versCount);
	for (unsigned i = 0; i < versCount; ++i)
		if (mainLabel == connectedLabels(i))
			mainLabelIdxes.push_back(i);
	mainLabelIdxes.shrink_to_fit();
	subFromIdxVec(versOut, vers, mainLabelIdxes);


	// 5. ��ȡ�����ͨ���������Ƭ��

	//		5.1 ������-������ӳ���
	std::vector<int> oldNewIdxInfo(versCount, -1);
	for (int i = 0; i < mainLabelIdxes.size(); ++i)
	{
		int oldIdx = mainLabelIdxes[i];
		oldNewIdxInfo[oldIdx] = i;
	}

	//		5.2 ����Ƭ�����еĶ�������ӳ�����������
	DerivedI trisCopy = tris;
	int* intPtr = trisCopy.data();
	for (int i = 0; i < trisCopy.size(); ++i)
	{
		int oldIdx = *intPtr;
		*intPtr = oldNewIdxInfo[oldIdx];
		intPtr++;
	}

	//		5.3 ��ȡ�����ͨ�����ڵ�����Ƭ��
	trisOut.setZero(trisCount, 3);
	int trisCountNew = 0;
	for (int i = 0; i < trisCount; ++i)
	{
		if (trisCopy(i, 0) >= 0 && trisCopy(i, 1) >= 0 && trisCopy(i, 2) >= 0)
		{
			trisOut.row(trisCountNew) = trisCopy.row(i);
			trisCountNew++;
		}
	}

	//		5.4 shrink_to_fit();
	trisOut.conservativeResize(trisCountNew, 3);


	return true;
}


// �����㵥��ͨ����ֽ�����
template <typename T>
bool simplyConnectedSplitMesh(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& compVers, std::vector<Eigen::MatrixXi>& compTris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
{
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	if (0 == vers.rows() || 0 == tris.rows())
		return false;

	compVers.clear();
	compTris.clear();

	// 0. Ԥ������ȥ���������㣬ȥ���Ƿ����ظ�����Ƭ��
	MatrixXT vers0;
	Eigen::MatrixXi tris0, trisTmp;
	trisTmp = tris;
	removeSickDupTris(vers, trisTmp);
	std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, tris);
	if (!isoVerIdxes.empty())
		removeIsoVers(vers0, tris0, vers, trisTmp, isoVerIdxes);
	else
	{
		vers0 = vers;
		tris0 = trisTmp;
	}

	// 1. ��ȡ���㵥��ͨ���ƣ�
	const int versCount = vers0.rows();
	int scCount = 0;							// ���㵥��ͨ���������
	int sctCount = 0;							// ����Ƭ����ͨ���������
	Eigen::VectorXi connectedLabels, connectedCount;

	Eigen::SparseMatrix<int> adjSM, adjSM_eCount, adjSM_eIdx;
	adjMatrix(adjSM_eCount, adjSM_eIdx, tris0);
	adjSM = adjSM_eCount;
	traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
		{
			if (iter.value() > 0)
				iter.valueRef() = 1;
		});
	scCount = simplyConnectedRegion(adjSM, connectedLabels, connectedCount);

	// 2. ���ɶ��㵥��ͨ����
	compVers.resize(scCount);
	compTris.resize(scCount);
	for (int i = 0; i < scCount; ++i)
	{
		Eigen::VectorXi flag = (connectedLabels.array() == i).select(Eigen::VectorXi::Ones(versCount), Eigen::VectorXi::Zero(versCount));
		std::vector<int> oldNewIdxInfo, newOldIdxInfo;
		subFromFlagVec(compVers[i], oldNewIdxInfo, newOldIdxInfo, vers0, flag);
		compTris[i] = tris0;
		int* ptrData = compTris[i].data();
		for (int k = 0; k < compTris[i].size(); ++k)
		{
			*ptrData = oldNewIdxInfo[*ptrData];					// ������ӳ�����������
			ptrData++;
		}
		removeSickDupTris(compVers[i], compTris[i]);		// ȥ���Ƿ�����Ƭ��
	}

	return true;
}



/////////////////////////////////////////////////////////////////////////////////////////////////// ����ȱ�ݼ����޸�


template <typename DerivedI>
void findRepTris(std::vector<int>& repIdxes, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	std::unordered_set<std::int64_t> codeSet;
	unsigned trisCount = tris.rows();
	repIdxes.reserve(trisCount);

	for (unsigned i = 0; i < trisCount; ++i)
	{
		std::int64_t code = encodeTriangle(tris(i, 0), tris(i, 1), tris(i, 2));
		auto retPair = codeSet.insert(code);
		if (!retPair.second)
			repIdxes.push_back(static_cast<int>(i));
	}
	repIdxes.shrink_to_fit();
}


// findOverlapTris()����Ѱ�������е��ص�����Ƭ����������Ƭ���������ķ�����
template <typename DerivedV, typename DerivedI>
int findOverLapTris(std::vector<std::pair<int, int>>& opTrisPairs, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx = 0, const double thetaThreshold = 0.99 * pi)
{ 
	int olCount = 0;

	const double dotThreshold1 = cos(thetaThreshold);
	const double dotThreshold2 = cos(pi - thetaThreshold);
	const int versCount = vers.rows();
	const int trisCount = tris.rows();
	const int edgesCount = 3 * trisCount;
	std::unordered_set<std::int64_t> opTriPairCodes;			// һ���ص�������Ƭ������ֵ��ǰС������У������std::int64_t

#ifdef  LOCAL_DEBUG
	std::set<int> olTriIdxes;
#endif 

	// 0. ���ʼ״̬��������������Ƭ�ķ���
	Eigen::MatrixXd triNorms;
	if (!getTriNorms(triNorms, vers, tris))				// �������˻�����Ƭ���ᵼ�¼�����Ƭ�������
		return false;

	// 1. ������Ƭ�ڽӹ�ϵ��
	Eigen::MatrixXi ttAdj_mnEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return false;

	// 2. ����Ƭ��������-Ѱ�����ڵ��ص�����Ƭ�Ե�ѭ����
	std::unordered_set<int> workingSet;								// ����Ѱ����������Ƭ�Ķ��У�
	std::vector<int> triTags(trisCount, 0);								// 0 : δ���ʣ�1: �ѷ���
	workingSet.insert(static_cast<int>(triIdx));						// ���в����һ������Ƭ�� 
	while (workingSet.size() > 0)
	{
		// 2.1 �׸�Ԫ�س��ӣ����뵽��ͨ����Ƭ�����У�
		int cTriIdx = *workingSet.begin();									// ��ǰ���������Ƭ����
		std::vector<int> adjTriIdxes;												// ��ǰ����Ƭ��δ�����ʵ���������Ƭ������
		Eigen::RowVector3d cTriNorm = triNorms.row(cTriIdx);
		unsigned adjTrisCount = 0;
		workingSet.erase(workingSet.begin());
		triTags[cTriIdx] = 1;											// ��ǰ����Ƭ���Ϊ�ѷ��ʣ�

		// 2.2 ȷ����ǰ����Ƭ��������������Ƭ
		for (int i = 0; i < 3; ++i)
		{
			int nbrTriIdx = ttAdj_mnEdge(cTriIdx, i);
			if (nbrTriIdx >= 0)
			{
				// a. ��ǰ��Ϊ���α�
				if (0 == triTags[nbrTriIdx])
					adjTriIdxes.push_back(nbrTriIdx);
			}
			else
			{
				// b. ��ǰ��Ϊ�����αߣ�
				auto& vec1 = ttAdj_nmnEdge[cTriIdx][i];						// ��ǰ�����α����ڵ���������Ƭ��
				auto& vec2 = ttAdj_nmnOppEdge[cTriIdx][i];				// ��ǰ�����αߵĶԱ����ڵ���������Ƭ��
				for (const auto& index : vec1)
				{
					if (index == cTriIdx)
						continue;
					else
						if (0 == triTags[index])
							adjTriIdxes.push_back(index);
				}

				for (const auto& index : vec2)
					if (0 == triTags[index])
						adjTriIdxes.push_back(index);
			}
		}

		// 2.3 ����-�ж���ǰ����Ƭ���ܱ�����Ƭ�Ƿ��ص���
		adjTrisCount = adjTriIdxes.size();
		if (0 == adjTrisCount)
			continue;
		for (unsigned k = 0; k < adjTrisCount; ++k)
		{
			int adjIdx = adjTriIdxes[k];
			double dotValue = cTriNorm.dot(triNorms.row(adjIdx));
			if (dotValue > dotThreshold2 || dotValue < dotThreshold1)
			{
				// ����ǰ����Ƭ����������Ƭ����ӽ���ͬ���෴�������������ƬͶӰ���䷨��ƽ�����Ƿ����ص����� 
				int A, B, C1, C2;				// AB�������������Ƭ����Ķ��㣬C1, C2�����Ǹ�������Ķ��㣻
				std::vector<int> ctv{ tris(cTriIdx, 0), tris(cTriIdx, 1), tris(cTriIdx, 2) };
				std::vector<int> adjtv{ tris(adjIdx, 0), tris(adjIdx, 1), tris(adjIdx, 2) };
				auto iterC = ctv.begin();
				for (; iterC != ctv.end(); iterC++)
				{
					int cIdx = *iterC;
					if (cIdx != adjtv[0] && cIdx != adjtv[1] && cIdx != adjtv[2])
					{
						C1 = cIdx;
						break;
					}
				}
				for (const auto& aIdx : adjtv)
				{
					if (aIdx != ctv[0] && aIdx != ctv[1] && aIdx != ctv[2])
					{
						C2 = aIdx;
						break;
					}
				}
				ctv.erase(iterC);
				A = ctv[0];
				B = ctv[1];

				// �ж�C1, C2�Ƿ���ͬ�ࣺ
				Eigen::RowVector3d AB = vers.row(B).array().cast<double>() - vers.row(A).array().cast<double>();
				Eigen::RowVector3d AC1 = vers.row(C1).array().cast<double>() - vers.row(A).array().cast<double>();
				Eigen::RowVector3d AC2 = vers.row(C2).array().cast<double>() - vers.row(A).array().cast<double>();
				Eigen::RowVector3d crossArrow1 = AB.cross(AC1);
				Eigen::RowVector3d crossArrow2 = AB.cross(AC2);
				if (crossArrow1.dot(cTriNorm) * crossArrow2.dot(cTriNorm) > 0)				// C1, C2�ڹ�����AB��ͬһ�ࣻ
				{
					auto retPair = opTriPairCodes.insert(encodeUedge(cTriIdx, adjIdx));
					if (retPair.second)
					{
						opTrisPairs.push_back(std::make_pair(cTriIdx, adjIdx));
						olCount++;								// �������ӣ�
#ifdef LOCAL_DEBUG
						olTriIdxes.insert(cTriIdx);
						olTriIdxes.insert(adjIdx);
#endif
					}
				}
			}
		}

		// 2.3 ���������в��뵱ǰ����Ƭδ�����ʵ���������Ƭ��
		for (const auto& index : adjTriIdxes)
			workingSet.insert(index);
	}

#ifdef LOCAL_DEBUG
	if (olTriIdxes.size() > 0)
	{
		Eigen::MatrixXd versCopy = vers.array().cast<double>();
		Eigen::MatrixXi olTris;
		subFromIdxCon(olTris, tris, olTriIdxes);
		objWriteMeshMat("E:/overlapTris.obj", versCopy, olTris);
	}
#endif 

	return olCount;
}



// �����������еķ���������ߣ���ĳ��������ظ���������1����ñ߼���Ա߶����ж�Ϊ�����αߣ�
template<typename DerivedI>
int nonManifoldEdges(Eigen::MatrixXi& nmnEdges, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	nmnEdges.resize(0, 0);
	Eigen::SparseMatrix<int> adjSM_eIdx, adjSM_eCount;
	if (!adjMatrix(adjSM_eCount, adjSM_eIdx, tris))
		return false;

	// ����adjSM_eCount�����������αߣ�
	std::unordered_set<std::int64_t> nmnEdgeIdxes;
	traverseSparseMatrix(adjSM_eCount, [&adjSM_eCount, &nmnEdgeIdxes](auto& iter)
		{
			if (iter.value() > 1)
			{
				std::int64_t code = encodeEdge(iter.row(), iter.col());
				nmnEdgeIdxes.insert(code);
				if (adjSM_eCount.coeff(iter.col(), iter.row()) > 0)						// ��ǰ�ߵĶԱ��п��ܲ����ڣ�
				{
					std::int64_t codeOpp = encodeEdge(iter.col(), iter.row());
					nmnEdgeIdxes.insert(codeOpp);
				}
			}
		});
	const int nmnCount = nmnEdgeIdxes.size();
	nmnEdges.resize(nmnCount, 2);
	int index = 0;
	for (const auto& code : nmnEdgeIdxes)
	{
		auto retPair = decodeEdge(code);
		nmnEdges(index, 0) = retPair.first;
		nmnEdges(index, 1) = retPair.second;
		index++;
	}

	return nmnCount;
}


// �����������еķ���������ߣ���ĳ��������ظ���������1����ñ߼���Ա߶����ж�Ϊ�����αߣ�
template<typename DerivedI>
int nonManifoldEdges(Eigen::MatrixXi& nmnEdges, std::vector<std::pair<std::pair<int, int>, int>>& nmnInfos, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	nmnEdges.resize(0, 0);
	Eigen::SparseMatrix<int> adjSM_eIdx, adjSM_eCount;
	if (!adjMatrix(adjSM_eCount, adjSM_eIdx, tris))
		return false;

	// ����adjSM_eCount�����������αߣ�
	std::unordered_map<std::int64_t, int> repCountMap;
	std::unordered_set<std::int64_t> nmnEdgeIdxes;
	traverseSparseMatrix(adjSM_eCount, [&adjSM_eCount, &nmnEdgeIdxes, &repCountMap](auto& iter)
		{
			if (iter.value() > 1)
			{
				std::int64_t code = encodeEdge(iter.row(), iter.col());
				nmnEdgeIdxes.insert(code);
				repCountMap.insert({ code, iter.value() });
				int oppEdgeCount = adjSM_eCount.coeff(iter.col(), iter.row());
				if (oppEdgeCount > 0)						// ��ǰ�ߵĶԱ��п��ܲ����ڣ�
				{
					std::int64_t codeOpp = encodeEdge(iter.col(), iter.row());
					nmnEdgeIdxes.insert(codeOpp);
					repCountMap.insert({ codeOpp, oppEdgeCount });
				}
			}
		});
	const int nmnCount = nmnEdgeIdxes.size();
	nmnEdges.resize(nmnCount, 2);
	int index = 0;
	for (const auto& code : nmnEdgeIdxes)
	{
		auto retPair = decodeEdge(code);
		nmnEdges(index, 0) = retPair.first;
		nmnEdges(index, 1) = retPair.second;
		index++;
	}

	// ���ɷ����α���Ϣ��
	nmnInfos.reserve(nmnCount);
	for (const auto& pair : repCountMap)
	{
		auto edgePair = decodeEdge(pair.first);
		nmnInfos.push_back(std::make_pair(std::make_pair(edgePair.first, edgePair.second), pair.second));
	}

	return nmnCount;
}


// �����������еķ���������ߣ�����1 
template <typename T>
int nonManifoldUEs(Eigen::MatrixXi& nmnUedges, const Eigen::SparseMatrix<T>& adjSM_eCount, \
	const Eigen::SparseMatrix<T>& adjSM_ueCount)
{
	nmnUedges.resize(0, 0);
	std::unordered_set<std::int64_t> nmnUeCodes;

	// ����adjSM_ueCount�����������αߣ�
	traverseSparseMatrix(adjSM_ueCount, [&adjSM_ueCount, &nmnUedges, &nmnUeCodes](auto& iter)
		{
			if (iter.valueRef() > 2)
			{
				std::int64_t code = encodeUedge(iter.row(), iter.col());
				nmnUeCodes.insert(code);
			}
		});
	const int nmnCount = nmnUeCodes.size();
	nmnUedges.resize(nmnCount, 2);
	int index = 0;
	for (const auto& code : nmnUeCodes)
	{
		auto retPair = decodeEdge(code);
		nmnUedges(index, 0) = retPair.first;
		nmnUedges(index, 1) = retPair.second;
		index++;
	}

	return nmnCount;
}


// �����������еķ���������ߣ�����2 
template <typename DerivedI>
int nonManifoldUEs(Eigen::MatrixXi& nmnUedges, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	Eigen::SparseMatrix<int> adjSM_eIdx, adjSM_eCount, adjSM_ueCount, tmpSpMat;
	if (!adjMatrix(adjSM_eCount, adjSM_eIdx, tris))
		return false;

	spMatTranspose(tmpSpMat, adjSM_eCount);
	adjSM_ueCount = adjSM_eCount + tmpSpMat;

	return nonManifoldUEs(nmnUedges, adjSM_eCount, adjSM_ueCount);
}


// �����������еķ���������ߣ�����3�����������������ߵ���ϸ��Ϣ��
template<typename DerivedI>
int nonManifoldUEs(Eigen::MatrixXi& nmnUedges, std::vector<std::tuple<std::pair<int, int>, int, int >>& nmnEdgeInfos, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		int nonManifoldUEs(																			���ؼ�⵽�ķ���������ߵ�������ʧ�ܷ���-1��
				Eigen::MatrixXi& nmnUedges,																nmnCount * 2, ���������������
				std::vector<std::pair<int, std::pair<int, int>>>& nmnEdgeInfos,			Ԫ�ص�first��������ظ��Ĵ�����second����������������ߵľ������ݣ�
				const Eigen::PlainObjectBase<DerivedI>& tris
				)
	*/
	nmnEdgeInfos.clear();
	Eigen::SparseMatrix<int> adjSM_eIdx, adjSM_eCount, adjSM_ueCount, tmpSpMat;
	if (!adjMatrix(adjSM_eCount, adjSM_eIdx, tris))
		return -1;

	spMatTranspose(tmpSpMat, adjSM_eCount);
	adjSM_ueCount = adjSM_eCount + tmpSpMat;

	int nmnUeCount = nonManifoldUEs(nmnUedges, adjSM_eCount, adjSM_ueCount);
	if (nmnUeCount < 0)
		return -1;

	nmnEdgeInfos.resize(nmnUeCount);
	for (int i = 0; i < nmnUeCount; ++i)
	{
		int vaIdx = nmnUedges(i, 0);
		int vbIdx = nmnUedges(i, 1);
		int repCount = adjSM_eCount.coeff(vaIdx, vbIdx);
		int repOppCount = adjSM_eCount.coeff(vbIdx, vaIdx);
		nmnEdgeInfos[i] = std::make_tuple(std::make_pair(vaIdx, vbIdx), repCount, repOppCount);
	}

	return nmnUeCount;
}


// �������ε㡪������������Ҫû���ظ�����Ƭ�ͷǷ�����Ƭ
template<typename DerivedV>
int nonManifoldVers(std::vector<int>& nmnVerIdxes, const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	nmnVerIdxes.clear();
	const int versCount = vers.rows();
	const int trisCount = tris.rows();
	std::unordered_multimap<int, std::int64_t> triMap;

	// 1. ���㶥�����Ե���������߱��룬�����ֵ�
	for (int i = 0; i < trisCount; ++i)
	{
		const int& vaIdx = tris(i, 0);
		const int& vbIdx = tris(i, 1);
		const int& vcIdx = tris(i, 2);
		triMap.insert({ vaIdx, encodeUedge(vbIdx, vcIdx) });
		triMap.insert({ vbIdx, encodeUedge(vcIdx, vaIdx) });
		triMap.insert({ vcIdx, encodeUedge(vaIdx, vbIdx) });
	}

	// 2. �������ж��㣬���������Ե���������ߣ�
	nmnVerIdxes.reserve(versCount);
	for (int i = 0; i < versCount; ++i)
	{
		// f1. ��ǰ�������Ե���������ߴ���˫������
		std::list<std::pair<int, int>> roundUes;
		auto iter = triMap.find(i);
		if (iter == triMap.end())
			continue;
		for (int k = 0; k < triMap.count(i); ++k)
			roundUes.push_back(decodeEdge((iter++)->second));

		// f2. �����������ѭ����
		std::unordered_set<int> connectedVers;
		std::pair<int, int> currentUe = roundUes.front();
		roundUes.pop_front();
		while (!roundUes.empty())
		{
			bool blFound = false;
			auto iter = roundUes.begin();
			connectedVers.insert(currentUe.first);
			connectedVers.insert(currentUe.second);

			// fw1. ����ǰ����ߺ�����ͨ����ߴ�����ͬ�Ķ��㣬����������ͨ�ģ�
			for (; iter != roundUes.end(); iter++)
			{
				auto iter1 = connectedVers.find(iter->first);
				auto iter2 = connectedVers.find(iter->second);
				if (iter1 != connectedVers.end() || iter2 != connectedVers.end())
				{
					blFound = true;
					currentUe = *iter;
					roundUes.erase(iter);
					break;
				}
			}

			// fw1. ����ǰ��ͨ����������ȫ����roundUes�����б�ıߣ�˵���˶�������������߸�������Σ���ǰ�����Ƿ����ε㣻
			if (!blFound && !roundUes.empty())
			{
				nmnVerIdxes.push_back(i);
				break;
			}
		}
	}
	nmnVerIdxes.shrink_to_fit();

	return nmnVerIdxes.size();
}



// correctTriDirs()�����������������Ƭ���򣬻�������Ƭ���������ķ�����ʹ��ǰ��Ҫ�ȴ��������е��ص�����Ƭ��
template <typename DerivedV, typename DerivedI>
int correctTriDirs(Eigen::PlainObjectBase<DerivedI>& trisOut, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris, const double thetaThreshold = 0.99 * pi)
{
	int corrCount = 0;
	int visitedCount = 0;
	const double dotThreshold = cos(thetaThreshold);
	const int versCount = vers.rows();
	const int trisCount = tris.rows();
	const int edgesCount = 3 * trisCount;
	trisOut = tris;

#ifdef  LOCAL_DEBUG
	std::vector<int> corrTriIdxes;							// �洢������������Ƭ��������
#endif 

	// 0. ���ʼ״̬��������������Ƭ�ķ���
	Eigen::MatrixXd triNorms;
	if (!getTriNorms(triNorms, vers, tris))				// �������˻�����Ƭ���ᵼ�¼�����Ƭ�������
		return false;

	// 1. ������Ƭ�ڽӹ�ϵ��
	Eigen::MatrixXi ttAdj_mnEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return false;

	// 2. ѭ�����������е�һ������Ƭt���� - tд��������� - ���t������������Ƭ������Ȼ�����
	int triIdx = 0;
	std::unordered_set<int> workingSet;								// ����Ѱ����������Ƭ�Ķ��У�
	std::vector<int> triTags(trisCount, 0);								// 0 : δ���ʣ�1: �ѷ��ʣ�-1: �ѷ��ʣ������Ϊ�����������Ƭ��
	while (visitedCount < trisCount)
	{
		for (int i = 0; i < trisCount; ++i)
		{
			if (0 == triTags[i])
			{
				triIdx = i;
				break;
			}
		}

		workingSet.insert(static_cast<int>(triIdx));						// ���в����һ������Ƭ�� 
		while (workingSet.size() > 0)
		{
			// w1. �׸�Ԫ�س��ӣ����뵽��ͨ����Ƭ�����У�
			int cTriIdx = *workingSet.begin();									// ��ǰ���������Ƭ����
			workingSet.erase(workingSet.begin());
			if (triTags[cTriIdx] != 0)
				continue;
			triTags[cTriIdx] = 1;													// ��ǰ����Ƭ���Ϊ�ѷ��ʣ�
			visitedCount++;

			std::vector<int> adjTriIdxes;									// ��ǰ����Ƭ��δ�����ʵ���������Ƭ������
			Eigen::RowVector3d cTriNorm = triNorms.row(cTriIdx);
			Eigen::MatrixXd adjTriNorms;
			unsigned adjTrisCount = 0;

			// w2. ȷ����ǰ����Ƭ��������������Ƭ
			for (int i = 0; i < 3; ++i)
			{
				int nbrTriIdx = ttAdj_mnEdge(cTriIdx, i);
				if (nbrTriIdx >= 0)
					if (0 == triTags[nbrTriIdx])								// a. ��ǰ��Ϊ���α�
						adjTriIdxes.push_back(nbrTriIdx);
					else
					{
						// b. ��ǰ��Ϊ�����αߣ�
						auto& vec1 = ttAdj_nmnEdge[cTriIdx][i];						// ��ǰ�����α����ڵ���������Ƭ��
						auto& vec2 = ttAdj_nmnOppEdge[cTriIdx][i];				// ��ǰ�����αߵĶԱ����ڵ���������Ƭ��
						for (const auto& index : vec1)
						{
							if (index == cTriIdx)
								continue;
							else
								if (0 == triTags[index])
									adjTriIdxes.push_back(index);
						}

						for (const auto& index : vec2)
							if (0 == triTags[index])
								adjTriIdxes.push_back(index);
					}
			}

			// 4.3 ����-������ǰ����Ƭ����������Ƭ��
			adjTrisCount = adjTriIdxes.size();
			if (0 == adjTrisCount)
				continue;
			subFromIdxVec(adjTriNorms, triNorms, adjTriIdxes);
			for (unsigned k = 0; k < adjTrisCount; ++k)
			{
				int adjIdx = adjTriIdxes[k];
				if (cTriNorm.dot(adjTriNorms.row(k)) < dotThreshold)
				{
					// ����Ƭ���Ϊ-1��
					triTags[adjIdx] = -1;
					visitedCount++;
					corrCount++;												// �������ӣ�
#ifdef LOCAL_DEBUG
					corrTriIdxes.push_back(adjIdx);
#endif

					////		��������Ƭ����
					//int tmp = trisOut(adjIdx, 1);
					//trisOut(adjIdx, 1) = trisOut(adjIdx, 2);
					//trisOut(adjIdx, 2) = tmp;
					//triNorms.row(adjIdx).array() *= -1;			//	���·���
				}
			}

			// 4.3 ���������в��뵱ǰ����Ƭδ�����ʵ���������Ƭ��
			for (const auto& index : adjTriIdxes)
				if (0 == index)
					workingSet.insert(index);
		}
	}



	// for debug:
#if 0
	Eigen::MatrixXi sickTris;
	Eigen::MatrixXd versCopy = vers.array().cast<double>();
	subFromIdxVec(sickTris, tris, corrTriIdxes);
	objWriteMeshMat("E:/sickTris.obj", versCopy, sickTris);
#endif

	return corrCount;
}

// ������������㣬ֱ�ӶԶ����������ʷ֣� 
template <typename IndexType>
bool fillSmallHoles(Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic>& newTris, \
	const std::vector<std::vector<int>>& holes)
{
	/*
	bool fillSmallHoles(
		Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic>& newTris,				// �������ɵ���triangle soup
		const std::vector<std::vector<int>>& holes					// ������Ϣ��ÿ��Ԫ��������Ļ�·��������������
		)

	*/
	newTris.resize(0, 0);
	const int holesCount = holes.size();
	for (int i = 0; i < holesCount; ++i)
	{
		if (holes[i].size() < 3)				// invalid input;
			return false;

		Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic> currentNewTris;
		if (!circuitGetTris(currentNewTris, holes[i], true))
			return false;
		matInsertRows(newTris, currentNewTris);
	}

	return true;
}


// ��������еĹ�������
template <typename DerivedV>
std::vector<unsigned> checkIsoVers(const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();

	Eigen::VectorXi verIdxVec = Eigen::VectorXi::LinSpaced(versCount, 0, versCount - 1);

	std::vector<unsigned> isoVerIdxes;
	const int* dataPtr = tris.data();
	for (unsigned i = 0; i < tris.size(); ++i)
	{
		verIdxVec(*dataPtr) = -1;
		dataPtr++;
	}
	for (unsigned i = 0; i < versCount; ++i)
		if (verIdxVec(i) >= 0)
			isoVerIdxes.push_back(i);
 
	return isoVerIdxes;
}


// ȥ�������еĹ������㣺
template <typename DerivedV>
bool removeIsoVers(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const std::vector<unsigned>& isoVerIdxes)
{
	versOut.resize(0, 0);
	trisOut.resize(0, 0);
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	unsigned isoVersCount = isoVerIdxes.size();
	unsigned versCountNew = versCount - isoVersCount;

	std::vector<int> oldNewIdxInfo, newOldIdxInfo, tmpIdxVec;
	tmpIdxVec.resize(versCount);
	for (int i = 0; i < versCount; ++i)
		tmpIdxVec[i] = i;
	for (const auto& index : isoVerIdxes)
		tmpIdxVec[index] = -1;

	newOldIdxInfo.reserve(versCountNew);
	for (const auto& index : tmpIdxVec)
		if (index >= 0)
			newOldIdxInfo.push_back(index);

	int tmpIdx = 0;
	oldNewIdxInfo = tmpIdxVec;
	for (auto& index : oldNewIdxInfo)
		if (index >= 0)
			index = tmpIdx++;

	// ������ƣ�
	subFromIdxVec(versOut, vers, newOldIdxInfo);

	// �������Ƭ��
	trisOut = tris;
	int* dataPtr = trisOut.data();
	for (int i = 0; i < trisCount * 3; ++i)
		*dataPtr++ = oldNewIdxInfo[*dataPtr];

	return true;
}


// ����������Ƿ��зǷ�����Ƭ������������������������ͬ�ģ�
template <typename DerivedI>
bool checkSickTris(std::vector<int>& sickIdxes, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	bool retFlag = false;
	for (int i = 0; i < tris.rows(); ++i)
	{
		if (tris(i, 0) == tris(i, 1) || tris(i, 0) == tris(i, 2) || tris(i, 1) == tris(i, 2))
		{
			sickIdxes.push_back(i);
			retFlag = true;
		}
	}
	return retFlag;
}




// ����˻��ߣ��߳����̵ıߣ�
template <typename DerivedV>
Eigen::VectorXi checkDegEdges(const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const double threshold = 1e-3)
{
	/*
		Eigen::VectorXi checkDegEdges(														// ���ر���˻��ߵ�flag����
				const Eigen::MatrixXi& edges,												// ���������
				const Eigen::PlainObjectBase<DerivedV>& edgeArrows,		// ����ߵı���������
				const Eigen::PlainObjectBase<DerivedV>& vers,
				const Eigen::MatrixXi& tris,
				const double threshold = 1e-3												// �ж��˻��ߵı߳���ֵ��
				)
	*/
	Eigen::VectorXd edgesLen = edgeArrows.rowwise().norm();
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	unsigned edgesCount = edges.rows();

	Eigen::VectorXi flags = (edgesLen.array() <= threshold).select(Eigen::VectorXi::Ones(edgesCount), Eigen::VectorXi::Zero(edgesCount));

	return flags;
}


// ȥ�������˻��ߣ��߳����̵ıߣ������������ںϣ�;
template <typename DerivedV>
int mergeDegEdges(Eigen::PlainObjectBase<DerivedV>& newVers, Eigen::MatrixXi& newTris, \
	const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const Eigen::VectorXi& degEdgeFlags)
{
	/*
		int mergeDegEdges(																				����ȥ�����Ķ���������ʧ���򷵻�-1��
				Eigen::PlainObjectBase<DerivedV>& newVers,
				Eigen::MatrixXi& newTris,
				const Eigen::MatrixXi& edges,
				const Eigen::PlainObjectBase<DerivedV>& edgeArrows,
				const Eigen::PlainObjectBase<DerivedV>& vers,
				const Eigen::MatrixXi& tris,
				const Eigen::VectorXi& degEdgeFlags								����˻��ߵ�flag����
				)

	*/
	int repVersCount = -1;
	int* dataPtr = nullptr;
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	unsigned edgesCount = edges.rows();

	// 1. ��ȡ�˻������
	int degCount = degEdgeFlags.sum();				// �˻�����ߵ�����
	Eigen::MatrixXi degEdges;
	subFromFlagVec(degEdges, edges, degEdgeFlags);

	// for debug:
	if (1)
	{
		Eigen::MatrixXd degEdgeVers;
		std::unordered_set<int> tmpSet;
		std::vector<int> degVerIdxes;

		objWriteEdgesMat("E:/degEdges.obj", degEdges, vers);
		for (unsigned i = 0; i < degCount; ++i)
		{
			tmpSet.insert(degEdges(i, 0));
			tmpSet.insert(degEdges(i, 1));
		}
		degVerIdxes.insert(degVerIdxes.end(), tmpSet.begin(), tmpSet.end());
		subFromIdxVec(degEdgeVers, vers, degVerIdxes);
		objWriteVerticesMat("E:/degEdgeVers.obj", degEdgeVers);
	}

	// �����˻����У�����������С�Ķ��㣬������Ķ����дΪ����С�Ķ��㣺
	Eigen::MatrixXi revDegEdges = degEdges;
	for (unsigned i = 0; i < degCount; ++i)
	{
		int vaIdx = degEdges(i, 0);
		int vbIdx = degEdges(i, 1);
		if (degEdges(i, 0) < degEdges(i, 1))
			revDegEdges.row(i) = Eigen::RowVector2i(vbIdx, vaIdx);
		else
			degEdges.row(i) = Eigen::RowVector2i(vbIdx, vaIdx);
	}

	std::list<std::set<int>> clusters;
	for (unsigned i = 0; i < degCount; ++i)
	{
		int vaIdx = degEdges(i, 0);
		int vbIdx = degEdges(i, 1);
		for (auto& clusterSet : clusters)
		{
			auto retIter1 = clusterSet.find(vaIdx);
			auto retIter2 = clusterSet.find(vbIdx);
			if (clusterSet.end() == retIter1 && clusterSet.end() == retIter2)
				continue;
			else       // ��ǰ�ߵĶ�����֮ǰ�Ĵ��еĶ�����ͬ�����뵽�ô��У�
			{
				clusterSet.insert(vaIdx);
				clusterSet.insert(vbIdx);
				break;
			}
		}

		// �����µĴأ�
		std::set<int> tmpCluster;
		tmpCluster.insert(vaIdx);
		tmpCluster.insert(vbIdx);
		clusters.push_back(tmpCluster);
	}

	// ���ɶ�������ӳ���ֵ䡪��һ�����������������㶼ӳ��Ϊ������С���Ǹ����㣺
	std::map<int, int> clusterMap;
	for (const auto clusterSet : clusters)
	{
		auto iter = clusterSet.begin();
		int headIdx = *iter;
		iter++;
		for (; iter != clusterSet.end(); ++iter)
			clusterMap.insert({ *iter, headIdx });
	}

	std::vector<int> repIdxes;							// ��Ҫȥ���Ķ����������
	repIdxes.reserve(degCount);
	for (const auto& pair : clusterMap)
		repIdxes.push_back(pair.second);
	repIdxes.shrink_to_fit();
	repVersCount = repIdxes.size();

	std::map<int, int> fullMap = clusterMap;
	for (int i = 0; i < versCount; ++i)
		fullMap.insert({ i, i });

	// ����ӳ�䣺
	Eigen::MatrixXi trisCopy = tris;
	dataPtr = trisCopy.data();
	for (unsigned i = 0; i < 3 * trisCount; ++i)
	{
		int index = *dataPtr;
		*dataPtr = fullMap[index];
		dataPtr++;
	}

	// ɾ���Ƿ�����Ƭ��
	std::vector<int> sickIdxes;
	checkSickTris(sickIdxes, trisCopy);
	removeTris(newTris, trisCopy, sickIdxes);

	// ɾ���������㣺
	trisCopy = newTris;
	std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, trisCopy);
	newVers.resize(0, 0);
	newTris.resize(0, 0);
	removeIsoVers(newVers, newTris, vers, trisCopy, isoVerIdxes);

	return degCount;
}



// ȥ���Ƿ�����Ƭ���ظ�����Ƭ��
template<typename T>
int removeSickDupTris(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris)
{
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	int versCount = vers.rows();
	int trisCount = tris.rows();

	Eigen::VectorXi flags{ Eigen::VectorXi::Ones(trisCount) };			// flag��������Ҫȥ��������Ƭ��0��
	std::unordered_set<std::uint64_t> triCodeSet;

	// ��������Ƭ
	for (int i = 0; i < trisCount; ++i)
	{
		// ����Ƭ�е�����ֵ����0~versCount-1֮����Ϊ�Ƿ�����Ƭ��ò����ЩSTL�ļ��л�����������Σ�
		if (tris(i, 0) < 0 || tris(i, 0) > versCount - 1 || tris(i, 1) < 0 || tris(i, 1) > versCount - 1 || tris(i, 2) < 0 || tris(i, 2) > versCount - 1)
		{
			flags(i) = 0;
			continue;
		}

		std::uint64_t code = encodeTriangle(tris(i, 0), tris(i, 1), tris(i, 2));
		if (0 == code)
		{
			flags(i) = 0;				// ����Ƭ����������������ͬ�ģ��Ƿ���
			continue;
		}

		auto retPair = triCodeSet.insert(code);
		if (!retPair.second)
			flags(i) = 0;				// �ظ�����Ƭ�Ƿ���
	}

	Eigen::MatrixXi tmpMat;
	subFromFlagVec(tmpMat, tris, flags);
	tris = tmpMat;

	int removeCount = trisCount - tris.rows();
	return removeCount;
}


// ������������ɱ�Ե��·(�����ı�Ե)���������Ƕ����ߵ����Ρ���������ǰ����̫���ӵĶ������صĶ���������˳����ܲ��ԡ�Ч������VBFindHole2
template <typename DerivedV, typename DerivedI>
int findHoles(std::vector<std::vector<int>>& holes, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	// ע: �������񲻿��Դ��ڷǷ�����Ƭ���ظ�����Ƭ��
	/*
		int findHoles(																			�ɹ����ض�����ʧ�ܷ���-1;
			std::vector<std::vector<int>>& holes,								�������պ�Ϊ��·�ı�Ե���㣬ÿȦ���㰴˳ʱ�����У�
			const Eigen::PlainObjectBase<DerivedV>& vers,
			const Eigen::PlainObjectBase<DerivedI>& tris
			)
	*/
	// lambda������������Ϣ��
	auto organizeHoles = [](std::vector<std::vector<int>>& holes)		// ����ÿ�����Ķ���˳��ʹ���׸�������������С���Ǹ�
	{
		for (auto& vec : holes)
		{
			auto iter = std::min_element(vec.begin(), vec.end());
			std::vector<int> tmpVec;
			tmpVec.insert(tmpVec.end(), iter, vec.end());
			tmpVec.insert(tmpVec.end(), vec.begin(), iter);
			vec = tmpVec;
		}
	};

	// ȷ�����ıߡ���ֻ����һ������Ƭ�ıߣ�
	const int  trisCount = tris.rows();
	const int  edgesCount = 3 * trisCount;
	const int  versCount = tris.maxCoeff() + 1;
	holes.clear();

	// 1. ��ȡ���б�Ե����ߣ�
	std::vector<Eigen::MatrixXi> circs;
	Eigen::MatrixXi bdrys;
	bdryEdges(bdrys, tris);
	const int bdrysCount = bdrys.rows();

	// 2. ���ɱ�Ե������
	std::list<std::pair<int, int>> bdryList;
	for (int i = 0; i < bdrysCount; ++i)
		bdryList.push_back({ bdrys(i, 0), bdrys(i, 1) });

	// 3. �����������Ѱ�����еı�Ե��·��
	std::list<std::list<int>> holeList;
	std::pair<int, int> currentEdge = bdryList.front();
	int currentHead = currentEdge.first;								// ��ǰ�����׶���������
	holeList.push_back(std::list<int>{currentEdge.first, currentEdge.second});
	bdryList.pop_front();
	while (!bdryList.empty())
	{
		int vaIdx1 = currentEdge.first;
		int vbIdx1 = currentEdge.second;
		std::list<int>& currentHole = holeList.back();
		auto iter = bdryList.begin();
		bool blFound = false;
		for (; iter != bdryList.end(); ++iter)
		{
			int vaIdx2 = iter->first;
			int vbIdx2 = iter->second;
			if (vbIdx1 == vaIdx2)
			{
				currentEdge = *iter;
				bdryList.erase(iter);
				currentHole.push_back(vbIdx2);
				blFound = true;
				break;
			}
		}

		if (!blFound)
		{
			currentEdge = bdryList.front();
			currentHead = currentEdge.first;
			holeList.push_back(std::list<int>{currentEdge.first, currentEdge.second});
			bdryList.pop_front();
		}
	}

	// 4. ���������
	const int holesCount = holeList.size();
	holes.reserve(holesCount);
	for (const auto& list : holeList)
	{
		const int holeVersCount = list.size() - 1;
		auto iterList = list.rbegin();									// ʹ�÷��������������������䣻
		std::vector<int> currentHole(holeVersCount);
		if (holeVersCount < 3)
			return -1;							// û���γɻ�·��
		if (*list.begin() != *list.rbegin())
			return -1;							// û���γɻ�·��
		for (int i = 0; i < holeVersCount; ++i)
			currentHole[i] = *iterList++;
		holes.push_back(currentHole);
	}
	organizeHoles(holes);

	return holesCount;
}





