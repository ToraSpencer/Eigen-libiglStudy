
/////////////////////////////////////////////////////////////////////////////////////////////////// 图形属性：



// 得到有向边编码-有向边索引的字典：
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
	// ！！！非流形网格不可使用；
	/*
		bool getUedges(
				Eigen::MatrixXi& uEdges,													无向边
				Eigen::VectorXi& edgeUeInfo,												edgeUeInfo(i)是索引为i的有向边对应的无向边的索引；
				Eigen::MatrixXi& ueEdgeInfo												ueCount * 2的矩阵；
																													ueEdgeInfo(i, 0)和ueEdgeInfo(i, 1)是索引为i的无向边对应的两条有向边的索引；
																													正序有向边在前，逆序有向边在后；
																													若当前无向边是边缘边，则只对应一条有向边，ueEdgeInfo(i, 1)==-1
				const Eigen::PlainObjectBase<DerivedI>& edges
				)
	*/
	edgeUeInfo.resize(0);
	if (!getUedges(uEdges, edges))
		return false;

	const unsigned uEdgesCount = uEdges.rows();
	const unsigned edgesCount = edges.rows();

	// 确定有向边索引→无向边索引的映射；
	std::unordered_map<std::int64_t, int> codeMap;
	for (unsigned i = 0; i < uEdgesCount; ++i)
		codeMap.insert({ encodeUedge(uEdges(i, 0), uEdges(i, 1)), i });
	edgeUeInfo.resize(edgesCount);
	for (unsigned i = 0; i < edgesCount; ++i)
	{
		std::int64_t code = encodeUedge(edges(i, 0), edges(i, 1));
		edgeUeInfo(i) = codeMap[code];
	}

	// 确定无向边索引→有向边索引的映射；
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

	// ueEdgeInfo中每行两条无向边排序——正序在前逆序在后；
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
				Eigen::MatrixXi& UeTrisInfo,						ueCount×2，无向边关联的三角片索引矩阵；
																					若索引为i的无向边为(A,B)，则UeTrisInfo(i, 0)为有向边(A,B)所在的三角片索引，
																							UeTrisInfo(i, 1)为有向边(B, A)所在的三角片索引。
				Eigen::MatrixXi& UeCornersInfo,				ueCount×2，无向边所对的顶点标记；
																					顶点标记：0-a, 1-b, 2-c
																					若索引为i的无向边为(A,B)，则UeCornersInfo(i, 0)为有向边(A,B)所对的顶点的标记，
																							则UeCornersInfo(i, 1)为有向边(B,A)所对的顶点的标记。

				const Eigen::VectorXi& edgeUeInfo,										有向边索引到无向边索引的映射；
				const Eigen::PlainObjectBase<DerivedI>& uEdges,				无向边数据
				const Eigen::PlainObjectBase<DerivedI>& tris						三角片数据
				)

		注意：无向边统一使用正序表示，即vaIdx < vbIdx;
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
				// 若当前有向边ab和无向边AB相同
				UeTrisInfo(ueIdx, 0) = i;
				UeCornersInfo(ueIdx, 0) = k;
			}
			else
			{
				// 若当前有向边ab和无向边AB相反
				UeTrisInfo(ueIdx, 1) = i;
				UeCornersInfo(ueIdx, 1) = k;
			}
		}
	}

	return true;
}


// 计算网格所有三角片法向：
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


// 得到坐标向量表示的有向边：
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


// 计算网格所有三角片的重心
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


// 计算网格所有三角片的归一化法向量
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

	// 法向量归一化，若存在退化三角片，则写为inf
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


// 计算网格所有三角片所在平面
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

		// 若存在退化三角片，则平面系数都写为NAN:
		if (std::isinf(norm(0)))
			planeCoeff(i, 3) = INFINITY;

		Eigen::RowVector3d va = vers.row(tris(i, 0)).array().cast<double>();

		// va点在平面上 → va点到平面的距离为0 → norm.dot(va) +d == 0; → d = -norm.dot(va);
		planeCoeff(i, 3) = -norm.dot(va);
	}

	return true;
}


// 计算三角网格每个三角片的面积：
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
	s = (lens1 + lens2 + lens3) / 2.0;				// 三角片周长的一半；

	// 使用Heron公式计算三角片面积:
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


// 计算三角网格的体积：
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

	// 一个三角片对应的四面体的符号体积（signed volume） V(t0) == (-x3y2z1 + x2y3z1 + x3y1z2 - x1y3z2 - x2y1z3 + x1y2z3) / 6.0;
	auto volumeVec = \
		(-versC.col(0).array() * versB.col(1).array() * versA.col(2).array() + versB.col(0).array() * versC.col(1).array() * versA.col(2).array()\
			+ versC.col(0).array() * versA.col(1).array() * versB.col(2).array() - versA.col(0).array() * versC.col(1).array() * versB.col(2).array()\
			- versB.col(0).array() * versA.col(1).array() * versC.col(2).array() + versA.col(0).array() * versB.col(1).array() * versC.col(2).array()) / 6.0;

	volume = volumeVec.sum();

	return volume;
}


// 求三角网格的不同权重的邻接矩阵 
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
	adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());										// 权重为该有向边重复的次数；
	adjSM_eIdx.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());		// 权重为该有向边的索引；

	return true;
}


// 边缘有向边，重载1：
template<typename DerivedI>
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
	adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
	Eigen::SparseMatrix<int> adjSM_tmp, adjSM_ueCount;

	// 1. 确定边缘无向边：
	spMatTranspose(adjSM_tmp, adjSM_eCount);
	adjSM_ueCount = adjSM_eCount + adjSM_tmp;				// 无向边邻接矩阵，权重是该边的重复次数；

	std::vector<std::pair<int, int>> bdryUeVec;						// 所有边缘无向边；约定无向边表示为前大后小的索引对；
	traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
		{
			if (1 == iter.value() && iter.row() > iter.col())
				bdryUeVec.push_back({ iter.row(), iter.col() });
		});

	// 2. 由边缘无向边确定边缘有向边：
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


// 边缘有向边， 重载2：
template<typename DerivedI>
bool bdryEdges(Eigen::PlainObjectBase<DerivedI>& bdrys, std::vector<int>& bdryTriIdxes, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		bool bdryEdges(
				Eigen::PlainObjectBase<DerivedI>& bdrys,					bdrysCount * 2，边缘边数据，是有向边；
				std::vector<int>& bdryTriIdxes,										包含边缘边的三角片索引；
				const Eigen::PlainObjectBase<DerivedI>& tris				trisCount * 3, 网格三角片数据；
				)

	*/

	// 1. 查找网格的边缘有向边；
	if (!bdryEdges(bdrys, tris))
		return false;

	if (0 == bdrys.rows())
		return true;

	const unsigned trisCount = tris.rows();

	// 2. 确认边缘非流形边所在的三角片索引：
	Eigen::MatrixXi edgeAs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi edgeBs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi edgeCs = Eigen::MatrixXi::Zero(trisCount, 2);
	Eigen::MatrixXi vaIdxes = tris.col(0);
	Eigen::MatrixXi vbIdxes = tris.col(1);
	Eigen::MatrixXi vcIdxes = tris.col(2);

	// 3. 生成有向边编码-三角片字典：
	std::unordered_multimap<std::int64_t, unsigned> edgesMap;			// 一条有向边正常情况下只关联一个三角片，若关联多个，则是非流形边；
	for (unsigned i = 0; i < trisCount; ++i)
	{
		std::int64_t codeA = encodeEdge(vbIdxes(i), vcIdxes(i));
		std::int64_t codeB = encodeEdge(vcIdxes(i), vaIdxes(i));
		std::int64_t codeC = encodeEdge(vaIdxes(i), vbIdxes(i));
		edgesMap.insert({ codeA, i });
		edgesMap.insert({ codeB, i });
		edgesMap.insert({ codeC, i });
	}

	// 4. 在字典中查找所有边缘边所在的三角片索引：
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


// 计算网格中三角片三个角的余切值；
template <typename Tv, typename Tl>
void trisCotValues(const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers,
	const Eigen::MatrixXi& tris, Eigen::Matrix<Tl, Eigen::Dynamic, Eigen::Dynamic>& cotValues)
{
	// 1/2*cotangents corresponding angles. for triangles, columns correspond to edges bc, ca, ab;
	const int trisCount = tris.rows();

	// 计算每个三角片的面积的两倍：
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



// buildAdjacency()——计算网格的三角片邻接信息：
using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool buildAdjacency(Eigen::MatrixXi& ttAdj_mnEdge, std::vector<tVec>& ttAdj_nmnEdge, std::vector<tVec>& ttAdj_nmnOppEdge, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		bool buildAdjacency(
					const Eigen::MatrixXi& tris,									输入的三角片数据

					Eigen::MatrixXi& ttAdj_mnEdge,							三角片的非边缘流形有向边邻接的三角片索引；
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

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_eIdx的triplet数据；
	std::unordered_multimap<std::int64_t, int> edgeMap;					// 边编码——64位整型数数表示的边数据（两个端点索引）；


	// 1. 求基本的边信息、邻接矩阵：
	Eigen::SparseMatrix<int> adjSM;								// 邻接矩阵；
	Eigen::SparseMatrix<int> adjSM_eCount;					// 邻接矩阵，索引为该有向边重复的次数；
	Eigen::SparseMatrix<int> adjSM_eIdx;						// 邻接矩阵，元素若对应流形边则为该边索引，若对应非流形边则为该边所有索引的和；
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;				// adjSM_eIdx的转置，表示对边的信息；
	Eigen::SparseMatrix<int> adjSM_ueCount;				// 权重为无向边ij关联的三角片数量；
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB;
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;
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
		getEdges(edges, tris);

		// 1.2 生成三种有向边邻接矩阵：adjSM, adjSM_eCount, adjSM_eCount;
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris, edges);
		adjSM = adjSM_eCount;																		// 有向边邻接矩阵；
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		若adjSM_eIdx(i, j)对应的是流形有向边，该权重值为该边的索引；若是非流形有向边，则该边对应两个边索引，权重为两个索引之和；
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

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
		spMatTranspose(adjSM_MN_NB_opp, adjSM_MN_NB);
	}


	// 2. 确定非边缘流形有向边、及其对边的索引；
	std::vector<int> edgesIdx_MN_NB;							// 非边缘流形有向边的索引；
	std::vector<int> edgesIdx_MN_NB_opp;					// 非边缘流形有向边的对边的索引；
	{
		unsigned edgesCount_MN_NB = adjSM_MN_NB.sum();
		edgesIdx_MN_NB.reserve(edgesCount_MN_NB);					// 非边缘流形有向边索引；
		edgesIdx_MN_NB_opp.reserve(edgesCount_MN_NB);			// 非边缘流形有向边的对边的索引；

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB.push_back(adjSM_eIdx.coeffRef(iter.row(), iter.col()));
			});

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB_opp.push_back(adjSM_eIdx_opp.coeffRef(iter.row(), iter.col()));
			});
	}


	// 3. 流形边-三角片邻接关系，三角片邻接关系：
	std::vector<int> etInfo;						// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
	Eigen::VectorXi etAdj_mnEdge;							// 边所在三角片的索引，非流形边或边缘边写-1；
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

			if (edgesIdxOpp < 0 || edgesIdxOpp >= edgesCount)
				std::cout << "pause" << std::endl;

			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//	3.3 求三角片邻接矩阵；ttAdj_mnEdge(i, :)是索引为i的三角片三条边邻接的三个三角片，若其中有非流形边或边缘边则写为-1；
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. 计算非流形边（关联两个三角片的有向边）信息
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_map;				// 非流形有向边的索引 - 该边对应的所有索引；
	std::unordered_map<int, std::vector<int>> edgeIdx_nmn_opp_map;		// 非流形有向边的索引 - 该边的对边对应的所有索引；
	{
		//	4.1 遍历邻接矩阵找出所有非流形边（关联两个三角片的有向边）；
		std::vector<int> edgeIdx_nmn;						// 非流形有向边索引；
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)) > 1)
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();

		int neCount = edgeIdx_nmn.size();						// 非流形有向边个数；
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
			// 注意！！存在两种以外情况： a. 当前边是非流形边，对边是流形边；	b. 当前边是非流形边，对边不存在；
			int neIdx = edgeIdx_nmn[i];		// 当前非流形边索引；
			int vaIdx = nmnEdges(i, 0);
			int vbIdx = nmnEdges(i, 1);
			std::int64_t eCode = encodeEdge(vaIdx, vbIdx);
			std::int64_t oppEcode = encodeEdge(vbIdx, vaIdx);
			auto iter = edgeMap.find(eCode);
			auto oppIter = edgeMap.find(oppEcode);

			int otherIdx = (iter->second == neIdx) ? ((++iter)->second) : (iter->second);
			edgeIdx_nmn_map.insert({ neIdx, std::vector<int>{neIdx, otherIdx} });

			if (edgeMap.end() == oppIter)			// 当前边是非流形边，对边不存在；
			{
				edgeIdx_nmn_opp_map.insert({ neIdx, std::vector<int>{} });
				continue;
			}
			int oppIdx1 = oppIter->second;
			int oppIdx2 = (++oppIter)->second;
			edgeIdx_nmn_opp_map.insert({ neIdx, std::vector<int>{oppIdx1, oppIdx2} });
		}
	}


	// 5. 求含有非流形边的三角片的邻接关系：
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


// buildAdjacency_new()——计算网格的三角片邻接信息：
using tVec = std::vector<std::vector<int>>;
template <typename DerivedI>
bool buildAdjacency_new(Eigen::MatrixXi& ttAdj_mnEdge, std::vector<tVec>& ttAdj_nmnEdge, std::vector<tVec>& ttAdj_nmnOppEdge, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		bool buildAdjacency_new(
					const Eigen::MatrixXi& tris,									输入的三角片数据

					Eigen::MatrixXi& ttAdj_mnEdge,							ttAdj_mnEdge(i, :)是索引为i的三角片三条边邻接的三个三角片
																									trisCount * 3
																									(i, 0)为索引为i的三角片的流形有向边(vb, vc)（非边缘）邻接的三角片的索引；
																									若(vb, vc)为非流形有向边，则(i, 0)为-1；
																									若(vb, vc)为边缘有向边，则(i, 0)为-2；

					std::vector<tVec>& ttAdj_nmnEdge,				三角片的非流形有向边所在的三角片索引；
																									size == trisCount;
																									每个tVec元素包含三个int向量；
																									ttAdj_nmnEdge[i][0]，为索引为i的三角片的非流形有向边(vb, vc)所在的所有三角片索引；


					std::vector<tVec>& ttAdj_nmnOppEdge		三角片的非流形有向边邻接的三角片索引（即对边所在的三角片的索引）；
																									size == trisCount;
																									ttAdj_nmnOppEdge[i][0]，为索引为i的三角片的非流形有向边(vb, vc)邻接的所有三角片索引；
					)

	*/
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	Eigen::MatrixXi edges;							// 有向边数据；
	Eigen::MatrixXi nmnEdges;					// 非流形有向边

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_eIdx的triplet数据；

	// 1. 求基本的边信息、邻接矩阵： 
	Eigen::SparseMatrix<int> adjSM_eCount;					// 邻接矩阵，索引为该有向边重复的次数；
	Eigen::SparseMatrix<int> adjSM_eCount_opp;
	Eigen::SparseMatrix<int> adjSM_eIdx;						// 邻接矩阵，元素若对应流形边则为该边索引，若对应非流形边则为该边所有索引的和；
	Eigen::SparseMatrix<int> adjSM_ueCount;				// 权重为无向边ij关联的三角片数量； 
	Eigen::SparseMatrix<int> adjSM_MN;						// 流形有向边（非边缘）邻接矩阵，权重为1； 
	std::unordered_set<int> bdryIdxes;									// 边缘有向边的索引；
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
		getEdges(edges, tris);

		// 1.2 生成三种有向边邻接矩阵：adjSM, adjSM_eIdx, adjSM_eCount;
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris, edges);
		spMatTranspose(adjSM_eCount_opp, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + adjSM_eCount_opp;						// 无向边邻接矩阵，权重为无向边ij关联的三角片数量； 						 

		//	1.3 生成流形有向边（非边缘）的邻接矩阵
		adjSM_MN = adjSM_ueCount;
		traverseSparseMatrix(adjSM_MN, [&adjSM_MN](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});
		adjSM_MN.prune([&adjSM_eCount](const Eigen::Index& row, const Eigen::Index& col, const int& value)->bool
			{
				if (1 == adjSM_eCount.coeffRef(row, col) && 1 == adjSM_eCount.coeffRef(col, row))
					return true;								// 返回true的元素被保留；
				else
					return false;
			});															// 是对称矩阵；

		// 1.4 确定边缘有向边的索引： 
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


	// 2. 确定流形有向边（非边缘）、及其对边的索引；
	std::vector<int> edgesIdx_MN;							// 流形有向边（非边缘）的索引；
	std::vector<int> edgesIdx_MN_opp;					// 流形有向边（非边缘）的对边的索引；
	unsigned edgesCount_MN = 0;
	{
		edgesCount_MN = adjSM_MN.sum();
		edgesIdx_MN.reserve(edgesCount_MN);					// 流形有向边（非边缘）索引；
		edgesIdx_MN_opp.reserve(edgesCount_MN);			// 流形有向边（非边缘）的对边的索引；

		traverseSparseMatrix(adjSM_MN, [&](auto& iter)
			{
				edgesIdx_MN.push_back(adjSM_eIdx.coeffRef(iter.row(), iter.col()));
				edgesIdx_MN_opp.push_back(adjSM_eIdx.coeffRef(iter.col(), iter.row()));
			});
	}

	// 3. 流形边-三角片邻接关系，三角片邻接关系：
	std::vector<int> etInfo;									// 有向边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
	Eigen::VectorXi etAdj_mnEdge;						// etAdj_mnEdge(i)是索引为i的有向边所对三角片的索引，没有的话写-1；
	{
		// 3.1 生成边索引 - 三角片索引映射表etInfo;
		etInfo.resize(edgesCount);							// etInfo(i)是索引为i的边所在的三角片的索引；
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

		// 3.2 求有向边所在三角片的索引
		etAdj_mnEdge = -Eigen::VectorXi::Ones(edgesCount);
		for (unsigned i = 0; i < edgesCount_MN; ++i)
		{
			const int& edgeIdx = edgesIdx_MN[i];
			const int& edgesIdxOpp = edgesIdx_MN_opp[i];
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		// 3.3 边缘边赋值为-2；
		for (const auto& index : bdryIdxes)
			etAdj_mnEdge(index) = -2;

		//	3.3 求三角片邻接矩阵；
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}

	// 4. 计算非流形边（关联两个三角片的有向边）信息
	std::unordered_map<std::int64_t, std::unordered_set<int>> edgeMap;							// 边编码——64位整型数数表示的边数据（两个端点索引）；
	std::unordered_map<int, std::unordered_set<int>> edgeIdx_nmn_map;						// 非流形有向边的索引 - 该边对应的所有索引；
	std::unordered_map<int, std::unordered_set<int>> edgeIdx_nmn_opp_map;				// 非流形有向边的索引 - 该边的对边对应的所有索引；
	{
		//	4.1 遍历邻接矩阵找出所有非流形有向边——若某条有向边或其对边重复次数大于1，则这两条边都被视为非流形有向边；
		std::vector<int> edgeIdx_nmn;							// 非流形有向边索引；
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)) > 1 || adjSM_eCount.coeffRef(edges(i, 1), edges(i, 0)) > 1)
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();

		int neCount = edgeIdx_nmn.size();						// 非流形有向边个数；
		subFromIdxVec(nmnEdges, edges, edgeIdx_nmn);

		// 4.2 建立边编码-边索引的哈希表；
		for (int i = 0; i < edges.rows(); ++i)
		{
			std::int64_t eCode = encodeEdge(edges(i, 0), edges(i, 1));
			auto iter = edgeMap.find(eCode);
			if (edgeMap.end() == iter)
				edgeMap.insert({ eCode, std::unordered_set<int>{} });
			edgeMap[eCode].insert(i);
		}

		// 4.3 在哈希表中搜索所有非流形边，建立非流形边及其对边的（边——边索引）映射关系；
		for (int i = 0; i < neCount; ++i)
		{
			int neIdx = edgeIdx_nmn[i];		// 当前非流形边索引；
			int vaIdx = nmnEdges(i, 0);
			int vbIdx = nmnEdges(i, 1);
			std::int64_t eCode = encodeEdge(vaIdx, vbIdx);
			std::int64_t oppEcode = encodeEdge(vbIdx, vaIdx);
			auto oppIter = edgeMap.find(oppEcode);

			// 搜索当前非流形有向边的所有索引： 
			edgeIdx_nmn_map.insert({ neIdx, edgeMap[eCode] });

			// 当前边的对边有可能不存在；
			if (edgeMap.end() != oppIter)
				edgeIdx_nmn_opp_map.insert({ neIdx, oppIter->second });
		}
	}

	// 5. 求含有非流形边的三角片的邻接关系：
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



/////////////////////////////////////////////////////////////////////////////////////////////////// 图形编辑：

// 斌杰写的环路光顺接口：
template <typename T>
bool smoothCircuit2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& circuit, const float param = 0.1)
{
	int versCount = circuit.rows();
	Eigen::MatrixXd circDouble(circuit.rows(), circuit.cols());
	circDouble.array() = circuit.array().cast<double>();

	// 1. 环路中心点移动到原点，如果环路太小则适当放大；
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

	// 2. 计算一堆数据；		A = [A1; p * A2];		B = [B1; p * B2];
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

	// 3. 最小二乘法解超定线性方程组A*X = B;
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

	// 4. 环路中心位移到原位，输出；
	if (scale > 0)
		circDouble.array() /= scale;
	circDouble.rowwise() += center;
	circuit.array() = circDouble.array().cast<T>();

	return true;
}

// 网格串联——合并两个孤立的网格到一个网格里
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

// 去除三角片：
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



/////////////////////////////////////////////////////////////////////////////////////////////////// 区域生长算法：

// triangleGrow()——从指定三角片开始区域生长；输入网格不可以有非流形边
template <typename DerivedV, typename DerivedI>
bool triangleGrow(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris, const int triIdx)
{
	/*

			bool  triangleGrow(
						versOut,							// 输出网格的点云
						trisOut,							// 输出网格的三角片
						vers,								// 输入网格的点云
						tris,								// 输入网格的三角片
						triIdx							// 三角片生长的起点三角片索引。
						)

			输入网格若存在非流形边，则会return false
	*/
	int versCount = vers.rows();
	int trisCount = tris.rows();
	int edgesCount = 3 * trisCount;

	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_eIdx的triplet数据；	

	// 1. 求基本的边信息、邻接矩阵：
	Eigen::SparseMatrix<int> adjSM;								// 邻接矩阵；
	Eigen::SparseMatrix<int> adjSM_eCount;					// 邻接矩阵，索引为该有向边重复的次数；
	Eigen::SparseMatrix<int> adjSM_eIdx;				// 邻接矩阵，元素若对应流形边则为该边索引，若对应非流形边则为该边所有索引的和；
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;		// adjSM_eIdx的转置，表示对边的信息；
	Eigen::SparseMatrix<int> adjSM_ueCount;				// 权重为无向边ij关联的三角片数量；
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;				// 非边缘流形无向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB;					// 非边缘流形有向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;			// 非边缘流形有向边对边的邻接矩阵
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

		// 1.1 生成三种有向边邻接矩阵：adjSM, adjSM_eCount, smElems_weighted
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
		adjSM = adjSM_eCount;																		// 有向边邻接矩阵；
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		若adjSM_eIdx(i, j)对应的是流形有向边，该权重值为该边的索引；若是非流形有向边，则该边对应两个边索引，权重为两个索引之和；
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

		//	1.2 生成无向边邻接矩阵：
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// 无向边邻接矩阵；
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.3 生成非边缘流形无向边邻接矩阵：
		smElems.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				if (2 == iter.value())
					smElems.push_back(Eigen::Triplet<int>(iter.row(), iter.col(), 1));
			});
		adjSM_ueMN_NB.resize(versCount, versCount);
		adjSM_ueMN_NB.setFromTriplets(smElems.begin(), smElems.end());
		smElems.clear();

		//	1.4 生成非边缘流形有向边、及其对边的邻接矩阵
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
	std::vector<int> edgesIdx_MN_NB;							// 非边缘流形有向边的索引；
	std::vector<int> edgesIdx_MN_NB_opp;					// 非边缘流形有向边的对边的索引；
	{
		unsigned edgesCount_MN_NB = adjSM_MN_NB.sum();
		edgesIdx_MN_NB.reserve(edgesCount_MN_NB);					// 非边缘流形有向边索引；
		edgesIdx_MN_NB_opp.reserve(edgesCount_MN_NB);			// 非边缘流形有向边的对边的索引；

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB.push_back(adjSM_eIdx.coeffRef(iter.row(), iter.col()));
			});

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB_opp.push_back(adjSM_eIdx_opp.coeffRef(iter.row(), iter.col()));
			});
	}


	// 3. 流形边-三角片邻接关系，三角片邻接关系：
	std::vector<int> etInfo;												// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
	Eigen::VectorXi etAdj_mnEdge;							// 边所在三角片的索引，非流形边或边缘边写-1；
	Eigen::MatrixXi ttAdj_mnEdge;					// 三角片的非边缘流形有向边邻接的三角片索引；	trisCount * 3(i, 0)为索引为i的三角片的非边缘流形有向边(vb, vc)邻接的三角片的索引；
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
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. 循环——队列中第一个三角片t出队 - t写入输出集合 - 求出t所有相邻三角片的索引然后入队
	std::unordered_set<int> finalTriIdx;						// 联通三角片集合；
	std::unordered_set<int> workingSet;						// 用于寻找相邻三角片的队列；
	workingSet.insert(static_cast<int>(triIdx));							// 队列插入第一个三角片；
	while (workingSet.size() > 0)
	{
		// 4.1 首个元素出队，插入到联通三角片集合中；
		int cTriIdx = *workingSet.begin();									// 当前处理的三角片索引
		std::vector<int> adjTriIdx;												// 当前三角片相邻三角片的索引
		workingSet.erase(workingSet.begin());
		finalTriIdx.insert(cTriIdx);

		// 4.2 确定当前三角片的相邻三角片
		for (int i = 0; i < 3; ++i)
			if (ttAdj_mnEdge(cTriIdx, i) >= 0)
				adjTriIdx.push_back(ttAdj_mnEdge(cTriIdx, i));

		// 4.3 若联通三角片集合中没有当前的三角片，插入该三角片，队列中也插入该三角片；
		for (const auto& index : adjTriIdx)
		{
			auto retPair = finalTriIdx.insert(index);
			if (retPair.second)
				workingSet.insert(index);
		}
	}


	// 5. 由选定的三角片取顶点：
	std::set<int> finalVerIdx;
	Eigen::VectorXi finalVerIdxVec;													// 输出顶点的索引，是原网格中的索引；
	Eigen::VectorXi oldNewIdxInfo = -Eigen::VectorXi::Ones(vers.rows());		// 新老索引表——下标为老索引，元素值为新索引，不在新网格中的顶点写为-1
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


	// 6. 原网格所有三角片中的顶点索引由老的换为新的。
	Eigen::MatrixXi tempMat = tris;
	int* intPtr = tempMat.data();
	for (int i = 0; i < tempMat.size(); ++i)
	{
		int oldIdx = *intPtr;
		*intPtr = oldNewIdxInfo[oldIdx];
		intPtr++;
	}

	// 7. 去除非法三角片，即含有-1元素的三角片；
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


// triangleGrowSplitMesh()——区域生长将输入网格分为多个单连通网格；输入网格不可以有非流形边
template <typename DerivedV>
bool triangleGrowSplitMesh(std::vector<DerivedV>& meshesVersOut, std::vector<Eigen::MatrixXi>& meshesTrisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	int versCount = vers.rows();
	int trisCount = tris.rows();
	int edgesCount = 3 * trisCount;

	// 输入参数中的triIdx为包含最高点的三角片索引；
	std::unordered_set<int> finalTriIdx;						// 联通三角片集合；
	std::unordered_set<int> workingSet;						// 用于寻找相邻三角片的队列；

	Eigen::MatrixXi ttAdj_mnEdge;		// 三角片的非边缘流形有向边邻接的三角片索引；	trisCount * 3(i, 0)为索引为i的三角片的非边缘流形有向边(vb, vc)邻接的三角片的索引；

	Eigen::SparseMatrix<int> adjSM;								// 邻接矩阵；
	Eigen::SparseMatrix<int> adjSM_eCount;					// 邻接矩阵，索引为该有向边重复的次数；
	Eigen::SparseMatrix<int> adjSM_eIdx;				// 邻接矩阵，元素若对应流形边则为该边索引，若对应非流形边则为该边所有索引的和；
	Eigen::SparseMatrix<int> adjSM_eIdx_opp;		// adjSM_eIdx的转置，表示对边的信息；
	Eigen::SparseMatrix<int> adjSM_ueCount;				// 权重为无向边ij关联的三角片数量；
	Eigen::SparseMatrix<int> adjSM_ueMN_NB;				// 非边缘流形无向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB;					// 非边缘流形有向边邻接矩阵
	Eigen::SparseMatrix<int> adjSM_MN_NB_opp;			// 非边缘流形有向边对边的邻接矩阵

	Eigen::VectorXi etAdj_mnEdge;							// 边所在三角片的索引，非流形边或边缘边写-1；

	std::vector<int> etInfo;												// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
	std::vector<int> edgesIdx_MN_NB;							// 非边缘流形有向边的索引；
	std::vector<int> edgesIdx_MN_NB_opp;					// 非边缘流形有向边的对边的索引；
	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_eIdx的triplet数据；	

	std::vector<bool> trisFlag(trisCount, true);		// true表示该三角片可收集，false表示已收集；


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

		// 1.1 生成三种有向边邻接矩阵：adjSM, adjSM_eCount, smElems_weighted
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
		adjSM = adjSM_eCount;																		// 有向边邻接矩阵；
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});

		//		若adjSM_eIdx(i, j)对应的是流形有向边，该权重值为该边的索引；若是非流形有向边，则该边对应两个边索引，权重为两个索引之和；
		spMatTranspose(adjSM_eIdx_opp, adjSM_eIdx);

		//	1.2 生成无向边邻接矩阵：
		Eigen::SparseMatrix<int> tmpSm;
		spMatTranspose(tmpSm, adjSM_eCount);
		adjSM_ueCount = adjSM_eCount + tmpSm;					// 无向边邻接矩阵；
		std::vector<int> tmpVec;
		tmpVec.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				tmpVec.push_back(iter.value());
			});

		//	1.3 生成非边缘流形无向边邻接矩阵：
		smElems.reserve(adjSM_ueCount.nonZeros());
		traverseSparseMatrix(adjSM_ueCount, [&](auto& iter)
			{
				if (2 == iter.value())
					smElems.push_back(Eigen::Triplet<int>(iter.row(), iter.col(), 1));
			});
		adjSM_ueMN_NB.resize(versCount, versCount);
		adjSM_ueMN_NB.setFromTriplets(smElems.begin(), smElems.end());
		smElems.clear();

		//	1.4 生成非边缘流形有向边、及其对边的邻接矩阵
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
				edgesIdx_MN_NB.push_back(adjSM_eIdx.coeffRef(iter.row(), iter.col()));
			});

		traverseSparseMatrix(adjSM_MN_NB, [&](auto& iter)
			{
				edgesIdx_MN_NB_opp.push_back(adjSM_eIdx_opp.coeffRef(iter.row(), iter.col()));
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
		ttAdj_mnEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);
	}


	// 4. 循环——队列中第一个三角片t出队 - t写入输出集合 - 求出t所有相邻三角片的索引然后入队
	std::vector<std::unordered_set<int>> connectedTriIdxes;
	int collectedTrisCount = 0;
	while (collectedTrisCount < trisCount)
	{
		finalTriIdx.clear();

		// 取当前首个未收集的三角片作为种子三角片：
		int seedTriIdx = 0;
		for (int i = 0; i < trisCount; ++i)
			if (trisFlag[i])
				seedTriIdx = i;

		workingSet.insert(static_cast<int>(seedTriIdx));					// 队列插入种子三角片；
		while (workingSet.size() > 0)
		{
			// 4.1 首个元素出队，插入到联通三角片集合中；
			int cTriIdx = *workingSet.begin();									// 当前处理的三角片索引
			trisFlag[cTriIdx] = false;

			std::vector<int> adjTriIdx;												// 当前三角片相邻三角片的索引
			workingSet.erase(workingSet.begin());
			finalTriIdx.insert(cTriIdx);

			// 4.2 确定当前三角片的相邻三角片
			for (int i = 0; i < 3; ++i)
				if (ttAdj_mnEdge(cTriIdx, i) >= 0)
					adjTriIdx.push_back(ttAdj_mnEdge(cTriIdx, i));

			// 4.3 若联通三角片集合中没有当前的三角片，插入该三角片，队列中也插入该三角片；
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


	// 5. 输出所有单连通网格
	int meshesCount = connectedTriIdxes.size();
	meshesVersOut.resize(meshesCount);
	meshesTrisOut.resize(meshesCount);
	for (int i = 0; i < meshesCount; ++i)
	{
		DerivedV& cMeshVers = meshesVersOut[i];
		Eigen::MatrixXi& cMeshTris = meshesTrisOut[i];

		// 5.1 由connectedTriIdxes中连通的三角片索引取顶点：
		std::set<int> cMeshVerIdx;																	// 当前的单连通网格中包含的顶点的索引；
		for (const auto& index : connectedTriIdxes[i])
		{
			cMeshVerIdx.insert(tris(static_cast<int>(index), 0));
			cMeshVerIdx.insert(tris(static_cast<int>(index), 1));
			cMeshVerIdx.insert(tris(static_cast<int>(index), 2));
		}

		Eigen::VectorXi cMeshVerIdxVec(cMeshVerIdx.size());									// 输出顶点的索引，是原网格中的索引；
		Eigen::VectorXi oldNewIdxInfo = -Eigen::VectorXi::Ones(vers.rows());			// 新老索引表——下标为老索引，元素值为新索引，不在新网格中的顶点写为-1
		auto iter = cMeshVerIdx.begin();
		for (int i = 0; i < cMeshVerIdx.size(); ++i)
			cMeshVerIdxVec(i) = static_cast<int>(*iter++);

		int setOrder = 0;
		for (const auto& index : cMeshVerIdx)
			oldNewIdxInfo[static_cast<int>(index)] = setOrder++;
		subFromIdxVec(cMeshVers, vers, cMeshVerIdxVec);

		// 5.2 确定当前的单连通网格的三角片：
		std::vector<int> cTriIdxes;
		cTriIdxes.insert(cTriIdxes.end(), connectedTriIdxes[i].begin(), connectedTriIdxes[i].end());
		subFromIdxVec(cMeshTris, tris, cTriIdxes);
		int* verIdxPtr = cMeshTris.data();
		for (int k = 0; k < 3 * cMeshTris.rows(); ++k)				// 三角片中的老的顶点索引改为新的： 
		{
			int oldIdx = *verIdxPtr;
			*verIdxPtr = oldNewIdxInfo[oldIdx];
			verIdxPtr++;
		}
	}

	return true;
}


// triangleGrowOuterSurf()——从指定三角片（必须是外部的三角片）开始区域生长，提取外层表面单连通流形网格
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
	std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;	// 用于生成稀疏矩阵adjSM_eCount, adjSM_eIdx的triplet数据；	

	// 0. 计算所有三角片法向，若存在退化三角片（面积为0），则去除；
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
		if (0 == areaDB)								// 若含有退化三角片，会导致计算面片法向出错；
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

	// 0+.处理异常的非流形边——有些情形下meshArrange会在自交处生成关联三个三角片的无向边，其中一条有向边关联一个，其对边关联两个；
	Eigen::MatrixXi nmnEdges, edges, trisTmp;
	std::vector<std::pair<std::pair<int, int>, int>> nmnInfos;
	std::unordered_map<std::int64_t, int> nmnMap;						// 非流形有向边code-边索引的字典；
	std::unordered_set<std::int64_t> sickCodes;							// 需要预处理的非流形有向边的code；
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
				int triIdx = i % trisCount;							// 当前流形边所在的三角片索引；
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

		// 非流形有向边关联的两个三角片，删除面积较小的那个：
		if (!removeTris(trisTmp, trisClean, sickTriIdxes))
		{
			std::cout << "error!!! removeTris() failed." << std::endl;
			return false;
		}
		trisClean = trisTmp;
	}

	// 1. 求三角片邻接关系：
	Eigen::MatrixXi ttAdj_mnEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency_new(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, trisClean))
	{
		std::cout << "error!!! buildAdjacency() failed." << std::endl;
		return false;
	}

	// 2. 若未指定种子三角片索引，则选取一个不包含非流形边的三角片作为扩散的种子三角片：
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

	// 3. 循环——队列中第一个三角片t出队 - t写入输出集合 - 求出t所有相邻三角片的索引然后入队 
	std::unordered_set<int> finalTriIdxes;
	std::unordered_set<int> workingSet;								// 用于寻找相邻三角片的队列；
	std::vector<int> triTags(trisCount, 0);								// 0 : 未访问，1: 已收录， -1: 已删除；

	// ！！！注意tag写1以后，是有可能被改写为-1的，目前考虑第一次访问到时的标签来决定是否保存三角片；
	workingSet.insert(static_cast<int>(seedIdx));						// 队列插入第一个三角片；
	triTags[seedIdx] = 1;
	while (workingSet.size() > 0)
	{
		// w.1 首个元素出队，插入到联通三角片集合中；
		int cTriIdx = *workingSet.begin();									// 当前处理的三角片索引
		workingSet.erase(workingSet.begin());
		if (triTags[cTriIdx] < 0)											// 有可能入队时是未访问状态，入队之后被标记为删除三角片；
			continue;
		finalTriIdxes.insert(cTriIdx);

		// w.2 确定当前三角片相邻的所有三角片，根据情况保留和删除
		std::vector<int> adjTriIdx;												// 当前三角片相邻三角片的索引
		for (int i = 0; i < 3; ++i)
		{
			int nmAdjTriIdx = ttAdj_mnEdge(cTriIdx, i);
			if (nmAdjTriIdx >= 0)										// 流形边（非边缘）；
				adjTriIdx.push_back(nmAdjTriIdx);						// wf.a. 流形边关联的相邻三角片都保留；
			else if (nmAdjTriIdx == -1)								// 非流形边
			{
				// wf.b.1. 删除当前非流形边关联的其他三角片：
				auto& vec1 = ttAdj_nmnEdge[cTriIdx][i];				// 当前非流形边所在的三角片，包括当前三角片自身；
				for (const auto& index : vec1)
				{
					if (index == cTriIdx)
						continue;
					else
						triTags[index] = -1;
				}

				// wf.b.2. 当前非流形边的对边所在的三角片中，选取法向夹角最大(即dot值最小)的那个；
				auto& vec2 = ttAdj_nmnOppEdge[cTriIdx][i];			// 当前非流形边的对边所在的所有三角片
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

					for (const auto& index : vec2)				// 未被选取的对面其他三角片全部删除 ； 
						if (selTriIdx != index)
							triTags[index] = -1;
				}

			}

			// nmAdjTriIdx == -2边缘边，直接跳过；
		}

		// w.3 finalTriIdxes中插入保留下来的相邻三角片；若该三角片之前未访问，则插入队列；
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

	// 4. 提取选中的三角片
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vers1 = vers;
	Eigen::MatrixXi tris1;
	subFromIdxCon(tris1, trisClean, finalTriIdxes);

	// 5. 去除孤立顶点；
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


// robustTriangleGrowOuterSurf()——提升了鲁棒性的三角片生长，避免了某些反向三角片的干扰；
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

	// 有时有反向三角片的干扰，导致生长出来的网格只是原网格很小的一部分，需要尝试多次避免这种错误；
	for (int i = 0; i < maxIter; ++i)
	{
		versOut.resize(0, 0);
		trisOut.resize(0, 0);
		startIdx = UID_i(e);
		triangleGrowOuterSurf(versOut, trisOut, vers, tris, blRemIso, startIdx);
		int trisCountNew = trisOut.rows();
		if (trisCountNew > trisCount / 2)						// 输出三角片过少则说明有错误，重来；
		{
			subFlag = true;
			break;
		}
	}

	return subFlag;
}



// 无向边的区域生长，提取环路边和边曲线；输入的边矩阵可以是有向边也可以是无向边；
template <typename DerivedI>
bool uEdgeGrow(std::vector<Eigen::MatrixXi>& circs, std::vector<Eigen::MatrixXi >& segs, \
	const Eigen::PlainObjectBase<DerivedI>& uEdges)
{
	/*
	bool uEdgeGrow(
			std::vector<Eigen::MatrixXi>& circs,									元素是组成一个环路的边缘边数据
			std::vector<Eigen::MatrixXi >& segs,									元素是组成一个非封闭曲线的边缘边数据
			const Eigen::PlainObjectBase<DerivedI>& uEdges			ueCount * 2, 无向边数据
			)

	*/
	std::list<std::list<std::pair<int, int>>> circList, segList;
	Eigen::SparseMatrix<int> adjSM;										// 邻接矩阵；
	std::vector<Eigen::Triplet<int>> smElems;

	// 1. 建立无向边邻接矩阵：
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

	// 2. 统计不重复的无向边数量；
	std::unordered_set<std::int64_t> tmpSet;
	for (unsigned i = 0; i < uEdges.rows(); ++i)
		tmpSet.insert(encodeUedge(uEdges(i, 0), uEdges(i, 1)));
	int remainCount = tmpSet.size();			// (3,4) (4,3)表示同一条无向边；

	// 3. 无向边生长的循环：
	while (remainCount > 0)
	{
		std::list<std::pair<int, int>> currentSeg;

		// w1. 提取邻接表中第一个无向边：
		bool breakFlag = false;
		for (unsigned i = 0; !breakFlag & i < adjSM.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, i); !breakFlag & iter; ++iter)
			{
				if (iter.valueRef() > 0)
				{
					currentSeg.push_back({ iter.row(), iter.col() });
					iter.valueRef() = -1;												// 表示删除该元素；
					adjSM.coeffRef(iter.col(), iter.row()) = -1;
					remainCount--;
					breakFlag = true;
				}
			}
		}

		// w2. 搜索第一条无向边关联的边，直到搜索不到，或者形成环路为止。
		int head = currentSeg.front().first;
		int tail = currentSeg.front().second;
		while (1)
		{
			int otherEnd = -1;
			breakFlag = false;
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, head); !breakFlag & iter; ++iter)	// 对第head列的遍历：
			{
				if (iter.valueRef() > 0)
				{
					otherEnd = iter.row();
					currentSeg.push_back({ head, otherEnd });
					iter.valueRef() = -1;												// 表示删除该元素；
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

	// 4. 输出：
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

// 确定三角网格中的所有单连通区域（顶点连通）——使用稀疏矩阵；
template <typename Index>
int simplyConnectedRegion(const Eigen::SparseMatrix<Index>& adjSM,
	Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount)
{
	/*
		int simplyConnectedRegion(												返回单连通区域的数量，出错时返回-1
					const Eigen::SparseMatrix<Index>& adjSM,			三角网格的邻接矩阵
					Eigen::VectorXi& connectedLabels,						同一连通区域内的顶点标签，依次为0, 1, 2,...；
																											索引为i的顶点的标签为connectedLabels(i);
					Eigen::VectorXi& connectedCount							每个标签对应的单连通区域内包含的顶点数
																											标签为i的单连通区域包含的顶点数为connectedCount(i)
					)
	*/
	if (adjSM.rows() == 0)
		return -1;

	if (adjSM.rows() != adjSM.cols())
		return -1;

	const unsigned versCount = adjSM.rows();

	// 1. 初始化：
	connectedLabels.setConstant(versCount, versCount);		// 最终结果中最大可能的标签为versCount-1;
	connectedCount.setZero(versCount);

	// 2. 遍历所有顶点，搜索与其联通的顶点：
	int currentLabel = 0;
	for (unsigned i = 0; i < versCount; ++i)
	{
		// 2.1 若当前顶点已被访问过(label < versCount)，则continue:
		if (connectedLabels(i) < versCount)
			continue;

		// 2.2 若当前顶点未被访问过，执行BFS收集其所有联通的顶点：
		std::queue<Index> workingQueue;
		workingQueue.push(i);
		while (!workingQueue.empty())
		{
			const Index curVerIdx = workingQueue.front();
			workingQueue.pop();

			if (connectedLabels(curVerIdx) < versCount)
				continue;

			// 2.2.1 队首顶点label赋值，label计数+1
			connectedLabels(curVerIdx) = currentLabel;
			connectedCount(currentLabel)++;

			// 2.2.2 在邻接矩阵中搜索当前顶点邻接的顶点：
			for (typename Eigen::SparseMatrix<Index>::InnerIterator iter(adjSM, curVerIdx); iter; ++iter)
			{
				const Index connectVerIdx = iter.row();					// 默认稀疏矩阵列优先存储；当前迭代器在第curVerIdx列；

				if (connectedLabels(connectVerIdx) < versCount)
					continue;

				workingQueue.push(connectVerIdx);
			}
		}

		// 2.3 上一个标签的顶点收集完毕，下一个循环收集下一个标签的顶点：
		currentLabel++;
	}

	// 3. shrink_to_fit()
	connectedCount.conservativeResize(currentLabel, 1);

	return currentLabel;
}


// 确定三角网格中的所有单连通区域（顶点连通）——使用std::unordered_set，不使用稀疏矩阵；
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
		} while (blCollectedNew);						// 本轮有新元素被收集，则继续执行一轮

		return versSet.size();										// 返回收集的顶点个数，即被标记为label的顶点个数；
	};

	if (vers.rows() == 0 || tris.rows() == 0)
		return -1;

	const int versCount = vers.rows();
	const int trisCount = tris.rows();

	// 1. 区域生长的循环：
	connectedLabels.resize(versCount, versCount);		// 最终结果中最大可能的标签为versCount-1;
	std::unordered_set<int> versSet;
	std::vector<bool> triCollected(trisCount, false);
	int currentLabel = 0;
	int collectedVersCount = 0;
	while (collectedVersCount < versCount)
	{
		// w1. 选取第一个未被标记的三角片，将其三个顶点作为区域生长的种子顶点；
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

		// w2. 区域生长的循环：
		int currentLabelVersCount = collectVers(connectedLabels, triCollected, versSet, tris, currentLabel);

		// w3. post procedure:
		versSet.clear();
		collectedVersCount += currentLabelVersCount;
		currentLabel++;
	}

	// 2. 统计：
	int scCount = currentLabel;
	connectedCount.resize(scCount, 0);
	for (int i = 0; i < versCount; ++i)
	{
		int label = connectedLabels[i];
		connectedCount[label]++;
	}

	return scCount;
}



// 确定三角网格中的所有单连通区域（三角片连通）——比triangleGrow更强，输入网格可以带有非流形元素；
template <typename DerivedI>
int simplyTrisConnectedRegion(Eigen::VectorXi& connectedLabels, Eigen::VectorXi& connectedCount, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
	int simplyTrisConnectedRegion(											返回三角片单连通区域的数量，出错时返回-1
				Eigen::VectorXi& connectedLabels,						同一连通区域内的三角片标签，依次为0, 1, 2,...；
																										索引为i的三角片的标签为connectedLabels(i);
				Eigen::VectorXi& connectedCount							每个标签对应的单连通区域内包含的三角片数
																										标签为i的单连通区域包含的三角片数为connectedCount(i)
				const Eigen::PlainObjectBase<DerivedI>& tris
				)
*/
	const unsigned trisCount = tris.rows();
	connectedLabels = Eigen::VectorXi::Zero(trisCount);
	unsigned visitedCount = 0;

	// 1. 计算网格三角片邻接关系：
	Eigen::MatrixXi ttAdj_mnEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency_new(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return -1;

	// 2. 邻接三角片生长的循环：
	int currentLabel = 1;
	while (visitedCount < trisCount)
	{
		// w1. 选取种子：
		int seedIdx = 0;
		for (unsigned i = 0; i < trisCount; ++i)
		{
			if (0 == connectedLabels(i))
			{
				seedIdx = i;
				break;
			}
		}

		// w2. 开始生长：
		std::unordered_set<int> workingQueue;
		workingQueue.insert(seedIdx);
		while (!workingQueue.empty())
		{
			// ww1. 队首
			int cIdx = *workingQueue.begin();
			visitedCount++;
			workingQueue.erase(cIdx);
			connectedLabels(cIdx) = currentLabel;

			// ww2. 搜索当前三角片的流形有向边相对的三角片
			for (int i = 0; i < 3; ++i)
				if (ttAdj_mnEdge(cIdx, i) >= 0)
					if (0 == connectedLabels(ttAdj_mnEdge(cIdx, i)))			// 标签为0表示该三角片未被访问过；
						workingQueue.insert(ttAdj_mnEdge(cIdx, i));					// 只插入未访问的三角片 

			// ww3. 搜索当前三角片的非流形有向边相对的三角片：
			for (int i = 0; i < 3; ++i)
				for (const auto& index : ttAdj_nmnOppEdge[cIdx][i])
					if (0 == index)
						workingQueue.insert(index);
		}

		// w3. 当前联通区域生长完成
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



// 提取三角网格中最大单连通区域（顶点连通）
template <typename DerivedV, typename DerivedI>
bool simplyConnectedLargest(Eigen::PlainObjectBase<DerivedV>& versOut, Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();

	// 1. 生成邻接矩阵：
	Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx, adjSM;
	adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
	adjSM = adjSM_eCount;
	traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
		{
			if (iter.value() > 0)
				iter.valueRef() = 1;
		});

	// 2. 确定网格中所有单连通区域：
	Eigen::VectorXi connectedLabels, connectedCount;
	int conCount = simplyConnectedRegion(adjSM, connectedLabels, connectedCount);
	if (conCount < 0)
		return false;

	// 3. 确定最大单连通区域（顶点数最多）；
	int mainLabel = 0;								// 最大单连通区域的标签；
	int mainLabelCount = 0;
	for (int i = 0; i < conCount; ++i)
	{
		if (connectedCount(i) > mainLabelCount)
		{
			mainLabel = i;
			mainLabelCount = connectedCount(i);
		}
	}

	// 4. 提取最大单连通区域的顶点：
	std::vector<int> mainLabelIdxes;
	mainLabelIdxes.reserve(versCount);
	for (unsigned i = 0; i < versCount; ++i)
		if (mainLabel == connectedLabels(i))
			mainLabelIdxes.push_back(i);
	mainLabelIdxes.shrink_to_fit();
	subFromIdxVec(versOut, vers, mainLabelIdxes);


	// 5. 提取最大单连通区域的三角片：

	//		5.1 生成老-新索引映射表
	std::vector<int> oldNewIdxInfo(versCount, -1);
	for (int i = 0; i < mainLabelIdxes.size(); ++i)
	{
		int oldIdx = mainLabelIdxes[i];
		oldNewIdxInfo[oldIdx] = i;
	}

	//		5.2 三角片数据中的顶点索引映射成新索引；
	DerivedI trisCopy = tris;
	int* intPtr = trisCopy.data();
	for (int i = 0; i < trisCopy.size(); ++i)
	{
		int oldIdx = *intPtr;
		*intPtr = oldNewIdxInfo[oldIdx];
		intPtr++;
	}

	//		5.3 提取最大联通区域内的三角片：
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


// 按顶点单连通区域分解网格：
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

	// 0. 预处理——去除孤立顶点，去除非法、重复三角片：
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

	// 1. 提取顶点单连通点云：
	const int versCount = vers0.rows();
	int scCount = 0;							// 顶点单连通区域个数；
	int sctCount = 0;							// 三角片单连通区域个数；
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

	// 2. 生成顶点单连通网格：
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
			*ptrData = oldNewIdxInfo[*ptrData];					// 老索引映射成新索引；
			ptrData++;
		}
		removeSickDupTris(compVers[i], compTris[i]);		// 去除非法三角片；
	}

	return true;
}



/////////////////////////////////////////////////////////////////////////////////////////////////// 网格缺陷检查和修复


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


// findOverlapTris()——寻找网格中的重叠三角片，基于三角片区域生长的方法；
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
	std::unordered_set<std::int64_t> opTriPairCodes;			// 一对重叠的三角片，索引值按前小后大排列，编码成std::int64_t

#ifdef  LOCAL_DEBUG
	std::set<int> olTriIdxes;
#endif 

	// 0. 求初始状态下网格所有三角片的法向：
	Eigen::MatrixXd triNorms;
	if (!getTriNorms(triNorms, vers, tris))				// 若含有退化三角片，会导致计算面片法向出错；
		return false;

	// 1. 求三角片邻接关系：
	Eigen::MatrixXi ttAdj_mnEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return false;

	// 2. 三角片区域生长-寻找相邻的重叠三角片对的循环；
	std::unordered_set<int> workingSet;								// 用于寻找相邻三角片的队列；
	std::vector<int> triTags(trisCount, 0);								// 0 : 未访问，1: 已访问
	workingSet.insert(static_cast<int>(triIdx));						// 队列插入第一个三角片； 
	while (workingSet.size() > 0)
	{
		// 2.1 首个元素出队，插入到联通三角片集合中；
		int cTriIdx = *workingSet.begin();									// 当前处理的三角片索引
		std::vector<int> adjTriIdxes;												// 当前三角片的未被访问的相邻三角片的索引
		Eigen::RowVector3d cTriNorm = triNorms.row(cTriIdx);
		unsigned adjTrisCount = 0;
		workingSet.erase(workingSet.begin());
		triTags[cTriIdx] = 1;											// 当前三角片标记为已访问；

		// 2.2 确定当前三角片的所有相邻三角片
		for (int i = 0; i < 3; ++i)
		{
			int nbrTriIdx = ttAdj_mnEdge(cTriIdx, i);
			if (nbrTriIdx >= 0)
			{
				// a. 当前边为流形边
				if (0 == triTags[nbrTriIdx])
					adjTriIdxes.push_back(nbrTriIdx);
			}
			else
			{
				// b. 当前边为非流形边；
				auto& vec1 = ttAdj_nmnEdge[cTriIdx][i];						// 当前非流形边所在的所有三角片；
				auto& vec2 = ttAdj_nmnOppEdge[cTriIdx][i];				// 当前非流形边的对边所在的所有三角片；
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

		// 2.3 检验-判定当前三角片和周边三角片是否重叠；
		adjTrisCount = adjTriIdxes.size();
		if (0 == adjTrisCount)
			continue;
		for (unsigned k = 0; k < adjTrisCount; ++k)
		{
			int adjIdx = adjTriIdxes[k];
			double dotValue = cTriNorm.dot(triNorms.row(adjIdx));
			if (dotValue > dotThreshold2 || dotValue < dotThreshold1)
			{
				// 若当前三角片和相邻三角片法向接近相同或相反，检测两个三角片投影到其法向平面上是否有重叠区域 
				int A, B, C1, C2;				// AB是这对相邻三角片共享的顶点，C1, C2是它们各自相异的顶点；
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

				// 判断C1, C2是否在同侧：
				Eigen::RowVector3d AB = vers.row(B).array().cast<double>() - vers.row(A).array().cast<double>();
				Eigen::RowVector3d AC1 = vers.row(C1).array().cast<double>() - vers.row(A).array().cast<double>();
				Eigen::RowVector3d AC2 = vers.row(C2).array().cast<double>() - vers.row(A).array().cast<double>();
				Eigen::RowVector3d crossArrow1 = AB.cross(AC1);
				Eigen::RowVector3d crossArrow2 = AB.cross(AC2);
				if (crossArrow1.dot(cTriNorm) * crossArrow2.dot(cTriNorm) > 0)				// C1, C2在公共边AB的同一侧；
				{
					auto retPair = opTriPairCodes.insert(encodeUedge(cTriIdx, adjIdx));
					if (retPair.second)
					{
						opTrisPairs.push_back(std::make_pair(cTriIdx, adjIdx));
						olCount++;								// 计数增加；
#ifdef LOCAL_DEBUG
						olTriIdxes.insert(cTriIdx);
						olTriIdxes.insert(adjIdx);
#endif
					}
				}
			}
		}

		// 2.3 工作队列中插入当前三角片未被访问的相邻三角片；
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



// 求三角网格中的非流形有向边；若某条有向边重复次数大于1，则该边及其对边都被判定为非流形边；
template<typename DerivedI>
int nonManifoldEdges(Eigen::MatrixXi& nmnEdges, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	nmnEdges.resize(0, 0);
	Eigen::SparseMatrix<int> adjSM_eIdx, adjSM_eCount;
	if (!adjMatrix(adjSM_eCount, adjSM_eIdx, tris))
		return false;

	// 遍历adjSM_eCount，搜索非流形边；
	std::unordered_set<std::int64_t> nmnEdgeIdxes;
	traverseSparseMatrix(adjSM_eCount, [&adjSM_eCount, &nmnEdgeIdxes](auto& iter)
		{
			if (iter.value() > 1)
			{
				std::int64_t code = encodeEdge(iter.row(), iter.col());
				nmnEdgeIdxes.insert(code);
				if (adjSM_eCount.coeff(iter.col(), iter.row()) > 0)						// 当前边的对边有可能不存在；
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


// 求三角网格中的非流形有向边；若某条有向边重复次数大于1，则该边及其对边都被判定为非流形边；
template<typename DerivedI>
int nonManifoldEdges(Eigen::MatrixXi& nmnEdges, std::vector<std::pair<std::pair<int, int>, int>>& nmnInfos, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	nmnEdges.resize(0, 0);
	Eigen::SparseMatrix<int> adjSM_eIdx, adjSM_eCount;
	if (!adjMatrix(adjSM_eCount, adjSM_eIdx, tris))
		return false;

	// 遍历adjSM_eCount，搜索非流形边；
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
				if (oppEdgeCount > 0)						// 当前边的对边有可能不存在；
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

	// 生成非流形边信息：
	nmnInfos.reserve(nmnCount);
	for (const auto& pair : repCountMap)
	{
		auto edgePair = decodeEdge(pair.first);
		nmnInfos.push_back(std::make_pair(std::make_pair(edgePair.first, edgePair.second), pair.second));
	}

	return nmnCount;
}


// 求三角网格中的非流形无向边，重载1 
template <typename T>
int nonManifoldUEs(Eigen::MatrixXi& nmnUedges, const Eigen::SparseMatrix<T>& adjSM_eCount, \
	const Eigen::SparseMatrix<T>& adjSM_ueCount)
{
	nmnUedges.resize(0, 0);
	std::unordered_set<std::int64_t> nmnUeCodes;

	// 遍历adjSM_ueCount，搜索非流形边；
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


// 求三角网格中的非流形无向边，重载2 
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


// 求三角网格中的非流形无向边，重载3——输出非流形无向边的详细信息：
template<typename DerivedI>
int nonManifoldUEs(Eigen::MatrixXi& nmnUedges, std::vector<std::tuple<std::pair<int, int>, int, int >>& nmnEdgeInfos, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		int nonManifoldUEs(																			返回检测到的非流形无向边的数量，失败返回-1；
				Eigen::MatrixXi& nmnUedges,																nmnCount * 2, 非流形无向边数据
				std::vector<std::pair<int, std::pair<int, int>>>& nmnEdgeInfos,			元素的first是有向边重复的次数，second是这条非流形有向边的具体数据；
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


// 检测非流形点——输入网格需要没有重复三角片和非法三角片
template<typename DerivedV>
int nonManifoldVers(std::vector<int>& nmnVerIdxes, const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	nmnVerIdxes.clear();
	const int versCount = vers.rows();
	const int trisCount = tris.rows();
	std::unordered_multimap<int, std::int64_t> triMap;

	// 1. 计算顶点所对的所有无向边编码，存入字典
	for (int i = 0; i < trisCount; ++i)
	{
		const int& vaIdx = tris(i, 0);
		const int& vbIdx = tris(i, 1);
		const int& vcIdx = tris(i, 2);
		triMap.insert({ vaIdx, encodeUedge(vbIdx, vcIdx) });
		triMap.insert({ vbIdx, encodeUedge(vcIdx, vaIdx) });
		triMap.insert({ vcIdx, encodeUedge(vaIdx, vbIdx) });
	}

	// 2. 遍历所有顶点，分析其所对的所有无向边：
	nmnVerIdxes.reserve(versCount);
	for (int i = 0; i < versCount; ++i)
	{
		// f1. 当前顶点所对的所有无向边存入双向链表：
		std::list<std::pair<int, int>> roundUes;
		auto iter = triMap.find(i);
		if (iter == triMap.end())
			continue;
		for (int k = 0; k < triMap.count(i); ++k)
			roundUes.push_back(decodeEdge((iter++)->second));

		// f2. 无向边生长的循环：
		std::unordered_set<int> connectedVers;
		std::pair<int, int> currentUe = roundUes.front();
		roundUes.pop_front();
		while (!roundUes.empty())
		{
			bool blFound = false;
			auto iter = roundUes.begin();
			connectedVers.insert(currentUe.first);
			connectedVers.insert(currentUe.second);

			// fw1. 若当前无向边和已连通无向边存在相同的顶点，则两者是连通的；
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

			// fw1. 若当前连通边已生长完全，但roundUes中仍有别的边，说明此顶点关联两个或者更多的扇形，当前顶点是非流形点；
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



// correctTriDirs()——矫正网格的三角片朝向，基于三角片区域生长的方法；使用前需要先处理网格中的重叠三角片；
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
	std::vector<int> corrTriIdxes;							// 存储方向错误的三角片的索引；
#endif 

	// 0. 求初始状态下网格所有三角片的法向：
	Eigen::MatrixXd triNorms;
	if (!getTriNorms(triNorms, vers, tris))				// 若含有退化三角片，会导致计算面片法向出错；
		return false;

	// 1. 求三角片邻接关系：
	Eigen::MatrixXi ttAdj_mnEdge;
	std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
	if (!buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris))
		return false;

	// 2. 循环——队列中第一个三角片t出队 - t写入输出集合 - 求出t所有相邻三角片的索引然后入队
	int triIdx = 0;
	std::unordered_set<int> workingSet;								// 用于寻找相邻三角片的队列；
	std::vector<int> triTags(trisCount, 0);								// 0 : 未访问，1: 已访问，-1: 已访问，被标记为朝向错误三角片；
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

		workingSet.insert(static_cast<int>(triIdx));						// 队列插入第一个三角片； 
		while (workingSet.size() > 0)
		{
			// w1. 首个元素出队，插入到联通三角片集合中；
			int cTriIdx = *workingSet.begin();									// 当前处理的三角片索引
			workingSet.erase(workingSet.begin());
			if (triTags[cTriIdx] != 0)
				continue;
			triTags[cTriIdx] = 1;													// 当前三角片标记为已访问；
			visitedCount++;

			std::vector<int> adjTriIdxes;									// 当前三角片的未被访问的相邻三角片的索引
			Eigen::RowVector3d cTriNorm = triNorms.row(cTriIdx);
			Eigen::MatrixXd adjTriNorms;
			unsigned adjTrisCount = 0;

			// w2. 确定当前三角片的所有相邻三角片
			for (int i = 0; i < 3; ++i)
			{
				int nbrTriIdx = ttAdj_mnEdge(cTriIdx, i);
				if (nbrTriIdx >= 0)
					if (0 == triTags[nbrTriIdx])								// a. 当前边为流形边
						adjTriIdxes.push_back(nbrTriIdx);
					else
					{
						// b. 当前边为非流形边；
						auto& vec1 = ttAdj_nmnEdge[cTriIdx][i];						// 当前非流形边所在的所有三角片；
						auto& vec2 = ttAdj_nmnOppEdge[cTriIdx][i];				// 当前非流形边的对边所在的所有三角片；
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

			// 4.3 检验-矫正当前三角片的相邻三角片：
			adjTrisCount = adjTriIdxes.size();
			if (0 == adjTrisCount)
				continue;
			subFromIdxVec(adjTriNorms, triNorms, adjTriIdxes);
			for (unsigned k = 0; k < adjTrisCount; ++k)
			{
				int adjIdx = adjTriIdxes[k];
				if (cTriNorm.dot(adjTriNorms.row(k)) < dotThreshold)
				{
					// 三角片标记为-1：
					triTags[adjIdx] = -1;
					visitedCount++;
					corrCount++;												// 计数增加；
#ifdef LOCAL_DEBUG
					corrTriIdxes.push_back(adjIdx);
#endif

					////		矫正三角片朝向
					//int tmp = trisOut(adjIdx, 1);
					//trisOut(adjIdx, 1) = trisOut(adjIdx, 2);
					//trisOut(adjIdx, 2) = tmp;
					//triNorms.row(adjIdx).array() *= -1;			//	更新法向：
				}
			}

			// 4.3 工作队列中插入当前三角片未被访问的相邻三角片；
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

// 补洞——不插点，直接对洞进行三角剖分： 
template <typename IndexType>
bool fillSmallHoles(Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic>& newTris, \
	const std::vector<std::vector<int>>& holes)
{
	/*
	bool fillSmallHoles(
		Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic>& newTris,				// 补洞生成的新triangle soup
		const std::vector<std::vector<int>>& holes					// 洞的信息，每个元素是有序的环路顶点索引向量；
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


// 检测网格中的孤立顶点
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


// 去除网格中的孤立顶点：
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

	// 输出点云：
	subFromIdxVec(versOut, vers, newOldIdxInfo);

	// 输出三角片：
	trisOut = tris;
	int* dataPtr = trisOut.data();
	for (int i = 0; i < trisCount * 3; ++i)
		*dataPtr++ = oldNewIdxInfo[*dataPtr];

	return true;
}


// 检测网格中是否有非法三角片（三个顶点索引中有两个相同的）
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




// 检测退化边（边长过短的边）
template <typename DerivedV>
Eigen::VectorXi checkDegEdges(const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const double threshold = 1e-3)
{
	/*
		Eigen::VectorXi checkDegEdges(														// 返回标记退化边的flag向量
				const Eigen::MatrixXi& edges,												// 网格有向边
				const Eigen::PlainObjectBase<DerivedV>& edgeArrows,		// 有向边的边向量数据
				const Eigen::PlainObjectBase<DerivedV>& vers,
				const Eigen::MatrixXi& tris,
				const double threshold = 1e-3												// 判定退化边的边长阈值；
				)
	*/
	Eigen::VectorXd edgesLen = edgeArrows.rowwise().norm();
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	unsigned edgesCount = edges.rows();

	Eigen::VectorXi flags = (edgesLen.array() <= threshold).select(Eigen::VectorXi::Ones(edgesCount), Eigen::VectorXi::Zero(edgesCount));

	return flags;
}


// 去除网格退化边（边长过短的边，将两个顶点融合）;
template <typename DerivedV>
int mergeDegEdges(Eigen::PlainObjectBase<DerivedV>& newVers, Eigen::MatrixXi& newTris, \
	const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& edgeArrows, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris, const Eigen::VectorXi& degEdgeFlags)
{
	/*
		int mergeDegEdges(																				返回去除掉的顶点数，若失败则返回-1；
				Eigen::PlainObjectBase<DerivedV>& newVers,
				Eigen::MatrixXi& newTris,
				const Eigen::MatrixXi& edges,
				const Eigen::PlainObjectBase<DerivedV>& edgeArrows,
				const Eigen::PlainObjectBase<DerivedV>& vers,
				const Eigen::MatrixXi& tris,
				const Eigen::VectorXi& degEdgeFlags								标记退化边的flag向量
				)

	*/
	int repVersCount = -1;
	int* dataPtr = nullptr;
	unsigned versCount = vers.rows();
	unsigned trisCount = tris.rows();
	unsigned edgesCount = edges.rows();

	// 1. 提取退化有向边
	int degCount = degEdgeFlags.sum();				// 退化有向边的数量
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

	// 所有退化边中，保留索引较小的顶点，索引大的顶点改写为索引小的顶点：
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
			else       // 当前边的顶点与之前的簇中的顶点相同，插入到该簇中；
			{
				clusterSet.insert(vaIdx);
				clusterSet.insert(vbIdx);
				break;
			}
		}

		// 生成新的簇；
		std::set<int> tmpCluster;
		tmpCluster.insert(vaIdx);
		tmpCluster.insert(vbIdx);
		clusters.push_back(tmpCluster);
	}

	// 生成顶点聚类的映射字典——一个簇中所有其他顶点都映射为索引最小的那个顶点：
	std::map<int, int> clusterMap;
	for (const auto clusterSet : clusters)
	{
		auto iter = clusterSet.begin();
		int headIdx = *iter;
		iter++;
		for (; iter != clusterSet.end(); ++iter)
			clusterMap.insert({ *iter, headIdx });
	}

	std::vector<int> repIdxes;							// 需要去除的顶点的索引；
	repIdxes.reserve(degCount);
	for (const auto& pair : clusterMap)
		repIdxes.push_back(pair.second);
	repIdxes.shrink_to_fit();
	repVersCount = repIdxes.size();

	std::map<int, int> fullMap = clusterMap;
	for (int i = 0; i < versCount; ++i)
		fullMap.insert({ i, i });

	// 聚类映射：
	Eigen::MatrixXi trisCopy = tris;
	dataPtr = trisCopy.data();
	for (unsigned i = 0; i < 3 * trisCount; ++i)
	{
		int index = *dataPtr;
		*dataPtr = fullMap[index];
		dataPtr++;
	}

	// 删除非法三角片：
	std::vector<int> sickIdxes;
	checkSickTris(sickIdxes, trisCopy);
	removeTris(newTris, trisCopy, sickIdxes);

	// 删除孤立顶点：
	trisCopy = newTris;
	std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, trisCopy);
	newVers.resize(0, 0);
	newTris.resize(0, 0);
	removeIsoVers(newVers, newTris, vers, trisCopy, isoVerIdxes);

	return degCount;
}



// 去除非法三角片，重复三角片：
template<typename T>
int removeSickDupTris(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris)
{
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	int versCount = vers.rows();
	int trisCount = tris.rows();

	Eigen::VectorXi flags{ Eigen::VectorXi::Ones(trisCount) };			// flag向量，需要去除的三角片标0；
	std::unordered_set<std::uint64_t> triCodeSet;

	// 遍历三角片
	for (int i = 0; i < trisCount; ++i)
	{
		// 三角片中的索引值不在0~versCount-1之间则为非法三角片；貌似有些STL文件中会出现这种情形；
		if (tris(i, 0) < 0 || tris(i, 0) > versCount - 1 || tris(i, 1) < 0 || tris(i, 1) > versCount - 1 || tris(i, 2) < 0 || tris(i, 2) > versCount - 1)
		{
			flags(i) = 0;
			continue;
		}

		std::uint64_t code = encodeTriangle(tris(i, 0), tris(i, 1), tris(i, 2));
		if (0 == code)
		{
			flags(i) = 0;				// 三角片的三个索引中有相同的，非法；
			continue;
		}

		auto retPair = triCodeSet.insert(code);
		if (!retPair.second)
			flags(i) = 0;				// 重复三角片非法；
	}

	Eigen::MatrixXi tmpMat;
	subFromFlagVec(tmpMat, tris, flags);
	tris = tmpMat;

	int removeCount = trisCount - tris.rows();
	return removeCount;
}


// 搜索网格中组成边缘回路(即洞的边缘)——不考虑洞共边的情形。！！！当前对于太复杂的洞，返回的洞顶点索引顺序可能不对。效果不如VBFindHole2
template <typename DerivedV, typename DerivedI>
int findHoles(std::vector<std::vector<int>>& holes, const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::PlainObjectBase<DerivedI>& tris)
{
	// 注: 输入网格不可以存在非法三角片、重复三角片。
	/*
		int findHoles(																			成功返回洞数，失败返回-1;
			std::vector<std::vector<int>>& holes,								洞，即闭合为环路的边缘顶点，每圈顶点按顺时针排列；
			const Eigen::PlainObjectBase<DerivedV>& vers,
			const Eigen::PlainObjectBase<DerivedI>& tris
			)
	*/
	// lambda——整理洞的信息：
	auto organizeHoles = [](std::vector<std::vector<int>>& holes)		// 整理每个洞的顶点顺序，使得首个顶点是索引最小的那个
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

	// 确定洞的边——只关联一个三角片的边：
	const int  trisCount = tris.rows();
	const int  edgesCount = 3 * trisCount;
	const int  versCount = tris.maxCoeff() + 1;
	holes.clear();

	// 1. 提取所有边缘有向边：
	std::vector<Eigen::MatrixXi> circs;
	Eigen::MatrixXi bdrys;
	bdryEdges(bdrys, tris);
	const int bdrysCount = bdrys.rows();

	// 2. 生成边缘边链表
	std::list<std::pair<int, int>> bdryList;
	for (int i = 0; i < bdrysCount; ++i)
		bdryList.push_back({ bdrys(i, 0), bdrys(i, 1) });

	// 3. 有向边生长，寻找所有的边缘回路；
	std::list<std::list<int>> holeList;
	std::pair<int, int> currentEdge = bdryList.front();
	int currentHead = currentEdge.first;								// 当前洞的首顶点索引；
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

	// 4. 整理，输出：
	const int holesCount = holeList.size();
	holes.reserve(holesCount);
	for (const auto& list : holeList)
	{
		const int holeVersCount = list.size() - 1;
		auto iterList = list.rbegin();									// 使用反向迭代器，按照逆序填充；
		std::vector<int> currentHole(holeVersCount);
		if (holeVersCount < 3)
			return -1;							// 没有形成回路；
		if (*list.begin() != *list.rbegin())
			return -1;							// 没有形成回路；
		for (int i = 0; i < holeVersCount; ++i)
			currentHole[i] = *iterList++;
		holes.push_back(currentHole);
	}
	organizeHoles(holes);

	return holesCount;
}





