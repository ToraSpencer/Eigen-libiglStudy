
/////////////////////////////////////////////////////////////////////////////////////////////////// general tools: 

// 传入函数子或函数指针遍历stl容器
template<typename T, typename F>
void traverseSTL(T& con, F f)
{
	std::for_each(con.begin(), con.end(), f);
	std::cout << std::endl;
}


// 反向遍历
template<typename T, typename F>
void revTraverseSTL(T& con, F f)
{
	std::for_each(con.rbegin(), con.rend(), f);
	std::cout << std::endl;
}



/////////////////////////////////////////////////////////////////////////////////////////////////// debug接口：
namespace MY_DEBUG 
{



	template <typename T1, typename T2>
	void dispPair(const std::pair<T1, T2>& pair)
	{
		std::cout << "(" << pair.first << ", " << pair.second << ")" << std::endl;
	}


	template<typename T>
	void dispQuat(const Eigen::Quaternion<T>& q)
	{
		std::cout << q.w() << std::endl;
		std::cout << q.vec() << std::endl;
	}
	 
}

  

//////////////////////////////////////////////////////////////////////////////////////////////// 三角网格处理：

// 返回网格中顶点verIdx0的1领域三角片索引；
template <typename DerivedV>
std::unordered_set<int> oneRingTriIdxes(const int verIdx0, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	std::unordered_set<int> nbrIdxes;
	const int versCount = vers.rows();
	const int trisCount = tris.rows();
	if (verIdx0 < 0 || verIdx0 >= versCount)
		return nbrIdxes;

	Eigen::MatrixXi tmpMat = (verIdx0 == tris.array()).select(Eigen::MatrixXi::Ones(trisCount, 3), Eigen::MatrixXi::Zero(trisCount, 3));
	Eigen::VectorXi flagVec = tmpMat.rowwise().sum();
	Eigen::VectorXi tmpVec{ Eigen::VectorXi::LinSpaced(trisCount, 0, trisCount - 1) };
	Eigen::VectorXi selected;
	subFromFlagVec(selected, tmpVec, flagVec);

	for (int i = 0; i < selected.size(); ++i)
		nbrIdxes.insert(selected(i));
	return nbrIdxes;
}


// 返回网格中顶点verIdx0的1领域顶点索引；
template <typename DerivedV>
std::unordered_set<int> oneRingVerIdxes(const int verIdx0, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	std::unordered_set<int> nbrIdxes;
	Eigen::MatrixXi seleTris;
	std::unordered_set<int> nbrIdxSet;
	const int versCount = vers.rows();
	const int trisCount = tris.rows();
	if (verIdx0 < 0 || verIdx0 >= versCount)
		return nbrIdxes;

	std::unordered_set<int> nbrTriIdxes = oneRingTriIdxes(verIdx0, vers, tris);
	if (nbrTriIdxes.empty())
		return nbrIdxes;

	subFromIdxCon(seleTris, tris, nbrTriIdxes);
	int* intPtr = seleTris.data();
	for (int i = 0; i < seleTris.size(); ++i)
	{
		int verIdx = *intPtr++;
		if (verIdx0 != verIdx)
			nbrIdxes.insert(verIdx);
	}
	return nbrIdxes;
}


// 返回网格中顶点verIdx0的1领域三角片索引——按面片法线方向逆时针排列； 
template <typename DerivedV>
std::vector<int> oneRingTriIdxesOrdered(const int verIdx0, \
	const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	std::vector<int> nbrIdxes = oneRingTriIdxes(verIdx0, vers, tris);
	if (nbrIdxes.empty())
		return nbrIdxes;
	const int trisCount = nbrIdxes.size();
	Eigen::MatrixXi nbrTris;
	subFromIdxVec(nbrTris, tris, nbrIdxes);
	std::unordered_map<int, std::pair<int, int>> tmpMap;			// 第一个int是vbIdx, 第二个是vcIdx，第三个是所在三角片的索引；
	for (int i = 0; i < trisCount; ++i)
	{
		int vbIdx, vcIdx;
		if (verIdx0 == tris(i, 0))
		{
			vbIdx = tris(i, 1);
			vcIdx = tris(i, 2);
		}
		else if (verIdx0 = tris(i, 1))
		{
			vbIdx = tris(i, 2);
			vcIdx = tris(i, 0);
		}
		else
		{
			vbIdx = tris(i, 0);
			vcIdx = tris(i, 1);
		}

		auto retPair = tmpMap.insert(std::make_pair(vbIdx, std::make_pair(vcIdx, nbrTris[i])));
		if (!retPair.second)
		{
			// 走到这里说明存在非流形边；
			nbrIdxes.clear();
			return nbrIdxes;
		}
	}
	 
	return nbrIdxes;
}
 

template <typename DerivedI>
bool sortLoopEdges(std::vector<int>& loopVerIdxes, \
	const Eigen::PlainObjectBase<DerivedI>& loopEdges)
{
	loopVerIdxes.clear();
	const int edgesCount = loopEdges.rows();
	loopVerIdxes.reserve(edgesCount);

	// for debug:
	Eigen::MatrixXi tmpMat = loopEdges;
	std::vector<doublet<int>> debugVec = mat2doublets(tmpMat);

	std::unordered_map<int, int> tmpMap;
	for (int i = 0; i < edgesCount; ++i)
	{
		auto retPair = tmpMap.insert(std::make_pair(loopEdges(i, 0), loopEdges(i, 1)));
		if (!retPair.second)
			return false;				// should not happen;
	}

	int head = tmpMap.begin()->first;
	loopVerIdxes.push_back(head);
	for (int i = 0; i < edgesCount - 1; ++i)
	{
		auto iter = tmpMap.find(head);
		if (iter == tmpMap.end())
			return false;										// should not happen;
		int tail = iter->second;
		tmpMap.erase(iter);
		loopVerIdxes.push_back(tail);
		head = tail;
	}

	if (*loopVerIdxes.begin() != tmpMap.begin()->second)
		return false;										// should not happen;

	return true;
}


namespace LAPLACIAN_MASS_MATRIX
{
	// MATLAB——repmat();
	template <typename DerivedA, typename DerivedB>
	void repmat(Eigen::PlainObjectBase<DerivedB>& B, \
		const Eigen::PlainObjectBase<DerivedA>& A, \
		const int repRows, const int repCols)
	{
		assert(repRows > 0);
		assert(repCols > 0);
		B.resize(repRows * A.rows(), repCols * A.cols());

		// copy tiled blocks
		for (int i = 0; i < repRows; i++)
			for (int j = 0; j < repCols; j++)
				B.block(i * A.rows(), j * A.cols(), A.rows(), A.cols()) = A;
	}

	// 生成质量矩阵：
	template <typename Tm, typename Tv>
	bool massMatrix_baryCentric(Eigen::SparseMatrix<Tm>& MM, \
		const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers, \
		const Eigen::MatrixXi& tris)
	{
		const int versCount = vers.rows();
		const int trisCount = tris.rows();
		Eigen::MatrixXd lenMat;
		edge_lengths(lenMat, vers, tris);

		Eigen::VectorXd dblA;
		trisArea(dblA, vers, tris);
		dblA.array() *= 2;
		Eigen::VectorXi MI, MJ;
		Eigen::VectorXd MV;

		// diagonal entries for each face corner
		MI.resize(trisCount * 3);
		MJ.resize(trisCount * 3);
		MV.resize(trisCount * 3);

		MI.segment(0, trisCount) = tris.col(0);
		MI.segment(trisCount, trisCount) = tris.col(1);
		MI.segment(2 * trisCount, trisCount) = tris.col(2);
		MJ = MI;
		repmat(MV, dblA, 3, 1);
		MV.array() /= 6.0;
		sparse(MM, MI, MJ, MV, versCount, versCount);

		return true;
	}

	// laplace光顺——使用质量矩阵；
	template <typename DerivedVo, typename DerivedVi>
	bool laplaceFaring_old(Eigen::PlainObjectBase<DerivedVo>& versOut, \
		const Eigen::PlainObjectBase<DerivedVi>& vers, \
		const Eigen::MatrixXi& tris, const float deltaLB, const unsigned loopCount)
	{
		using ScalarO = typename DerivedVo::Scalar;
		assert(3 == vers.cols() && "assert!!! Input mesh should be in 3-dimension space.");
		const int versCount = vers.rows();
		Eigen::SparseMatrix<ScalarO> I(versCount, versCount);
		I.setIdentity();

		// 1. 计算laplace矩阵：
		Eigen::SparseMatrix<ScalarO> Lmat;
		Eigen::Matrix<ScalarO, Eigen::Dynamic, Eigen::Dynamic> versCopy = vers.array().cast<ScalarO>();
		cotLaplacian(Lmat, versCopy, tris);

		// 2. 光顺的循环：
		for (unsigned i = 0; i < loopCount; ++i)
		{
			// f1. 计算当前质量矩阵：
			Eigen::SparseMatrix<ScalarO> mass;
			massMatrix_baryCentric(mass, versCopy, tris);

			// f2. 解线性方程组 (mass - delta*L) * newVers = mass * newVers
			const Eigen::SparseMatrix<ScalarO> A = (mass - deltaLB * Lmat);
			Eigen::SimplicialLLT<Eigen::SparseMatrix<ScalarO>> solver(A);
			assert(solver.info() == Eigen::Success);
			versOut = solver.solve(mass * versCopy).eval();
			versCopy = versOut;
		}

		return true;
	}
}



//////////////////////////////////////////////////////////////////////////////////////////////// 图像处理相关：

// 矩阵的空域线性滤波——！！！注：滤波mask尺寸必须为奇数！！！
template<typename T>
bool linearSpatialFilter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matOut, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matIn, \
	const Eigen::MatrixXd& mask)
{
	unsigned rows = matIn.rows();
	unsigned cols = matIn.cols();
	matOut.resize(rows, cols);
	matOut.setZero(rows, cols);

	unsigned sizeMask = mask.rows();
	unsigned im = (sizeMask + 1) / 2;					// 掩膜中心的下标；
	unsigned km = im;
	unsigned offset = (sizeMask - 1) / 2;

	// 1. 对输入矩阵进行边缘延拓，生成拓展矩阵：
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matExt = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows + 2 * offset, cols + 2 * offset);
	matExt.block(offset, offset, rows, cols) = matIn;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rowHead = matIn.block(0, 0, offset, cols);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rowTail = matIn.block(rows - offset, 0, offset, cols);

	matExt.block(0, offset, offset, cols) = rowHead;
	matExt.block(rows + offset, offset, offset, cols) = rowTail;

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> colHead = matExt.block(0, offset, rows + 2 * offset, offset);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> colTail = matExt.block(0, cols, rows + 2 * offset, offset);
	matExt.block(0, 0, rows + 2 * offset, offset) = colHead;
	matExt.block(0, cols + offset, rows + 2 * offset, offset) = colTail;

	// 2. 滑动领域操作：
	PARALLEL_FOR(0, rows, [&](int i)
		{
			for (unsigned k = 0; k < cols; ++k)
			{
				unsigned ie = i + offset;					// 此时mask中心元素在matExt中的位置下标；
				unsigned ke = k + offset;
				Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> coveredElems = matExt.block(ie - offset, ke - offset, sizeMask, sizeMask);
				Eigen::MatrixXd tmpMat = coveredElems.array().cast<double>().array() * mask.array();
				matOut(i, k) = static_cast<T>(tmpMat.sum());
			}
		});


	return true;
}


