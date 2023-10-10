
/////////////////////////////////////////////////////////////////////////////////////////////////// general tools: 

// ���뺯���ӻ���ָ�����stl����
template<typename T, typename F>
void traverseSTL(T& con, F f)
{
	std::for_each(con.begin(), con.end(), f);
	std::cout << std::endl;
}


// �������
template<typename T, typename F>
void revTraverseSTL(T& con, F f)
{
	std::for_each(con.rbegin(), con.rend(), f);
	std::cout << std::endl;
}



/////////////////////////////////////////////////////////////////////////////////////////////////// debug�ӿڣ�
namespace MY_DEBUG 
{
	// ���������Ƭ����ת��Ϊtriplet��������ʽ��
	template <typename T>
	std::vector<triplet<T>> mat2triplets(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat)
	{
		std::vector<triplet<T>> vec;
		if (0 == mat.rows() || 3 != mat.cols())
			return vec;

		vec.resize(mat.rows());
		for (unsigned i = 0; i < mat.rows(); ++i)
		{
			vec[i].x = mat(i, 0);
			vec[i].y = mat(i, 1);
			vec[i].z = mat(i, 2);
		}

		return vec;
	}


	template <typename T>
	std::vector<doublet<T>> mat2doublets(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat)
	{
		std::vector<doublet<T>> vec;
		if (0 == mat.rows() || 2 != mat.cols())
			return vec;

		vec.resize(mat.rows());
		for (unsigned i = 0; i < mat.rows(); ++i)
		{
			vec[i].x = mat(i, 0);
			vec[i].y = mat(i, 1);
		}

		return vec;
	}


	template <typename T>
	triplet<T> vec2triplet(const Eigen::Matrix<T, 1, 3>& vec)
	{
		triplet<T> t{ vec(0), vec(1), vec(2) };
		return t;
	}


	template <typename T>
	doublet<T> vec2doublet(const Eigen::Matrix<T, 1, 2>& vec)
	{
		doublet<T> d{ vec(0), vec(1) };
		return d;
	}


	// ��������triplet������������Ϊindex��tripletԪ��ʹ��ν��f����ֵΪtrue���򷵻�index; ���Ҳ���������򷵻�-1��
	template <typename T, typename F>
	int findTriplet(const std::vector<triplet<T>>& trips, F f)
	{
		// ν��F����ʽΪbool foo(const triplet<T>& t);
		if (trips.empty())
			return -1;
		for (unsigned i = 0; i < trips.size(); ++i)
			if (f(trips[i]))
				return static_cast<int>(i);

		return -1;
	}


	template <typename T, typename F>
	int findTriplet(const std::vector<doublet<T>>& doubs, F f)
	{
		// ν��F����ʽΪbool foo(const doublet<T>& d);
		if (doubs.empty())
			return -1;
		for (unsigned i = 0; i < doubs.size(); ++i)
			if (f(doubs[i]))
				return static_cast<int>(i);

		return -1;
	}


	// lambda������ӡstd::cout֧�ֵ����ͱ�����
	template <typename T>
	auto disp = [](const T& arg)
	{
		std::cout << arg << ", ";
	};


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





// MATLAB����repmat();
template <typename DerivedA, typename DerivedB>
void repmat(Eigen::PlainObjectBase<DerivedB>& B, const Eigen::PlainObjectBase<DerivedA>& A, const int repRows, const int repCols)
{
	assert(repRows > 0);
	assert(repCols > 0);
	B.resize(repRows * A.rows(), repCols * A.cols());

	// copy tiled blocks
	for (int i = 0; i < repRows; i++)
		for (int j = 0; j < repCols; j++)
			B.block(i * A.rows(), j * A.cols(), A.rows(), A.cols()) = A;
}

 


//////////////////////////////////////////////////////////////////////////////////////////////// ������������

// ���������ж���verIdx0��1��������Ƭ������
template <typename DerivedV>
std::unordered_set<int> oneRingTriIdxes(const int verIdx0, const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
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

// ���������ж���verIdx0��1���򶥵�������
template <typename DerivedV>
std::unordered_set<int> oneRingVerIdxes(const int verIdx0, const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
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


// ���������ж���verIdx0��1��������Ƭ������������Ƭ���߷�����ʱ�����У�
#if 0
template <typename DerivedV>
std::vector<int> oneRingTriIdxesOrdered(const int verIdx0, const Eigen::PlainObjectBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
{
	std::vector<int> nbrIdxes = oneRingTriIdxes(verIdx0, vers, tris);
	if (nbrIdxes.empty())
		return nbrIdxes;
	const int trisCount = nbrIdxes.size();
	Eigen::MatrixXi nbrTris;
	subFromIdxVec(nbrTris, tris, nbrIdxes);
	std::unordered_map<int, std::pair<int, int>> tmpMap;			// ��һ��int��vbIdx, �ڶ�����vcIdx������������������Ƭ��������
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
			// �ߵ�����˵�����ڷ����αߣ�
			nbrIdxes.clear();
			return nbrIdxes;
		}
	}




	return nbrIdxes;
}
#endif


template <typename DerivedI>
bool sortLoopEdges(std::vector<int>& loopVerIdxes, const Eigen::PlainObjectBase<DerivedI>& loopEdges)
{
	loopVerIdxes.clear();
	const int edgesCount = loopEdges.rows();
	loopVerIdxes.reserve(edgesCount);

	// for debug:
	Eigen::MatrixXi tmpMat = loopEdges;
	std::vector<doublet<int>> debugVec = mat2doublets(tmpMat);
	objWriteEdgesMat("E:/currentLoopEdge.obj", tmpMat, g_debugVers);

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



// ��������Ȩ��laplacian (ͬʱҲ�ǸնȾ���(stiffness matrix))�����������ģ�Polygon Laplacian Made Simple [Bunge et al. 2020]
template<typename Tv, typename Tl>
bool cotLaplacian(Eigen::SparseMatrix<Tl>& L, const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();
	L.resize(versCount, versCount);
	L.reserve(10 * versCount);

	Eigen::MatrixXi edges(3, 2);
	edges << 1, 2, 2, 0, 0, 1;

	// Gather cotangents
	Eigen::Matrix<Tl, Eigen::Dynamic, Eigen::Dynamic> C;
	trisCotValues(vers, tris, C);

	std::vector<Eigen::Triplet<Tl> > IJV;
	IJV.reserve(tris.rows() * edges.rows() * 4);

	// Loop over triangles
	for (int i = 0; i < tris.rows(); i++)
	{
		// loop over edges of element
		for (int e = 0; e < edges.rows(); e++)
		{
			int source = tris(i, edges(e, 0));
			int dest = tris(i, edges(e, 1));
			IJV.push_back(Eigen::Triplet<Tl>(source, dest, C(i, e)));
			IJV.push_back(Eigen::Triplet<Tl>(dest, source, C(i, e)));
			IJV.push_back(Eigen::Triplet<Tl>(source, source, -C(i, e)));
			IJV.push_back(Eigen::Triplet<Tl>(dest, dest, -C(i, e)));
		}
	}
	L.setFromTriplets(IJV.begin(), IJV.end());

	return true;
}


// ������������
template <typename Tm, typename Tv>
bool massMatrix_baryCentric(Eigen::SparseMatrix<Tm>& MM, const Eigen::Matrix<Tv, Eigen::Dynamic, Eigen::Dynamic>& vers, \
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


// laplace��˳ 
template <typename T>
bool laplaceFaring(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, \
	const float deltaLB, const unsigned loopCount)
{
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// 1. ����laplace����
	Eigen::SparseMatrix<T> Lmat;
	cotLaplacian(Lmat, vers, tris);

	// 2. ��˳��ѭ����
	MatrixXT versCopy = vers;
	for (unsigned i = 0; i < loopCount; ++i)
	{
		// f1. ���㵱ǰ��������
		Eigen::SparseMatrix<T> mass;
		massMatrix_baryCentric(mass, vers, tris);

		// f2. �����Է����� (mass - delta*L) * newVers = mass * newVers
		const auto& S = (mass - deltaLB * Lmat);
		Eigen::SimplicialLLT<Eigen::SparseMatrix<T>> solver(S);
		assert(solver.info() == Eigen::Success);
		versOut = solver.solve(mass * versCopy).eval();
		versCopy = versOut;
	}

	return true;
}







//////////////////////////////////////////////////////////////////////////////////////////////// ͼ������أ�

// ����Ŀ��������˲�����������ע���˲�mask�ߴ����Ϊ����������
template<typename T>
bool linearSpatialFilter(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matOut, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matIn, \
	const Eigen::MatrixXd& mask)
{
	unsigned rows = matIn.rows();
	unsigned cols = matIn.cols();
	matOut.resize(rows, cols);
	matOut.setZero(rows, cols);

	unsigned sizeMask = mask.rows();
	unsigned im = (sizeMask + 1) / 2;					// ��Ĥ���ĵ��±ꣻ
	unsigned km = im;
	unsigned offset = (sizeMask - 1) / 2;

	// 1. �����������б�Ե���أ�������չ����
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

	// 2. �������������
	PARALLEL_FOR(0, rows, [&](int i)
		{
			for (unsigned k = 0; k < cols; ++k)
			{
				unsigned ie = i + offset;					// ��ʱmask����Ԫ����matExt�е�λ���±ꣻ
				unsigned ke = k + offset;
				Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> coveredElems = matExt.block(ie - offset, ke - offset, sizeMask, sizeMask);
				Eigen::MatrixXd tmpMat = coveredElems.array().cast<double>().array() * mask.array();
				matOut(i, k) = static_cast<T>(tmpMat.sum());
			}
		});


	return true;
}

