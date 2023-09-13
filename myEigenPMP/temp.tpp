 

/////////////////////////////////////////////////////////////////////////////////////////////////// 网格编辑：


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


