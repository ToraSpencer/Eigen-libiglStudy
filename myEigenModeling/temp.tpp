 
///////////////////////////////////////////////////////////////////////////////////////////////////// modeling接口：
  
#ifdef USE_TRIANGLE_H

// genCylinder()重载1——生成（类）柱体：
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers, \
	const bool isCovered = true)
{
	/*
	bool genCylinder(
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers,
		Eigen::MatrixXi& tris,
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers,			轴线；
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers,			横截面回路顶点，必须要在XOY平面内；
		const bool isCovered																						是否封底
		)


	*/

	// lambda——柱体侧面的三角片生长，会循环调用，调用一次生长一层；
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// 会重复调用，tris容器不需要为空。
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// 待生成表面的圆柱体底圈的第一个顶点的索引。
		sideTris.conservativeResize(currentTrisCount + 2 * circVersCount, 3);

		int triIdx = currentTrisCount;
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx + circVersCount - 1, startIdx, startIdx + 2 * circVersCount - 1 };
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx, startIdx + circVersCount, startIdx + 2 * circVersCount - 1 };
		for (int i = startIdx + 1; i < startIdx + circVersCount; ++i)
		{
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i - 1, i, i + circVersCount - 1 };
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i, i + circVersCount, i + circVersCount - 1 };
		}

		return true;
	};

	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	unsigned circVersCount = btmVers.rows();					// 横截面一圈的顶点数；
	unsigned circCount = axisVers.rows();							// 圈数；
	unsigned versCount = circVersCount * circCount;
	std::vector<MatrixXT> circuitsVec(circCount);									// 每个横截面的一圈顶点；
	std::vector<RowVector3T> sectionNorms(circCount);					// 每个横截面的法向；

	// 1. 计算柱体circCount个横截面的法向、每个横截面的顶点；
	for (unsigned i = 0; i < circCount - 1; ++i)
	{
		sectionNorms[i] = axisVers.row(i + 1) - axisVers.row(i);
		sectionNorms[i].normalize();
		Matrix3T rotation = getRotationMat(RowVector3T{ 0, 0, 1 }, sectionNorms[i]);
		circuitsVec[i] = btmVers * rotation.transpose();
		circuitsVec[i].rowwise() += axisVers.row(i);
	}

	// 2. 计算最后一圈顶点：
	RowVector3T deltaNormAve{ RowVector3T::Zero() };
	for (unsigned i = 0; i < circCount - 2; ++i)
		deltaNormAve += (sectionNorms[i + 1] - sectionNorms[i]);
	deltaNormAve.array() /= (circCount - 2);
	sectionNorms[circCount - 1] = sectionNorms[circCount - 2] + deltaNormAve;
	Matrix3T rotation = getRotationMat(RowVector3T{ 0, 0, 1 }, sectionNorms[circCount - 1]);
	circuitsVec[circCount - 1] = btmVers * rotation.transpose();
	circuitsVec[circCount - 1].rowwise() += axisVers.row(circCount - 1);

	// 3. 生成柱体顶点：
	vers.resize(versCount, 3);
	for (unsigned i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.生成侧面三角片：
	for (unsigned i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);


	// 5. 加盖：
	if (isCovered)
	{
		MatrixXT capVers;
		Eigen::MatrixXi capTris1, capTris2;
		circuit2mesh(capVers, capTris1, btmVers);
		capTris2 = capTris1;
		for (int i = 0; i < capTris1.rows(); ++i)
		{
			int tmp = capTris1(i, 2);
			capTris1(i, 2) = capTris1(i, 1);
			capTris1(i, 1) = tmp;
		}
		capTris2.array() += versCount - circVersCount;
		matInsertRows(tris, capTris1);
		matInsertRows(tris, capTris2);
	}

	return true;
}


// genCylinder()重载2——生成圆柱体网格
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const float radius,\
	const double deltaTheta = 2 * pi / 30, const bool isCovered = true)
{
	// 生成XOY平面上采样角度步长为的圆圈顶点：
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	MatrixXT circuit(30, 3);
	circuit.setZero();
	for (unsigned i = 0; i < 30; ++i)
	{
		double theta = deltaTheta * i;
		circuit(i, 0) = radius * cos(theta);
		circuit(i, 1) = radius * sin(theta);
	}

	genCylinder(vers, tris, axisVers, circuit);

	return true;
}


// genCylinder()重载3——生成方柱体网格
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, \
	const std::pair<float, float> sizePair, const float SR = 0.5, const bool isCovered = true)
{ 
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// 1. 生成XOY平面内的方框顶点：
	float length = sizePair.first;
	float width = sizePair.second;
	std::vector<RowVector3T> corners(4);
	corners[0] = RowVector3T{ length / 2, width / 2, 0 };
	corners[1] = RowVector3T{ -length / 2, width / 2, 0 };
	corners[2] = RowVector3T{ -length / 2, -width / 2, 0 };
	corners[3] = RowVector3T{ length / 2, -width / 2, 0 };
	MatrixXT circuit, tmpVers1, tmpVers2, tmpVers3, tmpVers4;
	interpolateToLine(tmpVers1, corners[0], corners[1], SR, true);
	interpolateToLine(tmpVers2, corners[1], corners[2], SR, false);
	interpolateToLine(tmpVers3, corners[2], corners[3], SR, true);
	interpolateToLine(tmpVers4, corners[3], corners[0], SR, false);
	matInsertRows(circuit, tmpVers1);
	matInsertRows(circuit, tmpVers2);
	matInsertRows(circuit, tmpVers3);
	matInsertRows(circuit, tmpVers4);

	// for debug
	objWriteVerticesMat("E:/squareVers.obj", circuit);

	// 2. 
	genCylinder(vers, tris, axisVers, circuit);

	return true;
}


// genCylinder()重载4——输入上底面和下底面的边界环路 ，生成柱体：
template <typename T, typename DerivedVa, typename DerivedVt, typename DerivedVb>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedVa>& axisVers, const Eigen::PlainObjectBase<DerivedVt>& topLoop, \
	const Eigen::PlainObjectBase<DerivedVb>& btmLoop, const bool isCovered)
{
	assert(btmLoop.rows() == topLoop.rows(), "assert!!! topLoop and btmLoop should have the same amount of vertices.");
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// lambda——柱体侧面的三角片生长，会循环调用，调用一次生长一层；
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// 会重复调用，tris容器不需要为空。
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// 待生成表面的圆柱体底圈的第一个顶点的索引。
		sideTris.conservativeResize(currentTrisCount + 2 * circVersCount, 3);

		int triIdx = currentTrisCount;
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx + circVersCount - 1, startIdx, startIdx + 2 * circVersCount - 1 };
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx, startIdx + circVersCount, startIdx + 2 * circVersCount - 1 };
		for (int i = startIdx + 1; i < startIdx + circVersCount; ++i)
		{
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i - 1, i, i + circVersCount - 1 };
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i, i + circVersCount, i + circVersCount - 1 };
		}

		return true;
	};

	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	int circVersCount = topLoop.rows();					// 一圈横截面环路的顶点数；
	int circCount = axisVers.rows();							// 圈数；
	int versCount = circVersCount * circCount;
	std::vector<MatrixXT> circuitsVec(circCount);							// 每个横截面的一圈顶点；
	std::vector<RowVector3T> sectionNorms(circCount);				// 每个横截面的法向；

	// 1. 计算每个横截面的环路点集——从底向上，中间的环路顶点应该是从底环路到顶环路的过渡；
	MatrixXT topLoopCast, btmLoopCast;
	topLoopCast.array() = topLoop.array().cast < T>();
	btmLoopCast.array() = btmLoop.array().cast < T>();
	circuitsVec.begin()->array() = btmLoopCast;
	circuitsVec.rbegin()->array() = topLoopCast;
	
	//			插值生成中间的横截面环路：
	if (circCount > 2)
	{
		const int stepCount = circCount - 1;
		MatrixXT stepArrows(circVersCount, 3);
		for (int i = 0; i < circVersCount; ++i)
			stepArrows.row(i) = (topLoopCast.row(i) - btmLoopCast.row(i)) / stepCount;

		for (int i = 1; i < circCount - 1; ++i)
			circuitsVec[i] = btmLoopCast + i * stepArrows;
	}
	 

	// 2. 计算柱体circCount个横截面的法向、每个横截面的顶点；
	for (int i = 0; i < circCount - 1; ++i)
	{
		sectionNorms[i] = axisVers.row(i + 1) - axisVers.row(i);
		sectionNorms[i].normalize();
		Matrix3T rotation = getRotationMat(RowVector3T{ 0, 0, 1 }, sectionNorms[i]);

		//		仿射变换：
		circuitsVec[i] = (circuitsVec[i] * rotation.transpose()).eval();
		circuitsVec[i].rowwise() += axisVers.row(i);
	}

	// 3. 计算最后一圈顶点：
	RowVector3T deltaNormAve{ RowVector3T::Zero() };
	for (int i = 0; i < circCount - 2; ++i)
		deltaNormAve += (sectionNorms[i + 1] - sectionNorms[i]);
	deltaNormAve.array() /= (circCount - 2);
	sectionNorms[circCount - 1] = sectionNorms[circCount - 2] + deltaNormAve;
	Matrix3T rotation = getRotationMat(RowVector3T{ 0, 0, 1 }, sectionNorms[circCount - 1]);

	//		仿射变换：
	circuitsVec[circCount - 1] = (circuitsVec[circCount - 1] * rotation.transpose()).eval();
	circuitsVec[circCount - 1].rowwise() += axisVers.row(circCount - 1);

	// 4. 生成柱体顶点：
	vers.resize(versCount, 3);
	for (int i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.生成侧面三角片：
	for (int i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);

	// 5. 加盖： 
	if (isCovered)
	{
		MatrixXT capVersTop, capVersBtm;
		Eigen::MatrixXi capTrisTop, capTrisBtm;
		circuit2mesh(capVersTop, capTrisTop, topLoopCast);
		circuit2mesh(capVersBtm, capTrisBtm, btmLoopCast);
		for (int i = 0; i < capTrisBtm.rows(); ++i)
		{
			int tmp = capTrisBtm(i, 2);
			capTrisBtm(i, 2) = capTrisBtm(i, 1);
			capTrisBtm(i, 1) = tmp;
		}
		capTrisTop.array() += versCount - circVersCount;
		matInsertRows(tris, capTrisBtm);
		matInsertRows(tris, capTrisTop);
	} 

	return true;
}


 
// genCylinder()重载5——输入上底面和下底面的3D边界环路 ，生成柱体：
template <typename DerivedVo, typename DerivedVt, typename DerivedVb, typename ScalarN>
bool genCylinder(Eigen::PlainObjectBase<DerivedVo>& vers, 	Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedVt>& topLoop, \
	const Eigen::PlainObjectBase<DerivedVb>& btmLoop, \
	const Eigen::Matrix<ScalarN, 1, 3>& topNorm,\
	const Eigen::Matrix<ScalarN, 1, 3>& btmNorm, \
	const int layersCount, const float maxArea, const bool isCovered)
{
	assert(btmLoop.rows() == topLoop.rows() && "assert!!! topLoop and btmLoop should have the same amount of vertices.");
	assert(layersCount > 0 && "assert!!! layersCount should be a positive integer.");

	using ScalarO = typename DerivedVo::Scalar;
	using MatrixXO = Eigen::Matrix<ScalarO, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3O = Eigen::Matrix<ScalarO, 3, 3>;
	using RowVector3O = Eigen::Matrix<ScalarO, 1, 3>;

	// lambda——柱体侧面的三角片生长，会循环调用，调用一次生长一层；
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// 会重复调用，tris容器不需要为空。
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// 待生成表面的圆柱体底圈的第一个顶点的索引。
		sideTris.conservativeResize(currentTrisCount + 2 * circVersCount, 3);

		int triIdx = currentTrisCount;
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx + circVersCount - 1, startIdx, startIdx + 2 * circVersCount - 1 };
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx, startIdx + circVersCount, startIdx + 2 * circVersCount - 1 };
		for (int i = startIdx + 1; i < startIdx + circVersCount; ++i)
		{
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i - 1, i, i + circVersCount - 1 };
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i, i + circVersCount, i + circVersCount - 1 };
		}

		return true;
	}; 

	int circVersCount = topLoop.rows();							// 一圈横截面环路的顶点数；
	int circCount = 2 + layersCount - 1;							// 圈数；
	int versCount = circVersCount * circCount;
	std::vector<MatrixXO> circuitsVec(circCount);							// 每个横截面的一圈顶点；
	std::vector<RowVector3O> sectionNorms(circCount);				// 每个横截面的法向；

	// 1. 计算每个横截面的环路点集——从底向上，中间的环路顶点应该是从底环路到顶环路的过渡；
	MatrixXO topLoopCast, btmLoopCast;
	for (auto& mat : circuitsVec)
		mat.resize(circVersCount, 3);
	topLoopCast.array() = topLoop.array().cast<ScalarO>();
	btmLoopCast.array() = btmLoop.array().cast<ScalarO>();
	circuitsVec.begin()->array() = btmLoopCast;
	circuitsVec.rbegin()->array() = topLoopCast;

	// 2. 生成中间插入的回路顶点：
	if (circCount > 2)
	{
		MatrixXO arrowSegs(circVersCount, 3);
		for (int i = 0; i < circVersCount; ++i)
		{
			RowVector3O arrow = topLoopCast.row(i) - btmLoopCast.row(i);
			arrowSegs.row(i) = arrow / layersCount;
		}
		for (int i = 1; i < circCount - 1; ++i)
			circuitsVec[i] = btmLoopCast + i * arrowSegs;
	}

	// 3. 生成柱体顶点：
	vers.resize(versCount, 3);
	for (int i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.生成侧面三角片：
	for (int i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);

	// 5. 加盖： 
	if (isCovered)
	{
		// 5.1 上下底面环路进行插点的三角剖分，生成曲面网格；
		MatrixXO capVersTop, capVersBtm;
		Eigen::MatrixXi capTrisTop, capTrisBtm;
		circuit2mesh(capVersTop, capTrisTop, topLoopCast, topNorm, maxArea);
		circuit2mesh(capVersBtm, capTrisBtm, btmLoopCast, Eigen::Matrix<ScalarN, 1, 3>{-btmNorm}, maxArea);

		// 5.2 下表面三角片反向：
		for (int i = 0; i < capTrisBtm.rows(); ++i)
		{
			int tmp = capTrisBtm(i, 2);
			capTrisBtm(i, 2) = capTrisBtm(i, 1);
			capTrisBtm(i, 1) = tmp;
		}

		// 5.3 输出网格点云中插入新生成的顶点：
		const int newCircVersCount1 = capVersTop.rows() - circVersCount;		// 上表面三角剖分时新加入的顶点数；
		const int newCircVersCount2 = capVersBtm.rows() - circVersCount;			
		MatrixXO newVersTop = capVersTop.bottomRows(newCircVersCount1);
		MatrixXO newVersBtm = capVersBtm.bottomRows(newCircVersCount2);
		matInsertRows(vers, newVersTop);
		matInsertRows(vers, newVersBtm);

		// 5.4 新生成的三角片数据附加偏移量：
		int* idxPtr = nullptr;
		int offSet = 0;
		offSet = versCount - circVersCount;
		capTrisTop.array() += offSet;
		idxPtr = capTrisBtm.data();
		offSet = versCount + newCircVersCount1 - circVersCount;
		for (int i = 0; i < capTrisBtm.size(); ++i)
		{
			if ((*idxPtr) >= circVersCount)
				(*idxPtr) += offSet;
			idxPtr++;
		}

		// 5.5 插入上下底面的三角片：
		matInsertRows(tris, capTrisBtm);
		matInsertRows(tris, capTrisTop);
	}

	return true;
}
 

//		生成方柱，旋转分两次，以确保侧面和XOY平面平行或垂直；
template <typename T>
bool genAlignedCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const std::pair<float, float> sizePair, const float SR, \
	const bool isCovered = true)
{
	// 生成XOY平面上采样角度步长为的圆圈顶点：
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// 生成XOY平面内的方框顶点：
	float length = sizePair.first;
	float width = sizePair.second;
	std::vector<RowVector3T> corners(4);
	corners[0] = RowVector3T{ length / 2, width / 2, 0 };
	corners[1] = RowVector3T{ -length / 2, width / 2, 0 };
	corners[2] = RowVector3T{ -length / 2, -width / 2, 0 };
	corners[3] = RowVector3T{ length / 2, -width / 2, 0 };

	MatrixXT circuit, tmpVers1, tmpVers2, tmpVers3, tmpVers4;
	interpolateToLine(tmpVers1, corners[0], corners[1], SR, true);
	interpolateToLine(tmpVers2, corners[1], corners[2], SR, false);
	interpolateToLine(tmpVers3, corners[2], corners[3], SR, true);
	interpolateToLine(tmpVers4, corners[3], corners[0], SR, false);
	matInsertRows(circuit, tmpVers1);
	matInsertRows(circuit, tmpVers2);
	matInsertRows(circuit, tmpVers3);
	matInsertRows(circuit, tmpVers4);

	// lambda——柱体侧面的三角片生长，会循环调用，调用一次生长一层；
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// 会重复调用，tris容器不需要为空。
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// 待生成表面的圆柱体底圈的第一个顶点的索引。
		sideTris.conservativeResize(currentTrisCount + 2 * circVersCount, 3);

		int triIdx = currentTrisCount;
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx + circVersCount - 1, startIdx, startIdx + 2 * circVersCount - 1 };
		sideTris.row(triIdx++) = Eigen::RowVector3i{ startIdx, startIdx + circVersCount, startIdx + 2 * circVersCount - 1 };
		for (int i = startIdx + 1; i < startIdx + circVersCount; ++i)
		{
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i - 1, i, i + circVersCount - 1 };
			sideTris.row(triIdx++) = Eigen::RowVector3i{ i, i + circVersCount, i + circVersCount - 1 };
		}

		return true;
	};

	unsigned circVersCount = circuit.rows();					// 横截面一圈的顶点数；
	unsigned circCount = axisVers.rows();							// 圈数；
	unsigned versCount = circVersCount * circCount;
	std::vector<MatrixXT> circuitsVec(circCount);									// 每个横截面的一圈顶点；
	std::vector<RowVector3T> sectionNorms(circCount);					// 每个横截面的法向；

	// 1. 计算柱体circCount个横截面的法向、每个横截面的顶点；
	for (unsigned i = 0; i < circCount - 1; ++i)
	{
		sectionNorms[i] = axisVers.row(i + 1) - axisVers.row(i);
		sectionNorms[i].normalize();
		Matrix3T rotation1 = getRotationMat(RowVector3T{ 0, 0, 1 }, RowVector3T{ 0, 1, 0 });
		Matrix3T rotation2 = getRotationMat(RowVector3T{ 0, 1, 0 }, sectionNorms[i]);
		Matrix3T rotation = rotation2 * rotation1;
		circuitsVec[i] = circuit * rotation.transpose();
		circuitsVec[i].rowwise() += axisVers.row(i);
	}

	// 2. 计算最后一圈顶点：
	RowVector3T deltaNormAve{ RowVector3T::Zero() };
	for (unsigned i = 0; i < circCount - 2; ++i)
		deltaNormAve += (sectionNorms[i + 1] - sectionNorms[i]);
	deltaNormAve.array() /= (circCount - 2);
	sectionNorms[circCount - 1] = sectionNorms[circCount - 2] + deltaNormAve;
	Matrix3T rotation1 = getRotationMat(RowVector3T{ 0, 0, 1 }, RowVector3T{ 0, 1, 0 });
	Matrix3T rotation2 = getRotationMat(RowVector3T{ 0, 1, 0 }, sectionNorms[circCount - 1]);
	Matrix3T rotation = rotation2 * rotation1;
	circuitsVec[circCount - 1] = circuit * rotation.transpose();
	circuitsVec[circCount - 1].rowwise() += axisVers.row(circCount - 1);

	// 3. 生成柱体顶点：
	vers.resize(versCount, 3);
	for (unsigned i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.生成侧面三角片：
	for (unsigned i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);


	// 5. 加盖：
	if (isCovered)
	{
		MatrixXT capVers;
		Eigen::MatrixXi capTris1, capTris2;
		circuit2mesh(capVers, capTris1, circuit);
		capTris2 = capTris1;
		for (int i = 0; i < capTris1.rows(); ++i)
		{
			int tmp = capTris1(i, 2);
			capTris1(i, 2) = capTris1(i, 1);
			capTris1(i, 1) = tmp;
		}
		capTris2.array() += versCount - circVersCount;
		matInsertRows(tris, capTris1);
		matInsertRows(tris, capTris2);
	}

	return true;
}




#endif


 