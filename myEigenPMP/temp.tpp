

///////////////////////////////////////////////////////////////////////////////////////////////////// modeling接口：
// 生成栅格采样点云
template<typename T>
bool genGrids(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, const Eigen::Matrix<T, 1, 3>& origin, \
	const float step, const std::vector<unsigned>& gridCounts)
{
	/*
		bool genGrids(
			Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& gridCenters,			// 输出的栅格点云
			const Eigen::Matrix<T, 1, 3>& origin,														// 栅格原点，即三轴坐标最小的那个栅格点；
			const float step,																						// 采样步长
			const std::vector<unsigned>& gridCounts											// 三元数组，分别是XYZ三轴上的步数；
			)
	*/

	using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

	// 整个距离场的包围盒：
	RowVector3T minp = origin;
	RowVector3T maxp = origin + step * RowVector3T(gridCounts[0] - 1, gridCounts[1] - 1, gridCounts[2] - 1);

	// 生成栅格：
	/*
		按索引增大排列的栅格中心点为：
		gc(000), gc(100), gc(200), gc(300),...... gc(010), gc(110), gc(210), gc(310),...... gc(001), gc(101), gc(201).....

		x坐标：
		x0, x1, x2, x3......x0, x1, x2, x3......x0, x1, x2, x3......
		周期为xCount;
		重复次数为(yCount * zCount)

		y坐标：
		y0, y0, y0...y1, y1, y1...y2, y2, y2.........y0, y0, y0...y1, y1, y1...
		周期为(xCount * yCount);
		重复次数为zCount;
		单个元素重复次数为xCount

		z坐标：
		z0, z0, z0......z1, z1, z1......z2, z2, z2......
		单个元素重复次数为(xCount * yCount)
	*/
	VectorXT xPeriod = VectorXT::LinSpaced(gridCounts[0], minp(0), maxp(0));
	VectorXT yPeriod = VectorXT::LinSpaced(gridCounts[1], minp(1), maxp(1));
	VectorXT zPeriod = VectorXT::LinSpaced(gridCounts[2], minp(2), maxp(2));

	MatrixXT tmpVec0, tmpVec1, tmpVec2, tmpVec11;
	tmpVec0 = kron(VectorXT::Ones(gridCounts[1] * gridCounts[2]), xPeriod);
	tmpVec11 = kron(yPeriod, VectorXT::Ones(gridCounts[0]));
	tmpVec1 = kron(VectorXT::Ones(gridCounts[2]), tmpVec11);
	tmpVec2 = kron(zPeriod, VectorXT::Ones(gridCounts[0] * gridCounts[1]));

	gridCenters.resize(gridCounts[0] * gridCounts[1] * gridCounts[2], 3);
	gridCenters.col(0) = tmpVec0;
	gridCenters.col(1) = tmpVec1;
	gridCenters.col(2) = tmpVec2;

	return true;
}

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
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const float radius, const double deltaTheta = 2 * pi / 30, const bool isCovered = true)
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
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const std::pair<float, float> sizePair, const float SR = 0.5, \
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

	genCylinder(vers, tris, axisVers, circuit);

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



//			circuitToMesh()重载1：triangle库三角剖分——封闭边界线点集得到面网格，可以是平面也可以是曲面，三角片尺寸不可控，不会在网格内部插点。
template <typename T>
bool circuit2mesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& circVers)
{
	// ！！！貌似当前三角剖分有问题！！
	using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;

	unsigned circCount = circVers.rows();
	unsigned versCount = circVers.rows();

	// 0. 边缘环路顶点坐标数据拷入到输出网格中。
	vers = circVers;

	// 输入triangulate()的点集是投影到XOY平面的二维点，原点集应该旋转到合适的角度再投影。

	// 1. 取第一个点、1/3处的点、2/3处的点的所在平面的法向量作为原点集的法向量
	RowVector3T vers1 = circVers.row(0);
	RowVector3T vers2 = circVers.row(versCount / 3);
	RowVector3T vers3 = circVers.row(2 * versCount / 3);
	RowVector3T norm = (vers1 - vers2).cross(vers3 - vers2);
	norm.normalize();

	//  2. 旋转点集使得norm平行于z轴
	Matrix3T rotation = getRotationMat(norm, RowVector3T{ 0, 0, 1 });
	vers = (vers * rotation.transpose()).eval();

	// 3. 旋转后的点数据写入到triangulate()接口的输入结构体中。
	Eigen::MatrixXd vers2D;
	Eigen::MatrixXi edges2D;
	vers2D = vers.transpose().topRows(2).cast<double>();
	edges2D.resize(2, versCount);
	for (unsigned i = 0; i < versCount; i++)
	{
		edges2D(0, i) = i + 1;
		edges2D(1, i) = (i + 1) % versCount + 1;
	}

	triangulateio triIn, triOut;
	triIn.numberofpoints = versCount;
	triIn.pointlist = (double*)vers2D.data();
	triIn.numberofpointattributes = 0;
	triIn.pointattributelist = NULL;
	triIn.pointmarkerlist = NULL;

	triIn.numberofsegments = versCount;
	triIn.segmentlist = (int*)edges2D.data();
	triIn.segmentmarkerlist = NULL;

	triIn.numberoftriangles = 0;
	triIn.numberofcorners = 0;
	triIn.numberoftriangleattributes = 0;

	triIn.numberofholes = 0;
	triIn.holelist = NULL;

	triIn.numberofregions = 0;
	triIn.regionlist = NULL;
	memset(&triOut, 0, sizeof(triangulateio));

	// 4. 执行二维三角剖分，得到输出网格三角片数据，不取顶点坐标数据，顶点坐标数据使用旋转操作前的。
	char triStr[256] = "pY";
	triangulate(triStr, &triIn, &triOut, NULL);
	tris.resize(3, triOut.numberoftriangles);
	MatrixXT norms(versCount, 3);
	std::memcpy(tris.data(), triOut.trianglelist, sizeof(int) * 3 * triOut.numberoftriangles);
	tris.transposeInPlace();
	tris.array() -= 1;

	return true;
}
#endif





