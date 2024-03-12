 
///////////////////////////////////////////////////////////////////////////////////////////////////// modeling接口：
  
#ifdef USE_TRIANGLE_H

// 生成圆形面网格，用于展示三维空间中的一个平面：
template <typename DerivedVO, typename DerivedVC, typename DerivedVN>
bool genRoundSurfMesh(Eigen::PlainObjectBase<DerivedVO>& versOut, Eigen::MatrixXi& trisOut, \
	const Eigen::MatrixBase<DerivedVC>& planeCenter, \
	const Eigen::MatrixBase<DerivedVN>& planeNorm, \
	const double radius, const int versCount)
{
	const double eps = 1e-5;
	assert(std::abs(planeNorm.norm() - 1) <= eps && "assert!!! the planeNorm is not a unit vector.");

	using ScalarO = typename DerivedVO::Scalar;
	using Matrix3O = Eigen::Matrix<ScalarO, 3, 3>;
	using Matrix4O = Eigen::Matrix<ScalarO, 4, 4>;
	using MatrixXO = Eigen::Matrix<ScalarO, Eigen::Dynamic, Eigen::Dynamic>;

	versOut.resize(0, 0);
	trisOut.resize(0, 0);

	// 生成XOY平面内中心在原点的圆盘网格
	MatrixXO circVers;
	getCircleVers(circVers, radius,  versCount);
	circuit2mesh(versOut, trisOut, circVers);

	// 求一个仿射变换： 
	Matrix3O rotation;
	Matrix4O affineHomo{Matrix4O::Ones()};
	getRotationMat(rotation, Eigen::RowVector3f{0, 0, 1}, planeNorm);
	affineHomo.topLeftCorner(3, 3) = rotation;
	affineHomo(0, 3) = static_cast<ScalarO>(planeCenter(0));
	affineHomo(1, 3) = static_cast<ScalarO>(planeCenter(1));
	affineHomo(2, 3) = static_cast<ScalarO>(planeCenter(2));

	//
	MatrixXO versOutHomo;
	vers2HomoVers(versOutHomo, versOut);
	versOutHomo = (affineHomo * versOutHomo).eval();
	homoVers2Vers(versOut, versOutHomo);

	return true;
}


	// 重载1：2D点云三角剖分得到面网格——可以带洞也可以不带洞
	/*

		注：
			输入点云必须在XOY平面内，可以是2D的也可以是3D的；
			默认三角剖分模式为"pY"，表示三角剖分成不插点的平面直线图；


			switch string:
					-p			三角剖分生成一个平面直线图
					-r			(refine a previously generated mesh)对一个已有网格进行进一步的三角剖分；
					-q			(Quality mesh generation)后面跟一个数值，如"-q30"表示三角剖分结果中不可以存在小于30°的角；
					-a			后面跟一个数值，如"-a5"指定三角剖分结果中三角片面积不大于5mm^2;
					-Y			(Prohibits the insertion of Steiner points on the mesh boundary.)
								禁止在边缘边上插入新点；
					-YY		prohibits the insertion of Steiner points on any segment, including internal segments.
								禁止在任何原有边上插入新点；

	*/
template <typename DerivedVo, typename DerivedI, typename DerivedVi, typename DerivedVh>
bool triangulateVers2Mesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const std::vector<Eigen::VectorXi>& bdryLoops, \
	const Eigen::PlainObjectBase<DerivedVh>& holeCenters, \
	const char* strSwitcher)
{
	using ScalarO = typename DerivedVo::Scalar;
	using MatrixXR = Eigen::Matrix<TRI_REAL, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3R = Eigen::Matrix<TRI_REAL, 3, 3>;
	using Matrix4R = Eigen::Matrix<TRI_REAL, 4, 4>;
	using RowVector3R = Eigen::Matrix<TRI_REAL, 1, 3>;

	assert((2 == versIn.cols() || 3 == versIn.cols()) && "assert!!! input vertices should be in 2 or 3 dimension space.");

	// lambda——回路索引→回路边数据（索引从1开始）
	auto loop2edges = [](const Eigen::VectorXi& loop)->Eigen::MatrixXi
	{
		const int edgesCount = loop.size();
		Eigen::MatrixXi edges(edgesCount, 2);
		for (int i = 0; i < edgesCount - 1; ++i)
		{
			edges(i, 0) = loop(i) + 1;
			edges(i, 1) = loop(i + 1) + 1;
		}
		edges(edgesCount - 1, 0) = loop(edgesCount - 1) + 1;
		edges(edgesCount - 1, 1) = loop(0) + 1;
		return edges;
	};

	const int versCount = versIn.rows();
	versOut.resize(0, 0);
	trisOut.resize(0, 0);

	// 1. 生成输入顶点数据、边缘数据 
	Eigen::MatrixXi bdrys;													// 连成闭合回路的一系列的边
	MatrixXR vers2Dtrans(2, versCount);
	int bdryEdgesCount = 0;
	{
		for (const auto& vec : bdryLoops)
		{
			Eigen::MatrixXi loopEdges = loop2edges(vec);
			matInsertRows(bdrys, loopEdges);
		}
		bdrys.transposeInPlace();
		bdryEdgesCount = bdrys.cols();

		for (int i = 0; i < versCount; ++i)
		{
			vers2Dtrans(0, i) = static_cast<TRI_REAL>(versIn(i, 0));
			vers2Dtrans(1, i) = static_cast<TRI_REAL>(versIn(i, 1));
		}
	}

	// 2. 生成输入的洞的数据：
	const int holesCount = holeCenters.rows();
	MatrixXR  versHole2Dtrans;
	if (holesCount > 0)
	{
		versHole2Dtrans.resize(2, holesCount);
		for (int i = 0; i < holesCount; ++i)
		{
			versHole2Dtrans(0, i) = static_cast<TRI_REAL>(holeCenters(i, 0));
			versHole2Dtrans(1, i) = static_cast<TRI_REAL>(holeCenters(i, 1));
		}
	}

	// 3. 生成三角剖分信息：
	TRIANGLE_LIB::triangulateio inputTrig, outputTrig;
	{
		inputTrig.numberofpoints = versCount;
		inputTrig.pointlist = vers2Dtrans.data();
		inputTrig.numberofpointattributes = 0;
		inputTrig.pointattributelist = nullptr;
		inputTrig.pointmarkerlist = nullptr;

		inputTrig.numberofsegments = bdryEdgesCount;
		inputTrig.segmentlist = bdrys.data();
		inputTrig.segmentmarkerlist = nullptr;

		inputTrig.numberoftriangles = 0;
		inputTrig.trianglelist = nullptr;

		inputTrig.numberofcorners = 3;
		inputTrig.numberoftriangleattributes = 0;

		inputTrig.numberofholes = holesCount;					// 洞的数目
		inputTrig.holelist = holesCount > 0 ? versHole2Dtrans.data() : nullptr;		// 洞数据数组的首地址，用一个二维点表示一个洞，只要该点在洞内即可。

		inputTrig.numberofregions = 0;
		inputTrig.regionlist = nullptr;

		memset(&outputTrig, 0, sizeof(TRIANGLE_LIB::triangulateio));		// 初始化输出结构体。 
	}

	// 4. 执行三角剖分，拷贝结果：   
	char strCopy[512];
	strcpy(strCopy, strSwitcher);
	TRIANGLE_LIB::triangulate(strCopy, &inputTrig, &outputTrig, NULL);					// p——平面直线图，Y——不插点。  

	//		4.1 生成输出三角片
	int trisOutCount = outputTrig.numberoftriangles;
	{
		trisOut.resize(3, trisOutCount);
		std::memcpy(trisOut.data(), (int*)outputTrig.trianglelist, sizeof(int) * 3 * trisOutCount);
		trisOut.transposeInPlace();
		trisOut.array() -= 1;										// triangle库中顶点索引从1开始，所以三角片数据都要减去1
	}

	//		4.2 生成输出点云：
	const int versOutCount = outputTrig.numberofpoints;
	MatrixXR versOutTmp(2, versOutCount);
	std::memcpy(versOutTmp.data(), reinterpret_cast<TRI_REAL*>(\
		outputTrig.pointlist), sizeof(TRI_REAL) * 2 * versOutCount);
	versOutTmp.transposeInPlace();
	versOutTmp.conservativeResize(versOutCount, 3);
	versOutTmp.col(2).setZero();
	versOut = versOutTmp.array().cast<ScalarO>();

	return true;
}


// 重载2：2D点云三角剖分得到面网格——不带洞
template <typename DerivedVo, typename DerivedI, typename DerivedVi>
bool triangulateVers2Mesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const std::vector<Eigen::VectorXi>& bdryLoops, \
	const char* strSwitcher)
{
	return triangulateVers2Mesh(versOut, trisOut, versIn, bdryLoops, Eigen::MatrixXf{}, strSwitcher);
}


// 三角剖分提升网格质量：
template <typename DerivedVo, typename DerivedIt, typename DerivedVi, typename DerivedIe>
bool triangulateRefineMesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedIt>& trisOut, const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const Eigen::PlainObjectBase<DerivedIt>& trisIn, const Eigen::PlainObjectBase<DerivedIe>& edges, \
	const char* strSwitcher)
{
	using ScalarO = typename DerivedVo::Scalar;
	using MatrixXR = Eigen::Matrix<TRI_REAL, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3R = Eigen::Matrix<TRI_REAL, 3, 3>;
	using Matrix4R = Eigen::Matrix<TRI_REAL, 4, 4>;
	using RowVector3R = Eigen::Matrix<TRI_REAL, 1, 3>;

	assert((2 == versIn.cols() || 3 == versIn.cols()) && "assert!!! input vertices should be in 2 or 3 dimension space.");
	const int versCount = versIn.rows();
	versOut.resize(0, 0);
	trisOut.resize(0, 0);

	// 1. 生成边缘信息和洞的信息
	Eigen::MatrixXi edgesData;													// 连成闭合回路的一系列的边
	MatrixXR vers2Dtrans(2, versCount);
	int edgesCount = 0;
	{
		edgesData = edges.array() + 1;
		edgesData.transposeInPlace();
		edgesCount = edgesData.cols();
		for (int i = 0; i < versCount; ++i)
		{
			vers2Dtrans(0, i) = static_cast<TRI_REAL>(versIn(i, 0));
			vers2Dtrans(1, i) = static_cast<TRI_REAL>(versIn(i, 1));
		}
	}

	// 2. 整理已有三角片信息：
	const int trisCount = trisIn.rows();
	Eigen::MatrixXi trisData;
	if (trisCount > 0)
	{
		trisData = trisIn.transpose();
		trisData.array() += 1;
	}

	// 3. 生成三角剖分信息：
	TRIANGLE_LIB::triangulateio inputTrig, outputTrig;
	{
		inputTrig.numberofpoints = versCount;
		inputTrig.pointlist = vers2Dtrans.data();
		inputTrig.numberofpointattributes = 0;
		inputTrig.pointattributelist = NULL;
		inputTrig.pointmarkerlist = NULL;

		inputTrig.numberofsegments = edgesCount;
		inputTrig.segmentlist = edgesData.data();
		inputTrig.segmentmarkerlist = NULL;

		inputTrig.numberoftriangles = trisCount;
		inputTrig.trianglelist = trisCount > 0 ? trisData.data() : nullptr;

		inputTrig.numberofcorners = 3;
		inputTrig.numberoftriangleattributes = 0;

		inputTrig.numberofholes = 0;					// 洞的数目
		inputTrig.holelist = NULL;						// 洞数据数组的首地址，用一个二维点表示一个洞，只要该点在洞内即可。

		inputTrig.numberofregions = 0;
		inputTrig.regionlist = NULL;

		memset(&outputTrig, 0, sizeof(TRIANGLE_LIB::triangulateio));		// 初始化输出结构体。 
	}

	// 3. 执行三角剖分，拷贝结果：   
	char strCopy[512];
	strcpy(strCopy, strSwitcher);
	TRIANGLE_LIB::triangulate(strCopy, &inputTrig, &outputTrig, NULL);					// p——平面直线图，Y——不插点。  

	//		3.1 生成输出三角片
	int trisOutCount = outputTrig.numberoftriangles;
	{
		trisOut.resize(3, trisOutCount);
		std::memcpy(trisOut.data(), (int*)outputTrig.trianglelist, sizeof(int) * 3 * trisOutCount);
		trisOut.transposeInPlace();
		trisOut.array() -= 1;										// triangle库中顶点索引从1开始，所以三角片数据都要减去1
	}

	//		3.2 生成输出点云：
	const int versOutCount = outputTrig.numberofpoints;
	MatrixXR versOutTmp(2, versOutCount);
	std::memcpy(versOutTmp.data(), reinterpret_cast<TRI_REAL*>(
		outputTrig.pointlist), sizeof(TRI_REAL) * 2 * versOutCount);
	versOutTmp.transposeInPlace();
	versOutTmp.conservativeResize(versOutCount, 3);
	versOutTmp.col(2).setZero();
	versOut = versOutTmp.array().cast<ScalarO>();

	return true;
}
  

#endif


 