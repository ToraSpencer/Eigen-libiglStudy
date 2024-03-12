 
///////////////////////////////////////////////////////////////////////////////////////////////////// modeling�ӿڣ�
  
#ifdef USE_TRIANGLE_H

// ����Բ������������չʾ��ά�ռ��е�һ��ƽ�棺
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

	// ����XOYƽ����������ԭ���Բ������
	MatrixXO circVers;
	getCircleVers(circVers, radius,  versCount);
	circuit2mesh(versOut, trisOut, circVers);

	// ��һ������任�� 
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


	// ����1��2D���������ʷֵõ������񡪡����Դ���Ҳ���Բ�����
	/*

		ע��
			������Ʊ�����XOYƽ���ڣ�������2D��Ҳ������3D�ģ�
			Ĭ�������ʷ�ģʽΪ"pY"����ʾ�����ʷֳɲ�����ƽ��ֱ��ͼ��


			switch string:
					-p			�����ʷ�����һ��ƽ��ֱ��ͼ
					-r			(refine a previously generated mesh)��һ������������н�һ���������ʷ֣�
					-q			(Quality mesh generation)�����һ����ֵ����"-q30"��ʾ�����ʷֽ���в����Դ���С��30��Ľǣ�
					-a			�����һ����ֵ����"-a5"ָ�������ʷֽ��������Ƭ���������5mm^2;
					-Y			(Prohibits the insertion of Steiner points on the mesh boundary.)
								��ֹ�ڱ�Ե���ϲ����µ㣻
					-YY		prohibits the insertion of Steiner points on any segment, including internal segments.
								��ֹ���κ�ԭ�б��ϲ����µ㣻

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

	// lambda������·��������·�����ݣ�������1��ʼ��
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

	// 1. �������붥�����ݡ���Ե���� 
	Eigen::MatrixXi bdrys;													// ���ɱպϻ�·��һϵ�еı�
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

	// 2. ��������Ķ������ݣ�
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

	// 3. ���������ʷ���Ϣ��
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

		inputTrig.numberofholes = holesCount;					// ������Ŀ
		inputTrig.holelist = holesCount > 0 ? versHole2Dtrans.data() : nullptr;		// ������������׵�ַ����һ����ά���ʾһ������ֻҪ�õ��ڶ��ڼ��ɡ�

		inputTrig.numberofregions = 0;
		inputTrig.regionlist = nullptr;

		memset(&outputTrig, 0, sizeof(TRIANGLE_LIB::triangulateio));		// ��ʼ������ṹ�塣 
	}

	// 4. ִ�������ʷ֣����������   
	char strCopy[512];
	strcpy(strCopy, strSwitcher);
	TRIANGLE_LIB::triangulate(strCopy, &inputTrig, &outputTrig, NULL);					// p����ƽ��ֱ��ͼ��Y��������㡣  

	//		4.1 �����������Ƭ
	int trisOutCount = outputTrig.numberoftriangles;
	{
		trisOut.resize(3, trisOutCount);
		std::memcpy(trisOut.data(), (int*)outputTrig.trianglelist, sizeof(int) * 3 * trisOutCount);
		trisOut.transposeInPlace();
		trisOut.array() -= 1;										// triangle���ж���������1��ʼ����������Ƭ���ݶ�Ҫ��ȥ1
	}

	//		4.2 ����������ƣ�
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


// ����2��2D���������ʷֵõ������񡪡�������
template <typename DerivedVo, typename DerivedI, typename DerivedVi>
bool triangulateVers2Mesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const std::vector<Eigen::VectorXi>& bdryLoops, \
	const char* strSwitcher)
{
	return triangulateVers2Mesh(versOut, trisOut, versIn, bdryLoops, Eigen::MatrixXf{}, strSwitcher);
}


// �����ʷ���������������
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

	// 1. ���ɱ�Ե��Ϣ�Ͷ�����Ϣ
	Eigen::MatrixXi edgesData;													// ���ɱպϻ�·��һϵ�еı�
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

	// 2. ������������Ƭ��Ϣ��
	const int trisCount = trisIn.rows();
	Eigen::MatrixXi trisData;
	if (trisCount > 0)
	{
		trisData = trisIn.transpose();
		trisData.array() += 1;
	}

	// 3. ���������ʷ���Ϣ��
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

		inputTrig.numberofholes = 0;					// ������Ŀ
		inputTrig.holelist = NULL;						// ������������׵�ַ����һ����ά���ʾһ������ֻҪ�õ��ڶ��ڼ��ɡ�

		inputTrig.numberofregions = 0;
		inputTrig.regionlist = NULL;

		memset(&outputTrig, 0, sizeof(TRIANGLE_LIB::triangulateio));		// ��ʼ������ṹ�塣 
	}

	// 3. ִ�������ʷ֣����������   
	char strCopy[512];
	strcpy(strCopy, strSwitcher);
	TRIANGLE_LIB::triangulate(strCopy, &inputTrig, &outputTrig, NULL);					// p����ƽ��ֱ��ͼ��Y��������㡣  

	//		3.1 �����������Ƭ
	int trisOutCount = outputTrig.numberoftriangles;
	{
		trisOut.resize(3, trisOutCount);
		std::memcpy(trisOut.data(), (int*)outputTrig.trianglelist, sizeof(int) * 3 * trisOutCount);
		trisOut.transposeInPlace();
		trisOut.array() -= 1;										// triangle���ж���������1��ʼ����������Ƭ���ݶ�Ҫ��ȥ1
	}

	//		3.2 ����������ƣ�
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


 