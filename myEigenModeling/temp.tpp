 
///////////////////////////////////////////////////////////////////////////////////////////////////// modeling�ӿڣ�
  
#ifdef USE_TRIANGLE_H

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
	Eigen::MatrixXf vers2Dtrans(2, versCount);
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
			vers2Dtrans(0, i) = static_cast<float>(versIn(i, 0));
			vers2Dtrans(1, i) = static_cast<float>(versIn(i, 1));
		}
	}

	// 2. ��������Ķ������ݣ�
	const int holesCount = holeCenters.rows();
	Eigen::MatrixXf  versHole2Dtrans;
	if (holesCount > 0)
	{
		versHole2Dtrans.resize(2, holesCount);
		for (int i = 0; i < holesCount; ++i)
		{
			versHole2Dtrans(0, i) = static_cast<float>(holeCenters(i, 0));
			versHole2Dtrans(1, i) = static_cast<float>(holeCenters(i, 1));
		}
	}

	// 3. ���������ʷ���Ϣ��
	triangulateio inputTrig, outputTrig;
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

		memset(&outputTrig, 0, sizeof(triangulateio));		// ��ʼ������ṹ�塣 
	}

	// 4. ִ�������ʷ֣����������   
	char strCopy[512];
	strcpy(strCopy, strSwitcher);
	triangulate(strCopy, &inputTrig, &outputTrig, NULL);					// p����ƽ��ֱ��ͼ��Y��������㡣  

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
	Eigen::MatrixXf versOutF(2, versOutCount);
	std::memcpy(versOutF.data(), reinterpret_cast<float*>(outputTrig.pointlist), sizeof(float) * 2 * versOutCount);
	versOutF.transposeInPlace();
	versOutF.conservativeResize(versOutCount, 3);
	versOutF.col(2).setZero();
	versOut = versOutF.array().cast<ScalarO>();

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

	assert((2 == versIn.cols() || 3 == versIn.cols()) && "assert!!! input vertices should be in 2 or 3 dimension space.");
	const int versCount = versIn.rows();
	versOut.resize(0, 0);
	trisOut.resize(0, 0);

	// 1. ���ɱ�Ե��Ϣ�Ͷ�����Ϣ
	Eigen::MatrixXi edgesData;													// ���ɱպϻ�·��һϵ�еı�
	Eigen::MatrixXf vers2Dtrans(2, versCount);
	int edgesCount = 0;
	{
		edgesData = edges.array() + 1;
		edgesData.transposeInPlace();
		edgesCount = edgesData.cols();
		for (int i = 0; i < versCount; ++i)
		{
			vers2Dtrans(0, i) = static_cast<float>(versIn(i, 0));
			vers2Dtrans(1, i) = static_cast<float>(versIn(i, 1));
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
	triangulateio inputTrig, outputTrig;
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

		memset(&outputTrig, 0, sizeof(triangulateio));		// ��ʼ������ṹ�塣 
	}

	// 3. ִ�������ʷ֣����������   
	char strCopy[512];
	strcpy(strCopy, strSwitcher);
	triangulate(strCopy, &inputTrig, &outputTrig, NULL);					// p����ƽ��ֱ��ͼ��Y��������㡣  

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
	Eigen::MatrixXf versOutF(2, versOutCount);
	std::memcpy(versOutF.data(), reinterpret_cast<float*>(
		outputTrig.pointlist), sizeof(float) * 2 * versOutCount);
	versOutF.transposeInPlace();
	versOutF.conservativeResize(versOutCount, 3);
	versOutF.col(2).setZero();
	versOut = versOutF.array().cast<ScalarO>();

	return true;
}


// genCylinder()����1�������ɣ��ࣩ���壺
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers, \
	const bool isCovered = true)
{
	/*
	bool genCylinder(
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers,
		Eigen::MatrixXi& tris,
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers,			���ߣ�
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& btmVers,			������·���㣬����Ҫ��XOYƽ���ڣ�
		const bool isCovered																						�Ƿ���
		)


	*/

	// lambda����������������Ƭ��������ѭ�����ã�����һ������һ�㣻
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// ���ظ����ã�tris��������ҪΪ�ա�
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// �����ɱ����Բ�����Ȧ�ĵ�һ�������������
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
	unsigned circVersCount = btmVers.rows();					// �����һȦ�Ķ�������
	unsigned circCount = axisVers.rows();							// Ȧ����
	unsigned versCount = circVersCount * circCount;
	std::vector<MatrixXT> circuitsVec(circCount);									// ÿ��������һȦ���㣻
	std::vector<RowVector3T> sectionNorms(circCount);					// ÿ�������ķ���

	// 1. ��������circCount�������ķ���ÿ�������Ķ��㣻
	for (unsigned i = 0; i < circCount - 1; ++i)
	{
		sectionNorms[i] = axisVers.row(i + 1) - axisVers.row(i);
		sectionNorms[i].normalize();
		Matrix3T rotation = getRotationMat(RowVector3T{ 0, 0, 1 }, sectionNorms[i]);
		circuitsVec[i] = btmVers * rotation.transpose();
		circuitsVec[i].rowwise() += axisVers.row(i);
	}

	// 2. �������һȦ���㣺
	RowVector3T deltaNormAve{ RowVector3T::Zero() };
	for (unsigned i = 0; i < circCount - 2; ++i)
		deltaNormAve += (sectionNorms[i + 1] - sectionNorms[i]);
	deltaNormAve.array() /= (circCount - 2);
	sectionNorms[circCount - 1] = sectionNorms[circCount - 2] + deltaNormAve;
	Matrix3T rotation = getRotationMat(RowVector3T{ 0, 0, 1 }, sectionNorms[circCount - 1]);
	circuitsVec[circCount - 1] = btmVers * rotation.transpose();
	circuitsVec[circCount - 1].rowwise() += axisVers.row(circCount - 1);

	// 3. �������嶥�㣺
	vers.resize(versCount, 3);
	for (unsigned i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.���ɲ�������Ƭ��
	for (unsigned i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);


	// 5. �Ӹǣ�
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


// genCylinder()����2��������Բ��������
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const float radius,\
	const double deltaTheta = 2 * pi / 30, const bool isCovered = true)
{
	// ����XOYƽ���ϲ����ǶȲ���Ϊ��ԲȦ���㣺
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


// genCylinder()����3�������ɷ���������
template <typename T>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, \
	const std::pair<float, float> sizePair, const float SR = 0.5, const bool isCovered = true)
{ 
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// 1. ����XOYƽ���ڵķ��򶥵㣺
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


// genCylinder()����4���������ϵ�����µ���ı߽绷· ���������壺
template <typename T, typename DerivedVa, typename DerivedVt, typename DerivedVb>
bool genCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedVa>& axisVers, const Eigen::PlainObjectBase<DerivedVt>& topLoop, \
	const Eigen::PlainObjectBase<DerivedVb>& btmLoop, const bool isCovered)
{
	assert(btmLoop.rows() == topLoop.rows(), "assert!!! topLoop and btmLoop should have the same amount of vertices.");
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// lambda����������������Ƭ��������ѭ�����ã�����һ������һ�㣻
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// ���ظ����ã�tris��������ҪΪ�ա�
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// �����ɱ����Բ�����Ȧ�ĵ�һ�������������
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

	int circVersCount = topLoop.rows();					// һȦ����滷·�Ķ�������
	int circCount = axisVers.rows();							// Ȧ����
	int versCount = circVersCount * circCount;
	std::vector<MatrixXT> circuitsVec(circCount);							// ÿ��������һȦ���㣻
	std::vector<RowVector3T> sectionNorms(circCount);				// ÿ�������ķ���

	// 1. ����ÿ�������Ļ�·�㼯�����ӵ����ϣ��м�Ļ�·����Ӧ���Ǵӵ׻�·������·�Ĺ��ɣ�
	MatrixXT topLoopCast, btmLoopCast;
	topLoopCast.array() = topLoop.array().cast < T>();
	btmLoopCast.array() = btmLoop.array().cast < T>();
	circuitsVec.begin()->array() = btmLoopCast;
	circuitsVec.rbegin()->array() = topLoopCast;
	
	//			��ֵ�����м�ĺ���滷·��
	if (circCount > 2)
	{
		const int stepCount = circCount - 1;
		MatrixXT stepArrows(circVersCount, 3);
		for (int i = 0; i < circVersCount; ++i)
			stepArrows.row(i) = (topLoopCast.row(i) - btmLoopCast.row(i)) / stepCount;

		for (int i = 1; i < circCount - 1; ++i)
			circuitsVec[i] = btmLoopCast + i * stepArrows;
	}
	 

	// 2. ��������circCount�������ķ���ÿ�������Ķ��㣻
	for (int i = 0; i < circCount - 1; ++i)
	{
		sectionNorms[i] = axisVers.row(i + 1) - axisVers.row(i);
		sectionNorms[i].normalize();
		Matrix3T rotation = getRotationMat(RowVector3T{ 0, 0, 1 }, sectionNorms[i]);

		//		����任��
		circuitsVec[i] = (circuitsVec[i] * rotation.transpose()).eval();
		circuitsVec[i].rowwise() += axisVers.row(i);
	}

	// 3. �������һȦ���㣺
	RowVector3T deltaNormAve{ RowVector3T::Zero() };
	for (int i = 0; i < circCount - 2; ++i)
		deltaNormAve += (sectionNorms[i + 1] - sectionNorms[i]);
	deltaNormAve.array() /= (circCount - 2);
	sectionNorms[circCount - 1] = sectionNorms[circCount - 2] + deltaNormAve;
	Matrix3T rotation = getRotationMat(RowVector3T{ 0, 0, 1 }, sectionNorms[circCount - 1]);

	//		����任��
	circuitsVec[circCount - 1] = (circuitsVec[circCount - 1] * rotation.transpose()).eval();
	circuitsVec[circCount - 1].rowwise() += axisVers.row(circCount - 1);

	// 4. �������嶥�㣺
	vers.resize(versCount, 3);
	for (int i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.���ɲ�������Ƭ��
	for (int i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);

	// 5. �Ӹǣ� 
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


 
// genCylinder()����5���������ϵ�����µ����3D�߽绷· ���������壺
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

	// lambda����������������Ƭ��������ѭ�����ã�����һ������һ�㣻
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// ���ظ����ã�tris��������ҪΪ�ա�
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// �����ɱ����Բ�����Ȧ�ĵ�һ�������������
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

	int circVersCount = topLoop.rows();							// һȦ����滷·�Ķ�������
	int circCount = 2 + layersCount - 1;							// Ȧ����
	int versCount = circVersCount * circCount;
	std::vector<MatrixXO> circuitsVec(circCount);							// ÿ��������һȦ���㣻
	std::vector<RowVector3O> sectionNorms(circCount);				// ÿ�������ķ���

	// 1. ����ÿ�������Ļ�·�㼯�����ӵ����ϣ��м�Ļ�·����Ӧ���Ǵӵ׻�·������·�Ĺ��ɣ�
	MatrixXO topLoopCast, btmLoopCast;
	for (auto& mat : circuitsVec)
		mat.resize(circVersCount, 3);
	topLoopCast.array() = topLoop.array().cast<ScalarO>();
	btmLoopCast.array() = btmLoop.array().cast<ScalarO>();
	circuitsVec.begin()->array() = btmLoopCast;
	circuitsVec.rbegin()->array() = topLoopCast;

	// 2. �����м����Ļ�·���㣺
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

	// 3. �������嶥�㣺
	vers.resize(versCount, 3);
	for (int i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.���ɲ�������Ƭ��
	for (int i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);

	// 5. �Ӹǣ� 
	if (isCovered)
	{
		// 5.1 ���µ��滷·���в��������ʷ֣�������������
		MatrixXO capVersTop, capVersBtm;
		Eigen::MatrixXi capTrisTop, capTrisBtm;
		circuit2mesh(capVersTop, capTrisTop, topLoopCast, topNorm, maxArea);
		circuit2mesh(capVersBtm, capTrisBtm, btmLoopCast, Eigen::Matrix<ScalarN, 1, 3>{-btmNorm}, maxArea);

		// 5.2 �±�������Ƭ����
		for (int i = 0; i < capTrisBtm.rows(); ++i)
		{
			int tmp = capTrisBtm(i, 2);
			capTrisBtm(i, 2) = capTrisBtm(i, 1);
			capTrisBtm(i, 1) = tmp;
		}

		// 5.3 �����������в��������ɵĶ��㣺
		const int newCircVersCount1 = capVersTop.rows() - circVersCount;		// �ϱ��������ʷ�ʱ�¼���Ķ�������
		const int newCircVersCount2 = capVersBtm.rows() - circVersCount;			
		MatrixXO newVersTop = capVersTop.bottomRows(newCircVersCount1);
		MatrixXO newVersBtm = capVersBtm.bottomRows(newCircVersCount2);
		matInsertRows(vers, newVersTop);
		matInsertRows(vers, newVersBtm);

		// 5.4 �����ɵ�����Ƭ���ݸ���ƫ������
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

		// 5.5 �������µ��������Ƭ��
		matInsertRows(tris, capTrisBtm);
		matInsertRows(tris, capTrisTop);
	}

	return true;
}
 

//		���ɷ�������ת�����Σ���ȷ�������XOYƽ��ƽ�л�ֱ��
template <typename T>
bool genAlignedCylinder(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& axisVers, const std::pair<float, float> sizePair, const float SR, \
	const bool isCovered = true)
{
	// ����XOYƽ���ϲ����ǶȲ���Ϊ��ԲȦ���㣺
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// ����XOYƽ���ڵķ��򶥵㣺
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

	// lambda����������������Ƭ��������ѭ�����ã�����һ������һ�㣻
	auto growSurf = [](Eigen::MatrixXi& sideTris, const int circVersCount)->bool
	{
		// ���ظ����ã�tris��������ҪΪ�ա�
		if (circVersCount < 3)
			return false;

		int currentTrisCount = sideTris.rows();
		int currentVersCount = circVersCount + currentTrisCount / 2;
		int startIdx = currentTrisCount / 2;							// �����ɱ����Բ�����Ȧ�ĵ�һ�������������
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

	unsigned circVersCount = circuit.rows();					// �����һȦ�Ķ�������
	unsigned circCount = axisVers.rows();							// Ȧ����
	unsigned versCount = circVersCount * circCount;
	std::vector<MatrixXT> circuitsVec(circCount);									// ÿ��������һȦ���㣻
	std::vector<RowVector3T> sectionNorms(circCount);					// ÿ�������ķ���

	// 1. ��������circCount�������ķ���ÿ�������Ķ��㣻
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

	// 2. �������һȦ���㣺
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

	// 3. �������嶥�㣺
	vers.resize(versCount, 3);
	for (unsigned i = 0; i < circCount; ++i)
		vers.block(0 + circVersCount * i, 0, circVersCount, 3) = circuitsVec[i];

	// 4.���ɲ�������Ƭ��
	for (unsigned i = 1; i <= circCount - 1; ++i)
		growSurf(tris, circVersCount);


	// 5. �Ӹǣ�
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


 