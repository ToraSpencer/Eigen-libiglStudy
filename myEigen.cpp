#include "myEigen.h"
 

/// /////////////////////////////////////////////////////////////////////////////////////////// DEBUG �ӿ�
static std::string g_debugPath = "E:/";
 
static void debugDisp()			// �ݹ���ֹ
{						//		�ݹ���ֹ��Ϊ�޲λ���һ�����������ζ����ԡ�
	std::cout << std::endl;
	return;
}

template <typename T, typename... Types>
static void debugDisp(const T& firstArg, const Types&... args)
{
	std::cout << firstArg << " ";
	debugDisp(args...);
}


template <typename T, int M, int N>
static void dispData(const Eigen::Matrix<T, M, N>& m)
{
	auto dataPtr = m.data();
	unsigned elemsCount = m.size();

	for (unsigned i = 0; i < elemsCount; ++i)
		std::cout << dataPtr[i] << ", ";

	std::cout << std::endl;
}


template <typename Derived>
static void dispData(const Eigen::PlainObjectBase<Derived>& m)
{
	int m0 = m.RowsAtCompileTime;
	int n0 = m.ColsAtCompileTime;

	auto dataPtr = m.data();
	unsigned elemsCount = m.size();

	for (unsigned i = 0; i < elemsCount; ++i)
		std::cout << dataPtr[i] << ", ";

	std::cout << std::endl;
}


template <typename Derived>
static void dispElem(const Eigen::MatrixBase<Derived>& m)
{
	const Derived& mm = m.derived();
	std::cout << mm(1, 1) << std::endl;
}


template<typename DerivedV>
static void debugWriteVers(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteVerticesMat(path, vers);
}


template<typename T>
static void debugWriteMesh(const char* name, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteMeshMat(path, vers, tris);
}


static void debugWriteMesh(const char* name, T_MESH::Basic_TMesh& mesh)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	mesh.save(path);
}


template<typename DerivedV>
static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteEdgesMat(path, edges, vers);
}



///////////////////////////////////////////////////////////////////////////////////////////////// ʵ��

unsigned readNextData(char*& pszBuf, unsigned& nCount, char* validData, const unsigned nMaxSize) 
{
	unsigned nIndx = 0;

	while ((pszBuf[0] == ' ') ||
		(pszBuf[0] == '\n') ||
		(pszBuf[0] == '\t') ||
		(pszBuf[0] == '\r'))
	{
		pszBuf++;
		nCount++;
	}

	while ((pszBuf[0] != ' ') &&
		(pszBuf[0] != '\n') &&
		(pszBuf[0] != '\t') &&
		(pszBuf[0] != '\r') &&
		(pszBuf[0] != '\null') &&
		(pszBuf[0] != 0) &&
		(nIndx < nMaxSize))
	{
		validData[nIndx++] = pszBuf[0];
		pszBuf++;
		nCount++;
	}

	validData[nIndx] = 0;
	return nIndx;
};


// ��ȡ����OBJ�ļ��е����ݣ��洢���������ϵ��ʾ�ĵ��ƾ����У�ע�������ھ��������б�ʾ�ģ�����ά��Ԫ��ʼ��Ϊ1��
void objReadVerticesHomoMat(Eigen::MatrixXf& vers, const char* fileName)
{
	char* pTmp = NULL;
	std::ifstream ifs(fileName);			// cube bunny Eight
	if (false == ifs.is_open())
	{
		return;
	}

	std::streampos   pos = ifs.tellg();			//  save   current   position   
	ifs.seekg(0, std::ios::end);
	unsigned fileLen = (unsigned)ifs.tellg();
	if (0 == fileLen)
	{
		return;
	}

	ifs.seekg(pos);				  //   restore   saved   position   
	char* pFileBuf = new char[fileLen + 1];
	std::memset(pFileBuf, 0, fileLen + 1);
	ifs.read(pFileBuf, fileLen);
	char tmpBuffer[1024];
	unsigned nMaxSize = 1024;
	pTmp = pFileBuf;
	unsigned nReadLen = 0;
	unsigned nRet = 0;

	vers.resize(0, 0);
	int rows = vers.cols();
	while (nReadLen < fileLen)
	{
		nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
		if (0 == nRet)
			break;

		// ������Ϣ		
		if (std::strcmp(tmpBuffer, "v") == 0)
		{
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;

			Eigen::Vector4f ver(Eigen::Vector4f::Ones());
			ver(0) = (float)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(1) = (float)atof(tmpBuffer);
			nRet = readNextData(pTmp, nReadLen, tmpBuffer, nMaxSize);
			if (0 == nRet)
				break;
			ver(2) = (float)atof(tmpBuffer);

			vers.conservativeResize(++rows, 4);
			vers.row(rows - 1) = ver;
		}
		else
			break;
	}
	vers.transposeInPlace();
	delete[] pFileBuf;
};


// �������ϵ��������д�뵽OBJ�ļ��У�ע�������ھ��������б�ʾ�ģ�����ά��Ԫ��ʼ��Ϊ1��
void objWriteVerticesHomoMat(const char* fileName, const Eigen::MatrixXf& vers)
{
	std::ofstream dstFile(fileName);
	for (int i = 0; i < vers.cols(); i++)
	{
		dstFile << "v " << vers(0, i) << " " << vers(1, i) << " " << vers(2, i) << std::endl;
	}
	dstFile.close();

}


//	��ͨ���ƾ���ת��Ϊ�������ϵ�µĵ��ƾ���
void vers2homoVers(Eigen::MatrixXf& homoVers, const Eigen::MatrixXf& vers)
{
	homoVers = Eigen::MatrixXf::Ones(4, vers.rows());
	homoVers.topRows(3) = vers.transpose();
}


Eigen::MatrixXf vers2homoVers(const Eigen::MatrixXf& vers)
{
	Eigen::MatrixXf homoVers = Eigen::MatrixXf::Ones(4, vers.rows());
	homoVers.topRows(3) = vers.transpose();
	return homoVers;
}


// �������ϵ�µĵ��ƾ���任Ϊ��ͨ���ƾ���
void homoVers2vers(Eigen::MatrixXf& vers, const Eigen::MatrixXf& homoVers)
{
	Eigen::MatrixXf tempMat = homoVers.transpose();
	vers = tempMat.leftCols(3);
}


Eigen::MatrixXf homoVers2vers(const Eigen::MatrixXf& homoVers)
{
	Eigen::MatrixXf tempMat = homoVers.transpose();
	Eigen::MatrixXf vers = tempMat.leftCols(3);
	return vers;
}


void printDirEigen(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& dir)
{
	Eigen::MatrixXf line;

	const float SR = 0.5;			// �ռ������SR������������������ľ��루��λmm��
	const float length = 10;
	int versCount = std::round(length / SR);
	line.resize(versCount + 1, 3);
	line.row(0) = origin;
	for (int i = 1; i <= versCount; i++)
		line.row(i) = line.row(0) + SR * dir * i;

	objWriteVerticesMat(pathName, line);
};


void printCoordinateEigen(const char* pathName, const Eigen::RowVector3f& origin, const Eigen::RowVector3f& xdir, \
	const Eigen::RowVector3f& ydir, const Eigen::RowVector3f& zdir)
{
	const float SR = 0.5;			// �ռ������SR������������������ľ��루��λmm��
	const float length = 10;
	int versCount = std::round(length / SR);
	Eigen::MatrixXf line1(versCount, 3), line2(versCount, 3), line3(versCount, 3);
	for (int i = 0; i < versCount; i++)
		line1.row(i) = origin + SR * xdir * (i + 1);
	for (int i = 0; i < versCount; i++)
		line2.row(i) = origin + SR * ydir * (i + 1);
	for (int i = 0; i < versCount; i++)
		line3.row(i) = origin + SR * zdir * (i + 1);

	Eigen::MatrixXf line = origin;
	matInsertRows<float>(line, line1);
	matInsertRows<float>(line, line2);
	matInsertRows<float>(line, line3);
	objWriteVerticesMat(pathName, line);
}


// ��������ת��Ϊ����������
Eigen::VectorXi flagVec2IdxVec(const Eigen::VectorXi& flag)
{
	Eigen::VectorXi idxVec;

	for (unsigned i = 0; i < flag.rows(); ++i)
		if (0 != flag(i) && 1 != flag(i))
			return idxVec;

	std::vector<int> tempVec;
	tempVec.reserve(flag.rows());
	for (unsigned i = 0; i< flag.rows(); ++i) 
		if (flag(i) > 0)
			tempVec.push_back(i);
	idxVec = vec2Vec<int>(tempVec);
	
	return idxVec;
}


// ��������ת��Ϊ����������
Eigen::VectorXi IdxVec2FlagVec(const Eigen::VectorXi& idxVec, const unsigned size)
{
	Eigen::VectorXi flag;

	if (idxVec.rows() > size)
		return flag;

	for (unsigned i = 0; i < idxVec.rows(); ++i) 
		if (idxVec(i) >= size)
			return flag;

	flag = Eigen::VectorXi::Zero(size);
	for (unsigned i = 0; i < idxVec.rows(); ++i)
		flag(idxVec(i)) = 1;

	return flag;
}


// ����ʽ��ֵ
void polyInterpolation() 
{

}


// ��˹��ֵ
void gaussInterpolation() {}


// ��С���˶���ʽ������ߣ�
void leastSquarePolyFitting()
{}


std::pair<int, int> decodeEdge(const std::int64_t code)
{
	int a = static_cast<int>(code >> 32);
	int b = static_cast<int>(code - (static_cast<std::int64_t>(a) << 32));
	return std::make_pair(a, b);
}


std::vector<int> decodeTrianagle(const std::uint64_t code)
{
	std::uint64_t a = code >> 42;
	std::uint64_t resi = code - (a << 42);
	std::uint64_t b = resi >> 21;
	std::uint64_t c = resi - (b << 21);
	return std::vector<int>{static_cast<int>(a), static_cast<int>(b), static_cast<int>(c)};
}


// ����T_MESH::di_cell��Ӧ�������Χ������
void genAABBmesh(const T_MESH::di_cell& cell, Eigen::MatrixXd& vers, Eigen::MatrixXi& tris)
{
	vers.resize(0, 0);
	tris.resize(0, 0);
	T_MESH::Point arrow = cell.Mp - cell.mp;
	T_MESH::Point newOri = (cell.mp + cell.Mp) / 2.0;

	genCubeMesh(vers, tris);
	vers.col(0) *= arrow.x;
	vers.col(1) *= arrow.y;
	vers.col(2) *= arrow.z;
	vers.rowwise() += Eigen::RowVector3d{ newOri.x, newOri.y, newOri.z };
}




// ����myEigen�еĽӿ�
namespace TEST_MYEIGEN 
{
	// ���Ա�����룺
	void test0() 
	{
		std::vector<int> retVec = decodeTrianagle(encodeTriangle(65535, 65534, 65533));
		traverseSTL(retVec, disp<int>);
		std::cout << std::hex << encodeTriangle(65535, 65534, 65533) << std::endl;
		std::cout << std::dec << encodeTriangle(1, 1, 2) << std::endl;;

		retVec.clear();
		retVec = decodeTrianagle(encodeTriangle(651135, 522534, 653533));
		traverseSTL(retVec, disp<int>);

		std::cout << "finished." << std::endl;
	}


	// ������������triangleGrow();
	void test1() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/����/fatTeeth.obj");
		objWriteMeshMat("E:/triangleGrowInput.obj", vers, tris);

		std::vector<Eigen::MatrixXd> meshesVersOut;
		std::vector<Eigen::MatrixXi> meshesTrisOut;
		triangleGrowSplitMesh(meshesVersOut, meshesTrisOut, vers, tris);

		for (int i = 0; i<meshesVersOut.size(); ++i) 
		{
			char str[256];
			sprintf_s(str, "E:/triangleGrowSplitOut%d.obj", i);
			objWriteMeshMat(str, meshesVersOut[i], meshesTrisOut[i]);
		}

		triangleGrow(versOut, trisOut, vers, tris, 0);
		objWriteMeshMat("E:/triangleGrowOut0.obj", versOut, trisOut);

		versOut.resize(0, 0);
		trisOut.resize(0, 0);
		triangleGrow(versOut, trisOut, vers, tris, 16666);
		objWriteMeshMat("E:/triangleGrowOut16666.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// ������������ʵ�ֵ��ص�����Ƭ��⣺
	void test11()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/����/jawMeshDense_algSimp_60000.obj");
		objWriteMeshMat("E:/triangleGrowInput.obj", vers, tris);

		std::vector<std::pair<int, int>> opTrisPairs;
		int olCount = findOverLapTris(opTrisPairs, vers, tris);


		debugDisp("finished.");
	}


	// ������������ʵ�ֵ���������Ƭ���������
	void test111() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/����/jawMeshDense_algSimp_60000.obj");
		objWriteMeshMat("E:/triangleGrowInput.obj", vers, tris);

		int corrCount = correctTriDirs(trisOut, vers, tris, 0);
		debugDisp("correctTriDirs == ", corrCount);
		debugWriteMesh("meshsOut", vers, trisOut);
		
		debugDisp("finished.");
	}


	//		���԰��������αߵ�����Ƭ����������
	void test1111() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/����/meshRepairInput_d_meshArranged.obj");
		objWriteMeshMat("E:/triangleGrowInput.obj", vers, tris);
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();

		Eigen::MatrixXi nmnEdges;
		std::vector<std::pair<int, std::pair<int, int>>> nmnInfos;
		if (nonManifoldEdges(nmnEdges, nmnInfos, tris) > 0)
			debugWriteEdges("nmnEdgesOrigin", nmnEdges, vers);
 
		if (!triangleGrowOuterSurf(versOut, trisOut, vers, tris, 0))
			debugDisp("triangleGrowOuterSurf() failed!");
		debugWriteMesh("triangleGrowOuterSurf", versOut, trisOut);
		debugDisp("ȥ������Ƭ������", trisCount - trisOut.rows());


		// �������������Ƿ��з����αߣ�
		nmnEdges.resize(0, 0);
		nmnInfos.clear();
		if (nonManifoldEdges(nmnEdges, nmnInfos, trisOut) > 0)
			debugWriteEdges("nmnEdges", nmnEdges, versOut);


		std::cout << "finished." << std::endl;
	}


	// ����ظ�����Ƭ���ظ����㣬
	void test2() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/����/shrinkedMeshDirty.obj");

		Eigen::MatrixXi bdrys;
		bdryEdges(bdrys, tris);
		std::cout << "boundary edges count is " << bdrys.rows() << std::endl;
		objWriteEdgesMat("E:/bdrys.obj", bdrys, vers);

		std::vector<int> repIdxes;
		findRepTris(repIdxes, tris);


		std::cout << "finished." << std::endl;
	}


	// ���Ի�ȡ��������������ԵĽӿڣ�
	void test3() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		bool retFlag = true;
		objReadMeshMat(vers, tris, "E:/����/jawCore.obj");

		const unsigned versCount = vers.rows();
		const unsigned trisCount = tris.rows();

		// ��������ÿ������Ƭ�������
		Eigen::VectorXd trisAreaVec;
		if (!trisArea(trisAreaVec, vers, tris))
			return;

		// �ֶ������һ������Ƭ�������
		Eigen::RowVector3d arrow1, arrow2, dir1, dir2, va, vb, vc;
		va = vers.row(tris(0,0));
		vb = vers.row(tris(0, 1));
		vc = vers.row(tris(0, 2));
		arrow1 = vb - va;
		arrow2 = vc - va;
		dir1 = arrow1.normalized();
		dir2 = arrow2.normalized();
		double cosTheta = dir1.dot(dir2);
		double sinTheta = std::sqrt(1 - cosTheta*cosTheta);
		double area0 = 0.5 * arrow1.norm() * arrow2.norm() * sinTheta;
		std::cout << "area0 == " << area0 << std::endl;
		dispVecSeg(trisAreaVec, 0, 9);

		std::cout << "finished." << std::endl;
	}


	// ���Ի�����ѧ�ӿ�
	void test4()
	{
		dispPair(cart2polar(5.0, -5.0));
		dispPair(polar2cart(20.0, 3.14159/3));

		std::cout << "finished." << std::endl;
	}


	// �������ɻ���ͼ�εĽӿڣ�
	void test5() 
	{
		Eigen::MatrixXf axis(15, 3);
		axis.setZero();
		double deltaTheta = pi / 10;
		for (unsigned i = 0; i < 15; ++i)
		{
			double theta = deltaTheta * i;
			axis(i, 0) = 50 * cos(theta);
			axis(i, 1) = 50 * sin(theta);
		}
 
		debugWriteVers("axis", axis);
		std::cout << "finished." << std::endl;

		// �����ʷ֣�
		Eigen::MatrixXd surfVers, circleVers;
		Eigen::MatrixXi surfTris;
		objReadVerticesMat(circleVers, "E:/����/circleVers.obj");
		circuit2mesh(surfVers, surfTris, circleVers);
		debugWriteMesh("surfMesh", surfVers, surfTris);

		// �������壺
		Eigen::MatrixXf cylinderVers;
		Eigen::MatrixXi cylinderTris;
		genCylinder(cylinderVers, cylinderTris, axis, 10);				// ����Բ��
		debugWriteMesh("cylinder", cylinderVers, cylinderTris);

		cylinderVers.resize(0, 0);
		cylinderTris.resize(0, 0);
		axis.resize(0, 0);
		matInsertRows<float , 3>(axis, Eigen::RowVector3f{0, 0, 0});
		matInsertRows<float , 3>(axis, Eigen::RowVector3f{ 0, 0, 1 });

		genCylinder(cylinderVers, cylinderTris, axis, std::make_pair(10.0, 20.0));
		debugWriteMesh("pillar", cylinderVers, cylinderTris);

		cylinderVers.resize(0, 0);
		cylinderTris.resize(0, 0);
		genAlignedCylinder(cylinderVers, cylinderTris, axis, std::make_pair(1.5, 1.5), 0.5); 
		debugWriteMesh("AlignedPillar", cylinderVers, cylinderTris);

		std::cout << "finished." << std::endl;
	}

	
	// ���Կռ�任��صĽӿڣ�
	void test6()
	{
		Eigen::Matrix3d rotation;
		Eigen::RowVector3d axisArrow;
		Eigen::RowVector3d arrow1, arrow2, arrowRotted;
		float theta = pi;

		// ����������ת��
		theta = pi / 3;
		axisArrow = Eigen::RowVector3d{20, 0, 0};
		rotation = getRotationMat(axisArrow, theta);
		arrowRotted = arrow1 * rotation.transpose();
		dispVec<double, 3>(arrowRotted);

		// ��ת��Ŀ��������
		arrow1.setRandom();
		arrow2.setRandom();
		arrow1.normalize();
		arrow2.normalize();
		dispVec<double, 3>(arrow1);
		dispVec<double, 3>(arrow2);
		rotation = getRotationMat(arrow1, arrow2);
		arrowRotted = arrow1 * rotation.transpose();
		dispVec<double, 3>(arrowRotted);

		std::cout << "finished." << std::endl;
	}


	// ���Ҷ��ͱ�Ե���ߣ�
	void test7() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		std::vector<std::vector<int>> holes, bdrySegs;

		objReadMeshMat(vers, tris, "E:/����/meshRepairInput_c_noOpTris.obj");
		debugWriteMesh("meshInput", vers, tris);

		findHolesBdrySegs(holes, bdrySegs, vers, tris);


		debugDisp("finished.");
	}
}



// ����ͼ����
namespace TEST_DIP
{
	void test0()
	{
		Eigen::MatrixXi mat1(9, 9);
		mat1.setRandom();
		
		Eigen::MatrixXd mask(3, 3);
		mask.setOnes();

		Eigen::MatrixXi matOut;
		linearSpatialFilter(matOut, mat1, mask);

		dispMat(mat1);
		dispMat(matOut);

		std::cout << "finished." << std::endl;
	}

}

 

// ��������������tmesh
namespace TEST_TMESH
{
	// TMESH��IO����������
	void test0()
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();												// ������This is mandatory
		T_MESH::Basic_TMesh tmesh1;
		T_MESH::Node* nodePtr = nullptr;
		tmesh1.load("E:/����/tooth.obj");
		tmesh1.save("E:/meshInput.obj");

		const int versCount = tmesh1.V.numels();
		const int edgesCount = tmesh1.E.numels();
		const int trisCount = tmesh1.T.numels();

		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		TMesh2MeshMat(vers, tris, tmesh1);
		debugWriteMesh("meshInputCopy", vers, tris);

		// �������ݣ�
		nodePtr = tmesh1.V.head();
		nodePtr = nodePtr->next();					// �ڶ�������
		T_MESH::Vertex* vPtr = reinterpret_cast<T_MESH::Vertex*>(nodePtr->data);
		T_MESH::Vertex& ver0 = *vPtr;

		// �����ݣ�
		nodePtr =tmesh1.E.head();
		T_MESH::Edge* ePtr = reinterpret_cast<T_MESH::Edge*>(nodePtr->data);
		T_MESH::Edge& e0 = *ePtr;

		// ����Ƭ���ݣ�
		nodePtr = tmesh1.T.head();
		T_MESH::Triangle* tPtr = reinterpret_cast<T_MESH::Triangle*>(nodePtr->data);
		T_MESH::Triangle& t0 = *tPtr;

		// ����ת����
		T_MESH::Basic_TMesh tMesh2;
		meshMat2tMesh(tMesh2, vers, tris);
		debugWriteMesh("tMesh2", tmesh1);

		std::cout << "finished." << std::endl;
	}


	// ����T_MESH::di_cell��
	void test1() 
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh mesh;
		T_MESH::Node* nodePtr = nullptr;
		mesh.load("E:/����/tooth.obj");				// ��������������
		debugWriteMesh("meshInput", mesh);

		T_MESH::di_cell cell(&mesh);			// ������������di_cell����
		Eigen::MatrixXd aabbVers;
		Eigen::MatrixXi aabbTris;

		// ��ӡcell��Ӧ�İ�Χ�У�
		genAABBmesh(cell, aabbVers, aabbTris);
		debugWriteMesh("aabb", aabbVers, aabbTris);

		// fork()������
		T_MESH::di_cell cell2 = *cell.fork();
		genAABBmesh(cell, aabbVers, aabbTris);
		debugWriteMesh("aabb11", aabbVers, aabbTris);
		genAABBmesh(cell2, aabbVers, aabbTris);
		debugWriteMesh("aabb12", aabbVers, aabbTris);

		std::cout << "finished." << std::endl;
	}


	// ����TMESH�еļ���Ԫ�ر�ǹ��ܣ�
	void test2() 
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/����/tooth.obj");
		debugWriteMesh("meshInput", tMesh);

		// 1.�������㣬���x����С��0�Ķ��㣻
		traverseVersList(tMesh.V, [&](T_MESH::Vertex* verPtr) 
			{
				if (verPtr->x < 0)
					MARK_VISIT(verPtr);
			});

		// 2. ѡ�����а���ѡ�ж��������Ƭ��
		tMesh.growSelection();

		// 3. ɾ������ѡ������Ƭ����������ĵ�ͱߣ�
		tMesh.removeSelectedTriangles();
		debugWriteMesh("meshOut", tMesh);

		std::cout << "finished." << std::endl;
	}


	// �������ʡ������޸�
	void test3() 
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/����/tooth.obj");
		debugWriteMesh("meshInput", tMesh);

		// ��������������ʣ�
		const char* retStr = tMesh.checkConnectivity();				// ��������������᷵���ַ�����
		if (NULL != retStr)
			std::cout << retStr << std::endl;

		std::cout << "finished." << std::endl;
	}


	// ȥ���˻�����Ƭ��
	void test4()
	{
		int max_iters = 5;
		int inner_loops = 6;							// ÿ�δ�ѭ����ȥ���˻�����Ƭ��ȥ���Խ��ĵ���������
		int holesCount = 0;
		bool flagIsct = false;
		bool flagDeg = false;

		T_MESH::TMesh::init();
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/����/meshRepairInput.obj");
		debugWriteMesh("meshInput", tMesh);

		int removedCount = tMesh.removeSmallestComponents();						// d_boundaries, d_handles, d_shells��ֵ
		if (removedCount > 0)
			std::cout << "����������������" << removedCount + 1 << "������ͨ����" << std::endl;

		tMesh.deselectTriangles();
		tMesh.invertSelection();

		flagDeg = tMesh.strongDegeneracyRemoval(inner_loops);			// ȫ������ɹ�����true�� ���򷵻�false

		debugWriteMesh("meshOut", tMesh);
		if(!flagDeg)
			debugDisp("�������˻�����Ƭû��ȫ������ɹ���");
		
		debugDisp("finished.");
	}


	// ѭ������ȥ���˻�����Ƭ������
	void test44() 
	{
		int max_iters = 5;			
		int inner_loops = 6;							// ÿ�δ�ѭ����ȥ���˻�����Ƭ��ȥ���Խ��ĵ���������
		int holesCount = 0;
		bool flagIsct = false;
		bool flagDeg = false;

		T_MESH::TMesh::init();				 
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/����/jawMeshDense_qslim_150000_noDegOpTris.obj");
		debugWriteMesh("meshInput", tMesh);	

		int removedCount = tMesh.removeSmallestComponents();						// d_boundaries, d_handles, d_shells��ֵ
		if (removedCount > 1)
			std::cout << "����������������" << removedCount << "������ͨ����" << std::endl;
		
		tMesh.deselectTriangles();
		tMesh.invertSelection();

		int loopCount = 1;
		for(int i = 0; i < max_iters; ++i)
		{
			debugDisp("loopCount == ", loopCount++);
			holesCount = tMesh.boundaries();
			if (holesCount)
				tMesh.fillSmallBoundaries(0, true);
			flagDeg = tMesh.strongDegeneracyRemoval(inner_loops);			// ȫ������ɹ�����true�� ���򷵻�false
			holesCount = tMesh.boundaries();
		}

		debugWriteMesh("meshOut", tMesh);
		debugDisp("holesCount == ", holesCount);
		debugDisp("flagDeg == ", flagDeg);
		debugDisp("finished.");
	}

}