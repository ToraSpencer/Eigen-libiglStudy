#include "myEigen.h"
 

////////////////////////////////////////////////////////////////////////////////////////////// DEBUG 接口
static std::string g_debugPath = "E:/";

 
static void debugDisp()			// 递归终止
{						//		递归终止设为无参或者一个参数的情形都可以。
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

 

template<typename DerivedV>
static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteEdgesMat(path, edges, vers);
}



///////////////////////////////////////////////////////////////////////////////////////////////// 实现

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


// 读取点云OBJ文件中的数据，存储到齐次坐标系表示的点云矩阵中；注：顶点在矩阵中是列表示的，第四维的元素始终为1；
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

		// 顶点信息		
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


// 齐次坐标系点云数据写入到OBJ文件中；注：顶点在矩阵中是列表示的，第四维的元素始终为1；
void objWriteVerticesHomoMat(const char* fileName, const Eigen::MatrixXf& vers)
{
	std::ofstream dstFile(fileName);
	for (int i = 0; i < vers.cols(); i++)
	{
		dstFile << "v " << vers(0, i) << " " << vers(1, i) << " " << vers(2, i) << std::endl;
	}
	dstFile.close();

}


//	普通点云矩阵转换为齐次坐标系下的点云矩阵：
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


// 齐次坐标系下的点云矩阵变换为普通点云矩阵
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

	const float SR = 0.5;			// 空间采样率SR——相邻两个采样点的距离（单位mm）
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
	const float SR = 0.5;			// 空间采样率SR——相邻两个采样点的距离（单位mm）
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


// 布尔向量转化为索引向量；
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


// 索引向量转化为布尔向量；
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


// 多项式插值
void polyInterpolation() 
{

}


// 高斯插值
void gaussInterpolation() {}


// 最小二乘多项式拟合曲线：
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


// 测试myEigen中的接口
namespace TEST_MYEIGEN 
{
	// 测试编码解码：
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


	// 测试区域生长triangleGrow();
	void test1() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/材料/fatTeeth.obj");
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


	// 测试区域生长实现的重叠三角片检测：
	void test11()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/材料/jawMeshDense_algSimp_60000.obj");
		objWriteMeshMat("E:/triangleGrowInput.obj", vers, tris);

		std::vector<std::pair<int, int>> opTrisPairs;
		int olCount = findOverLapTris(opTrisPairs, vers, tris);


		debugDisp("finished.");
	}


	// 测试区域生长实现的网格三角片朝向矫正：
	void test111() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/材料/jawMeshDense_algSimp_60000.obj");
		objWriteMeshMat("E:/triangleGrowInput.obj", vers, tris);

		int corrCount = correctTriDirs(trisOut, vers, tris, 0);
		debugDisp("correctTriDirs == ", corrCount);
		debugWriteMesh("meshsOut", vers, trisOut);
		
		debugDisp("finished.");
	}


	//		测试triangleGrowOuterSurf提取包含非流行三角片的网格提取外层表面：
	void test1111() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/材料/tmpSickArrResult.obj");
		objWriteMeshMat("E:/triangleGrowInput.obj", vers, tris);
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();

		// 1. 确定所有非流形无向边：
		Eigen::MatrixXi nmnUedges; 
		std::vector<std::tuple<std::pair<int, int>, int, int>> nmnInfos;
		int nmnCount = nonManifoldUEs(nmnUedges, nmnInfos, tris);
		if (nmnCount < 0) 
		{
			debugDisp("error! nonManifoldUEs() run failed.");
			return;
		}
		if (nmnUedges.rows() > 0)
		{
			debugDisp("非流形无向边数量：nmnUedges.rows() == ", nmnUedges.rows());
			debugWriteEdges("nmnEdgesOrigin", nmnUedges, vers);
		}
 
		// 3. 一次提取外表面可能不能完全去除非流形边，需要多次调用：
		while (nmnCount > 0) 
		{
			debugDisp("当前非流形有向边数：", nmnCount);
			debugDisp("执行triangleGrowOuterSurf()...");
			if (!triangleGrowOuterSurf(versOut, trisOut, vers, tris, true))
				debugDisp("triangleGrowOuterSurf() failed!");
			debugDisp("去除三角片个数：", trisCount - trisOut.rows());

			// 检查输出网格中是否有非流形边：
			nmnUedges.resize(0, 0);
			nmnCount = nonManifoldUEs(nmnUedges, trisOut);
			if (nmnCount < 0)
			{
				debugDisp("error! nonManifoldUEs() run failed.");
				return;
			}
			if (nmnUedges.rows() > 0)
				debugWriteEdges("nmnUedges", nmnUedges, versOut);
			nmnCount = nmnUedges.rows();

			vers = versOut;
			tris = trisOut;
			trisCount = tris.rows();
		}

		debugWriteMesh("triangleGrowOuterSurf", vers, tris);
		std::cout << "finished." << std::endl;
	}


	// 检测重复三角片，重复顶点，
	void test2() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/材料/shrinkedMeshDirty.obj");

		Eigen::MatrixXi bdrys;
		bdryEdges(bdrys, tris);
		std::cout << "boundary edges count is " << bdrys.rows() << std::endl;
		objWriteEdgesMat("E:/bdrys.obj", bdrys, vers);

		std::vector<int> repIdxes;
		findRepTris(repIdxes, tris);


		std::cout << "finished." << std::endl;
	}


	// 测试获取三角网格基础属性的接口：
	void test3() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		bool retFlag = true;
		objReadMeshMat(vers, tris, "E:/材料/jawCore.obj");

		const unsigned versCount = vers.rows();
		const unsigned trisCount = tris.rows();

		// 计算网格每个三角片的面积：
		Eigen::VectorXd trisAreaVec;
		if (!trisArea(trisAreaVec, vers, tris))
			return;

		// 手动计算第一个三角片的面积：
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


	// 测试基本数学接口
	void test4()
	{
		dispPair(cart2polar(5.0, -5.0));
		dispPair(polar2cart(20.0, 3.14159/3));

		std::cout << "finished." << std::endl;
	}

#ifdef USE_TRIANGLE_H
	// 测试生成基础图形的接口：
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

		// 三角剖分：
		Eigen::MatrixXd surfVers, circleVers;
		Eigen::MatrixXi surfTris;
		objReadVerticesMat(circleVers, "E:/材料/circleVers.obj");
		circuit2mesh(surfVers, surfTris, circleVers);
		debugWriteMesh("surfMesh", surfVers, surfTris);

		// 生成柱体：
		Eigen::MatrixXf cylinderVers;
		Eigen::MatrixXi cylinderTris;
		genCylinder(cylinderVers, cylinderTris, axis, 10);				// 生成圆柱
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
#endif
	
	// 测试空间变换相关的接口：
	void test6()
	{
		Eigen::Matrix3d rotation;
		Eigen::RowVector3d axisArrow;
		Eigen::RowVector3d arrow1, arrow2, arrowRotted;
		float theta = pi;

		// 绕任意轴旋转：
		theta = pi / 3;
		axisArrow = Eigen::RowVector3d{20, 0, 0};
		rotation = getRotationMat(axisArrow, theta);
		arrowRotted = arrow1 * rotation.transpose();
		dispVec<double, 3>(arrowRotted);

		// 旋转到目标向量：
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


	// 查找洞和边缘曲线：
	void test7() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		std::vector<std::vector<int>> holes, bdrySegs;

		objReadMeshMat(vers, tris, "E:/材料/meshRepairInput_c_noOpTris.obj");
		debugWriteMesh("meshInput", vers, tris);

		findHolesBdrySegs(holes, bdrySegs, vers, tris);


		debugDisp("finished.");
	}


	// 网格IO:
	void test8() 
	{
		std::string file1{"E:/材料/tooth.obj"};							// 原文件是单精度
		std::string file2{ "E:/材料/meshRepairInput.obj" };			// 原文件是双精度
		Eigen::MatrixXd vers;
		Eigen::MatrixXf versF;
		Eigen::MatrixXi tris;
 
	}
}



// 测试图像处理：
namespace TEST_DIP
{
	void test0()
	{
		Eigen::MatrixXi mat1(9, 9);
		mat1.setRandom();
		
		Eigen::MatrixXd mask(5, 5);
		mask.setOnes();

		Eigen::MatrixXi matOut;
		linearSpatialFilter(matOut, mat1, mask);

		dispMat(mat1);
		dispMat(matOut);

		std::cout << "finished." << std::endl;
	}

}

 

