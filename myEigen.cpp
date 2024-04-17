#include "myEigen.h"
 

////////////////////////////////////////////////////////////////////////////////////////////// DEBUG 接口
namespace MY_DEBUG
{
	static std::string g_debugPath = "E:/";

	// lambda——打印std::cout支持的类型变量。
	template <typename T>
	static auto disp = [](const T& arg)
	{
		std::cout << arg << ", ";
	};

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


	template <typename Derived>
	static void dispData(const Eigen::MatrixBase<Derived>& m)
	{
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
	static void debugWriteVers(const char* name, const Eigen::MatrixBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat(path, vers);
	}


	template<typename DerivedV>
	static void debugWriteVers2D(const char* name, const Eigen::MatrixBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat2D(path, vers);
	}


	template<typename DerivedV>
	static void debugWriteMesh(const char* name, \
		const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteMeshMat(path, vers, tris);
	}


	template<typename DerivedV>
	static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, \
		const Eigen::MatrixBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteEdgesMat(path, edges, vers);
	}

}
using namespace MY_DEBUG;


 
///////////////////////////////////////////////////////////////////////////////////////////////// 测试函数

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
		getRotationMat(rotation, axisArrow, theta);
		arrowRotted = arrow1 * rotation.transpose();
		dispVec<double, 3>(arrowRotted);

		// 旋转到目标向量：
		arrow1.setRandom();
		arrow2.setRandom();
		arrow1.normalize();
		arrow2.normalize();
		dispVec<double, 3>(arrow1);
		dispVec<double, 3>(arrow2);
		getRotationMat(rotation, arrow1, arrow2);
		arrowRotted = arrow1 * rotation.transpose();
		dispVec<double, 3>(arrowRotted);

		std::cout << "finished." << std::endl;
	}


	// 找洞、补洞；
	void test7() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		std::vector<std::vector<int>> holes;

		objReadMeshMat(vers, tris, "E:/材料/meshHoles.obj");
		debugWriteMesh("meshInput", vers, tris);
 
		// 1. 找洞：
		if (findHoles(holes, vers, tris) < 0)
		{
			debugDisp("error!!! findHoles() failed!");
			return;
		}
		
		// 2. 补洞：
		Eigen::MatrixXi newTris;
		if (!fillSmallHoles(newTris, holes))
		{
			debugDisp("error!!! fillSmallHoles() failed!");
			return;
		}
		matInsertRows(tris, newTris);
		debugWriteMesh("meshHoleFilled", vers, tris);
		 
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


	//		测试triangleGrowOuterSurf提取包含非流行三角片的网格提取外层表面：
	void test9()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/材料/tmpArrangeResult.obj");
		objWriteMeshMat("E:/triangleGrowInput.obj", vers, tris);
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();

		// 1. 确定所有非流形无向边：
		Eigen::MatrixXi nmnEdges;
		int nmnCount = nonManifoldEdges(nmnEdges, tris);
		if (nmnCount < 0)
		{
			debugDisp("error! nonManifoldUEs() run failed.");
			return;
		}
		if (nmnEdges.rows() > 0)
		{
			debugDisp("非流形无向边数量：nmnEdges.rows() == ", nmnEdges.rows());
			debugWriteEdges("nmnEdgesOrigin", nmnEdges, vers);
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
			nmnEdges.resize(0, 0);
			nmnCount = nonManifoldEdges(nmnEdges, trisOut);
			if (nmnCount < 0)
			{
				debugDisp("error! nonManifoldUEs() run failed.");
				return;
			}
			if (nmnEdges.rows() > 0)
				debugWriteEdges("nmnEdges", nmnEdges, versOut);
			nmnCount = nmnEdges.rows();

			vers = versOut;
			tris = trisOut;
			trisCount = tris.rows();
		}

		debugWriteMesh("triangleGrowOuterSurf", vers, tris);
		std::cout << "finished." << std::endl;
	}


	//		测试triangleGrowOuterSurf()——从指定三角片开始生长，只调用一次；
	void test99()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/材料/tmpArrangeResult.obj");
		objWriteMeshMat("E:/triangleGrowInput.obj", vers, tris);
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		int startIdx = 17801;

		// 1. 确定所有非流形无向边：
		Eigen::MatrixXi nmnEdges;
		int nmnCount = nonManifoldEdges(nmnEdges, tris);
		if (nmnCount < 0)
		{
			debugDisp("error! nonManifoldUEs() run failed.");
			return;
		}
		if (nmnEdges.rows() > 0)
		{
			debugDisp("非流形无向边数量：nmnEdges.rows() == ", nmnEdges.rows());
			debugWriteEdges("nmnEdgesOrigin", nmnEdges, vers);
		}

		// 3.  
		debugDisp("当前非流形有向边数：", nmnCount);
		debugDisp("执行triangleGrowOuterSurf()...");
		if (!triangleGrowOuterSurf(versOut, trisOut, vers, tris, true))
			debugDisp("triangleGrowOuterSurf() failed!");
		debugDisp("去除三角片个数：", trisCount - trisOut.rows());

		// 检查输出网格中是否有非流形边：
		nmnEdges.resize(0, 0);
		nmnCount = nonManifoldEdges(nmnEdges, trisOut);
		if (nmnCount < 0)
		{
			debugDisp("error! nonManifoldUEs() run failed.");
			return;
		}
		if (nmnEdges.rows() > 0)
			debugWriteEdges("nmnEdges", nmnEdges, versOut);
		nmnCount = nmnEdges.rows(); 
 

		debugWriteMesh("triangleGrowOuterSurf", versOut, trisOut);
		std::cout << "finished." << std::endl;
	}


	// 轴向包围盒类Eigen::AlignedBox
	void test10() 
	{
		Eigen::MatrixXd vers, vers1, versBox;
		Eigen::MatrixXi tris, tris1, trisBox;
		objReadMeshMat(vers, tris, "E:/材料/tooth.obj");
		objReadMeshMat(vers1, tris1, "E:/材料/jawMesh.obj");
		const int versCount = vers.rows();
		const int versCount1 = vers1.rows();
		const int trisCount = tris.rows();
		const int trisCount1 = tris1.rows();

		// extend()——包围盒扩张至将列向量表示的顶点纳入其中；
		Eigen::AlignedBox3d aabb;
		for (int i = 0; i < versCount; ++i) 
		{ 
			Eigen::Vector3d v = vers.row(i).transpose();
			aabb.extend(v);			
		}
		genAABBmesh(versBox, trisBox, aabb);
		debugWriteMesh("boxMesh", versBox, trisBox);

		// max(), min()——返回包围盒最大、最小向量；
		Eigen::RowVector3d arrowMax = aabb.max().transpose();				// 注意返回的是列向量；
		Eigen::RowVector3d arrowMin = aabb.min().transpose();
		debugWriteVers("arrowMax", arrowMax);
		debugWriteVers("arrowMin", arrowMin);
		double boxLen = (arrowMax - arrowMin).norm();
		debugDisp("boxLen == ", boxLen); 

		// setEmpty()——初始化AABB
		aabb.setEmpty();
		for (int i = 0; i < versCount1; ++i)
		{
			Eigen::Vector3d v = vers1.row(i).transpose();
			aabb.extend(v);
		}
		genAABBmesh(versBox, trisBox, aabb);
		debugWriteMesh("boxMesh1", versBox, trisBox);

		debugDisp("finished.");
	}
}


namespace TEST_MYEIGEN_IO
{
	void test0()
	{
		Eigen::MatrixXf versF1, versF2;
		Eigen::MatrixXd versD1, versD2;
		Eigen::MatrixXi tris1, tris2;

		objReadMeshMat(versF1, tris1, "E:/材料/tooth.obj");
		debugWriteMesh("meshInput", versF1, tris1);
		debugWriteMesh("meshInput2", 2.0f * versF1, tris1); 
		stlWriteMeshMat("E:/tooth.stl", versF1, tris1); 

		stlReadMeshMat(versD2, tris2, "E:/材料/jawMeshSimplified.stl");
		debugWriteMesh("meshInput22", versD2, tris2);

		debugDisp("test0() finished.");
	}
}



namespace TEST_MYEIGEN_BASIC_MATH
{
	// 测试矩阵的增删查改和运算相关的轮子：
	void test0()
	{
		// 1. 矩阵串联、并联
		debugDisp("1. 矩阵串联、并联");
		{
			Eigen::MatrixXi m1(4, 3);
			Eigen::MatrixXi m2{ Eigen::Matrix3i::Identity() };
			m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
			debugDisp("m1 == \n", m1, "\n");

			matInsertRows(m1, m2);									// 尾后插入矩阵，串联；
			debugDisp("m1 == \n", m1, "\n");

			m1 = m1.topRows(3).eval();							// 赋值的时候等号两个有同一个对象，最好用eval()来确保安全；
			matInsertRows(m1, Eigen::RowVector3i{ 999, 999, 999 });		// 尾后插入行向量
			matInsertRows(m1, -Eigen::MatrixXi::Ones(4, 3));				// 尾后插入return type对象；
			debugDisp("m1 == \n", m1, "\n");

			m1 = m1.topRows(4).eval();
			matInsertCols(m1, 88 * Eigen::Vector4i::Ones());		// 并联；
			debugDisp("m1 == \n", m1, "\n");
		}


		// 2. subFrom.... 
		debugDisp("2. subFrom...");
		{
			Eigen::MatrixXi m1(4, 3);
			Eigen::MatrixXi m2;
			std::vector<int> oldNewIdxInfo, newOldIdxInfo;
			m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
			debugDisp("m1 == \n", m1, "\n");

			subFromIdxVec(m2, m1, Eigen::Vector2i{0, 2});
			debugDisp("m2 == \n", m2, "\n");

			subFromIdxVec(m2, m1, std::vector{ 0, 2, 3 });
			debugDisp("m2 == \n", m2, "\n");

			subFromIdxCon(m2, m1, std::list{ 1, 2, 3 });
			debugDisp("m2 == \n", m2, "\n");

			subFromFlagVec(m2, m1, eigenVec2Vec(Eigen::Vector4i{0, 1, 0, 1}));
			debugDisp("m2 == \n", m2, "\n");

			subFromFlagVec(m2, oldNewIdxInfo, newOldIdxInfo, m1, eigenVec2Vec(Eigen::Vector4i{0, 1, 1, 1}));
			debugDisp("m2 == \n", m2, "\n");

			debugDisp("oldNewIdxInfo: ");
			traverseSTL(oldNewIdxInfo, disp<int>);
			debugDisp("newOldIdxInfo: ");
			traverseSTL(newOldIdxInfo, disp<int>);

			flagVec2oldNewIdxInfo(oldNewIdxInfo, newOldIdxInfo, eigenVec2Vec(Eigen::Vector4i{ 0, 1, 1, 1 }));
			debugDisp("oldNewIdxInfo: ");
			traverseSTL(oldNewIdxInfo, disp<int>);
			debugDisp("newOldIdxInfo: ");
			traverseSTL(newOldIdxInfo, disp<int>);
		}

		debugDisp("finished.");
	}


	// 测试最小二乘拟合：
	void test1() 
	{
		Eigen::MatrixXf versLoop;
		objReadVerticesMat(versLoop, "E:\\颌板\\颌板示例\\2\\curveFitLower.obj");
		debugWriteVers("versLoop", versLoop);

		Eigen::RowVector3f centerF = versLoop.colwise().mean();
		Eigen::MatrixXf versPlane;
		Eigen::MatrixXi trisPlane;
		Eigen::Vector4f x;
		Eigen::VectorXf beta;
		Eigen::RowVector3f pos0, norm;
		float a, b, c, d, x0, y0, z0, lambda;
		x0 = centerF(0);
		y0 = centerF(1);
		z0 = centerF(2);
		 
		// 拟合平面：
		std::vector<float> coeffs;
		versFitPlane(coeffs, versLoop);
		a = coeffs[0];
		b = coeffs[1];
		c = coeffs[2];
		d = coeffs[3];

		// 打印拟合出的平面；
		norm = Eigen::RowVector3f{ a, b, c };
		lambda = -d - centerF.dot(norm);
		pos0 = centerF + lambda * norm;
		genRoundSurfMesh(versPlane, trisPlane, pos0, norm, 30, 100);
		debugWriteMesh("plane", versPlane, trisPlane);

		debugDisp("finished.");
	}


	// 测试线性变换：
	void test2() 
	{
		Eigen::MatrixXd versGlobal, versLocal, versResult;
		Eigen::MatrixXi trisGlobal, trisLocal;
		Eigen::Matrix4d affine;
		objReadMeshMat(versGlobal, trisGlobal, "E:/托槽环/upper_comps_fixed/托槽21.obj");
		objReadMeshMat(versLocal, trisLocal, "E:/托槽环/upper_comps_fixed/local托槽21.obj");

		if (!getAffineL2G(affine, versLocal, versGlobal))
		{
			debugDisp("error!!! getAffineL2G() failed.");
			return;
		}
		versResult = homoVers2VersD(affine * vers2HomoVersD(versLocal));
		debugWriteMesh("global", versGlobal, trisGlobal);
		debugWriteMesh("local", versLocal, trisLocal);
		debugWriteMesh("result", versResult, trisGlobal);

		debugDisp("test2 finished.");

	}
}



// 测试myEigenPMP中的接口
namespace TEST_MYEIGEN_PMP
{
	void test0()
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi edges;
		objReadVerticesMat(vers, "E:/材料/loop2D.obj"); 
		debugWriteVers("versInput", vers);

		const int versCount = vers.rows();
		bool ret = getSortedLoopEdges(edges, versCount);
		debugWriteEdges("edges", edges, vers);

		debugDisp("finished.");
	}


	// 测试laplace光顺：
	void test1() 
	{
		Eigen::MatrixXd versIn, versOut;
		Eigen::MatrixXi trisIn;
		objReadMeshMat(versIn, trisIn, "E:/材料/surfMeshNotSmooth.obj");
		debugWriteMesh("meshIn", versIn, trisIn);

		const float lambda = 0.1;
		const int iterCount = 50;
		std::vector<int> fixedVerIdxes(300);			// 前300个点是边缘点；
		for (int i = 0; i < 300; ++i)
			fixedVerIdxes[i] = i;
		bool ret = laplaceFaring(versOut, versIn, trisIn, lambda,  iterCount, fixedVerIdxes);
		if (!ret)
		{
			debugDisp("error!!! laplace faring failed.");
			return;
		}

		debugWriteMesh("meshOut", versOut, trisIn);

		debugDisp("finished.");
	}


	// 测试对回路点云的处理：
	void test2() 
	{
#if 0
		Eigen::MatrixXf versLoop;
		objReadVerticesMat(versLoop, "E:/材料/LoopBdryUnsorted.obj");
		debugWriteVers("versInput", versLoop);

		// 1. 测试回路有序化：
		sortLoop2D(versLoop);
		debugWriteVers("loopSorted", versLoop);

		// 2. 测试回路均匀化：
		Eigen::MatrixXd versLoopArr;
		arrangeLoop(versLoopArr, versLoop, 500);
		debugWriteVers("loopArr", versLoopArr);
#else
		Eigen::MatrixXf versInput, versOut;
		objReadVerticesMat(versInput, "G:/合板测试结果/失败/1/bdrySortedAS.obj");
		debugWriteVers("versInput", versInput);

		arrangeLoop(versOut, versInput, 500);
		debugWriteVers("versOut", versOut);

#endif


		debugDisp("finished.");
	}


	// 测试计算网格的领域性质
	void test3() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/材料/tooth.obj");
		debugWriteMesh("meshInput", vers, tris);

		// get1ring()——得到网格顶点1领域的信息；
		const int index0 = 0;
		std::vector<std::unordered_set<int>> vIdx1ring, tIdx1ring;
		Eigen::MatrixXf vers1ring;
		Eigen::MatrixXi tris1ring;
		get1ring(vIdx1ring, tIdx1ring, vers, tris);
		subFromIdxCon(vers1ring, vers, vIdx1ring[index0]);
		subFromIdxCon(tris1ring, tris, tIdx1ring[index0]);
		debugWriteVers("ver0", vers.row(0));
		debugWriteVers("vers1ring0", vers1ring);
		debugWriteMesh("tris1ring0", vers, tris1ring);


		debugDisp("finished.");
	}


	// 测试区域生长算法：
	void test4() 
	{
		Eigen::MatrixXf vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/材料/meshArrangeResult.obj");
		debugWriteMesh("meshInput", vers, tris);

		// 1. 区域生长提取 最外层网格；
		Eigen::MatrixXd versOutter;
		Eigen::MatrixXi trisOutter;
		robustTriangleGrowOuterSurf(versOutter, trisOutter, vers, tris);
		debugWriteMesh("meshOutter", versOutter, trisOutter);

		// 2. 测试区域生长实现的网格三角片朝向矫正：
#if 0				// 当前有错误；
		objReadMeshMat(vers, tris, "E:/材料/meshWrongTriDir.obj");
		int corrCount = correctTriDirs(trisOut, vers, tris);
		debugDisp("correctTriDirs == ", corrCount);
		debugWriteMesh("meshsOut", vers, trisOut);
#endif

		// 3. 区域生长逐个分解出单连通网格
		Eigen::MatrixXd versTmp;
		Eigen::MatrixXi trisTmp;
		std::vector<Eigen::MatrixXd> compVers;
		std::vector<Eigen::MatrixXi> compTris;
		objReadMeshMat(vers, tris, "E:/材料/threeTeeth.obj");
		simplyConnectedSplitMesh(compVers, compTris, vers, tris);
		for (int k = 0; k < compVers.size(); ++k)
			debugWriteMesh((std::string{ "compMesh_" } + std::to_string(k)).c_str(), compVers[k], compTris[k]); 

		debugDisp("finished.");
	} 
	    

	// 测试获取网格属性：
	void test5() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/材料/holeMesh.obj"); 
		debugWriteMesh("meshInput", vers, tris);

		// 1. 获取网格边缘：
		Eigen::MatrixXi bdrys;
		bdryEdges(bdrys, tris);
		debugWriteEdges("bdrys", bdrys, vers);

		debugDisp("finished.");
	}


	// 有向边的区域生长：
	/*
		bool edgeGrow(
			std::vector<std::vector<int>>& loops,									
			const Eigen::MatrixBase<DerivedI>& edges
			) 
	*/
	template <typename IndexType, typename DerivedI>
	bool edgeGrow(std::vector<std::vector<IndexType>>& loops, \
		const Eigen::MatrixBase<DerivedI>& edges) 
	{   
		assert(edges.rows() > 0, "Assert!!! input edges is a empty matrix.");
		assert(edges.cols() == 2, "Assert!!! input edges is not a edge matrix.");
		loops.clear();

		// 1. 通过有向边得到邻接矩阵：
		Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
		edges2adjMat(adjSM_eCount, adjSM_eIdx, edges);
		Eigen::SparseMatrix<int>& adjSM = adjSM_eCount;

		// 2. 检测是否有非流形边（即重复的有向边）存在；
		bool blNmnEdge = false;
		traverseSparseMatrix(adjSM, [&adjSM, &blNmnEdge](auto& iter) 
			{
				if (iter.value() > 1)
				{
					blNmnEdge = true;
					return;
				}
			});
		if (blNmnEdge)
		{
			// debugDisp("输入边数据中存在非流形边");
			return false;
		}

		// 3. 建立有向边数据的顶点索引字典； 
		std::unordered_set<IndexType> isctIdxes;				 
		for (int i = 0; i < edges.rows(); ++i)
		{
			isctIdxes.insert(static_cast<IndexType>(edges(i, 0)));
			isctIdxes.insert(static_cast<IndexType>(edges(i, 1)));
		}

		// 4. 有向边生长的循环 
		while (!isctIdxes.empty())			  
		{
			std::vector<IndexType> currentLoop;						// 当前收集顶点的loop
			std::unordered_set<IndexType> auxiSet;
			IndexType tail = *isctIdxes.begin();
			IndexType head = -1;
			isctIdxes.erase(isctIdxes.begin());				// 收集到一个顶点则从字典中删除；
			auxiSet.insert(tail);
			currentLoop.push_back(tail);

			while (1)
			{
				bool blFindHead = false;
				for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM, tail); it; ++it)			// 列优先存储时，InnerIterator即列内迭代器；
				{
					head = static_cast<IndexType>(it.row());
					auto iter = auxiSet.insert(head);
					if (!iter.second)			 // 若插入失败，则跳过这条边，继续找下一条；否则存入该边起点，break；
						continue;
					else
					{
						blFindHead = true;
						currentLoop.push_back(head);
						isctIdxes.erase(head);
						break;
					}
				}

				if (!blFindHead)
					break;

				// 在邻接矩阵中删除本次循环中收集到的有向边：
				adjSM.prune([&](const Eigen::Index& row, const Eigen::Index& col, const float& value)->bool
					{
						if (head == static_cast<IndexType>(row) && \
							tail == static_cast<IndexType>(col))
							return false;
						else
							return true;            // 返回true的元素被保留；
					});

				// 以本次循环收集到的有向边的head顶点为下一次循环的tail顶点；
				tail = head;
			}

			loops.push_back(currentLoop);
		}

		return true;
	}


	// 测试有向边区域生长：
	void test6() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/材料/surfsTwo.obj");
		debugWriteMesh("meshInput", vers, tris);

		// 1. 找出所有边缘有向边：
		Eigen::MatrixXi bdrys;
		bdryEdges(bdrys, tris);
		debugWriteEdges("bdrys", bdrys, vers);

		// 2. 有向边生长：
		std::vector<std::vector<int>> loops;
		edgeGrow(loops, bdrys);

		// 3. 输出：
		for (int i = 0; i < loops.size(); ++i)
		{
			Eigen::MatrixXd loopVers;
			subFromIdxVec(loopVers, vers, loops[i]);
			char str[256];
			sprintf_s(str, "loop%d", i);
			debugWriteVers(str, loopVers);
		}

		debugDisp("test6() finished.");
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

  
