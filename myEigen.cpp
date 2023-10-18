#include "myEigen.h"
 

////////////////////////////////////////////////////////////////////////////////////////////// DEBUG 接口
namespace MY_DEBUG 
{
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


	template<typename DerivedV>
	static void debugWriteVers2D(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat2D(path, vers);
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
		objReadMeshMat(vers, tris, "E:/材料/meshWrongTriDir.obj");
		
		debugWriteMesh("meshInput", vers, tris);

		int corrCount = correctTriDirs(trisOut, vers, tris);
		debugDisp("correctTriDirs == ", corrCount);
		debugWriteMesh("meshsOut", vers, trisOut);
		
		debugDisp("finished.");
	}


	// 顶点区域生长提取单连通网格：
	void test1111() 
	{
		Eigen::MatrixXd vers, versTmp;
		Eigen::MatrixXi tris, trisTmp;
		objReadMeshMat(vers, tris, "E:/材料/arrangeResult.obj");

		std::vector<Eigen::MatrixXd> compVers;
		std::vector<Eigen::MatrixXi> compTris;
		simplyConnectedSplitMesh(compVers, compTris, vers, tris);

		for(int k = 0; k<compVers.size(); ++k)
			debugWriteMesh((std::string{ "compMesh_" } + std::to_string(k)).c_str(), compVers[k], compTris[k]);

		debugDisp("finished.");
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


// 测试myEigenModeling中的接口
namespace TEST_MYEIGEN_MODELING
{
#ifdef USE_TRIANGLE_H

	// 测试生成柱体
	void test1()
	{
		Eigen::MatrixXf axis(15, 3);
		Eigen::MatrixXd surfVers, circleVers;
		Eigen::MatrixXi surfTris;
		Eigen::MatrixXf cylinderVers;
		Eigen::MatrixXi cylinderTris;

		// 1. 生成柱体轴线；
		axis.setZero();
		double deltaTheta = pi / 10;
		for (unsigned i = 0; i < 15; ++i)
		{
			double theta = deltaTheta * i;
			axis(i, 0) = 50 * cos(theta);
			axis(i, 1) = 50 * sin(theta);
		}
		debugWriteVers("axis", axis);

		// 2. 输入回路顶点三角剖分，得到圆形底面网格：
		objReadVerticesMat(circleVers, "E:/材料/circleVers.obj");
		circuit2mesh(surfVers, surfTris, circleVers);
		debugWriteMesh("surfMesh", surfVers, surfTris);

		// 3. 生成圆柱体：
		genCylinder(cylinderVers, cylinderTris, axis, 10);				// 生成圆柱
		debugWriteMesh("cylinder", cylinderVers, cylinderTris);

		// 4. 生成方柱体；
		axis.resize(0, 0);
		interpolateToLine(axis, Eigen::RowVector3f{ 0, 0, 0 }, Eigen::RowVector3f{ 0, 0, 5 }, 0.5);
		cylinderVers.resize(0, 0);
		cylinderTris.resize(0, 0);
		genCylinder(cylinderVers, cylinderTris, axis, std::make_pair(10.0, 20.0));
		debugWriteMesh("pillar", cylinderVers, cylinderTris);

		cylinderVers.resize(0, 0);
		cylinderTris.resize(0, 0);
		genAlignedCylinder(cylinderVers, cylinderTris, axis, std::make_pair(1.5, 1.5), 0.5);
		debugWriteMesh("AlignedPillar", cylinderVers, cylinderTris);

		// 5. 读取上下底面2D边界环路，生成柱体：
		Eigen::MatrixXd topLoop, btmLoop;
		cylinderVers.resize(0, 0);
		cylinderTris.resize(0, 0);
		axis.resize(0, 0);
		objReadVerticesMat(topLoop, "E:/颌板/bdryUpper_final.obj");
		objReadVerticesMat(btmLoop, "E:/颌板/bdryLower_final.obj");
		interpolateToLine(axis, Eigen::RowVector3f{ 0, 0, 0 }, Eigen::RowVector3f{ 0, 0, 5 }, 1.0);
		genCylinder(cylinderVers, cylinderTris, axis, topLoop, btmLoop, true);
		debugWriteMesh("pillar2", cylinderVers, cylinderTris);

		// 5. 读取上下底面3D边界环路，生成柱体：
		topLoop.resize(0, 0);
		btmLoop.resize(0, 0);
		cylinderVers.resize(0, 0);
		cylinderTris.resize(0, 0);
		axis.resize(0, 0);
		objReadVerticesMat(topLoop, "G:\\gitRepositories\\matlabCode\\颌板\\data/curveFitUpper.obj");
		objReadVerticesMat(btmLoop, "G:\\gitRepositories\\matlabCode\\颌板\\data/curveFitLower.obj");
		genCylinder(cylinderVers, cylinderTris, topLoop, btmLoop, 5, true);
		debugWriteMesh("pillar3", cylinderVers, cylinderTris);

		std::cout << "finished." << std::endl;
	}


#endif
	 
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
		bool ret = getLoopEdges(edges, versCount);
		debugWriteEdges("edges", edges, vers);

		debugDisp("finished.");
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

  
