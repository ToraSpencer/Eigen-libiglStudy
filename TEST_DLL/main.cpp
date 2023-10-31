#include "smarteeDLL.h"
#include "my_cgal_dll.h"

// 动态库的导入lib:
#pragma comment(lib, "smarteeDLL.lib")
#pragma comment(lib, "MY_CGAL_DLL.lib")

// myEigen静态库：
#include "myEigenIO/myEigenIO.h"
#pragma comment(lib, "myEigenIO.lib")
#include "myEigenBasicMath/myEigenBasicMath.h"
#pragma comment(lib, "myEigenBasicMath.lib")
#include "myEigenModeling/myEigenModeling.h"
#pragma comment(lib, "myEigenModeling.lib")
#include "myEigenPMP/myEigenPMP.h"
#pragma comment(lib, "myEigenPMP.lib")
#include "myEigenHomoCoor/myEigenHomoCoor.h"
#pragma comment(lib, "myEigenHomoCoor.lib")

// 已禁用vcpkg;


// ！！！注：在项目属性→调试→环境中给变量PATH赋值为DLL文件所在目录；
/*
	若不设置此项的话，程序只会在项目所在文件夹以及输出文件夹中搜索DLL文件；
	不仅是本程序测试的dll文件，其依赖的所有dll文件都要放入到可以被搜索到的目录中；
*/


/////////////////////////////////////////////////////////////////////////////////////////////////////////// DEBUG接口；
namespace MY_DEBUG
{
	static const std::string g_debugPath = "E:/";

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

	template<typename DerivedV>
	static void debugWriteVers(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		writeOBJ(path, vers);
	}

	template<typename DerivedV>
	static void debugWriteMesh(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers, \
		const Eigen::MatrixXi& tris)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		writeOBJ(path, vers, tris);
	}
}
using namespace MY_DEBUG;


/////////////////////////////////////////////////////////////////////////////////////////////////////////// 测试smarteeDLL.dll中的导出函数
namespace TEST_MYDLL
{
	// 测试封装的IO功能：
	void test0()
	{
		Eigen::MatrixXf versF1, versF2;
		Eigen::MatrixXd versD1, versD2;
		Eigen::MatrixXi tris21, tris22;

		readOBJ(versF1, "E:/材料/circleVers.obj");
		readOBJ(versD1, "E:/材料/circleVers.obj");
		readOBJ(versF2, tris21, "E:/材料/tooth.obj");
		readOBJ(versD2, tris22, "E:/材料/tooth.obj");

		debugWriteVers("versOut1", versD1);
		debugWriteMesh("meshOut1", versF2, tris21);
		debugWriteMesh("meshOut2", versD2, tris22);

		debugDisp("finished.");
	}


	// 测试老版B样条插值：
	void test1()
	{
		Eigen::MatrixXf versIn, versOut;
		Eigen::RowVector3f planeVer{ 0, 0, 0 };
		Eigen::RowVector3f planeNorm{ 0, 0, 1 };
		readOBJ(versIn, "E:/材料/circleVers.obj");
		BsplineInterp(versOut, versIn, planeVer, planeNorm);
		writeOBJ("E:/versBsplineInterp.obj", versOut);

		debugDisp("finished.");
	}


	// 测试新版三次B样条插值拟合曲线：
	void test2()
	{
		Eigen::MatrixXf versIn, versOut;
		readOBJ(versIn, "E:/颌板/sampleVersUpperBdry.obj");
		interpltnCubicBSpline(versOut, versIn, 300, true);

		debugWriteVers("versOut", versOut);
		debugDisp("finished.");
	}


	// 测试射线求交功能：
	void test3()
	{
		Eigen::MatrixXd srcVers, meshVers;
		Eigen::MatrixXi srcTris, meshTris;
		Eigen::RowVector3d dir{ 0, 0, -1 };
		meshRayOut result;
		readOBJ(srcVers, srcTris, "E:/材料/fatTeeth.obj");
		readOBJ(meshVers, meshTris, "E:/材料/jawMesh.obj");
		debugWriteMesh("meshInput1", srcVers, srcTris);
		debugWriteMesh("meshInput2", meshVers, meshTris);

		meshRayIntersect(result, srcVers, dir, meshVers, meshTris);

		std::vector<Eigen::RowVector3d> tmpVec;
		tmpVec.reserve(srcVers.rows());
		for (int i = 0; i < srcVers.rows(); ++i)
		{
			if (result.rayLen[i].size() > 0)
			{
				double projDepth = result.rayLen[i][0];
				Eigen::RowVector3d tmpVer = srcVers.row(i) + projDepth * dir;
				tmpVec.push_back(tmpVer);
			}
		}
		tmpVec.shrink_to_fit();

		Eigen::MatrixXd positiveProVers(tmpVec.size(), 3);
		for (int i = 0; i < positiveProVers.rows(); ++i)
			positiveProVers.row(i) = tmpVec[i];

		debugWriteVers("投影点", positiveProVers);

		debugDisp("finished.");
	}


	// alpha shapes
	void test4() 
	{
		Eigen::MatrixXf versIn, versOut;
		readOBJ(versIn, "E:/颌板/2DbdryUpper.obj");
		alphaShapes2D(versOut, versIn, 50);
		debugWriteVers("versOut", versOut);
		debugDisp("finished.");
	}


	// SDF, marching cubes
	void test5()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		SDF_RESULT sdfResult;

		readOBJ(vers, tris, "E:/材料/tooth.obj");
		genSDF(sdfResult, vers, tris, 0.5, 3);
		marchingCubes(versOut, trisOut, sdfResult, 0.6);

		debugWriteMesh("meshOut", versOut, trisOut);
		debugDisp("finished.");
	}


	void test55() 
	{
		Eigen::MatrixXf vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		SDF_RESULT sdfResult;

		readSTL(vers, tris, "E:/颌板/zhangweiliang-S.stl");
		genSDF(sdfResult, vers, tris, 0.3, 3);
		marchingCubes(versOut, trisOut, sdfResult, 0.1);
		debugWriteMesh("zhangweiliang-S_0.1", versOut, trisOut);

		vers.resize(0, 0);
		versOut.resize(0, 0);
		tris.resize(0, 0);
		trisOut.resize(0, 0);
		readSTL(vers, tris, "E:/颌板/zhangweiliang-X.stl");
		genSDF(sdfResult, vers, tris, 0.3, 3);
		marchingCubes(versOut, trisOut, sdfResult, 0.1);
		debugWriteMesh("zhangweiliang-X_0.1", versOut, trisOut);

		debugDisp("finished.");		 
	}


	// 网格布尔运算
	void test6()
	{
		Eigen::MatrixXd versOut, vers1, vers2;
		Eigen::MatrixXi trisOut, tris1, tris2;
		readOBJ(vers1, tris1, "E:/材料/tooth.obj");
		readOBJ(vers2, tris2, "E:/材料/cylinder1.obj");
		meshCross(versOut, trisOut, vers1, tris1, vers2, tris2);

		debugWriteMesh("meshOut", versOut, trisOut);

		debugDisp("finished.");
	}


	void test66() 
	{
		Eigen::MatrixXd versOut, vers1, vers2;
		Eigen::MatrixXi trisOut, tris1, tris2;

		readOBJ(vers1, tris1, "E:/颌板/palate.obj");
		readOBJ(vers2, tris2, "E:/颌板/zhangweiliang-S_0.1.obj");
		meshCross(versOut, trisOut, vers1, tris1, vers2, tris2);
		debugWriteMesh("meshCross1", versOut, trisOut);

		vers1 = versOut;
		tris1 = trisOut;
		versOut.resize(0, 0);
		trisOut.resize(0, 0);
		vers2.resize(0, 0);
		tris2.resize(0, 0);
		readOBJ(vers2, tris2, "E:/颌板/zhangweiliang-X_0.1.obj");
		meshCross(versOut, trisOut, vers1, tris1, vers2, tris2);
		debugWriteMesh("meshCross2", versOut, trisOut);

		debugDisp("finished.");
	}
}



namespace TEST_JAW_PALATE 
{
	// step3——生成胚子网格然后刻牙印，侧面使用genCylinder()生成 
	void test0()
	{
		// 1. 读取上下底面拟合曲线，生成胚子网格；
		const int layersCount = 3;
		const float maxTriArea = 3.0;
		Eigen::MatrixXf versJawUpper, versJawLower;
		Eigen::MatrixXi trisJawUpper, trisJawLower;
		Eigen::MatrixXd versLoopUpper, versLoopLower, versPalate;
		Eigen::MatrixXi trisPalate;
		Eigen::RowVector3f upperPlaneNorm{ 0, 0, -1 };
		Eigen::RowVector3f lowerPlaneNorm{ -0.0987428597222625, -0.123428574652828, 0.987428597222625 };
		Eigen::RowVector3f topNorm = -upperPlaneNorm;
		Eigen::RowVector3f btmNorm = -lowerPlaneNorm;
		readSTL(versJawUpper, trisJawUpper, "E:/颌板/颌板示例/2/upper.stl");
		readSTL(versJawLower, trisJawLower, "E:/颌板/颌板示例/2/lower.stl");
		readOBJ(versLoopUpper, "G:\\gitRepositories\\matlabCode\\颌板\\data\\curveFitUpper2.obj");
		readOBJ(versLoopLower, "G:\\gitRepositories\\matlabCode\\颌板\\data\\curveFitLower2.obj");
		genCylinder(versPalate, trisPalate, versLoopUpper, versLoopLower,\
				topNorm, btmNorm, layersCount, maxTriArea, true);
		debugWriteMesh("meshPalate", versPalate, trisPalate);


		// for debug——生成不插点三角剖分的老版本：
		{
			Eigen::MatrixXd versTmp;
			Eigen::MatrixXi trisTmp;
			circuit2mesh(versTmp, trisTmp, versLoopUpper);
			debugWriteMesh("meshUpperOld", versTmp, trisTmp);
		}

		// for debug——输出上下底面网格
		{
			Eigen::MatrixXd versUpper, versLower;
			Eigen::MatrixXi trisUpper, trisLower;
			circuit2mesh(versUpper, trisUpper, versLoopUpper, topNorm, maxTriArea);
			circuit2mesh(versLower, trisLower, versLoopLower, btmNorm, maxTriArea);
			debugWriteMesh("meshUpper", versUpper, trisUpper);
			debugWriteMesh("meshLower", versLower, trisLower);
		}


		// 2. 读取颌网格，膨胀0.1mm, meshcross
		{
			Eigen::MatrixXd versCrossed;
			Eigen::MatrixXf versJawNew;
			Eigen::MatrixXi trisJawNew, trisCrossed;
			SDF_RESULT sdfResult;
			genSDF(sdfResult, versJawUpper, trisJawUpper, 0.3, 3);
			marchingCubes(versJawNew, trisJawNew, sdfResult, 0.1);
			meshCross(versCrossed, trisCrossed, versPalate, trisPalate, versJawNew, trisJawNew);

			versPalate = versCrossed;
			trisPalate = trisCrossed;
			debugWriteMesh("meshCrossed1", versCrossed, trisCrossed);
		}

		{
			Eigen::MatrixXd versCrossed;
			Eigen::MatrixXf versJawNew;
			Eigen::MatrixXi trisJawNew, trisCrossed;
			SDF_RESULT sdfResult;
			genSDF(sdfResult, versJawLower, trisJawLower, 0.3, 3);
			marchingCubes(versJawNew, trisJawNew, sdfResult, 0.1);
			meshCross(versCrossed, trisCrossed, versPalate, trisPalate, versJawNew, trisJawNew);
			 
			debugWriteMesh("meshOut", versCrossed, trisCrossed);
		}


		debugDisp("finished.");
	}


	// step3——生成胚子网格然后刻牙印，侧面使用圆环三角剖分生成
	void test00()
	{
		// 0
		const int layersCount = 3;
		const float maxTriArea = 3.0;
		Eigen::MatrixXf versJawUpper, versJawLower;
		Eigen::MatrixXi trisJawUpper, trisJawLower;
		Eigen::MatrixXd versLoopUpper, versLoopLower;
		Eigen::RowVector3f upperPlaneNorm{ 0, 0, -1 };
		Eigen::RowVector3f lowerPlaneNorm{ -0.0987428597222625, -0.123428574652828, 0.987428597222625 };
		Eigen::RowVector3f topNorm = -upperPlaneNorm;
		Eigen::RowVector3f btmNorm = -lowerPlaneNorm;
		readSTL(versJawUpper, trisJawUpper, "E:/颌板/颌板示例/2/upper.stl");
		readSTL(versJawLower, trisJawLower, "E:/颌板/颌板示例/2/lower.stl");
		readOBJ(versLoopUpper, "G:\\gitRepositories\\matlabCode\\颌板\\data\\curveFitUpper2.obj");
		readOBJ(versLoopLower, "G:\\gitRepositories\\matlabCode\\颌板\\data\\curveFitLower2.obj");
		const int versCount1 = versLoopUpper.rows();
		const int versCount2 = versLoopUpper.rows();

		// 1. 生成上下底面网格；
		Eigen::MatrixXd versUpper, versLower;
		Eigen::MatrixXi trisUpper, trisLower;
		{
			circuit2mesh(versUpper, trisUpper, versLoopUpper, topNorm, maxTriArea);
			circuit2mesh(versLower, trisLower, versLoopLower, btmNorm, maxTriArea);
			debugWriteMesh("meshUpper", versUpper, trisUpper);
			debugWriteMesh("meshLower", versLower, trisLower);
		}

		// 2. 生成侧面网格：
		Eigen::MatrixXd versSide;
		Eigen::MatrixXi trisSide;
		{
			Eigen::MatrixXd versInner, versOuter, versTG, tmpVers;
			Eigen::MatrixXi trisTG;
			Eigen::VectorXi bdryOuter, bdryInner;

			getCircleVers(versInner, 10, versCount1);
			getCircleVers(versOuter, 20, versCount2);
			matInsertRows(tmpVers, versInner);
			matInsertRows(tmpVers, versOuter);
			debugWriteVers("versInput", tmpVers);

			Eigen::RowVector3d innerCenter = versInner.colwise().mean();
			bdryInner = Eigen::VectorXi::LinSpaced(versCount1, 0, versCount1 - 1);
			bdryOuter = Eigen::VectorXi::LinSpaced(versCount1, versCount1 + versCount2 - 1, versCount1);
			triangulateVers2Mesh(versTG, trisTG, tmpVers, std::vector<Eigen::VectorXi>{bdryInner, bdryOuter}, innerCenter);
			debugWriteMesh("meshTG", versTG, trisTG);

			matInsertRows(versSide, versLoopUpper);
			matInsertRows(versSide, versLoopLower);
			trisSide = trisTG;
			debugWriteMesh("meshSide", versSide, trisSide);
		}


		// 3. 生成整体的网格：
		Eigen::MatrixXd versPalate;
		Eigen::MatrixXi trisPalate;
		{
			const int versUpperCount = versUpper.rows();
			const int versLowerCount = versLower.rows();
			matInsertRows(versPalate, versUpper);
			matInsertRows(versPalate, versLower);
			trisLower.array() += versUpperCount;

			// 侧面网格三角片中，大于等于versCount1的顶点索引都要施加偏移量：
			const int offset = versUpperCount - versCount1;
			int* ptrInt = trisSide.data();
			for (int i = 0; i < trisSide.size(); ++i)
			{
				int& index = *ptrInt;
				if (index >= versCount1) 
					index += offset;	 
				ptrInt++;
			}

			matInsertRows(trisPalate, trisUpper);
			matInsertRows(trisPalate, trisLower);
			matInsertRows(trisPalate, trisSide);

			debugWriteMesh("meshPalate", versPalate, trisPalate);
		}

		// 4. 读取颌网格，膨胀0.1mm, meshcross
		{
			Eigen::MatrixXd versCrossed;
			Eigen::MatrixXf versJawNew;
			Eigen::MatrixXi trisJawNew, trisCrossed;
			SDF_RESULT sdfResult;
			genSDF(sdfResult, versJawUpper, trisJawUpper, 0.3, 3);
			marchingCubes(versJawNew, trisJawNew, sdfResult, 0.1);
			meshCross(versCrossed, trisCrossed, versPalate, trisPalate, versJawNew, trisJawNew);

			versPalate = versCrossed;
			trisPalate = trisCrossed;
			debugWriteMesh("meshCrossed1", versCrossed, trisCrossed);
		}
		{
			Eigen::MatrixXd versCrossed;
			Eigen::MatrixXf versJawNew;
			Eigen::MatrixXi trisJawNew, trisCrossed;
			SDF_RESULT sdfResult;
			genSDF(sdfResult, versJawLower, trisJawLower, 0.3, 3);
			marchingCubes(versJawNew, trisJawNew, sdfResult, 0.1);
			meshCross(versCrossed, trisCrossed, versPalate, trisPalate, versJawNew, trisJawNew);

			debugWriteMesh("meshOut", versCrossed, trisCrossed);
		}


		debugDisp("finished.");
	}
}


// 测试triangle库的三角剖分：
namespace TEST_TRIANGULATION 
{
	 
	// 测试点云三角剖分成网格；
	void test1()
	{
		Eigen::MatrixXd versIn;
		Eigen::MatrixXf versOut;
		Eigen::MatrixXi trisIn, trisOut;

		// 1. 点云剖分成不带洞的网格：
		readOBJ(versIn, "E:/材料/circleVers.obj");
		const int versCount = versIn.rows();
		Eigen::VectorXi loopVec = Eigen::VectorXi::LinSpaced(versCount, 0, versCount - 1);
		triangulateVers2Mesh(versOut, trisOut, versIn, std::vector<Eigen::VectorXi>{loopVec});
		debugWriteVers("versInput1", versIn);
		debugWriteMesh("meshOut1", versOut, trisOut);

		// 2. 点云剖分成带洞的网格：

		//			2.1 ...
		Eigen::MatrixXd vers0, vers10, vers11, vers12, vers2, allVers, holeCenters;
		Eigen::RowVector3d arrow10, arrow11, arrow12, arrow2, holeCenter10, holeCenter11, holeCenter12;
		const float radius0 = 5;
		const float radius1 = 1;
		const float radius2 = 2;
		getCircleVers(vers0, radius0, 100);
		getCircleVers(vers10, radius1, 20);
		getCircleVers(vers2, radius2, 30);
		vers11 = vers10;
		vers12 = vers10;
		arrow10 << 0, 2, 0;
		arrow11 << -2, 0, 0;
		arrow12 << 3, 0, 0;
		arrow2 << 10, 0, 0;
		vers10.rowwise() += arrow10;
		vers11.rowwise() += arrow11;
		vers12.rowwise() += arrow12;
		vers2.rowwise() += arrow2;
		holeCenter10 = vers10.colwise().mean();
		holeCenter11 = vers11.colwise().mean();
		holeCenter12 = vers12.colwise().mean();
		matInsertRows(allVers, vers0);
		matInsertRows(allVers, vers10);
		matInsertRows(allVers, vers11);
		matInsertRows(allVers, vers12);
		matInsertRows(allVers, vers2);
		matInsertRows(holeCenters, holeCenter10);
		matInsertRows(holeCenters, holeCenter11);
		matInsertRows(holeCenters, holeCenter12);
		debugWriteVers("allVers", allVers);

		//		2.2. 生成边缘信息：
		const int versCount0 = vers0.rows();
		const int versCount10 = vers10.rows();
		const int versCount11 = vers11.rows();
		const int versCount12 = vers12.rows();
		const int versCount2 = vers2.rows();
		std::vector<int> loops0(versCount0), loops10(versCount10), loops11(versCount11), loops12(versCount12), loops2(versCount2);
		int startIdx = 0;
		for (int i = 0; i < versCount0; ++i)
			loops0[i] = startIdx + i;
		startIdx += versCount0;
		for (int i = 0; i < versCount10; ++i)
			loops10[i] = startIdx + i;
		startIdx += versCount10;
		for (int i = 0; i < versCount11; ++i)
			loops11[i] = startIdx + i;
		startIdx += versCount11;
		for (int i = 0; i < versCount12; ++i)
			loops12[i] = startIdx + i;
		startIdx += versCount12;
		for (int i = 0; i < versCount2; ++i)
			loops2[i] = startIdx + i;
		startIdx += versCount2;

		//   
		versIn = allVers;
		std::vector<Eigen::VectorXi> loops{ vec2EigenVec(loops0), vec2EigenVec(loops10), \
				vec2EigenVec(loops11), vec2EigenVec(loops12), vec2EigenVec(loops2) };
		triangulateVers2Mesh(versOut, trisOut, versIn, loops, holeCenters);
		debugWriteVers("versInput2", versIn);
		debugWriteMesh("meshOut2", versOut, trisOut);

		debugDisp("main() finished.");

	}


	// 测试网格三角剖分成质量更好的网格：
	void test2() 
	{
		Eigen::MatrixXd versIn, versOut;
		Eigen::MatrixXi edgesIn, trisIn, trisOut;

		readOBJ(versIn, trisIn, "E:/材料/roundSurfNoInterp.obj");		
		getEdges(edgesIn, trisIn);
		triangulateRefineMesh(versOut, trisOut, versIn, trisIn, edgesIn, "pra2.0");
		debugWriteMesh("meshOut1", versOut, trisOut);
 


		debugDisp("finished.");
	}


	// 测试上下底面环路三角剖分分别生成侧面、上底面、下底面网格：
	void test3()
	{
		const int versCount1 = 300;
		const int versCount2 = 300;
		Eigen::MatrixXd versInner, versOuter, allVers;
		Eigen::MatrixXd versLoopUpper, versLoopLower, versOut;
		Eigen::RowVector3f upperPlaneNorm{ 0, 0, -1 };
		Eigen::RowVector3f lowerPlaneNorm{ -0.0987428597222625, -0.123428574652828, 0.987428597222625 };
		Eigen::RowVector3f topNorm = -upperPlaneNorm;
		Eigen::RowVector3f btmNorm = -lowerPlaneNorm;
		objReadVerticesMat(versLoopUpper, "E:/颌板/curveFitUpper2.obj");
		objReadVerticesMat(versLoopLower, "E:/颌板/curveFitLower2.obj");
	
		getCircleVers(versInner, 10, versCount1);
		getCircleVers(versOuter, 20, versCount2);
		matInsertRows(allVers, versInner);
		matInsertRows(allVers, versOuter);
		debugWriteVers("versInput", allVers);

		// 1. 三角剖分生成侧面：
		Eigen::MatrixXd versTG;
		Eigen::MatrixXi trisTG;
		Eigen::VectorXi bdryOuter, bdryInner;
		Eigen::RowVector3d innerCenter = versInner.colwise().mean();
		bdryInner = Eigen::VectorXi::LinSpaced(versCount1, 0, versCount1-1);
		bdryOuter = Eigen::VectorXi::LinSpaced(versCount1, versCount1 + versCount2 - 1, versCount1);
		triangulateVers2Mesh(versTG, trisTG, allVers, std::vector<Eigen::VectorXi>{bdryInner, bdryOuter}, innerCenter);
		debugWriteMesh("meshTG", versTG, trisTG);

		// 2. 三角剖分生成上下底面：
		Eigen::MatrixXd tmpVersOut;
		Eigen::MatrixXi trisUpper, trisLower;
		Eigen::MatrixXd versLoopUpperProj, versLoopLowerProj;
		versLoopUpperProj = versLoopUpper;
		versLoopLowerProj = versLoopLower;

		versLoopUpperProj.col(2).setZero(); 
		versLoopLowerProj.col(2).setZero();

		circuit2mesh(tmpVersOut, trisUpper, versLoopUpperProj);
		circuit2mesh(tmpVersOut, trisLower, versLoopLowerProj);
		Eigen::VectorXi tmpVec = trisLower.col(2);
		trisLower.col(2) = trisLower.col(1);
		trisLower.col(1) = tmpVec;

		debugWriteMesh("meshUpper", versLoopUpper, trisUpper);
		debugWriteMesh("meshLower", versLoopLower, trisLower);

		// 3. 生成最终网格：
		matInsertRows(versOut, versLoopUpper);
		matInsertRows(versOut, versLoopLower);

		//		3.1 
		Eigen::MatrixXi trisOut = trisTG;
		trisLower.array() += versCount1;
		matInsertRows(trisOut, trisUpper);
		matInsertRows(trisOut, trisLower);

		debugWriteMesh("meshOut", versOut, trisOut);



		debugDisp("finished.");
	}


	// 参数化之后的侧面网格重新三角剖分成质量更高的网格；
	void test4() 
	{
		Eigen::MatrixXd versIn, versOut, versInner;
		Eigen::MatrixXi edgesIn, trisIn, trisOut;
		Eigen::VectorXi bdryOuter, bdryInner;
		Eigen::RowVector3d innerCenter;
		const int versCount1 = 300;
		const int versCount2 = 300;

		readOBJ(versIn, trisIn, "E:/颌板/mesh_side_ini_param_2.obj");
		bdryInner = Eigen::VectorXi::LinSpaced(versCount1, 0, versCount1 - 1);
		bdryOuter = Eigen::VectorXi::LinSpaced(versCount1, versCount1 + versCount2 - 1, versCount1);
		subFromFlagVec(versInner, versIn, bdryInner);
		innerCenter = versInner.colwise().mean();

		triangulateVers2Mesh(versOut, trisOut, versIn, std::vector<Eigen::VectorXi>{bdryInner, bdryOuter}, innerCenter, "pra3.0");
		debugWriteMesh("mesh_side_ini_param_refined_2", versOut, trisOut);

		debugDisp("finished.");
	}

}



int main(int argc, char** argv)
{
	TEST_JAW_PALATE::test00();

	return 0;
}