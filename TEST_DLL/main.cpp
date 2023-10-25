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

	void step3() 
	{
		// 1. 读取上下底面拟合曲线，生成胚子网格；
		const int layersCount = 3;
		const float maxTriArea = 1.0;
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

		// for debug
		{
			Eigen::MatrixXd upperVers, lowerVers;
			Eigen::MatrixXi upperTris, lowerTris;
			circuit2mesh(upperVers, upperTris, versLoopUpper, topNorm, maxTriArea);
			circuit2mesh(lowerVers, lowerTris, versLoopLower, btmNorm, maxTriArea);
			debugWriteMesh("meshUpper", upperVers, upperTris);
			debugWriteMesh("meshLower", lowerVers, lowerTris);
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


}




int main(int argc, char** argv)
{
	// TEST_MYDLL::test66();

	TEST_JAW_PALATE::step3();

	debugDisp("main() finished.");

	return 0;
}