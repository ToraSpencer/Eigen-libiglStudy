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
	const Eigen::PlainObjectBase<DerivedVi>& versIn,\
	const std::vector<Eigen::VectorXi>& bdryLoops, \
	const Eigen::PlainObjectBase<DerivedVh>& holeCenters, \
	const char* strSwitcher = "pY")
{
	using ScalarO = typename DerivedVo::Scalar;

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

	// 2. 生成输入的洞的数据：
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

	// 3. 生成三角剖分信息：
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

		inputTrig.numberofholes = holesCount;					// 洞的数目
		inputTrig.holelist = holesCount > 0 ? versHole2Dtrans.data() : nullptr;		// 洞数据数组的首地址，用一个二维点表示一个洞，只要该点在洞内即可。

		inputTrig.numberofregions = 0;
		inputTrig.regionlist = nullptr;

		memset(&outputTrig, 0, sizeof(triangulateio));		// 初始化输出结构体。 
	}

	// 4. 执行三角剖分，拷贝结果：   
	char strCopy[512];
	strcpy(strCopy, strSwitcher);
	triangulate(strCopy, &inputTrig, &outputTrig, NULL);					// p——平面直线图，Y——不插点。  

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
	Eigen::MatrixXf versOutF(2, versOutCount);
	std::memcpy(versOutF.data(), reinterpret_cast<float*>(outputTrig.pointlist), sizeof(float) * 2 * versOutCount);
	versOutF.transposeInPlace();
	versOutF.conservativeResize(versOutCount, 3);
	versOutF.col(2).setZero();
	versOut = versOutF.array().cast<ScalarO>();

	return true;
}


// 重载2：2D点云三角剖分得到面网格——不带洞
template <typename DerivedVo, typename DerivedI, typename DerivedVi>
bool triangulateVers2Mesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedI>& trisOut, \
	const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const std::vector<Eigen::VectorXi>& bdryLoops, \
	const char* strSwitcher = "pY")
{
	return triangulateVers2Mesh(versOut, trisOut, versIn, bdryLoops, Eigen::MatrixXf{}, strSwitcher);
}


// 三角剖分提升网格质量：
template <typename DerivedVo, typename DerivedIt, typename DerivedVi, typename DerivedIe>
bool triangulateRefineMesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::PlainObjectBase<DerivedIt>& trisOut, const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const Eigen::PlainObjectBase<DerivedIt>& trisIn, const Eigen::PlainObjectBase<DerivedIe>& edges, \
	const char* strSwitcher = "pr")
{
	using ScalarO = typename DerivedVo::Scalar;

	assert((2 == versIn.cols() || 3 == versIn.cols()) && "assert!!! input vertices should be in 2 or 3 dimension space.");
	const int versCount = versIn.rows();
	versOut.resize(0, 0);
	trisOut.resize(0, 0);

	// 1. 生成边缘信息和洞的信息
	Eigen::MatrixXi edgesData;													// 连成闭合回路的一系列的边
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

	// 2. 整理已有三角片信息：
	const int trisCount = trisIn.rows();
	Eigen::MatrixXi trisData;
	if (trisCount > 0)
	{
		trisData = trisIn.transpose();
		trisData.array() += 1;
	}

	// 3. 生成三角剖分信息：
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

		inputTrig.numberofholes = 0;					// 洞的数目
		inputTrig.holelist = NULL;						// 洞数据数组的首地址，用一个二维点表示一个洞，只要该点在洞内即可。

		inputTrig.numberofregions = 0;
		inputTrig.regionlist = NULL;

		memset(&outputTrig, 0, sizeof(triangulateio));		// 初始化输出结构体。 
	}

	// 3. 执行三角剖分，拷贝结果：   
	char strCopy[512];
	strcpy(strCopy, strSwitcher);
	triangulate(strCopy, &inputTrig, &outputTrig, NULL);					// p——平面直线图，Y——不插点。  

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
	Eigen::MatrixXf versOutF(2, versOutCount);
	std::memcpy(versOutF.data(), reinterpret_cast<float*>(outputTrig.pointlist), sizeof(float) * 2 * versOutCount);
	versOutF.transposeInPlace();
	versOutF.conservativeResize(versOutCount, 3);
	versOutF.col(2).setZero();
	versOut = versOutF.array().cast<ScalarO>();

	return true;
}


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


int main(int argc, char** argv)
{


	return 0;
}