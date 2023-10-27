#include "smarteeDLL.h"
#include "my_cgal_dll.h"

// ��̬��ĵ���lib:
#pragma comment(lib, "smarteeDLL.lib")
#pragma comment(lib, "MY_CGAL_DLL.lib")

// myEigen��̬�⣺
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

// �ѽ���vcpkg;


// ������ע������Ŀ���ԡ����ԡ������и�����PATH��ֵΪDLL�ļ�����Ŀ¼��
/*
	�������ô���Ļ�������ֻ������Ŀ�����ļ����Լ�����ļ���������DLL�ļ���
	�����Ǳ�������Ե�dll�ļ���������������dll�ļ���Ҫ���뵽���Ա���������Ŀ¼�У�
*/


/////////////////////////////////////////////////////////////////////////////////////////////////////////// DEBUG�ӿڣ�
namespace MY_DEBUG
{
	static const std::string g_debugPath = "E:/";

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


/////////////////////////////////////////////////////////////////////////////////////////////////////////// ����smarteeDLL.dll�еĵ�������
namespace TEST_MYDLL
{
	// ���Է�װ��IO���ܣ�
	void test0()
	{
		Eigen::MatrixXf versF1, versF2;
		Eigen::MatrixXd versD1, versD2;
		Eigen::MatrixXi tris21, tris22;

		readOBJ(versF1, "E:/����/circleVers.obj");
		readOBJ(versD1, "E:/����/circleVers.obj");
		readOBJ(versF2, tris21, "E:/����/tooth.obj");
		readOBJ(versD2, tris22, "E:/����/tooth.obj");

		debugWriteVers("versOut1", versD1);
		debugWriteMesh("meshOut1", versF2, tris21);
		debugWriteMesh("meshOut2", versD2, tris22);

		debugDisp("finished.");
	}


	// �����ϰ�B������ֵ��
	void test1()
	{
		Eigen::MatrixXf versIn, versOut;
		Eigen::RowVector3f planeVer{ 0, 0, 0 };
		Eigen::RowVector3f planeNorm{ 0, 0, 1 };
		readOBJ(versIn, "E:/����/circleVers.obj");
		BsplineInterp(versOut, versIn, planeVer, planeNorm);
		writeOBJ("E:/versBsplineInterp.obj", versOut);

		debugDisp("finished.");
	}


	// �����°�����B������ֵ������ߣ�
	void test2()
	{
		Eigen::MatrixXf versIn, versOut;
		readOBJ(versIn, "E:/��/sampleVersUpperBdry.obj");
		interpltnCubicBSpline(versOut, versIn, 300, true);

		debugWriteVers("versOut", versOut);
		debugDisp("finished.");
	}


	// ���������󽻹��ܣ�
	void test3()
	{
		Eigen::MatrixXd srcVers, meshVers;
		Eigen::MatrixXi srcTris, meshTris;
		Eigen::RowVector3d dir{ 0, 0, -1 };
		meshRayOut result;
		readOBJ(srcVers, srcTris, "E:/����/fatTeeth.obj");
		readOBJ(meshVers, meshTris, "E:/����/jawMesh.obj");
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

		debugWriteVers("ͶӰ��", positiveProVers);

		debugDisp("finished.");
	}


	// alpha shapes
	void test4() 
	{
		Eigen::MatrixXf versIn, versOut;
		readOBJ(versIn, "E:/��/2DbdryUpper.obj");
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

		readOBJ(vers, tris, "E:/����/tooth.obj");
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

		readSTL(vers, tris, "E:/��/zhangweiliang-S.stl");
		genSDF(sdfResult, vers, tris, 0.3, 3);
		marchingCubes(versOut, trisOut, sdfResult, 0.1);
		debugWriteMesh("zhangweiliang-S_0.1", versOut, trisOut);

		vers.resize(0, 0);
		versOut.resize(0, 0);
		tris.resize(0, 0);
		trisOut.resize(0, 0);
		readSTL(vers, tris, "E:/��/zhangweiliang-X.stl");
		genSDF(sdfResult, vers, tris, 0.3, 3);
		marchingCubes(versOut, trisOut, sdfResult, 0.1);
		debugWriteMesh("zhangweiliang-X_0.1", versOut, trisOut);

		debugDisp("finished.");		 
	}


	// ���񲼶�����
	void test6()
	{
		Eigen::MatrixXd versOut, vers1, vers2;
		Eigen::MatrixXi trisOut, tris1, tris2;
		readOBJ(vers1, tris1, "E:/����/tooth.obj");
		readOBJ(vers2, tris2, "E:/����/cylinder1.obj");
		meshCross(versOut, trisOut, vers1, tris1, vers2, tris2);

		debugWriteMesh("meshOut", versOut, trisOut);

		debugDisp("finished.");
	}


	void test66() 
	{
		Eigen::MatrixXd versOut, vers1, vers2;
		Eigen::MatrixXi trisOut, tris1, tris2;

		readOBJ(vers1, tris1, "E:/��/palate.obj");
		readOBJ(vers2, tris2, "E:/��/zhangweiliang-S_0.1.obj");
		meshCross(versOut, trisOut, vers1, tris1, vers2, tris2);
		debugWriteMesh("meshCross1", versOut, trisOut);

		vers1 = versOut;
		tris1 = trisOut;
		versOut.resize(0, 0);
		trisOut.resize(0, 0);
		vers2.resize(0, 0);
		tris2.resize(0, 0);
		readOBJ(vers2, tris2, "E:/��/zhangweiliang-X_0.1.obj");
		meshCross(versOut, trisOut, vers1, tris1, vers2, tris2);
		debugWriteMesh("meshCross2", versOut, trisOut);

		debugDisp("finished.");
	}
}


namespace TEST_JAW_PALATE 
{

	void step3() 
	{
		// 1. ��ȡ���µ���������ߣ�������������
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
		readSTL(versJawUpper, trisJawUpper, "E:/��/��ʾ��/2/upper.stl");
		readSTL(versJawLower, trisJawLower, "E:/��/��ʾ��/2/lower.stl");
		readOBJ(versLoopUpper, "G:\\gitRepositories\\matlabCode\\��\\data\\curveFitUpper2.obj");
		readOBJ(versLoopLower, "G:\\gitRepositories\\matlabCode\\��\\data\\curveFitLower2.obj");
		genCylinder(versPalate, trisPalate, versLoopUpper, versLoopLower,\
				topNorm, btmNorm, layersCount, maxTriArea, true);
		debugWriteMesh("meshPalate", versPalate, trisPalate);


		// for debug�������ɲ���������ʷֵ��ϰ汾��
		{
			Eigen::MatrixXd versTmp;
			Eigen::MatrixXi trisTmp;
			circuit2mesh(versTmp, trisTmp, versLoopUpper);
			debugWriteMesh("meshUpperOld", versTmp, trisTmp);
		}

		// for debug����������µ�������
		{
			Eigen::MatrixXd upperVers, lowerVers;
			Eigen::MatrixXi upperTris, lowerTris;
			circuit2mesh(upperVers, upperTris, versLoopUpper, topNorm, maxTriArea);
			circuit2mesh(lowerVers, lowerTris, versLoopLower, btmNorm, maxTriArea);
			debugWriteMesh("meshUpper", upperVers, upperTris);
			debugWriteMesh("meshLower", lowerVers, lowerTris);
		}


		// 2. ��ȡ���������0.1mm, meshcross
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
	const Eigen::PlainObjectBase<DerivedVi>& versIn,\
	const std::vector<Eigen::VectorXi>& bdryLoops, \
	const Eigen::PlainObjectBase<DerivedVh>& holeCenters, \
	const char* strSwitcher = "pY")
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
	const char* strSwitcher = "pY")
{
	return triangulateVers2Mesh(versOut, trisOut, versIn, bdryLoops, Eigen::MatrixXf{}, strSwitcher);
}


// �����ʷ���������������
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
	std::memcpy(versOutF.data(), reinterpret_cast<float*>(outputTrig.pointlist), sizeof(float) * 2 * versOutCount);
	versOutF.transposeInPlace();
	versOutF.conservativeResize(versOutCount, 3);
	versOutF.col(2).setZero();
	versOut = versOutF.array().cast<ScalarO>();

	return true;
}


// ���Ե��������ʷֳ�����
void test1() 
{
	Eigen::MatrixXd versIn;
	Eigen::MatrixXf versOut;
	Eigen::MatrixXi trisIn, trisOut;

	// 1. �����ʷֳɲ�����������
	readOBJ(versIn, "E:/����/circleVers.obj");
	const int versCount = versIn.rows();
	Eigen::VectorXi loopVec = Eigen::VectorXi::LinSpaced(versCount, 0, versCount - 1);
	triangulateVers2Mesh(versOut, trisOut, versIn, std::vector<Eigen::VectorXi>{loopVec});
	debugWriteVers("versInput1", versIn);
	debugWriteMesh("meshOut1", versOut, trisOut);

	// 2. �����ʷֳɴ���������

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

	//		2.2. ���ɱ�Ե��Ϣ��
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