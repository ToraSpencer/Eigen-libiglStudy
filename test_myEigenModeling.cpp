#include "test_myEigenModeling.h"
#include "representations.h"

// SDFGEN_DLL��̬�⣺
#include "SDFGEN_DLL.h"
#pragma comment(lib,"SDFGEN_DLL.lib")	


////////////////////////////////////////////////////////////////////////////////////////////// DEBUG �ӿ�
namespace MY_DEBUG
{
	static std::string g_debugPath = "E:/";

	// lambda������ӡstd::cout֧�ֵ����ͱ�����
	template <typename T>
	static auto disp = [](const T& arg)
	{
		std::cout << arg << ", ";
	};

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


	static void debugWriteVers(const char* fileName, const std::vector<verF>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), fileName);
		writeOBJ(path, vers);
	}


	static void debugWriteVers(const char* fileName, const std::vector<verD>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), fileName);
		writeOBJ(path, vers);
	}


	static void debugWriteMesh(const char* fileName, const triMeshD& triMesh)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), fileName);
		writeOBJ(path, triMesh);
	}


	static void debugWriteMesh(const char* fileName, const triMeshF& triMesh)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), fileName);
		writeOBJ(path, triMesh);
	}
}
using namespace MY_DEBUG;


///////////////////////////////////////////////////////////////////////////////////////////////// ���Ժ���

// ����myEigenModeling�еĽӿ�
namespace TEST_MYEIGEN_MODELING
{
#ifdef USE_TRIANGLE_H

	// ������������
	void test1()
	{
		Eigen::MatrixXf axis(15, 3);
		Eigen::MatrixXd surfVers, circleVers;
		Eigen::MatrixXi surfTris;
		Eigen::MatrixXf cylinderVers;
		Eigen::MatrixXi cylinderTris;

		// 1. �����������ߣ�
		axis.setZero();
		double deltaTheta = pi / 10;
		for (unsigned i = 0; i < 15; ++i)
		{
			double theta = deltaTheta * i;
			axis(i, 0) = 50 * cos(theta);
			axis(i, 1) = 50 * sin(theta);
		}
		debugWriteVers("axis", axis);

		// 2. �����·���������ʷ֣��õ�Բ�ε�������
		objReadVerticesMat(circleVers, "E:/����/circleVers.obj");
		circuit2mesh(surfVers, surfTris, circleVers);
		debugWriteMesh("surfMesh", surfVers, surfTris);

		// 3. ����Բ���壺
		genCylinder(cylinderVers, cylinderTris, axis, 10);				// ����Բ��
		debugWriteMesh("cylinder", cylinderVers, cylinderTris);

		// 4. ���ɷ����壻
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

		// 5. ��ȡ���µ���2D�߽绷·���������壺
		Eigen::MatrixXd topLoop, btmLoop;
		cylinderVers.resize(0, 0);
		cylinderTris.resize(0, 0);
		axis.resize(0, 0);
		objReadVerticesMat(topLoop, "E:/��/bdryUpper_final.obj");
		objReadVerticesMat(btmLoop, "E:/��/bdryLower_final.obj");
		interpolateToLine(axis, Eigen::RowVector3f{ 0, 0, 0 }, Eigen::RowVector3f{ 0, 0, 5 }, 1.0);
		genCylinder(cylinderVers, cylinderTris, axis, topLoop, btmLoop, true);
		debugWriteMesh("pillar2", cylinderVers, cylinderTris);

		// 5. ��ȡ���µ���3D�߽绷·���������壺
		topLoop.resize(0, 0);
		btmLoop.resize(0, 0);
		cylinderVers.resize(0, 0);
		cylinderTris.resize(0, 0);
		axis.resize(0, 0);
		objReadVerticesMat(topLoop, "G:\\gitRepositories\\matlabCode\\��\\data/curveFitUpper.obj");
		objReadVerticesMat(btmLoop, "G:\\gitRepositories\\matlabCode\\��\\data/curveFitLower.obj");
		genCylinder(cylinderVers, cylinderTris, topLoop, btmLoop, Eigen::RowVector3f(0, 0, 1), Eigen::RowVector3f(0, 0, -1), 3, 1, true);
		debugWriteMesh("pillar3", cylinderVers, cylinderTris);

		std::cout << "finished." << std::endl;
	}


	// ���Ա߽��·������������
	void test2()
	{
		Eigen::MatrixXf vers;
		Eigen::MatrixXd versLoop, versLoop2D;
		Eigen::MatrixXi tris;

		objReadVerticesMat(versLoop, "G:\\gitRepositories\\matlabCode\\��\\data\\curveFitUpper2.obj");

		Eigen::RowVector3d center = versLoop.colwise().mean();
		versLoop2D = versLoop.rowwise() - center;
		versLoop2D = versLoop2D.leftCols(2).eval();

		circuit2mesh(vers, tris, versLoop, Eigen::RowVector3d{ 0, 0, -1 }, 1.0f);

		debugWriteMesh("meshOut", vers, tris);

		debugDisp("finished.");
	}


	// �������ɹ���ͼ�Σ�
	void test3()
	{
		Eigen::MatrixXd versRound;
		Eigen::MatrixXi trisRound;
		genRoundSurfMesh(versRound, trisRound, Eigen::RowVector3f{ 1, 2, 3 }, Eigen::RowVector3f{ 0, 0, -1 });
		debugWriteMesh("meshRound", versRound, trisRound);

		debugDisp("finished.");
	}

#endif

}


namespace TEST_SDF
{
	// ��ȡ��������SDF���ݣ�Ȼ��ִ��MC
	void test1() 
	{
		triMeshF mesh;
		readOBJ(mesh, "E:/����/tooth.obj");
		debugWriteMesh("meshInput", mesh);

		// 1. 
		SDF_RESULT sdfData;
		genSDF3D(sdfData, mesh, 0.3, 10);

		// 2. 
		Eigen::MatrixXf grids;
		Eigen::RowVector3f origin{ sdfData.origin.x, sdfData.origin.y, sdfData.origin.z };
		genGrids(grids, origin, sdfData.step, sdfData.stepsCount);

		// 2. 
		const float isoValue = 0.5;
		Eigen::MatrixXd versOut;
		Eigen::MatrixXi trisOut;
		Eigen::VectorXf sdfVec;
		vec2EigenVec(sdfVec, sdfData.SDFvalues);
		marchingCubes(versOut, trisOut, sdfVec, grids, sdfData.stepsCount[0], sdfData.stepsCount[1],\
			sdfData.stepsCount[2], isoValue, true);
		debugWriteMesh("meshOut", versOut, trisOut);

		debugDisp("finished.");
	}

}
