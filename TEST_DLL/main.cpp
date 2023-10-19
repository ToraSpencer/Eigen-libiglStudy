#include "smarteeDLL.h"
#include "my_cgal_dll.h"

// ��̬��ĵ���lib:
#pragma comment(lib, "smarteeDLL.lib")
#pragma comment(lib, "MY_CGAL_DLL.lib")

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
}



int main(int argc, char** argv)
{
	TEST_MYDLL::test5();

	debugDisp("main() finished.");

	return 0;
}