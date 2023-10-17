#include "smarteeDLL.h"

// ��̬��ĵ���lib:
#pragma comment(lib, "smarteeDLL.lib")


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
		Eigen::MatrixXf srcVers, meshVers;
		Eigen::MatrixXi srcTris, meshTris;
		Eigen::RowVector3f dir{ 0, 0, -1 };
		meshRayOut result;
		readOBJ(srcVers, srcTris, "E:/����/fatTeeth.obj");
		readOBJ(meshVers, meshTris, "E:/����/jawMesh.obj");
		debugWriteMesh("meshInput1", srcVers, srcTris);
		debugWriteMesh("meshInput2", meshVers, meshTris);

		meshRayIntersect(result, srcVers, dir, meshVers, meshTris);

		std::vector<Eigen::RowVector3f> tmpVec;
		tmpVec.reserve(srcVers.rows());
		for (int i = 0; i < srcVers.rows(); ++i)
		{
			if (result.rayLen[i].size() > 0)
			{
				float projDepth = result.rayLen[i][0];
				Eigen::RowVector3f tmpVer = srcVers.row(i) + projDepth * dir;
				tmpVec.push_back(tmpVer);
			}
		}
		tmpVec.shrink_to_fit();

		Eigen::MatrixXf positiveProVers(tmpVec.size(), 3);
		for (int i = 0; i < positiveProVers.rows(); ++i)
			positiveProVers.row(i) = tmpVec[i];

		debugWriteVers("ͶӰ��", positiveProVers);

		debugDisp("finished.");
	}
}



int main(int argc, char** argv)
{
	TEST_MYDLL::test3();

	debugDisp("main() finished.");

	return 0;
}