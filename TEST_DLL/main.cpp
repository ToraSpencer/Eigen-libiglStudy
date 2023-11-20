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
	static void debugWriteMesh(const char* name, const Eigen::MatrixBase<DerivedV>& vers, \
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


	// CGAL:: alpha shapes
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


	// CGAL:: remsh()
	void test7() 
	{
		Eigen::MatrixXf versIn;
		Eigen::MatrixXd versOut;
		Eigen::MatrixXi trisIn, trisOut;
		readOBJ(versIn, trisIn, "E:/��/mesh_side_ini_2.obj");
		debugWriteMesh("meshInput", versIn, trisIn);

		remeshSurfMesh(versOut, trisOut, versIn, trisIn, 2.0);
		debugWriteMesh("meshOut", versOut, trisOut);

		debugDisp("finished.");
	}

}


namespace TEST_JAW_PALATE 
{
	// step0����Ӧ���ڳ���һ��ʼ��ʱ��Ͷ�ȡ����������������ڿ���������������д���ڴ棻
	void test0()
	{
		Eigen::MatrixXf versJawUpper, versJawLower;
		Eigen::MatrixXi trisJawUpper, trisJawLower;
		readSTL(versJawUpper, trisJawUpper, "E:/��/��ʾ��/2/upper.stl");
		readSTL(versJawLower, trisJawLower, "E:/��/��ʾ��/2/lower.stl");

		// 1. ��ȡ���������0.1mm
		Eigen::MatrixXf versJawUpperExt, versJawLowerExt;
		Eigen::MatrixXi trisJawUpperExt, trisJawLowerExt;
		{ 
			SDF_RESULT sdfResult;
			genSDF(sdfResult, versJawUpper, trisJawUpper, 0.3, 3);
			marchingCubes(versJawUpperExt, trisJawUpperExt, sdfResult, 0.1); 
		}
		{
			SDF_RESULT sdfResult;
			genSDF(sdfResult, versJawLower, trisJawLower, 0.3, 3);
			marchingCubes(versJawLowerExt, trisJawLowerExt, sdfResult, 0.1);
		}

		debugWriteMesh("meshUpperEdited", versJawUpperExt, trisJawUpperExt);
		debugWriteMesh("meshLowerEdited", versJawLowerExt, trisJawLowerExt);

		debugDisp("finished.");
	}


	// step1�������������µ���߽��ߣ�
	/*
		���룺
				jawMeshUpper, jawMeshLower						���������
				{planeVerUpper, planeNormUpper}				�ϻ�׼ƽ��������һ�㡢ƽ�淨��
																								�ȼ��ڿռ����������������{va, vb, vc}��
																								�ȼ��ڿռ���ƽ�淽�̵�4��ϵ��{a0, b0, c0, d0}, ����[a0, b0, c0]��ƽ�淨��  
				{planeVerLower, planeNormLower}
								
		�����
				versUpperBdry, versLowerBdry						����֮ǰ�����µ���߽��ߣ�
				affineUpperHomo, affineLowerHomo				��������任����	
	*/
	void test1()
	{
		const float thresholdDis = 5.0;
		const int alphaValue = 10;
		Eigen::MatrixXf versJawUpper, versJawLower;
		Eigen::MatrixXi trisJawUpper, trisJawLower;
		Eigen::RowVector3f	planeCenterUpper{ 0, -15, -25 };
		Eigen::RowVector3f	planeNormUpper{ 0.0, 0.0, -1 };
		Eigen::RowVector3f	planeCenterLower{ 8, -20, -30 };
		Eigen::RowVector3f	planeNormLower{ -0.1, -0.5, 4 };

		readSTL(versJawUpper, trisJawUpper, "E:/��/��ʾ��/2/upper.stl");
		readSTL(versJawLower, trisJawLower, "E:/��/��ʾ��/2/lower.stl");
		planeNormUpper.normalize();
		planeNormLower.normalize();

		// 1. �������»�׼ƽ�棻
		{
			Eigen::MatrixXd versPlaneUpper, versPlaneLower;
			Eigen::MatrixXi trisPlaneUpper, trisPlaneLower;
			genRoundSurfMesh(versPlaneUpper, trisPlaneUpper, planeCenterUpper, planeNormUpper, 50);
			genRoundSurfMesh(versPlaneLower, trisPlaneLower, planeCenterLower, planeNormLower, 50);

			debugWriteMesh("meshPlaneUpper", versPlaneUpper, trisPlaneUpper);
			debugWriteMesh("meshPlaneLower", versPlaneLower, trisPlaneLower);
		}
		
		// 2. �����任������任���õ������ͶӰ�㣺
		Eigen::Matrix4f affineHomoUpper{ Eigen::Matrix4f::Identity() };
		Eigen::Matrix4f affineHomoInvUpper{ Eigen::Matrix4f::Identity() };
		Eigen::Matrix4f affineHomoLower{ Eigen::Matrix4f::Identity() };
		Eigen::Matrix4f affineHomoInvLower{ Eigen::Matrix4f::Identity() };
		Eigen::MatrixXf versProjUpper, versProjLower;
		{
			const int versCountUpper = versJawUpper.rows();
			const int versCountLower = versJawLower.rows();
			Eigen::Matrix3f rotationUpper, rotationLower;
			Eigen::MatrixXf versUpperHomo, versLowerHomo, versUpperSelectedHomo, versLowerSelectedHomo;
			getRotationMat(rotationUpper, Eigen::RowVector3f{ 0, 0, -1}, planeNormUpper);
			getRotationMat(rotationLower, Eigen::RowVector3f{ 0, 0, 1 }, planeNormLower);
			affineHomoUpper.topLeftCorner(3, 3) = rotationUpper;
			affineHomoUpper.topRightCorner(3, 1) = planeCenterUpper.transpose().array().cast<float>();
			affineHomoLower.topLeftCorner(3, 3) = rotationLower;
			affineHomoLower.topRightCorner(3, 1) = planeCenterLower.transpose().array().cast<float>();
			affineHomoInvUpper = affineHomoUpper.inverse();
			affineHomoInvLower = affineHomoLower.inverse();
			vers2HomoVers(versUpperHomo, versJawUpper);
			vers2HomoVers(versLowerHomo, versJawLower);
			versUpperHomo = (affineHomoInvUpper * versUpperHomo).eval();
			versLowerHomo = (affineHomoInvLower * versLowerHomo).eval(); 
			for (int i = 0; i < versCountUpper; ++i) 
				if (versUpperHomo(2, i) <= thresholdDis)
					matInsertCols(versUpperSelectedHomo, Eigen::Vector4f{ versUpperHomo.col(i) });
			for (int i = 0; i < versCountLower; ++i)
				if (versLowerHomo(2, i) >= -thresholdDis)
					matInsertCols(versLowerSelectedHomo, Eigen::Vector4f{ versLowerHomo.col(i) }); 

			versUpperSelectedHomo = (affineHomoUpper * versUpperSelectedHomo).eval();
			versLowerSelectedHomo = (affineHomoLower * versLowerSelectedHomo).eval();
			homoVers2Vers(versProjUpper, versUpperSelectedHomo);
			homoVers2Vers(versProjLower, versLowerSelectedHomo);
			versProjUpper.col(2).setZero();
			versProjLower.col(2).setZero();
			debugWriteVers("versProjUpper", versProjUpper);
			debugWriteVers("versProjLower", versProjLower);
		}

		// 3. alphashapes��ȡ�߽��ߣ�
		Eigen::MatrixXf bdryUpper, bdryLower;
		{
			// ����ͶӰ��ϲ���ʹ��ALPHA SHAPES��ȡ�߽��ߣ�
			Eigen::MatrixXf versProj, versBdry;
			matInsertRows(versProj, versProjUpper);
			matInsertRows(versProj, versProjLower);
			alphaShapes2D(versBdry, versProj, alphaValue);
			bdryUpper = versBdry;
			bdryLower = versBdry;

			debugWriteVers("versProj", versProj);
			debugWriteVers("versBdry", versBdry);
		}

		// 4. �߽����������ʱ�������Ļ�·��
		{
			// 4.1 �任�ر�������ϵ
			Eigen::MatrixXf bdryUpperHomo, bdryLowerHomo;
			vers2HomoVers(bdryUpperHomo, bdryUpper);
			homoVers2Vers(bdryUpper, affineHomoInvUpper * bdryUpperHomo);
			vers2HomoVers(bdryLowerHomo, bdryLower);
			homoVers2Vers(bdryLower, affineHomoInvLower * bdryLowerHomo);
			
			// 4.2 �����·
			sortLoop2D(bdryUpper);
			sortLoop2D(bdryLower);

			// 4. 3 �任��ԭ����ϵ��
			vers2HomoVers(bdryUpperHomo, bdryUpper);
			homoVers2Vers(bdryUpper, affineHomoUpper * bdryUpperHomo);
			vers2HomoVers(bdryLowerHomo, bdryLower);
			homoVers2Vers(bdryLower, affineHomoLower * bdryLowerHomo);

			// for debug
			debugWriteVers("bdryUpperSorted", bdryUpper);
			debugWriteVers("bdryLowerSorted", bdryLower);
		}



		debugDisp("finished.");
	}


	// step2�����������µ���߽��ߣ�
	/*
		���룺
			versUpperBdry, versLowerBdry						����֮ǰ�����µ���߽��ߣ�
			affineUpperHomo, affineLowerHomo				��������任����

		�����
			versUpperBdry, versLowerBdry						����֮������µ���߽��ߣ�

	*/
	void test2()
	{
		// 1. ���µ���߽��������任�������ռ�

		// 2. XOYƽ���ڵľ���

		// 3. Z�����ϵľ�����

		// 4. ������ı߽��߱任��ԭ�ռ䣻

	}


	// step3����������������Ȼ�����ӡ������ʹ��Բ�������ʷ�Ȼ��remesh����
	/*
		���룺
			versUpperBdry, versLowerBdry						����֮������µ���߽��ߣ�
			upperPlaneNorm, lowerPlaneNorm				���»�׼ƽ��ķ���

		�����
			meshPalate														��������						
	*/
	void test3()
	{
		// 0
		const int layersCount = 3;
		const float maxTriArea = 1.0;

		Eigen::MatrixXd versLoopUpper, versLoopLower;
		Eigen::RowVector3f upperPlaneNorm{ 0, 0, -1 };
		Eigen::RowVector3f lowerPlaneNorm{ -0.4, -0.5, 4 };
		upperPlaneNorm.normalize();
		lowerPlaneNorm.normalize();

		Eigen::RowVector3f topNorm = -upperPlaneNorm;
		Eigen::RowVector3f btmNorm = -lowerPlaneNorm;
		readOBJ(versLoopUpper, "E:\\��\\��ʾ��\\2\\curveFitUpper.obj");
		readOBJ(versLoopLower, "E:\\��\\��ʾ��\\2\\curveFitLower.obj");
		const int versCount1 = versLoopUpper.rows();
		const int versCount2 = versLoopUpper.rows();

		// 1. �������µ�������
		Eigen::MatrixXd versUpper, versLower;
		Eigen::MatrixXi trisUpper, trisLower;
		{
			circuit2mesh(versUpper, trisUpper, versLoopUpper, topNorm, maxTriArea);
			circuit2mesh(versLower, trisLower, versLoopLower, btmNorm, maxTriArea);
			debugWriteMesh("meshUpper", versUpper, trisUpper);
			debugWriteMesh("meshLower", versLower, trisLower);
		}

		// 2. ���ɲ�������
		Eigen::MatrixXd versSide, versSideExtra;
		Eigen::MatrixXi trisSide;
		{
			Eigen::MatrixXd versInner, versOuter, versTG, tmpVers;
			Eigen::MatrixXi trisTG, tmpTris;
			Eigen::VectorXi bdryOuter, bdryInner;

			getCircleVers(versInner, 10, versCount1);
			getCircleVers(versOuter, 20, versCount2);
			matInsertRows(tmpVers, versInner);
			matInsertRows(tmpVers, versOuter);
			debugWriteVers("versInner", versInner);
			debugWriteVers("versOuter", versOuter);
			debugWriteVers("versInput", tmpVers);

			Eigen::RowVector3d innerCenter = versInner.colwise().mean();
			bdryInner = Eigen::VectorXi::LinSpaced(versCount1, 0, versCount1 - 1);
			bdryOuter = Eigen::VectorXi::LinSpaced(versCount1, versCount1 + versCount2 - 1, versCount1);
			triangulateVers2Mesh(versTG, trisTG, tmpVers, std::vector<Eigen::VectorXi>{bdryInner, bdryOuter}, innerCenter);
			debugWriteMesh("meshTG", versTG, trisTG);

			matInsertRows(versSide, versLoopUpper);
			matInsertRows(versSide, versLoopLower);
			trisSide = trisTG;
			debugWriteMesh("meshSide0", versSide, trisSide);

			remeshSurfMesh(tmpVers, tmpTris, versSide, trisSide, 1.0);
			versSide = tmpVers;
			trisSide = tmpTris;
			const int versSideCount = versSide.rows();
			versSideExtra = versSide.bottomRows(versSideCount - versCount1 - versCount2);

			debugWriteMesh("meshSide", versSide, trisSide);
		}

		// 3. �������������
		Eigen::MatrixXd versPalate;
		Eigen::MatrixXi trisPalate;
		const int versUpperCount = versUpper.rows();
		const int versLowerCount = versLower.rows();
		const int versUpperExtraCount = versUpperCount - versCount1;
		const int versLowerExtraCount = versLowerCount - versCount2;
		{
			matInsertRows(versPalate, versUpper);
			matInsertRows(versPalate, versLower);
			matInsertRows(versPalate, versSideExtra);
			trisLower.array() += versUpperCount;

			// ������������Ƭ�У����ڵ���versCount1�Ķ���������Ҫʩ��ƫ������
			const int offset1 = versUpperExtraCount;
			const int offset2 = versUpperExtraCount + versLowerExtraCount;
			int* ptrInt = trisSide.data();
			for (int i = 0; i < trisSide.size(); ++i)
			{
				int& index = *ptrInt;
				if ((index >= versCount1) && (index < versCount1 + versCount2))
					index += offset1;
				else if (index >= versCount1 + versCount2)
					index += offset2;
				ptrInt++;
			}

			matInsertRows(trisPalate, trisUpper);
			matInsertRows(trisPalate, trisLower);
			matInsertRows(trisPalate, trisSide);

			debugWriteMesh("meshPalate", versPalate, trisPalate);
		} 


		// 4. �̶���Ȧ�߽�㣬��������й�˳��
		int versCount = versPalate.rows();
		int trisCount = trisPalate.rows();
		{
			Eigen::MatrixXd versTmp; 
			Eigen::MatrixXd versFixed;
			Eigen::MatrixXi trisTmp; 
			const int iterCount = 30;
			std::vector<int> fixedVerIdxes;

			// �̶���Ȧ�߽�㣺
			std::unordered_set<int> bdryVerIdxes;
			fixedVerIdxes.reserve(versCount);
			for (int i = 0; i < versCount; ++i)
			{
				if (i < versCount1)
					fixedVerIdxes.push_back(i);
				if ((i >= versCount1 + versUpperExtraCount) && (i < versCount1 + versCount2 + versUpperExtraCount))
					fixedVerIdxes.push_back(i); 
			}
			for (const auto& index : fixedVerIdxes)
				bdryVerIdxes.insert(index);
			subFromIdxVec(versFixed, versPalate, fixedVerIdxes);
			debugWriteVers("versFixed1", versFixed);

			smoothSurfMesh(versTmp, trisTmp, versPalate, trisPalate, fixedVerIdxes, 0.01);

			debugWriteMesh("meshLaplace1", versTmp, trisPalate);
			versPalate = versTmp;
		}


		// 5. �̶��߽��2��������Ķ��㣬��������й�˳��
		{
			Eigen::MatrixXd versTmp;
			Eigen::MatrixXd versFixed;
			Eigen::MatrixXi trisTmp;
			const float lambda = 0.1;
			const int iterCount = 30;
			std::vector<int> fixedVerIdxes;
			std::unordered_set<int> bdryVerIdxes, tmpSet;
			std::vector<int> triIdxesTaged;
			std::unordered_set<int> verIdxes1ring, verIdxes2ring;

			bdryVerIdxes.reserve(versCount);
			for (int i = 0; i < versCount; ++i)
			{
				if (i < versCount1)
					bdryVerIdxes.insert(i);
				if ((i >= versCount1 + versUpperExtraCount) && (i < versCount1 + versCount2 + versUpperExtraCount))
					bdryVerIdxes.insert(i);
			} 

			for (int i = 0; i < trisCount; ++i)
			{
				int vaIdx = trisPalate(i, 0);
				int vbIdx = trisPalate(i, 1);
				int vcIdx = trisPalate(i, 2);
				if (bdryVerIdxes.find(vaIdx) != bdryVerIdxes.end() || \
					bdryVerIdxes.find(vbIdx) != bdryVerIdxes.end() || \
					bdryVerIdxes.find(vcIdx) != bdryVerIdxes.end())
				{
					verIdxes1ring.insert(vaIdx);
					verIdxes1ring.insert(vbIdx);
					verIdxes1ring.insert(vcIdx);
				}
			}

			for (int i = 0; i < trisCount; ++i)
			{
				int vaIdx = trisPalate(i, 0);
				int vbIdx = trisPalate(i, 1);
				int vcIdx = trisPalate(i, 2);
				if (verIdxes1ring.find(vaIdx) != verIdxes1ring.end() || \
					verIdxes1ring.find(vbIdx) != verIdxes1ring.end() || \
					verIdxes1ring.find(vcIdx) != verIdxes1ring.end())
				{
					verIdxes2ring.insert(vaIdx);
					verIdxes2ring.insert(vbIdx);
					verIdxes2ring.insert(vcIdx);
				}
			}

			for (int i = 0; i < versCount; ++i)
				if (verIdxes2ring.end() == verIdxes2ring.find(i))
					fixedVerIdxes.push_back(i);

			for (const auto& index : fixedVerIdxes)
				tmpSet.insert(index);
			fixedVerIdxes.clear();
			fixedVerIdxes.insert(fixedVerIdxes.end(), tmpSet.begin(), tmpSet.end());

			versFixed.resize(0, 0);
			subFromIdxVec(versFixed, versPalate, fixedVerIdxes);
			debugWriteVers("versFixed2", versFixed);

			smoothSurfMesh(versTmp, trisTmp, versPalate, trisPalate, fixedVerIdxes, 0.01);

			debugWriteMesh("meshLaplace2", versTmp, trisPalate);

			versPalate = versTmp;
		}


		debugDisp("finished.");
	}


	// step4��������ӡ��
	/*
		���룺
			meshPalate														��������	
			meshUpperEdit, meshLowerEdit						���͡����������������
			
		�����
			meshJawPalate												���յ�������
	*/
	void test4() 
	{
		Eigen::MatrixXd versOut, versPalate;
		Eigen::MatrixXi trisOut, trisPalate;
		readOBJ(versPalate, trisPalate, "E:\\��\\��ʾ��\\2\\meshLaplace2.obj");
		 
		{
			Eigen::MatrixXd versUpper, versTmp;
			Eigen::MatrixXi trisUpper, trisTmp;
			readOBJ(versUpper, trisUpper, "E:\\��\\��ʾ��\\2\\meshUpperEdited.obj");
			meshCross(versTmp, trisTmp, versPalate, trisPalate, versUpper, trisUpper);
			versPalate = versTmp;
			trisPalate = trisTmp;

			debugWriteMesh("meshCross1", versTmp, trisTmp);
		}

		{
			Eigen::MatrixXd versLower, versTmp;
			Eigen::MatrixXi trisLower, trisTmp;
			readOBJ(versLower, trisLower, "E:\\��\\��ʾ��\\2\\meshLowerEdited.obj");

			meshCross(versTmp, trisTmp, versPalate, trisPalate, versLower, trisLower);
			versPalate = versTmp;
			trisPalate = trisTmp;
		}

		debugWriteMesh("meshJawPalate", versPalate, trisPalate);
		debugDisp("finished.");
	}


	// step5������ף�
	void test5() 
	{
		// 1. ���ɴ�׵�Բ����


	}
}


// ����triangle��������ʷ֣�
namespace TEST_TRIANGULATION 
{
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


	// �������������ʷֳ��������õ�����
	void test2() 
	{
		Eigen::MatrixXd versIn, versOut;
		Eigen::MatrixXi edgesIn, trisIn, trisOut;

		readOBJ(versIn, trisIn, "E:/����/roundSurfNoInterp.obj");		
		getEdges(edgesIn, trisIn);
		triangulateRefineMesh(versOut, trisOut, versIn, trisIn, edgesIn, "pra2.0");
		debugWriteMesh("meshOut1", versOut, trisOut);
 


		debugDisp("finished.");
	}


	// �������µ��滷·�����ʷֱַ����ɲ��桢�ϵ��桢�µ�������
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
		objReadVerticesMat(versLoopUpper, "E:/��/curveFitUpper2.obj");
		objReadVerticesMat(versLoopLower, "E:/��/curveFitLower2.obj");
	
		getCircleVers(versInner, 10, versCount1);
		getCircleVers(versOuter, 20, versCount2);
		matInsertRows(allVers, versInner);
		matInsertRows(allVers, versOuter);
		debugWriteVers("versInput", allVers);

		// 1. �����ʷ����ɲ��棺
		Eigen::MatrixXd versTG;
		Eigen::MatrixXi trisTG;
		Eigen::VectorXi bdryOuter, bdryInner;
		Eigen::RowVector3d innerCenter = versInner.colwise().mean();
		bdryInner = Eigen::VectorXi::LinSpaced(versCount1, 0, versCount1-1);
		bdryOuter = Eigen::VectorXi::LinSpaced(versCount1, versCount1 + versCount2 - 1, versCount1);
		triangulateVers2Mesh(versTG, trisTG, allVers, std::vector<Eigen::VectorXi>{bdryInner, bdryOuter}, innerCenter);
		debugWriteMesh("meshTG", versTG, trisTG);

		// 2. �����ʷ��������µ��棺
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

		// 3. ������������
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


	// ������֮��Ĳ����������������ʷֳ��������ߵ�����
	void test4() 
	{
		Eigen::MatrixXd versIn, versOut, versInner;
		Eigen::MatrixXi edgesIn, trisIn, trisOut;
		Eigen::VectorXi bdryOuter, bdryInner;
		Eigen::RowVector3d innerCenter;
		const int versCount1 = 300;
		const int versCount2 = 300;

		readOBJ(versIn, trisIn, "E:/��/mesh_side_ini_param_2.obj");
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
	// TEST_MYDLL::test7();

	TEST_JAW_PALATE::test1();

	debugDisp("main finished.");

	return 0;
}