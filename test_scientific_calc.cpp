#include "test_scientific_calc.h"


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

}
using namespace MY_DEBUG;


 
namespace TEST_SCIENTIFIC_CALC
{
	// ��������ʾ�µ�����任��
	void test1()
	{
		Eigen::MatrixXf vers, versTmpHomo;
		Eigen::MatrixXd versTmp;

		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/����/bunny.obj"); 
		debugWriteMesh("meshInput", vers, tris);

		// 1. ����ʹ�������������ʾ����һ������һ����Ԫ������(w*x, w*y, w*z, w)��wͨ��Ϊ1��Ϊ0ʱ��ʾ����Զ���ĵ㡣
		Eigen::MatrixXf versHomo;
		vers2HomoVers(versHomo, vers); 
		dispMatBlock(vers, 0, 3, 0, 2);
		dispMatBlock(versHomo, 0, 3, 0, 3); 

		versTmp.resize(4, 3);
		versTmp << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		dispMat(versTmp);
		vers2HomoVers(versTmpHomo, versTmp);
		dispMat(versTmpHomo);
		homoVers2Vers(versTmp, versTmpHomo);
		dispMat(versTmp);

		versTmp.resize(0, 0);
		versTmpHomo.resize(0, 0);
		versTmp.resize(4, 2);
		versTmp << 1, 2, 3, 4, 5, 6, 7, 8;
		dispMat(versTmp);
		vers2HomoVers(versTmpHomo, versTmp);
		dispMat(versTmpHomo);
		homoVers2Vers(versTmp, versTmpHomo);
		dispMat(versTmp);

		// 2. �ѿ���������ʩ�ӷ���任������ת�����������ƽ�ƣ�
		Eigen::Matrix3f scale = Eigen::Matrix3f::Identity();
		Eigen::Matrix3f rotation;
		Eigen::Vector3f moveVec(10, 0, 0);															// ����x������ƽ��20mm;
		Eigen::MatrixXf versDesc = vers.transpose();
		
		// scale(1, 1) = 2;																							// y��������������Ϊ2��
		rotation << cos(pi / 6), 0, sin(pi / 6), 0, 1, 0, -sin(pi / 6), 0, cos(pi / 6);		// ��z����ʱ����תpi/6�ȡ�
		versDesc = (rotation * scale * versDesc).eval();
		versDesc = (versDesc.colwise() + moveVec).eval();
		versDesc = versDesc.transpose().eval();
		objWriteMeshMat("E:/�ѿ��������·���任���bunny.obj", versDesc, tris);

		// 3. ���������ʩ�ӷ���任������ת�������� ���ƽ�ƣ�
		Eigen::Matrix4f scaleHomo = Eigen::Matrix4f::Identity();
		Eigen::Matrix4f rotationHomo = Eigen::Matrix4f::Identity();
		Eigen::Matrix4f moveMat = Eigen::Matrix4f::Identity();
		Eigen::MatrixXf finalVersHomo;
		scaleHomo.block<3, 3>(0, 0) = scale;
		rotationHomo.block<3, 3>(0, 0) = rotation;
		moveMat.topRightCorner(3, 1) = moveVec;
		 finalVersHomo = moveMat * rotationHomo * scaleHomo * versHomo;
		 homoVers2Vers(vers, finalVersHomo);
		 debugWriteMesh("��������·���任���bunny", vers, tris); 

		Eigen::MatrixXf transMat = moveMat * rotationHomo * scaleHomo;
		Eigen::MatrixXf originVersHomo = transMat.inverse() * finalVersHomo;
		homoVers2Vers(vers, originVersHomo);
		debugWriteMesh("��������»�ԭ��bunny", vers, tris); 

		std::cout << "finished" << std::endl;
	}


	// ���Ի��ɷ������ؾ����㷨���������ʽ�����ٳ˷����ô���������������ʣ�
	void test2()
	{
		// x^4 + 3x^3 + 2x^2 + x + c = c + x(1 + x(2*x + x(3 + x)));
		double x = 3.0;
		Eigen::Vector4d coeff(9, -1, 3, 1);
		std::cout << "���ɷ������㣺result == " << hornersPoly(coeff, x) << std::endl;

		Eigen::RowVector4d coeffTrans = coeff.transpose();
		Eigen::Vector4d X(1, x, std::powf(x, 2), std::powf(x, 3));
		std::cout << "��ͳ�����������ʽ��result == " << coeffTrans * X << std::endl;

		x = 5.2;
		Eigen::VectorXd coeff2(6);
		coeff2 << 2, 4, 5, 6, 8, 9;
		std::cout << "���ɷ������㣺result == " << hornersPoly(coeff2, x) << std::endl;

		Eigen::RowVectorXd coeffTrans2 = coeff2.transpose();
		Eigen::VectorXd X2(6);
		X2 << 1, x, std::powf(x, 2), std::powf(x, 3), std::powf(x, 4), std::powf(x, 5);
		std::cout << "��ͳ�����������ʽ��result == " << coeffTrans2 * X2 << std::endl;

		debugDisp("finished.");
	}


	// ���Լ��ֲ�ֵ������
	void test3()
	{

	}


	// ���Լ��ֺ�����Ϸ�����
	void test4()
	{


	}


	// ���ɷַ�����
	void test5()
	{
		Eigen::MatrixXf vers;
		Eigen::MatrixXi tris;
		//objReadMeshMat(vers, tris, "E:/����/jawMesh1.obj");
		//debugWriteMesh("meshInput", vers, tris); 
		objReadVerticesMat(vers, "E:/temp/versProjIntr2.obj");
		debugWriteVers("versInput", vers);

		std::vector<Eigen::RowVector3d> resultVecs;
		std::vector<double> resultVals;
		PCA3d(resultVecs, resultVals, vers);
		Eigen::RowVector3d vec_normal, vec_pca1, vec_pca2;
		vec_normal = resultVecs[0];
		vec_pca1 = resultVecs[1];
		vec_pca2 = resultVecs[2]; 

		// 
		Eigen::RowVector3f center = vers.colwise().mean();
		objWriteDirection("E:/normDir.obj", center, vec_normal);
		objWriteDirection("E:/pcaDir1.obj", center, vec_pca1);
		objWriteDirection("E:/pcaDir2.obj", center, vec_pca2);

		debugDisp("finished.");
	}


	// ���Լ�������ڿ˻�(Kronecker product);
	void test6()
	{
		Eigen::MatrixXd result, m1, m2;
		Eigen::Matrix2d m22;
		Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(9, 1, 9);
		m1 = Eigen::Map<Eigen::MatrixXd>(vec.data(), 3, 3);

		m2.resize(2, 2);
		m2 << 1, 1.1, 1, 1;
		m22 << 1, 1.1, 1, 1;
		
		result = kron(m1, m2);
		debugDisp("m1 == \n", m1);
		debugDisp("m2 == \n", m2);
		debugDisp("kron(m1, m2) == \n", result, "\n");

		result = kron(m1, m22); 
		debugDisp("kron(m1, m22) == \n", result, "\n");

		Eigen::Vector3f v1(1, 2, 3);
		Eigen::Vector4f v2{ Eigen::Vector4f::Ones() };
		result = kron(v1, v2); 
		debugDisp("kron(v1, v2) == \n", result, "\n");

		std::cout << "finished." << std::endl;
	}


	// ���������Բ��
	void test7()
	{
		// 1. ��ȡһȦ��������
		const unsigned teethCount = 14;
		std::vector<Eigen::MatrixXd> teethVers(teethCount);
		std::vector<Eigen::MatrixXi> teethTris(teethCount);
		Eigen::MatrixXd teethCenters(teethCount, 3);
		for (unsigned i = 0; i<teethCount; ++i) 
		{
			char fileName[256];
			sprintf_s(fileName, "E:/����/teeth/compMesh%d.obj", i);
			objReadMeshMat(teethVers[i], teethTris[i], fileName);
			teethCenters.row(i) = teethVers[i].colwise().mean();
		}

		// 2. �������ĵ�ͶӰ��XOY���ϣ���Ϊ�����㣬�����Բ��
		teethCenters.col(2).array() = 0;
		Eigen::VectorXd coeffs = fittingStandardEllipse(teethCenters);
		double a = coeffs(0);
		double c = coeffs(1);
		double d = coeffs(2);
		double e = coeffs(3);
		double f = coeffs(4);
		double x0 = -d / (2 * a);		// ��ԲԲ�ģ��볤�����᣺
		double y0 = -e / (2 * c);
		double p = d * d / (4 * a) + e * e / (4 * c) - f;
		double a0 = std::sqrt(p / a);
		double b0 = std::sqrt(p / c);

		// 3. ������Բ�㼯��
		unsigned versCount = 361;
		Eigen::VectorXd theta = Eigen::VectorXd::LinSpaced(versCount, 0, 2 * pi);
		Eigen::MatrixXd elliVers(Eigen::MatrixXd::Zero(versCount, 3));
		elliVers.col(0).array() = a0 * Eigen::cos(theta.array()) + x0;
		elliVers.col(1).array() = b0 * Eigen::sin(theta.array()) + y0;

		debugWriteVers("elliVers", elliVers);

		debugDisp("finished.");
	}
	 
}

 

 
namespace TEST_PMP
{

	// ��������������������
	void test2()
	{
		Eigen::MatrixXf vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/����/jawMesh.obj");

		double volume = meshVolume(vers, tris);
		std::cout << "volume == " << volume << std::endl;

		std::cout << "finished." << std::endl;
	}


	// ���������������������:
	void test3() 
	{
		// �����������еķ����ΰ�ߣ�
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		Eigen::MatrixXi nmnEdges;

		objReadMeshMat(vers, tris, "E:/����/jawMeshDense4.obj");
		nonManifoldUEs(nmnEdges, tris);
		objWriteEdgesMat("E:/nmnEdges.obj", nmnEdges, vers);

		std::cout << "trisCount == " << tris.rows() << std::endl;
		std::cout << "���������������" << nmnEdges.rows() << std::endl;

		unsigned trisSum = 0;
		for (unsigned i = 0; i < 40; ++i) 
		{
			char str[512];
			sprintf_s(str, 512, "E:/���񾫼�/splittedData4/splitedMesh%d.obj", i);
			Eigen::MatrixXd vers;
			Eigen::MatrixXi tris;
			Eigen::MatrixXi nmnEdges;
			objReadMeshMat(vers, tris, str);
			nonManifoldUEs(nmnEdges, tris);

			trisSum += tris.rows();
			if (nmnEdges.rows() > 0) 
				std::cout << "splitedMesh" << i << "�з����������" << nmnEdges.rows() << "����" << std::endl;
		}
		std::cout << "trisSum == " << trisSum << std::endl;
		std::cout << "finished." << std::endl;


		objReadMeshMat(vers, tris, "E:/����/meshNmnEdges.obj");
		nonManifoldUEs(nmnEdges, tris);
		objWriteEdgesMat("E:/nmnEdges.obj", nmnEdges, vers);
	}

}
 
 





