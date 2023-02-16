#include "scientific_calc.h"


static igl::opengl::glfw::Viewer viewer;				// libigl�еĻ���glfw����ʾ���ڣ�

namespace SCIENTIFICCALC
{
	// ��������ʾ�µ�����任��
	void test1()
	{
		MatrixXf vers;
		MatrixXi tris;
		objReadMeshMat(vers, tris, "./data/bunny.obj");

		// 1. ����ʹ�������������ʾ����һ������һ����Ԫ������(w*x, w*y, w*z, w)��wͨ��Ϊ1��Ϊ0ʱ��ʾ����Զ���ĵ㡣
		MatrixXf versHomo;
		vers2homoVers(versHomo, vers);
		dispMatBlock(vers, 0, 3, 0, 2);
		dispMatBlock(versHomo, 0, 3, 0, 3);

		// 2. �ѿ���������ʩ�ӷ���任������ת�����������ƽ�ƣ�
		Matrix3f scale = Matrix3f::Identity();
		scale(1, 1) = 2;													// y��������������Ϊ2��
		Matrix3f rotation;
		rotation << cos(pi / 6), 0, sin(pi / 6), 0, 1, 0, -sin(pi / 6), 0, cos(pi / 6);	// ��z����ʱ����תpi/6�ȡ�
		Vector3f moveVec(10, 0, 0);			// ����x������ƽ��20mm;
		MatrixXf versDesc = vers.transpose();
		versDesc = (rotation * scale * versDesc).eval();
		versDesc = (versDesc.colwise() + moveVec).eval();
		objWriteMeshMat<float>("E:/�ѿ��������·���任���bunny.obj", versDesc.transpose(), tris);

		// 3. ���������ʩ�ӷ���任������ת�������� ���ƽ�ƣ�
		Matrix4f scaleHomo = Matrix4f::Identity();
		scaleHomo.block<3, 3>(0, 0) = scale;
		Matrix4f rotationHomo = Matrix4f::Identity();
		rotationHomo.block<3, 3>(0, 0) = rotation;
		Matrix4f moveMat = Matrix4f::Identity();
		moveMat.topRightCorner(3, 1) = moveVec;

		MatrixXf finalVersHomo = moveMat * rotationHomo * scaleHomo * versHomo;
		objWriteMeshMat("E:/��������·���任���bunny.obj", homoVers2vers(finalVersHomo), tris);
 
		MatrixXf transMat = moveMat * rotationHomo * scaleHomo;
		MatrixXf originVersHomo = transMat.inverse() * finalVersHomo;
		objWriteMeshMat("E:/��������»�ԭ��bunny.obj", homoVers2vers(originVersHomo), tris);


		std::cout << "finished" << std::endl;
	}


	// ���Ի��ɷ������ؾ����㷨���������ʽ��
	void test2()
	{
		Eigen::Vector4f coeff(9, -1, 3, 1);
		Eigen::RowVector4f coeffTrans = coeff.transpose();
		float x = 3.0;
		Eigen::Vector4f X(1, x, std::powf(x, 2), std::powf(x, 3));

		std::cout << "result == " << hornersPoly(coeff, x) << std::endl;
		std::cout << "real result == " << coeffTrans * X << std::endl;
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


	}


	// ���Լ�������ڿ˻�(Kronecker product);
	void test6() 
	{
		Eigen::MatrixXd result, m1, m2;
		Eigen::VectorXd vec = Eigen::VectorXd::LinSpaced(100, 1, 100);

		m1 = Eigen::Map<MatrixXd>(vec.data(), 10, 10);
		m2.resize(2, 2);
		m2 << 1, 1.1, 1, 1;
		Eigen::Matrix2d m22;
		m22 << 1, 1.1, 1, 1;

		kron(result, m1, m2);
		dispMat(m1);
		dispMat(result);

		kron(result, m1, m22);
		dispMat(result);

		Eigen::Vector3f v1(1,2,3);
		Eigen::Vector4f v2{Eigen::Vector4f::Ones()};
		kron(result, v1, v2);
		dispMat(result);

		std::cout << "finished." << std::endl;
	}
}


namespace TEST_PMP
{
	// ��������Ľӿ�
	void test1() 
	{
		Eigen::MatrixXd vers, triNorms, barycenters, edgeArrows, vas, vbs;
		Eigen::MatrixXi tris, edges;
		objReadMeshMat(vers, tris, "E:/����/tooth.obj");

		getEdges(edges, vers, tris);
		trianglesBarycenter(barycenters, vers, tris);
		trianglesNorm(triNorms, vers, tris);
		objWriteVerticesMat("E:/barycenters.obj", barycenters);

		int edgesCount = edges.rows();
		Eigen::VectorXi vaIdxes = edges.col(0);
		Eigen::VectorXi vbIdxes = edges.col(1);
		subFromIdxVec(vas, vers, vaIdxes);
		subFromIdxVec(vbs, vers, vbIdxes);
		edgeArrows = vbs - vas;
		Eigen::VectorXd edgesLen = edgeArrows.rowwise().norm();

		viewer.data().set_mesh(vers, tris);
		viewer.data().show_lines = false;                 // ����������
		viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);    // �趨����������ת

		// ����Ƭ����ָʾ���ñ����ݵ���ʽ��Ⱦ������
		const RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);    // RGBɫ��������
		double aveLen = edgesLen.mean();							// ���б߳���ƽ��ֵ��
		viewer.data().add_edges(barycenters - aveLen * triNorms, barycenters + aveLen * triNorms, red);         // ������ʷ����ú�ɫָʾ�߱�ʶ

		viewer.launch();

		std::cout << "finished." << std::endl;
	}


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
		nonManifoldEdges(tris, nmnEdges);
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
			nonManifoldEdges(tris, nmnEdges);

			trisSum += tris.rows();
			if (nmnEdges.rows() > 0) 
				std::cout << "splitedMesh" << i << "�з����������" << nmnEdges.rows() << "����" << std::endl;
		}
		std::cout << "trisSum == " << trisSum << std::endl;
		std::cout << "finished." << std::endl;


		objReadMeshMat(vers, tris, "E:/����/meshNmnEdges.obj");
		nonManifoldEdges(tris, nmnEdges);
		objWriteEdgesMat("E:/nmnEdges.obj", nmnEdges, vers);
	}

}


// libigl�е���ѧ����
namespace IGL_MATH 
{
	// ��������ֵ������������
	void test0() 
	{

	}

	// �����ڿ˻���
	void test1() 
	{
		Eigen::MatrixXi m1(2, 3);
		Eigen::Matrix2d m2;
		Eigen::MatrixXd result;
		m1 << 1, 2, 3, 4, 5, 6;
		m2 << 100, 10, 1, 0.1;
		dispMat(m1);
		dispMat(m2);
		dispMat(kron(m1, m2));

		std::cout << "finished." << std::endl;
	}
}





