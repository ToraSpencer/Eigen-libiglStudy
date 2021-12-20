#include "scientific_calc.h"


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
		dispMatBlock<float>(vers, 0, 3, 0, 2);
		dispMatBlock<float>(versHomo, 0, 3, 0, 3);

		// 2. �ѿ���������ʩ�ӷ���任������ת�����������ƽ�ƣ�
		Matrix3f scale = Matrix3f::Identity();
		scale(1, 1) = 2;													// y��������������Ϊ2��
		Matrix3f rotation;
		rotation << cos(pi / 6), 0, sin(pi / 6), 0, 1, 0, -sin(pi / 6), 0, cos(pi / 6);	// ��z����ʱ����תpi/6�ȡ�
		Vector3f moveVec(10, 0, 0);			// ����x������ƽ��20mm;
		MatrixXf versDesc = vers.transpose();
		versDesc = (rotation * scale * versDesc).eval();
		versDesc = (versDesc.colwise() + moveVec).eval();
		objWriteMeshMat("E:/�ѿ��������·���任���bunny.obj", versDesc.transpose(), tris);

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
}

