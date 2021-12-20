#include "scientific_calc.h"


namespace SCIENTIFICCALC
{
	// 齐次坐标表示下的坐标变换；
	void test1()
	{
		MatrixXf vers;
		MatrixXi tris;
		objReadMeshMat(vers, tris, "./data/bunny.obj");

		// 1. 顶点使用齐次坐标来表示——一个点是一个四元列向量(w*x, w*y, w*z, w)，w通常为1，为0时表示无穷远处的点。
		MatrixXf versHomo;
		vers2homoVers(versHomo, vers);
		dispMatBlock<float>(vers, 0, 3, 0, 2);
		dispMatBlock<float>(versHomo, 0, 3, 0, 3);

		// 2. 笛卡尔坐标下施加仿射变换——旋转、放缩，最后平移：
		Matrix3f scale = Matrix3f::Identity();
		scale(1, 1) = 2;													// y方向上缩放因子为2；
		Matrix3f rotation;
		rotation << cos(pi / 6), 0, sin(pi / 6), 0, 1, 0, -sin(pi / 6), 0, cos(pi / 6);	// 绕z轴逆时针旋转pi/6度。
		Vector3f moveVec(10, 0, 0);			// 朝着x正方向平移20mm;
		MatrixXf versDesc = vers.transpose();
		versDesc = (rotation * scale * versDesc).eval();
		versDesc = (versDesc.colwise() + moveVec).eval();
		objWriteMeshMat("E:/笛卡尔坐标下仿射变换后的bunny.obj", versDesc.transpose(), tris);

		// 3. 齐次坐标下施加仿射变换——旋转、放缩， 最后平移：
		Matrix4f scaleHomo = Matrix4f::Identity();
		scaleHomo.block<3, 3>(0, 0) = scale;
		Matrix4f rotationHomo = Matrix4f::Identity();
		rotationHomo.block<3, 3>(0, 0) = rotation;
		Matrix4f moveMat = Matrix4f::Identity();
		moveMat.topRightCorner(3, 1) = moveVec;

		MatrixXf finalVersHomo = moveMat * rotationHomo * scaleHomo * versHomo;
		objWriteMeshMat("E:/齐次坐标下仿射变换后的bunny.obj", homoVers2vers(finalVersHomo), tris);
 
		MatrixXf transMat = moveMat * rotationHomo * scaleHomo;
		MatrixXf originVersHomo = transMat.inverse() * finalVersHomo;
		objWriteMeshMat("E:/齐次坐标下还原的bunny.obj", homoVers2vers(originVersHomo), tris);


		std::cout << "finished" << std::endl;
	}


	// 测试霍纳方法（秦九昭算法）计算多项式：
	void test2()
	{
		Eigen::Vector4f coeff(9, -1, 3, 1);
		Eigen::RowVector4f coeffTrans = coeff.transpose();
		float x = 3.0;
		Eigen::Vector4f X(1, x, std::powf(x, 2), std::powf(x, 3));

		std::cout << "result == " << hornersPoly(coeff, x) << std::endl;
		std::cout << "real result == " << coeffTrans * X << std::endl;
	}


	// 测试几种插值方法：
	void test3() 
	{

	}


	// 测试几种函数拟合方法：
	void test4() 
	{


	}
}

