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
		MatrixXf versHomo(4, vers.rows());
		versHomo.topRows(3) = vers.transpose().eval();
		versHomo.bottomRows(1) = RowVectorXf::Ones(vers.rows());

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

		objWriteMeshMat("./笛卡尔坐标下仿射变换后的bunny.obj", versDesc.transpose(), tris);


		// 3. 齐次坐标下施加仿射变换——旋转、放缩， 最后平移：
		Matrix4f scaleHomo = Matrix4f::Identity();
		scaleHomo.block<3, 3>(0, 0) = scale;
		Matrix4f rotationHomo = Matrix4f::Identity();
		rotationHomo.block<3, 3>(0, 0) = rotation;
		Matrix4f moveMat = Matrix4f::Identity();
		moveMat.topRightCorner(3, 1) = moveVec;
		versHomo = (moveMat * rotationHomo * scaleHomo * versHomo).eval();
		dispMat<float>(scaleHomo);
		dispMat<float>(rotationHomo);
		dispMat<float>(moveMat);
		dispMat<float>(moveMat * rotationHomo * scaleHomo);
		dispMatBlock<float>(versHomo, 0, 3, 0, 3);

		MatrixXf finalVers = versHomo.transpose().eval().leftCols(3);
		objWriteMeshMat("./齐次坐标下仿射变换后的bunny.obj", finalVers, tris);
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

}

