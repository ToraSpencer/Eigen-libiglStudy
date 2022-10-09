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
		dispMatBlock(vers, 0, 3, 0, 2);
		dispMatBlock(versHomo, 0, 3, 0, 3);

		// 2. 笛卡尔坐标下施加仿射变换——旋转、放缩，最后平移：
		Matrix3f scale = Matrix3f::Identity();
		scale(1, 1) = 2;													// y方向上缩放因子为2；
		Matrix3f rotation;
		rotation << cos(pi / 6), 0, sin(pi / 6), 0, 1, 0, -sin(pi / 6), 0, cos(pi / 6);	// 绕z轴逆时针旋转pi/6度。
		Vector3f moveVec(10, 0, 0);			// 朝着x正方向平移20mm;
		MatrixXf versDesc = vers.transpose();
		versDesc = (rotation * scale * versDesc).eval();
		versDesc = (versDesc.colwise() + moveVec).eval();
		objWriteMeshMat<float>("E:/笛卡尔坐标下仿射变换后的bunny.obj", versDesc.transpose(), tris);

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


	// 主成分分析：
	void test5() 
	{


	}


	// 计算三角网格的体积——
	void test6() 
	{
		Eigen::MatrixXf vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/材料/jawMesh.obj");
 
		double volume = meshVolume(vers, tris);
		std::cout << "volume == " << volume << std::endl;

		std::cout << "finished." << std::endl;
	}


	// 测试计算克罗内克积(Kronecker product);
	void test7() 
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



// libigl中的数学工具
namespace IGL_MATH 
{
	// 计算特征值、特征向量：
	void test0() 
	{

	}
}

