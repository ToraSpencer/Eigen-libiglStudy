#include "scientific_calc.h"


static igl::opengl::glfw::Viewer viewer;				// libigl中的基于glfw的显示窗口；

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


	// 测试计算克罗内克积(Kronecker product);
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
	// 测试求法向的接口
	void test1() 
	{
		Eigen::MatrixXd vers, triNorms, barycenters, edgeArrows, vas, vbs;
		Eigen::MatrixXi tris, edges;
		objReadMeshMat(vers, tris, "E:/材料/tooth.obj");

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
		viewer.data().show_lines = false;                 // 隐藏网格线
		viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);    // 设定可以三轴旋转

		// 三角片法向指示线用边数据的形式渲染出来；
		const RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);    // RGB色彩向量；
		double aveLen = edgesLen.mean();							// 所有边长的平均值；
		viewer.data().add_edges(barycenters - aveLen * triNorms, barycenters + aveLen * triNorms, red);         // 最大曲率方向用红色指示线标识

		viewer.launch();

		std::cout << "finished." << std::endl;
	}


	// 计算三角网格的体积——
	void test2()
	{
		Eigen::MatrixXf vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/材料/jawMesh.obj");

		double volume = meshVolume(vers, tris);
		std::cout << "volume == " << volume << std::endl;

		std::cout << "finished." << std::endl;
	}


	// 计算三角网格的拓扑性质:
	void test3() 
	{
		// 求三角网格中的非流形半边：
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		Eigen::MatrixXi nmnEdges;

		objReadMeshMat(vers, tris, "E:/材料/jawMeshDense4.obj");
		nonManifoldEdges(tris, nmnEdges);
		objWriteEdgesMat("E:/nmnEdges.obj", nmnEdges, vers);

		std::cout << "trisCount == " << tris.rows() << std::endl;
		std::cout << "非流形有向边数：" << nmnEdges.rows() << std::endl;

		unsigned trisSum = 0;
		for (unsigned i = 0; i < 40; ++i) 
		{
			char str[512];
			sprintf_s(str, 512, "E:/网格精简/splittedData4/splitedMesh%d.obj", i);
			Eigen::MatrixXd vers;
			Eigen::MatrixXi tris;
			Eigen::MatrixXi nmnEdges;
			objReadMeshMat(vers, tris, str);
			nonManifoldEdges(tris, nmnEdges);

			trisSum += tris.rows();
			if (nmnEdges.rows() > 0) 
				std::cout << "splitedMesh" << i << "有非流形有向边" << nmnEdges.rows() << "条。" << std::endl;
		}
		std::cout << "trisSum == " << trisSum << std::endl;
		std::cout << "finished." << std::endl;


		objReadMeshMat(vers, tris, "E:/材料/meshNmnEdges.obj");
		nonManifoldEdges(tris, nmnEdges);
		objWriteEdgesMat("E:/nmnEdges.obj", nmnEdges, vers);
	}

}


// libigl中的数学工具
namespace IGL_MATH 
{
	// 计算特征值、特征向量：
	void test0() 
	{

	}

	// 克罗内克积：
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





