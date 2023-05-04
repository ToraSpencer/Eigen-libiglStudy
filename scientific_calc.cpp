#include "scientific_calc.h"


static igl::opengl::glfw::Viewer viewer;				// libigl中的基于glfw的显示窗口；
static std::string g_debugPath = "E:/";

/////////////////////////////////////////////////////////////////////////////////////////////////DEBUG接口：

static void debugDisp()			// 递归终止
{						//		递归终止设为无参或者一个参数的情形都可以。
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
	objWriteVerticesMat(path, vers);
}

 
template<typename T>
static void debugWriteMesh(const char* name, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteMeshMat(path, vers, tris);
}



namespace SCIENTIFICCALC
{
	// 齐次坐标表示下的坐标变换；
	void test1()
	{
		Eigen::MatrixXf vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "./data/bunny.obj");

		// 1. 顶点使用齐次坐标来表示——一个点是一个四元列向量(w*x, w*y, w*z, w)，w通常为1，为0时表示无穷远处的点。
		Eigen::MatrixXf versHomo;
		vers2homoVers(versHomo, vers);
		dispMatBlock(vers, 0, 3, 0, 2);
		dispMatBlock(versHomo, 0, 3, 0, 3);

		// 2. 笛卡尔坐标下施加仿射变换——旋转、放缩，最后平移：
		Eigen::Matrix3f scale = Eigen::Matrix3f::Identity();
		scale(1, 1) = 2;													// y方向上缩放因子为2；
		Eigen::Matrix3f rotation;
		rotation << cos(pi / 6), 0, sin(pi / 6), 0, 1, 0, -sin(pi / 6), 0, cos(pi / 6);	// 绕z轴逆时针旋转pi/6度。
		Eigen::Vector3f moveVec(10, 0, 0);			// 朝着x正方向平移20mm;
		Eigen::MatrixXf versDesc = vers.transpose();
		versDesc = (rotation * scale * versDesc).eval();
		versDesc = (versDesc.colwise() + moveVec).eval();
		objWriteMeshMat<float>("E:/笛卡尔坐标下仿射变换后的bunny.obj", versDesc.transpose(), tris);

		// 3. 齐次坐标下施加仿射变换——旋转、放缩， 最后平移：
		Eigen::Matrix4f scaleHomo = Eigen::Matrix4f::Identity();
		scaleHomo.block<3, 3>(0, 0) = scale;
		Eigen::Matrix4f rotationHomo = Eigen::Matrix4f::Identity();
		rotationHomo.block<3, 3>(0, 0) = rotation;
		Eigen::Matrix4f moveMat = Eigen::Matrix4f::Identity();
		moveMat.topRightCorner(3, 1) = moveVec;

		Eigen::MatrixXf finalVersHomo = moveMat * rotationHomo * scaleHomo * versHomo;
		objWriteMeshMat("E:/齐次坐标下仿射变换后的bunny.obj", homoVers2vers(finalVersHomo), tris);

		Eigen::MatrixXf transMat = moveMat * rotationHomo * scaleHomo;
		Eigen::MatrixXf originVersHomo = transMat.inverse() * finalVersHomo;
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

		m1 = Eigen::Map<Eigen::MatrixXd>(vec.data(), 10, 10);
		m2.resize(2, 2);
		m2 << 1, 1.1, 1, 1;
		Eigen::Matrix2d m22;
		m22 << 1, 1.1, 1, 1;

		kron(result, m1, m2);
		dispMat(m1);
		dispMat(result);

		kron(result, m1, m22);
		dispMat(result);

		Eigen::Vector3f v1(1, 2, 3);
		Eigen::Vector4f v2{ Eigen::Vector4f::Ones() };
		kron(result, v1, v2);
		dispMat(result);

		std::cout << "finished." << std::endl;
	}


	// 测试拟合椭圆：
	void test7()
	{
		// 1. 读取一圈牙齿网格：
		const unsigned teethCount = 14;
		std::vector<Eigen::MatrixXd> teethVers(teethCount);
		std::vector<Eigen::MatrixXi> teethTris(teethCount);
		Eigen::MatrixXd teethCenters(teethCount, 3);
		for (unsigned i = 0; i<teethCount; ++i) 
		{
			char fileName[256];
			sprintf_s(fileName, "E:/材料/teeth/compMesh%d.obj", i);
			objReadMeshMat(teethVers[i], teethTris[i], fileName);
			teethCenters.row(i) = teethVers[i].colwise().mean();
		}

		// 2. 牙齿中心点投影到XOY面上，作为样本点，拟合椭圆；
		teethCenters.col(2).array() = 0;
		Eigen::VectorXd coeffs = fittingStandardEllipse(teethCenters);
		double a = coeffs(0);
		double c = coeffs(1);
		double d = coeffs(2);
		double e = coeffs(3);
		double f = coeffs(4);
		double x0 = -d / (2 * a);		// 椭圆圆心，半长轴半短轴：
		double y0 = -e / (2 * c);
		double p = d * d / (4 * a) + e * e / (4 * c) - f;
		double a0 = std::sqrt(p / a);
		double b0 = std::sqrt(p / c);

		// 3. 生成椭圆点集：
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
	// 测试求法向的接口
	void test1() 
	{
		Eigen::MatrixXd vers, triNorms, barycenters, edgeArrows, vas, vbs;
		Eigen::MatrixXi tris, edges;
		objReadMeshMat(vers, tris, "E:/材料/tooth.obj");

		getEdges(edges, tris);
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
		const Eigen::RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);    // RGB色彩向量；
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
		nonManifoldUEs(nmnEdges, tris);
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
			nonManifoldUEs(nmnEdges, tris);

			trisSum += tris.rows();
			if (nmnEdges.rows() > 0) 
				std::cout << "splitedMesh" << i << "有非流形有向边" << nmnEdges.rows() << "条。" << std::endl;
		}
		std::cout << "trisSum == " << trisSum << std::endl;
		std::cout << "finished." << std::endl;


		objReadMeshMat(vers, tris, "E:/材料/meshNmnEdges.obj");
		nonManifoldUEs(nmnEdges, tris);
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





