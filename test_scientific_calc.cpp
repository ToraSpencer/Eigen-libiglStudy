#include "test_scientific_calc.h"

 
static std::string g_debugPath = "E:/";



/////////////////////////////////////////////////////////////////////////////////////////////////DEBUG接口：
namespace MY_DEBUG 
{
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

	template<typename DerivedV>
	static void debugWriteVers2D(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat2D(path, vers);
	}

	template<typename T>
	static void debugWriteMesh(const char* name, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteMeshMat(path, vers, tris);
	}
}
using namespace MY_DEBUG;


 
namespace TEST_SCIENTIFIC_CALC
{
	// 齐次坐标表示下的坐标变换；
	void test1()
	{
		Eigen::MatrixXf vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/材料/bunny.obj"); 

		// 1. 顶点使用齐次坐标来表示——一个点是一个四元列向量(w*x, w*y, w*z, w)，w通常为1，为0时表示无穷远处的点。
		Eigen::MatrixXf versHomo = vers2homoVersF(vers);
		dispMatBlock(vers, 0, 3, 0, 2);
		dispMatBlock(versHomo, 0, 3, 0, 3);
		objWriteHomoMeshMat("E:/meshInput.obj", versHomo, tris);


		// 2. 笛卡尔坐标下施加仿射变换——旋转、放缩，最后平移：
		Eigen::Matrix3f scale = Eigen::Matrix3f::Identity();
		Eigen::Matrix3f rotation;
		Eigen::Vector3f moveVec(10, 0, 0);															// 朝着x正方向平移20mm;
		Eigen::MatrixXf versDesc = vers.transpose();
		
		// scale(1, 1) = 2;																							// y方向上缩放因子为2；
		rotation << cos(pi / 6), 0, sin(pi / 6), 0, 1, 0, -sin(pi / 6), 0, cos(pi / 6);		// 绕z轴逆时针旋转pi/6度。
		versDesc = (rotation * scale * versDesc).eval();
		versDesc = (versDesc.colwise() + moveVec).eval();
		versDesc = versDesc.transpose().eval();
		objWriteMeshMat("E:/笛卡尔坐标下仿射变换后的bunny.obj", versDesc, tris);

		// 3. 齐次坐标下施加仿射变换——旋转、放缩， 最后平移：
		Eigen::Matrix4f scaleHomo = Eigen::Matrix4f::Identity();
		Eigen::Matrix4f rotationHomo = Eigen::Matrix4f::Identity();
		Eigen::Matrix4f moveMat = Eigen::Matrix4f::Identity();
		Eigen::MatrixXf finalVersHomo;
		scaleHomo.block<3, 3>(0, 0) = scale;
		rotationHomo.block<3, 3>(0, 0) = rotation;
		moveMat.topRightCorner(3, 1) = moveVec;
		 finalVersHomo = moveMat * rotationHomo * scaleHomo * versHomo;
		 objWriteHomoMeshMat("E:/齐次坐标下仿射变换后的bunny.obj", finalVersHomo, tris); 

		Eigen::MatrixXf transMat = moveMat * rotationHomo * scaleHomo;
		Eigen::MatrixXf originVersHomo = transMat.inverse() * finalVersHomo;
		objWriteHomoMeshMat("E:/齐次坐标下还原的bunny.obj", originVersHomo, tris); 

		std::cout << "finished" << std::endl;
	}


	// 测试霍纳方法（秦九昭算法）计算多项式：减少乘法调用次数，提高运行速率；
	void test2()
	{
		// x^4 + 3x^3 + 2x^2 + x + c = c + x(1 + x(2*x + x(3 + x)));
		double x = 3.0;
		Eigen::Vector4d coeff(9, -1, 3, 1);
		std::cout << "霍纳方法计算：result == " << hornersPoly(coeff, x) << std::endl;

		Eigen::RowVector4d coeffTrans = coeff.transpose();
		Eigen::Vector4d X(1, x, std::powf(x, 2), std::powf(x, 3));
		std::cout << "传统方法计算多项式：result == " << coeffTrans * X << std::endl;

		x = 5.2;
		Eigen::VectorXd coeff2(6);
		coeff2 << 2, 4, 5, 6, 8, 9;
		std::cout << "霍纳方法计算：result == " << hornersPoly(coeff2, x) << std::endl;

		Eigen::RowVectorXd coeffTrans2 = coeff2.transpose();
		Eigen::VectorXd X2(6);
		X2 << 1, x, std::powf(x, 2), std::powf(x, 3), std::powf(x, 4), std::powf(x, 5);
		std::cout << "传统方法计算多项式：result == " << coeffTrans2 * X2 << std::endl;

		debugDisp("finished.");
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
 
 





