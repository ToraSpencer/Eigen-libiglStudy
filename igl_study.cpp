#include "igl_study.h"

#define DATA_PATH "./data/"

namespace IGLSTUDY
{
	//// libigl中的文件IO
	//void test0()
	//{
	//	Eigen::MatrixXd vers;
	//	Eigen::MatrixXi tris;
	//	igl::readOFF("./data/bunny.off", vers, tris);
	//	igl::writeOBJ("./data/bunny.obj", vers, tris);

	//	dispMatBlock<double>(vers, 0, 10, 0, 2);
	//}

	void test0()
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXd tris;
		igl::readOBJ("./data/bunny.obj", vers, tris);

		return;
	}


	
	// libigl中的显示窗口类Viewer
	void test1() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("./data/bunny.obj", vers, tris);

		igl::opengl::glfw::Viewer viewer;				// libigl中的基于glfw的显示窗口；

		// Viewer::data()——返回viewer中数据对象的引用；

		// ViewerData::set_mesh()——输入顶点和三角片矩阵生成网格数据，写入到ViewerData的成员变量中；
		viewer.data().set_mesh(vers, tris);		// 窗口中装载数据；

		 // 默认是两轴旋转，加上下面这一句设定为三轴旋转。
		viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

		viewer.launch();
	}

#if 0
	void test2()
	{
		using namespace Eigen;
		using namespace std;


		MatrixXd V;
		MatrixXi F;
		igl::readOBJ("./data/rootTooth1.obj", V, F);

		VectorXd U;
		igl::readDMAT("./data/cheburashka-scalar.dmat", U);		// 一系列的函数值

		std::cout << "U == \n" << U << std::endl;

		SparseMatrix<double> G;			// 梯度算子
		igl::grad(V, F, G);


		// Compute gradient of U
		MatrixXd GU = Map<const MatrixXd>((G * U).eval().data(), F.rows(), 3);

		// Compute gradient magnitude
		const VectorXd GU_mag = GU.rowwise().norm();


		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(V, F);

		//viewer.data().set_data(U);

		//// Average edge length divided by average gradient (for scaling)
		//const double max_size = igl::avg_edge_length(V, F) / GU_mag.mean();


		//// 每个三角片重心上画一根指示线，方向为梯度方向。 
		//MatrixXd BC;
		//igl::barycenter(V, F, BC);
		//const RowVector3d black(0, 0, 0);
		//viewer.data().add_edges(BC, BC + max_size * GU, black);

		viewer.data().show_lines = false;	  // 隐藏网格线

		viewer.launch();

	}
#endif


	// 使用Laplacian平滑网格
	void test3() 
	{
		Eigen::MatrixXd vers, newVers;
		Eigen::MatrixXi tris;
		igl::readOBJ( "./data/mesh1.obj", vers, tris);
		newVers = vers;
		Eigen::SparseMatrix<double> L;
		igl::cotmatrix(vers, tris, L);

		int loopCount = 10;
		for (int i = 0; i<loopCount; ++i) 
		{
			// 重新计算质量矩阵
			SparseMatrix<double> mass;
			igl::massmatrix(newVers, tris, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

			// 解线性方程组 (mass - delta*L) * newVers = mass * newVers
			float delta = 0.01;
			const auto& S = (mass - delta * L);
			Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
			assert(solver.info() == Eigen::Success);
			newVers = solver.solve(mass * newVers).eval();

#if 0
			// Compute centroid and subtract (also important for numerics)
			VectorXd dblA;                                     // 每个三角片面积的两倍；
			igl::doublearea(newVers, tris, dblA);
			double areaSum = 0.5 * dblA.sum();
			MatrixXd centers;
			igl::barycenter(newVers, tris, centers);
			RowVector3d centroid(0, 0, 0);
			for (int i = 0; i < centers.rows(); i++)
			{
				centroid += 0.5 * dblA(i) / areaSum * centers.row(i);
			}
			newVers.rowwise() -= centroid;

			// 面积归一化
			newVers.array() /= sqrt(areaSum);
#endif
			vers = newVers;
		}


		igl::writeOBJ("./data/mesh1_平滑后.obj", vers, tris);

		std::cout << "finished." << std::endl;
	}
}
