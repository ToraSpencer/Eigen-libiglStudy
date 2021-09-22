#include "igl_study.h"


namespace IGLSTUDY
{
	// libigl中的文件IO
	void test0()
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOFF("./data/bunny.off", vers, tris);
		igl::writeOBJ("./data/bunny.obj", vers, tris);

	}

	
	// libigl中的显示窗口
	void test1() 
	{



	}


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

}
