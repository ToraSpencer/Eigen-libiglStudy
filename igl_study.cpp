#include "igl_study.h"

#define DATA_PATH "./data/"

// libigl基本功能
namespace IGL_BASIC
{
	Eigen::MatrixXd vers, newVers;
	Eigen::MatrixXi tris;
	Eigen::SparseMatrix<double> L;
	igl::opengl::glfw::Viewer viewer;		// libigl中的基于glfw的显示窗口；

	// 文件IO
	void test0()
	{
		// readOBJ(), writeObj()——OBJ文件的IO，有多个重载，至少有路径、点云、三角片这三个数据。
		igl::readOBJ("./data/bunny.obj", vers, tris);
		igl::writeOBJ("./data/bunny_export.obj", vers, tris);
		igl::writeOBJ("./data/bunnyVers.obj", vers, Eigen::MatrixXi{});			// 只要点云不需要三角片的话，传入空矩阵；

		return;
	}


	// libigl中的显示窗口类Viewer
	void test1() 
	{
		igl::readOBJ("bunny.obj", vers, tris);

		// Viewer::data()——返回viewer中数据对象的引用；

		// ViewerData::set_mesh()——输入顶点和三角片矩阵生成网格数据，写入到ViewerData的成员变量中；
		viewer.data().set_mesh(vers, tris);		// 窗口中装载数据；


		 // 默认是两轴旋转，set_rotation_type()方法可以指定其他旋转风格
		viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);				// 三轴旋转

		viewer.data().show_lines = 0;			// show_lines指定是否画出网格线；

		viewer.launch();
	}

 
	// 计算梯度
	void test2()
	{
		using namespace Eigen;
		using namespace std;

		igl::readOBJ("./data/rootTooth1.obj", vers, tris);

		VectorXd U;
		igl::readDMAT("./data/cheburashka-scalar.dmat", U);		// 一系列的函数值

		SparseMatrix<double> G;			// 梯度算子
		igl::grad(vers, tris, G);

		// 计算函数值U的梯度
		MatrixXd GU = Map<const MatrixXd>((G * U).eval().data(), tris.rows(), 3);

		// 函数值U的梯度GU的模长
		const VectorXd GU_mag = GU.rowwise().norm();

		viewer.data().set_mesh(vers, tris);
		viewer.data().set_data(U);

		// Average edge length divided by average gradient (for scaling)
		const double max_size = igl::avg_edge_length(vers, tris) / GU_mag.mean();

		// 每个三角片重心上画一根指示线，方向为梯度方向。 
		MatrixXd BC;
		igl::barycenter(vers, tris, BC);
		const RowVector3d black(0, 0, 0);
		viewer.data().add_edges(BC, BC + max_size * GU, black);
		viewer.data().show_lines = false;	  // 隐藏网格线

		viewer.launch();
	}
 

	// lambda——键盘事件：使用laplacian光顺网格
	const auto& key_down = [](igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)->bool
	{
		switch (key)
		{
		case 'r':

		case 'R':			// 复位程序
			newVers = vers;
			break;

		case ' ':				// 空格键，执行一次laplace光顺
		{
			// 重新计算质量矩阵
			Eigen::SparseMatrix<double> mass;
			igl::massmatrix(newVers, tris, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

			// 解线性方程组 (mass - delta*L) * newVers = mass * newVers
			float delta = 0.001;
			const auto& S = (mass - delta * L);
			Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
			assert(solver.info() == Eigen::Success);
			newVers = solver.solve(mass * newVers).eval();
 
			break;
		}

		default:
			return false;
		}

		viewer.data().set_vertices(newVers);
		viewer.data().compute_normals();
		viewer.core().align_camera_center(newVers, tris);
		return true;
	};


	// 使用Laplacian光顺网格
	void test3() 
	{
		 igl::readOBJ( "./data/bunny.obj", vers, tris);
		newVers = vers;

		// 1.a 直接构造laplacian——Compute Laplace-Beltrami operator: 
		igl::cotmatrix(vers, tris, L);

		// 1.b 分步构造laplacian
		{
			SparseMatrix<double> Gradient, L2;

			igl::grad(vers, tris, Gradient);      // 离散梯度

			// Diagonal per-triangle "mass matrix"
			VectorXd dblA;
			igl::doublearea(vers, tris, dblA);             // 每个三角片面积的两倍

			// Place areas along diagonal  #dim times
			const auto& T = 1. * (dblA.replicate(3, 1) * 0.5).asDiagonal();

			L2 = -Gradient.transpose() * T * Gradient;         // discrete Dirichelet energy Hessian 离散狄利克雷能量海塞矩阵
			std::cout << "两种方法得到的laplacian的差的范数：" << std::endl;
			cout << "(L2 - L).norm() == " << (L2 - L).norm() << endl;
		}

		// 2. 根据原始的法向量，使用伪色
		MatrixXd norms;
		igl::per_vertex_normals(vers, tris, norms);
		MatrixXd colors = norms.rowwise().normalized().array() * 0.5 + 0.5;

		// 3. viewr填充初始数据
		newVers = vers;
		viewer.data().set_mesh(newVers, tris);
		viewer.data().set_colors(colors);
		viewer.callback_key_down = key_down;

		// 4. 运行
		cout << "Press [space] to smooth." << endl;;
		cout << "Press [r] to reset." << endl;;
		viewer.launch();
 
	}
 
}


// libigl中的微分几何相关
namespace IGL_DIF_GEO 
{

	void test1() 
	{



	}

}