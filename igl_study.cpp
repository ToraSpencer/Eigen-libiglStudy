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
		//// readOBJ(), writeObj()——OBJ文件的IO，有多个重载，至少有路径、点云、三角片这三个数据。
		igl::readOBJ("./data/bunny.obj", vers, tris);
		igl::writeOBJ("./data/bunny_export.obj", vers, tris);
		igl::writeOBJ("./data/bunnyVers.obj", vers, Eigen::MatrixXi{});			// 只要点云不需要三角片的话，传入空矩阵；
	
		vers.resize(0, 0);
		tris.resize(0, 0);
		igl::readOBJ("E:/fatTeeth1_预处理后.obj", vers, tris);
		igl::writeOFF("E:/fatTeeth1_预处理后.off", vers, tris);

		vers.resize(0, 0);
		tris.resize(0, 0);
		igl::readOBJ("E:/fatTeeth2_预处理后.obj", vers, tris);
		igl::writeOFF("E:/fatTeeth2_预处理后.off", vers, tris);

		vers.resize(0, 0);
		tris.resize(0, 0);
		igl::readOBJ("E:/s9block_预处理后.obj", vers, tris);
		igl::writeOFF("E:/s9block_预处理后.off", vers, tris);


		std::cout << "finished." << std::endl;
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

	// libigl中的网格布尔操作；！！！需要配置CGAL库，当前未完成；
#if 0
	// https://stackoverflow.com/questions/37898084/how-to-perform-boolean-operation-on-thin-surface-using-libigl
	Eigen::MatrixXd VA, VB, VC;
	Eigen::VectorXi J, I;
	Eigen::MatrixXi FA, FB, FC;
	igl::MeshBooleanType boolean_type(
		igl::MESH_BOOLEAN_TYPE_UNION);

	const char* MESH_BOOLEAN_TYPE_NAMES[] =
	{
	  "Union",
	  "Intersect",
	  "Minus",
	  "XOR",
	  "Resolve",
	};

	const auto& key_down = [](igl::opengl::glfw::Viewer& viewer, unsigned char key, int mods) -> bool
	{
		switch (key)
		{
		default:
			return false;
		case 'A':
			viewer.data().clear();
			std::cout << "Loading A" << std::endl;
			viewer.data().set_mesh(VA, FA);
			break;
		case 'B':
			viewer.data().clear();
			std::cout << "Loading B" << std::endl;
			viewer.data().set_mesh(VB, FB);
			break;
		case 'C':
			viewer.data().clear();
			std::cout << "Loading C" << std::endl;
			viewer.data().set_mesh(VC, FC);
			return true;
		}

		return true;
	};

	void test1()
	{
		using namespace Eigen;
		using namespace std;

		double prismSize = 150;
		double Heigh = 300;
		VA.resize(6, 3);
		VA << -prismSize, prismSize, 0,
			prismSize, prismSize, 0,
			0, 2 * prismSize, 0,
			-prismSize, prismSize, Heigh,
			prismSize, prismSize, Heigh,
			0, 2 * prismSize, Heigh;
		FA.resize(8, 3);
		FA << 1, 0, 2,
			5, 3, 4,
			4, 1, 2,
			2, 5, 4,
			3, 5, 2,
			2, 0, 3,
			0, 1, 4,
			4, 3, 0;

		double tetsize = 300;
		VB.resize(4, 3);
		VB << 0, 0, tetsize,
			-tetsize, 0, 0,
			tetsize, 0, 0,
			0, tetsize * 2, 0;
		FB.resize(4, 3);
		FB << 2, 1, 3,
			2, 0, 1,
			3, 0, 2,
			1, 0, 3;


		igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, igl::MESH_BOOLEAN_TYPE_INTERSECT, VC, FC);

		std::cout
			<< "VA:" << std::endl << VA << std::endl << "==============" << std::endl
			<< "FA:" << std::endl << FA << std::endl << "==============" << std::endl
			<< "VB:" << std::endl << VB << std::endl << "==============" << std::endl
			<< "FB:" << std::endl << FB << std::endl << "==============" << std::endl
			<< "VC:" << std::endl << VC << std::endl << "==============" << std::endl
			<< "FC:" << std::endl << FC << std::endl << "==============" << std::endl;

		// Plot the mesh with pseudocolors
		igl::opengl::glfw::Viewer viewer;

		viewer.data().set_mesh(VA, FA);
 

		viewer.data().show_lines = 1;
		viewer.callback_key_down = key_down;
		viewer.core().camera_dnear = 3.9;
		cout <<
			"Press '.' to switch to next boolean operation type." << endl <<
			"Press ',' to switch to previous boolean operation type." << endl <<
			"Press ']' to push near cutting plane away from camera." << endl <<
			"Press '[' to pull near cutting plane closer to camera." << endl <<
			"Hint: investigate _inside_ the model to see orientation changes." << endl;
		viewer.launch();
	}
#endif

}