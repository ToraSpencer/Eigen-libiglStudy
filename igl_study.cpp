#include "igl_study.h"

#define DATA_PATH "./data/"

// libigl��������
namespace IGL_BASIC
{
	Eigen::MatrixXd vers, newVers;
	Eigen::MatrixXi tris;
	Eigen::SparseMatrix<double> L;
	igl::opengl::glfw::Viewer viewer;		// libigl�еĻ���glfw����ʾ���ڣ�

	// �ļ�IO
	void test0()
	{
		// readOBJ(), writeObj()����OBJ�ļ���IO���ж�����أ�������·�������ơ�����Ƭ���������ݡ�
		igl::readOBJ("./data/bunny.obj", vers, tris);
		igl::writeOBJ("./data/bunny_export.obj", vers, tris);
		igl::writeOBJ("./data/bunnyVers.obj", vers, Eigen::MatrixXi{});			// ֻҪ���Ʋ���Ҫ����Ƭ�Ļ�������վ���

		return;
	}


	// libigl�е���ʾ������Viewer
	void test1() 
	{
		igl::readOBJ("bunny.obj", vers, tris);

		// Viewer::data()��������viewer�����ݶ�������ã�

		// ViewerData::set_mesh()�������붥�������Ƭ���������������ݣ�д�뵽ViewerData�ĳ�Ա�����У�
		viewer.data().set_mesh(vers, tris);		// ������װ�����ݣ�


		 // Ĭ����������ת��set_rotation_type()��������ָ��������ת���
		viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);				// ������ת

		viewer.data().show_lines = 0;			// show_linesָ���Ƿ񻭳������ߣ�

		viewer.launch();
	}

 
	// �����ݶ�
	void test2()
	{
		using namespace Eigen;
		using namespace std;

		igl::readOBJ("./data/rootTooth1.obj", vers, tris);

		VectorXd U;
		igl::readDMAT("./data/cheburashka-scalar.dmat", U);		// һϵ�еĺ���ֵ

		SparseMatrix<double> G;			// �ݶ�����
		igl::grad(vers, tris, G);

		// ���㺯��ֵU���ݶ�
		MatrixXd GU = Map<const MatrixXd>((G * U).eval().data(), tris.rows(), 3);

		// ����ֵU���ݶ�GU��ģ��
		const VectorXd GU_mag = GU.rowwise().norm();

		viewer.data().set_mesh(vers, tris);
		viewer.data().set_data(U);

		// Average edge length divided by average gradient (for scaling)
		const double max_size = igl::avg_edge_length(vers, tris) / GU_mag.mean();

		// ÿ������Ƭ�����ϻ�һ��ָʾ�ߣ�����Ϊ�ݶȷ��� 
		MatrixXd BC;
		igl::barycenter(vers, tris, BC);
		const RowVector3d black(0, 0, 0);
		viewer.data().add_edges(BC, BC + max_size * GU, black);
		viewer.data().show_lines = false;	  // ����������

		viewer.launch();
	}
 

	// lambda���������¼���ʹ��laplacian��˳����
	const auto& key_down = [](igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)->bool
	{
		switch (key)
		{
		case 'r':

		case 'R':			// ��λ����
			newVers = vers;
			break;

		case ' ':				// �ո����ִ��һ��laplace��˳
		{
			// ���¼�����������
			Eigen::SparseMatrix<double> mass;
			igl::massmatrix(newVers, tris, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

			// �����Է����� (mass - delta*L) * newVers = mass * newVers
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


	// ʹ��Laplacian��˳����
	void test3() 
	{
		 igl::readOBJ( "./data/bunny.obj", vers, tris);
		newVers = vers;

		// 1.a ֱ�ӹ���laplacian����Compute Laplace-Beltrami operator: 
		igl::cotmatrix(vers, tris, L);

		// 1.b �ֲ�����laplacian
		{
			SparseMatrix<double> Gradient, L2;

			igl::grad(vers, tris, Gradient);      // ��ɢ�ݶ�

			// Diagonal per-triangle "mass matrix"
			VectorXd dblA;
			igl::doublearea(vers, tris, dblA);             // ÿ������Ƭ���������

			// Place areas along diagonal  #dim times
			const auto& T = 1. * (dblA.replicate(3, 1) * 0.5).asDiagonal();

			L2 = -Gradient.transpose() * T * Gradient;         // discrete Dirichelet energy Hessian ��ɢ��������������������
			std::cout << "���ַ����õ���laplacian�Ĳ�ķ�����" << std::endl;
			cout << "(L2 - L).norm() == " << (L2 - L).norm() << endl;
		}

		// 2. ����ԭʼ�ķ�������ʹ��αɫ
		MatrixXd norms;
		igl::per_vertex_normals(vers, tris, norms);
		MatrixXd colors = norms.rowwise().normalized().array() * 0.5 + 0.5;

		// 3. viewr����ʼ����
		newVers = vers;
		viewer.data().set_mesh(newVers, tris);
		viewer.data().set_colors(colors);
		viewer.callback_key_down = key_down;

		// 4. ����
		cout << "Press [space] to smooth." << endl;;
		cout << "Press [r] to reset." << endl;;
		viewer.launch();
 
	}
 
}


// libigl�е�΢�ּ������
namespace IGL_DIF_GEO 
{

	void test1() 
	{



	}

}