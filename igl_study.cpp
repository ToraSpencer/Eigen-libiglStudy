#include "igl_study.h"

#define DATA_PATH "./data/"

namespace IGLSTUDY
{
	//// libigl�е��ļ�IO
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


	
	// libigl�е���ʾ������Viewer
	void test1() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("./data/bunny.obj", vers, tris);

		igl::opengl::glfw::Viewer viewer;				// libigl�еĻ���glfw����ʾ���ڣ�

		// Viewer::data()��������viewer�����ݶ�������ã�

		// ViewerData::set_mesh()�������붥�������Ƭ���������������ݣ�д�뵽ViewerData�ĳ�Ա�����У�
		viewer.data().set_mesh(vers, tris);		// ������װ�����ݣ�

		 // Ĭ����������ת������������һ���趨Ϊ������ת��
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
		igl::readDMAT("./data/cheburashka-scalar.dmat", U);		// һϵ�еĺ���ֵ

		std::cout << "U == \n" << U << std::endl;

		SparseMatrix<double> G;			// �ݶ�����
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


		//// ÿ������Ƭ�����ϻ�һ��ָʾ�ߣ�����Ϊ�ݶȷ��� 
		//MatrixXd BC;
		//igl::barycenter(V, F, BC);
		//const RowVector3d black(0, 0, 0);
		//viewer.data().add_edges(BC, BC + max_size * GU, black);

		viewer.data().show_lines = false;	  // ����������

		viewer.launch();

	}
#endif


	// ʹ��Laplacianƽ������
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
			// ���¼�����������
			SparseMatrix<double> mass;
			igl::massmatrix(newVers, tris, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

			// �����Է����� (mass - delta*L) * newVers = mass * newVers
			float delta = 0.01;
			const auto& S = (mass - delta * L);
			Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
			assert(solver.info() == Eigen::Success);
			newVers = solver.solve(mass * newVers).eval();

#if 0
			// Compute centroid and subtract (also important for numerics)
			VectorXd dblA;                                     // ÿ������Ƭ�����������
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

			// �����һ��
			newVers.array() /= sqrt(areaSum);
#endif
			vers = newVers;
		}


		igl::writeOBJ("./data/mesh1_ƽ����.obj", vers, tris);

		std::cout << "finished." << std::endl;
	}
}
