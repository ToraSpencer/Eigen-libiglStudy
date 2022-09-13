#include "igl_study.h"

#define DATA_PATH "./data/"


// libigl��������
namespace IGL_BASIC
{
	Eigen::MatrixXd vers, newVers, normals;
	Eigen::MatrixXi tris;
	Eigen::SparseMatrix<double> L;
	igl::opengl::glfw::Viewer viewer;		// libigl�еĻ���glfw����ʾ���ڣ�

	// �ļ�IO
	void test0()
	{
		//// readOBJ(), writeObj()����OBJ�ļ���IO���ж�����أ�������·�������ơ�����Ƭ���������ݡ�
		igl::readOBJ("./data/bunny.obj", vers, tris);
		igl::writeOBJ("./data/bunny_export.obj", vers, tris);
		igl::writeOBJ("./data/bunnyVers.obj", vers, Eigen::MatrixXi{});			// ֻҪ���Ʋ���Ҫ����Ƭ�Ļ�������վ���
	
		//vers.resize(0, 0);
		//tris.resize(0, 0);
		//igl::readOBJ("E:/fatTeeth1_Ԥ�����.obj", vers, tris);
		//igl::writeOFF("E:/fatTeeth1_Ԥ�����.off", vers, tris);

		// ��ȡstl�ļ���
		std::string fileName{"E:/����/jawMesh_noHP_noPro"};
		vers.resize(0, 0);
		tris.resize(0, 0);
		std::ifstream fileIn((fileName + std::string{ ".stl" }).c_str(), std::ios::binary);			// stl�ļ��Ƕ������ļ���
		igl::readSTL(fileIn, vers, tris, normals);
		fileIn.close();
		igl::writeOBJ((fileName + std::string{".obj"}).c_str(), vers, tris);
 

		std::cout << "finished." << std::endl;
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


	// ������ά�ռ��е�դ��
	void test4() 
	{
		Eigen::MatrixXd vers, gridCenters;
		Eigen::MatrixXi tris;	
		Eigen::RowVector3i gridCounts;			// ����ά����դ�����Ŀ��
		unsigned num = 10;				// ��������Ǹ�ά��(xyz�е�һ��)��դ������

		igl::readOBJ("E:/����/tooth.obj", vers, tris);
		igl::voxel_grid(vers, 0, num, 1, gridCenters, gridCounts);

		igl::writeOBJ("E:/gridCenters.obj", gridCenters, Eigen::MatrixXi{});

		std::cout << "finished." << std::endl;
	}


	// ������ž��볡��marching cubes��ȡ��ֵ������
	void test5()
	{
		MatrixXi tris;
		MatrixXd vers;
		igl::readOBJ("E:/����/jawMesh.obj", vers, tris);

		tiktok& tt = tiktok::getInstance();
		double gridStep = 0.5;			
		double range[3];
		for (unsigned i = 0; i < 3; ++i)
		{
			Eigen::VectorXd coors = vers.col(i);
			range[i] = coors.maxCoeff() - coors.minCoeff();
		}

		// 1. ����դ��
		double largestRange = std::max({ range[0], range[1], range[2] });
		int num = std::ceil((largestRange/gridStep));              // ��������Ǹ�ά��(xyz�е�һ��)��դ������
		double gridsOffset = (gridStep * num - largestRange) / 2.;
		
		MatrixXd gridCenters;
		Eigen::RowVector3i gridCounts;
		tt.start();
		igl::voxel_grid(vers, gridsOffset, num, 1, gridCenters, gridCounts);
		tt.endCout("Elapsed time of igl::voxel_grid() is ");

		VectorXd SDF, signValues;

		{
			VectorXi I;                     // useless
			MatrixXd C, N;              // useless

			// 2. ������ž��볡
			tt.start();
			igl::signed_distance(gridCenters, vers, tris, igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, SDF, I, C, N);
			tt.endCout("Elapsed time of igl::signed_distance() is ");

			// 3. ���ž��볡��дΪ���ų�����������Ϊ-1����������Ϊ0������Ϊ1��
			signValues = SDF;
			for_each(signValues.data(), signValues.data() + signValues.size(), [](double& b)\
			{
				b = (b > 0 ? 1 : (b < 0 ? -1 : 0));
			});
		}

		// 4. marching cubes�㷨�����������棺
		MatrixXd versResult_SDF, versResults_signs;
		MatrixXi trisResult_SDF, trisResults_signs;
		double selectedSDF = -1.;
		tt.start();
		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), selectedSDF, versResult_SDF, trisResult_SDF);
		tt.endCout("Elapsed time of igl::marching_cubes() is ");
		
		// igl::marching_cubes(signValues, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), 0, versResults_signs, trisResults_signs);
		igl::writeOBJ("E:/shrinkedMesh.obj", versResult_SDF, trisResult_SDF);


		std::cout << "finished." << std::endl;
	}



	// ��ȡSDFGen.exe���ɵľ��볡���ݣ�ʹ��igl::marching_cubes()��ȡ��ֵ������
	void test55() 
	{
		std::vector<int> stepCouts(3);				// xyz����ά����դ����
		Eigen::RowVector3d gridsOri;					// դ��ԭ�㣺
		std::ifstream sdfFile("E:/inputMesh.sdf");
		std::string readStr(1024, '\0');

		// ��һ�У�
		{
			std::string tmpStr;
			sdfFile.getline(&readStr[0], 1024);

			unsigned index = 0;
			for (const auto& ch : readStr)
			{
				if (ch >= '0' && ch <= '9' || ch == '.')
					tmpStr.push_back(ch);
				else
				{
					if (tmpStr.size() > 0)
					{
						stepCouts[index] = std::stoi(tmpStr);
						index++;
						tmpStr.clear();
					}
				}
			}
		}

		// �ڶ��У�
		{
			std::string tmpStr;
			sdfFile.getline(&readStr[0], 1024);

			unsigned index = 0;
			for (const auto& ch : readStr)
			{
				if (ch >= '0' && ch <= '9' || ch == '.')
					tmpStr.push_back(ch);
				else
				{
					if (tmpStr.size() > 0)
					{
						gridsOri(index) = std::stod(tmpStr);
						index++;
						tmpStr.clear();
					}
				}
			}
		}

		// �ڶ���֮�󣬾��볡���ݣ�
		


		sdfFile.close();


		traverseSTL(stepCouts, disp<int>);
		dispVec<double, 3>(gridsOri);
		std::cout << "finished." << std::endl;
	}
}


// libigl�е�΢�ּ������
namespace IGL_DIF_GEO 
{

	// libigl�е����񲼶���������������Ҫ����CGAL�⣬��ǰδ��ɣ�
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


// ͼ�㷨
namespace IGL_GRAPH 
{
	// ͼ���ݽṹ��ת����
	void test0() 
	{
		MatrixXd vers;
		MatrixXi tris;
		igl::readOBJ("E:/����/cube�ظ�����Ƭ.obj", vers, tris);

		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		
		// ����Ƭ��Ϣ�õ�ͼ���ڽӱ��ڽӾ���
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);		// ����ߴ�����Ԫ��Ϊ1������Ϊ0��

		// �ֶ������ڽӾ���
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		unsigned edgesCount = 3 * trisCount;
		Eigen::SparseMatrix<int> adjSM_eCount;				// ������ڽӾ���Ȩ��Ϊ������ߵ�������

		MatrixXi edges = MatrixXi::Zero(edgesCount, 2);
		edges.block(0, 0, trisCount, 1) = tris.col(0);
		edges.block(0, 1, trisCount, 1) = tris.col(1);
		edges.block(trisCount, 0, trisCount, 1) = tris.col(1);
		edges.block(trisCount, 1, trisCount, 1) = tris.col(2);
		edges.block(2 * trisCount, 0, trisCount, 1) = tris.col(2);
		edges.block(2 * trisCount, 1, trisCount, 1) = tris.col(0);

		// FOR DEBUG:
		dispMatBlock<int>(edges, 0, 10, 0, 1);
		dispMatBlock<int>(edges, 3 * trisCount - 10, 3 * trisCount - 1, 0, 1);

		std::vector<Eigen::Triplet<int>> elems;
		elems.reserve(2 * (versCount + trisCount + 100));				// euler formula: E = V + S - 2 + 2*g;
		for (unsigned i = 0; i < edgesCount; ++i)
			elems.push_back(Eigen::Triplet<int>(edges(i, 0), edges(i, 1), 1));

		adjSM_eCount.resize(versCount, versCount);
		adjSM_eCount.reserve(elems.size());
		adjSM_eCount.setFromTriplets(elems.begin(), elems.end());
		adjSM_eCount.makeCompressed();

		std::cout << "����Ϊ0�Ķ����1���򶥵㣺" << std::endl;
		traverseSTL(adjList[0], disp<int>);
		std::cout << std::endl << std::endl;

		std::cout << "����Ϊ1�Ķ������ߣ����ڽӾ����е�1�������в�Ϊ���Ԫ���±꣩��" << std::endl;
		for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM, 0); it; ++it)
		{
			std::cout << "value == " << it.value() << std::endl;
			std::cout << "row == " << it.row() << std::endl;			 // row index
			std::cout << "col == " << it.col() << std::endl;			 // col index (here it is equal to k)
			std::cout << std::endl;
		} 

		std::cout << "����Ϊ2�Ķ���ĳ��ߣ����ڽӾ����е�2�������в�Ϊ���Ԫ���±꣩��" << std::endl;
		for (int i = 0; i < adjSM_eCount.outerSize(); ++i)			// ���еı�����
		{
			for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM_eCount, i); it; ++it)		// ��ǰ�е����ڵ�������
			{
				if (2 == it.row())
				{
					std::cout << "value == " << it.value() << std::endl;
					std::cout << "row == " << it.row() << std::endl;			 // row index
					std::cout << "col == " << it.col() << std::endl;			 // col index (here it is equal to k)
				}
			}
		}

		// �ж��Ƿ����ظ�����Ƭ���������ظ�����Ƭ��������ظ�����ߣ�
		bool hasDupTri = false;;
		for (int i = 0; i < adjSM_eCount.outerSize(); ++i)			// ���еı�����
		{
			if (hasDupTri)
				break;
			for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM_eCount, i); it; ++it)		// ��ǰ�е����ڵ�������
			{
				if (it.value() > 1)
				{
					std::cout << "�����ظ�����Ƭ" << std::endl;
					hasDupTri = true;
					break;
				}
			}
		}
		if (!hasDupTri)
			std::cout << "û���ظ�����Ƭ" << std::endl;

		std::cout << "finished." << std::endl;
	}


	// ������ͼ�㷨��
	void test1() 
	{
		MatrixXd vers;
		MatrixXi tris;
		igl::readOBJ("E:/����/roundSurf.obj", vers, tris);
 
		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);

		// dfs:
		Eigen::VectorXi disCoveredIdx, bfsTreeVec, dfsTreeVec, closedIdx;
		size_t startIdx = 165;			// �ӽ�Բ�ĵĶ��㣻
		igl::dfs(adjList, startIdx, disCoveredIdx, dfsTreeVec, closedIdx);
		objWriteTreePath("E:/dfsTree.obj", bfsTreeVec, vers);
		auto retCrev = closedIdx.reverse();
		auto flag = (retCrev == disCoveredIdx);

		std::vector<int> vec1, vec2, vec3;
		vec1 = vec2Vec<int>(disCoveredIdx);
		vec2 = vec2Vec<int>(closedIdx);
		vec3 = vec2Vec<int>(retCrev);

		// bfs:
		igl::bfs(adjList, startIdx, disCoveredIdx, bfsTreeVec);
		objWriteTreePath("E:/bfsTree.obj", bfsTreeVec, vers);

		// myDfs:
		std::vector<bool> visited(vers.rows(), false);
		std::list<int> discoveredVersIdx;				// �ѱ����ʵĶ��㣺
		std::list<int> closedVersIdx;

		std::function<void(const int, const int)> myDfs = [&](const int index, const int parentIdx)
		{
			// �ݹ���ֹ��
			if (visited[index])			
				return;

			visited[index] = true;
			discoveredVersIdx.push_back(index);
			const std::vector<int>& adjVersIdx = adjList[index];

			// �ݹ���ƣ�
			for (const auto& adjIdx : adjVersIdx)			
				myDfs(adjIdx, index);

			closedVersIdx.push_back(index);
		};


 
		std::cout << "finished." << std::endl;
	}


	// prime, dijkstra
	void test2() 
	{
		// dijkstra
		MatrixXd vers;
		MatrixXi tris;
		igl::readOBJ("E:/����/roundSurf.obj", vers, tris);

		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);

		Eigen::VectorXd min_distance;				// ͼ�����ж��㵽ָ����������·�����ȣ�
		Eigen::VectorXi mst;						// ��С��������
		int verIdx0 = 165;						// �����Ϊ�ӽ�Բ�ĵĶ��㣻

		// ��ָ���ص�ʱ�� igl::dijkstra()����ͼ����verIdx0Ϊ������С��������
		int retIdx = igl::dijkstra(vers, adjList, verIdx0, std::set<int>{}, min_distance, mst);
		std::cout << "retIdx == " << retIdx << std::endl;
		objWriteTreePath("E:/mst.obj", mst, vers);


		std::cout << "finished." << std::endl;
	}
 
}


// �ռ仮��
namespace IGL_SPACE_PARTITION
{
	// aabb��ʵ�ֵ�BVH(��ΰ�Χ��)
	void test0()
	{
		Eigen::MatrixXd vers, vers0, minDisVers;
		Eigen::MatrixXi tris;
		Eigen::VectorXd minSqrDis;
		Eigen::VectorXi minDisIdx;
		igl::readOBJ("E:/����/tooth.obj", vers, tris);

		vers0.resize(2, 3);
		vers0.row(0) = vers.row(0);
		vers0.row(1) = vers.row(99);

		// igl::point_mesh_squared_distance()����ʹ��BVH�󶥵㵽������С���룻
		igl::point_mesh_squared_distance(vers0, vers, tris, minSqrDis, minDisIdx, minDisVers);
		igl::writeOBJ("E:/inputMesh.obj", vers, tris);
		igl::writeOBJ("E:/vers0.obj", vers0, Eigen::MatrixXi{});
		igl::writeOBJ("E:/minDisVers.obj", minDisVers, Eigen::MatrixXi{});

		//// ���������BVH����
		//igl::AABB<double, 3> bvh;

		

		std::cout << "finished." << std::endl;
	}

}


// IGLʵ�ֵĻ��������������㷨��
namespace IGL_BASIC_PMP 
{
	// ��std::pair<int, int>���Զ����ϣ������
	struct edgeHash
	{
		bool operator()(const std::pair<int, int>& edge) const
		{
			return (std::hash<int>()(edge.first) + std::hash<int>()(edge.second));
		}
	};

	// ��std::pair<int, int>�ĵȼ۱Ƚ�����
	struct edgeComparator 
	{
		bool operator()(const std::pair<int, int>& edge1, const std::pair<int, int>& edge2) const
		{
			return (edge1.first == edge2.first && edge1.second == edge2.second);
		}
	};


	// �����������еı�-����Ƭ�ڽӹ�ϵ��
	void test1() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris, edges;
		igl::readOBJ("E:/����/meshArranged_hole.obj", vers, tris);
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		unsigned edgesCount = 3 * trisCount;		// ������������edgesCount����ʵ���������Ҫ����Ϊ�����ظ���ͬһ����߽������ظ�������

		// edges == [ea; eb; ec] == [vbIdxes, vcIdxes; vcIdxes, vaIdxes; vaIdxes, vbIdxes];
		/*
			������������������edges�б���ÿ���߶���unique�ģ�
			�����ڷ����αߣ�������α���edges����ظ��洢��ʵ�ʱ�����Ҳ������edges������
			����˵ʹ������representation�Ļ��������α��в�ֹһ��������ȡ�������÷����αߵ�����Ƭ������
		*/
		edges = Eigen::MatrixXi::Zero(edgesCount, 2);
		Eigen::MatrixXi vaIdxes = tris.col(0);
		Eigen::MatrixXi vbIdxes = tris.col(1);
		Eigen::MatrixXi vcIdxes = tris.col(2);
		edges.block(0, 0, trisCount, 1) = vbIdxes;
		edges.block(trisCount, 0, trisCount, 1) = vcIdxes;
		edges.block(trisCount * 2, 0, trisCount, 1) = vaIdxes;
		edges.block(0, 1, trisCount, 1) = vcIdxes;
		edges.block(trisCount, 1, trisCount, 1) = vaIdxes;
		edges.block(trisCount * 2, 1, trisCount, 1) = vbIdxes;

		// 1. �ڽӾ���
		std::vector<Eigen::Triplet<int>> smElems;
		smElems.reserve(edgesCount);
		for (unsigned i = 0; i < edgesCount; ++i)
			smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});

		Eigen::SparseMatrix<int> adjSM_eCount;
		adjSM_eCount.resize(versCount, versCount);
		adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());		// Ȩ��Ϊ��������ظ��Ĵ�����

		// 2. ������-����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
		std::vector<int> etInfo(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

#if __cplusplus < 202002L
		// 3. ȷ�������αߣ�
		std::unordered_set<std::pair<int, int>, edgeHash, edgeComparator> edgesNM;
#else
		// lambda����std::pair<int, int>��ʾ�ıߵ��Զ����ϣ������
		auto edgeHashLamb = [](const std::pair<int, int>& edge)->std::size_t
		{
			return (std::hash<int>()(edge.first) + std::hash<int>()(edge.second));
		};

		// lambda����std::pair<int, int>��ʾ�ıߵĵȼ۱Ƚ�����
		auto edgeComLamb = [](const std::pair<int, int>& edge1, const std::pair<int, int>& edge2)->bool
		{
			return (edge1.first == edge2.first && edge1.second == edge2.second);
		};

		// 3. ȷ�������αߣ�
		std::unordered_set<std::pair<int, int>, decltype(edgeHashLamb), decltype(edgeComLamb)> edgesNM;				// ��ҪC++20;�����αߣ�ʹ��unordered_setȥ�أ�
#endif
		 for (unsigned i = 0; i < adjSM_eCount.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_eCount, i); iter; ++iter)	// ��i�е����ڵ�������
			{
				if (iter.value() > 1)
				{
					int num1 = iter.row();
					int num2 = iter.col();
					edgesNM.insert({num1, num2});
				}
			}
		}

		
#if __cplusplus < 202002L
		 //		�����������-�������ļ�ֵ�ԣ�һ�������α߶�Ӧ�Ŷ����������
		 std::unordered_map<std::pair<int, int>, std::vector<int>, edgeHash, edgeComparator> edgesNMmap;

		 //		�����������-�����ڵ�����Ƭ�������ļ�ֵ�ԣ�һ�������α߶�Ӧ�Ŷ������Ƭ������
		 std::unordered_map<std::pair<int, int>, std::vector<int>, edgeHash, edgeComparator> etNMmap;
#else
		 //		�����������-�������ļ�ֵ�ԣ�һ�������α߶�Ӧ�Ŷ����������
		 std::unordered_map<std::pair<int, int>, std::vector<int>, decltype(edgeHashLamb), decltype(edgeComLamb)> edgesNMmap;
		 //		�����������-�����ڵ�����Ƭ�������ļ�ֵ�ԣ�һ�������α߶�Ӧ�Ŷ������Ƭ������
		 std::unordered_map<std::pair<int, int>, std::vector<int>, decltype(edgeHashLamb), decltype(edgeComLamb)> etNMmap;
#endif

		for (auto& nmEdge : edgesNM)
		{
			for (int i = 0; i < edgesCount; ++i)
			{
				if (edges(i, 0) == nmEdge.first && edges(i, 1) == nmEdge.second)
				{
					auto retPair = edgesNMmap.insert({nmEdge, std::vector<int>{i} });
					if (!retPair.second)			// ������ʧ�ܣ���˵�����д˼���
					{
						auto iter = edgesNMmap.find(nmEdge);
						iter->second.push_back(i);
					}
				}
			}
		}
		for (const auto& pair : edgesNMmap)
		{
			auto copyPair = pair;
			for (auto& index : copyPair.second)
				index = etInfo[index];
			etNMmap.insert(copyPair);
		}

		// 4. �ҳ����з����α���������Ƭ��
		std::vector<int> trisIdxNM;
		for (auto& pair : etNMmap)
		{
			for (auto& index : pair.second)
				trisIdxNM.push_back(index);
		}
		Eigen::MatrixXi trisNM(Eigen::MatrixXi::Zero(trisIdxNM.size(), 3));
		for (int i = 0; i < trisIdxNM.size(); ++i)
			trisNM.row(i) = tris.row(trisIdxNM[i]);

		// for debug: ��ӡtrisNM:
		objWriteMeshMat("E:/trisNM.obj", vers, trisNM);
		std::cout << "finished." << std::endl;

		// 5. ȷ����������Ƭ������������Ƭ���ڽӹ�ϵ��


	}

	using ttTuple = std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>;

	bool buildAdjacency(const Eigen::MatrixXi& tris, Eigen::MatrixXi& ttAdj_nmEdge, \
		std::vector<ttTuple>& ttAdj_nmnEdge, std::vector<ttTuple>& ttAdj_nmnOppEdge)
	{
		/*
			bool buildAdjacency(
						const Eigen::MatrixXi& tris,												���������Ƭ����
						Eigen::MatrixXi& ttAdj_nmEdge,										����Ƭ�ķǱ�Ե����������ڽӵ�����Ƭ������
						std::vector<ttTuple>& ttAdj_nmnEdge,							����Ƭ�ķ�������������ڵ�����Ƭ������
						std::vector<ttTuple>& ttAdj_nmnOppEdge					����Ƭ�ķ�����������ڽӵ�����Ƭ���������Ա����ڵ�����Ƭ����������
						)

		*/
		const unsigned trisCount = tris.rows();
		const unsigned edgesCount = 3 * trisCount;
		const unsigned versCount = tris.maxCoeff() + 1;

		// 1. �󶥵��ڽӹ�ϵ��

		// edges == [ea; eb; ec] == [vbIdxes, vcIdxes; vcIdxes, vaIdxes; vaIdxes, vbIdxes];
		/*
			������������������edges�б���ÿ���߶���unique�ģ�
			�����ڷ����αߣ�������α���edges����ظ��洢��ʵ�ʱ�����Ҳ������edges������
			����˵ʹ������representation�Ļ��������α��в�ֹһ��������ȡ�������÷����αߵ�����Ƭ������
		*/

		// ����Ƭ�����ߵ�������teIdx == [eaIdx, ebIdx, ecIdx] == [(0: trisCount-1)', (trisCount: 2*trisCount-1)', (2*trisCount, 3*trisCount-1)'];
		/*
			teIdx(i, j) = trisCount *j + i;
		*/
		Eigen::MatrixXi edges = Eigen::MatrixXi::Zero(edgesCount, 2);
		Eigen::MatrixXi vaIdxes = tris.col(0);
		Eigen::MatrixXi vbIdxes = tris.col(1);
		Eigen::MatrixXi vcIdxes = tris.col(2);
		edges.block(0, 0, trisCount, 1) = vbIdxes;
		edges.block(trisCount, 0, trisCount, 1) = vcIdxes;
		edges.block(trisCount * 2, 0, trisCount, 1) = vaIdxes;
		edges.block(0, 1, trisCount, 1) = vcIdxes;
		edges.block(trisCount, 1, trisCount, 1) = vaIdxes;
		edges.block(trisCount * 2, 1, trisCount, 1) = vbIdxes;

		std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;
		smElems.reserve(edgesCount);
		smElems_weighted.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
		{
			smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});
			smElems_weighted.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), i});
		}

		Eigen::SparseMatrix<int> adjSM_eCount, adjSM_weighted;
		adjSM_eCount.resize(versCount, versCount);
		adjSM_weighted.resize(versCount, versCount);
		adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());		// Ȩ��Ϊ��������ظ��Ĵ�����
		adjSM_weighted.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());		// Ȩ��Ϊ��������ظ��Ĵ���(�����αߵ�ֵ��Ч)
		Eigen::SparseMatrix<int> adjSM_weighted_opp = adjSM_weighted.transpose();

		Eigen::SparseMatrix<int> adjSM = adjSM_eCount;		// ������ڽӾ���
		for (unsigned i = 0; i < adjSM.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, i); iter; ++iter)
			{
				if (iter.value() > 0)
					adjSM.coeffRef(iter.row(), iter.col()) = 1;
			}
		}


		// 2. ȷ�����зǱ�Ե����������ߣ�
		Eigen::SparseMatrix<int> adjSM_eCount_ND = adjSM_eCount + Eigen::SparseMatrix<int>(adjSM_eCount.transpose());		// Ȩ��Ϊ�����ij����������Ƭ������

		//		�Ǳ�Ե��������ߣ�
		Eigen::SparseMatrix<int> adjSM_MNnonBdry_ND = adjSM_eCount_ND;
		for (unsigned i = 0; i < adjSM_MNnonBdry_ND.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_MNnonBdry_ND, i); iter; ++iter)
			{
				if (iter.value() > 0)
				{
					if (iter.value() == 2)
						adjSM_MNnonBdry_ND.coeffRef(iter.row(), iter.col()) = 1;
					else
						adjSM_MNnonBdry_ND.coeffRef(iter.row(), iter.col()) = 0;
				}
			}
		}

		//		�Ǳ�Ե����������ڽӾ���
		Eigen::SparseMatrix<int> adjSM_MNnonBdry = adjSM + adjSM_MNnonBdry_ND;		// adjSM & adjSM_MNnonBdry_ND
		adjSM_MNnonBdry.prune([](const Index& row, const Index& col, const float& value)->bool
			{
				if (2 == value)
					return true;
				else
					return false;
			});
		adjSM_MNnonBdry /= 2;

		unsigned edgesCount_MNnonBdry = adjSM_MNnonBdry.sum();
		Eigen::MatrixXi edges_MNnonBdry(edgesCount_MNnonBdry, 2);		// �Ǳ�Ե���������
		std::vector<int> edgesIdx_MNnonBdry;												// �Ǳ�Ե��������ߵ�������
		edgesIdx_MNnonBdry.reserve(edgesCount_MNnonBdry);
		unsigned index = 0;
		for (unsigned i = 0; i < adjSM_MNnonBdry.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_MNnonBdry, i); iter; ++iter)
			{
				edges_MNnonBdry(index, 0) = iter.row();
				edges_MNnonBdry(index, 1) = iter.col();
				edgesIdx_MNnonBdry.push_back(adjSM_weighted.coeffRef(iter.row(), iter.col()));
				index++;
			}
		}

		Eigen::SparseMatrix<int> adjSM_MNnonBdry_opp = adjSM_MNnonBdry.transpose();
		Eigen::MatrixXi edges_MNnonBdry_opp(edgesCount_MNnonBdry, 2);		// �Ǳ�Ե��������ߵĶԱߣ�
		std::vector<int> edgesIdx_MNnonBdry_opp;								// �Ǳ�Ե��������ߵĶԱߵ�������
		edgesIdx_MNnonBdry_opp.reserve(edgesCount_MNnonBdry);
		index = 0;
		for (unsigned i = 0; i < adjSM_MNnonBdry_opp.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_MNnonBdry_opp, i); iter; ++iter)
			{
				edges_MNnonBdry_opp(index, 0) = iter.row();
				edges_MNnonBdry_opp(index, 1) = iter.col();
				edgesIdx_MNnonBdry_opp.push_back(adjSM_weighted_opp.coeffRef(iter.row(), iter.col()));
				index++;
			}
		}

		// 3. ���α�-����Ƭ�ڽӹ�ϵ������Ƭ�ڽӹ�ϵ��
		std::vector<int> etInfo(edgesCount);				// ������ - ����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

		Eigen::VectorXi etAdj_mnEdge(-Eigen::VectorXi::Ones(edgesCount));	// ����������Ƭ�������������α߻��Ե��д-1��
		for (unsigned i = 0; i < edgesIdx_MNnonBdry.size(); ++i)
		{
			const int& edgeIdx = edgesIdx_MNnonBdry[i];
			const int& edgesIdxOpp = edgesIdx_MNnonBdry_opp[i];
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//			����Ƭ�ڽӾ���ttAdj_mnEdge(i, :)������Ϊi������Ƭ�������ڽӵ���������Ƭ���������з����α߻��Ե����дΪ-1��
		ttAdj_nmEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);

		// 4. �����αߣ�������������Ƭ������ߣ���Ϣ
		index = 0;
		std::vector<int> edgeIdx_nmn;
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (2 == adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)))
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();


		//			�ҳ�ͬһ������������߶�Ӧ�Ķ����������map<pair��ʾ�ķ����α����ݣ� ������vector>
		std::unordered_map<std::pair<int, int>, std::vector<int>, edgeHash, edgeComparator> edges_nmn_map;
		for (const auto& eIdx : edgeIdx_nmn)
		{
			std::pair<int, int> edge{ edges(eIdx, 0), edges(eIdx, 1) };
			auto retPair = edges_nmn_map.insert({ edge , std::vector<int>{eIdx} });
			if (!retPair.second)		// ������ʧ�ܣ���˵�����д˼���
			{
				auto iter = edges_nmn_map.find(edge);
				iter->second.push_back(eIdx);
			}
		}

		//			map<�����α�Idx, �ñߵĶԱ߶�Ӧ�Ķ����������vector>
		std::unordered_map<int, std::vector<int>> edgeIdx_nmn_map;
		for (const auto& eIdx : edgeIdx_nmn)
		{
			std::pair<int, int> edge{ edges(eIdx, 0), edges(eIdx, 1) };
			const std::vector<int>& commonEidxes = edges_nmn_map.find(edge)->second;
			edgeIdx_nmn_map.insert({ eIdx, commonEidxes });
		}
		std::unordered_map<int, std::vector<int>> edgeIdx_nmn_opp_map;
		for (const auto& pair : edgeIdx_nmn_map)
		{
			int eIdx = pair.first;
			int vaIdx = edges(eIdx, 0);
			int vbIdx = edges(eIdx, 1);
			int eOppIdx = -1;
			for (int i = 0; i < edgesCount; ++i)
				if (edges(i, 0) == vbIdx && edges(i, 1) == vaIdx)
				{
					eOppIdx = i;
					break;
				}
			auto iter = edgeIdx_nmn_map.find(eOppIdx);
			edgeIdx_nmn_opp_map.insert({eIdx, iter->second});
		}


		// 5. ���з����αߵ�����Ƭ���ڽӹ�ϵ��
		ttAdj_nmnEdge.resize(trisCount);
		for (const auto& pair : edgeIdx_nmn_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;

			std::vector<int> trisIdx_nmn = pair.second;
			for (auto& index : trisIdx_nmn)
				index = etInfo[index];

			switch (col)
			{
			case 0:	std::get<0>(ttAdj_nmnEdge[row]) = trisIdx_nmn;
				break;
			case 1:	std::get<1>(ttAdj_nmnEdge[row]) = trisIdx_nmn;
				break;
			case 2:	std::get<2>(ttAdj_nmnEdge[row]) = trisIdx_nmn;
				break;
			default:
				return false;
			}
		}

		ttAdj_nmnOppEdge.resize(trisCount);
		for (const auto& pair : edgeIdx_nmn_opp_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;

			std::vector<int> trisIdx_nmn_opp = pair.second;
			for (auto& index : trisIdx_nmn_opp)
				index = etInfo[index];

			switch (col)
			{
			case 0:	std::get<0>(ttAdj_nmnOppEdge[row]) = trisIdx_nmn_opp;
				break;
			case 1:	std::get<1>(ttAdj_nmnOppEdge[row]) = trisIdx_nmn_opp;
				break;
			case 2:	std::get<2>(ttAdj_nmnOppEdge[row]) = trisIdx_nmn_opp;
				break;
			default:
				return false;
			}
		}


		// for debug
		{
			Eigen::MatrixXd vers;
			Eigen::MatrixXi tris_useless;
			igl::readOBJ("E:/����/meshArranged_hole.obj", vers, tris_useless);

			for (int i = 0; i < ttAdj_nmnEdge.size(); ++i)
			{
				const auto& tuple = ttAdj_nmnEdge[i];
				const auto& tupleOpp = ttAdj_nmnOppEdge[i];
				if (!std::get<0>(tuple).empty() || !std::get<1>(tuple).empty() || !std::get<2>(tuple).empty())
				{
					const std::vector<int>& tmpVec = (std::get<0>(tuple).empty() && std::get<1>(tuple).empty()) ? std::get<2>(tuple)\
								: (std::get<0>(tuple).empty() ? std::get<1>(tuple): std::get<0>(tuple));
					Eigen::MatrixXi tris0 = tris.row(i);
					Eigen::MatrixXi tris00(tmpVec.size(), 3);
					for (int k = 0; k < tmpVec.size(); ++k)
						tris00.row(k) = tris.row(tmpVec[k]);
					igl::writeOBJ("E:/tris0.obj", vers, tris0);
					igl::writeOBJ("E:/tris00.obj", vers, tris00);

					const std::vector<int>& tmpVecOpp = (std::get<0>(tupleOpp).empty() && std::get<1>(tupleOpp).empty()) ? std::get<2>(tupleOpp)\
						: (std::get<0>(tupleOpp).empty() ? std::get<1>(tupleOpp) : std::get<0>(tupleOpp));
					Eigen::MatrixXi tris000(tmpVecOpp.size(), 3);
					for (int k = 0; k < tmpVecOpp.size(); ++k)
						tris000.row(k) = tris.row(tmpVecOpp[k]);
					igl::writeOBJ("E:/tris000.obj", vers, tris000);

					break;
				}
			}
		}
 
		return true;
	}


	// test buildAdjacency();
	void test2() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris, edges;
		igl::readOBJ("E:/����/meshArranged_hole.obj", vers, tris);

		Eigen::MatrixXi ttAdj_nmEdge;
		std::vector<ttTuple> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		bool retFlag = buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);
	}


	// boolean����select tris:
	enum BOOLEAN_TYPE
	{
		UNION = 0,
		DIFF = 1,
	};

	bool booleanSelectTris() 
	{



		return true;
	}



}
