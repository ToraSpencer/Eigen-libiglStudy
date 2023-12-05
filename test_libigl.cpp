#include "test_libigl.h"
#define DATA_PATH "./data/"

static igl::opengl::glfw::Viewer viewer;				// libigl�еĻ���glfw����ʾ���ڣ�
static std::mutex g_mutex;


////////////////////////////////////////////////////////////////////////////////////////////// DEBUG �ӿ�
namespace MY_DEBUG
{
	static std::string g_debugPath = "E:/";

	// lambda������ӡstd::cout֧�ֵ����ͱ�����
	template <typename T>
	static auto disp = [](const T& arg)
	{
		std::cout << arg << ", ";
	};

	static void debugDisp()			// �ݹ���ֹ
	{						//		�ݹ���ֹ��Ϊ�޲λ���һ�����������ζ����ԡ�
		std::cout << std::endl;
		return;
	}


	template <typename T, typename... Types>
	static void debugDisp(const T& firstArg, const Types&... args)
	{
		std::cout << firstArg << " ";
		debugDisp(args...);
	}


	template <typename Derived>
	static void dispData(const Eigen::MatrixBase<Derived>& m)
	{
		auto dataPtr = m.data();
		unsigned elemsCount = m.size();

		for (unsigned i = 0; i < elemsCount; ++i)
			std::cout << dataPtr[i] << ", ";

		std::cout << std::endl;
	}


	template <typename Derived>
	static void dispElem(const Eigen::MatrixBase<Derived>& m)
	{
		const Derived& mm = m.derived();
		std::cout << mm(1, 1) << std::endl;
	}


	template<typename DerivedV>
	static void debugWriteVers(const char* name, const Eigen::MatrixBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat(path, vers);
	}


	template<typename DerivedV>
	static void debugWriteVers2D(const char* name, const Eigen::MatrixBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat2D(path, vers);
	}


	template<typename DerivedV>
	static void debugWriteMesh(const char* name, \
		const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteMeshMat(path, vers, tris);
	}


	template<typename DerivedV>
	static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, \
		const Eigen::MatrixBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteEdgesMat(path, edges, vers);
	}

}
using namespace MY_DEBUG;


////////////////////////////////////////////////////////////////////////////// libigl��������
namespace IGL_BASIC
{
	Eigen::MatrixXd vers, newVers, normals;
	Eigen::MatrixXi tris, tets;
	Eigen::SparseMatrix<double> L;

	// igl�л����ľ������ӣ�
	void test00() 
	{
		Eigen::MatrixXi m1(3, 3);
		Eigen::MatrixXi m11;
		Eigen::VectorXi vec1(2);

		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		dispMat(m1);

		// slice()����ʹ��������������������ȡ�����е�Ԫ�أ�
		/*
			  template <typename MatX,  typename DerivedR,  typename MatY>
			  IGL_INLINE void slice(
						const MatX& X,																	�������
						const Eigen::DenseBase<DerivedR> & R,							��������
						const int dim,															ά�ȡ���1-�У�2-�У�
						MatY& Y																				�������
					)
		*/

		vec1 << 0, 2;
		igl::slice(m1, vec1, 1, m11);				// ��ȡ��0, 2�У�
		dispMat(m11);

		igl::slice(m1, vec1, 2, m11);				// ��ȡ��0, 2�У�
		dispMat(m11);

		Eigen::VectorXi vec2(2);
		vec2 << 1, 2;
		igl::slice(m1, vec1, vec2, m11);		// �ȼ��ڣ�m11 == m1(0: 2, 1: 2);
		dispMat(m11);


		std::cout << "finished." << std::endl;
	}


	// igl�м�������������ԵĽӿڣ�
	void test000()
	{
		bool retFlag = igl::readOBJ("E:/����/tooth.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		Eigen::MatrixXd verNorms, triNorms, conerNorms;			// ���ַ���ÿ����һ����������
		Eigen::MatrixXi edges;									// ���
		Eigen::MatrixXi uEdges;								// ����ߣ�����undirected������unique�ģ�
		Eigen::VectorXi edgeUeInfo;						// �����������ж�Ӧ��������edgeUeInfo(i)������Ϊi�İ�߶�Ӧ������ߵ�������
		Eigen::MatrixXi UeTrisInfo;							// �����-��������Ƭӳ�������i�е�Ԫ�أ�����Ϊi������߹���������Ƭ��������
		Eigen::MatrixXi UeCornersInfo;					// �����-�Զ���ӳ�������i�е�Ԫ�أ�����Ϊi����������Ե������������������Ϊ�Ǳ�Ե���α����������Զ��㣻
		std::vector<std::vector<int>> UeEdgeInfo;
		Eigen::VectorXi uEC, uEE;
		Eigen::VectorXd dbArea;								// ������ÿ������Ƭ�����������
		std::vector<int> nbrVersIdx, nbrTrisIdx;

		// per_vertex_normals
		igl::per_vertex_normals(vers, tris, verNorms);
		igl::per_face_normals(vers, tris, triNorms);
		igl::per_corner_normals(vers, tris, conerNorms);

		// unique_edge_map()������������İ�ߡ�����ߡ������Ӧ��ϵ��
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, UeEdgeInfo);
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, uEC, uEE);

		// edge_flaps()������������߹���������Ƭ��
		igl::edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
		objWriteEdgesMat("E:/uEdge0.obj", Eigen::RowVector2i(uEdges.row(0)), vers);
		Eigen::RowVector3i tri00 = tris.row(UeTrisInfo(0, 0));
		Eigen::RowVector3i tri01 = tris.row(UeTrisInfo(0, 1));
		Eigen::MatrixXi tmpTris;
		matInsertRows(tmpTris, tri00);
		matInsertRows(tmpTris, tri01);
		objWriteMeshMat("E:/uEdge0relaTris.obj", vers, tmpTris);

		Eigen::MatrixXd oppVers(2, 3);
		oppVers.row(0) = vers.row(tri00(UeCornersInfo(0, 0)));
		oppVers.row(1) = vers.row(tri01(UeCornersInfo(0, 1)));
		objWriteVerticesMat("E:/oppVers.obj", oppVers);


		// circulation()����Ѱ������ߵĶ˵����������Ƭ������1���򶥵㣻
		igl::circulation(0, true, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx, nbrTrisIdx);

		Eigen::MatrixXd vers0;
		Eigen::MatrixXi tris0;
		objWriteVerticesMat("E:/edgeHead0.obj", Eigen::MatrixXd{ vers.row(uEdges(0, 0)) });
		subFromIdxVec(vers0, vers, nbrVersIdx);
		objWriteVerticesMat("E:/vers0.obj", vers0);
		subFromIdxVec(tris0, tris, nbrTrisIdx);
		igl::writeOBJ("E:/tris0.obj", vers, tris0);

		igl::circulation(0, false, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx, nbrTrisIdx);

		Eigen::MatrixXd vers1;
		Eigen::MatrixXi tris1;
		objWriteVerticesMat("E:/edgeTail0.obj", Eigen::MatrixXd{ vers.row(uEdges(0, 1)) });
		subFromIdxVec(vers1, vers, nbrVersIdx);
		objWriteVerticesMat("E:/vers1.obj", vers1);
		subFromIdxVec(tris1, tris, nbrTrisIdx);
		igl::writeOBJ("E:/tris1.obj", vers, tris1);


		// igl::doubleArea()����������������ÿ������Ƭ�����������
		igl::doublearea(vers, tris, dbArea);
 
		std::vector<double> dbAreaVec = eigenVec2Vec(dbArea);
		std::sort(dbAreaVec.begin(), dbAreaVec.end());
		std::cout << "finished." << std::endl;
	}


	// �ļ�IO
	void test0()
	{
		//// readOBJ(), writeObj()����OBJ�ļ���IO���ж�����أ�������·�������ơ�����Ƭ���������ݡ�
		igl::readOBJ("./data/bunny.obj", vers, tris);
		igl::writeOBJ("./data/bunny_export.obj", vers, tris);
		igl::writeOBJ("./data/bunnyVers.obj", vers, Eigen::MatrixXi{});			// ֻҪ���Ʋ���Ҫ����Ƭ�Ļ�������վ���
	
		vers.resize(0, 0);
		tris.resize(0, 0);
		igl::readOBJ("E:/����/jawMeshDense.obj", vers, tris);
		igl::writeOFF("E:/����/jawMeshDense.off", vers, tris);

		// ��ȡstl�ļ���
		std::string fileName{"E:/����/jawMeshSimplified"};
		vers.resize(0, 0);
		tris.resize(0, 0);
		std::ifstream fileIn((fileName + std::string{ ".stl" }).c_str(), std::ios::binary);			// stl�ļ��Ƕ������ļ���
		igl::readSTL(fileIn, vers, tris, normals);
		fileIn.close();
		igl::writeOBJ((fileName + std::string{".obj"}).c_str(), vers, tris);

		// ��.mesh�ļ�
		fileName = "E:/����/bunny";
		vers.resize(0, 0);
		tris.resize(0, 0);
		Eigen::MatrixXi tets;
		igl::readMESH((fileName + std::string{ ".mesh" }).c_str(), vers, tets, tris);
		igl::writeOBJ((fileName + std::string{ ".obj" }).c_str(), vers, tris);
		igl::writeOBJ((fileName + std::string{ "_tets.obj" }).c_str(), vers, tets);
		std::cout << "finished." << std::endl;
	}


	// libigl�е���ʾ������Viewer���������÷���
	void test1() 
	{
		igl::readMESH("E:/����/hand.mesh", vers, tets, tris);

		// Viewer::data()��������viewer�����ݶ�������ã�

		// ViewerData::set_mesh()�������붥�������Ƭ���������������ݣ�д�뵽ViewerData�ĳ�Ա�����У�

		// 1. ������װ�����ݣ�
		viewer.data().set_mesh(vers, tris);		

		 // 2. �趨������ת��Ĭ����������ת��set_rotation_type()��������ָ��������ת���
		viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);			 

		// 3. show_linesָ���Ƿ񻭳������ߣ�
		viewer.data().show_lines = 0;			

		// 4. ������Ⱦѭ����
		viewer.launch();
	}


	// libigl�е���ʾ������Viewer����������
	void test11() 
	{



		debugDisp("finished.");
	}
 

	// �����ݶ�
	void test2()
	{
		using namespace Eigen;
		using namespace std;

		igl::readOBJ("./data/rootTooth1.obj", vers, tris);

		Eigen::VectorXd U;
		igl::readDMAT("./data/cheburashka-scalar.dmat", U);		// һϵ�еĺ���ֵ

		Eigen::SparseMatrix<double> G;			// �ݶ�����
		igl::grad(vers, tris, G);

		// ���㺯��ֵU���ݶ�
		Eigen::MatrixXd GU = Map<const Eigen::MatrixXd>((G * U).eval().data(), tris.rows(), 3);

		// ����ֵU���ݶ�GU��ģ��
		const Eigen::VectorXd GU_mag = GU.rowwise().norm();

		viewer.data().set_mesh(vers, tris);
		viewer.data().set_data(U);

		// Average edge length divided by average gradient (for scaling)
		const double max_size = igl::avg_edge_length(vers, tris) / GU_mag.mean();

		// ÿ������Ƭ�����ϻ�һ��ָʾ�ߣ�����Ϊ�ݶȷ��� 
		Eigen::MatrixXd BC;
		igl::barycenter(vers, tris, BC);
		const Eigen::RowVector3d black(0, 0, 0);
		viewer.data().add_edges(BC, BC + max_size * GU, black);
		viewer.data().show_lines = false;	  // ����������

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


	// �����Χ���࣬�����Χ����
	void test6() 
	{
		// ���ɰ�Χ������
		Eigen::MatrixXd aabbVers, obbVers;
		Eigen::MatrixXi aabbTris, obbTris;
		Eigen::RowVector3d minp = Eigen::RowVector3d(-1, 2, -3);
		Eigen::RowVector3d maxp = Eigen::RowVector3d(4, 5, 6);
		Eigen::AlignedBox<double, 3> aabb(minp, maxp);
		genAABBmesh(aabbVers, aabbTris, aabb);
		igl::writeOBJ("E:/aabbMesh.obj", aabbVers, aabbTris);
 
		OBB<double> obb(aabb, Eigen::RowVector3d(3,4,5).normalized(), Eigen::RowVector3d(1, 2, 3) );
		genOBBmesh(obb, obbVers, obbTris);
		igl::writeOBJ("E:/obbMesh.obj", obbVers, obbTris);
 
		// contains()�����ж϶������Χ�еĹ�ϵ��������Ҫ����������ʾ���ڰ�Χ�б�����Ҳ��Ϊ���ڲ���
		std::cout << "aabb.contains (3,4,5) ? " << aabb.contains(Eigen::RowVector3d(3, 4, 5).transpose()) << std::endl;
		std::cout << "aabb.contains (4,5,5) ? " << aabb.contains(Eigen::RowVector3d(4, 5, 5).transpose()) << std::endl;
		std::cout << "aabb.contains (9,9,9) ? " << aabb.contains(Eigen::RowVector3d(9, 9, 9).transpose()) << std::endl;

		// OBB����д��contains()������
		std::cout << "obb.contains(0, 0, 0) ? " << obb.contains(Eigen::RowVector3d(0, 0, 0)) << std::endl;
		std::cout << "obb.contains(4, 2, 6) ? " << obb.contains(Eigen::RowVector3d(4, 2, 6)) << std::endl;
		std::cout << "obb.contains(-5, -5, -5) ? " << obb.contains(Eigen::RowVector3d(-5, -5, -5)) << std::endl;
 
		igl::writeOBJ("E:/426.obj", Eigen::MatrixXd{ Eigen::RowVector3d(4, 2, 6) }, Eigen::MatrixXi{});

		std::cout << "finished." << std::endl;
	}

}



/////////////////////////////////////////////////////////////////////////////// ͼ�㷨
namespace IGL_GRAPH 
{
	// ͼ���ݽṹ��ת����
	void test0() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("E:/����/cube�ظ�����Ƭ.obj", vers, tris);

		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		
		// ����Ƭ��Ϣ�õ�ͼ���ڽӱ��ڽӾ���
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);				// ����ߴ�����Ԫ��Ϊ1������Ϊ0��

		// �ֶ������ڽӾ���
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		unsigned edgesCount = 3 * trisCount;
		Eigen::SparseMatrix<int> adjSM_eCount;				// ������ڽӾ���Ȩ��Ϊ������ߵ�������

		Eigen::MatrixXi edges = Eigen::MatrixXi::Zero(edgesCount, 2);
		edges.block(0, 0, trisCount, 1) = tris.col(0);
		edges.block(0, 1, trisCount, 1) = tris.col(1);
		edges.block(trisCount, 0, trisCount, 1) = tris.col(1);
		edges.block(trisCount, 1, trisCount, 1) = tris.col(2);
		edges.block(2 * trisCount, 0, trisCount, 1) = tris.col(2);
		edges.block(2 * trisCount, 1, trisCount, 1) = tris.col(0);

		// FOR DEBUG:
		dispMatBlock(edges, 0, 10, 0, 1);
		dispMatBlock(edges, 3 * trisCount - 10, 3 * trisCount - 1, 0, 1);

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
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("E:/����/roundSurf.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);
 
		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);

		size_t startIdx = 165;			// �ӽ�Բ�ĵĶ��㣻
		Eigen::RowVector3d startVer = vers.row(startIdx);
		objWriteVerticesMat("E:/startVer.obj", startVer);

		// dfs:
		Eigen::VectorXi disCoveredIdx, bfsTreeVec, dfsTreeVec, closedIdx;
		igl::dfs(adjList, startIdx, disCoveredIdx, dfsTreeVec, closedIdx);
		objWriteTreePath("E:/dfsTree.obj", dfsTreeVec, vers);
		Eigen::VectorXi retCrev = closedIdx.reverse();
		auto cornerState = (retCrev == disCoveredIdx);

		std::vector<int> vec1, vec2, vec3;
		vec1 = eigenVec2Vec(disCoveredIdx);
		vec2 = eigenVec2Vec(closedIdx);
		vec3 = eigenVec2Vec(retCrev);

		// bfs:
		igl::bfs(adjList, startIdx, disCoveredIdx, bfsTreeVec);
		objWriteTreePath("E:/bfsTree.obj", bfsTreeVec, vers);


		// �Լ�д��DFS�㷨����Ŀǰ������
		std::vector<bool> visited(vers.rows(), false);
		std::list<int> discoveredVersIdx;						// �ѱ����ʵĶ��㣺
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


		myDfs(startIdx, startIdx);
		std::vector<int> tmpVec;
		tmpVec.insert(tmpVec.end(), closedVersIdx.begin(), closedVersIdx.end());
		Eigen::VectorXi myDfsTreeVec = vec2EigenVec(tmpVec);
		objWriteTreePath("E:/myDfsTree.obj", myDfsTreeVec, vers);

 
		std::cout << "finished." << std::endl;
	}


	// prime, dijkstra
	void test2() 
	{
		// dijkstra
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
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



/////////////////////////////////////////////////////////////////////////////// �ռ仮��
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



/////////////////////////////////////////////////////////////////////////////// IGLʵ�ֵĻ��������������㷨��
namespace IGL_BASIC_PMP 
{ 
	// ����������ͨ������ȡconnected_components()
	void test5() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/����/ԭ������/originalMesh.obj");
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();


		// 1. �����ڽӾ���
		Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
		tris2adjMat(adjSM_eCount, adjSM_eIdx, tris);

		Eigen::SparseMatrix<int> adjSM = adjSM_eCount;
		traverseSparseMatrix(adjSM, [&](auto& iter)
			{
				iter.valueRef() = 1;
			});

		// 2. ȷ������ͨ���򡪡�connected_components()
		Eigen::VectorXi connectedLabels, connectedCount;
		int conCount = igl::connected_components(adjSM, connectedLabels, connectedCount);

		// 3. ��ȡ���ĵ���ͨ�����еĶ��㣺 
		std::vector<int> retVec1 = eigenVec2Vec(connectedLabels);
		std::vector<int> retVec2 = eigenVec2Vec(connectedCount);

		int mainLabel = 0;						// ���������ͨ�������ı�ǩ��
		int mainLabelCount = 0;
		for (int i = 0; i < conCount; ++i)
		{
			if (connectedCount(i) > mainLabelCount)
			{
				mainLabel = i;
				mainLabelCount = connectedCount(i);
			}
		}

		std::unordered_set<int> indexSet;
		for (int i = 0; i < versCount; ++i)
		{
			if (mainLabel == connectedLabels(i))
				indexSet.insert(i);
		}

		std::vector<int> indexVec;
		indexVec.insert(indexVec.end(), indexSet.begin(), indexSet.end());
		Eigen::MatrixXd versOut;
		subFromIdxVec(versOut, vers, indexVec);

		std::vector<int> oldNewIdxInfo(versCount, -1);
		for (int i = 0; i < indexVec.size(); ++i)
		{
			int oldIdx = indexVec[i];
			oldNewIdxInfo[oldIdx] = i;
		}

		// 4. ��ȡ�����ͨ�����е�����Ƭ��
		Eigen::MatrixXi trisCopy = tris;
		int* intPtr = trisCopy.data();
		for (int i = 0; i < trisCopy.size(); ++i)
		{
			int oldIdx = *intPtr;
			*intPtr = oldNewIdxInfo[oldIdx];
			intPtr++;
		}

		Eigen::MatrixXi trisOut = Eigen::MatrixXi::Zero(trisCount, 3);
		int trisCountNew = 0;
		for (int i = 0; i < trisCount; ++i)
		{
			if (trisCopy(i, 0) >= 0 && trisCopy(i, 1) >= 0 && trisCopy(i, 2) >= 0)
			{
				trisOut.row(trisCountNew) = trisCopy.row(i);
				trisCountNew++;
			}
		}

		trisOut.conservativeResize(trisCountNew, 3);
		objWriteMeshMat("E:/mainConnectedMesh.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// �Լ�ʵ�ֵ�connected_conponents()
	void test55()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/����/ԭ������/originalMesh.obj");
		objWriteMeshMat("E:/meshInput.obj", vers, tris);

		bool retFlag = simplyConnectedLargest(versOut, trisOut, vers, tris);
		if (!retFlag)
			std::cout << "function failed!!!" << std::endl;

		objWriteMeshMat("E:/meshOut.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// ������������֮���hausdorff distance
	void test6() 
	{
		Eigen::MatrixXd vers0, vers00, vers1, vers11, vers2, vers22;
		Eigen::MatrixXi tris0, tris00, tris1, tris11, tris2, tris22;
		double hd0, hd00, hd1, hd2, hd3, hd4;
		hd00 = hd0 = hd1 = hd2 = hd3 = hd4 = 0;

		tiktok& tt = tiktok::getInstance();
		tt.start();
		objReadMeshMat(vers0, tris0, "E:/����/tooth.obj");
		objReadMeshMat(vers00, tris00, "E:/����/tooth1.obj");
		objReadMeshMat(vers1, tris1, "E:/����/rootTooth1.obj");
		objReadMeshMat(vers11, tris11, "E:/����/rootTooth2.obj"); 
		tt.endCout("meshes loading finished:");

		igl::hausdorff(vers0, tris0, vers00, tris00, hd0);
		igl::hausdorff(vers00, tris00, vers0, tris0, hd00);
		igl::hausdorff(vers1, tris1, vers11, tris11, hd1); 
		debugDisp("Hausdorff distance0 == ", hd0);
		debugDisp("Hausdorff distance00 == ", hd00);
		debugDisp("Hausdorff distance1 == ", hd1); 

		std::cout << "finished." << std::endl;
	}
 

	// LIBIGL������Ľӿ�
	void test7()
	{
		Eigen::MatrixXd vers, triNorms, barycenters, edgeArrows, vas, vbs;
		Eigen::MatrixXi tris, edges;
		objReadMeshMat(vers, tris, "E:/����/tooth.obj");

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
		viewer.data().show_lines = false;                 // ����������
		viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);    // �趨����������ת

		// ����Ƭ����ָʾ���ñ����ݵ���ʽ��Ⱦ������
		const Eigen::RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);    // RGBɫ��������
		double aveLen = edgesLen.mean();							// ���б߳���ƽ��ֵ��
		viewer.data().add_edges(barycenters - aveLen * triNorms, barycenters + aveLen * triNorms, red);         // ������ʷ����ú�ɫָʾ�߱�ʶ

		viewer.launch();

		std::cout << "finished." << std::endl;
	}
}
