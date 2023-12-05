#include "test_libigl.h"
#define DATA_PATH "./data/"

static igl::opengl::glfw::Viewer viewer;				// libigl中的基于glfw的显示窗口；
static std::mutex g_mutex;


////////////////////////////////////////////////////////////////////////////////////////////// DEBUG 接口
namespace MY_DEBUG
{
	static std::string g_debugPath = "E:/";

	// lambda——打印std::cout支持的类型变量。
	template <typename T>
	static auto disp = [](const T& arg)
	{
		std::cout << arg << ", ";
	};

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


////////////////////////////////////////////////////////////////////////////// libigl基本功能
namespace IGL_BASIC
{
	Eigen::MatrixXd vers, newVers, normals;
	Eigen::MatrixXi tris, tets;
	Eigen::SparseMatrix<double> L;

	// igl中基础的矩阵算子：
	void test00() 
	{
		Eigen::MatrixXi m1(3, 3);
		Eigen::MatrixXi m11;
		Eigen::VectorXi vec1(2);

		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		dispMat(m1);

		// slice()——使用索引矩阵、索引向量提取矩阵中的元素；
		/*
			  template <typename MatX,  typename DerivedR,  typename MatY>
			  IGL_INLINE void slice(
						const MatX& X,																	输入矩阵
						const Eigen::DenseBase<DerivedR> & R,							索引向量
						const int dim,															维度——1-行，2-列；
						MatY& Y																				输出矩阵；
					)
		*/

		vec1 << 0, 2;
		igl::slice(m1, vec1, 1, m11);				// 提取第0, 2行；
		dispMat(m11);

		igl::slice(m1, vec1, 2, m11);				// 提取第0, 2列；
		dispMat(m11);

		Eigen::VectorXi vec2(2);
		vec2 << 1, 2;
		igl::slice(m1, vec1, vec2, m11);		// 等价于：m11 == m1(0: 2, 1: 2);
		dispMat(m11);


		std::cout << "finished." << std::endl;
	}


	// igl中计算基础几何属性的接口：
	void test000()
	{
		bool retFlag = igl::readOBJ("E:/材料/tooth.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		Eigen::MatrixXd verNorms, triNorms, conerNorms;			// 各种法向；每行是一个法向量；
		Eigen::MatrixXi edges;									// 半边
		Eigen::MatrixXi uEdges;								// 无向边，既是undirected的又是unique的；
		Eigen::VectorXi edgeUeInfo;						// 半边在无向边中对应的索引；edgeUeInfo(i)是索引为i的半边对应的无向边的索引；
		Eigen::MatrixXi UeTrisInfo;							// 无向边-关联三角片映射表——第i行的元素：索引为i的无向边关联的三角片的索引；
		Eigen::MatrixXi UeCornersInfo;					// 无向边-对顶点映射表——第i行的元素：索引为i的无向边所对的两个顶点的索引；若为非边缘流形边则有两个对顶点；
		std::vector<std::vector<int>> UeEdgeInfo;
		Eigen::VectorXi uEC, uEE;
		Eigen::VectorXd dbArea;								// 网格中每个三角片的两倍面积；
		std::vector<int> nbrVersIdx, nbrTrisIdx;

		// per_vertex_normals
		igl::per_vertex_normals(vers, tris, verNorms);
		igl::per_face_normals(vers, tris, triNorms);
		igl::per_corner_normals(vers, tris, conerNorms);

		// unique_edge_map()——计算网格的半边、无向边、及其对应关系；
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, UeEdgeInfo);
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, uEC, uEE);

		// edge_flaps()——计算无向边关联的三角片；
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


		// circulation()——寻找输入边的端点关联的三角片，及其1领域顶点；
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


		// igl::doubleArea()——计算输入网格每个三角片面积的两倍：
		igl::doublearea(vers, tris, dbArea);
 
		std::vector<double> dbAreaVec = eigenVec2Vec(dbArea);
		std::sort(dbAreaVec.begin(), dbAreaVec.end());
		std::cout << "finished." << std::endl;
	}


	// 文件IO
	void test0()
	{
		//// readOBJ(), writeObj()——OBJ文件的IO，有多个重载，至少有路径、点云、三角片这三个数据。
		igl::readOBJ("./data/bunny.obj", vers, tris);
		igl::writeOBJ("./data/bunny_export.obj", vers, tris);
		igl::writeOBJ("./data/bunnyVers.obj", vers, Eigen::MatrixXi{});			// 只要点云不需要三角片的话，传入空矩阵；
	
		vers.resize(0, 0);
		tris.resize(0, 0);
		igl::readOBJ("E:/材料/jawMeshDense.obj", vers, tris);
		igl::writeOFF("E:/材料/jawMeshDense.off", vers, tris);

		// 读取stl文件；
		std::string fileName{"E:/材料/jawMeshSimplified"};
		vers.resize(0, 0);
		tris.resize(0, 0);
		std::ifstream fileIn((fileName + std::string{ ".stl" }).c_str(), std::ios::binary);			// stl文件是二进制文件；
		igl::readSTL(fileIn, vers, tris, normals);
		fileIn.close();
		igl::writeOBJ((fileName + std::string{".obj"}).c_str(), vers, tris);

		// 读.mesh文件
		fileName = "E:/材料/bunny";
		vers.resize(0, 0);
		tris.resize(0, 0);
		Eigen::MatrixXi tets;
		igl::readMESH((fileName + std::string{ ".mesh" }).c_str(), vers, tets, tris);
		igl::writeOBJ((fileName + std::string{ ".obj" }).c_str(), vers, tris);
		igl::writeOBJ((fileName + std::string{ "_tets.obj" }).c_str(), vers, tets);
		std::cout << "finished." << std::endl;
	}


	// libigl中的显示窗口类Viewer——基本用法：
	void test1() 
	{
		igl::readMESH("E:/材料/hand.mesh", vers, tets, tris);

		// Viewer::data()——返回viewer中数据对象的引用；

		// ViewerData::set_mesh()——输入顶点和三角片矩阵生成网格数据，写入到ViewerData的成员变量中；

		// 1. 窗口中装载数据；
		viewer.data().set_mesh(vers, tris);		

		 // 2. 设定三轴旋转；默认是两轴旋转，set_rotation_type()方法可以指定其他旋转风格
		viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);			 

		// 3. show_lines指定是否画出网格线；
		viewer.data().show_lines = 0;			

		// 4. 进入渲染循环：
		viewer.launch();
	}


	// libigl中的显示窗口类Viewer——动画：
	void test11() 
	{



		debugDisp("finished.");
	}
 

	// 计算梯度
	void test2()
	{
		using namespace Eigen;
		using namespace std;

		igl::readOBJ("./data/rootTooth1.obj", vers, tris);

		Eigen::VectorXd U;
		igl::readDMAT("./data/cheburashka-scalar.dmat", U);		// 一系列的函数值

		Eigen::SparseMatrix<double> G;			// 梯度算子
		igl::grad(vers, tris, G);

		// 计算函数值U的梯度
		Eigen::MatrixXd GU = Map<const Eigen::MatrixXd>((G * U).eval().data(), tris.rows(), 3);

		// 函数值U的梯度GU的模长
		const Eigen::VectorXd GU_mag = GU.rowwise().norm();

		viewer.data().set_mesh(vers, tris);
		viewer.data().set_data(U);

		// Average edge length divided by average gradient (for scaling)
		const double max_size = igl::avg_edge_length(vers, tris) / GU_mag.mean();

		// 每个三角片重心上画一根指示线，方向为梯度方向。 
		Eigen::MatrixXd BC;
		igl::barycenter(vers, tris, BC);
		const Eigen::RowVector3d black(0, 0, 0);
		viewer.data().add_edges(BC, BC + max_size * GU, black);
		viewer.data().show_lines = false;	  // 隐藏网格线

		viewer.launch();
	}
 

	// 生成三维空间中的栅格：
	void test4() 
	{
		Eigen::MatrixXd vers, gridCenters;
		Eigen::MatrixXi tris;	
		Eigen::RowVector3i gridCounts;			// 三个维度上栅格的数目；
		unsigned num = 10;				// 跨度最大的那个维度(xyz中的一个)的栅格数；

		igl::readOBJ("E:/材料/tooth.obj", vers, tris);
		igl::voxel_grid(vers, 0, num, 1, gridCenters, gridCounts);

		igl::writeOBJ("E:/gridCenters.obj", gridCenters, Eigen::MatrixXi{});

		std::cout << "finished." << std::endl;
	}


	// 轴向包围盒类，方向包围盒类
	void test6() 
	{
		// 生成包围盒网格：
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
 
		// contains()——判断顶点与包围盒的关系，顶点需要用列向量表示；在包围盒表面上也视为在内部；
		std::cout << "aabb.contains (3,4,5) ? " << aabb.contains(Eigen::RowVector3d(3, 4, 5).transpose()) << std::endl;
		std::cout << "aabb.contains (4,5,5) ? " << aabb.contains(Eigen::RowVector3d(4, 5, 5).transpose()) << std::endl;
		std::cout << "aabb.contains (9,9,9) ? " << aabb.contains(Eigen::RowVector3d(9, 9, 9).transpose()) << std::endl;

		// OBB类重写的contains()方法：
		std::cout << "obb.contains(0, 0, 0) ? " << obb.contains(Eigen::RowVector3d(0, 0, 0)) << std::endl;
		std::cout << "obb.contains(4, 2, 6) ? " << obb.contains(Eigen::RowVector3d(4, 2, 6)) << std::endl;
		std::cout << "obb.contains(-5, -5, -5) ? " << obb.contains(Eigen::RowVector3d(-5, -5, -5)) << std::endl;
 
		igl::writeOBJ("E:/426.obj", Eigen::MatrixXd{ Eigen::RowVector3d(4, 2, 6) }, Eigen::MatrixXi{});

		std::cout << "finished." << std::endl;
	}

}



/////////////////////////////////////////////////////////////////////////////// 图算法
namespace IGL_GRAPH 
{
	// 图数据结构的转换：
	void test0() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("E:/材料/cube重复三角片.obj", vers, tris);

		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		
		// 由面片信息得到图的邻接表、邻接矩阵
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);				// 有向边存在则元素为1，否则为0；

		// 手动计算邻接矩阵：
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		unsigned edgesCount = 3 * trisCount;
		Eigen::SparseMatrix<int> adjSM_eCount;				// 有向边邻接矩阵，权重为该有向边的数量；

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

		std::cout << "索引为0的顶点的1邻域顶点：" << std::endl;
		traverseSTL(adjList[0], disp<int>);
		std::cout << std::endl << std::endl;

		std::cout << "索引为1的顶点的入边（即邻接矩阵中第1列中所有不为零的元素下标）：" << std::endl;
		for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM, 0); it; ++it)
		{
			std::cout << "value == " << it.value() << std::endl;
			std::cout << "row == " << it.row() << std::endl;			 // row index
			std::cout << "col == " << it.col() << std::endl;			 // col index (here it is equal to k)
			std::cout << std::endl;
		} 

		std::cout << "索引为2的顶点的出边（即邻接矩阵中第2行中所有不为零的元素下标）：" << std::endl;
		for (int i = 0; i < adjSM_eCount.outerSize(); ++i)			// 对列的遍历：
		{
			for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM_eCount, i); it; ++it)		// 当前列的列内迭代器；
			{
				if (2 == it.row())
				{
					std::cout << "value == " << it.value() << std::endl;
					std::cout << "row == " << it.row() << std::endl;			 // row index
					std::cout << "col == " << it.col() << std::endl;			 // col index (here it is equal to k)
				}
			}
		}

		// 判断是否有重复三角片——若有重复三角片，则存在重复有向边：
		bool hasDupTri = false;;
		for (int i = 0; i < adjSM_eCount.outerSize(); ++i)			// 对列的遍历：
		{
			if (hasDupTri)
				break;
			for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM_eCount, i); it; ++it)		// 当前列的列内迭代器；
			{
				if (it.value() > 1)
				{
					std::cout << "存在重复三角片" << std::endl;
					hasDupTri = true;
					break;
				}
			}
		}
		if (!hasDupTri)
			std::cout << "没有重复三角片" << std::endl;

		std::cout << "finished." << std::endl;
	}


	// 基本的图算法；
	void test1() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("E:/材料/roundSurf.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);
 
		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);

		size_t startIdx = 165;			// 接近圆心的顶点；
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


		// 自己写的DFS算法——目前有问题
		std::vector<bool> visited(vers.rows(), false);
		std::list<int> discoveredVersIdx;						// 已被访问的顶点：
		std::list<int> closedVersIdx;

		std::function<void(const int, const int)> myDfs = [&](const int index, const int parentIdx)
		{
			// 递归终止；
			if (visited[index])			
				return;

			visited[index] = true;
			discoveredVersIdx.push_back(index);
			const std::vector<int>& adjVersIdx = adjList[index];

			// 递归递推；
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
		igl::readOBJ("E:/材料/roundSurf.obj", vers, tris);

		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);

		Eigen::VectorXd min_distance;				// 图中所有顶点到指定顶点的最短路径长度；
		Eigen::VectorXi mst;						// 最小生成树；
		int verIdx0 = 165;						// 设起点为接近圆心的顶点；

		// 不指定重点时， igl::dijkstra()返回图中以verIdx0为起点的最小生成树；
		int retIdx = igl::dijkstra(vers, adjList, verIdx0, std::set<int>{}, min_distance, mst);
		std::cout << "retIdx == " << retIdx << std::endl;
		objWriteTreePath("E:/mst.obj", mst, vers);

		std::cout << "finished." << std::endl;
	}
 
}



/////////////////////////////////////////////////////////////////////////////// 空间划分
namespace IGL_SPACE_PARTITION
{
	// aabb树实现的BVH(层次包围体)
	void test0()
	{
		Eigen::MatrixXd vers, vers0, minDisVers;
		Eigen::MatrixXi tris;
		Eigen::VectorXd minSqrDis;
		Eigen::VectorXi minDisIdx;
		igl::readOBJ("E:/材料/tooth.obj", vers, tris);

		vers0.resize(2, 3);
		vers0.row(0) = vers.row(0);
		vers0.row(1) = vers.row(99);

		// igl::point_mesh_squared_distance()——使用BVH求顶点到网格最小距离；
		igl::point_mesh_squared_distance(vers0, vers, tris, minSqrDis, minDisIdx, minDisVers);
		igl::writeOBJ("E:/inputMesh.obj", vers, tris);
		igl::writeOBJ("E:/vers0.obj", vers0, Eigen::MatrixXi{});
		igl::writeOBJ("E:/minDisVers.obj", minDisVers, Eigen::MatrixXi{});

		//// 生成网格的BVH对象：
		//igl::AABB<double, 3> bvh;

		

		std::cout << "finished." << std::endl;
	}

}



/////////////////////////////////////////////////////////////////////////////// IGL实现的基础三角网格处理算法；
namespace IGL_BASIC_PMP 
{ 
	// 测试网格单连通区域提取connected_components()
	void test5() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/材料/原型数据/originalMesh.obj");
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();


		// 1. 生成邻接矩阵：
		Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
		tris2adjMat(adjSM_eCount, adjSM_eIdx, tris);

		Eigen::SparseMatrix<int> adjSM = adjSM_eCount;
		traverseSparseMatrix(adjSM, [&](auto& iter)
			{
				iter.valueRef() = 1;
			});

		// 2. 确定单连通区域——connected_components()
		Eigen::VectorXi connectedLabels, connectedCount;
		int conCount = igl::connected_components(adjSM, connectedLabels, connectedCount);

		// 3. 提取最大的单连通区域中的顶点： 
		std::vector<int> retVec1 = eigenVec2Vec(connectedLabels);
		std::vector<int> retVec2 = eigenVec2Vec(connectedCount);

		int mainLabel = 0;						// 包含最大联通顶点数的标签；
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

		// 4. 提取最大单连通区域中的三角片：
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


	// 自己实现的connected_conponents()
	void test55()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/材料/原型数据/originalMesh.obj");
		objWriteMeshMat("E:/meshInput.obj", vers, tris);

		bool retFlag = simplyConnectedLargest(versOut, trisOut, vers, tris);
		if (!retFlag)
			std::cout << "function failed!!!" << std::endl;

		objWriteMeshMat("E:/meshOut.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// 测量两个网格之间的hausdorff distance
	void test6() 
	{
		Eigen::MatrixXd vers0, vers00, vers1, vers11, vers2, vers22;
		Eigen::MatrixXi tris0, tris00, tris1, tris11, tris2, tris22;
		double hd0, hd00, hd1, hd2, hd3, hd4;
		hd00 = hd0 = hd1 = hd2 = hd3 = hd4 = 0;

		tiktok& tt = tiktok::getInstance();
		tt.start();
		objReadMeshMat(vers0, tris0, "E:/材料/tooth.obj");
		objReadMeshMat(vers00, tris00, "E:/材料/tooth1.obj");
		objReadMeshMat(vers1, tris1, "E:/材料/rootTooth1.obj");
		objReadMeshMat(vers11, tris11, "E:/材料/rootTooth2.obj"); 
		tt.endCout("meshes loading finished:");

		igl::hausdorff(vers0, tris0, vers00, tris00, hd0);
		igl::hausdorff(vers00, tris00, vers0, tris0, hd00);
		igl::hausdorff(vers1, tris1, vers11, tris11, hd1); 
		debugDisp("Hausdorff distance0 == ", hd0);
		debugDisp("Hausdorff distance00 == ", hd00);
		debugDisp("Hausdorff distance1 == ", hd1); 

		std::cout << "finished." << std::endl;
	}
 

	// LIBIGL中求法向的接口
	void test7()
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
}
