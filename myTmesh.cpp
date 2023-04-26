#include "myTmesh.h"



// 生成T_MESH::di_cell对应的轴向包围盒网格；
void genAABBmesh(const T_MESH::di_cell& cell, Eigen::MatrixXd& vers, Eigen::MatrixXi& tris)
{
	vers.resize(0, 0);
	tris.resize(0, 0);
	T_MESH::Point arrow = cell.Mp - cell.mp;
	T_MESH::Point newOri = (cell.mp + cell.Mp) / 2.0;

	genCubeMesh(vers, tris);
	vers.col(0) *= arrow.x;
	vers.col(1) *= arrow.y;
	vers.col(2) *= arrow.z;
	vers.rowwise() += Eigen::RowVector3d{ newOri.x, newOri.y, newOri.z };
}


// 测试拓扑网格类tmesh
namespace TEST_TMESH
{
	// TMESH的IO，基本功能
	void test0()
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();												// ？？？This is mandatory
		T_MESH::Basic_TMesh tmesh1;
		T_MESH::Node* nodePtr = nullptr;
		tmesh1.load("E:/材料/tooth.obj");
		tmesh1.save("E:/meshInput.obj");

		const int versCount = tmesh1.V.numels();
		const int edgesCount = tmesh1.E.numels();
		const int trisCount = tmesh1.T.numels();

		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		TMesh2MeshMat(vers, tris, tmesh1);
		debugWriteMesh("meshInputCopy", vers, tris);

		// 顶点数据：
		nodePtr = tmesh1.V.head();
		nodePtr = nodePtr->next();					// 第二个顶点
		T_MESH::Vertex* vPtr = reinterpret_cast<T_MESH::Vertex*>(nodePtr->data);
		T_MESH::Vertex& ver0 = *vPtr;

		// 边数据；
		nodePtr = tmesh1.E.head();
		T_MESH::Edge* ePtr = reinterpret_cast<T_MESH::Edge*>(nodePtr->data);
		T_MESH::Edge& e0 = *ePtr;

		// 三角片数据：
		nodePtr = tmesh1.T.head();
		T_MESH::Triangle* tPtr = reinterpret_cast<T_MESH::Triangle*>(nodePtr->data);
		T_MESH::Triangle& t0 = *tPtr;

		// 表象转换：
		T_MESH::Basic_TMesh tMesh2;
		meshMat2tMesh(tMesh2, vers, tris);
		debugWriteMesh("tMesh2", tmesh1);

		std::cout << "finished." << std::endl;
	}


	// 测试T_MESH::di_cell类
	void test1()
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh mesh;
		T_MESH::Node* nodePtr = nullptr;
		mesh.load("E:/材料/tooth.obj");				// 必须是流形网格；
		debugWriteMesh("meshInput", mesh);

		T_MESH::di_cell cell(&mesh);			// 输入网格，生成di_cell对象；
		Eigen::MatrixXd aabbVers;
		Eigen::MatrixXi aabbTris;

		// 打印cell对应的包围盒；
		genAABBmesh(cell, aabbVers, aabbTris);
		debugWriteMesh("aabb", aabbVers, aabbTris);

		// fork()方法：
		T_MESH::di_cell cell2 = *cell.fork();
		genAABBmesh(cell, aabbVers, aabbTris);
		debugWriteMesh("aabb11", aabbVers, aabbTris);
		genAABBmesh(cell2, aabbVers, aabbTris);
		debugWriteMesh("aabb12", aabbVers, aabbTris);

		std::cout << "finished." << std::endl;
	}


	// 测试TMESH中的几何元素标记功能：
	void test2()
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/材料/tooth.obj");
		debugWriteMesh("meshInput", tMesh);

		// 1.遍历顶点，标记x坐标小于0的顶点；
		traverseVersList(tMesh.V, [&](T_MESH::Vertex* verPtr)
			{
				if (verPtr->x < 0)
					MARK_VISIT(verPtr);
			});

		// 2. 选中所有包含选中顶点的三角片；
		tMesh.growSelection();

		// 3. 删除所有选中三角片，及其包含的点和边；
		tMesh.removeSelectedTriangles();
		debugWriteMesh("meshOut", tMesh);

		std::cout << "finished." << std::endl;
	}


	// 拓扑性质、拓扑修复
	void test3()
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/材料/tooth.obj");
		debugWriteMesh("meshInput", tMesh);

		// 检测网格拓扑性质：
		const char* retStr = tMesh.checkConnectivity();				// 若有拓扑问题则会返回字符串；
		if (NULL != retStr)
			std::cout << retStr << std::endl;

		std::cout << "finished." << std::endl;
	}


	// 去除退化三角片；
	void test4()
	{
		int max_iters = 5;
		int inner_loops = 6;							// 每次大循环中去除退化三角片、去除自交的迭代次数；
		int holesCount = 0;
		bool flagIsct = false;
		bool flagDeg = false;

		T_MESH::TMesh::init();
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/材料/meshRepairInput.obj");
		debugWriteMesh("meshInput", tMesh);

		int removedCount = tMesh.removeSmallestComponents();						// d_boundaries, d_handles, d_shells赋值
		if (removedCount > 0)
			std::cout << "！！！输入网格有" << removedCount + 1 << "个单连通区域。" << std::endl;

		tMesh.deselectTriangles();
		tMesh.invertSelection();

		flagDeg = tMesh.strongDegeneracyRemoval(inner_loops);			// 全部清除成功返回true， 否则返回false

		debugWriteMesh("meshOut", tMesh);
		if (!flagDeg)
			debugDisp("！！！退化三角片没有全部清除成功。");

		debugDisp("finished.");
	}


	// 循环调用去除退化三角片，补洞
	void test44()
	{
		int max_iters = 5;
		int inner_loops = 6;							// 每次大循环中去除退化三角片、去除自交的迭代次数；
		int holesCount = 0;
		bool flagIsct = false;
		bool flagDeg = false;

		T_MESH::TMesh::init();
		T_MESH::Basic_TMesh tMesh;

		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/材料/holeMeshIsctFree2.obj");
		debugWriteMesh("meshInput", vers, tris);
		meshMat2tMesh(tMesh, vers, tris);
 
		int removedCount = tMesh.removeSmallestComponents();						// d_boundaries, d_handles, d_shells赋值
		if (removedCount > 1)
			std::cout << "！！！输入网格有" << removedCount << "个单连通区域。" << std::endl;

		tMesh.deselectTriangles();
		tMesh.invertSelection();

		int loopCount = 1;
		for (int i = 0; i < max_iters; ++i)
		{
			debugDisp("loopCount == ", loopCount++);
			holesCount = tMesh.boundaries();
			if (holesCount)
				tMesh.fillSmallBoundaries(0, true);
			flagDeg = tMesh.strongDegeneracyRemoval(inner_loops);			// 全部清除成功返回true， 否则返回false
			holesCount = tMesh.boundaries();
			if (flagDeg && 0 == holesCount)
				break;
		}

		debugWriteMesh("meshOut", tMesh);

		vers.resize(0, 0);
		tris.resize(0, 0);
		TMesh2MeshMat(vers, tris, tMesh);
		debugWriteMesh("meshMatOut", vers, tris);


		debugDisp("holesCount == ", holesCount);
		debugDisp("输出结果是否没有退化三角片： ", flagDeg);
		debugDisp("finished.");
	}

}