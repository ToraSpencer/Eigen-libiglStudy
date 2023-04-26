#include "myTmesh.h"



// ����T_MESH::di_cell��Ӧ�������Χ������
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


// ��������������tmesh
namespace TEST_TMESH
{
	// TMESH��IO����������
	void test0()
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();												// ������This is mandatory
		T_MESH::Basic_TMesh tmesh1;
		T_MESH::Node* nodePtr = nullptr;
		tmesh1.load("E:/����/tooth.obj");
		tmesh1.save("E:/meshInput.obj");

		const int versCount = tmesh1.V.numels();
		const int edgesCount = tmesh1.E.numels();
		const int trisCount = tmesh1.T.numels();

		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		TMesh2MeshMat(vers, tris, tmesh1);
		debugWriteMesh("meshInputCopy", vers, tris);

		// �������ݣ�
		nodePtr = tmesh1.V.head();
		nodePtr = nodePtr->next();					// �ڶ�������
		T_MESH::Vertex* vPtr = reinterpret_cast<T_MESH::Vertex*>(nodePtr->data);
		T_MESH::Vertex& ver0 = *vPtr;

		// �����ݣ�
		nodePtr = tmesh1.E.head();
		T_MESH::Edge* ePtr = reinterpret_cast<T_MESH::Edge*>(nodePtr->data);
		T_MESH::Edge& e0 = *ePtr;

		// ����Ƭ���ݣ�
		nodePtr = tmesh1.T.head();
		T_MESH::Triangle* tPtr = reinterpret_cast<T_MESH::Triangle*>(nodePtr->data);
		T_MESH::Triangle& t0 = *tPtr;

		// ����ת����
		T_MESH::Basic_TMesh tMesh2;
		meshMat2tMesh(tMesh2, vers, tris);
		debugWriteMesh("tMesh2", tmesh1);

		std::cout << "finished." << std::endl;
	}


	// ����T_MESH::di_cell��
	void test1()
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh mesh;
		T_MESH::Node* nodePtr = nullptr;
		mesh.load("E:/����/tooth.obj");				// ��������������
		debugWriteMesh("meshInput", mesh);

		T_MESH::di_cell cell(&mesh);			// ������������di_cell����
		Eigen::MatrixXd aabbVers;
		Eigen::MatrixXi aabbTris;

		// ��ӡcell��Ӧ�İ�Χ�У�
		genAABBmesh(cell, aabbVers, aabbTris);
		debugWriteMesh("aabb", aabbVers, aabbTris);

		// fork()������
		T_MESH::di_cell cell2 = *cell.fork();
		genAABBmesh(cell, aabbVers, aabbTris);
		debugWriteMesh("aabb11", aabbVers, aabbTris);
		genAABBmesh(cell2, aabbVers, aabbTris);
		debugWriteMesh("aabb12", aabbVers, aabbTris);

		std::cout << "finished." << std::endl;
	}


	// ����TMESH�еļ���Ԫ�ر�ǹ��ܣ�
	void test2()
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/����/tooth.obj");
		debugWriteMesh("meshInput", tMesh);

		// 1.�������㣬���x����С��0�Ķ��㣻
		traverseVersList(tMesh.V, [&](T_MESH::Vertex* verPtr)
			{
				if (verPtr->x < 0)
					MARK_VISIT(verPtr);
			});

		// 2. ѡ�����а���ѡ�ж��������Ƭ��
		tMesh.growSelection();

		// 3. ɾ������ѡ������Ƭ����������ĵ�ͱߣ�
		tMesh.removeSelectedTriangles();
		debugWriteMesh("meshOut", tMesh);

		std::cout << "finished." << std::endl;
	}


	// �������ʡ������޸�
	void test3()
	{
		g_debugPath = "E:/";

		T_MESH::TMesh::init();				// This is mandatory
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/����/tooth.obj");
		debugWriteMesh("meshInput", tMesh);

		// ��������������ʣ�
		const char* retStr = tMesh.checkConnectivity();				// ��������������᷵���ַ�����
		if (NULL != retStr)
			std::cout << retStr << std::endl;

		std::cout << "finished." << std::endl;
	}


	// ȥ���˻�����Ƭ��
	void test4()
	{
		int max_iters = 5;
		int inner_loops = 6;							// ÿ�δ�ѭ����ȥ���˻�����Ƭ��ȥ���Խ��ĵ���������
		int holesCount = 0;
		bool flagIsct = false;
		bool flagDeg = false;

		T_MESH::TMesh::init();
		T_MESH::Basic_TMesh tMesh;
		tMesh.load("E:/����/meshRepairInput.obj");
		debugWriteMesh("meshInput", tMesh);

		int removedCount = tMesh.removeSmallestComponents();						// d_boundaries, d_handles, d_shells��ֵ
		if (removedCount > 0)
			std::cout << "����������������" << removedCount + 1 << "������ͨ����" << std::endl;

		tMesh.deselectTriangles();
		tMesh.invertSelection();

		flagDeg = tMesh.strongDegeneracyRemoval(inner_loops);			// ȫ������ɹ�����true�� ���򷵻�false

		debugWriteMesh("meshOut", tMesh);
		if (!flagDeg)
			debugDisp("�������˻�����Ƭû��ȫ������ɹ���");

		debugDisp("finished.");
	}


	// ѭ������ȥ���˻�����Ƭ������
	void test44()
	{
		int max_iters = 5;
		int inner_loops = 6;							// ÿ�δ�ѭ����ȥ���˻�����Ƭ��ȥ���Խ��ĵ���������
		int holesCount = 0;
		bool flagIsct = false;
		bool flagDeg = false;

		T_MESH::TMesh::init();
		T_MESH::Basic_TMesh tMesh;

		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/����/holeMeshIsctFree2.obj");
		debugWriteMesh("meshInput", vers, tris);
		meshMat2tMesh(tMesh, vers, tris);
 
		int removedCount = tMesh.removeSmallestComponents();						// d_boundaries, d_handles, d_shells��ֵ
		if (removedCount > 1)
			std::cout << "����������������" << removedCount << "������ͨ����" << std::endl;

		tMesh.deselectTriangles();
		tMesh.invertSelection();

		int loopCount = 1;
		for (int i = 0; i < max_iters; ++i)
		{
			debugDisp("loopCount == ", loopCount++);
			holesCount = tMesh.boundaries();
			if (holesCount)
				tMesh.fillSmallBoundaries(0, true);
			flagDeg = tMesh.strongDegeneracyRemoval(inner_loops);			// ȫ������ɹ�����true�� ���򷵻�false
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
		debugDisp("�������Ƿ�û���˻�����Ƭ�� ", flagDeg);
		debugDisp("finished.");
	}

}