#include "dense_mat.h"
#include "sparse_mat.h"
#include "scientific_calc.h"
#include "igl_study.h"
#include <iomanip>

#define DATA_PATH "./data/"

// 当前问题-easy
/*
	
*/


// 当前问题-hard
/*
	1. 给定流形网格上两个点，找出它们之间的最短路径。
	
	2. 使用OPENGL实现上面的问题，可以用鼠标点击确定网格上的任意两个点作为输入。 


*/



// 项目信息
/*
	编译环境： x64 Relase
 
	使用的第三方库:
		eigen
		libigl				 
		glad
		glfw
*/


// 项目中几个预定义宏：
/*
	IGL_STATIC_LIBRARY
	
	NOMINMAX
			解决windows.h和<algorithm>中同时定义std::max()，std::min()的冲突。		
	
	TUTORIAL_SHARED_PATH="G:/gitRepositories/ligIgl/libigl_CGAL_openGL/cmake/../external/../tutorial/data"
	
	CMAKE_INTDIR="Release"
*/


// 测试并行for循环PARALLEL_FOR()
void test0() 
{
	tiktok& tt = tiktok::getInstance();
	const unsigned colsCount = 500000;
	Eigen::MatrixXd m1{ Eigen::MatrixXd::Random(100, colsCount) };
	Eigen::VectorXd sumVec(colsCount);
	auto myPlus = [&](const unsigned colIdx)
	{
		Eigen::VectorXd colVec = m1.col(colIdx);
		sumVec(colIdx) = colVec.sum();
	};

	tt.start();
	for (unsigned i = 0; i < colsCount; ++i)
	{
		Eigen::VectorXd colVec = m1.col(i);
		sumVec(i) = colVec.sum();
	}
	tt.endCout("regular for-loop time comsumed: ");

	tt.start();
	PARALLEL_FOR(0, colsCount, myPlus);
	tt.endCout("PARALLEL_FOR loop time comsumed: ");

	std::cout << "finished." << std::endl;
}


template <typename T, int M, int N>
void dispData(const Eigen::Matrix<T, M, N>& m)
{
	auto dataPtr = m.data();
	unsigned elemsCount = m.size();

	for (unsigned i = 0; i < elemsCount; ++i)
		std::cout << dataPtr[i] << ", ";

	std::cout << std::endl;
}


template <typename Derived>
void dispData(const Eigen::PlainObjectBase<Derived>& m)
{
	int m0 = m.RowsAtCompileTime;
	int n0 = m.ColsAtCompileTime;

	auto dataPtr = m.data();
	unsigned elemsCount = m.size();

	for (unsigned i = 0; i < elemsCount; ++i)
		std::cout << dataPtr[i] << ", ";

	std::cout << std::endl;
}


template <typename Derived>
void dispElem(const Eigen::MatrixBase<Derived>& m)
{
	const Derived& mm = m.derived();
	std::cout << mm(1, 1) << std::endl;
}


// 网格精简：
namespace DECIMATION 
{
	// 循环调用igl::qslim()精简一批网格：
	void test0()
	{
		const unsigned meshesCount = 40;
		for (unsigned i = 0; i < meshesCount; ++i)
		{
			Eigen::MatrixXd vers, versOut;
			Eigen::MatrixXi tris, trisOut;
			char fileName[256];
			sprintf_s(fileName, 256, "E:/网格精简/splittedData4/splitedMesh%d.obj", i);
			igl::readOBJ(fileName, vers, tris);

			int trisCount = tris.rows();
			int tarTrisCount = std::round(0.6397 * trisCount);
			Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
			Eigen::VectorXi newOldVersInfo;
			igl::qslim(vers, tris, tarTrisCount, versOut, trisOut, newOldTrisInfo, newOldVersInfo);

			sprintf_s(fileName, 256, "E:/qslimOutput_%d.obj", i);
			igl::writeOBJ(fileName, versOut, trisOut);
			std::cout << "Loop " << i << " finished." << std::endl;
		}

		std::cout << "finished." << std::endl;
	}
} 


// 网格缺陷的检测和修复：
namespace MESH_REPAIR 
{
	// 检测网格边缘边、非流形边、孤立点、重复点等缺陷：
	void test0() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		bool retFlag = true;
		objReadMeshMat(vers, tris, "E:/材料/jawCoreNoDupVers.obj");

		const unsigned versCount = vers.rows();
		const unsigned trisCount = tris.rows();

		// 检测边缘有向边：
		Eigen::MatrixXi bdrys;
		if (!bdryEdges(bdrys, tris))
			return;

		if (bdrys.rows() > 0)
			std::cout << "bdrys.rows() == " << bdrys.rows() << std::endl;
		objWriteEdgesMat("E:/bdry.obj", bdrys, vers);

		// 检测非流形有向边：
		Eigen::MatrixXi nmnEdges;
		if (!nonManifoldEdges(tris, nmnEdges))
			return;
		if(nmnEdges.rows() > 0)
			std::cout << "nmnEdges.rows() == " << nmnEdges.rows() << std::endl;
		objWriteEdgesMat("E:/nmnEdges.obj", nmnEdges, vers);

		
		// 检测孤立顶点：
		Eigen::VectorXi verIdxVec = Eigen::VectorXi::LinSpaced(versCount, 0, versCount - 1);
		Eigen::MatrixXd isoVers;
		std::vector<unsigned> isoVerIdxes;
		int* dataPtr = tris.data();
		for (unsigned i = 0; i<tris.size(); ++i) 
		{
			verIdxVec(*dataPtr) = -1;
			dataPtr++;
		}
		for (unsigned i = 0; i < versCount; ++i)
			if (verIdxVec(i) >= 0)
				isoVerIdxes.push_back(i);
		if (!isoVerIdxes.empty())
		{
			std::cout << "isoVerIdxes.size() == " << isoVerIdxes.size() << std::endl;
			subFromIdxVec(isoVers, vers, isoVerIdxes);
			objWriteVerticesMat("E:/isoVers.obj", isoVers);
		}

		// 检测退化三角片：
		Eigen::VectorXd trisAreaVec;
		std::vector<unsigned> degenTriIdxes;
		const double eps = 10e-9;
		if (!trisArea(trisAreaVec, vers, tris))
			return;
		for (unsigned i = 0; i < trisCount; ++i)
			if (trisAreaVec(i) < eps)
				degenTriIdxes.push_back(i);
		if (!degenTriIdxes.empty())
			std::cout << "degenTriIdxes.size() == " << degenTriIdxes.size() << std::endl;

		// 提取所有单连通区域：
		Eigen::SparseMatrix<int> adjSM, adjSM_eCount, adjSM_eIdx;
		Eigen::VectorXi connectedLabels, connectedCount;
		adjMatrix(tris, adjSM_eCount, adjSM_eIdx);
		adjSM = adjSM_eCount;
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter) 
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});
		int scCount = simplyConnectedRegion(adjSM, connectedLabels, connectedCount);
		if (scCount < 0)
			return;
		if (scCount > 1)
			std::cout << "scCount == " << scCount << std::endl;
		
		std::cout << "finished." << std::endl;
	}


	// 去除网格重复顶点(duplicate vertices)，！！！当前有问题；
	void test1()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut, trisCopy;
		Eigen::VectorXi selectedIdxes, oldNewIdxInfo;

		igl::readOBJ("E:/meshInnerRev.obj", vers, tris);

		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		// 打印初始信息：
		Eigen::MatrixXi bdrys, bdryTris;
		std::vector<int> bdryTriIdxes;
		bdryEdges(bdrys, bdryTriIdxes, tris);
		subFromIdxVec(bdryTris, tris, bdryTriIdxes);
		igl::writeOBJ("E:/bdryTris.obj", vers, bdryTris);
		std::cout << "versCount == " << versCount << std::endl;
		std::cout << "trisCount == " << trisCount << std::endl;
		std::cout << "bdrysCount == " << bdrys.rows() << std::endl;
		std::cout << std::endl;

		// 1. 去除duplicated vertices——！！！当前有问题；
		igl::remove_duplicate_vertices(vers, 0, versOut, selectedIdxes, oldNewIdxInfo);
		objWriteVerticesMat("E:/versCleaned.obj", versOut);
		std::cout << "重复顶点数：" << versCount - versOut.rows() << std::endl;

		trisCopy = tris;
		int* ptr = trisCopy.data();
		for (unsigned i = 0; i < 3 * trisCount; ++i)
		{
			int oldIdx = *ptr;
			*ptr = oldNewIdxInfo(oldIdx);
			ptr++;
		}

		//		去除非法三角片：
		std::vector<unsigned> sickTriIdxes;
		checkSickTris(sickTriIdxes, trisCopy);
		trisOut = trisCopy;
		removeTris(trisOut, trisCopy, sickTriIdxes);

		bdrys.resize(0, 0);
		bdryTriIdxes.clear();
		bdryEdges(bdrys, bdryTriIdxes, trisOut);
		std::cout << "去除重复顶点后bdrysCount == " << bdrys.rows() << std::endl;
		std::cout << std::endl;
		igl::writeOBJ("E:/mesh去除重复顶点之后.obj", versOut, trisOut);

		trisCopy = trisOut;
		removeTris(trisOut, trisCopy, bdryTriIdxes);
		bdrys.resize(0, 0);
		bdryTriIdxes.clear();

		igl::writeOBJ("E:/meshOut.obj", versOut, trisOut);

		// 打印最终信息：
		versCount = versOut.rows();
		trisCount = trisOut.rows();
		std::cout << "final versCount == " << versCount << std::endl;
		std::cout << "final trisCount == " << trisCount << std::endl;
		std::cout << "final bdrysCount == " << bdrys.rows() << std::endl;

		std::cout << "finished." << std::endl;
	}
}


namespace TEMP_TEST
{
	// 处理STL原型网格——转换为OBJ文件，网格重心移动到原点，然后底面移动到和XOY平面平行
	void test1() 
	{
		Eigen::MatrixXd vers, normals;
		Eigen::MatrixXi tris;

		std::string fileName = "E:/胡凤棋668994L_13_22123111021";
		std::ifstream fileIn((fileName + std::string{ ".stl" }).c_str(), std::ios::binary);			// stl文件是二进制文件；
		igl::readSTL(fileIn, vers, tris, normals);
		Eigen::RowVector3d bary = vers.colwise().mean();
		vers = (vers.rowwise() - bary).eval();
		Eigen::VectorXd zValues = vers.col(2);
		double zmin = zValues.minCoeff();
		vers = (vers.rowwise() - Eigen::RowVector3d(0, 0, zmin)).eval();
		igl::writeOBJ(fileName + std::string{".obj"}, vers, tris);

		std::cout << "finished." << std::endl;
	}

}


int main()
{
	// DENSEMAT::test7();
	// SPARSEMAT::test0();

	// IGL_BASIC::test55();
	// IGL_DIF_GEO::test0();
	// IGL_GRAPH::test1();
	// IGL_SPACE_PARTITION::test0();
	// IGL_BASIC_PMP::test4();

	// SCIENTIFICCALC::test7();
	// TEST_PMP::test3();
	// IGL_MATH::test1();

	// DECIMATION::test0();

	// TEST_MYEIGEN::test3();

	// TEMP_TEST::test1();

	MESH_REPAIR::test0();



	std::cout << "main() finished." << std::endl;
}
