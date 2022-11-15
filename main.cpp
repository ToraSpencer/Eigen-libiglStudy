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


int main()
{
	// DENSEMAT::test7();
	// SPARSEMAT::test1();

	
	// IGL_BASIC::test55();
	// IGL_DIF_GEO::test0();
	// IGL_GRAPH::test2();
	// IGL_SPACE_PARTITION::test0();
	IGL_BASIC_PMP::test3();


	// SCIENTIFICCALC::test7();
	// TEST_PMP::test3();
	// IGL_MATH::test1();

	// DECIMATION::test0();

	std::cout << "main() finished." << std::endl;
}
