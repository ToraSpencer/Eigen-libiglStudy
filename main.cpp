#include "dense_mat.h"
#include "sparse_mat.h"
#include "scientific_calc.h"
#include "igl_study.h"
#include <iomanip>
#include <winuser.h>
 
#include<stdio.h>
#include<assert.h>

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
	// igl::qslim()
	void test0() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		std::string fileName = "E:/材料/jawMeshDense.obj";
		igl::readOBJ(fileName, vers, tris);

		int trisCount = tris.rows();
		int tarTrisCount = std::round(0.6397 * trisCount);
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
		Eigen::VectorXi newOldVersInfo;
		igl::qslim(vers, tris, tarTrisCount, versOut, trisOut, newOldTrisInfo, newOldVersInfo);

		igl::writeOBJ("E:/meshSimplified_qslim.obj", versOut, trisOut);
 

		std::cout << "finished." << std::endl;
	}

	// 循环调用igl::qslim()精简一批网格：
	void test00()
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
		Eigen::MatrixXd vers, edgeArrows;
		Eigen::MatrixXi tris, edges;
		bool retFlag = true;
		objReadMeshMat(vers, tris, "E:/材料/jawMeshDense_qslim_meshFixed.obj");
		objWriteMeshMat("E:/meshInput.obj", vers, tris);

		const unsigned versCount = vers.rows();
		const unsigned trisCount = tris.rows();
 
		// 0. 计算边数据：
		getEdges(edges, vers, tris);
		getEdgeArrows(edgeArrows, edges, vers);
		Eigen::VectorXd edgesLen = edgeArrows.rowwise().norm();
		double minLen = edgesLen.minCoeff();
		std::cout << "minimum edge len is " << minLen << std::endl;

		// 1. 检测边缘有向边：
		Eigen::MatrixXi bdrys, bdryTris;
		std::vector<int> bdryTriIdxes;
		if (!bdryEdges(bdrys, bdryTriIdxes, tris))
			return;
		if (bdrys.rows() > 0)
		{
			std::cout << "bdrys.rows() == " << bdrys.rows() << std::endl;
			std::cout << "bdry data: " << std::endl;
			dispMat(bdrys);
			std::cout << std::endl;
			subFromIdxVec(bdryTris, tris, bdryTriIdxes);
			objWriteEdgesMat("E:/bdry.obj", bdrys, vers);
			objWriteMeshMat("E:/bdryTris.obj", vers, bdryTris);
		}

		// 2. 检测非流形有向边：
		Eigen::MatrixXi nmnEdges;
		if (!nonManifoldEdges(tris, nmnEdges))
			return;
		if(nmnEdges.rows() > 0)
			std::cout << "nmnEdges.rows() == " << nmnEdges.rows() << std::endl;
		objWriteEdgesMat("E:/nmnEdges.obj", nmnEdges, vers);


		// 3. 检测孤立顶点：
		std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, tris);
		if (!isoVerIdxes.empty())
		{
			std::cout << "isoVerIdxes.size() == " << isoVerIdxes.size() << std::endl;
			Eigen::MatrixXd isoVers;
			subFromIdxVec(isoVers, vers, isoVerIdxes);
			objWriteVerticesMat("E:/isoVers.obj", isoVers);
		}


		// 4. 检测退化边
		double degEdgeThreshold = 1e-3;
		Eigen::VectorXi degEdgeFlags = checkDegEdges(edges, edgeArrows, vers, tris, degEdgeThreshold);
		unsigned degEdgesCount = degEdgeFlags.sum();
		if (degEdgesCount > 0)
			std::cout << "degenerate edges count == " << degEdgesCount << std::endl;


		// 5. 检测退化三角片：
		Eigen::VectorXd trisAreaVec;
		std::vector<unsigned> degenTriIdxes;
		const double eps = 10e-9;
		if (!trisArea(trisAreaVec, vers, tris))
			return;
		for (unsigned i = 0; i < trisCount; ++i)
			if (trisAreaVec(i) < eps)
				degenTriIdxes.push_back(i);
		if (!degenTriIdxes.empty())
		{
			std::cout << "degenTriIdxes.size() == " << degenTriIdxes.size() << std::endl;

			Eigen::MatrixXi deTris;
			subFromIdxVec(deTris, tris, degenTriIdxes);
			std::cout << "degenerate tris data: " << std::endl;
			dispMat(deTris);
			std::cout << std::endl;
		}


		// 6. 提取所有单连通区域：
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

		igl::readOBJ("E:/meshNoDeTris.obj", vers, tris);

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

		//	2. 去除非法三角片：
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


	// 去除网格退化边（边长过短的边，将两个顶点融合）;
	void test11() 
	{
		Eigen::MatrixXd vers, versOut, edgeArrows;
		Eigen::MatrixXi tris, trisOut, edges;
		Eigen::VectorXi selectedIdxes, oldNewIdxInfo;
		igl::readOBJ("E:/材料/meshDegTris.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);
 
		// 1. 计算边数据：
		getEdges(edges, vers, tris);
		getEdgeArrows(edgeArrows, edges, vers);		
		Eigen::VectorXd edgesLen = edgeArrows.rowwise().norm();
		double minLen = edgesLen.minCoeff();
		std::cout << "minimum edge len is " << minLen << std::endl;

		// 2. 检测退化边：
		double degEdgeThreshold = 1e-3;			// 判定退化边的边长阈值；
		Eigen::VectorXi degEdgeFlags = checkDegEdges(edges, edgeArrows, vers, tris, degEdgeThreshold);
		int degEdgesCount = degEdgeFlags.sum();
		std::cout << "degEdgesCount == " << degEdgesCount << std::endl;

		// 3. 融合退化边——！！！若退化边阈值设置得太高，融合后的结果可能会有缺陷！！！
		int repVersCount = mergeDegEdges(versOut, trisOut, edges, edgeArrows, vers, tris, degEdgeFlags);
		std::cout << "repVersCount == " << repVersCount << std::endl;

		// 4. 融合之后可能有单联通区域分裂，取最大单连通区域：
		vers = versOut; 
		tris = trisOut;
		versOut.resize(0, 0);
		trisOut.resize(0, 0);
		simplyConnectedLargest(vers, tris, versOut, trisOut);
		std::cout << "remove isolated mesh: versCount == " << vers.rows() - versOut.rows() << ", trisCount == "\
			<< tris.rows() - trisOut.rows() << std::endl;

		// 5. 融合之后的检测：
		edges.resize(0, 0);
		edgeArrows.resize(0, 0);
		vers.resize(0, 0);
		tris.resize(0, 0);
		getEdges(edges, versOut, trisOut);
		getEdgeArrows(edgeArrows, edges, versOut);
		degEdgeFlags = checkDegEdges(edges, edgeArrows, versOut, trisOut, degEdgeThreshold);
		degEdgesCount = degEdgeFlags.sum();
		std::cout << "degenerate edges count == " << degEdgesCount << " after mergeDegEdges procedure." << std::endl;

		objWriteMeshMat("E:/meshOut.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// 检测、去除网格孤立顶点：
	void test2() 
	{
		Eigen::MatrixXd vers, versOut, edgeArrows;
		Eigen::MatrixXi tris, trisOut, edges;
		Eigen::VectorXi selectedIdxes, oldNewIdxInfo;
		igl::readOBJ("E:/材料/meshIsoVers.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, tris);
		if (isoVerIdxes.size() > 0)
			std::cout << isoVerIdxes.size() << " isolated vertices detected." << std::endl;

		removeIsoVers(versOut, trisOut, vers, tris, isoVerIdxes);
		igl::writeOBJ("E:/meshOut.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// 去除退化三角片：
	void test3() 
	{
		/*
			退化三角形的三种情形：
			A. 至少存在一条退化边，即三角形的两个或者三个顶点几乎重叠在一起；
					在融合duplicate vertices过程中，可以消除此种三角形；
			B. 不存在退化边，即三个顶点是接近共线的关系；		
		*/
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		bool retFlag = true;
		objReadMeshMat(vers, tris, "E:/材料/meshDegTris.obj");
		objWriteMeshMat("E:/meshInput.obj", vers, tris);

		const unsigned versCount = vers.rows();
		const unsigned trisCount = tris.rows();

		// 生成有向边编码-三角片字典：
		Eigen::MatrixXi edgeAs = Eigen::MatrixXi::Zero(trisCount, 2);
		Eigen::MatrixXi edgeBs = Eigen::MatrixXi::Zero(trisCount, 2);
		Eigen::MatrixXi edgeCs = Eigen::MatrixXi::Zero(trisCount, 2);
		Eigen::MatrixXi vaIdxes = tris.col(0);
		Eigen::MatrixXi vbIdxes = tris.col(1);
		Eigen::MatrixXi vcIdxes = tris.col(2);
		std::unordered_multimap<std::int64_t, unsigned> edgesMap;
		for (unsigned i = 0; i < trisCount; ++i)
		{
			std::int64_t codeA = encodeEdge(vbIdxes(i), vcIdxes(i));
			std::int64_t codeB = encodeEdge(vcIdxes(i), vaIdxes(i));
			std::int64_t codeC = encodeEdge(vaIdxes(i), vbIdxes(i));
			edgesMap.insert({ codeA, i });
			edgesMap.insert({ codeB, i });
			edgesMap.insert({ codeC, i });
		}


		// 检测退化三角片：
		Eigen::VectorXd trisAreaVec;
		std::vector<unsigned> degenTriIdxes;
		const double eps = 10e-9;								// 退化三角片面积阈值；
		if (!trisArea(trisAreaVec, vers, tris))
			return;
		for (unsigned i = 0; i < trisCount; ++i)
			if (trisAreaVec(i) < eps)
				degenTriIdxes.push_back(i);
		if (!degenTriIdxes.empty())
			std::cout << "degenTriIdxes.size() == " << degenTriIdxes.size() << std::endl;
		unsigned degCount = degenTriIdxes.size();
  
		// 打印退化三角片：
		Eigen::MatrixXi deTris;
		subFromIdxVec(deTris, tris, degenTriIdxes);
		objWriteMeshMat("E:/deTris.obj", vers, deTris);

		Eigen::MatrixXd vas, vbs, vcs, arrows1, arrows2, arrows3;
		Eigen::MatrixXd deTrisEdgeLen;
		vaIdxes = deTris.col(0);
		vbIdxes = deTris.col(1);
		vcIdxes = deTris.col(2);
		subFromIdxVec(vas, vers, vaIdxes);
		subFromIdxVec(vbs, vers, vbIdxes);
		subFromIdxVec(vcs, vers, vcIdxes);
		arrows1 = vbs - vas;
		arrows2 = vcs - vbs;
		arrows3 = vas - vcs;
		deTrisEdgeLen.resize(deTris.rows(), 3);
		deTrisEdgeLen.col(0) = arrows1.rowwise().norm();				// 边ab
		deTrisEdgeLen.col(1) = arrows2.rowwise().norm();				// 边bc
		deTrisEdgeLen.col(2) = arrows3.rowwise().norm();				// 边ca

		std::cout << "degenerate tris data: " << std::endl;
		dispMat(deTris);
		std::cout << std::endl;

		std::cout << "退化三角形的边长：" << std::endl;
		dispMat(deTrisEdgeLen);
		std::cout << std::endl;

		// 假设当前已不存在A型退化三角形，只需要处理B型退化三角形：

		// 1. 删除所有退化三角形：
		Eigen::MatrixXi trisCopy = tris;
		for (const auto& index : degenTriIdxes)
			trisCopy.row(index) = -Eigen::RowVector3i::Ones();				// 退化三角形标记为(-1, -1, -1)

		// 2. 提取退化三角形中最长的那条边关联的所有三角片：
		std::vector<std::int64_t> longEdgeOppCodes(degCount);
		std::vector<int> longEdgeOppVerIdx(degCount);
		for (unsigned i = 0; i < degCount; ++i)
		{
			int vaIdx0 = deTris(i, 0);
			int vbIdx0 = deTris(i, 1);
			int vcIdx0 = deTris(i, 2);
			if (deTrisEdgeLen(i, 1) >= deTrisEdgeLen(i, 0) && deTrisEdgeLen(i, 1) >= deTrisEdgeLen(i, 2))		// bc最长；
			{
				longEdgeOppCodes[i] = encodeEdge(vcIdx0, vbIdx0);
				longEdgeOppVerIdx[i] = vaIdx0;
			}
			else if (deTrisEdgeLen(i, 2) >= deTrisEdgeLen(i, 0) && deTrisEdgeLen(i, 2) >= deTrisEdgeLen(i, 1))		// ca边最长；
			{
				longEdgeOppCodes[i] = encodeEdge(vaIdx0, vcIdx0);
				longEdgeOppVerIdx[i] = vbIdx0;
			}
			else									      // ab边最长；
			{				
				longEdgeOppCodes[i] = encodeEdge(vbIdx0, vaIdx0);
				longEdgeOppVerIdx[i] = vcIdx0;
			}
		}

		// for debug
		std::vector<int> oppTriIdxes;


		// 3. 退化三角片ABC，若最长边是AB, 对边所在的三角片为BAX, 则BAX分解为BCX和CAX，
		for (unsigned i = 0; i < degCount; ++i) 
		{
			unsigned oppTriIdx = edgesMap.find(longEdgeOppCodes[i])->second;		

			// for debug
			oppTriIdxes.push_back(oppTriIdx);

			Eigen::RowVector3i oppTri = tris.row(oppTriIdx);
			std::pair<int, int> retPair = decodeEdge(longEdgeOppCodes[i]);
			int vbIdx = retPair.first;
			int vaIdx = retPair.second;
			int vcIdx = longEdgeOppVerIdx[i];
			int vxIdx = 0;
			for (unsigned k = 0; k < 3; ++k)
			{
				if (oppTri(k) != vaIdx && oppTri(k) != vbIdx)
				{
					vxIdx = oppTri(k);
					break;
				}
			}
			trisCopy.row(oppTriIdx) = Eigen::RowVector3i{vbIdx, vcIdx, vxIdx};
			matInsertRows<int , 3>(trisCopy, Eigen::RowVector3i{vcIdx, vaIdx, vxIdx});
		}

		// for debug;
		Eigen::MatrixXi oppTris;
		subFromIdxVec(oppTris, tris, oppTriIdxes);
		objWriteMeshMat("E:/tris2modify.obj", vers, oppTris);
 
		// 4. 删除被标记的退化三角片：
		tris.resize(trisCopy.rows(), 3);
		unsigned index = 0;
		for (unsigned i = 0; i < trisCopy.rows(); ++i) 
			if (trisCopy(i, 0) >= 0)
				tris.row(index++) = trisCopy.row(i);
		tris.conservativeResize(index, 3);

		objWriteMeshMat("E:/meshOut.obj", vers, tris);

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
	// DENSEMAT::test8();
	// SPARSEMAT::test0();

	// IGL_BASIC::test55();
	// IGL_DIF_GEO::test1();
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
 

	//TEST_DIP::test0();

	std::cout << "main() finished." << std::endl;
}
