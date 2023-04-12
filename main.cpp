#include "dense_mat.h"
#include "sparse_mat.h"
#include "scientific_calc.h"
#include "igl_study.h"
#include "myTmesh.h"

#include<stdio.h>
#include<assert.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <winuser.h>
#include <string>

#include <windows.h>
#include <atlstr.h>			// 包含CString类。属于microsoft ATL(活动模板库avtive template library)
#include <atlconv.h>
#include <io.h>




#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <functional>
#include <algorithm>
#include <numeric>
#include <thread>
#include <windows.h>
#include <atlstr.h>		// 包含CString类。属于microsoft ATL(活动模板库avtive template library)
#include <io.h>
#include <SDKDDKVer.h>	
#include <windows.h>		
#include <iostream>
#include <locale>


#define DATA_PATH "./data/"

static igl::opengl::glfw::Viewer viewer;				// 全局的网格查看器对象；
static std::mutex g_mutex;						// 全局的互斥锁；


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


////////////////////////////////////////////////////////////////////////////// 基于WINAPI的一些接口：

// 读取某个目录下所有文件名、目录名；
void getFileNames(std::string path, std::vector<std::string>& files, bool blRecur = true)
{
	std::string str;
	struct _finddata_t fileinfo;			// 文件信息
	intptr_t hFile = _findfirst(str.assign(path).append("/*").c_str(), &fileinfo);							// 文件句柄	
	bool blFileValid = (hFile != -1);

	if (blFileValid)
	{
		do
		{
			bool isSubDir = (fileinfo.attrib & _A_SUBDIR);
			//如果是目录,递归查找；如果不是,把文件绝对路径存入vector中
			if (isSubDir & blRecur)
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					getFileNames(str.assign(path).append("/").append(fileinfo.name), files);
			}
			else
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					files.push_back(str.assign(path).append("/").append(fileinfo.name));
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}


// config文件路径字符串字面量 → std::wstring字符串；
static void GetConfigDir(std::wstring& confDir, const TCHAR* pszConfFile)
{
#ifdef WIN32
	std::wstring strRet(pszConfFile);
	if ((std::string::npos != strRet.find(L"./")) ||
		(std::string::npos != strRet.find(L".\\")))
	{
		TCHAR szCurDir[256] = { 0 };
		GetCurrentDirectory(256, szCurDir);
		confDir = szCurDir;
		confDir += strRet.substr(1, strRet.length());

		return;
	}
	confDir = strRet;
#endif
	return;
}


// 读取config文件中的某一个类型为float的键值
float INIGetFloat(const TCHAR* pszKey, const TCHAR* pszConfFile)
{
#ifdef WIN32
	std::wstring strFileName;
	GetConfigDir(strFileName, pszConfFile);
	TCHAR szValue[256] = { 0 };
	GetPrivateProfileString(L"ZHENGYADENTALCONFIG", pszKey, L"", szValue, 256, strFileName.c_str());

#ifdef UNICODE
	return _wtof(szValue);
#else
	return atof(szValue);
#endif
#else
	return 0.0f;
#endif

}


// 读取config文件中的某一个类型为int的键值
DWORD INIGetInt(const TCHAR* pszKey, const TCHAR* pszConfFile)
{
#ifdef WIN32
	std::wstring strFileName;
	GetConfigDir(strFileName, pszConfFile);
	return GetPrivateProfileInt(L"ZHENGYADENTALCONFIG", pszKey, 0, strFileName.c_str());
#else
	return 0;
#endif
}



////////////////////////////////////////////////////////////////////////////// TEST: 网格精简：
namespace DECIMATION 
{
	using Quadric = std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>;

	// q矩阵重载运算符；
	Quadric operator+(const Quadric& q1, const Quadric& q2)
	{
		Quadric result;
		std::get<0>(result) = (std::get<0>(q1) + std::get<0>(q2)).eval();
		std::get<1>(result) = (std::get<1>(q1) + std::get<1>(q2)).eval();
		std::get<2>(result) = (std::get<2>(q1) + std::get<2>(q2));
		return result;
	}


	// 判断边是否可以折叠：
	bool edge_collapse_is_valid(std::vector<int>& srcNbrIdx, std::vector<int>& desNbrIdx)
	{
		// Do we really need to check if edge is IGL_COLLAPSE_EDGE_NULL ?
		if (srcNbrIdx.size() < 2 || desNbrIdx.size() < 2)
		{
			// Bogus data
			assert(false);
			return false;
		}

		// determine if the first two vertices are the same before reordering.
		// If they are and there are 3 each, then (I claim) this is an edge on a single tet.
		const bool first_two_same = (srcNbrIdx[0] == desNbrIdx[0]) && (srcNbrIdx[1] == desNbrIdx[1]);
		if (srcNbrIdx.size() == 3 && desNbrIdx.size() == 3 && first_two_same)
			return false;           // single tet


		  // https://stackoverflow.com/a/19483741/148668
		std::sort(srcNbrIdx.begin(), srcNbrIdx.end());
		std::sort(desNbrIdx.begin(), desNbrIdx.end());
		std::vector<int> Nint;
		std::set_intersection(srcNbrIdx.begin(), srcNbrIdx.end(), desNbrIdx.begin(), desNbrIdx.end(), std::back_inserter(Nint));

		// check if edge collapse is valid: intersection of vertex neighbors of s and d should be exactly 2+(s,d) = 4

		// http://stackoverflow.com/a/27049418/148668
		if (Nint.size() != 2)
			return false;


		return true;
	}


	// 折叠索引为uEdgeIdx的无向边；
	bool collapseSingleEdge(const int uEdgeIdx, const Eigen::RowVectorXd& collapsedVer,
		std::vector<int>& nbrVersIdx_src,	const std::vector<int>& nbrTrisIdx_src,
		std::vector<int>& nbrVersIdx_des,	const std::vector<int>& nbrTrisIdx_des,
		Eigen::MatrixXd& vers, Eigen::MatrixXi& tris,	Eigen::MatrixXi& uEdges,
		Eigen::VectorXi& edgeUeInfo, Eigen::MatrixXi& UeTrisInfo,	Eigen::MatrixXi& UeCornersInfo,
		int& a_e1, int& a_e2,	int& a_f1, int& a_f2)
	{
		/*
		   Assign this to 0 rather than,  say,  -1 so that deleted elements will get draw as degenerate elements at vertex 0
				  (which should always exist and never get collapsed to anything else since it is the smallest index)
	   */
		const int eFlipFlag = uEdges(uEdgeIdx, 0) > uEdges(uEdgeIdx, 1);

		// lambda——某一条边的所有信息设置为无效信息；
		const auto& kill_edge = [&uEdges, &UeCornersInfo, &UeTrisInfo](const int uEdgeIdx)
		{
			uEdges(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			uEdges(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
			UeTrisInfo(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			UeTrisInfo(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
			UeCornersInfo(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			UeCornersInfo(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
		};

		// 定义无向边的起点和终点——索引小的为起点，索引大的为终点；
		const int srcIdx = eFlipFlag ? uEdges(uEdgeIdx, 1) : uEdges(uEdgeIdx, 0);
		const int desIdx = eFlipFlag ? uEdges(uEdgeIdx, 0) : uEdges(uEdgeIdx, 1);

		// 0. 判断该边是否可折叠
		if (!edge_collapse_is_valid(nbrVersIdx_src, nbrVersIdx_des))
			return false;

		// Important to grab neighbors of desIdx before monkeying with uEdges
		const std::vector<int>& nV2Fd = (!eFlipFlag ? nbrTrisIdx_src : nbrTrisIdx_des);

		assert(srcIdx < desIdx && "srcIdx should be less than desIdx");			 // The following implementation strongly relies on srcIdx<desIdx

		// 1. 边两端点都替换为collapsedVer
		vers.row(srcIdx) = collapsedVer;
		vers.row(desIdx) = collapsedVer;

		// update edge info for each flap
		const int trisCount = tris.rows();
		for (int side = 0; side < 2; side++)
		{
			const int f = UeTrisInfo(uEdgeIdx, side);
			const int v = UeCornersInfo(uEdgeIdx, side);
			const int sign = (eFlipFlag == 0 ? 1 : -1) * (1 - 2 * side);

			// next edge emanating from desIdx
			const int e1 = edgeUeInfo(f + trisCount * ((v + sign * 1 + 3) % 3));

			// prev edge pointing to srcIdx
			const int e2 = edgeUeInfo(f + trisCount * ((v + sign * 2 + 3) % 3));
			assert(uEdges(e1, 0) == desIdx || uEdges(e1, 1) == desIdx);
			assert(uEdges(e2, 0) == srcIdx || uEdges(e2, 1) == srcIdx);

			// face adjacent to f on e1,  also incident on desIdx
			const bool flip1 = UeTrisInfo(e1, 1) == f;
			const int f1 = flip1 ? UeTrisInfo(e1, 0) : UeTrisInfo(e1, 1);
			assert(f1 != f);
			assert(tris(f1, 0) == desIdx || tris(f1, 1) == desIdx || tris(f1, 2) == desIdx);

			// across from which vertex of f1 does e1 appear?
			const int vaIdx = flip1 ? UeCornersInfo(e1, 0) : UeCornersInfo(e1, 1);

			// Kill e1
			kill_edge(e1);

			// Kill f
			tris(f, 0) = IGL_COLLAPSE_EDGE_NULL;
			tris(f, 1) = IGL_COLLAPSE_EDGE_NULL;
			tris(f, 2) = IGL_COLLAPSE_EDGE_NULL;

			// map f1'srcIdx edge on e1 to e2
			assert(edgeUeInfo(f1 + trisCount * vaIdx) == e1);
			edgeUeInfo(f1 + trisCount * vaIdx) = e2;

			// side opposite f2,  the face adjacent to f on e2,  also incident on srcIdx
			const int opp2 = (UeTrisInfo(e2, 0) == f ? 0 : 1);
			assert(UeTrisInfo(e2, opp2) == f);
			UeTrisInfo(e2, opp2) = f1;
			UeCornersInfo(e2, opp2) = vaIdx;

			// remap e2 from desIdx to srcIdx
			uEdges(e2, 0) = uEdges(e2, 0) == desIdx ? srcIdx : uEdges(e2, 0);
			uEdges(e2, 1) = uEdges(e2, 1) == desIdx ? srcIdx : uEdges(e2, 1);
			if (side == 0)
			{
				a_e1 = e1;
				a_f1 = f;
			}
			else
			{
				a_e2 = e1;
				a_f2 = f;
			}
		}

		/*
			 finally,  reindex faces and uEdges incident on desIdx. Do this last so asserts make sense.
			 Could actually skip first and last,  since those are always the two collpased faces.
			 Nah,  this is handled by (tris(f, v) == desIdx)
			 Don't attempt to use Nde, Nse here because edgeUeInfo has changed
		*/
		{
			int p1 = -1;
			for (auto f : nV2Fd)
			{
				for (int v = 0; v < 3; v++)
				{
					if (tris(f, v) == desIdx)
					{
						const int e1 = edgeUeInfo(f + trisCount * ((v + 1) % 3));
						const int flip1 = (UeTrisInfo(e1, 0) == f) ? 1 : 0;
						assert(uEdges(e1, flip1) == desIdx || uEdges(e1, flip1) == srcIdx);
						uEdges(e1, flip1) = srcIdx;
						const int e2 = edgeUeInfo(f + trisCount * ((v + 2) % 3));

						// Skip if we just handled this edge (claim: this will be all except for the first non-trivial face)
						if (e2 != p1)
						{
							const int flip2 = (UeTrisInfo(e2, 0) == f) ? 0 : 1;
							assert(uEdges(e2, flip2) == desIdx || uEdges(e2, flip2) == srcIdx);
							uEdges(e2, flip2) = srcIdx;
						}

						tris(f, v) = srcIdx;
						p1 = e1;
						break;
					}
				}
			}
		}

		// Finally,  "remove" this edge and its information
		kill_edge(uEdgeIdx);
		return true;
	}


	// 打印优先队列队首的若干元素：
	void dispCosts(const igl::min_heap<std::tuple<double, int, int>>& pQueue, const unsigned num)
	{
		igl::min_heap<std::tuple<double, int, int>> pCopy = pQueue;
		std::vector<double> costValues;
		std::vector<std::vector<int>> edges;
		costValues.reserve(num);
		edges.reserve(num);

		for (unsigned i = 0; i < num; ++i)
		{
			auto& t = pCopy.top();
			costValues.push_back(std::get<0>(t));
			edges.push_back(std::vector<int>{std::get<1>(t), std::get<2>(t)});
			pCopy.pop();
		}

		for (unsigned i = 0; i < num; ++i)
			std::cout << costValues[i] << ", ";

		std::cout << std::endl << std::endl;
	}


	// igl::qslim()
	void test0()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		std::string filePath = "E:/材料/";
		std::string fileName = "gumMeshDense";
		igl::readOBJ(filePath + fileName + std::string{ ".obj" }, vers, tris);

		int trisCount = tris.rows();
		int tarTrisCount = 10000;						// 精简目标三角片 
		Eigen::VectorXi newOldTrisInfo;				// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
		Eigen::VectorXi newOldVersInfo;
		if (!igl::qslim(vers, tris, tarTrisCount, versOut, trisOut, newOldTrisInfo, newOldVersInfo))
			std::cout << "igl::qslim() failed!!!" << std::endl;

		std::stringstream ss;
		ss << "E:/" << fileName << "_qslim_" << tarTrisCount << ".obj";
		igl::writeOBJ(ss.str(), versOut, trisOut);


		std::cout << "finished." << std::endl;
	}


	// 循环调用igl::qslim()精简一批网格：
	void test00()
	{
		for (unsigned i = 3; i <= 7; ++i)
		{
			Eigen::MatrixXd vers, versOut;
			Eigen::MatrixXi tris, trisOut;
			char fileName[256];
			sprintf_s(fileName, 256, "E:/材料/jawMeshDense%d.obj", i);
			igl::readOBJ(fileName, vers, tris);

			int trisCount = tris.rows();
			int tarTrisCount = 150000;
			Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
			Eigen::VectorXi newOldVersInfo;
			igl::qslim(vers, tris, tarTrisCount, versOut, trisOut, newOldTrisInfo, newOldVersInfo);

			sprintf_s(fileName, 256, "E:/jawMeshDense%d_qslim.obj", i);
			igl::writeOBJ(fileName, versOut, trisOut);
			std::cout << "Loop " << i << " finished." << std::endl;
		}

		std::cout << "finished." << std::endl;
	}


	// 使用viewer对象的动画功能逐步查看精简过程中的边折叠：
	void test000()
	{
		using namespace igl;
		Eigen::MatrixXd vers;
	Eigen:Eigen::MatrixXi tris;
		igl::readOBJ("E:/材料/tmpMesh.obj", vers, tris);

		bool keepWorkingFlag = true;
		bool simplest_decimate_flag = false;
		bool qslim_decimate_flag = ~simplest_decimate_flag;

		// 1. 准备数据
		Eigen::MatrixXd versCopy = vers;
		Eigen::MatrixXi trisCopy = tris;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// 优先队列； 

		int vaIdx = -1;					// State variables keeping track of edge we just collapsed
		int vbIdx = -1;
		Eigen::VectorXi timeStamps;
		Eigen::MatrixXd collapsedVers;
		typedef std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double> Quadric;
		std::vector<Quadric> quadrics;				// 每个顶点的Q矩阵；

		edge_flaps(trisCopy, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());						// Could reserve with https://stackoverflow.com/a/29236236/148668        
		collapsedVers.resize(uEdges.rows(), versCopy.cols());				// If an edge were collapsed, we'd collapse it to these points:

		// 检测是否有非流形边：
		{
			Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> BF;
			Eigen::Array<bool, Eigen::Dynamic, 1> BE;
			if (!is_edge_manifold(trisCopy, uEdges.rows(), edgeUeInfo, BF, BE))
				return;
		}

		decimate_cost_and_placement_callback cost_and_placement;
		decimate_pre_collapse_callback pre_collapse;
		decimate_post_collapse_callback post_collapse;

		if (simplest_decimate_flag)
		{
			// 获取最简边折叠精简算法的函数子
			cost_and_placement = shortest_edge_and_midpoint;
			decimate_trivial_callbacks(pre_collapse, post_collapse);    // 生成最简边折叠中的pre_collapse和post_collapse函数子——什么都不做；
		}
		else
		{
			//  计算每个顶点的Q矩阵：
			per_vertex_point_to_plane_quadrics(vers, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, quadrics);

			// 获取qslim算法中的函数子：
			igl::qslim_optimal_collapse_edge_callbacks(uEdges, quadrics, vaIdx, vbIdx, cost_and_placement, pre_collapse, post_collapse);
		}

		decimate_stopping_condition_callback stopping_condition;

		int loopCount = 0;
		int num_collapsed;                                // 边收缩循环的计数；

		const auto& exportCurrentMesh = [&]()
		{
			Eigen::MatrixXi tris0(trisCopy.rows(), 3);
			Eigen::VectorXi _1;
			Eigen::VectorXi newOldTrisInfo;
			Eigen::VectorXi newOldVersInfo;
			newOldTrisInfo.resize(trisCopy.rows());
			int m = 0;
			for (int i = 0; i < trisCopy.rows(); i++)
			{
				if (trisCopy(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
					trisCopy(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
					trisCopy(i, 2) != IGL_COLLAPSE_EDGE_NULL)
				{
					tris0.row(m) = trisCopy.row(i);
					newOldTrisInfo(m) = i;
					m++;
				}
			}
			tris0.conservativeResize(m, tris0.cols());              // 这里相当于shrink_to_fit();
			newOldTrisInfo.conservativeResize(m);

			// 3.2. 删除网格中的孤立顶点：
			Eigen::MatrixXd versOut;
			Eigen::MatrixXi trisOut;
			igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);
			char str[256];
			sprintf_s(str, 256, "E:/decimating_output_%d.obj", loopCount);
			igl::writeOBJ(str, versOut, trisOut);
		};

		// lambda——所有数据复位
		const auto& reset = [&]()
		{
			tris = trisCopy;
			vers = versCopy;
			edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
			collapsedVers.resize(uEdges.rows(), vers.cols());
			Eigen::VectorXd costs(uEdges.rows());

			pQueue = {};            // https://stackoverflow.com/questions/2852140/priority-queue-clear-method
			timeStamps = Eigen::VectorXi::Zero(uEdges.rows());

			// r1. 计算每条边的折叠cost值，以及折叠之后的顶点坐标，存入优先队列pQueue；
			{
				Eigen::VectorXd costs(uEdges.rows());
				igl::parallel_for(uEdges.rows(), \
					[&](const int i)
					{
						double cost = i;
						Eigen::RowVectorXd edgeCenter(1, 3);           // 取边的中点作为边折叠之后的顶点；
						cost_and_placement(i, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, cost, edgeCenter);
						collapsedVers.row(i) = edgeCenter;
						costs(i) = cost;
					}, \
					10000);

				for (int i = 0; i < uEdges.rows(); i++)
					pQueue.emplace(costs(i), i, 0);
			}

			num_collapsed = 0;
			viewer.data().clear();
			viewer.data().set_mesh(vers, tris);
			viewer.data().set_face_based(true);
		};


		// lambda——执行网格精简，准备渲染的数据；
		const auto& pre_draw = [&](igl::opengl::glfw::Viewer& viewer)->bool
		{
			// p1. 每一次动画循环中，收缩1%的边；
			if (viewer.core().is_animating && !pQueue.empty())
			{
				bool FlagCollapsed = false;           // 本次循环中边收缩是否执行成功；

				// p1.1 执行边收缩——collapse edge
				const int max_iter = std::ceil(0.01 * pQueue.size());
				for (int j = 0; j < max_iter; j++)
				{
					if (!collapse_edge(cost_and_placement, vers, tris, uEdges, \
						edgeUeInfo, UeTrisInfo, UeCornersInfo, pQueue, timeStamps, collapsedVers))              // collapse_edge()重载2.1
						break;
					FlagCollapsed = true;
					num_collapsed++;
				}

				// p1.2 
				if (FlagCollapsed)
				{
					viewer.data().clear();
					viewer.data().set_mesh(vers, tris);
					viewer.data().set_face_based(true);
				}
			}

			// p2. 计算折叠后的网格三角片数（三个点索引都没有被标记为IGL_COLLAPSE_EDGE_NULL的三角片数）
			unsigned currentTrisCount = 0;
			for (int i = 0; i < tris.rows(); i++)
			{
				if (tris(i, 0) != IGL_COLLAPSE_EDGE_NULL || tris(i, 1) != IGL_COLLAPSE_EDGE_NULL || tris(i, 2) != IGL_COLLAPSE_EDGE_NULL)
					currentTrisCount++;
			}

			if (loopCount < 10)
				exportCurrentMesh();
			std::cout << "loop : " << loopCount++ << ", current trisCount == " << currentTrisCount << std::endl;
			return false;
		};


		// lambda——键盘事件响应；
		const auto& key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)->bool
		{
			switch (key)
			{
			case ' ':
				viewer.core().is_animating ^= 1;
				break;
			case 'R':
			case 'r':
				reset();
				break;
			case '0':
				keepWorkingFlag = false;
				viewer.core().is_animating = 0;
				viewer.launch_shut();
				viewer.shutdown_plugins();
				break;
			default:
				return false;
			}
			return true;
		};


		// 1. 初始化（复位）；
		reset();

		// 2. 打开窗口，动画循环
		viewer.core().background_color.setConstant(1);
		viewer.core().is_animating = true;                                 // 动画
		viewer.core().animation_max_fps = 1.0;						// 指定最大动画帧率；
		viewer.callback_key_down = key_down;
		viewer.callback_pre_draw = pre_draw;
		viewer.launch(keepWorkingFlag);

		// 3. 动画终止后，输出网格：
		Eigen::VectorXi newOldTrisInfo;
		Eigen::VectorXi newOldVersInfo;

		// 3.1. 删除所有含有标记为IGL_COLLAPSE_EDGE_NULL边的三角片：
		Eigen::MatrixXi tris0(trisCopy.rows(), 3);
		Eigen::VectorXi _1;
		newOldTrisInfo.resize(trisCopy.rows());
		int m = 0;
		for (int i = 0; i < trisCopy.rows(); i++)
		{
			if (trisCopy(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
				trisCopy(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
				trisCopy(i, 2) != IGL_COLLAPSE_EDGE_NULL)
			{
				tris0.row(m) = trisCopy.row(i);
				newOldTrisInfo(m) = i;
				m++;
			}
		}
		tris0.conservativeResize(m, tris0.cols());              // 这里相当于shrink_to_fit();
		newOldTrisInfo.conservativeResize(m);

		// 3.2. 删除网格中的孤立顶点：
		Eigen::MatrixXd versOut;
		Eigen::MatrixXi trisOut;
		igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);
		igl::writeOBJ("E:/decimating_output.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// 将网格处理成一个封闭网格——若存在边缘有向边，则将其和一个无限远点连接生成新三角片
	template <typename T>
	bool connect_boundary_to_infinity(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut,	Eigen::MatrixXi& trisOut,
			const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
	{
		// 1. 点云末尾加上一个无限远顶点：
		versOut.resize(vers.rows() + 1, vers.cols());
		versOut.topLeftCorner(vers.rows(), vers.cols()) = vers;
		auto inf = std::numeric_limits<T>::infinity();
		versOut.row(vers.rows()).setConstant(inf);

		// 2. 将所有边缘有向边和无限远顶点连接形成新的三角片；
		Eigen::MatrixXi bdrys, newTris;
		bdryEdges(bdrys, tris);
		if (0 == bdrys.size())
			return true;
		newTris.resize(bdrys.rows(), 3);
		newTris.col(0) = bdrys.col(1);
		newTris.col(1) = bdrys.col(0);
		newTris.col(2).array() = vers.rows();

		//// for debug:
		//auto newTrisVec = mat2triplets(bdrys);

		trisOut = tris;
		matInsertRows(trisOut, newTris);

		return true;
	}
 

	// 解析qslim()
	void test0000()
	{
		Eigen::MatrixXd vers, versOut, vers0, versOri;
		Eigen::MatrixXi tris, trisOut, tris0, trisOri;
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
		Eigen::VectorXi newOldVersInfo;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi edges, uEdges, UeTrisInfo, UeCornersInfo;
		Eigen::VectorXi timeStamps;
		Eigen::MatrixXd collapsedVers;
		Eigen::VectorXd costs;
		std::vector<Quadric> quadrics;													// 每个顶点的Q矩阵；
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// 优先队列；
		Eigen::VectorXi _1, I2;										// 最后删除孤立点和内补三角片时用到；
		int vaIdx = -1;											// State variables keeping track of edge we just collapsed
		int vbIdx = -1;
		bool clean_finish = false;
		tiktok& tt = tiktok::getInstance();

		// for debug:
		std::vector<triplet<double>> versMoni;
		std::vector<triplet<int>>	trisMoni;
		std::vector<doublet<int>> edgesMoni;
		std::vector<int>			vecMoni;
		int retIdx = -1;

		tt.start();
		igl::readOBJ("E:/材料/roundSurf.obj", versOri, trisOri);
		tt.endCout("elapsed time of loading mesh is: ");
		igl::writeOBJ("E:/meshIn.obj", versOri, trisOri);

		unsigned versCount = versOri.rows();
		unsigned trisCount = trisOri.rows();
		int trisCountNew = trisCount;
		int tarTrisCount = std::round(trisCount / 3);

		// 0. 检测是否有非流形有向边，有则直接退出；
		Eigen::MatrixXi nmnEdges;
		nonManifoldEdges(nmnEdges, trisOri);
		if (nmnEdges.size() > 0)
			return;

		// 1. 将网格处理成一个封闭网格——若存在边缘有向边，则将其和一个无限远点连接生成新三角片
		connect_boundary_to_infinity(vers, tris, versOri, trisOri);

		// 1.1 计算无向边信息：
		getEdges(edges, tris);
		getUedges(uEdges, edgeUeInfo, edges);
		getUeInfos(UeTrisInfo, UeCornersInfo, edgeUeInfo, uEdges, tris); 

		// 2. 计算每个顶点的Q矩阵：
		igl::per_vertex_point_to_plane_quadrics(vers, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, quadrics);

		// 3. 计算每条无向边的cost值，以此为优先级存入优先队列
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());
		collapsedVers.resize(uEdges.rows(), vers.cols());
		costs.resize(uEdges.rows());

		igl::parallel_for(uEdges.rows(), [&](const int ueIdx)
			{
				// 以下是cost_and_placement()的内容：
				double cost = ueIdx;
				Eigen::RowVectorXd p(1, 3);

				// Combined quadric
				Quadric quadric_p;
				quadric_p = quadrics[uEdges(ueIdx, 0)] + quadrics[uEdges(ueIdx, 1)];

				// Quadric: p'Ap + 2b'p + c,  optimal point: Ap = -b, or rather because we have row vectors: pA=-b
				const auto& A = std::get<0>(quadric_p);
				const auto& b = std::get<1>(quadric_p);
				const auto& c = std::get<2>(quadric_p);
				p = -b * A.inverse();
				cost = p.dot(p * A) + 2 * p.dot(b) + c;

				// Force infs and nans to infinity
				if (std::isinf(cost) || cost != cost)
				{
					cost = std::numeric_limits<double>::infinity();
					p.setConstant(0);
				}

				collapsedVers.row(ueIdx) = p;
				costs(ueIdx) = cost;
			},
			10000);

		for (int i = 0; i < uEdges.rows(); i++)
			pQueue.emplace(costs(i), i, 0);

		// 4. 边折叠的循环：
		int uEdgeIdx0;							// 优先队列首部的边；
		int e1, e2, f1, f2;
		while (true)
		{
			bool collapsed = true;
			std::tuple<double, int, int> edgeTuple;							// 优先队列首部的元素；
			std::vector<int> nbrTrisIdx_src, nbrVersIdx_src;
			std::vector<int>  nbrTrisIdx_des, nbrVersIdx_des;

			// w1. 取队首元素，队首元素出队：
			while (true)
			{
				// 1.1 若队列为空，退出循环；
				if (pQueue.empty())   // no uEdges to collapse
					assert("边折叠光了");

				// 1.2取队首元素，队首元素出队：
				edgeTuple = pQueue.top();
				if (std::get<0>(edgeTuple) == std::numeric_limits<double>::infinity())
					assert("队首边的cost是无穷大");

				pQueue.pop();
				uEdgeIdx0 = std::get<1>(edgeTuple);          // 队首的无向边索引;

				// 1.3 Check if matches timestamp
				if (std::get<2>(edgeTuple) == timeStamps(uEdgeIdx0))
					break;

				// 1.4 错误处理
				assert(std::get<2>(edgeTuple) < timeStamps(uEdgeIdx0) || timeStamps(uEdgeIdx0) == -1);          // must be stale or dead.
			}

			// w2. 计算当前边两端点1领域的顶点、三角片：
			igl::circulation(uEdgeIdx0, true, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_src, nbrTrisIdx_src);
			igl::circulation(uEdgeIdx0, false, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_des, nbrTrisIdx_des);

			//	w3. pre_collapse； 折叠队首的边——collapse_edge()重载1
			vaIdx = uEdges(uEdgeIdx0, 0);
			vbIdx = uEdges(uEdgeIdx0, 1);
			collapsed = collapseSingleEdge(uEdgeIdx0, collapsedVers.row(uEdgeIdx0), nbrVersIdx_src, nbrTrisIdx_src, \
				nbrVersIdx_des, nbrTrisIdx_des, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, e1, e2, f1, f2);

			//	w4. post_collapses:
			if (collapsed)
				quadrics[vaIdx < vbIdx ? vaIdx : vbIdx] = quadrics[vaIdx] + quadrics[vbIdx];

			// w5. 折叠操作之后，更新相关timeStamp， 更新相关边的cost值
			if (collapsed)
			{
				// 4.1 Erase the two,  other collapsed uEdges by marking their timestamps as -1
				timeStamps(e1) = -1;
				timeStamps(e2) = -1;

				// 4.2
				std::vector<int> nbrTrisIdx;
				nbrTrisIdx.reserve(nbrTrisIdx_src.size() + nbrTrisIdx_des.size());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_src.begin(), nbrTrisIdx_src.end());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_des.begin(), nbrTrisIdx_des.end());
				std::sort(nbrTrisIdx.begin(), nbrTrisIdx.end());
				nbrTrisIdx.erase(std::unique(nbrTrisIdx.begin(), nbrTrisIdx.end()), nbrTrisIdx.end());

				// 4.3 Collect all uEdges that must be updated
				std::vector<int> Ne;
				Ne.reserve(3 * nbrTrisIdx.size());
				for (auto& triIdx : nbrTrisIdx)
				{
					if (tris(triIdx, 0) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 1) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 2) != IGL_COLLAPSE_EDGE_NULL)
					{
						for (int i = 0; i < 3; i++)
						{
							const int ueIdx = edgeUeInfo(i * tris.rows() + triIdx);
							Ne.push_back(ueIdx);
						}
					}
				}

				// Only process edge once
				std::sort(Ne.begin(), Ne.end());
				Ne.erase(std::unique(Ne.begin(), Ne.end()), Ne.end());             // 去除重复元素：
				for (auto& ueIdx : Ne)
				{
					// 计算边折叠的cost值，及折叠后的顶点坐标：
					double cost;
					Eigen::RowVectorXd place;

					// 以下是cost_and_placement

					 // Combined quadric
					Quadric quadric_p;
					quadric_p = quadrics[uEdges(ueIdx, 0)] + quadrics[uEdges(ueIdx, 1)];

					// Quadric: place'Ap + 2b'place + c,  optimal point: Ap = -b, or rather because we have row vectors: pA=-b
					const auto& A = std::get<0>(quadric_p);
					const auto& b = std::get<1>(quadric_p);
					const auto& c = std::get<2>(quadric_p);
					place = -b * A.inverse();
					cost = place.dot(place * A) + 2 * place.dot(b) + c;

					// Force infs and nans to infinity
					if (std::isinf(cost) || cost != cost)
					{
						cost = std::numeric_limits<double>::infinity();						
						place.setConstant(0);
					}

					// Increment timestamp
					timeStamps(ueIdx)++;

					// Replace in queue
					pQueue.emplace(cost, ueIdx, timeStamps(ueIdx));
					collapsedVers.row(ueIdx) = place;
				}
			}
			else
				assert("edge collapse failed.");

			// w6. 判断是否满足终止条件stopping_condition，若满足则跳出while循环；
			if (collapsed)
			{
				if (f1 < trisCount)
					trisCountNew -= 1;
				if (f2 < trisCount)
					trisCountNew -= 1;

				bool stopConditionFlag = (trisCountNew <= (int)tarTrisCount);
				if (stopConditionFlag)
				{
					clean_finish = true;
					break;
				}
			}
			else				 // 边折叠失败，退出循环：
				assert("edge collapse failed.");
		}

		// 5. 删除所有含有标记为IGL_COLLAPSE_EDGE_NULL边的三角片：
		tris0.resize(tris.rows(), 3);
		newOldTrisInfo.resize(tris.rows());
		int index = 0;
		for (int i = 0; i < tris.rows(); i++)
		{
			if (tris(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
				tris(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
				tris(i, 2) != IGL_COLLAPSE_EDGE_NULL)
			{
				tris0.row(index) = tris.row(i);
				newOldTrisInfo(index) = i;
				index++;
			}
		}
		tris0.conservativeResize(index, tris0.cols());              // 这里相当于shrink_to_fit();
		newOldTrisInfo.conservativeResize(index);

		// 6. 删除网格中的孤立顶点：
		igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);

		//		for debug:
		igl::writeOBJ("E:/删除内部三角片前.obj", versOut, trisOut);

		// 7. ？？？删除内部三角片
		const Eigen::Array<bool, Eigen::Dynamic, 1> keep = (newOldTrisInfo.array() < trisCount);
		igl::slice_mask(Eigen::MatrixXi(trisOut), keep, 1, trisOut);
		igl::slice_mask(Eigen::VectorXi(newOldTrisInfo), keep, 1, newOldTrisInfo);
		igl::remove_unreferenced(Eigen::MatrixXd(versOut), Eigen::MatrixXi(trisOut), versOut, trisOut, _1, I2);
		igl::slice(Eigen::VectorXi(newOldVersInfo), I2, 1, newOldVersInfo);

		igl::writeOBJ("E:/qslimOutput.obj", versOut, trisOut);
		std::cout << "finished." << std::endl;
	}


	// 最简单的边折叠精简算法——使用边长来度量折叠损耗；
	void test1() 
	{
		Eigen::MatrixXd vers, versOut, vers0, versOri;
		Eigen::MatrixXi tris, trisOut, tris0, trisOri;
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
		Eigen::VectorXi newOldVersInfo;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi edges, uEdges, UeTrisInfo, UeCornersInfo;
		Eigen::VectorXi timeStamps;
		Eigen::MatrixXd collapsedVers;
		Eigen::VectorXd costs;
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// 优先队列；(无向边的cost, 无向边索引, 时间戳)
		Eigen::VectorXi _1, I2;										// 最后删除孤立点和内补三角片时用到；
		int vaIdx = -1;											// State variables keeping track of edge we just collapsed
		int vbIdx = -1;
		bool clean_finish = false;
		tiktok& tt = tiktok::getInstance();

		tt.start();
		igl::readOBJ("E:/材料/roundSurf.obj", versOri, trisOri);
		tt.endCout("elapsed time of loading mesh is: ");
		igl::writeOBJ("E:/meshIn.obj", versOri, trisOri);

		unsigned trisCount = trisOri.rows();
		int trisCountNew = trisCount;
		int tarTrisCount = std::round(trisCount/3);

		tt.start();

		// 0. 检测是否有非流形有向边，有则直接退出；
		Eigen::MatrixXi nmnEdges;
		nonManifoldEdges(nmnEdges, trisOri);
		if (nmnEdges.size() > 0)
			return;

		// 1. 将网格处理成一个封闭网格——若存在边缘有向边，则将其和一个无限远点连接生成新三角片
		connect_boundary_to_infinity(vers, tris, versOri, trisOri);

		// 1.1 计算无向边信息：
		getEdges(edges, tris);
		getUedges(uEdges, edgeUeInfo, edges);
		getUeInfos(UeTrisInfo, UeCornersInfo, edgeUeInfo, uEdges, tris);

		// 3. 计算每条无向边的cost值（边长），以此为优先级存入优先队列
		const unsigned ueCount = uEdges.rows();
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());
		collapsedVers.resize(uEdges.rows(), vers.cols());
		costs.resize(uEdges.rows());
		igl::parallel_for(ueCount, [&](const int ueIdx)
			{
				// 以下是cost_and_placement()的内容：
				int vaIdx = uEdges(ueIdx, 0);
				int vbIdx = uEdges(ueIdx, 1);
				double cost = (vers.row(vaIdx) - vers.row(vbIdx)).norm();
				Eigen::RowVector3d midPoint = (vers.row(vaIdx) + vers.row(vbIdx)) / 2;

				// Force infs and nans to infinity
				if (std::isinf(cost) || cost != cost)
				{
					cost = std::numeric_limits<double>::infinity();
					midPoint.setConstant(0);
				}

				collapsedVers.row(ueIdx) = midPoint;
				costs(ueIdx) = cost;
			},
			10000);
		for (int i = 0; i < ueCount; i++)
			pQueue.emplace(costs(i), i, 0);

		// 4. 边折叠的循环：
		int uEdgeIdx0;							// 优先队列首部的边；
		int e1, e2, f1, f2;
		while (true)
		{
			bool collapsed = true;
			std::tuple<double, int, int> edgeTuple;							// 优先队列首部的元素；
			std::vector<int> nbrTrisIdx_src, nbrVersIdx_src;				// vaIdx的1领域的三角片索引、顶点索引
			std::vector<int>  nbrTrisIdx_des, nbrVersIdx_des;			// vbIdx的1领域的三角片索引、顶点索引

			// w1. 取队首元素，队首元素出队：
			while (true)
			{
				// 1.1 若队列为空，退出循环；
				if (pQueue.empty())						 
					assert("边折叠光了");

				// 1.2取队首元素，队首元素出队：
				edgeTuple = pQueue.top();
				if (std::get<0>(edgeTuple) == std::numeric_limits<double>::infinity())
					assert("队首边的cost是无穷大");

				pQueue.pop();
				uEdgeIdx0 = std::get<1>(edgeTuple);          // 队首的无向边索引;

				// 1.3 Check if matches timestamp
				if (std::get<2>(edgeTuple) == timeStamps(uEdgeIdx0))
					break;
			}

			// w2. 计算当前边两端点1领域的顶点、三角片：
			vaIdx = uEdges(uEdgeIdx0, 0);
			vbIdx = uEdges(uEdgeIdx0, 1);
			igl::circulation(uEdgeIdx0, true, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_src, nbrTrisIdx_src);
			igl::circulation(uEdgeIdx0, false, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_des, nbrTrisIdx_des);
			//nbrVersIdx_des = oneRingVerIdxes(vaIdx, vers, tris);
			//nbrTrisIdx_des = oneRingTriIdxes(vaIdx, vers, tris);
			//nbrVersIdx_src = oneRingVerIdxes(vbIdx, vers, tris);
			//nbrTrisIdx_src = oneRingTriIdxes(vbIdx, vers, tris);

			//	w3. 折叠队首的边：
			collapsed = collapseSingleEdge(uEdgeIdx0, collapsedVers.row(uEdgeIdx0), nbrVersIdx_src, nbrTrisIdx_src, \
				nbrVersIdx_des, nbrTrisIdx_des, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, e1, e2, f1, f2);

			//	w4. post_collapses; 默认情形下什么都不做；

			// w5. 折叠操作之后，更新相关timeStamp， 更新相关边的cost值
			if (collapsed)
			{
				// w5.1. Erase the two,  other collapsed uEdges by marking their timestamps as -1
				timeStamps(e1) = -1;
				timeStamps(e2) = -1;

				// w5.2.
				std::vector<int> nbrTrisIdx;
				nbrTrisIdx.reserve(nbrTrisIdx_src.size() + nbrTrisIdx_des.size());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_src.begin(), nbrTrisIdx_src.end());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_des.begin(), nbrTrisIdx_des.end());
				std::sort(nbrTrisIdx.begin(), nbrTrisIdx.end());
				nbrTrisIdx.erase(std::unique(nbrTrisIdx.begin(), nbrTrisIdx.end()), nbrTrisIdx.end());

				// w5.3. Collect all uEdges that must be updated
				std::vector<int> Ne;
				Ne.reserve(3 * nbrTrisIdx.size());
				for (auto& triIdx : nbrTrisIdx)
				{
					if (tris(triIdx, 0) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 1) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 2) != IGL_COLLAPSE_EDGE_NULL)
					{
						for (int i = 0; i < 3; i++)
						{
							const int ueIdx = edgeUeInfo(i * tris.rows() + triIdx);
							Ne.push_back(ueIdx);
						}
					}
				}

				// w5.4. Only process edge once
				std::sort(Ne.begin(), Ne.end());
				Ne.erase(std::unique(Ne.begin(), Ne.end()), Ne.end());             // 去除重复元素：
				for (auto& ueIdx : Ne)			// 更新无向边cost值，折叠后位置
				{
					// 计算边折叠的cost值，及折叠后的顶点坐标：
					int vaIdx = uEdges(ueIdx, 0);
					int vbIdx = uEdges(ueIdx, 1);
					double cost = (vers.row(vaIdx) - vers.row(vbIdx)).norm();
					Eigen::RowVector3d midPoint = (vers.row(vaIdx) + vers.row(vbIdx)) / 2;

					// Force infs and nans to infinity
					if (std::isinf(cost) || cost != cost)
					{
						cost = std::numeric_limits<double>::infinity();
						midPoint.setConstant(0);
					}

					// Increment timestamp
					timeStamps(ueIdx)++;

					// Replace in queue
					pQueue.emplace(cost, ueIdx, timeStamps(ueIdx));
					collapsedVers.row(ueIdx) = midPoint;
				}
			}
			else
				assert("edge collapse failed.");

			// w6. 判断是否满足终止条件stopping_condition，若满足则跳出while循环；
			if (collapsed)
			{
				// 若stopping_condition函数子返回true，则满足终止条件，跳出折叠循环
				if (f1 < trisCount)
					trisCountNew -= 1;
				if (f2 < trisCount)
					trisCountNew -= 1;

				bool stopConditionFlag = (trisCountNew <= (int)tarTrisCount);
				if (stopConditionFlag)
				{
					clean_finish = true;
					break;
				}
			}
			else			 // 边折叠失败，退出循环：
				assert("edge collapse failed.");
		}

		// 5. 删除所有含有标记为IGL_COLLAPSE_EDGE_NULL边的三角片：
		tris0.resize(tris.rows(), 3);
		newOldTrisInfo.resize(tris.rows());
		int index = 0;
		for (int i = 0; i < tris.rows(); i++)
		{
			if (tris(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
				tris(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
				tris(i, 2) != IGL_COLLAPSE_EDGE_NULL)
			{
				tris0.row(index) = tris.row(i);
				newOldTrisInfo(index) = i;
				index++;
			}
		}
		tris0.conservativeResize(index, tris0.cols());              // 这里相当于shrink_to_fit();
		newOldTrisInfo.conservativeResize(index);

		// 6. 删除网格中的孤立顶点：
		igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);

		// 7. ？？？删除内部三角片
		const Eigen::Array<bool, Eigen::Dynamic, 1> keep = (newOldTrisInfo.array() < trisCount);
		igl::slice_mask(Eigen::MatrixXi(trisOut), keep, 1, trisOut);
		igl::slice_mask(Eigen::VectorXi(newOldTrisInfo), keep, 1, newOldTrisInfo);
		igl::remove_unreferenced(Eigen::MatrixXd(versOut), Eigen::MatrixXi(trisOut), versOut, trisOut, _1, I2);
		igl::slice(Eigen::VectorXi(newOldVersInfo), I2, 1, newOldVersInfo);

		tt.endCout("精简耗时：");


		igl::writeOBJ("E:/decimateOutput.obj", versOut, trisOut);
		std::cout << "finished." << std::endl;
	}


	// 求精简后的网格和原始网格的近似误差（来自QEM的paper）——速度太慢
	double calcSimpApproxError(const Eigen::MatrixXd& versSimp, const Eigen::MatrixXi& trisSimp, \
		const Eigen::MatrixXd& versOri, const Eigen::MatrixXi& trisOri)
	{
		// 注意不能出现太大的矩阵，否则容易bad alloc

		double appErr = 0;
		std::pair<unsigned, float> pair0;
		int versCount0 = versOri.rows();
		int versCount1 = versSimp.rows();
		const Eigen::MatrixXd& versMat0 = versOri;
		const Eigen::MatrixXd& versMat1 = versSimp;
		const Eigen::MatrixXi& trisMat0 = trisOri;
		const Eigen::MatrixXi& trisMat1 = trisSimp;
		Eigen::MatrixXd planeCoeff0, planeCoeff1;
		bool ret0 = trianglesPlane(planeCoeff0, versMat0, trisMat0);
		bool ret1 = trianglesPlane(planeCoeff1, versMat1, trisMat1);

		std::vector<long double> squaredDis0, squaredDis1;
		squaredDis0.reserve(versCount0);
		squaredDis1.reserve(versCount1);

		// 点到平面距离： dis == p.dot(v)，注意是有符号的距离； p == (a,b,c,d), v = (x0, y0, z0, 1), (a,b,c)是平面的归一化法向量；
		Eigen::MatrixXd versExt0{ Eigen::MatrixXd::Ones(versCount0, 4) };
		Eigen::MatrixXd versExt1{ Eigen::MatrixXd::Ones(versCount1, 4) };
		versExt0.leftCols(3) = versMat0.array().cast<double>();
		versExt1.leftCols(3) = versMat1.array().cast<double>();

		// for debug
		Eigen::Vector4d verColVec = versExt0.row(3).transpose();
		// 索引为verIdx的simp网格中的顶点，到ori网格中所有三角片平面的符号距离：
		Eigen::VectorXd signedDises = planeCoeff0 * verColVec;
		Eigen::VectorXd sqrDises = signedDises.array() * signedDises.array();

		// for new method:
		PARALLEL_FOR(0, versCount1, [&](int verIdx)
			{
				std::lock_guard<std::mutex> guard(g_mutex);
				Eigen::Vector4d verColVec = versExt1.row(verIdx).transpose();
				// 索引为verIdx的simp网格中的顶点，到ori网格中所有三角片平面的符号距离：
				Eigen::VectorXd signedDises = planeCoeff0 * verColVec;
				Eigen::VectorXd sqrDises = signedDises.array() * signedDises.array();
				squaredDis1.push_back(sqrDises.minCoeff());
			});

		PARALLEL_FOR(0, versCount0, [&](int verIdx)
			{
				std::lock_guard<std::mutex> guard(g_mutex);
				Eigen::Vector4d verColVec = versExt0.row(verIdx).transpose();
				// 索引为verIdx的simp网格中的顶点，到ori网格中所有三角片平面的符号距离：
				Eigen::VectorXd signedDises = planeCoeff1 * verColVec;
				Eigen::VectorXd sqrDises = signedDises.array() * signedDises.array();
				squaredDis0.push_back(sqrDises.minCoeff());
			});

		long double s0 = 0;
		long double s1 = 0;
		for (const auto& num : squaredDis0)
			s0 += num;
		for (const auto& num : squaredDis1)
			s1 += num;
		appErr = (s0 + s1) / (static_cast<long>(versCount0) + static_cast<long>(versCount1));

		return appErr;
	}


	// 测量精简网格的近似误差：
	void test2()
	{
		Eigen::MatrixXd vers1, vers2;
		Eigen::MatrixXi tris1, tris2;
		double appErr = 0;

		igl::readOBJ("E:/材料/jawMeshDense.obj", vers1, tris1);
		igl::readOBJ("E:/材料/jawMeshDense.obj", vers2, tris2);
		appErr = calcSimpApproxError(vers2, tris2, vers1, tris1);
		std::cout << "appErr == " << appErr << std::endl;

		vers2.resize(0, 0);
		tris2.resize(0, 0);
		igl::readOBJ("E:/材料/jawMeshDense_geoSimp_150000.obj", vers2, tris2);
		appErr = calcSimpApproxError(vers2, tris2, vers1, tris1);
		std::cout << "appErr == " << appErr << std::endl;

		vers2.resize(0, 0);
		tris2.resize(0, 0);
		igl::readOBJ("E:/材料/jawMeshDense_geoSimp_120000.obj", vers2, tris2);
		appErr = calcSimpApproxError(vers2, tris2, vers1, tris1);
		std::cout << "appErr == " << appErr << std::endl;

		vers2.resize(0, 0);
		tris2.resize(0, 0);
		igl::readOBJ("E:/材料/jawMeshDense_geoSimp_90000.obj", vers2, tris2);
		appErr = calcSimpApproxError(vers2, tris2, vers1, tris1);
		std::cout << "appErr == " << appErr << std::endl;

		vers2.resize(0, 0);
		tris2.resize(0, 0);
		igl::readOBJ("E:/材料/jawMeshDense_geoSimp_60000.obj", vers2, tris2);
		appErr = calcSimpApproxError(vers2, tris2, vers1, tris1);
		std::cout << "appErr == " << appErr << std::endl;

		std::cout << "finished." << std::endl;
	}

} 



////////////////////////////////////////////////////////////////////////////// TEST: 网格缺陷的检测和修复：
namespace MESH_REPAIR 
{
	// 检测网格边缘边、非流形边、孤立点、重复点、洞、重叠三角片等缺陷：
	void test0() 
	{
		Eigen::MatrixXd vers, edgeArrows, normals;
		Eigen::MatrixXi tris, edges;

		bool retFlag = true;
		//retFlag = igl::readSTL(fileIn, vers, tris, normals);							// 貌似igl::readSTL读取的网格顶点数据不对；

		objReadMeshMat(vers, tris, "E:/meshOut0.obj");
		objWriteMeshMat("E:/meshInput.obj", vers, tris);

		// 0. 去除重复三角片：
		int retNum = removeSickDupTris(vers, tris);
		if (retNum > 0)
			std::cout << "去除重复三角片：" << retNum << std::endl;

		const unsigned versCount = vers.rows();
		const unsigned trisCount = tris.rows();

		// 0.1. 计算边数据：
		Eigen::Index minEdgeIdx = 0;
		getEdges(edges, tris);
		getEdgeArrows(edgeArrows, edges, vers);
		Eigen::VectorXd edgesLen = edgeArrows.rowwise().norm();
		double minLen = edgesLen.minCoeff(&minEdgeIdx);
		std::cout << "minimum edge len is " << minLen << std::endl;

		// 1. 检测边缘有向边：
		Eigen::MatrixXi bdrys, bdryTris;
		std::vector<int> bdryTriIdxes;
		if (!bdryEdges(bdrys, bdryTriIdxes, tris))
			return;
		if (bdrys.rows() > 0)
		{
			std::cout << "bdrys.rows() == " << bdrys.rows() << std::endl;
			//std::cout << "bdry data: " << std::endl;
			//dispMat(bdrys);
			//std::cout << std::endl;
			subFromIdxVec(bdryTris, tris, bdryTriIdxes);
			objWriteEdgesMat("E:/bdry.obj", bdrys, vers);
			objWriteMeshMat("E:/bdryTris.obj", vers, bdryTris);

			// 1.1 若存在边缘，找洞：

		}

		// 2. 检测非流形有向边：
		Eigen::MatrixXi nmnEdges;
		std::vector<std::pair<int, std::pair<int, int>>> nmnEdgeInfos;
		if (nonManifoldEdges(nmnEdges, nmnEdgeInfos, tris) < 0)
			return;
		if (nmnEdges.rows() > 0)
			std::cout << "！！！存在非流形边，nmnEdges.rows() == " << nmnEdges.rows() << std::endl;
		objWriteEdgesMat("E:/nmnEdges.obj", nmnEdges, vers);

		// 3. 检测孤立顶点：
		std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, tris);
		if (!isoVerIdxes.empty())
		{
			std::cout << "isoVerIdxes.size() == " << isoVerIdxes.size() << std::endl;
			Eigen::MatrixXd isoVers;
			subFromIdxVec(isoVers, vers, isoVerIdxes);
			objWriteVerticesMat("E:/isoVers.obj", isoVers);
			Eigen::MatrixXd vers1;
			Eigen::MatrixXi tris1;
			removeIsoVers(vers1, tris1, vers, tris, isoVerIdxes);
			vers = vers1;
			tris = tris1;
			debugWriteMesh("meshNoIsoVers", vers, tris);
		}

		// 4. 检测重叠三角片：
		std::vector<std::pair<int, int>> opTrisPairs;
		int olCount = findOverLapTris(opTrisPairs, vers, tris);
		debugDisp("重叠三角片对数：", olCount);

		// 5. 检测退化边
		double degEdgeThreshold = 1e-3;
		Eigen::VectorXi degEdgeFlags = checkDegEdges(edges, edgeArrows, vers, tris, degEdgeThreshold);
		unsigned degEdgesCount = degEdgeFlags.sum();
		if (degEdgesCount > 0)
			std::cout << "degenerate edges count == " << degEdgesCount << std::endl;

		// 6. 检测退化三角片：
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

		// 7. 提取所有单连通区域：
		Eigen::SparseMatrix<int> adjSM, adjSM_eCount, adjSM_eIdx;
		Eigen::VectorXi connectedLabels, connectedCount, connectedTriLabels, connectedTriCount;
		Eigen::VectorXi connectedLabelsSCT, connectedCountSCT;
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
		unsigned nzCount = adjSM_eCount.nonZeros();
		adjSM = adjSM_eCount;
		traverseSparseMatrix(adjSM, [&adjSM](auto& iter) 
			{
				if (iter.value() > 0)
					iter.valueRef() = 1;
			});
		int scCount = simplyConnectedRegion(adjSM, connectedLabels, connectedCount);
		int sctCount = simplyTrisConnectedRegion(connectedLabelsSCT, connectedCountSCT, tris);
		if (scCount > 1 || sctCount > 1)
		{
			if (scCount != sctCount)
			{
				debugDisp("存在奇异点。");
				debugDisp("scCount == ", scCount);
			}

			if (sctCount > 1) 
			{
				std::vector<int> tmpVec = vec2Vec(connectedLabelsSCT);
				Eigen::VectorXi flagVec = connectedLabelsSCT;
				
				for (int i = 0; i < sctCount; ++i) 
				{
					flagVec.setZero();
					for (int k = 0; k < trisCount; ++k)
						if (i == connectedLabelsSCT[k])
							flagVec[k] = 1;
					Eigen::MatrixXi subTris;
					subFromFlagVec(subTris, tris, flagVec);
					debugWriteMesh((std::string{"meshSplitted"} + std::to_string(i) + std::string{".obj"}).c_str(), vers, subTris);
				}

			}
			debugDisp("sctCount == ", sctCount);
		}

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
		getEdges(edges, tris);
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
		simplyConnectedLargest(versOut, trisOut, vers, tris);
		std::cout << "remove isolated mesh: versCount == " << vers.rows() - versOut.rows() << ", trisCount == "\
			<< tris.rows() - trisOut.rows() << std::endl;

		// 5. 融合之后的检测：
		edges.resize(0, 0);
		edgeArrows.resize(0, 0);
		vers.resize(0, 0);
		tris.resize(0, 0);
		getEdges(edges, trisOut);
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


	// MeshFix——不能输入非流形网格，否则修复阶段可能会死循环；
	void test4()
	{
		T_MESH::TMesh::init();												// ？？？This is mandatory
		T_MESH::Basic_TMesh mesh;
		mesh.load("E:/材料/jawMeshDense_qslim_150000_noOpTris.obj");				// 网格载入时会计算壳体数n_shell;
		mesh.save("E:/meshFixInput.obj");

		tiktok& tt = tiktok::getInstance();
		unsigned versCount = mesh.V.numels();
		unsigned trisCount = mesh.T.numels();

		// 1. 提取最大单连通网格；
		int removedCount = mesh.removeSmallestComponents();					// d_boundaries, d_handles, d_shells赋值
		if (removedCount > 1)
			std::cout << "！！！输入网格有" << removedCount << "个单连通区域。" << std::endl;

		// 2. 补洞
		int holesCount = mesh.boundaries();				// ？？？
		int patchedCount = 0;
		if (holesCount)
		{
			std::cout << "！！！输入网格环形边界数+环柄数 == " << holesCount << "。" << std::endl;
			T_MESH::TMesh::warning("Patching holes\n");
			patchedCount = mesh.fillSmallBoundaries(0, true);
		}

		// 3. meshclean前计算环形边界数、环柄数；
		holesCount = mesh.boundaries();
		if (holesCount > 0)
			std::cout << "网格补洞后边界数+环柄数 == " << holesCount << "。" << std::endl;

		// 4. meshClean()——去除退化结构和三角片自交——默认max_iters == 10, inner_loops == 3；
		bool flagClean = false;
		int max_iters = 20;
		int inner_loops = 6;					// 每次大循环中去除退化三角片、去除自交的迭代次数；
		bool flagIsct, flagDeg;
		T_MESH::Triangle* t;
		T_MESH::Node* m;

		//		4.1. 
		mesh.deselectTriangles();
		mesh.invertSelection();

		//		4.2.修复流程的大循环 
		for (int i = 0; i < max_iters; i++)
		{
			//		f1. 去除退化三角片；
			flagDeg = mesh.strongDegeneracyRemoval(inner_loops);			// 全部清除成功返回true， 否则返回false

			//		f2. 
			mesh.deselectTriangles();
			mesh.invertSelection();

			//		f3. 去除自交三角片，补洞；
			flagIsct = mesh.strongIntersectionRemoval(inner_loops);			// 自交全部清除返回true，否则返回false;

			//		f4. 若前两项全部清除成功，进一步检查确认：
			if (flagIsct && flagDeg)
			{
				// 遍历三角片检测是否有退化；
				for (m = mesh.T.head(), t = (m) ? ((T_MESH::Triangle*)m->data) : NULL; m != NULL; m = m->next(), t = (m) ? ((T_MESH::Triangle*)m->data) : NULL)
					if (t->isExactlyDegenerate())
						flagIsct = false;

				// 若进一步检查没有问题，退出大循环；
				if (flagIsct)
				{
					flagClean = true;
					break;
				}
			}
		}

#ifdef LOCAL_DEBUG
		if (flagClean)
			std::cout << "meshclean() succeeded." << std::endl;
		else
			std::cout << "!!!meshclean() is not completed!!!" << std::endl;
#endif

		// 5. meshclean最后计算环形边界数、环柄数；
		holesCount = mesh.boundaries();
		if (holesCount > 0)
			std::cout << "输出网格边界数+环柄数 == " << holesCount << "。" << std::endl;

		// 6. 输出：
		mesh.save("E:/meshFixOutput.obj");


		std::cout << "finished." << std::endl;
	}


	// 修复边折叠精简之后的网格——寻找重叠三角片 
	void test5() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		bool retFlag = false;
		objReadMeshMat(vers, tris, "E:/材料/meshRepairInput_b_noDegTris.obj");
		objWriteMeshMat("E:/meshInput.obj", vers, tris);
		unsigned trisCount = tris.rows();

		// 1. 搜索重叠三角片：
		std::vector<std::pair<int, int>> opTrisPairs;
		std::vector<int> opTriIdxes;
		Eigen::MatrixXi opTris;
		int olCount = findOverLapTris(opTrisPairs, vers, tris);

		if (olCount > 0)
		{
			debugDisp("重叠三角片对数：", olCount);
			for (const auto& pair : opTrisPairs)
			{
				opTriIdxes.push_back(pair.first);
				opTriIdxes.push_back(pair.second);
			}
			subFromIdxVec(opTris, tris, opTriIdxes);
			debugWriteMesh("opTris", vers, opTris);

			// 2. 处理重叠三角片——删掉面积小的，保留面积大的：
			Eigen::VectorXi triFlags{ Eigen::VectorXi::Ones(trisCount) };
			for (const auto& pair : opTrisPairs)
			{
				Eigen::MatrixXi opTris(2, 3);
				Eigen::VectorXd areas;
				opTris.row(0) = tris.row(pair.first);
				opTris.row(1) = tris.row(pair.second);
				trisArea(areas, vers, opTris);
				int smallIdx = areas(0) < areas(1) ? pair.first : pair.second;
				triFlags(smallIdx) = 0;
			}
			subFromFlagVec(trisOut, tris, triFlags);
			tris = trisOut;
			debugWriteMesh("noOpTris", vers, tris);
		}
		else
			debugDisp("输入网格没有重叠三角片。");

		// 2. 打印结果网格中的洞、边缘曲线：


		debugDisp("finished.");
	}
}



////////////////////////////////////////////////////////////////////////////// TEST: 暂时无法分类的测试：
namespace TEMP_TEST
{
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


////////////////////////////////////////////////////////////////////////////// 生成控制台程序工具：

// 批量读取inputOBJ文件夹中的网格，调用MESHFIX修复：
int testCmd_meshFix(int argc, char** argv)
{
	tiktok& tt = tiktok::getInstance();
	float deciRatio = 0.5;							// 精简率

	// 生成路径：
	int   nPos;
	CString   cPath;
	GetModuleFileName(NULL, cPath.GetBufferSetLength(MAX_PATH + 1), MAX_PATH);		// 获取当前进程加载的模块的路径。
	nPos = cPath.ReverseFind('\\');
	cPath = cPath.Left(nPos);
	std::string path{ CT2CA{cPath} };
	std::string pathOBJ = path + "\\inputOBJ";
	std::string pathOutput = path + "\\outputData";

	// 1. 读取部件网格 
	tt.start();
	std::cout << "读取输入网格..." << std::endl;
	std::vector<std::string> fileNames, tmpStrVec, OBJfileNames;
	std::vector<Eigen::MatrixXd> meshesVers, outVers;
	std::vector<Eigen::MatrixXi> meshesTris, outTris;

	getFileNames(pathOBJ.c_str(), tmpStrVec, false);
	const unsigned meshesCount = tmpStrVec.size();
	OBJfileNames.reserve(meshesCount);
	for (const auto& str : tmpStrVec)
	{
		std::string tailStr = str.substr(str.size() - 4, 4);				//	".obj"
		if (".obj" == tailStr)
		{
			fileNames.push_back(str);
			unsigned index = str.find_last_of("/");
			std::string OBJfileName = str.substr(index, str.size() - index - 4);			// "/" + obj文件名，不含路径和.obj后缀；
			OBJfileNames.push_back(OBJfileName);
		}
	}

	meshesVers.resize(meshesCount);
	meshesTris.resize(meshesCount);
	outVers.resize(meshesCount);
	outTris.resize(meshesCount);
	for (unsigned i = 0; i < meshesCount; ++i)
		objReadMeshMat(meshesVers[i], meshesTris[i], fileNames[i].c_str());
	tt.endCout("读取输入网格耗时：");

	// 2. 执行meshFix：
	tt.start();
	debugDisp("执行meshFix：");
	for (unsigned i = 0; i < meshesCount; ++i)
	{
		T_MESH::TMesh::init();												// ？？？This is mandatory
		T_MESH::Basic_TMesh tMesh;
		meshMat2tMesh(tMesh, meshesVers[i], meshesTris[i]);

		unsigned versCount = tMesh.V.numels();
		unsigned trisCount = tMesh.T.numels();

		// f1. 提取最大单连通网格；
		int removedCount = tMesh.removeSmallestComponents();					// d_boundaries, d_handles, d_shells赋值

		// f2. 补洞
		int holesCount = tMesh.boundaries();				// ？？？
		int patchedCount = 0;
		if (holesCount)
			patchedCount = tMesh.fillSmallBoundaries(0, true);

		// f3. meshclean前计算环形边界数、环柄数；
		holesCount = tMesh.boundaries();

		// f4. meshClean()——去除退化结构和三角片自交——默认max_iters == 10, inner_loops == 3；
		bool flagClean = false;
		int max_iters = 20;
		int inner_loops = 6;					// 每次大循环中去除退化三角片、去除自交的迭代次数；
		bool flagIsct, flagDeg;
		T_MESH::Triangle* t;
		T_MESH::Node* m;

		//		f4.1. 
		tMesh.deselectTriangles();
		tMesh.invertSelection();

		//		f4.2.修复流程的大循环 
		for (int i = 0; i < max_iters; i++)
		{
			//		ff1. 去除退化三角片；
			flagDeg = tMesh.strongDegeneracyRemoval(inner_loops);			// 全部清除成功返回true， 否则返回false

			//		ff2. 
			tMesh.deselectTriangles();
			tMesh.invertSelection();

			//		ff3. 去除自交三角片，补洞；
			flagIsct = tMesh.strongIntersectionRemoval(inner_loops);			// 自交全部清除返回true，否则返回false;

			//		ff4. 若前两项全部清除成功，进一步检查确认：
			if (flagIsct && flagDeg)
			{
				// 遍历三角片检测是否有退化；
				for (m = tMesh.T.head(), t = (m) ? ((T_MESH::Triangle*)m->data) : NULL; m != NULL; m = m->next(), t = (m) ? ((T_MESH::Triangle*)m->data) : NULL)
					if (t->isExactlyDegenerate())
						flagIsct = false;

				// 若进一步检查没有问题，退出大循环；
				if (flagIsct)
				{
					flagClean = true;
					break;
				}
			}
		}

#ifdef LOCAL_DEBUG
		if (!flagClean)
			std::cout << "!!!meshclean() is not completed!!!" << std::endl;
#endif

		TMesh2MeshMat(outVers[i], outTris[i], tMesh);
	}
	tt.endCout("meshFix耗时：");

	// 3. 输出：
	for (unsigned i = 0; i < meshesCount; ++i)
	{
		std::string str = pathOutput + OBJfileNames[i] + std::string{ ".obj" };
		objWriteMeshMat(str.c_str(), outVers[i], outTris[i]);
	}

	debugDisp("finished.");
	getchar();

	return 0;
}


// 批量读取本地网格执行qslim边折叠精简
int testCmd_qslimDecimation(int argc, char** argv)
{
	tiktok& tt = tiktok::getInstance();
	float deciRatio = 0.5;							// 精简率

	// 生成路径：
	int   nPos;
	CString   cPath;
	GetModuleFileName(NULL, cPath.GetBufferSetLength(MAX_PATH + 1), MAX_PATH);		// 获取当前进程加载的模块的路径。
	nPos = cPath.ReverseFind('\\');
	cPath = cPath.Left(nPos);
	std::string path{ CT2CA{cPath} };
	std::string pathOBJ = path + "\\inputOBJ";
	std::string pathOutput = path + "\\outputData";

	// 1. 读取部件网格 
	tt.start();
	std::cout << "读取输入网格..." << std::endl;
	std::vector<std::string> fileNames, tmpStrVec, OBJfileNames;
	std::vector<Eigen::MatrixXd> meshesVers, outVers;
	std::vector<Eigen::MatrixXi> meshesTris, outTris;

	getFileNames(pathOBJ.c_str(), tmpStrVec, false);
	const unsigned meshesCount = tmpStrVec.size();
	OBJfileNames.reserve(meshesCount);
	for (const auto& str : tmpStrVec)
	{
		std::string tailStr = str.substr(str.size() - 4, 4);				//	".obj"
		if (".obj" == tailStr)
		{
			fileNames.push_back(str);
			unsigned index = str.find_last_of("/");
			std::string OBJfileName = str.substr(index, str.size() - index - 4);			// "/" + obj文件名，不含路径和.obj后缀；
			OBJfileNames.push_back(OBJfileName);
		}
	}

	meshesVers.resize(meshesCount);
	meshesTris.resize(meshesCount);
	outVers.resize(meshesCount);
	outTris.resize(meshesCount);
	for (unsigned i = 0; i < meshesCount; ++i)
		objReadMeshMat(meshesVers[i], meshesTris[i], fileNames[i].c_str());
	tt.endCout("读取输入网格耗时：");

	// 2. 执行qslim精简：
	tt.start();
	debugDisp("执行qslim精简：");
	for (int i = 0; i < meshesCount; ++i)
	{
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
		Eigen::VectorXi newOldVersInfo;
		int trisCount = meshesTris[i].rows();
		int tarTrisCount = std::round(trisCount * deciRatio);

		igl::qslim(meshesVers[i], meshesTris[i], tarTrisCount, outVers[i], outTris[i], newOldTrisInfo, newOldVersInfo);
	}
	tt.endCout("qslim精简耗时：");

	// 3. 输出：
	for (unsigned i = 0; i < meshesCount; ++i)
	{
		std::string str = pathOutput + OBJfileNames[i] + std::string{ ".obj" };
		objWriteMeshMat(str.c_str(), outVers[i], outTris[i]);
	}

	debugDisp("finished.");
	getchar();

	return 0;
}

 
// 批量读取本地网格执行laplace光顺：
int testCmd_laplaceFaring(int argc, char** argv)
{
	tiktok& tt = tiktok::getInstance();
	float deciRatio = 0.5;							// 精简率

	// 生成路径：
	int   nPos;
	CString   cPath;
	GetModuleFileName(NULL, cPath.GetBufferSetLength(MAX_PATH + 1), MAX_PATH);		// 获取当前进程加载的模块的路径。
	nPos = cPath.ReverseFind('\\');
	cPath = cPath.Left(nPos);
	std::string path{ CT2CA{cPath} };
	std::string pathOBJ = path + "\\inputOBJ";
	std::string pathOutput = path + "\\outputData";
	CString fileConfig = cPath + "\\config.ini";

	// 1. 读取部件网格 
	tt.start();
	std::cout << "读取输入网格..." << std::endl;
	std::vector<std::string> fileNames, tmpStrVec, OBJfileNames;
	std::vector<Eigen::MatrixXd> meshesVers, outVers;
	std::vector<Eigen::MatrixXi> meshesTris, outTris;

	getFileNames(pathOBJ.c_str(), tmpStrVec, false);
	const unsigned meshesCount = tmpStrVec.size();
	OBJfileNames.reserve(meshesCount);
	for (const auto& str : tmpStrVec)
	{
		std::string tailStr = str.substr(str.size() - 4, 4);				//	".obj"
		if (".obj" == tailStr)
		{
			fileNames.push_back(str);
			unsigned index = str.find_last_of("/");
			std::string OBJfileName = str.substr(index, str.size() - index - 4);			// "/" + obj文件名，不含路径和.obj后缀；
			OBJfileNames.push_back(OBJfileName);
		}
	}
	meshesVers.resize(meshesCount);
	meshesTris.resize(meshesCount);
	outVers.resize(meshesCount);
	outTris.resize(meshesCount);
	for (unsigned i = 0; i < meshesCount; ++i)
		objReadMeshMat(meshesVers[i], meshesTris[i], fileNames[i].c_str());
	tt.endCout("读取输入网格耗时：");

	// 2. 网格逐个执行laplace光顺：
	tt.start();
	debugDisp("执行laplace光顺：");
	float deltaLB = INIGetFloat(TEXT("deltaLB"), fileConfig);
	unsigned loopCount = INIGetInt(TEXT("laplaceFaringLoopCount"), fileConfig);
	for (int i = 0; i < meshesCount; ++i)
		laplaceFaring(outVers[i], meshesVers[i], meshesTris[i], deltaLB, loopCount);
	outTris = std::move(meshesTris);
	tt.endCout("laplace光顺：");

	// 3. 输出：
	for (unsigned i = 0; i < meshesCount; ++i)
	{
		std::string str = pathOutput + OBJfileNames[i] + std::string{ ".obj" };
		objWriteMeshMat(str.c_str(), outVers[i], outTris[i]);
	}

	debugDisp("finished.");

	return 0;
}
 


int main(int argc, char** argv)
{
	// DENSEMAT::test7();
	
	// SPARSEMAT::test0();

	 DECIMATION::test1();

	// IGL_DIF_GEO::test1();
	// IGL_GRAPH::test1();
	// IGL_SPACE_PARTITION::test0();
	// IGL_BASIC_PMP::test44();
 
	// SCIENTIFICCALC::test7();
	
	// TEST_PMP::test3();
	
	// IGL_MATH::test1();

	// DECIMATION::test1();

	// TEST_MYEIGEN::test5();

	// TEMP_TEST::test1();

	// MESH_REPAIR::test0();
 
	// TEST_DIP::test0();

	// TEST_TMESH::test44();

	//DECIMATION::testCmd_qslimDecimation(argc, argv);

	// MESH_REPAIR::testCmd_meshFix(argc, argv);

	// testCmd_laplaceFaring(argc, argv);
 
	std::cout << "main() finished." << std::endl;
}
