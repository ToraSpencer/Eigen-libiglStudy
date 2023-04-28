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
#include <string>

#include <windows.h>
#include <atlstr.h>			// 包含CString类。属于microsoft ATL(活动模板库avtive template library)
#include <atlconv.h>
#include <io.h>
#include <winuser.h>


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

#define IGL_COLLAPSE_EDGE_NULL -1

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


	void circulation(const int uEdgeIdx,	const bool ccw, const Eigen::MatrixXi& tris,	const Eigen::VectorXi& edgeUeInfo,
		const Eigen::MatrixXi& UeTrisInfo,	const Eigen::MatrixXi& UeCornersInfo,	std::vector<int>& nbrVersIdx, std::vector<int>& nbrTrisIdx)
	{
		nbrVersIdx.clear(); nbrVersIdx.reserve(10);
		nbrTrisIdx.clear(); nbrTrisIdx.reserve(10);
		const int m = edgeUeInfo.size() / 3;
		assert(m * 3 == edgeUeInfo.size());

		const auto& step = [&](const int uEdgeIdx, const int ff, int& ne, int& rv, int& nf)
		{
			assert((UeTrisInfo(uEdgeIdx, 1) == ff || UeTrisInfo(uEdgeIdx, 0) == ff) && "uEdgeIdx should touch ff");
			const int nside = UeTrisInfo(uEdgeIdx, 0) == ff ? 1 : 0;
			const int nv = UeCornersInfo(uEdgeIdx, nside);
			nf = UeTrisInfo(uEdgeIdx, nside);
			const int dir = ccw ? -1 : 1;
			rv = tris(nf, nv);
			ne = edgeUeInfo(nf + m * ((nv + dir + 3) % 3));
		};

		const int triIdx = UeTrisInfo(uEdgeIdx, 0);
		int cTriIdx = triIdx;
		int cUeIdx = uEdgeIdx;

		while (true)
		{
			int re, rv;
			step(cUeIdx, cTriIdx, cUeIdx, rv, cTriIdx);
			nbrTrisIdx.push_back(cTriIdx);
			nbrVersIdx.push_back(rv);
 
			if (cTriIdx == triIdx)
			{
				assert(cUeIdx == uEdgeIdx);
				break;
			}
		}
	}


	void getAllTheShit(const int uEdgeIdx, const Eigen::MatrixXi& uEdges, const Eigen::MatrixXd& vers, const Eigen::MatrixXi& tris, \
				std::vector<int>& nbrVersIdx_src, std::vector<int>& nbrTrisIdx_src, std::vector<int>& nbrVersIdx_des, std::vector<int>& nbrTrisIdx_des)
	{
		auto blEdgeInTri = [&](const int vaIdx0, const int vbIdx0, const int triIdx0) -> bool 
		{
			if (triIdx0 < 0)
				return false;
			const int A = tris(triIdx0, 0);
			const int B = tris(triIdx0, 1);
			const int C = tris(triIdx0, 2);
			if (A == vaIdx0 && B == vbIdx0)
				return true;
			if (B == vaIdx0 && C == vbIdx0)
				return true;
			if (C == vaIdx0 && A == vbIdx0)
				return true;

			return false;
		};

		auto getEdgeArrow = [&](const int triIdx, const int centerIdx)->std::pair<int, int>
		{
			if (triIdx < 0 || centerIdx < 0)
				return std::make_pair(-1, -1);
			
			const int A = tris(triIdx, 0);
			const int B = tris(triIdx, 1);
			const int C = tris(triIdx, 2);
			if(centerIdx == A)
				return std::make_pair(B, C);
			if (centerIdx == B)
				return std::make_pair(C, A);
			if (centerIdx == C)
				return std::make_pair(A, B);

			return std::make_pair(-1, -1);
		};

		const int vaIdx = uEdges(uEdgeIdx, 0);
		const int vbIdx = uEdges(uEdgeIdx, 1);
		std::unordered_set<int> nbrTriIdxesA = oneRingTriIdxes(vaIdx, vers, tris);
		std::unordered_set<int> nbrTriIdxesB = oneRingTriIdxes(vbIdx, vers, tris);

		int arrowStart;
		std::unordered_map<int, std::pair<int, int>> srcTrisMap;				// 以B为中心点的一圈三角片；
		std::unordered_map<int, std::pair<int, int>> desTrisMap;				// 以A为中心点的一圈三角片；

		arrowStart = vaIdx;		
		for (const auto& triIdx : nbrTriIdxesB)
		{
			std::pair<int, int> edgeArrow = getEdgeArrow(triIdx, vbIdx);
			if (edgeArrow.first < 0)
				continue;
			srcTrisMap.insert(std::make_pair(edgeArrow.first, std::make_pair(edgeArrow.second, triIdx)));
		}
		while (true) 
		{
			auto retIter = srcTrisMap.find(arrowStart);
			if (retIter == srcTrisMap.end())
				break;
			std::pair<int, int> triInfo = retIter->second;
			nbrTrisIdx_src.push_back(triInfo.second);
			arrowStart = triInfo.first;
			nbrVersIdx_src.push_back(arrowStart);
			srcTrisMap.erase(retIter);
		}

		arrowStart = vbIdx;
		for (const auto& triIdx : nbrTriIdxesA)
		{
			std::pair<int, int> edgeArrow = getEdgeArrow(triIdx, vaIdx);
			if (edgeArrow.first < 0)
				continue;
			desTrisMap.insert(std::make_pair(edgeArrow.first, std::make_pair(edgeArrow.second, triIdx)));
		}
		while (true)
		{
			auto retIter = desTrisMap.find(arrowStart);
			if (retIter == desTrisMap.end())
				break;
			nbrVersIdx_des.push_back(arrowStart);
			std::pair<int, int> triInfo = retIter->second;
			nbrTrisIdx_des.push_back(triInfo.second);
			arrowStart = triInfo.first;
			desTrisMap.erase(retIter);
		}
		std::reverse(nbrTrisIdx_des.begin(), nbrTrisIdx_des.end());
		std::reverse(nbrVersIdx_des.begin(), nbrVersIdx_des.end());
	}


	// 计算边的折叠cost，以及折叠后的新顶点——最简单的策略：cost为边长，折叠后的新顶点为中点；
	template <typename T>
	double costAndNewPos_simplest(Eigen::Matrix<T, 1, 3>& newPos, const int ueIdx, const Eigen::MatrixXi& uEdges, \
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
	{
		int vaIdx = uEdges(ueIdx, 0);
		int vbIdx = uEdges(ueIdx, 1);
		double cost = (vers.row(vaIdx) - vers.row(vbIdx)).norm();
		newPos = (vers.row(vaIdx) + vers.row(vbIdx)) / 2;

		// Force infs and nans to infinity
		if (std::isinf(cost) || cost != cost)
		{
			cost = std::numeric_limits<double>::infinity();
			newPos.setConstant(0);
		}

		return cost;
	}


	// 计算边的折叠cost，以及折叠后的新顶点——qslim策略；
	template <typename T>
	double costAndNewPos_qslim(Eigen::Matrix<T, 1, 3>& newPos, const int ueIdx, const Eigen::MatrixXi& uEdges, \
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, const std::vector<Quadric>& quadrics)
	{
		// Quadric: newPos'Ap + 2b'newPos + c,  optimal point: Ap = -b, or rather because we have row vectors: pA=-b
		Quadric quadric_p = quadrics[uEdges(ueIdx, 0)] + quadrics[uEdges(ueIdx, 1)];		// Combined quadric
		const auto& A = std::get<0>(quadric_p);
		const auto& b = std::get<1>(quadric_p);
		const auto& c = std::get<2>(quadric_p);
		newPos = -b * A.inverse();
		double cost = newPos.dot(newPos * A) + 2 * newPos.dot(b) + c;

		// Force infs and nans to infinity
		if (std::isinf(cost) || cost != cost)
		{
			cost = std::numeric_limits<double>::infinity();
			newPos.setConstant(0);
		}

		return cost;
	}


	// 判断边是否可以折叠：
	bool edge_collapse_is_valid(std::vector<int>& srcNbrIdx, std::vector<int>& desNbrIdx)
	{
		// Do we really need to check if edge is IGL_COLLAPSE_EDGE_NULL ?
		if (srcNbrIdx.size() < 2 || desNbrIdx.size() < 2)
			return false;

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

	 
	// 将网格处理成一个封闭网格——若存在边缘有向边，则将其和一个无限远点连接生成新三角片
	template <typename T>
	bool connect_boundary_to_infinity(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, Eigen::MatrixXi& trisOut,
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
		{
			versOut = vers;
			trisOut = tris;
			return true;
		}

		newTris.resize(bdrys.rows(), 3);
		newTris.col(0) = bdrys.col(1);
		newTris.col(1) = bdrys.col(0);
		newTris.col(2).array() = vers.rows();

		trisOut = tris;
		matInsertRows(trisOut, newTris);

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
 

	// 折叠索引为uEdgeIdx的无向边；
	bool collapseSingleEdge(const int uEdgeIdx, const Eigen::RowVectorXd& collapsedVer,
		std::vector<int>& nbrVersIdx_src, const std::vector<int>& nbrTrisIdx_src,
		std::vector<int>& nbrVersIdx_des, const std::vector<int>& nbrTrisIdx_des,
		Eigen::MatrixXd& vers, Eigen::MatrixXi& tris, Eigen::MatrixXi& uEdges, Eigen::MatrixXi& edges, \
		Eigen::VectorXi& edgeUeInfo, Eigen::MatrixXi& UeTrisInfo, Eigen::MatrixXi& UeCornersInfo,
		int& a_e1, int& a_e2, int& a_f1, int& a_f2)
	{
		// lambda——某一条无向边的所有信息设置为无效信息；
		const auto& kill_edge = [&uEdges, &UeCornersInfo, &UeTrisInfo](const int uEdgeIdx)
		{
			uEdges(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			uEdges(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
			UeTrisInfo(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			UeTrisInfo(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
			UeCornersInfo(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			UeCornersInfo(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
		};

		// 注意！！当前所有有向边都是正序的，即vaIdx < vbIdx; 且所有无向边都对应着两条有向边，因为边缘边已被消除；

	   // 定义无向边的起点和终点——索引小的为起点，索引大的为终点；
		const int trisCount = tris.rows();
		const int versCount = vers.rows();
		const int vaIdx = uEdges(uEdgeIdx, 0);
		const int vbIdx = uEdges(uEdgeIdx, 1);
		const std::vector<int>& nbrTriIdxes = nbrTrisIdx_src;

		// 0. 判断该边是否可折叠
		if (!edge_collapse_is_valid(nbrVersIdx_src, nbrVersIdx_des))
			return false;

		// 1. 要折叠的边两端点都替换为collapsedVer
		vers.row(vaIdx) = collapsedVer;
		vers.row(vbIdx) = collapsedVer;

		int corner = UeCornersInfo(uEdgeIdx, 0);
		int triIdx = UeTrisInfo(uEdgeIdx, 0);
		int eIdx = corner * trisCount + triIdx;
		int corner1 = (corner + 1) % 3;
		int corner2 = (corner + 2) % 3;
		int eIdx1 = triIdx + trisCount * corner1;
		int eIdx2 = triIdx + trisCount * corner2;
		const int ueIdx1 = edgeUeInfo(eIdx1);
		const int ueIdx2 = edgeUeInfo(eIdx2);

		{
			const bool flip1 = UeTrisInfo(ueIdx1, 1) == triIdx;
			const int triOppIdx1 = flip1 ? UeTrisInfo(ueIdx1, 0) : UeTrisInfo(ueIdx1, 1);
			// const int triOppIdx1 = eOppIdx1 % trisCount;
			const int cornerOpp1 = flip1 ? UeCornersInfo(ueIdx1, 0) : UeCornersInfo(ueIdx1, 1);
			const int eOppIdx1 = triOppIdx1 + trisCount * cornerOpp1;
			//const int eOppIdx1 = eOppIdx1;

			// kill edge, kill tri:
			kill_edge(ueIdx1);
			tris(triIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			tris(triIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
			tris(triIdx, 2) = IGL_COLLAPSE_EDGE_NULL;

			// 更新有向边无向边的索引映射；
			edgeUeInfo(eOppIdx1) = ueIdx2;

			// 更新其他：
			const int opp2 = (UeTrisInfo(ueIdx2, 0) == triIdx ? 0 : 1);
			UeTrisInfo(ueIdx2, opp2) = triOppIdx1;
			UeCornersInfo(ueIdx2, opp2) = cornerOpp1;

			a_e1 = ueIdx1;
			a_f1 = triIdx;
		}

		int CORNER = UeCornersInfo(uEdgeIdx, 1);
		int TRIIDX = UeTrisInfo(uEdgeIdx, 1);
		int EIDX = CORNER * trisCount + TRIIDX;
		int CORNER1 = (CORNER + 2) % 3;
		int CORNER2 = (CORNER + 1) % 3;
		int EIDX1 = TRIIDX + trisCount * CORNER1;
		int EIDX2 = TRIIDX + trisCount * CORNER2;
		const int UEIDX1 = edgeUeInfo(EIDX1);
		const int UEIDX2 = edgeUeInfo(EIDX2);

		{
			const bool flip1 = UeTrisInfo(UEIDX1, 1) == TRIIDX;
			const int TRIOPPIDX1 = flip1 ? UeTrisInfo(UEIDX1, 0) : UeTrisInfo(UEIDX1, 1);
			const int CORNEROPP1 = flip1 ? UeCornersInfo(UEIDX1, 0) : UeCornersInfo(UEIDX1, 1);
			const int EOPPIDX1 = TRIOPPIDX1 + trisCount * CORNEROPP1;

			kill_edge(UEIDX1);
			tris(TRIIDX, 0) = IGL_COLLAPSE_EDGE_NULL;
			tris(TRIIDX, 1) = IGL_COLLAPSE_EDGE_NULL;
			tris(TRIIDX, 2) = IGL_COLLAPSE_EDGE_NULL;

			edgeUeInfo(EOPPIDX1) = UEIDX2;

			const int opp2 = (UeTrisInfo(UEIDX2, 0) == TRIIDX ? 0 : 1);
			UeTrisInfo(UEIDX2, opp2) = TRIOPPIDX1;
			UeCornersInfo(UEIDX2, opp2) = CORNEROPP1;

			a_e2 = UEIDX1;
			a_f2 = TRIIDX;
		}

		// 3. 被折叠的无向边uEdgeIdx中, vaIdx被保留，vbIdx称为孤立点；
		for (auto i : nbrTriIdxes)
		{
			for (int k = 0; k < 3; k++)
			{
				if (tris(i, k) == vbIdx)
				{
					const int eIdx1 = i + trisCount * ((k + 1) % 3);
					const int eIdx2 = i + trisCount * ((k + 2) % 3);
					const int ueIdx2 = edgeUeInfo(eIdx2);
					const int ueIdx1 = edgeUeInfo(eIdx1);
					if (vbIdx == uEdges(ueIdx1, 0))
						uEdges(ueIdx1, 0) = vaIdx;
					if (vbIdx == uEdges(ueIdx1, 1))
						uEdges(ueIdx1, 1) = vaIdx;

					// Skip if we just handled this edge (claim: this will be all except for the first non-trivial face)
					if (ueIdx2 != -1)
					{
						if (vbIdx == uEdges(ueIdx2, 0))
							uEdges(ueIdx2, 0) = vaIdx;
						if (vbIdx == uEdges(ueIdx2, 1))
							uEdges(ueIdx2, 1) = vaIdx;
					}

					tris(i, k) = vaIdx;
					break;
				}
			}
		}

		// 4. Finally,  "remove" this edge and its information
		kill_edge(uEdgeIdx);

		return true;
	}
 


	// 最简单的边折叠精简算法——使用边长来度量折叠损耗；
	void test1()
	{
		Eigen::MatrixXd vers, versOut, versOri;
		Eigen::MatrixXi tris, trisOut, trisTmp, trisOri;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi edges, uEdges, UeTrisInfo, UeCornersInfo, ueEdgeInfo;
		Eigen::VectorXi timeStamps;
		Eigen::MatrixXd collapsedVers;
		Eigen::VectorXd costs;
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// 优先队列；(无向边的cost, 无向边索引, 时间戳)
		bool clean_finish = false;
		tiktok& tt = tiktok::getInstance();

		igl::readOBJ("E:/材料/tooth.obj", versOri, trisOri);
		igl::writeOBJ("E:/meshIn.obj", versOri, trisOri);

		tt.start();
		const unsigned trisCountOri = trisOri.rows();
		int trisCountNew = trisCountOri;
		int tarTrisCount = std::round(trisCountOri / 3);

		// 0. 检测是否有非流形有向边，有则直接退出；
		Eigen::MatrixXi nmnUedges;
		nonManifoldEdges(nmnUedges, trisOri);
		if (nmnUedges.size() > 0)
			return;

		// 1. 将网格处理成一个封闭网格——若存在边缘有向边，则将其和一个无限远点连接生成新三角片
		if (!connect_boundary_to_infinity(vers, tris, versOri, trisOri))
			return;

		// 2. 计算无向边信息：
		if (!getEdges(edges, tris))
			return;
		if (!getUedges(uEdges, edgeUeInfo, ueEdgeInfo, edges))
			return;
		if (!getUeInfos(UeTrisInfo, UeCornersInfo, edgeUeInfo, uEdges, tris))
			return;

		// 3. 计算每条无向边的cost值（边长），以此为优先级存入优先队列
		const unsigned ueCount = uEdges.rows();
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());
		collapsedVers.resize(uEdges.rows(), vers.cols());
		costs.resize(uEdges.rows());
		igl::parallel_for(ueCount, [&](const int ueIdx)
			{
				Eigen::RowVector3d newPos;
				double cost = costAndNewPos_simplest(newPos, ueIdx, uEdges, vers, tris);
				collapsedVers.row(ueIdx) = newPos;
				costs(ueIdx) = cost;
			},
			10000);
		for (int i = 0; i < ueCount; i++)
			pQueue.emplace(costs(i), i, 0);
 
		// 4. 边折叠的循环：
		int uEdgeIdx0;							// 优先队列首部的边；
		int ueIdx1, ueIdx02, triIdx1, triIdx02;
		unsigned loopCount = 0;
		while (true)
		{
			bool collapsed = true;
			std::tuple<double, int, int> edgeTuple;							// 优先队列首部的元素；
			std::vector<int> nbrTrisIdx_src, nbrVersIdx_src;				// vaIdx的1领域的三角片索引、顶点索引
			std::vector<int>  nbrTrisIdx_des, nbrVersIdx_des;			// vbIdx的1领域的三角片索引、顶点索引

			// w1. 取队首元素，队首元素出队：
			while (true)
			{
				// ww1. 若队列为空，退出循环；
				if (pQueue.empty())
					assert("边折叠光了");

				// ww2. 取队首元素，队首元素出队：
				edgeTuple = pQueue.top();
				if (std::get<0>(edgeTuple) == std::numeric_limits<double>::infinity())
					assert("队首边的cost是无穷大");
				pQueue.pop();
				uEdgeIdx0 = std::get<1>(edgeTuple);          // 队首的无向边索引;

				// ww3. Check if matches timestamp
				if (std::get<2>(edgeTuple) == timeStamps(uEdgeIdx0))
					break;
			}

			// w2. 计算当前边两端点1领域的顶点、三角片：
			getAllTheShit(uEdgeIdx0, uEdges, vers, tris, nbrVersIdx_src, nbrTrisIdx_src, nbrVersIdx_des, nbrTrisIdx_des);

			//	w3. 折叠队首的边：
			collapsed = collapseSingleEdge(uEdgeIdx0, collapsedVers.row(uEdgeIdx0), nbrVersIdx_src, nbrTrisIdx_src, \
				nbrVersIdx_des, nbrTrisIdx_des, vers, tris, uEdges, edges, edgeUeInfo, UeTrisInfo, UeCornersInfo, \
				ueIdx1, ueIdx02, triIdx1, triIdx02);
			loopCount++;

			//	w4. post_collapses; 默认情形下什么都不做；

			// w5. 折叠操作之后，更新相关timeStamp， 更新相关边的cost值
			if (collapsed)
			{
				// w5.1. Erase the two,  other collapsed uEdges by marking their timestamps as -1
				timeStamps(ueIdx1) = -1;
				timeStamps(ueIdx02) = -1;

				// w5.2.确定被折叠的边1领域的所有三角片
				std::vector<int> nbrTrisIdx;							// 被折叠的边1领域的所有三角片索引；
				nbrTrisIdx.reserve(nbrTrisIdx_src.size() + nbrTrisIdx_des.size());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_src.begin(), nbrTrisIdx_src.end());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_des.begin(), nbrTrisIdx_des.end());
				std::sort(nbrTrisIdx.begin(), nbrTrisIdx.end());
				nbrTrisIdx.erase(std::unique(nbrTrisIdx.begin(), nbrTrisIdx.end()), nbrTrisIdx.end());

				// w5.3. 确定被折叠的边1领域的所有无向边 
				std::vector<int> nbrUeIdxes;						// 被折叠的边1领域的所有无向边索引；
				nbrUeIdxes.reserve(3 * nbrTrisIdx.size());
				for (auto& triIdx : nbrTrisIdx)
				{
					if (tris(triIdx, 0) == IGL_COLLAPSE_EDGE_NULL || tris(triIdx, 1) == IGL_COLLAPSE_EDGE_NULL || tris(triIdx, 2) == IGL_COLLAPSE_EDGE_NULL)
						continue;

					for (int i = 0; i < 3; i++)
					{
						const int eIdx = i * tris.rows() + triIdx;
						const int ueIdx = edgeUeInfo(eIdx);
						nbrUeIdxes.push_back(ueIdx);
					}
				}
				std::sort(nbrUeIdxes.begin(), nbrUeIdxes.end());
				nbrUeIdxes.erase(std::unique(nbrUeIdxes.begin(), nbrUeIdxes.end()), nbrUeIdxes.end());             // 去除重复元素：
 
				// w5.4. 更新被折叠的边1领域的所有无向边的信息——cost值、折叠后顶点、时间戳；
				for (auto& ueIdx : nbrUeIdxes)
				{
					// wf1. 计算边折叠的cost值，及折叠后的顶点坐标：
					Eigen::RowVector3d newPos;
					double cost = costAndNewPos_simplest(newPos, ueIdx, uEdges, vers, tris);

					// wf2. 更新的无向边的时间戳+1
					timeStamps(ueIdx)++;

					// wf3. 更新后的无向边被重新插入队列，之前插入的由于时间戳不匹配所以已经无效；
					pQueue.emplace(cost, ueIdx, timeStamps(ueIdx));
					collapsedVers.row(ueIdx) = newPos;
				}
			}
			else
				assert("edge collapse failed.");

			// w6. 判断是否满足终止条件stopping_condition，若满足则跳出边折叠循环；
			if (collapsed)
			{
				// w6.1. 更新三角片数
				if (triIdx1 < trisCountOri)
					trisCountNew -= 1;
				if (triIdx02 < trisCountOri)
					trisCountNew -= 1;

				// w6.2. 若三角片数达到目标，则退出边折叠循环；
				bool stopConditionFlag = (trisCountNew <= (int)tarTrisCount);
				if (stopConditionFlag)
				{
					clean_finish = true;
					break;
				}
			}
			else
				assert("edge collapse failed.");		 // 边折叠失败，退出循环：
		}


		// 5. 删除所有含有IGL_COLLAPSE_EDGE_NULL的三角片，以及包含无穷远点的三角片
		trisTmp.resize(tris.rows(), 3);
		int index = 0;
		for (int i = 0; i < tris.rows(); i++)
		{
			int vaIdx = tris(i, 0);
			int vbIdx = tris(i, 1);
			int vcIdx = tris(i, 2);
			bool continueFlag = false;
			Eigen::Matrix3d triVers;					// 当前三角片的三个顶点；
			double* coorPtr = nullptr;
			if (vaIdx == vbIdx || vbIdx == vcIdx || vcIdx == vaIdx)
				continue;
			triVers.row(0) = vers.row(vaIdx);
			triVers.row(1) = vers.row(vbIdx);
			triVers.row(2) = vers.row(vcIdx);
			coorPtr = triVers.data();
			for (int k = 0; k < 9; ++k)
			{
				if (std::isinf(*coorPtr++))
				{
					continueFlag = true;
					break;
				}
			}
			if (continueFlag)
				continue;
			trisTmp.row(index) = tris.row(i);
			index++;
		}
		trisTmp.conservativeResize(index, trisTmp.cols());					  // 这里相当于shrink_to_fit(); 


		// 6. 删除孤立顶点，并更新三角片：
		std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, trisTmp);
		if (!removeIsoVers(versOut, trisOut, vers, trisTmp, isoVerIdxes))
			return;

		debugDisp("loopCount == ", loopCount);
		tt.endCout("精简耗时：");
		igl::writeOBJ("E:/decimateOutput.obj", versOut, trisOut);
		std::cout << "finished." << std::endl;
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

		igl::edge_flaps(trisCopy, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());						// Could reserve with https://stackoverflow.com/a/29236236/148668        
		collapsedVers.resize(uEdges.rows(), versCopy.cols());				// If an edge were collapsed, we'd collapse it to these points:

		// 检测是否有非流形边：
		{
			Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> BF;
			Eigen::Array<bool, Eigen::Dynamic, 1> BE;
			if (!igl::is_edge_manifold(trisCopy, uEdges.rows(), edgeUeInfo, BF, BE))
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
			Eigen::MatrixXi trisTmp(trisCopy.rows(), 3);
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
					trisTmp.row(m) = trisCopy.row(i);
					newOldTrisInfo(m) = i;
					m++;
				}
			}
			trisTmp.conservativeResize(m, trisTmp.cols());              // 这里相当于shrink_to_fit();
			newOldTrisInfo.conservativeResize(m);

			// 3.2. 删除网格中的孤立顶点：
			Eigen::MatrixXd versOut;
			Eigen::MatrixXi trisOut;
			igl::remove_unreferenced(vers, trisTmp, versOut, trisOut, _1, newOldVersInfo);
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
		Eigen::MatrixXi trisTmp(trisCopy.rows(), 3);
		Eigen::VectorXi _1;
		newOldTrisInfo.resize(trisCopy.rows());
		int m = 0;
		for (int i = 0; i < trisCopy.rows(); i++)
		{
			if (trisCopy(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
				trisCopy(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
				trisCopy(i, 2) != IGL_COLLAPSE_EDGE_NULL)
			{
				trisTmp.row(m) = trisCopy.row(i);
				newOldTrisInfo(m) = i;
				m++;
			}
		}
		trisTmp.conservativeResize(m, trisTmp.cols());              // 这里相当于shrink_to_fit();
		newOldTrisInfo.conservativeResize(m);

		// 3.2. 删除网格中的孤立顶点：
		Eigen::MatrixXd versOut;
		Eigen::MatrixXi trisOut;
		igl::remove_unreferenced(vers, trisTmp, versOut, trisOut, _1, newOldVersInfo);
		igl::writeOBJ("E:/decimating_output.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}

 
	// 解析qslim()
	void test0000()
	{
		Eigen::MatrixXd vers, versOut, vers0, versOri;
		Eigen::MatrixXi tris, trisOut, trisTmp, trisOri;
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
		Eigen::VectorXi newOldVersInfo;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi edges, uEdges, ueEdgeInfo, UeTrisInfo, UeCornersInfo;
		Eigen::VectorXi timeStamps;
		Eigen::MatrixXd collapsedVers;
		Eigen::VectorXd costs;
		std::vector<Quadric> quadrics;													// 每个顶点的Q矩阵；
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// 优先队列；
		Eigen::VectorXi _1, I2;											// 最后删除孤立点和内补三角片时用到；
		bool clean_finish = false;
		tiktok& tt = tiktok::getInstance();

		// for debug:
		std::vector<triplet<double>> versMoni;
		std::vector<triplet<int>>	trisMoni;
		std::vector<doublet<int>> edgesMoni;
		std::vector<int>			vecMoni;
		int retIdx = -1;

		tt.start();
		igl::readOBJ("E:/材料/jawMeshDense.obj", versOri, trisOri);
		tt.endCout("elapsed time of loading mesh is: ");
		igl::writeOBJ("E:/meshIn.obj", versOri, trisOri);

		unsigned versCount = versOri.rows();
		unsigned trisCount = trisOri.rows();
		int trisCountNew = trisCount;
		int tarTrisCount = 30000;

		// 0. 检测是否有非流形有向边，有则直接退出；
		Eigen::MatrixXi nmnUedges;
		nonManifoldEdges(nmnUedges, trisOri);
		if (nmnUedges.size() > 0)
			return;

		// 1. 将网格处理成一个封闭网格——若存在边缘有向边，则将其和一个无限远点连接生成新三角片
		connect_boundary_to_infinity(vers, tris, versOri, trisOri);

		// 1.1 计算无向边信息：
		getEdges(edges, tris);
		getUedges(uEdges, edgeUeInfo, ueEdgeInfo, edges);
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
		int ueIdx1, ueIdx02, triIdx1, triIdx02;
		while (true)
		{
			bool collapsed = true;
			std::tuple<double, int, int> edgeTuple;							// 优先队列首部的元素；
			std::vector<int> nbrTrisIdx_src, nbrVersIdx_src;				// vaIdx的1领域的三角片索引、顶点索引
			std::vector<int>  nbrTrisIdx_des, nbrVersIdx_des;			// vbIdx的1领域的三角片索引、顶点索引

			// w1. 取队首元素，队首元素出队：
			while (true)
			{
				// ww1. 若队列为空，退出循环；
				if (pQueue.empty())
					assert("边折叠光了");

				// ww2. 取队首元素，队首元素出队：
				edgeTuple = pQueue.top();
				if (std::get<0>(edgeTuple) == std::numeric_limits<double>::infinity())
					assert("队首边的cost是无穷大");
				pQueue.pop();
				uEdgeIdx0 = std::get<1>(edgeTuple);          // 队首的无向边索引;

				// ww3. Check if matches timestamp
				if (std::get<2>(edgeTuple) == timeStamps(uEdgeIdx0))
					break;
			}

			// w2. 计算当前边两端点1领域的顶点、三角片：
			int vaIdx = uEdges(uEdgeIdx0, 0);
			int vbIdx = uEdges(uEdgeIdx0, 1);
			int smallIdx = (vaIdx < vbIdx ? vaIdx: vbIdx);
			igl::circulation(uEdgeIdx0, true, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_src, nbrTrisIdx_src);
			igl::circulation(uEdgeIdx0, false, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_des, nbrTrisIdx_des);

			//	w3. 折叠队首的边：
			collapsed = collapseSingleEdge(uEdgeIdx0, collapsedVers.row(uEdgeIdx0), nbrVersIdx_src, nbrTrisIdx_src, \
				nbrVersIdx_des, nbrTrisIdx_des, vers, tris, uEdges, edges, edgeUeInfo, UeTrisInfo, UeCornersInfo, ueIdx1, ueIdx02, triIdx1, triIdx02);

			//	w4. post_collapses——更新q矩阵
			if (collapsed)
				quadrics[smallIdx] = quadrics[vaIdx] + quadrics[vbIdx];

			// w5. 折叠操作之后，更新相关timeStamp， 更新相关边的cost值
			if (collapsed)
			{
				// w5.1. Erase the two,  other collapsed uEdges by marking their timestamps as -1
				timeStamps(ueIdx1) = -1;
				timeStamps(ueIdx02) = -1;

				// w5.2.确定被折叠的边1领域的所有三角片
				std::vector<int> nbrTrisIdx;							// 被折叠的边1领域的所有三角片索引；
				nbrTrisIdx.reserve(nbrTrisIdx_src.size() + nbrTrisIdx_des.size());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_src.begin(), nbrTrisIdx_src.end());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_des.begin(), nbrTrisIdx_des.end());
				std::sort(nbrTrisIdx.begin(), nbrTrisIdx.end());
				nbrTrisIdx.erase(std::unique(nbrTrisIdx.begin(), nbrTrisIdx.end()), nbrTrisIdx.end());

				// w5.3. 确定被折叠的边1领域的所有无向边 
				std::vector<int> nbrUeIdxes;						// 被折叠的边1领域的所有无向边索引；
				nbrUeIdxes.reserve(3 * nbrTrisIdx.size());
				for (auto& triIdx : nbrTrisIdx)
				{
					if (tris(triIdx, 0) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 1) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 2) != IGL_COLLAPSE_EDGE_NULL)
					{
						for (int i = 0; i < 3; i++)
						{
							const int eIdx = i * tris.rows() + triIdx;
							const int ueIdx = edgeUeInfo(eIdx);
							nbrUeIdxes.push_back(ueIdx);
						}
					}
				}
				std::sort(nbrUeIdxes.begin(), nbrUeIdxes.end());
				nbrUeIdxes.erase(std::unique(nbrUeIdxes.begin(), nbrUeIdxes.end()), nbrUeIdxes.end());             // 去除重复元素：

				// w5.4. 更新被折叠的边1领域的所有无向边的信息——cost值、折叠后顶点、时间戳；
				for (auto& ueIdx : nbrUeIdxes)
				{
					// wf1. 计算边折叠的cost值，及折叠后的顶点坐标：
					Eigen::RowVector3d newPos;
					double cost = costAndNewPos_qslim(newPos, ueIdx, uEdges, vers, tris, quadrics);

					// wf2. 更新的无向边的时间戳+1
					timeStamps(ueIdx)++;

					// wf3. 更新后的无向边被重新插入队列，之前插入的由于时间戳不匹配所以已经无效；
					pQueue.emplace(cost, ueIdx, timeStamps(ueIdx));
					collapsedVers.row(ueIdx) = newPos;
				}
			}
			else
				assert("edge collapse failed.");

			// w6. 判断是否满足终止条件stopping_condition，若满足则跳出边折叠循环；
			if (collapsed)
			{
				// w6.1. 更新三角片数
				if (triIdx1 < trisCount)
					trisCountNew -= 1;
				if (triIdx02 < trisCount)
					trisCountNew -= 1;

				// w6.2. 若三角片数达到目标，则退出边折叠循环；
				bool stopConditionFlag = (trisCountNew <= (int)tarTrisCount);
				if (stopConditionFlag)
				{
					clean_finish = true;
					break;
				}
			}
			else
				assert("edge collapse failed.");		 // 边折叠失败，退出循环：
		}

		// 5. 删除所有含有IGL_COLLAPSE_EDGE_NULL的三角片，以及包含无穷远点的三角片
		trisTmp.resize(tris.rows(), 3);
		int index = 0;
		for (int i = 0; i < tris.rows(); i++)
		{
			int vaIdx = tris(i, 0);
			int vbIdx = tris(i, 1);
			int vcIdx = tris(i, 2);
			bool continueFlag = false;
			Eigen::Matrix3d triVers;					// 当前三角片的三个顶点；
			double* coorPtr = nullptr;
			if (vaIdx == IGL_COLLAPSE_EDGE_NULL || vbIdx == IGL_COLLAPSE_EDGE_NULL || vcIdx == IGL_COLLAPSE_EDGE_NULL)
				continue;
			triVers.row(0) = vers.row(vaIdx);
			triVers.row(1) = vers.row(vbIdx);
			triVers.row(2) = vers.row(vcIdx);
			coorPtr = triVers.data();
			for (int k = 0; k < 9; ++k)
			{
				if (std::isinf(*coorPtr++))
				{
					continueFlag = true;
					break;
				}
			}
			if (continueFlag)
				continue;
			trisTmp.row(index) = tris.row(i);
			index++;
		}
		trisTmp.conservativeResize(index, trisTmp.cols());					  // 这里相当于shrink_to_fit(); 

		// 6. 删除孤立顶点，并更新三角片：
		std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, trisTmp);
		if (!removeIsoVers(versOut, trisOut, vers, trisTmp, isoVerIdxes))
			return;

		igl::writeOBJ("E:/qslimOutput.obj", versOut, trisOut);
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

		//// for debug
		//Eigen::Vector4d verColVec = versExt0.row(3).transpose();		
		//Eigen::VectorXd signedDises = planeCoeff0 * verColVec;			// 索引为verIdx的simp网格中的顶点，到ori网格中所有三角片平面的符号距离：
		//Eigen::VectorXd sqrDises = signedDises.array() * signedDises.array();

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


} 



////////////////////////////////////////////////////////////////////////////// TEST: 网格缺陷的检测和修复：
namespace MESH_REPAIR 
{
	// 检测inputOBJ文件夹中网格的边缘边、非流形边、孤立点、重复点、洞、重叠三角片等缺陷：
	int testCmd_meshDefectsDetect(int argc, char** argv)
	{
		tiktok& tt = tiktok::getInstance();
		CString   cPath, fileConfig;
		std::string path, pathOBJ, pathOutput;
		std::vector<std::string> fileNames, tmpStrVec, OBJfileNames;
		bool debugFlag = false;
		int meshesCount = 0;
		std::stringstream ss;

		// 00. 读取路径、参数；
		{
			GetModuleFileName(NULL, cPath.GetBufferSetLength(MAX_PATH + 1), MAX_PATH);		// 获取当前进程加载的模块的路径。
			int nPos = cPath.ReverseFind('\\');
			cPath = cPath.Left(nPos);
			path = CT2CA{ cPath };
			pathOBJ = path + "\\inputOBJ";
			pathOutput = path + "\\outputData";
			fileConfig = cPath + "\\config.ini";

			// 读取配置文件中的参数：
			unsigned debugInt = INIGetInt(TEXT("debugFlag"), fileConfig);
			debugFlag = debugInt > 0 ? true : false;
			if (debugFlag)
				debugDisp("Debug mode: ");

			getFileNames(pathOBJ.c_str(), tmpStrVec, false);
			meshesCount = tmpStrVec.size();
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
		}

		// 0. 读取输入网格
		tt.start();
		std::cout << "读取输入网格..." << std::endl; 
		std::vector<Eigen::MatrixXd> meshesVers;
		std::vector<Eigen::MatrixXi> meshesTris;
		{ 
			meshesVers.resize(meshesCount);
			meshesTris.resize(meshesCount);
			for (unsigned i = 0; i < meshesCount; ++i)
				objReadMeshMat(meshesVers[i], meshesTris[i], fileNames[i].c_str());
			tt.endCout("读取输入网格耗时：");
		}

		// 2. 检测网格的循环：
		for (int i = 0; i < meshesCount; ++i) 
		{
			std::stringstream ss;
			Eigen::MatrixXd& vers = meshesVers[i];
			Eigen::MatrixXi& tris = meshesTris[i];

			// f0. 检测是否存在非法三角片、重复三角片；
			if (removeSickDupTris(vers, tris) > 0)
			{
				debugDisp(OBJfileNames[i], ".obj ！！！存在非法三角片或重复三角片。");
				ss.str("");
				ss << OBJfileNames[i] << "_meshNoSickDupTris";
				debugWriteMesh(ss.str().c_str(), vers, tris);
			}

			// f1. 检测孤立顶点：
			int versCount = vers.rows();
			int trisCount = tris.rows();
			std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, tris);
			if (!isoVerIdxes.empty())
			{
				Eigen::MatrixXd isoVers;
				Eigen::MatrixXd vers1;
				Eigen::MatrixXi tris1;
				debugDisp(OBJfileNames[i], ".obj ！！！存在孤立顶点，isoVerIdxes.size() == ", isoVerIdxes.size());
				subFromIdxVec(isoVers, vers, isoVerIdxes);

				ss.str("");
				ss << OBJfileNames[i] << "_isoVers";
				debugWriteVers(ss.str().c_str(), isoVers);
				removeIsoVers(vers1, tris1, vers, tris, isoVerIdxes);
				vers = vers1;
				tris = tris1;
				versCount = vers.rows();
				trisCount = tris.rows();
				ss.str("");
				ss << OBJfileNames[i] << "_meshNoIsoVers";
				debugWriteMesh(ss.str().c_str(), vers, tris);
			}

			// f2. 计算边数据， 检测边缘有向边：
			Eigen::Index minEdgeIdx = 0;
			double minLen = 0;
			Eigen::MatrixXd edgeArrows, normals;
			Eigen::MatrixXi edges;
			Eigen::MatrixXi bdrys, bdryTris;
			Eigen::VectorXd edgesLen;
			std::vector<int> bdryTriIdxes;
			getEdges(edges, tris);
			getEdgeArrows(edgeArrows, edges, vers);
			edgesLen = edgeArrows.rowwise().norm();
			minLen = edgesLen.minCoeff(&minEdgeIdx);
			debugDisp(OBJfileNames[i], ".obj minimum edge len is ", minLen);
			if (!bdryEdges(bdrys, bdryTriIdxes, tris))
			{
				debugDisp("error! bdryEdges() run failed.");
				return -1;
			}
			if (bdrys.rows() > 0)
			{
				debugDisp(OBJfileNames[i], ".obj ！！！存在边缘边；bdrys.rows() == ", bdrys.rows()); 
				subFromIdxVec(bdryTris, tris, bdryTriIdxes);

				ss.str("");
				ss << OBJfileNames[i] << "_bdrys";
				debugWriteEdges(ss.str().c_str(), bdrys, vers);
				ss.str("");
				ss << OBJfileNames[i] << "_bdryTris";
				debugWriteMesh(ss.str().c_str(), vers, bdryTris);

				// 1.1 若存在边缘，找洞：
			}

			// f3. 检测非流形有向边：
			Eigen::MatrixXi nmnUedges;
			int nmnCount = nonManifoldEdges(nmnUedges, tris);
			if (nmnCount < 0)
			{
				debugDisp("error! nonManifoldEdges() run failed.");
				return -1;
			}
			if (nmnUedges.rows() > 0)
			{
				Eigen::MatrixXd	nmnVers;
				std::unordered_set<int> nmnVerIdxes;
				int* ptrInt = nullptr;
				debugDisp(OBJfileNames[i], ".obj ！！！存在非流形边，nmnUedges.rows() == ", nmnUedges.rows());

				ss.str("");
				ss << OBJfileNames[i] << "_nmnEdges";
				debugWriteEdges(ss.str().c_str(), nmnUedges, vers);
				ptrInt = nmnUedges.data();
				for (int i = 0; i < nmnUedges.size(); ++i)
					nmnVerIdxes.insert(*ptrInt++);
				subFromIdxCon(nmnVers, vers, nmnVerIdxes);

				ss.str("");
				ss << OBJfileNames[i] << "_nmnVers";
				debugWriteVers(ss.str().c_str(), nmnVers);
			}
 
			// f4. 检测重叠三角片：
			std::vector<std::pair<int, int>> opTrisPairs;
			int olCount = findOverLapTris(opTrisPairs, vers, tris);
			if(olCount > 0)
				debugDisp(OBJfileNames[i], ".obj 重叠三角片对数：", olCount);

			// f5. 检测退化边
			double degEdgeThreshold = 1e-3;
			Eigen::VectorXi degEdgeFlags = checkDegEdges(edges, edgeArrows, vers, tris, degEdgeThreshold);
			unsigned degEdgesCount = degEdgeFlags.sum();
			if (degEdgesCount > 0)
				debugDisp(OBJfileNames[i], ".obj degenerate edges count == ", degEdgesCount);

			// f6. 检测退化三角片：
			Eigen::VectorXd trisAreaVec;
			std::vector<unsigned> degenTriIdxes;
			const double eps = 10e-9;
			if (!trisArea(trisAreaVec, vers, tris))
			{
				debugDisp("error! trisArea() run failed.");
				return -1;
			}
			for (unsigned i = 0; i < trisCount; ++i)
				if (trisAreaVec(i) < eps)
					degenTriIdxes.push_back(i);
			if (!degenTriIdxes.empty())
			{
				Eigen::MatrixXi deTris;
				debugDisp(OBJfileNames[i], ".obj degenTriIdxes.size() == ", degenTriIdxes.size());
				subFromIdxVec(deTris, tris, degenTriIdxes);
			}

			// f7. 提取所有单连通区域：
			int scCount = 0;							// 顶点单连通区域个数；
			int sctCount = 0;							// 三角片单连通区域个数；
			Eigen::SparseMatrix<int> adjSM, adjSM_eCount, adjSM_eIdx;
			Eigen::VectorXi connectedLabels, connectedCount, connectedTriLabels, connectedTriCount;
			Eigen::VectorXi connectedLabelsSCT, connectedCountSCT;
			adjMatrix(adjSM_eCount, adjSM_eIdx, tris);
			adjSM = adjSM_eCount;
			traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
				{
					if (iter.value() > 0)
						iter.valueRef() = 1;
				});
			scCount = simplyConnectedRegion(adjSM, connectedLabels, connectedCount);
			sctCount = simplyTrisConnectedRegion(connectedLabelsSCT, connectedCountSCT, tris);
			if (scCount > 1 || sctCount > 1)
			{
				if (scCount != sctCount)
				{
					debugDisp(OBJfileNames[i], ".obj ！！！存在奇异点。");
					debugDisp(OBJfileNames[i], ".obj scCount == ", scCount);
				}
				debugDisp(OBJfileNames[i], ".obj sctCount == ", sctCount);
			}
		}

		debugDisp("finished.");
		getchar();

		return 0;
	}

	// 查找hole, gap，并尝试修补：
	void test0() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris, bdrys;
		std::vector<int> bdryTriIdxes;

		objReadMeshMat(vers, tris, "E:/材料/holeMeshIsctFree2.obj");
 

		//std::vector<std::vector<int>> holes;
		//std::vector<std::vector<int>> bdrySegs;
		//if (findHolesBdrySegs(holes, bdrySegs, vers, tris) < 0)
		//{		
		//	debugDisp("error!!! findHolesBdrySegs() failed.");
		//	return;
		//}

		debugDisp("finished.");
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
		std::vector<int> sickTriIdxes;
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


		// 3. 退化三角片ABC，若最长边是AB, 对边所在的三角片为BAX, 则BAX分解为BCX和CAX，
		for (unsigned i = 0; i < degCount; ++i) 
		{
			unsigned oppTriIdx = edgesMap.find(longEdgeOppCodes[i])->second;		

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


	// 批量读取inputOBJ文件夹中的网格，调用MESHFIX修复——不能输入非流形网格，否则修复阶段可能会死循环；
	int testCmd_meshFix(int argc, char** argv)
	{
		tiktok& tt = tiktok::getInstance(); 
		CString   cPath, fileConfig;
		std::string path, pathOBJ, pathOutput;
		std::vector<std::string> fileNames, tmpStrVec, OBJfileNames;  
		bool debugFlag = false;
		int meshesCount = 0;
		std::stringstream ss;

		// 0. 读取路径、参数；
		{
			GetModuleFileName(NULL, cPath.GetBufferSetLength(MAX_PATH + 1), MAX_PATH);		// 获取当前进程加载的模块的路径。
			int nPos = cPath.ReverseFind('\\');
			cPath = cPath.Left(nPos);
			path = CT2CA{ cPath };
			pathOBJ = path + "\\inputOBJ";
			pathOutput = path + "\\outputData";
			fileConfig = cPath + "\\config.ini";

			// 读取配置文件中的参数：
			unsigned debugInt = INIGetInt(TEXT("debugFlag"), fileConfig);
			debugFlag = debugInt > 0 ? true : false;
			if (debugFlag)
				debugDisp("Debug mode: ");

			getFileNames(pathOBJ.c_str(), tmpStrVec, false);
			meshesCount = tmpStrVec.size();
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
		}

		// 1. 读取网格 
		std::vector<Eigen::MatrixXd> meshesVers;
		std::vector<Eigen::MatrixXi> meshesTris;
		{
			tt.start();
			std::cout << "读取输入网格..." << std::endl; 
			meshesVers.resize(meshesCount);
			meshesTris.resize(meshesCount); 
			for (int i = 0; i < meshesCount; ++i)
				objReadMeshMat(meshesVers[i], meshesTris[i], fileNames[i].c_str());
			tt.endCout("读取输入网格耗时：");
		}

		// 2. 执行meshFix，输出： 
		debugDisp("执行meshFix：");
		bool flagClean = false;
		int max_iters = 20;					// 大循环次数；
		int inner_loops = 6;					// 每次大循环中去除退化三角片、去除自交的迭代次数；
		bool flagIsct = false;
		bool flagDeg = false;
		for (int i = 0; i < meshesCount; ++i)
		{
			// f0. 表象转换：
			T_MESH::TMesh::init();												// ？？？This is mandatory
			T_MESH::Basic_TMesh tMesh;
			T_MESH::TMesh::quiet = false;						// true——不要在控制台上打印修复信息；

			Eigen::MatrixXd versOut;
			Eigen::MatrixXi trisOut;
			meshMat2tMesh(tMesh, meshesVers[i], meshesTris[i]);
			int versCount = tMesh.V.numels();
			int trisCount = tMesh.T.numels();

			// f1. 提取最大单连通网格；
			int removedCount = tMesh.removeSmallestComponents();					// d_boundaries, d_handles, d_shells赋值

			// f2. 补洞
			int holesCount = tMesh.boundaries();				
			int patchedCount = 0;
			if (holesCount)
				patchedCount = tMesh.fillSmallBoundaries(0, true);

			// f3. meshclean前计算环形边界数、环柄数；
			holesCount = tMesh.boundaries();

			// f4. meshClean()——去除退化结构和三角片自交——默认max_iters == 10, inner_loops == 3；
			T_MESH::Triangle* t;
			T_MESH::Node* m;

			//		f4.1. 
			tMesh.deselectTriangles();
			tMesh.invertSelection();

			//		f4.2.修复流程的大循环 
			flagDeg = false;
			flagClean = false;
			for (int k = 0; k < max_iters; k++)
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
 
			//		f4.3 检测是否已经修复完全——无自交、无退化三角片；
			if (!flagClean)
				debugDisp("warning!!! ", OBJfileNames[i], ".obj meshclean() is not completed.");

			// f5. 写输出数据
			TMesh2MeshMat(versOut, trisOut, tMesh);
			ss.str("");
			ss << pathOutput << OBJfileNames[i] << "_meshFixed.obj";
			objWriteMeshMat(ss.str().c_str(), versOut, trisOut);
			if (debugFlag)
			{
				ss.str("");
				ss << OBJfileNames[i] << "_meshFixed";
				debugWriteMesh(ss.str().c_str(), versOut, trisOut);
			}
		} 

		debugDisp("finished.");
		getchar();

		return 0;
	}


	// 使用meshFix去除退化三角片、补洞；
	void test4() 
	{
		T_MESH::TMesh::init();																// ？？？This is mandatory
		T_MESH::Basic_TMesh mesh;
		Eigen::MatrixXd versMat;
		Eigen::MatrixXi trisMat;
		objReadMeshMat(versMat, trisMat, "E:/材料/meshDegEdges.obj");				// TMesh的load方法担心会有精度问题，反正save方法肯定有；
		meshMat2tMesh(mesh, versMat, trisMat);		 
		mesh.save("E:/meshFixInput.obj"); 
		debugWriteMesh("meshFixInput_matRepr", versMat, trisMat); 

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

		// 4. 
		bool flagClean = false;
		int max_iters = 20;
		int inner_loops = 6;			 
		bool flagDeg = false;
		bool flagNoHoles = false;
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
 
			//		f4. 若前两项全部清除成功，进一步检查确认：
			if (flagDeg)
			{
				// 检查是否有洞
				holesCount = mesh.boundaries();				 
				if (holesCount)
				{ 
					T_MESH::TMesh::warning("Patching holes\n");
					mesh.fillSmallBoundaries(0, true);
					continue;
				}

				// 遍历三角片检测是否有退化；
				for (m = mesh.T.head(), t = (m) ? ((T_MESH::Triangle*)m->data) : NULL; m != NULL; m = m->next(), t = (m) ? ((T_MESH::Triangle*)m->data) : NULL)
					if (t->isExactlyDegenerate())
						flagDeg = false;

				// 若进一步检查没有问题，退出大循环；
				if (flagDeg)
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

		versMat.resize(0, 0);
		trisMat.resize(0, 0);
		TMesh2MeshMat(versMat, trisMat, mesh);
		debugWriteMesh("meshFixOutput_matRepr", versMat, trisMat);

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
	MESH_REPAIR::testCmd_meshDefectsDetect(argc, argv);

	std::cout << "main() finished." << std::endl;
}

