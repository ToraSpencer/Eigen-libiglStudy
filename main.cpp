#include <windows.h>
#include <atlstr.h>					// 包含CString类。属于microsoft ATL(活动模板库avtive template library)
#include <atlconv.h>
#include <io.h>
#include <winuser.h>

#include "triMesh.h"
#include "representations.h"

#include "test_dense_mat.h"
#include "test_sparse_mat.h"
#include "test_scientific_calc.h" 
#include "test_myEigenModeling.h"
#include "test_imgui.h"
 
 
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
	已禁用vcpkg
	  
	使用的Eigen库版本为3.4.0
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
namespace MY_WIN_API 
{
	// 读取某个目录下所有文件名、目录名；
	void getFileNames(std::string path, std::vector<std::string>& files, bool blRecur = true)
	{
		std::string str;
		struct _finddata_t fileinfo;			// 文件信息。宏_finddata_t来自头文件<io.h>
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
}
using namespace MY_WIN_API;



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




////////////////////////////////////////////////////////////////////////////// TEST: 不同的几何表象 
namespace TEST_REPRESENTATIONS
{
	void test1() 
	{
		triMeshF mesh;
		readOBJ(mesh, "E:/材料/tooth.obj");

		// 1. 元素索引
		std::vector<verF>& vers = mesh.vertices;
		debugDisp("vers[0][0] == ", vers[0][0]);

		vers[0][0] = INFINITY;
		debugDisp("vers[0][0] == ", vers[0][0]);

		// debugDisp("vers[0][3] == ", vers[0][3]);			// 会抛出std::out_of_range异常；

		// 2. 简单的运算
		verF ver1 = vers[1] - vers[3];
		verF ver2 = vers[2] - vers[3];
		const verF& ver3 = vers[3];
		debugDisp("ver1 == ", ver1);
		debugDisp("ver2 == ", ver2);
		debugDisp("ver3 == ", ver3);
		debugDisp("ver2 - ver1 == ", ver2 - ver1);
		debugDisp("ver3 - ver1 == ", ver3 - ver1);
		debugDisp("0.5 * ver1 == ", 0.5 * ver1);
		debugDisp("ver1 * 0.5 == ", ver1 * 0.5);
		debugDisp("ver1 / 0.5 == ", ver1 / 0.5);
		debugDisp("\n\n");

		ver1 *= 2;
		debugDisp("ver1 *= 2; ver1 == ", ver1);

		ver2 = ver2 + ver3; 
		debugDisp("ver2 = ver2 + ver3; ver2 == ", ver2);

		ver1 = ver1 - ver3;
		auto tmp = -ver3;
		debugDisp("ver1 = ver1 - ver3; ver1 == ", ver1);
		debugDisp("tmp == ", tmp);

		// 3. cast()方法
		auto retCast = ver1.cast<int>();
		debugDisp("typeid(ver1).name() == ", typeid(ver1).name());
		debugDisp("typeid(retCast).name() == ", typeid(retCast).name());

		debugDisp("finished.");
	}

	// test IO:
	void test2() 
	{
		triMeshF meshF;
		triMeshD meshD;
		readOBJ(meshF, "E:/材料/tooth.obj");

		meshD.triangles = meshF.triangles;
		meshD.vertices.reserve(meshF.vertices.size());
		for (const auto& ver : meshF.vertices)
			meshD.vertices.push_back(ver.cast<double>());
		writeSTL("E:/toothF.stl", meshF);
		writeSTL("E:/toothD.stl", meshD);

		debugDisp("test2() finished.");
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

			// 检测非流形点
			std::vector<int> nmnVerIdxes;
			int nmnVersCount = nonManifoldVers(nmnVerIdxes, vers, tris);
			if (nmnVersCount > 0)
			{
				debugDisp(OBJfileNames[i], ".obj ！！！存在非流形点，nmnVersCount == ", nmnVersCount);
				Eigen::MatrixXd nmnVers;
				subFromIdxCon(nmnVers, vers, nmnVerIdxes);
				ss.str("");
				ss << OBJfileNames[i] << "_nmnVers";
				debugWriteVers(ss.str().c_str(), nmnVers);
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
			int genus = 0;
			int ueCount = 0;
			Eigen::Index minEdgeIdx = 0;
			double minLen = 0;
			Eigen::MatrixXd edgeArrows, normals;
			Eigen::MatrixXi edges, uEdges;
			Eigen::MatrixXi bdrys, bdryTris;
			Eigen::VectorXd edgesLen;
			std::vector<int> bdryTriIdxes;
			getEdges(edges, tris);
			getEdgeArrows(edgeArrows, edges, vers);
			getUedges(uEdges, edges);
			ueCount = uEdges.rows();
			genus = (ueCount - versCount - trisCount) / 2 + 1;				//  g == (e - v- f )/2+1;
			if (genus > 0)
				debugDisp(OBJfileNames[i], ".obj genus == ", genus);

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
			Eigen::MatrixXi nmnEdges;
			std::vector<std::pair<std::pair<int, int>, int>> nmnInfos;
			int nmnCount = nonManifoldEdges(nmnEdges, nmnInfos, tris);
			if (nmnCount < 0)
			{
				debugDisp("error! nonManifoldUEs() run failed.");
				return -1;
			}

			// for debug;
			std::unordered_map<std::int64_t, int> nmnMap;
			for (const auto& pair: nmnInfos) 
			{
				auto edgePair = pair.first;
				std::int64_t code = encodeUedge(edgePair.first, edgePair.second);
				auto iter = nmnMap.find(code);
				if (nmnMap.end() == iter)
					nmnMap.insert({ code, pair.second });
				else
					nmnMap[code] += pair.second;
			}
			
			std::vector<std::int64_t> sickCodes;
			for (const auto& pair : nmnMap)
				if (3 == pair.second)
					sickCodes.push_back(pair.first);
			int sickCount = sickCodes.size();
			Eigen::MatrixXi sickNMNedges(sickCount, 2);
			for (int k = 0; k < sickCount; ++k) 
			{
				auto pair = decodeEdge(sickCodes[k]);
				sickNMNedges(k, 0) = pair.first;
				sickNMNedges(k, 1) = pair.second;
			}
			
			debugWriteEdges("sickNMNedges", sickNMNedges, vers); 
 

			if (nmnEdges.rows() > 0)
			{
				Eigen::MatrixXd	nmnVers;
				std::unordered_set<int> nmnVerIdxes;
				int* ptrInt = nullptr;
				debugDisp(OBJfileNames[i], ".obj ！！！存在非流形边，nmnEdges.rows() == ", nmnEdges.rows());

				ss.str("");
				ss << OBJfileNames[i] << "_nmnEdges";
				debugWriteEdges(ss.str().c_str(), nmnEdges, vers);
				ptrInt = nmnEdges.data();
				for (int i = 0; i < nmnEdges.size(); ++i)
					nmnVerIdxes.insert(*ptrInt++);
				subFromIdxCon(nmnVers, vers, nmnVerIdxes);

				ss.str("");
				ss << OBJfileNames[i] << "_nmnVers";
				debugWriteVers(ss.str().c_str(), nmnVers);
			}
 
#if 0
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
#endif

			// f6. 检测退化三角片：
			Eigen::MatrixXd triNorms;
			if (!getTriNorms(triNorms, vers, tris))
				debugDisp("error!!! degenerate tris detected.");

			// f7. 提取所有单连通区域：
			int scCount = 0;							// 顶点单连通区域个数；
			int sctCount = 0;							// 三角片单连通区域个数；
			Eigen::SparseMatrix<int> adjSM, adjSM_eCount, adjSM_eIdx;
			Eigen::VectorXi connectedLabels, connectedCount, connectedTriLabels, connectedTriCount;
			Eigen::VectorXi connectedLabelsSCT, connectedCountSCT;
			tris2adjMat(adjSM_eCount, adjSM_eIdx, tris);
			adjSM = adjSM_eCount;
			traverseSparseMatrix(adjSM, [&adjSM](auto& iter)
				{
					if (iter.value() > 0)
						iter.valueRef() = 1;
				});
			scCount = simplyConnectedRegion(connectedLabels, connectedCount, adjSM);
			if (scCount > 1)
				debugDisp(OBJfileNames[i], ".obj ！！！存在", scCount, "个顶点联通区域。");

			sctCount = simplyTrisConnectedRegion(connectedLabelsSCT, connectedCountSCT, tris);
			if (scCount != sctCount)
			{
				debugDisp(OBJfileNames[i], ".obj ！！！三角片连通区域数目大于顶点联通区域数目。");
				debugDisp(OBJfileNames[i], ".obj scCount == ", scCount);
				debugDisp(OBJfileNames[i], ".obj sctCount == ", sctCount);
			}
			if (sctCount > 0)
			{
				// 打印各个联通区域的三角片：
				for (int i = 0; i < sctCount; ++i)
				{
					Eigen::MatrixXi selectedTris;
					Eigen::VectorXi flags = (i == connectedLabelsSCT.array()).select(\
						Eigen::VectorXi::Ones(trisCount), Eigen::VectorXi::Zero(trisCount));
					subFromFlagVec(selectedTris, tris, eigenVec2Vec(flags));

					ss.str("");
					ss << "meshSCT_" << i;
					debugWriteMesh(ss.str().c_str(), vers, selectedTris);
				}
			}
			
		}

		debugDisp("finished.");
		getchar();

		return 0;
	}


	// 三角片生长提取最外层网格：
	int testCmd_triangleGrowOutterSurf(int argc, char** argv)
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
			Eigen::MatrixXd& vers = meshesVers[i];
			Eigen::MatrixXi& tris = meshesTris[i];
			Eigen::MatrixXd versOut;
			Eigen::MatrixXi trisOut; 
			unsigned versCount = vers.rows();
			unsigned trisCount = tris.rows();

			// f1. 确定所有非流形无向边：
			Eigen::MatrixXi nmnEdges;
			int nmnCount = nonManifoldEdges(nmnEdges, tris);
			if (nmnCount < 0)
			{
				debugDisp("error!!!", OBJfileNames[i],  ".obj nonManifoldUEs() run failed.");
				return -1;
			}
			if (nmnEdges.rows() > 0)
			{
				debugDisp(OBJfileNames[i], ".obj 非流形无向边数量：nmnEdges.rows() == ", nmnEdges.rows());
				ss.str("");
				ss << OBJfileNames[i] << "_nmnEdgesOrigin";
				debugWriteEdges(ss.str().c_str(), nmnEdges, vers);
			}

			// f2. 一次提取外表面可能不能完全去除非流形边，需要多次调用：
			while (nmnCount > 0)
			{ 
				if (!triangleGrowOuterSurf(versOut, trisOut, vers, tris, true))
					debugDisp("error!!!", OBJfileNames[i], ".obj triangleGrowOuterSurf() failed!"); 

				// 检查输出网格中是否有非流形边：
				nmnEdges.resize(0, 0);
				nmnCount = nonManifoldEdges(nmnEdges, trisOut);
				if (nmnCount < 0)
				{
					debugDisp("error!!!", OBJfileNames[i], ".obj nonManifoldUEs() run failed.");
					return -1;
				}
				if (nmnEdges.rows() > 0)
					debugWriteEdges("nmnEdges", nmnEdges, versOut);
				nmnCount = nmnEdges.rows();

				vers = versOut;
				tris = trisOut;
				trisCount = tris.rows();
			}

			// f3. 输出结果：
			ss.str("");
			ss << pathOutput << OBJfileNames[i] << "_noNMNedges.obj";
			objWriteMeshMat(ss.str().c_str(), vers, tris);
			if (debugFlag)
			{
				ss.str("");
				ss << OBJfileNames[i] << "_noNMNedges";
				debugWriteMesh(ss.str().c_str(), vers, tris);
			}
		}

		debugDisp("finished.");
		getchar();

		return 0;
	}


	// 对inputOBJ文件夹中的网格不插点补洞
	int testCmd_fillSmallHoles(int argc, char** argv)
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

		// 2. 网格找洞、补洞的循环：
		for (int i = 0; i < meshesCount; ++i)
		{ 
			Eigen::MatrixXd& vers = meshesVers[i];
			Eigen::MatrixXi& tris = meshesTris[i];
			std::vector<std::vector<int>> holes;
			
			// f1. 找洞
			int holesCount = findHoles(holes, vers, tris);
			if (holesCount <= 0)
			{
				std::string tmpStr = (0 == holesCount) ? (OBJfileNames[i] + std::string{".obj 没有洞，不需要补。"}) : \
					(std::string{ "error!!! " } + OBJfileNames[i] + std::string{ ".obj调用findHoles()失败。" });
				debugDisp(tmpStr);
				continue;
			}

			// f2. 补洞
			Eigen::MatrixXi newTris;
			if (!fillSmallHoles(newTris, holes))
			{
				debugDisp("error!!! ", OBJfileNames[i], ".obj调用fillSmallHoles()失败。");
				continue;
			}
			else
				debugDisp(OBJfileNames[i], ".obj补了", holesCount, "个洞。");
			matInsertRows(tris, newTris);
			 
			// f3. 输出结果：
			ss.str("");
			ss << pathOutput << OBJfileNames[i] << "_holeFilled.obj";
			objWriteMeshMat(ss.str().c_str(), vers, tris);
			if (debugFlag)
			{
				ss.str("");
				ss << OBJfileNames[i] << "_holeFilled";
				debugWriteMesh(ss.str().c_str(), vers, tris);
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
	 

	// 去除网格退化边（边长过短的边，将两个顶点融合）;
	void test11() 
	{
		Eigen::MatrixXd vers, versOut, edgeArrows;
		Eigen::MatrixXi tris, trisOut, edges;
		Eigen::VectorXi selectedIdxes, oldNewIdxInfo;
		objReadMeshMat(vers, tris, "E:/材料/meshDegTris.obj"); 
		debugWriteMesh("meshInput", vers, tris);
 
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
		objReadMeshMat(vers,tris,"E:/材料/meshIsoVers.obj"); 

		std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, tris);
		if (isoVerIdxes.size() > 0)
			std::cout << isoVerIdxes.size() << " isolated vertices detected." << std::endl;

		removeIsoVers(versOut, trisOut, vers, tris, isoVerIdxes);
		debugWriteMesh("meshOut", versOut, trisOut);

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
			matInsertRows(trisCopy, Eigen::RowVector3i{vcIdx, vaIdx, vxIdx});
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
			subFromFlagVec(trisOut, tris, eigenVec2Vec(triFlags));
			tris = trisOut;
			debugWriteMesh("noOpTris", vers, tris);
		}
		else
			debugDisp("输入网格没有重叠三角片。");

		// 2. 打印结果网格中的洞、边缘曲线：


		debugDisp("finished.");
	}


	// 删除包含短边的三角片；
	void test6()
	{
		const double lenThreshold = 1e-5;
		Eigen::MatrixXd vers, arrows, vas, vbs, triEdgeLens;
		Eigen::MatrixXi tris, edges;
		Eigen::VectorXi vaIdxes, vbIdxes;
		Eigen::VectorXd edgeLens;
		objReadMeshMat(vers, tris, "E:/材料/tmpTriangleGrowOuterSurf.obj");
		int versCount = vers.rows();
		int trisCount = tris.rows();

		// 生成三角片边长信息：
		if (!getEdges(edges, tris))
		{
			debugDisp("error!!! getEdges() failed.");
			return;
		}
		vaIdxes = edges.col(0);
		vbIdxes = edges.col(1);
		subFromIdxVec(vas, vers, vaIdxes);
		subFromIdxVec(vbs, vers, vbIdxes);
		arrows = vbs - vas;
		edgeLens = arrows.rowwise().norm();
		triEdgeLens = Eigen::Map<Eigen::MatrixXd>(edgeLens.data(), trisCount, 3);

		// 搜索边长都小于阈值的三角片：
		std::vector<int> sickTriIdxes;
		Eigen::MatrixXi sickTris;
		int sickCount = 0;
		sickTriIdxes.reserve(trisCount);
		for (int i = 0; i < trisCount; ++i)
			if (triEdgeLens(i, 0) < lenThreshold || triEdgeLens(i, 1) < lenThreshold || triEdgeLens(i, 2) < lenThreshold)
				sickTriIdxes.push_back(i);
		sickCount = sickTriIdxes.size();
		sickTriIdxes.shrink_to_fit();
		subFromIdxVec(sickTris, tris, sickTriIdxes);
		debugDisp("sickCount == ", sickCount);

		// 去除被标记的三角片：
		Eigen::VectorXi flags = Eigen::VectorXi::Ones(trisCount);
		for (const auto& index : sickTriIdxes)
			flags(index) = 0;
		Eigen::MatrixXi tmpTris;
		subFromFlagVec(tmpTris, tris, eigenVec2Vec(flags));
		tris = tmpTris;
			
		// 
		Eigen::MatrixXd versOut;
		Eigen::MatrixXi trisOut;
		if (removeSickDupTris(vers, tris) < 0)
		{
			debugDisp("error!!! removeSickDupTris() failed");
			return;
		}
		std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, tris);
		removeIsoVers(versOut, trisOut, vers, tris, isoVerIdxes);
		debugWriteMesh("meshNoSmallTris", versOut, trisOut);

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

}



////////////////////////////////////////////////////////////////////////////// 生成控制台程序工具：
namespace TEST_CMD 
{
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


	//
	template<typename DerivedVp, typename DerivedV>
	double calcSolidAngle(const Eigen::MatrixBase<DerivedVp>& pos, \
		const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
	{
		const int versCount = vers.rows();
		const int trisCount = tris.rows();

		// 1. 计算每个三角片的重心： 
		Eigen::MatrixXd barys(trisCount, 3), normals(trisCount, 3);
		Eigen::Matrix3d triVers;
		Eigen::RowVector3d va, vb, vc, normDir;
		Eigen::VectorXd areas(trisCount);
		for (int i = 0; i < trisCount; ++i)
		{
			va = vers.row(tris(i, 0)).array().cast<double>();
			vb = vers.row(tris(i, 1)).array().cast<double>();
			vc = vers.row(tris(i, 2)).array().cast<double>();
			triVers.row(0) = va;
			triVers.row(1) = vb;
			triVers.row(2) = vc;
			barys.row(i) = triVers.colwise().mean();
			normDir = (vc - vb).cross(va - vb);
			double crossNorm = normDir.norm();
			areas(i) = 0.5 * crossNorm;
			normals.row(i) = normDir.normalized();
		}

		// 2. 
		Eigen::RowVector3d posD = pos.array().cast<double>();
		Eigen::MatrixXd arrows = barys.rowwise() - posD;
		Eigen::VectorXd weights(trisCount), arrowLenDb(trisCount);
		Eigen::VectorXd resultVec(trisCount);
		for (int i = 0; i < trisCount; ++i)
		{
			Eigen::RowVector3d arr = arrows.row(i);
			arrowLenDb(i) = arr.dot(arr);
			arr.normalize();
			Eigen::RowVector3d norm = normals.row(i);
			weights(i) = arr.dot(norm);
			// resultVec(i) = weights(i) * areas(i) / arrowLenDb(i);
		}

		resultVec = weights.array() * areas.array() / arrowLenDb.array();

		return resultVec.sum();
	}


}


// 生成满秩的随机矩阵
void test1() 
{



}


int main(int argc, char** argv)
{  
	TEST_IMGUI::test2();
	
	//TEST_PMP::test2(); 

	// TEST_MYEIGEN_MODELING::test4();

	// TEST_MYEIGEN_IO::test0();

	// TEST_REPRESENTATIONS::test2();
	
	// TEST_SPARSE_MAT::test6();
	
	// TEST_DENSE_MAT::test3();

	debugDisp("main() finished."); 
}

