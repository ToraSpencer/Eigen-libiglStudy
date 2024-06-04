#include "test_libigl.h" 

#include<stdio.h>
#include<assert.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <typeinfo>

#include <windows.h>
#include <atlstr.h>					// 包含CString类。属于microsoft ATL(活动模板库avtive template library)
#include <atlconv.h>
#include <io.h>
#include <winuser.h>


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



// 读取inputMESH文件夹中的.mesh文件，在viewer中展示：
int testCmd_showMESHfile(int argc, char** argv)
{  
	// 生成路径：
	int   nPos;
	CString   cPath;
	GetModuleFileName(NULL, cPath.GetBufferSetLength(MAX_PATH + 1), MAX_PATH);		// 获取当前进程加载的模块的路径。
	nPos = cPath.ReverseFind('\\');
	cPath = cPath.Left(nPos);
	std::string path{ CT2CA{cPath} };
	std::string pathMESH = path + "\\inputMESH"; 
	CString fileConfig = cPath + "\\config.ini";

	// 1. 读取部件网格  
	std::cout << "读取输入网格..." << std::endl;
	std::vector<std::string> fileNames, tmpStrVec, MESHfileNames;
	std::vector<Eigen::MatrixXd> meshesVers;
	std::vector<Eigen::MatrixXi> meshesTris;
	std::vector<Eigen::MatrixXi> meshesTets;

	getFileNames(pathMESH.c_str(), tmpStrVec, false);
	const unsigned meshesCount = tmpStrVec.size();

	if (0 == meshesCount)
	{
		debugDisp("输入文件夹为空。");
		return -1;
	}

	if (meshesCount > 1)
	{
		debugDisp("只能输入一个.mesh文件");
		return -1;
	}

	MESHfileNames.reserve(meshesCount);
	for (const auto& str : tmpStrVec)
	{
		std::string tailStr = str.substr(str.size() - 5, 5);				//	".mesh"
		if (".mesh" == tailStr)
		{
			fileNames.push_back(str);
			unsigned index = str.find_last_of("/");
			std::string MESHfileName = str.substr(index, str.size() - index - 5);			// "/" + obj文件名，不含路径和.obj后缀；
			MESHfileNames.push_back(MESHfileName);
		}
	}
	meshesVers.resize(meshesCount);
	meshesTris.resize(meshesCount);
	meshesTets.resize(meshesCount); 
	for (unsigned i = 0; i < meshesCount; ++i)
		igl::readMESH(fileNames[i].c_str(), meshesVers[i], meshesTets[i], meshesTris[i]); 
 

	// 2. 使用viewer渲染网格：

	// 1. 窗口中装载数据；
	viewer.data().set_mesh(meshesVers[0], meshesTris[0]);

	// 2. 设定三轴旋转；默认是两轴旋转，set_rotation_type()方法可以指定其他旋转风格
	viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

	// 3. show_lines指定是否画出网格线；
	viewer.data().show_lines = 0;

	// 4. 进入渲染循环：
	viewer.launch();
 
	debugDisp("finished.");

	return 0;
}

  

int main(int argc, char** argv)
{
	// IGL_DEFORMATION::test0();

	// IGL_MODELLING::test1();

	// IGL_GRAPH::test1();

	// testCmd_showMESHfile(argc, argv);

	std::cout << "main() finished." << std::endl;
}

