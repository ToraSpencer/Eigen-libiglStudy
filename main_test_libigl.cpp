#include "test_libigl.h" 

#include<stdio.h>
#include<assert.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <typeinfo>

#include <windows.h>
#include <atlstr.h>					// ����CString�ࡣ����microsoft ATL(�ģ���avtive template library)
#include <atlconv.h>
#include <io.h>
#include <winuser.h>


static igl::opengl::glfw::Viewer viewer;				// ȫ�ֵ�����鿴������
static std::mutex g_mutex;						// ȫ�ֵĻ�������


// ��ǰ����-easy
/*

*/


// ��ǰ����-hard
/*
	1. �������������������㣬�ҳ�����֮������·����

	2. ʹ��OPENGLʵ����������⣬�����������ȷ�������ϵ�������������Ϊ���롣


*/


// ��Ŀ��Ϣ
/*
	���뻷���� x64 Relase

	ʹ�õĵ�������:
		eigen
		libigl
		glad
		glfw
*/


// ��Ŀ�м���Ԥ����꣺
/*
	IGL_STATIC_LIBRARY

	NOMINMAX
			���windows.h��<algorithm>��ͬʱ����std::max()��std::min()�ĳ�ͻ��

	TUTORIAL_SHARED_PATH="G:/gitRepositories/ligIgl/libigl_CGAL_openGL/cmake/../external/../tutorial/data"

	CMAKE_INTDIR="Release"
*/


////////////////////////////////////////////////////////////////////////////// ����WINAPI��һЩ�ӿڣ�
namespace MY_WIN_API
{
	// ��ȡĳ��Ŀ¼�������ļ�����Ŀ¼����
	void getFileNames(std::string path, std::vector<std::string>& files, bool blRecur = true)
	{
		std::string str;
		struct _finddata_t fileinfo;			// �ļ���Ϣ����_finddata_t����ͷ�ļ�<io.h>
		intptr_t hFile = _findfirst(str.assign(path).append("/*").c_str(), &fileinfo);							// �ļ����	
		bool blFileValid = (hFile != -1);

		if (blFileValid)
		{
			do
			{
				bool isSubDir = (fileinfo.attrib & _A_SUBDIR);
				//�����Ŀ¼,�ݹ���ң��������,���ļ�����·������vector��
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


	// config�ļ�·���ַ��������� �� std::wstring�ַ�����
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


	// ��ȡconfig�ļ��е�ĳһ������Ϊfloat�ļ�ֵ
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


	// ��ȡconfig�ļ��е�ĳһ������Ϊint�ļ�ֵ
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



////////////////////////////////////////////////////////////////////////////////////////////// DEBUG �ӿ�
namespace MY_DEBUG
{
	static std::string g_debugPath = "E:/";

	// lambda������ӡstd::cout֧�ֵ����ͱ�����
	template <typename T>
	static auto disp = [](const T& arg)
	{
		std::cout << arg << ", ";
	};

	static void debugDisp()			// �ݹ���ֹ
	{						//		�ݹ���ֹ��Ϊ�޲λ���һ�����������ζ����ԡ�
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



// ��ȡinputMESH�ļ����е�.mesh�ļ�����viewer��չʾ��
int testCmd_showMESHfile(int argc, char** argv)
{  
	// ����·����
	int   nPos;
	CString   cPath;
	GetModuleFileName(NULL, cPath.GetBufferSetLength(MAX_PATH + 1), MAX_PATH);		// ��ȡ��ǰ���̼��ص�ģ���·����
	nPos = cPath.ReverseFind('\\');
	cPath = cPath.Left(nPos);
	std::string path{ CT2CA{cPath} };
	std::string pathMESH = path + "\\inputMESH"; 
	CString fileConfig = cPath + "\\config.ini";

	// 1. ��ȡ��������  
	std::cout << "��ȡ��������..." << std::endl;
	std::vector<std::string> fileNames, tmpStrVec, MESHfileNames;
	std::vector<Eigen::MatrixXd> meshesVers;
	std::vector<Eigen::MatrixXi> meshesTris;
	std::vector<Eigen::MatrixXi> meshesTets;

	getFileNames(pathMESH.c_str(), tmpStrVec, false);
	const unsigned meshesCount = tmpStrVec.size();

	if (0 == meshesCount)
	{
		debugDisp("�����ļ���Ϊ�ա�");
		return -1;
	}

	if (meshesCount > 1)
	{
		debugDisp("ֻ������һ��.mesh�ļ�");
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
			std::string MESHfileName = str.substr(index, str.size() - index - 5);			// "/" + obj�ļ���������·����.obj��׺��
			MESHfileNames.push_back(MESHfileName);
		}
	}
	meshesVers.resize(meshesCount);
	meshesTris.resize(meshesCount);
	meshesTets.resize(meshesCount); 
	for (unsigned i = 0; i < meshesCount; ++i)
		igl::readMESH(fileNames[i].c_str(), meshesVers[i], meshesTets[i], meshesTris[i]); 
 

	// 2. ʹ��viewer��Ⱦ����

	// 1. ������װ�����ݣ�
	viewer.data().set_mesh(meshesVers[0], meshesTris[0]);

	// 2. �趨������ת��Ĭ����������ת��set_rotation_type()��������ָ��������ת���
	viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

	// 3. show_linesָ���Ƿ񻭳������ߣ�
	viewer.data().show_lines = 0;

	// 4. ������Ⱦѭ����
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

