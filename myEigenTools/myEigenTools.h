#ifndef _DYNAMIC_MY_EIGEN_TOOLS_H
#define _DYNAMIC_MY_EIGEN_TOOLS_H

// 包括 SDKDDKVer.h 将定义可用的最高版本的 Windows 平台。
/*
		如果要为以前的 Windows 平台生成应用程序，请包括 WinSDKVer.h，并将
		将 _WIN32_WINNT 宏设置为要支持的平台，然后再包括 SDKDDKVer.h。
*/
#include <SDKDDKVer.h>	

#define WIN32_LEAN_AND_MEAN        // 从 Windows 头文件中排除极少使用的内容
#include <windows.h>							// Windows 头文件

#include "Eigen/Dense"
#include "Eigen/Sparse"


// 宏MYEIGENTOOLS_EXPORTS
/*
	本头文件在生成动态库、使用动态库的项目中都要使用
	生成动态库的项目中应该定义预编译宏MYEIGENTOOLS_EXPORTS，从而执行输出动态库的分支。
		 方法1：项目属性->C/C++ -> 预处理器 -> 预处理器定义中加入宏MYEIGENTOOLS_EXPORTS
		 方法2：stdafx.h中#define MYEIGENTOOLS_EXPORTS
	使用动态库的项目中执行导入动态库的分支，使用__declspec(dllimport)
*/
#ifdef MYEIGENTOOLS_EXPORTS
#define DLL_API __declspec(dllexport) 
#else
#define DLL_API __declspec(dllimport)
#endif


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 导出函数
DLL_API void readOBJ(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::Matrix<int, Eigen::Dynamic, 	Eigen::Dynamic>& tris, const char* fileName);
DLL_API void readOBJ(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& tris, const char* fileName);
 
 
#endif
 
