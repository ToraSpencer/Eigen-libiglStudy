#ifndef _MY_CGAL_DLL_H
#define _MY_CGAL_DLL_H

#include <SDKDDKVer.h>	
#include <windows.h>					// Windows 头文件: 

#include <iostream>


// 使用eigen库中的向量、矩阵表示顶点、点云、面片等信息；
#include "Eigen/Dense"


#ifdef MY_CGAL_DLL_EXPORTS
#define MY_CGAL_DLL_API __declspec(dllexport) 
#else
#define MY_CGAL_DLL_API __declspec(dllimport)
#endif


// 导出函数；注！！！——函数实现开头也要加上宏MY_CGAL_DLL_API：
template <typename DerivedVo, typename DerivedVi>
MY_CGAL_DLL_API bool alphaShapes2D(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	const Eigen::PlainObjectBase<DerivedVi>& versIn, const int alphaValue);



#endif

