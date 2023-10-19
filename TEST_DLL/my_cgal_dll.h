#ifndef _MY_CGAL_DLL_H
#define _MY_CGAL_DLL_H

#include <SDKDDKVer.h>	
#include <windows.h>					// Windows ͷ�ļ�: 

#include <iostream>


// ʹ��eigen���е������������ʾ���㡢���ơ���Ƭ����Ϣ��
#include "Eigen/Dense"


#ifdef MY_CGAL_DLL_EXPORTS
#define MY_CGAL_DLL_API __declspec(dllexport) 
#else
#define MY_CGAL_DLL_API __declspec(dllimport)
#endif


// ����������ע��������������ʵ�ֿ�ͷҲҪ���Ϻ�MY_CGAL_DLL_API��
template <typename DerivedVo, typename DerivedVi>
MY_CGAL_DLL_API bool alphaShapes2D(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	const Eigen::PlainObjectBase<DerivedVi>& versIn, const int alphaValue);



#endif

