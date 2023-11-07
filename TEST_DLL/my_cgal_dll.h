#ifndef _MY_CGAL_DLL_H
#define _MY_CGAL_DLL_H

#include <SDKDDKVer.h>	
#include <windows.h>					// Windows ͷ�ļ�:  

#include "Eigen/Dense"


#ifdef MY_CGAL_DLL_EXPORTS
#define MY_CGAL_DLL_API __declspec(dllexport) 
#else
#define MY_CGAL_DLL_API __declspec(dllimport)
#endif


// ����������ע��������������ʵ�ֿ�ͷҲҪ���Ϻ�MY_CGAL_DLL_API��


// 
/*
	bool alphaShapes2D(
			Eigen::PlainObjectBase<DerivedVo>& versOut,						����߽��ߵ���
			const Eigen::PlainObjectBase<DerivedVi>& versIn,				������ƣ�������XOYƽ���ڣ�
			const int alphaValue
			);

*/
template <typename DerivedVo, typename DerivedVi>
MY_CGAL_DLL_API bool alphaShapes2D(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	const Eigen::PlainObjectBase<DerivedVi>& versIn, const int alphaValue);


// ����ͬ�Ե������ػ��֣�
template <typename DerivedVo, typename DerivedVi>
MY_CGAL_DLL_API bool remeshSurfMesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::MatrixXi& trisOut, const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const Eigen::MatrixXi& trisIn, const double tarEdgeLen, const int iterCount = 3, const bool blKeepBdry = true);


// ����ƽ���������������˳��
template <typename DerivedVo, typename DerivedVi, typename STLcontainer>
MY_CGAL_DLL_API bool smoothSurfMesh(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	Eigen::MatrixXi& trisOut, const Eigen::PlainObjectBase<DerivedVi>& versIn, \
	const Eigen::MatrixXi& trisIn, const STLcontainer& versConstrained, \
	const double time = 0.0001, const unsigned iterCount = 10);


/////////////////////////////////////////////////////////////////////////////////////////////////////////// temp, for test:
namespace TEST_MYDLL
{
	MY_CGAL_DLL_API void testRemesh();
	MY_CGAL_DLL_API void testSmoothing();
}


#endif

