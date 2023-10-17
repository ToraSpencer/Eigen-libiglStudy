#include "myEigenTools.h"
 
// µ¼Èë¾²Ì¬¿â£º
#include "myEigenIO/myEigenIO.h"
#pragma comment(lib, "myEigenIO.lib")


DLL_API void readOBJ(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& tris, const char* fileName)
{
	objReadMeshMat(vers, tris, fileName);
}

DLL_API void readOBJ(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& tris, const char* fileName)
{
	objReadMeshMat(vers, tris, fileName);
}
 

