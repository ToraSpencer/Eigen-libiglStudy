#include "myEigenIO.h"

// 模板函数需要特化之后才能在静态库中输出： 


/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化： 
template void vecWriteToFile<int>(const char* fileName, const std::vector<int>& vec);
template void vecWriteToFile<float>(const char* fileName, const std::vector<float>& vec);
template void vecWriteToFile<double>(const char* fileName, const std::vector<double>& vec);

template void vecReadFromFile<int>(std::vector<int>& vec, const char* fileName, const unsigned elemCount);
template void vecReadFromFile<float>(std::vector<float>& vec, const char* fileName, const unsigned elemCount);
template void vecReadFromFile<double>(std::vector<double>& vec, const char* fileName, const unsigned elemCount);

template bool matWriteToFile<int>(const char* fileName, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& mat);
template bool matWriteToFile<float>(const char* fileName, const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& mat);
template bool matWriteToFile<double>(const char* fileName, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mat);

template bool matReadFromFile<int>(Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName);
template bool matReadFromFile<float>(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName);
template bool matReadFromFile<double>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mat, const char* fileName);

template void objReadMeshMat<float, int>(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::Matrix<int, Eigen::Dynamic, \
	Eigen::Dynamic>& tris, const char* fileName);
template void objReadMeshMat<double, int>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::Matrix<int, Eigen::Dynamic, \
	Eigen::Dynamic>& tris, const char* fileName);

template void objReadVerticesMat<float>(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, const char* fileName);
template void objReadVerticesMat<double>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, const char* fileName);
