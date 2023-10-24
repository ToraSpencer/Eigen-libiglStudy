#include "myEigenModeling.h"

// 模板函数需要特化之后才能在静态库中输出： 


/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化：



template void genCubeMesh(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);
template void genCubeMesh(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);
 
template void genAABBmesh(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<float, 3>& aabb);
template void genAABBmesh(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<double, 3>& aabb);
 
template bool genGrids(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, \
	const Eigen::Matrix<float, 1, 3>& origin, const float step, const std::vector<unsigned>& gridCounts);
template bool genGrids(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, \
	const Eigen::Matrix<double, 1, 3>& origin, const float step, const std::vector<unsigned>& gridCounts);
template bool genGrids(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, \
	const Eigen::Matrix<float, 1, 3>& origin, const float step, const std::vector<unsigned>& gridCounts);
template bool genGrids(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, \
	const Eigen::Matrix<double, 1, 3>& origin, const float step, const std::vector<unsigned>& gridCounts);