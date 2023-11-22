#include "myEigenModeling.h"

// ģ�庯����Ҫ�ػ�֮������ھ�̬��������� 


/////////////////////////////////////////////////////////////////////////////////////////////////////////// ģ���ػ���

template void genCubeMesh(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);
template void genCubeMesh(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);
 
template void genAABBmesh(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<float, 3>& aabb);
template void genAABBmesh(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<double, 3>& aabb);
 