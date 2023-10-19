#include "myEigenHomoCoor.h"


// 模板函数需要特化之后才能在静态库中输出： 

/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化：
template void objWriteHomoMeshMat<double>(const char*, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&, \
	const Eigen::MatrixXi& tris);
template void objWriteHomoMeshMat<float>(const char*, const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>&, \
	const Eigen::MatrixXi& tris);

template Eigen::MatrixXf vers2homoVersF<Eigen::Matrix<double, -1, -1, 0, -1, -1>>(\
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>& vers);
template Eigen::MatrixXf vers2homoVersF<Eigen::Matrix<float, -1, -1, 0, -1, -1>>(\
	const Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>& vers);

template Eigen::MatrixXd vers2homoVersD<Eigen::Matrix<double, -1, -1, 0, -1, -1>>(\
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>& vers);
template Eigen::MatrixXd vers2homoVersD<Eigen::Matrix<float, -1, -1, 0, -1, -1>>(\
	const Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>& vers);