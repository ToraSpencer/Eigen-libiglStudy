#include "myEigenPMP.h"

// 模板函数需要特化之后才能在静态库中输出： 


/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化：

template void edges2mat<Eigen::MatrixXi>(Eigen::PlainObjectBase<Eigen::MatrixXi>&,\
	const std::vector<std::pair<int, int>>&);

template std::int64_t encodeEdge<int>(const int, const int);
template std::int64_t encodeEdge<unsigned>(const unsigned, const unsigned);
template std::int64_t encodeEdge<std::int64_t>(const std::int64_t, const std::int64_t);

template std::int64_t encodeUedge<int>(const int, const int);
template std::int64_t encodeUedge<unsigned>(const unsigned, const unsigned);
template std::int64_t encodeUedge<std::int64_t>(const std::int64_t, const std::int64_t);
 
template std::uint64_t encodeTriangle<int>(const int , const int, const int);
template std::uint64_t encodeTriangle<unsigned>(const unsigned, const unsigned, const unsigned);
template std::uint64_t encodeTriangle<std::int64_t>(const std::int64_t, const std::int64_t, const std::int64_t);

template bool interpolateToLine<Eigen::MatrixXf, float>(Eigen::PlainObjectBase<Eigen::MatrixXf>&, \
	const Eigen::Matrix<float, 1, 3>&, const Eigen::Matrix<float, 1, 3>&, const float, const bool);
template bool interpolateToLine<Eigen::MatrixXd, double>(Eigen::PlainObjectBase<Eigen::MatrixXd>&, \
	const Eigen::Matrix<double, 1, 3>&, const Eigen::Matrix<double, 1, 3>&, const float, const bool);
template bool interpolateToLine<Eigen::MatrixXd, float>(Eigen::PlainObjectBase<Eigen::MatrixXd>&, \
	const Eigen::Matrix<float, 1, 3>&, const Eigen::Matrix<float, 1, 3>&, const float, const bool);
template bool interpolateToLine<Eigen::MatrixXf, double>(Eigen::PlainObjectBase<Eigen::MatrixXf>&, \
	const Eigen::Matrix<double, 1, 3>&, const Eigen::Matrix<double, 1, 3>&, const float, const bool);

template bool getCircleVers(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, const unsigned versCount);
template bool getCircleVers(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, const unsigned versCount);

template bool getEdges<Eigen::MatrixXi>(Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXi>&);

template bool getLoopEdges<Eigen::MatrixXi, int>(Eigen::PlainObjectBase<Eigen::MatrixXi>&, const int);
template bool getLoopEdges<Eigen::MatrixXi, unsigned>(Eigen::PlainObjectBase<Eigen::MatrixXi>&, const unsigned);

template bool cotLaplacian<float, Eigen::MatrixXf>(Eigen::SparseMatrix<float>&, \
	const Eigen::PlainObjectBase<Eigen::MatrixXf>&, const Eigen::MatrixXi&);
template bool cotLaplacian<double, Eigen::MatrixXd>(Eigen::SparseMatrix<double>&, \
	const Eigen::PlainObjectBase<Eigen::MatrixXd>&, const Eigen::MatrixXi&);
template bool cotLaplacian<float, Eigen::MatrixXd>(Eigen::SparseMatrix<float>&, \
	const Eigen::PlainObjectBase<Eigen::MatrixXd>&, const Eigen::MatrixXi&);
template bool cotLaplacian<double, Eigen::MatrixXf>(Eigen::SparseMatrix<double>&, \
	const Eigen::PlainObjectBase<Eigen::MatrixXf>&, const Eigen::MatrixXi&);

template bool laplaceFaring<Eigen::MatrixXf, Eigen::MatrixXf>(\
	Eigen::PlainObjectBase<Eigen::MatrixXf>&, const Eigen::PlainObjectBase<Eigen::MatrixXf>&, \
	const Eigen::MatrixXi&, const float, const unsigned, const std::vector<int>&);
template bool laplaceFaring<Eigen::MatrixXd, Eigen::MatrixXd>(\
	Eigen::PlainObjectBase<Eigen::MatrixXd>&, const Eigen::PlainObjectBase<Eigen::MatrixXd>&, \
	const Eigen::MatrixXi&, const float, const unsigned, const std::vector<int>&);
template bool laplaceFaring<Eigen::MatrixXd, Eigen::MatrixXf>(\
	Eigen::PlainObjectBase<Eigen::MatrixXd>&, const Eigen::PlainObjectBase<Eigen::MatrixXf>&, \
	const Eigen::MatrixXi&, const float, const unsigned, const std::vector<int>&);
template bool laplaceFaring<Eigen::MatrixXf, Eigen::MatrixXd>(\
	Eigen::PlainObjectBase<Eigen::MatrixXf>&, const Eigen::PlainObjectBase<Eigen::MatrixXd>&, \
	const Eigen::MatrixXi&, const float, const unsigned, const std::vector<int>&);
 