#include "myEigenPMP.h"

// 模板函数需要特化之后才能在静态库中输出： 


/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化：

template bool interpolateToLine(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::Matrix<float, 1, 3>& start, \
	const Eigen::Matrix<float, 1, 3>& end, const float SR, const bool SE);
template bool interpolateToLine(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::Matrix<double, 1, 3>& start, \
	const Eigen::Matrix<double, 1, 3>& end, const float SR, const bool SE);

template bool getCircleVers(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, const unsigned versCount);
template bool getCircleVers(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, const unsigned versCount);

template void genCubeMesh(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);
template void genCubeMesh(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);
 
template bool circuitGetTris(Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& tris, \
	const std::vector<int>& indexes, const bool regularTri);

template void genAABBmesh(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<float, 3>& aabb);
template void genAABBmesh(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<double, 3>& aabb);

template void edges2mat(Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, const std::vector<std::pair<int, int>>& edges);

template std::int64_t encodeEdge(const int vaIdx, const int vbIdx);
template std::int64_t encodeEdge(const unsigned int vaIdx, const unsigned int vbIdx);
template std::int64_t encodeEdge(const std::int64_t vaIdx, const std::int64_t vbIdx);

template std::int64_t encodeUedge(const int vaIdx, const int vbIdx);
template std::int64_t encodeUedge(const unsigned int vaIdx, const unsigned int vbIdx);
template std::int64_t encodeUedge(const std::int64_t vaIdx, const std::int64_t vbIdx);
 
template std::uint64_t encodeTriangle(const int vaIdx, const int vbIdx, const int vcIdx);
template std::uint64_t encodeTriangle(const unsigned int vaIdx, const unsigned int vbIdx, const unsigned int vcIdx);
template std::uint64_t encodeTriangle(const std::int64_t vaIdx, const std::int64_t vbIdx, const std::int64_t vcIdx);

 