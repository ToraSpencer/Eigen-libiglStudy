#include "myEigenPMP.h"

// 模板函数需要特化之后才能在静态库中输出： 


/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化：

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

template bool getEdges(Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<Eigen::MatrixXi>& tris);

template bool getLoopEdges(Eigen::PlainObjectBase<Eigen::MatrixXi>& edges, const int versCount);
template bool getLoopEdges(Eigen::PlainObjectBase<Eigen::MatrixXi>& edges, const unsigned versCount);
 