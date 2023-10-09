#include "myEigenPMP.h"
 

/////////////////////////////////////////////////////////////////////////////////////////////////// 表象转换接口：

template <typename DerivedI>
void edges2mat(Eigen::PlainObjectBase<DerivedI>& mat, const std::vector<std::pair<int, int>>& edges)
{
	mat.resize(edges.size(), 2);
	for (unsigned i = 0; i < edges.size(); ++i)
	{
		mat(i, 0) = edges[i].first;
		mat(i, 1) = edges[i].second;
	}
}


// 生成边编码——边两个端点的索引对映射成一个64位整型数
template <typename Index>
std::int64_t encodeEdge(const Index vaIdx, const Index vbIdx)
{
	std::int64_t a = static_cast<std::int64_t>(vaIdx);
	std::int64_t b = static_cast<std::int64_t>(vbIdx);
	std::int64_t code = 0;
	code |= (a << 32);
	code |= b;
	return code;
}


// 生成无向边编码——边两个端点的索引对整理成前小后大的顺序，映射成一个64位整型数
template <typename Index>
std::int64_t encodeUedge(const Index vaIdx, const Index vbIdx)
{
	Index vaIdx0 = (vaIdx < vbIdx ? vaIdx : vbIdx);
	Index vbIdx0 = (vaIdx < vbIdx ? vbIdx : vaIdx);
	std::int64_t a = static_cast<std::int64_t>(vaIdx0);
	std::int64_t b = static_cast<std::int64_t>(vbIdx0);
	std::int64_t code = 0;
	code |= (a << 32);
	code |= b;
	return code;
}


// 生成三角片编码——三个顶点索引排序后映射成64位无符号整型数
template <typename Index>
std::uint64_t encodeTriangle(const Index vaIdx, const Index vbIdx, const Index vcIdx)
{
	std::uint64_t triIdxLimit = 0x1FFFFF;								// 最多为21位全1, 0x1FFFFF ==  2097181, 两百多万三角片；
	std::uint64_t vaIdx0 = static_cast<std::uint64_t>(vaIdx);
	std::uint64_t vbIdx0 = static_cast<std::uint64_t>(vbIdx);
	std::uint64_t vcIdx0 = static_cast<std::uint64_t>(vcIdx);

	assert(vaIdx0 <= triIdxLimit && vbIdx0 <= triIdxLimit && vcIdx0 <= triIdxLimit, "triangle index out of range !!!");
	assert(vaIdx0 != vbIdx0 && vaIdx0 != vcIdx0 && vbIdx0 != vcIdx0, "invalid triangle index !!!");

	// 区分正反面——即(a, b, c)和(a, c, b)会映射成不同的编码；但是(a,b,c)(b,c,a)(c,a,b)映射成相同的编码；
	std::uint64_t A, B, C;
	if (vaIdx0 < vbIdx0 && vaIdx0 < vcIdx0)
	{
		A = vaIdx0; B = vbIdx0; C = vcIdx0;						// abc
	}
	else if (vbIdx0 < vaIdx0 && vbIdx0 < vcIdx0)
	{
		A = vbIdx0; B = vcIdx0; C = vaIdx0;						// bca
	}
	else
	{
		A = vcIdx0; B = vaIdx0; C = vbIdx0;						// cab;
	}
	std::uint64_t code = 0;
	code |= (A << 42);
	code |= (B << 21);
	code |= C;
	return code;
}


// 解码边编码；
std::pair<int, int> decodeEdge(const std::int64_t code)
{
	int a = static_cast<int>(code >> 32);
	int b = static_cast<int>(code - (static_cast<std::int64_t>(a) << 32));
	return std::make_pair(a, b);
}
 

// 解码三角片编码：
std::vector<int> decodeTrianagle(const std::uint64_t code)
{
	std::uint64_t a = code >> 42;
	std::uint64_t resi = code - (a << 42);
	std::uint64_t b = resi >> 21;
	std::uint64_t c = resi - (b << 21);
	return std::vector<int>{static_cast<int>(a), static_cast<int>(b), static_cast<int>(c)};
}



/////////////////////////////////////////////////////////////////////////////////////////////////// 图形属性：
 
// 得到三角网格的有向边数据
template <typename DerivedI>
bool getEdges(Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedI>& tris)
{
	/*
		bool getEdges(
				Eigen::MatrixXi& edges,
				const Eigen::PlainObjectBase<DerivedI>& tris
				)

		三条边在edges中的排列顺序为[bc; ca; ab]；其所对的顶点标记corner分别为0, 1, 2，即a,b,c;
		边索引到三角片索引的映射——边eIdx0所在的三角片triIdx0 == eIdx0 % trisCount;
		三角片triIdx0和corner0到边索引的映射——eIdx0 = corner0 * trisCount + triIdx0;
		边索引到corner的映射——corner0 = eIdx0 / trisCount;
	*/
	const unsigned trisCount = tris.rows();
	const unsigned edgesCount = 3 * trisCount;
	const unsigned versCount = tris.maxCoeff() + 1;

	edges = Eigen::MatrixXi::Zero(edgesCount, 2);
	Eigen::MatrixXi vaIdxes = tris.col(0).array().cast<int>();
	Eigen::MatrixXi vbIdxes = tris.col(1).array().cast<int>();
	Eigen::MatrixXi vcIdxes = tris.col(2).array().cast<int>();
	edges.block(0, 0, trisCount, 1) = vbIdxes;
	edges.block(trisCount, 0, trisCount, 1) = vcIdxes;
	edges.block(trisCount * 2, 0, trisCount, 1) = vaIdxes;
	edges.block(0, 1, trisCount, 1) = vcIdxes;
	edges.block(trisCount, 1, trisCount, 1) = vaIdxes;
	edges.block(trisCount * 2, 1, trisCount, 1) = vbIdxes;

	return true;
}


// 生成环路的边数据；
template <typename IndexType, typename DerivedI>
bool getLoopEdges(Eigen::PlainObjectBase<DerivedI>& edges, const IndexType versCount)
{
	using Index = typename DerivedI::Scalar;							// 使用从属名称，需要加上"typename"前缀；
	const int versCount0 = static_cast<int>(versCount);
	if (versCount0 < 2)
		return false;

	edges.resize(0, 0);
	edges.resize(versCount0, 2);

	edges.col(0) = Eigen::Matrix<Index, Eigen::Dynamic, 1>::LinSpaced(versCount, 0, versCount0 - 1);
	edges.col(1) = Eigen::Matrix<Index, Eigen::Dynamic, 1>::LinSpaced(versCount, 1, versCount0);
	edges(versCount0 - 1, 1) = 0;

	return true;
}


// 模板特化输出：
#include "templateSpecialization.cpp"