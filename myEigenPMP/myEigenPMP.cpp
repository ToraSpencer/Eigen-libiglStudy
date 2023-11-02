#include "myEigenPMP.h"
 

/////////////////////////////////////////////////////////////////////////////////////////////////// 表象转换接口：

template <typename DerivedI>
void edges2mat(Eigen::PlainObjectBase<DerivedI>& mat, \
	const std::vector<std::pair<int, int>>& edges)
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


/////////////////////////////////////////////////////////////////////////////////////////////////// 生成基础图元：

// 输入起点、终点、空间采样率，插值生成一条直线点云；
template <typename DerivedVo, typename ScalarVi>
bool interpolateToLine(Eigen::PlainObjectBase<DerivedVo>& vers, \
	const Eigen::Matrix<ScalarVi, 1, 3>& start, const Eigen::Matrix<ScalarVi, 1, 3>& end, \
	const float SR, const bool SE)
{
	if (vers.rows() > 0)
		return false;

	Eigen::RowVector3d startD = start.array().cast<double>();
	Eigen::RowVector3d endD = end.array().cast<double>();
	Eigen::RowVector3d dir = endD - startD;
	float length = dir.norm();
	dir.normalize();
	if (length <= SR)
		return true;

	if (SE)
		matInsertRows(vers, startD);

	float lenth0 = 0;
	for (unsigned i = 1; (length - lenth0) > 0.8 * SR; i++)			// 确保最后一个点距离终点不能太近。
	{
		Eigen::RowVector3d temp = startD + SR * i * dir;
		matInsertRows(vers, temp);
		lenth0 = SR * (i + 1);		// 下一个temp的长度。
	}

	if (SE)
		matInsertRows(vers, endD);

	return true;
};


// 生成XOY平面内的圆圈点集：
template <typename T>
bool getCircleVers(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const float radius, \
	const unsigned versCount)
{
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	double deltaTheta = 2 * pi / versCount;
	vers.resize(versCount, 3);
	vers.setZero();
	for (unsigned i = 0; i < versCount; ++i)
	{
		double theta = deltaTheta * i;
		vers(i, 0) = radius * cos(theta);
		vers(i, 1) = radius * sin(theta);
	}

	return true;
}


// 对环路点集进行不插点三角剖分——手动缝合方法； 
bool circuitGetTris(Eigen::MatrixXi& tris, const std::vector<int>& indexes, const bool regularTri)
{
	/*
		bool circuitGetTris(
			Eigen::MatrixXi& tris								三角剖分生成的triangle soup
			const std::vector<int>& indexes,			单个环路点集，顶点需要有序排列。
			const bool regularTri								true——右手螺旋方向为生成面片方向；false——反向；
			)

	*/

	int versCount = static_cast<int>(indexes.size());
	if (versCount < 3)
		return false;						// invalid input

	// lambda——环路中的顶点索引转换；
	auto getIndex = [versCount](const int index0) -> int
	{
		int index = index0;
		if (index >= versCount)
			while (index >= versCount)
				index -= versCount;
		else if (index < 0)
			while (index < 0)
				index += versCount;
		return index;
	};

	tris.resize(versCount - 2, 3);
	int count = 0;
	int k = 0;
	if (regularTri)
	{
		while (count < versCount - 2)
		{
			tris.row(count++) = Eigen::RowVector3i{ indexes[getIndex(-k - 1)], indexes[getIndex(-k)], indexes[getIndex(k + 1)] };
			if (count == versCount - 2)
				break;

			tris.row(count++) = Eigen::RowVector3i{ indexes[getIndex(-k - 1)], indexes[getIndex(k + 1)], indexes[getIndex(k + 2)] };
			k++;
		}
	}
	else
	{
		while (count < versCount - 2)
		{
			tris.row(count++) = Eigen::RowVector3i{ indexes[getIndex(-k - 1)], indexes[getIndex(k + 1)], indexes[getIndex(-k)] };
			if (count == versCount - 2)
				break;

			tris.row(count++) = Eigen::RowVector3i{ indexes[getIndex(-k - 1)], indexes[getIndex(k + 2)], indexes[getIndex(k + 1)] };
			k++;
		}
	}

	return true;
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
template <typename DerivedI, typename IndexType>
bool getLoopEdges(Eigen::PlainObjectBase<DerivedI>& edges, \
	const IndexType versCount)
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


// 生成余切权的laplacian (同时也是刚度矩阵(stiffness matrix))——基于论文：Polygon Laplacian Made Simple [Bunge et al. 2020]
template<typename Tl, typename DerivedV>
bool cotLaplacian(Eigen::SparseMatrix<Tl>& L, \
	const Eigen::PlainObjectBase<DerivedV>& vers, \
	const Eigen::MatrixXi& tris)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();
	L.resize(versCount, versCount);
	L.reserve(10 * versCount);

	Eigen::MatrixXi edges(3, 2);
	edges << 1, 2, 2, 0, 0, 1;

	// Gather cotangents
	Eigen::Matrix<Tl, Eigen::Dynamic, Eigen::Dynamic> C;
	trisCotValues(C, vers, tris);

	std::vector<Eigen::Triplet<Tl> > IJV;
	IJV.reserve(tris.rows() * edges.rows() * 4);

	// Loop over triangles
	for (int i = 0; i < tris.rows(); i++)
	{
		// loop over edges of element
		for (int e = 0; e < edges.rows(); e++)
		{
			int source = tris(i, edges(e, 0));
			int dest = tris(i, edges(e, 1));
			IJV.push_back(Eigen::Triplet<Tl>(source, dest, C(i, e)));
			IJV.push_back(Eigen::Triplet<Tl>(dest, source, C(i, e)));
			IJV.push_back(Eigen::Triplet<Tl>(source, source, -C(i, e)));
			IJV.push_back(Eigen::Triplet<Tl>(dest, dest, -C(i, e)));
		}
	}
	L.setFromTriplets(IJV.begin(), IJV.end());

	return true;
}


/////////////////////////////////////////////////////////////////////////////////////////////////// 图形编辑：

// laplace光顺 
template <typename DerivedVo, typename DerivedVi>
bool laplaceFaring(Eigen::PlainObjectBase<DerivedVo>& versOut, \
	const Eigen::PlainObjectBase<DerivedVi>& vers, \
	const Eigen::MatrixXi& tris, const float deltaLB, const unsigned loopCount, \
	const std::vector<int>& fixedVerIdxes)
{
	assert(3 == vers.cols() && "assert!!! Input mesh should be in 3-dimension space.");

	using ScalarO = typename DerivedVo::Scalar;
	const int versCount = vers.rows();
	Eigen::SparseMatrix<ScalarO> I(versCount, versCount);
	I.setIdentity();
	versOut.resize(0, 0); 

	// 1. 计算laplace矩阵：
	Eigen::SparseMatrix<ScalarO> Lmat;
	Eigen::Matrix<ScalarO, Eigen::Dynamic, Eigen::Dynamic> versCopy = vers.array().cast<ScalarO>();
	cotLaplacian(Lmat, versCopy, tris);

	// 
	Eigen::Matrix<ScalarO, Eigen::Dynamic, Eigen::Dynamic> versFixed;
	if (!fixedVerIdxes.empty())
	{
		const int versFixedCount = fixedVerIdxes.size();
		int idxNew = 0;
		versFixed.resize(versFixedCount, 3);
		for (const auto& index : fixedVerIdxes)
			versFixed.row(idxNew++) = versCopy.row(index);
	}

	// 2. 光顺的循环：
	for (unsigned i = 0; i < loopCount; ++i)
	{
		// 2.1 ???
		versOut = (I + deltaLB * Lmat) * versCopy;

		// 2.2 固定不动的顶点拷贝回去：
		if (!fixedVerIdxes.empty())
		{
			int idxNew = 0;
			for (const auto& index : fixedVerIdxes)
				versOut.row(index) = versFixed.row(idxNew++);
		}

		// 
		versCopy = versOut;
	}

	return true;
}


/////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化输出：
#include "templateSpecialization.cpp"