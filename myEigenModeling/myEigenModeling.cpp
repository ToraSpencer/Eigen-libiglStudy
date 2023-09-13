#include "myEigenModeling.h"

static const double pi = 3.14159265359;
 

///////////////////////////////////////////////////////////////////////////////////////////////////// modeling接口：

// 输入起点、终点、空间采样率，插值生成一条直线点云；
template <typename T>
bool interpolateToLine(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::Matrix<T, 1, 3>& start, \
	const Eigen::Matrix<T, 1, 3>& end, const float SR, const bool SE)
{
	if (vers.rows() > 0)
		return false;

	Eigen::Matrix<T, 1, 3> dir = end - start;
	float length = dir.norm();
	dir.normalize();

	if (length <= SR)
		return true;

	if (SE)
		matInsertRows(vers, start);

	float lenth0 = 0;
	for (unsigned i = 1; (length - lenth0) > 0.8 * SR; i++)			// 确保最后一个点距离终点不能太近。
	{
		Eigen::Matrix<T, 1, 3> temp = start + SR * i * dir;
		matInsertRows(vers, temp);
		lenth0 = SR * (i + 1);		// 下一个temp的长度。
	}

	if (SE)
		matInsertRows(vers, end);

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
template <typename IndexType>
bool circuitGetTris(Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic>& tris, \
	const std::vector<int>& indexes, const bool regularTri)
{
	/*
		bool circuitGetTris(
			Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic>& tris,		三角剖分生成的triangle soup
			const std::vector<int>& indexes,				单个环路点集，顶点需要有序排列。
			const bool regularTri																				true——右手螺旋方向为生成面片方向；false——反向；
			)

	*/
	using MatrixXI = Eigen::Matrix<IndexType, Eigen::Dynamic, Eigen::Dynamic>;
	using RowVector3I = Eigen::Matrix<IndexType, 1, 3>;

	int versCount = static_cast<int>(indexes.size());
	if (versCount < 3)
		return false;						// invalid input

	// lambda——环路中的顶点索引转换；
	auto getIndex = [versCount](const int index0) -> IndexType
	{
		int index = index0;
		if (index >= versCount)
			while (index >= versCount)
				index -= versCount;
		else if (index < 0)
			while (index < 0)
				index += versCount;
		return static_cast<IndexType>(index);
	};

	tris.resize(versCount - 2, 3);
	int count = 0;
	int k = 0;
	if (regularTri)
	{
		while (count < versCount - 2)
		{
			tris.row(count++) = RowVector3I{ indexes[getIndex(-k - 1)], indexes[getIndex(-k)], indexes[getIndex(k + 1)] };
			if (count == versCount - 2)
				break;

			tris.row(count++) = RowVector3I{ indexes[getIndex(-k - 1)], indexes[getIndex(k + 1)], indexes[getIndex(k + 2)] };
			k++;
		}
	}
	else
	{
		while (count < versCount - 2)
		{
			tris.row(count++) = RowVector3I{ indexes[getIndex(-k - 1)], indexes[getIndex(k + 1)], indexes[getIndex(-k)] };
			if (count == versCount - 2)
				break;

			tris.row(count++) = RowVector3I{ indexes[getIndex(-k - 1)], indexes[getIndex(k + 2)], indexes[getIndex(k + 1)] };
			k++;
		}
	}

	return true;
}


// 生成中心在原点，边长为1，三角片数为12的正方体网格； 
template	<typename T>
void genCubeMesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris)
{
	vers.resize(0, 0);
	tris.resize(0, 0);
	vers.resize(8, 3);
	vers << -0.5000000, -0.5000000, -0.5000000, -0.5000000, 0.5000000, -0.5000000, \
		0.5000000, -0.5000000, -0.5000000, 0.5000000, 0.5000000, -0.5000000, \
		0.5000000, -0.5000000, 0.5000000, 0.5000000, 0.5000000, 0.5000000, \
		- 0.5000000, -0.5000000, 0.5000000, -0.5000000, 0.5000000, 0.5000000;

	tris.resize(12, 3);
	tris << 1, 2, 0, 1, 3, 2, 3, 4, 2, 3, 5, 4, 0, 4, 6, 0, 2, 4, 7, 3, 1, 7, 5, 3, 7, 0, 6, 7, 1, 0, 5, 6, 4, 5, 7, 6;
}


// 生成轴向包围盒的三角网格；
template <typename T>
void genAABBmesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<T, 3>& aabb)
{
	vers.resize(0, 0);
	tris.resize(0, 0);
	Eigen::Matrix<T, 3, 1> minp = aabb.min();
	Eigen::Matrix<T, 3, 1> maxp = aabb.max();
	Eigen::Matrix<T, 3, 1> newOri = (minp + maxp) / 2.0;
	Eigen::Matrix<T, 3, 1> sizeVec = maxp - minp;

	genCubeMesh(vers, tris);
	vers.col(0) *= sizeVec(0);
	vers.col(1) *= sizeVec(1);
	vers.col(2) *= sizeVec(2);
	vers.rowwise() += newOri.transpose();
}
 

// 生成栅格采样点云
template<typename Tg, typename To>
bool genGrids(Eigen::Matrix<Tg, Eigen::Dynamic, Eigen::Dynamic>& gridCenters, const Eigen::Matrix<To, 1, 3>& origin, \
	const float step, const std::vector<unsigned>& gridCounts)
{

	/*
		bool genGrids(
			Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& gridCenters,			// 输出的栅格点云
			const Eigen::Matrix<T, 1, 3>& origin,														// 栅格原点，即三轴坐标最小的那个栅格点；
			const float step,																						// 采样步长
			const std::vector<unsigned>& gridCounts											// 三元数组，分别是XYZ三轴上的步数；
			)
	*/

	using VectorXT = Eigen::Matrix<Tg, Eigen::Dynamic, 1>;
	using RowVector3T = Eigen::Matrix<Tg, 1, 3>;
	using MatrixXT = Eigen::Matrix<Tg, Eigen::Dynamic, Eigen::Dynamic>;	

	// 整个距离场的包围盒：
	RowVector3T minp = origin.array().cast<Tg>();
	RowVector3T maxp = minp + static_cast<Tg>(step) * RowVector3T(gridCounts[0] - 1, gridCounts[1] - 1, gridCounts[2] - 1);


	// 生成栅格：
	/*
		按索引增大排列的栅格中心点为：
		gc(000), gc(100), gc(200), gc(300),...... gc(010), gc(110), gc(210), gc(310),...... gc(001), gc(101), gc(201).....

		x坐标：
		x0, x1, x2, x3......x0, x1, x2, x3......x0, x1, x2, x3......
		周期为xCount;
		重复次数为(yCount * zCount)

		y坐标：
		y0, y0, y0...y1, y1, y1...y2, y2, y2.........y0, y0, y0...y1, y1, y1...
		周期为(xCount * yCount);
		重复次数为zCount;
		单个元素重复次数为xCount

		z坐标：
		z0, z0, z0......z1, z1, z1......z2, z2, z2......
		单个元素重复次数为(xCount * yCount)
	*/
	VectorXT xPeriod = VectorXT::LinSpaced(gridCounts[0], minp(0), maxp(0));
	VectorXT yPeriod = VectorXT::LinSpaced(gridCounts[1], minp(1), maxp(1));
	VectorXT zPeriod = VectorXT::LinSpaced(gridCounts[2], minp(2), maxp(2));

	MatrixXT tmpVec0, tmpVec1, tmpVec2, tmpVec11;
	//tmpVec0 = kron(VectorXT::Ones(gridCounts[1] * gridCounts[2]), xPeriod);
	//tmpVec11 = kron(yPeriod, VectorXT::Ones(gridCounts[0]));
	//tmpVec1 = kron(VectorXT::Ones(gridCounts[2]), tmpVec11);
	//tmpVec2 = kron(zPeriod, VectorXT::Ones(gridCounts[0] * gridCounts[1]));

	kron(tmpVec0, VectorXT::Ones(gridCounts[1] * gridCounts[2]), xPeriod);
	kron(tmpVec11, yPeriod, VectorXT::Ones(gridCounts[0]));
	kron(tmpVec1, VectorXT::Ones(gridCounts[2]), tmpVec11);
	kron(tmpVec2, zPeriod, VectorXT::Ones(gridCounts[0] * gridCounts[1]));
 
	gridCenters.resize(gridCounts[0] * gridCounts[1] * gridCounts[2], 3);
	gridCenters.col(0) = tmpVec0;
	gridCenters.col(1) = tmpVec1;
	gridCenters.col(2) = tmpVec2;

 
	return true;
}


// 模板特化输出：
#include "templateSpecialization.cpp"