#include "myEigenModeling.h"
 
///////////////////////////////////////////////////////////////////////////////////////////////////// modeling接口：

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


#ifdef USE_TRIANGLE_H

//	重载1：封闭边界线点集得到面网格，可以是平面也可以是曲面，不在网格内部插点，三角片尺寸不可控。
template <typename DerivedO, typename DerivedC>
bool circuit2mesh(Eigen::PlainObjectBase<DerivedO>& vers, Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedC>& circVers)
{
	using ScalarO = typename DerivedO::Scalar;
	using RowVector3O = Eigen::Matrix<ScalarO, 1, 3>;
	using Matrix3O = Eigen::Matrix<ScalarO, 3, 3>;
	using MatrixXO = Eigen::Matrix<ScalarO, Eigen::Dynamic, Eigen::Dynamic>;

	unsigned circCount = circVers.rows();
	unsigned versCount = circVers.rows();

	// 0. 边缘环路顶点坐标数据拷入到输出网格中。
	vers = circVers.array().cast<ScalarO>();

	// 输入triangulate()的点集是投影到XOY平面的二维点，原点集应该旋转到合适的角度再投影。

	// 1. 取第一个点、1/3处的点、2/3处的点的所在平面的法向量作为原点集的法向量
	RowVector3O vers1 = circVers.row(0).array().cast<ScalarO>();
	RowVector3O vers2 = circVers.row(versCount / 3).array().cast<ScalarO>();
	RowVector3O vers3 = circVers.row(2 * versCount / 3).array().cast<ScalarO>();
	RowVector3O norm = (vers1 - vers2).cross(vers3 - vers2);
	norm.normalize();

	//  2. 旋转点集使得norm平行于z轴
	Matrix3O rotation = getRotationMat(norm, RowVector3O{ 0, 0, 1 });
	vers = (vers * rotation.transpose()).eval();

	// 3. 旋转后的点数据写入到triangulate()接口的输入结构体中。
	Eigen::MatrixXf vers2D;
	Eigen::MatrixXi edges2D;
	vers2D = vers.transpose().topRows(2).cast<float>();
	edges2D.resize(2, versCount);
	for (unsigned i = 0; i < versCount; i++)
	{
		edges2D(0, i) = i + 1;
		edges2D(1, i) = (i + 1) % versCount + 1;
	}

	triangulateio triIn;
	triangulateio triOut;
	triIn.numberofpoints = versCount;
	triIn.pointlist = vers2D.data();
	triIn.numberofpointattributes = 0;
	triIn.pointattributelist = NULL;
	triIn.pointmarkerlist = NULL;

	triIn.numberofsegments = versCount;
	triIn.segmentlist = (int*)edges2D.data();
	triIn.segmentmarkerlist = NULL;

	triIn.numberoftriangles = 0;
	triIn.numberofcorners = 0;
	triIn.numberoftriangleattributes = 0;

	triIn.numberofholes = 0;
	triIn.holelist = NULL;

	triIn.numberofregions = 0;
	triIn.regionlist = NULL;
	memset(&triOut, 0, sizeof(triangulateio));

	// 4. 执行二维三角剖分，得到输出网格三角片数据，不取顶点坐标数据，顶点坐标数据使用旋转操作前的。
	char triStr[256] = "pY";
	triangulate(triStr, &triIn, &triOut, NULL);
	tris.resize(3, triOut.numberoftriangles);
	std::memcpy(tris.data(), triOut.trianglelist, sizeof(int) * 3 * triOut.numberoftriangles);
	tris.transposeInPlace();
	tris.array() -= 1;

	return true;
}


// 重载2：会在回路内部插点，三角片面积上限由参数maxTriArea决定；
template <typename DerivedV, typename DerivedL, typename DerivedN>
bool circuit2mesh(Eigen::PlainObjectBase<DerivedV>& vers, Eigen::MatrixXi& tris, \
	const Eigen::PlainObjectBase<DerivedL>& versLoop, \
	const Eigen::PlainObjectBase<DerivedN>& normDir, const float maxTriArea)
{
	using ScalarV = typename DerivedV::Scalar;
	using ScalarN = typename DerivedN::Scalar;
	const int versLoopCount = versLoop.rows();

	// 1. 求一个仿射变换：
	Eigen::RowVector3f loopCenter = versLoop.colwise().mean().array().cast<float>();
	Eigen::Matrix4f rotationHomo{ Eigen::Matrix4f::Identity() };
	Eigen::Matrix4f translationHomo{ Eigen::Matrix4f::Identity() };
	Eigen::Matrix4f affineHomo, affineInvHomo;
	Eigen::Matrix3f rotation = getRotationMat(Eigen::RowVector3f{ 0, 0, 1 }, Eigen::RowVector3f{ normDir.array().cast<float>() });
	rotationHomo.topLeftCorner(3, 3) = rotation;
	translationHomo(0, 3) = loopCenter(0);
	translationHomo(1, 3) = loopCenter(1);
	translationHomo(2, 3) = loopCenter(2);
	affineHomo = translationHomo * rotationHomo;
	affineInvHomo = affineHomo.inverse();

	// 1. 插点的三角剖分：
	Eigen::MatrixXf versLoopHomo0, versLoopHomo, versLoop0, vers2D;
	Eigen::MatrixXi edges2D(2, versLoopCount);
	vers2HomoVers(versLoopHomo, versLoop);
	versLoopHomo0 = affineInvHomo * versLoopHomo;
	homoVers2Vers(versLoop0, versLoopHomo0);
	vers2D = versLoop0.transpose().topRows(2);
	for (unsigned i = 0; i < versLoopCount; i++)
	{
		edges2D(0, i) = i + 1;
		edges2D(1, i) = (i + 1) % versLoopCount + 1;
	}

	triangulateio triIn;
	triangulateio triOut;
	triIn.numberofpoints = versLoopCount;
	triIn.pointlist = vers2D.data();
	triIn.numberofpointattributes = 0;
	triIn.pointattributelist = NULL;
	triIn.pointmarkerlist = NULL;

	triIn.numberofsegments = versLoopCount;
	triIn.segmentlist = (int*)edges2D.data();
	triIn.segmentmarkerlist = NULL;

	triIn.numberoftriangles = 0;
	triIn.numberofcorners = 0;
	triIn.numberoftriangleattributes = 0;

	triIn.numberofholes = 0;
	triIn.holelist = NULL;

	triIn.numberofregions = 0;
	triIn.regionlist = NULL;
	memset(&triOut, 0, sizeof(triangulateio));

	char triStr[256];
	//		"p" - 折线段; "a" - 最大三角片面积; "q" - Quality; "Y" - 边界上不插入点;
	sprintf_s(triStr, "pq30.0a%fYY", maxTriArea);
	triangulate(triStr, &triIn, &triOut, NULL);

	// 2. 处理三角剖分后的结果；
	Eigen::MatrixXf versTmp, versTmpHomo;
	const int versCount = triOut.numberofpoints;						// 插点后的网格的点数； 
	versTmp.resize(versCount, 3);
	versTmp.setZero();
	for (size_t i = 0; i < versCount; i++)
	{
		versTmp(i, 0) = triOut.pointlist[i * 2];
		versTmp(i, 1) = triOut.pointlist[i * 2 + 1];
	}
	vers2HomoVers(versTmpHomo, versTmp);
	versTmpHomo = (affineHomo * versTmpHomo).eval();
	homoVers2Vers(versTmp, versTmpHomo);
	versTmp.topRows(versLoopCount) = versLoop.array().cast<float>();
	tris.resize(3, triOut.numberoftriangles);
	std::memcpy(tris.data(), triOut.trianglelist, sizeof(int) * 3 * triOut.numberoftriangles);
	tris.transposeInPlace();
	tris.array() -= 1;

	// for debug——打印光顺前的网格；
#if 0
	{
		objWriteMeshMat("E:/meshOutNoFaring.obj", versTmp, tris);
	}
#endif

	// 3. 此时网格内部是个平面，边界和输入边界线相同，需要在保持边界线的情形下光顺网格： 
	std::vector<int> fixedVerIdxes(versLoopCount);
	for (int i = 0; i < versLoopCount; ++i)
		fixedVerIdxes[i] = i;
	if (!laplaceFaring(vers, versTmp, tris, 0.1, 50, fixedVerIdxes))
		return false;

	return true;
}



#endif

// 模板特化输出：
#include "templateSpecialization.cpp"