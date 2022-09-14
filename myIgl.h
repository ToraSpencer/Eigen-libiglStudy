#pragma once

#include "myEigen.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <functional>

#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/readSTL.h>
#include <igl/writeSTL.h>

#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/readDMAT.h>
#include <igl/repdiag.h>
#include <igl/opengl/glfw/Viewer.h>
 
#include "tutorial_shared_path.h"


// class——方向包围盒类OBB（本质上是中点在原点的AABB加上一个仿射变换）；
template <typename _Scalar>
class OBB : public AlignedBox<_Scalar, 3>
{
	// 成员数据：
public:
	Eigen::Matrix<_Scalar, 1, 3> m_dir;
	Eigen::Matrix<_Scalar, 1, 3> m_center;
	Eigen::Matrix3d m_rotation;

	// constructor:
public:
	inline OBB() :m_dir(0, 0, 1), m_center(Eigen::Matrix<_Scalar, 1, 3>::Zero(3)), m_rotation(Eigen::Matrix3d::Identity()) {}

	inline OBB(const AlignedBox<_Scalar, 3>& aabb)
	{
		auto minp = aabb.m_min;
		auto maxp = aabb.m_max;
		this->center = (minp + maxp) / 2.0;
		this->m_dir = Eigen::Matrix<_Scalar, 1, 3>(0, 0, 1);
		this->m_rotation = Eigen::Matrix3d::Identity();

		this->m_min = minp - this->center;
		this->m_max = maxp - this->center;
	}

	inline OBB(const AlignedBox<_Scalar, 3>& aabb, const Eigen::Matrix<_Scalar, 1, 3>& dir, \
		const Eigen::Matrix<_Scalar, 1, 3>& center): m_dir(dir), m_center(center)
	{
		const double epsilon = 1e-6;
		assert((dir.norm() - 1) < epsilon && "input dir is invalid.");

		auto minp = aabb.min();
		auto maxp = aabb.max();
		Eigen::Matrix<_Scalar, 1, 3> centerOri = (minp + maxp) / 2.0;
		this->m_min = minp - centerOri;
		this->m_max = maxp - centerOri;

		Eigen::Matrix<_Scalar, 1, 3> crossValue = Eigen::Matrix<_Scalar, 1, 3>(0, 0, 1).cross(dir);		// 旋转轴；
		double dotValue = Eigen::Matrix<_Scalar, 1, 3>(0, 0, 1).dot(dir);

		// 若指定方向和z轴平行：
		if (crossValue.norm() < epsilon)
		{
			this->m_rotation = Eigen::Matrix3d::Identity();
			if (dotValue < 0)
				this->m_rotation(2, 2) = -1;
			return;
		}

		Eigen::Matrix<_Scalar, 1, 3> zNew = dir;
		Eigen::Matrix<_Scalar, 1, 3> xNew = -crossValue.normalized();
		Eigen::Matrix<_Scalar, 1, 3> yNew = zNew.cross(xNew);

		for (int i = 0; i < 3; ++i)
		{
			this->m_rotation(i, 0) = static_cast<double>(xNew(i));
			this->m_rotation(i, 1) = static_cast<double>(yNew(i));
			this->m_rotation(i, 2) = static_cast<double>(zNew(i));
		}
	}

};


// 生成中心在原点，边长为1，三角片数为12的正方体网格；
template	<typename DerivedV, typename DerivedI>
void genCubeMesh(Eigen::Matrix<DerivedV, Dynamic, Dynamic>& vers, Eigen::Matrix<DerivedI, Dynamic, Dynamic>& tris) 
{
	vers.resize(8, 3);
	vers << -0.5000000, -0.5000000, -0.5000000, -0.5000000, 0.5000000, -0.5000000, \
		0.5000000, -0.5000000, -0.5000000, 0.5000000, 0.5000000, -0.5000000, \
		0.5000000, -0.5000000, 0.5000000, 0.5000000, 0.5000000, 0.5000000, \
		-0.5000000, -0.5000000, 0.5000000, -0.5000000 ,0.5000000, 0.5000000;

	tris.resize(12 ,3);
	tris << 1,2,0, 1,3,2, 3,4,2, 3,5,4, 0,4,6, 0,2,4, 7,3,1, 7,5,3, 7,0,6, 7,1,0, 5,6,4, 5,7,6;
}

 
// 生成轴向包围盒的三角网格；
template <typename _Scalar, int _AmbientDim>
void genAABBmesh(const Eigen::AlignedBox<_Scalar, _AmbientDim>& aabb, Eigen::Matrix<_Scalar, Dynamic, Dynamic>& vers, \
		Eigen::MatrixXi& tris)
{
	Matrix<_Scalar, _AmbientDim, 1> minp = aabb.min();
	Matrix<_Scalar, _AmbientDim, 1> maxp = aabb.max();
	Matrix<_Scalar, _AmbientDim, 1> newOri = (minp + maxp) / 2.0;
	Matrix<_Scalar, _AmbientDim, 1> sizeVec = maxp - minp;
	
	genCubeMesh(vers, tris);
	vers.col(0) *= sizeVec(0);
	vers.col(1) *= sizeVec(1);
	vers.col(2) *= sizeVec(2);
	vers.rowwise() += newOri.transpose();
}


template <typename _Scalar>
void genOBBmesh(const OBB<_Scalar>& obb, Eigen::Matrix<_Scalar, Dynamic, Dynamic>& vers, \
	Eigen::MatrixXi& tris)
{
	Matrix<_Scalar, 3, 1> minp = obb.min();
	Matrix<_Scalar, 3, 1> maxp = obb.max();
	Matrix<_Scalar, 3, 1> sizeVec = maxp - minp;

	genCubeMesh(vers, tris);
	vers.col(0) *= sizeVec(0);
	vers.col(1) *= sizeVec(1);
	vers.col(2) *= sizeVec(2);

	vers = ((obb.m_rotation * vers.transpose()).eval()).transpose().eval();
	vers.rowwise() += obb.m_center;
}

