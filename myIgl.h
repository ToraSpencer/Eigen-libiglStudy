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
	// 成员数据：注！！！继承自AlignedBox中的m_min, m_max是列向量，下面两个是行向量；
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

		this->m_min = minp - this->center.transpose();
		this->m_max = maxp - this->center.transpose();
	}


	inline OBB(const AlignedBox<_Scalar, 3>& aabb, const Eigen::Matrix<_Scalar, 1, 3>& dir, \
		const Eigen::Matrix<_Scalar, 1, 3>& center): m_dir(dir), m_center(center)
	{
		// 注：方向包围盒的方向取dir和-dir没有区别，所以控制dir和z轴正向夹角不大于pi/2;
		const double epsilon = 1e-6;
		assert((dir.norm() - 1) < epsilon && "input dir is invalid.");

		double dotValue = Eigen::Matrix<_Scalar, 1, 3>(0, 0, 1).dot(dir);
		if (dotValue < 0)
			this->m_dir = -dir;
		Eigen::Matrix<_Scalar, 1, 3> crossValue = Eigen::Matrix<_Scalar, 1, 3>(0, 0, 1).cross(this->m_dir);		// 旋转轴；

		auto minp = aabb.min();
		auto maxp = aabb.max();
		Eigen::Matrix<_Scalar, 1, 3> centerOri = (minp + maxp) / 2.0;
		this->m_min = minp - centerOri.transpose();
		this->m_max = maxp - centerOri.transpose();

		// 若指定方向和z轴平行：
		if (crossValue.norm() < epsilon)
		{
			this->m_rotation = Eigen::Matrix3d::Identity();
			return;
		}

		Eigen::Matrix<_Scalar, 1, 3> zNew = this->m_dir;
		Eigen::Matrix<_Scalar, 1, 3> xNew = -crossValue.normalized();
		Eigen::Matrix<_Scalar, 1, 3> yNew = zNew.cross(xNew);

		for (int i = 0; i < 3; ++i)
		{
			this->m_rotation(i, 0) = static_cast<double>(xNew(i));
			this->m_rotation(i, 1) = static_cast<double>(yNew(i));
			this->m_rotation(i, 2) = static_cast<double>(zNew(i));
		}
	}


	inline OBB(const AlignedBox<_Scalar, 3>& aabb, const Eigen::Matrix<_Scalar, 1, 3>& center, \
		const Eigen::Matrix3d& rotation) : m_center(center), m_rotation(rotation)
	{
		const double epsilon = 1e-6;
		auto minp = aabb.min();
		auto maxp = aabb.max();
		Eigen::Matrix<_Scalar, 1, 3> centerOri = (minp + maxp) / 2.0;
		this->m_min = minp - centerOri.transpose();
		this->m_max = maxp - centerOri.transpose();

		this->m_dir = (rotation * Eigen::Matrix<_Scalar, 3, 1>(0, 0, 1)).transpose(); 
	}

	
	//// 输入点云，生成方向包围盒，方向为点云的长轴：
	//inline OBB(const Eigen::Matrix<_Scalar, Dynamic, 3>& vers) 
	//{
	//	// 


	//}


	template<typename Derived>
	inline bool contains(const Eigen::MatrixBase<Derived>& ver) const
	{
		Eigen::Vector3d minp = this->m_min.cast<double>();
		Eigen::Vector3d maxp = this->m_max.cast<double>();

		const Derived& ver0 = ver.derived();
		assert((ver0.rows() == 1 && ver0.cols() == 3) || (ver0.cols() == 1 && ver0.rows() == 3) && "Input vertice must be a 3D vector!");

		// 顶点逆仿射变换到AABB的空间：
		auto dataPtr = ver0.data();
		Eigen::Vector3d p00(static_cast<double>(dataPtr[0]), static_cast<double>(dataPtr[1]), static_cast<double>(dataPtr[2]));
		Eigen::Vector3d p0 = this->m_rotation.inverse() * (p00 - this->m_center.transpose()); 
		return (minp.array() <= p0.array()).all() && (p0.array() <= maxp.array()).all();
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

