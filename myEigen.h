#pragma once
#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <fstream>
#include <queue>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
#include <initializer_list>
#include <memory>
#include <thread>
#include <limits>
#include <chrono>
#include <windows.h>

#include "Eigen/Dense"
#include "Eigen/Sparse"


/////////////////////////////////////////////////////////////////////////////////////////////////// 自己写的各个静态库：
#include "myEigenIO/myEigenIO.h"
#pragma comment(lib,"myEigenIO.lib")	

#include "myEigenBasicMath/myEigenBasicMath.h"
#pragma comment(lib, "myEigenBasicMath.lib")

#include "myEigenPMP/myEigenPMP.h"
#pragma comment(lib, "myEigenPMP.lib")

#include "myEigenModeling/myEigenModeling.h"
#pragma comment(lib, "myEigenModeling.lib")

 
/////////////////////////////////////////////////////////////////////////////////////////////////// debug tools

// 传入函数子或函数指针遍历stl容器
template<typename T, typename F>
void traverseSTL(T& con, F f)
{
	std::for_each(con.begin(), con.end(), f);
	std::cout << std::endl;
}


// 反向遍历
template<typename T, typename F>
void revTraverseSTL(T& con, F f)
{
	std::for_each(con.rbegin(), con.rend(), f);
	std::cout << std::endl;
}
 


/////////////////////////////////////////////////////////////////////////////////////////////////// 未分类自定义类型：
	
// 基于std::chrono的自定义计时器
using namespace std::chrono;
class tiktok
{ 
// 禁用constructor和destructor
private:
	tiktok() = default;
	tiktok(const tiktok&) {}
	~tiktok() = default;

// 成员数据：
public:
	time_point<steady_clock> startTik;
	time_point<steady_clock> endTik;
	unsigned recordCount;
	std::vector<time_point<steady_clock>> records;
	double lastDur = 0;									// 上一次结束计时记录下的时间间隔，单位为秒；

// 方法：
public:
	static tiktok& getInstance()
	{
		static tiktok tt_instance;
		return tt_instance;
	}

	// 开始计时
	void start()
	{
		this->startTik = steady_clock::now();
		this->recordCount = 0;
		this->records.clear();
	}


	// 结束计时
	void end()
	{
		this->endTik = steady_clock::now();
		microseconds duration = duration_cast<microseconds>(this->endTik - this->startTik);
		this->lastDur = static_cast<double>(duration.count()) * \
			microseconds::period::num / microseconds::period::den;
	}


	// 结束计时，控制台上打印时间间隔，单位为秒
	void endCout(const char* str)
	{
		end();
		std::cout << str << this->lastDur << " s." << std::endl;
	}


	// 结束计时，时间间隔写入到fileName的文本文件中，单位为秒；
	bool endWrite(const char* fileName, const char* str)
	{
		end();
		std::ofstream file(fileName, std::ios_base::out | std::ios_base::app);
		if (!file)
			return false; 
		file << str << this->lastDur << std::endl;
		file.close();
		return true;
	}


	// 结束计时，返回std::chrono::microseconds类型的时间间隔；
	microseconds endGetCount()
	{
		this->endTik = steady_clock::now();
		microseconds duration = duration_cast<microseconds>(this->endTik - this->startTik);
		this->lastDur = static_cast<double>(duration.count()) * \
			microseconds::period::num / microseconds::period::den;
		return duration;
	}


	// 按下秒表，记录此刻的时刻，保存在this->records向量中；
	void takeArecord()
	{
		this->records.push_back(steady_clock::now());
		recordCount++;
	} 
};
  
  
// class——方向包围盒类OBB（本质上是中点在原点的AABB加上一个仿射变换）；
template <typename _Scalar>
class OBB : public Eigen::AlignedBox<_Scalar, 3>
{
	// 成员数据：注！！！继承自Eigen::AlignedBox中的m_min, m_max是列向量，下面两个是行向量；
public:
	Eigen::Matrix<_Scalar, 1, 3> m_dir;
	Eigen::Matrix<_Scalar, 1, 3> m_center;
	Eigen::Matrix3d m_rotation;

	// constructor:
public:
	inline OBB() :m_dir(0, 0, 1), m_center(Eigen::Matrix<_Scalar, 1, 3>::Zero(3)), m_rotation(Eigen::Matrix3d::Identity()) {}


	inline OBB(const Eigen::AlignedBox<_Scalar, 3>& aabb)
	{
		auto minp = aabb.m_min;
		auto maxp = aabb.m_max;
		this->center = (minp + maxp) / 2.0;
		this->m_dir = Eigen::Matrix<_Scalar, 1, 3>(0, 0, 1);
		this->m_rotation = Eigen::Matrix3d::Identity();

		this->m_min = minp - this->center.transpose();
		this->m_max = maxp - this->center.transpose();
	}


	inline OBB(const Eigen::AlignedBox<_Scalar, 3>& aabb, const Eigen::Matrix<_Scalar, 1, 3>& dir, \
		const Eigen::Matrix<_Scalar, 1, 3>& center) : m_dir(dir), m_center(center)
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


	inline OBB(const Eigen::AlignedBox<_Scalar, 3>& aabb, const Eigen::Matrix<_Scalar, 1, 3>& center, \
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


template <typename _Scalar>
void genOBBmesh(const OBB<_Scalar>& obb, Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	Eigen::MatrixXi& tris)
{
	Eigen::Matrix<_Scalar, 3, 1> minp = obb.min();
	Eigen::Matrix<_Scalar, 3, 1> maxp = obb.max();
	Eigen::Matrix<_Scalar, 3, 1> sizeVec = maxp - minp;

	genCubeMesh(vers, tris);
	vers.col(0) *= sizeVec(0);
	vers.col(1) *= sizeVec(1);
	vers.col(2) *= sizeVec(2);

	vers = ((obb.m_rotation * vers.transpose()).eval()).transpose().eval();
	vers.rowwise() += obb.m_center;
}



/////////////////////////////////////////////////////////////////////////////////////////////////// 测试函数：
namespace TEST_MYEIGEN
{
	void test0();
	void test1();
	void test11();
	void test111();
	void test1111();
	void test2();
	void test3();
	void test4();
	void test5();
	void test6();
	void test7();
	void test8();
	void test9();
	void test99();
	void test10();

	void test12();
	void test13();
	void test14();
	void test15();
	void test16();
}


namespace TEST_MYEIGEN_IO 
{
	void test0();
	void test1(); 
	void test2();
	void test3();
	void test4();
	void test5();
	void test6();
	void test7();
	void test8();
	void test9(); 
}


namespace TEST_MYEIGEN_BASIC_MATH
{
	void test0();
	void test1();
	void test2();
	void test3();
	void test4();
	void test5();
	void test6();
	void test7();
	void test8();
	void test9();
}

 
// 测试myEigenPMP中的接口
namespace TEST_MYEIGEN_PMP 
{
	void test0();
	void test1();
	void test2();
	void test3();
	void test4();
	void test5();
	void test6();
	void test7();
	void test8();
	void test9();
	void test10();
}


// 数字图像处理
namespace TEST_DIP 
{
	void test0();

}


#include "myEigen.tpp"

 
