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


/////////////////////////////////////////////////////////////////////////////////////////////////// �Լ�д�ĸ�����̬�⣺
#include "myEigenIO/myEigenIO.h"
#pragma comment(lib,"myEigenIO.lib")	

#include "myEigenBasicMath/myEigenBasicMath.h"
#pragma comment(lib, "myEigenBasicMath.lib")

#include "myEigenPMP/myEigenPMP.h"
#pragma comment(lib, "myEigenPMP.lib")

#include "myEigenModeling/myEigenModeling.h"
#pragma comment(lib, "myEigenModeling.lib")

 
/////////////////////////////////////////////////////////////////////////////////////////////////// debug tools

// ���뺯���ӻ���ָ�����stl����
template<typename T, typename F>
void traverseSTL(T& con, F f)
{
	std::for_each(con.begin(), con.end(), f);
	std::cout << std::endl;
}


// �������
template<typename T, typename F>
void revTraverseSTL(T& con, F f)
{
	std::for_each(con.rbegin(), con.rend(), f);
	std::cout << std::endl;
}
 


/////////////////////////////////////////////////////////////////////////////////////////////////// δ�����Զ������ͣ�
	
// ����std::chrono���Զ����ʱ��
using namespace std::chrono;
class tiktok
{ 
// ����constructor��destructor
private:
	tiktok() = default;
	tiktok(const tiktok&) {}
	~tiktok() = default;

// ��Ա���ݣ�
public:
	time_point<steady_clock> startTik;
	time_point<steady_clock> endTik;
	unsigned recordCount;
	std::vector<time_point<steady_clock>> records;
	double lastDur = 0;									// ��һ�ν�����ʱ��¼�µ�ʱ��������λΪ�룻

// ������
public:
	static tiktok& getInstance()
	{
		static tiktok tt_instance;
		return tt_instance;
	}

	// ��ʼ��ʱ
	void start()
	{
		this->startTik = steady_clock::now();
		this->recordCount = 0;
		this->records.clear();
	}


	// ������ʱ
	void end()
	{
		this->endTik = steady_clock::now();
		microseconds duration = duration_cast<microseconds>(this->endTik - this->startTik);
		this->lastDur = static_cast<double>(duration.count()) * \
			microseconds::period::num / microseconds::period::den;
	}


	// ������ʱ������̨�ϴ�ӡʱ��������λΪ��
	void endCout(const char* str)
	{
		end();
		std::cout << str << this->lastDur << " s." << std::endl;
	}


	// ������ʱ��ʱ����д�뵽fileName���ı��ļ��У���λΪ�룻
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


	// ������ʱ������std::chrono::microseconds���͵�ʱ������
	microseconds endGetCount()
	{
		this->endTik = steady_clock::now();
		microseconds duration = duration_cast<microseconds>(this->endTik - this->startTik);
		this->lastDur = static_cast<double>(duration.count()) * \
			microseconds::period::num / microseconds::period::den;
		return duration;
	}


	// ���������¼�˿̵�ʱ�̣�������this->records�����У�
	void takeArecord()
	{
		this->records.push_back(steady_clock::now());
		recordCount++;
	} 
};
  
  
// class���������Χ����OBB�����������е���ԭ���AABB����һ������任����
template <typename _Scalar>
class OBB : public Eigen::AlignedBox<_Scalar, 3>
{
	// ��Ա���ݣ�ע�������̳���Eigen::AlignedBox�е�m_min, m_max����������������������������
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
		// ע�������Χ�еķ���ȡdir��-dirû���������Կ���dir��z������нǲ�����pi/2;
		const double epsilon = 1e-6;
		assert((dir.norm() - 1) < epsilon && "input dir is invalid.");

		double dotValue = Eigen::Matrix<_Scalar, 1, 3>(0, 0, 1).dot(dir);
		if (dotValue < 0)
			this->m_dir = -dir;
		Eigen::Matrix<_Scalar, 1, 3> crossValue = Eigen::Matrix<_Scalar, 1, 3>(0, 0, 1).cross(this->m_dir);		// ��ת�᣻

		auto minp = aabb.min();
		auto maxp = aabb.max();
		Eigen::Matrix<_Scalar, 1, 3> centerOri = (minp + maxp) / 2.0;
		this->m_min = minp - centerOri.transpose();
		this->m_max = maxp - centerOri.transpose();

		// ��ָ�������z��ƽ�У�
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

		// ���������任��AABB�Ŀռ䣺
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



/////////////////////////////////////////////////////////////////////////////////////////////////// ���Ժ�����
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

 
// ����myEigenPMP�еĽӿ�
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


// ����ͼ����
namespace TEST_DIP 
{
	void test0();

}


#include "myEigen.tpp"

 
