#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>

namespace TRIANGLE_MESH
{
	////////////////////////////////////////////////////////////////////////////////////////////// triplet<>类
	template <typename T>
	struct triplet
	{
		T x;
		T y;
		T z;

		triplet<T>() : x(0), y(0), z(0)
		{}

		triplet<T>(const T x0, const T y0, const T z0) : x(x0), y(y0), z(z0)
		{}

		template <typename Ti>
		triplet<T>(const triplet<Ti>& t): x(static_cast<T>(t.x)), y(static_cast<T>(t.y)), z(static_cast<T>(t.z))
		{}

		T& at(const size_t idx)
		{
			switch (idx)
			{
			case 0:
				return x;
			case 1:
				return y;
			case 2:
				return z;
			default:
				throw(std::out_of_range("exception!!! idx should be 0, 1, or 2;"));
			}
		}

		T& operator[](const size_t idx)
		{
			return at(idx);
		}

		const T& at(const size_t idx) const
		{
			return at(idx);
		}

		const T& operator[](const size_t idx) const
		{
			switch (idx)
			{
			case 0:
				return x;
			case 1:
				return y;
			case 2:
				return z;
			default:
				throw(std::out_of_range("exception!!! idx should be 0, 1, or 2;"));
			}
		}

		const triplet<T>& operator=(const triplet<T>& t)
		{
			this->x = t.x;
			this->y = t.y;
			this->z = t.z;
			return *this;
		}

		template <typename Scalar>
		const triplet<T>& operator+=(const Scalar num)
		{
			this->x = static_cast<T>(x + num);
			this->y = static_cast<T>(y + num);
			this->z = static_cast<T>(z + num);
			return *this;
		}

		template <typename Ti>
		const triplet<T>& operator+=(const triplet<Ti>& t)
		{
			this->x = static_cast<T>(x + t.x);
			this->y = static_cast<T>(y + t.y);
			this->z = static_cast<T>(z + t.z);
			return *this;
		}

		template <typename Scalar>
		const triplet<T>& operator-=(const Scalar num)
		{
			this->x = static_cast<T>(x - num);
			this->y = static_cast<T>(y - num);
			this->z = static_cast<T>(z - num);
			return *this;
		}

		template <typename Ti>
		const triplet<T>& operator-=(const triplet<Ti>& t)
		{
			this->x = static_cast<T>(x - t.x);
			this->y = static_cast<T>(y - t.y);
			this->z = static_cast<T>(z - t.z);
			return *this;
		}

		template <typename Scalar>
		const triplet<T>& operator*=(const Scalar num)
		{
			this->x = static_cast<T>(x * num);
			this->y = static_cast<T>(y * num);
			this->z = static_cast<T>(z * num);
			return *this;
		}

		template <typename Scalar>
		const triplet<T>& operator/=(const Scalar num)
		{
			this->x = static_cast<T>(x / num);
			this->y = static_cast<T>(y / num);
			this->z = static_cast<T>(z / num);
			return *this;
		}

		double length() const 
		{
			return static_cast<double>(std::sqrt(x*x + y*y + z*z));
		}

		template <typename Ti>
		double distance(const triplet<Ti>& t) const 
		{
			triplet<double> arrow{static_cast<double>(x - t.x), \
				static_cast<double>(y - t.y), static_cast<double>(z - t.z), };

			return arrow.length();
		}

		template <typename Ti>
		double dot(const triplet<Ti>& t) const 
		{
			return static_cast<double>(x * t.x + y * t.y + z * t.z);
		}

		template <typename Ti>
		triplet<T> cross(const triplet<Ti>& t) const
		{
			return triplet<T>{static_cast<T>(y* t.z - z * t.y),\
				static_cast<T>(z * t.x - x*t.z), \
				static_cast<T>(x* t.y - y * t.x)}; 
		}

		template <typename Ti>
		triplet<Ti> cast() const 
		{
			return triplet<Ti>{static_cast<Ti>(x), static_cast<Ti>(y), static_cast<Ti>(z)};
		}
	};


	////////////////////////////////////////////////////////////////////////////////////////////// triMesh()三角网格类：
	template <typename TV, typename TI>
	struct triMesh
	{
		std::vector<triplet<TV>> vertices;
		std::vector<triplet<TI>> triangles;
	};


	////////////////////////////////////////////////////////////////////////////////////////////// 重载运算符

	template <typename T>
	std::ostream& operator<<(std::ostream& os, const triplet<T>& t)
	{
		os << "(" << t.x << ", " << t.y << ", " << t.z << ")";
		return os;
	}


	template <typename T1, typename T2>
	triplet<T1> operator+(const triplet<T1>& t1, const triplet<T2>& t2) 
	{
		return triplet<T1>{t1.x + t2.x, t1.y + t2.y, t1.z + t2.z};
	}


	template <typename T1, typename T2>
	triplet<T1> operator-(const triplet<T1>& t1, const triplet<T2>& t2)
	{
		return triplet<T1>{t1.x - t2.x, t1.y - t2.y, t1.z - t2.z};
	}

	// 负号；
	template <typename T>
	triplet<T> operator-(const triplet<T>& t)				
	{
		return triplet<T>{-t.x, -t.y, -t.z};
	}

	template <typename T, typename Scalar>
	triplet<T> operator*(const triplet<T>& t, const Scalar num)
	{
		return triplet<T>{static_cast<T>(num* t.x), static_cast<T>(num* t.y), static_cast<T>(num* t.z)};
	}


	template <typename Scalar, typename T>
	triplet<T> operator*(const Scalar num, const triplet<T>& t)
	{
		return t * num;
	} 

	template <typename T, typename Scalar>
	triplet<T> operator/(const triplet<T>& t, const Scalar num)
	{
		return triplet<T>{static_cast<T>(t.x / num), static_cast<T>(t.y / num), static_cast<T>(t.z / num)};
	}
}

using verF = TRIANGLE_MESH::triplet<float>;									// 单精度顶点；
using verD = TRIANGLE_MESH::triplet<double>;								// 双精度顶点；
using triMeshF = TRIANGLE_MESH::triMesh<float, int>;					// 单精度顶点网格；
using triMeshD = TRIANGLE_MESH::triMesh<double, int>;				// 双精度顶点网格；

bool readOBJ(std::vector<verF>& vers, const char* fileName);
bool readOBJ(std::vector<verD>& vers, const char* fileName);
bool readOBJ(triMeshF& mesh, const char* fileName);
bool readOBJ(triMeshD& mesh, const char* fileName);
bool readSTL(triMeshF& mesh, const char* fileName, const bool blIsAscii = false);
bool readSTL(triMeshD& mesh, const char* fileName, const bool blIsAscii = false);
bool writeOBJ(const char* fileName, const std::vector<verF>& vers);
bool writeOBJ(const char* fileName, const std::vector<verD>& vers);
bool writeOBJ(const char* fileName, const triMeshF& mesh);
bool writeOBJ(const char* fileName, const triMeshD& mesh);