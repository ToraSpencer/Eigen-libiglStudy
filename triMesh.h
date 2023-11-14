#pragma once
#include <vector>

namespace TRIANGLE_MESH
{
	template <typename T>
	struct triplet
	{
		T x;
		T y;
		T z;
	};

	template <typename TV, typename TI>
	struct triMesh
	{
		std::vector<triplet<TV>> vertices;
		std::vector<triplet<TI>> triangles;
	};
}

using verF = TRIANGLE_MESH::triplet<float>;									// 单精度顶点；
using verD = TRIANGLE_MESH::triplet<double>;								// 双精度顶点；
using triMeshF = TRIANGLE_MESH::triMesh<float, int>;					// 单精度顶点网格；
using triMeshD = TRIANGLE_MESH::triMesh<double, int>;				// 双精度顶点网格；

bool readOBJ(std::vector<verF>& vers, const char* fileName);
bool readOBJ(std::vector<verD>& vers, const char* fileName);
bool readOBJ(triMeshF& mesh, const char* fileName);
bool readOBJ(triMeshD& mesh, const char* fileName);
bool writeOBJ(const char* fileName, const std::vector<verF>& vers);
bool writeOBJ(const char* fileName, const std::vector<verD>& vers);
bool writeOBJ(const char* fileName, const triMeshF& mesh);
bool writeOBJ(const char* fileName, const triMeshD& mesh);