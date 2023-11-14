#include "Eigen/Dense"
#include "triMesh.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////// 表象之间的转换：

template <typename DerivedV, typename TV>
void vers2mat(Eigen::PlainObjectBase<DerivedV>& versMat, \
	const std::vector<TRIANGLE_MESH::triplet<TV>>& vers)
{
	using ScalarV = typename DerivedV::Scalar;
	const int versCount = static_cast<int>(vers.size());
	versMat.resize(0, 0);
	versMat.resize(versCount, 3);
	for (int i = 0; i < versCount; ++i)
	{
		versMat(i, 0) = static_cast<ScalarV>(vers[i].x);
		versMat(i, 1) = static_cast<ScalarV>(vers[i].y);
		versMat(i, 2) = static_cast<ScalarV>(vers[i].z);
	}
}


template <typename TV, typename DerivedV>
void mat2vers(std::vector<TRIANGLE_MESH::triplet<TV>>& vers, \
		const Eigen::MatrixBase<DerivedV>& versMat)
{
	const int versCount = static_cast<int>(versMat.rows());
	vers.clear();
	vers.resize(versCount);
	for (int i = 0; i < versCount; ++i)
	{
		vers[i].x = static_cast<TV>(versMat(i, 0));
		vers[i].y = static_cast<TV>(versMat(i, 1));
		vers[i].z = static_cast<TV>(versMat(i, 2));
	}
}


template <typename DerivedV, typename TV, typename TI>
void triMesh2mat(Eigen::PlainObjectBase<DerivedV>& versMat,\
	Eigen::MatrixXi& trisMat, const TRIANGLE_MESH::triMesh<TV, TI>& mesh)
{
	using ScalarV = typename DerivedV::Scalar;
	const int versCount = static_cast<int>(mesh.vertices.size());
	const int trisCount = static_cast<int>(mesh.triangles.size());
	versMat.resize(0, 0);
	trisMat.resize(0, 0);
	versMat.resize(versCount, 3);
	for (int i = 0; i < versCount; ++i)
	{
		versMat(i, 0) = static_cast<ScalarV>(mesh.vertices[i].x);
		versMat(i, 1) = static_cast<ScalarV>(mesh.vertices[i].y);
		versMat(i, 2) = static_cast<ScalarV>(mesh.vertices[i].z);
	}
	trisMat.resize(trisCount, 3);
	for (int i = 0; i < trisCount; ++i)
	{
		trisMat(i, 0) = static_cast<int>(mesh.triangles[i].x);
		trisMat(i, 1) = static_cast<int>(mesh.triangles[i].y);
		trisMat(i, 2) = static_cast<int>(mesh.triangles[i].z);
	}
}


template<typename TV, typename TI, typename DerivedV, typename DerivedT>
void mat2triMesh(TRIANGLE_MESH::triMesh<TV, TI>& mesh, \
		const Eigen::MatrixBase<DerivedV>& versMat, \
		const Eigen::MatrixBase<DerivedT>& trisMat)
{  
	const int trisCount = static_cast<int>(trisMat.rows());
	mat2vers(mesh.vertices, versMat); 
	mesh.triangles.clear();
	mesh.triangles.resize(trisCount);
	for (int i = 0; i < trisCount; ++i)
	{
		mesh.triangles[i].x = static_cast<TI>(trisMat(i, 0));
		mesh.triangles[i].y = static_cast<TI>(trisMat(i, 1));
		mesh.triangles[i].z = static_cast<TI>(trisMat(i, 2));
	}


}