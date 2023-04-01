#pragma once

#include "myEigen.h"
#include "myIgl.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <queue>
#include <deque>
#include <forward_list>
#include <algorithm>
#include <functional>
#include <mutex>

#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/readSTL.h>
#include <igl/writeSTL.h>
#include <igl/readMESH.h>
#include <igl/writeMESH.h>

#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/barycenter.h>
#include <igl/edge_flaps.h>
#include <igl/unique_edge_map.h>
#include <igl/doublearea.h>
#include <igl/slice.h>
#include <igl/slice_mask.h>
#include <igl/circulation.h>
#include <igl/remove_unreferenced.h>
#include <igl/connect_boundary_to_infinity.h>
#include <igl/is_edge_manifold.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/remove_unreferenced.h>
#include <igl/connected_components.h>

#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/readDMAT.h>
#include <igl/repdiag.h>
#include <igl/opengl/glfw/Viewer.h>

#include <igl/adjacency_list.h>
#include <igl/adjacency_matrix.h>
#include <igl/dfs.h>
#include <igl/bfs.h>
#include <igl/dijkstra.h>

#include <igl/AABB.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/voxel_grid.h>
#include <igl/signed_distance.h>
#include <igl/marching_cubes.h>
#include <igl/march_cube.h>
#include "MC_tables.h"

#include <igl/decimate.h>
#include <igl/collapse_edge.h>
#include <igl/collapse_small_triangles.h>
#include <igl/qslim.h>
#include <igl/per_vertex_point_to_plane_quadrics.h>
#include <igl/quadric_binary_plus_operator.h>
#include <igl/decimate_trivial_callbacks.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/qslim_optimal_collapse_edge_callbacks.h>
#include <igl/hausdorff.h>

#include <igl/winding_number.h>

#include <igl/copyleft/cgal/mesh_boolean.h>

#include <igl/topological_hole_fill.h>

#include "tutorial_shared_path.h"


// 封装的基于libigl的laplace光顺接口：
template <typename T>
bool laplaceFaring(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, const float deltaLB, const unsigned loopCount)
{
	using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix3T = Eigen::Matrix<T, 3, 3>;
	using RowVector3T = Eigen::Matrix<T, 1, 3>;

	// 1. 计算laplace矩阵：
	Eigen::SparseMatrix<T> Lmat;
	//igl::cotmatrix(vers, tris, Lmat);
	cotLaplacian(Lmat, vers, tris);

	// 2. 光顺的循环：
	MatrixXT versCopy = vers;
	for (unsigned i = 0; i < loopCount; ++i)
	{
		// f1. 计算当前质量矩阵：
		Eigen::SparseMatrix<T> mass;
		igl::massmatrix(versCopy, tris, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

		// f2. 解线性方程组 (mass - delta*L) * newVers = mass * newVers
		const auto& S = (mass - deltaLB * Lmat);
		Eigen::SimplicialLLT<Eigen::SparseMatrix<T>> solver(S);
		assert(solver.info() == Eigen::Success);
		versOut = solver.solve(mass * versCopy).eval();
		versCopy = versOut;
	}

	return true;
}



namespace IGL_BASIC
{
	const double pi = 3.14159;
	void test00();
	void test000();
	void test0();
	void test1();
	void test2();
	void test3();
	void test4();
	void test5();
	void test55();
	void test6();
	void test7();
	void test77();
	void test777();
	void test7777();
	void test77777();
	void test777777();
	void test8();
	void test9();
	void test10();
	void test11();
	void test12();
	void test13();
	void test14();
	void test15();
	void test16();

}


namespace IGL_DIF_GEO
{
	const double pi = 3.14159;
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
	void test11();
	void test12();
	void test13();
	void test14();
	void test15();
	void test16();

}


namespace IGL_GRAPH 
{
	const double pi = 3.14159;
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
	void test11();
	void test12();
	void test13();
	void test14();
	void test15();
	void test16();
}


namespace IGL_SPACE_PARTITION
{
	const double pi = 3.14159;
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
	void test11();
	void test12();
	void test13();
	void test14();
	void test15();
	void test16();
}


// IGL实现的基础三角网格处理算法；
namespace IGL_BASIC_PMP
{
	const double pi = 3.14159;
	void test0();
	void test1();
	void test11();
	void test2();
	void test3();
	void test33();
	void test4();
	void test44();
	void test5();
	void test6();
	void test7();
	void test77();
	void test8();
	void test9();
	void test10();
	void test11();
	void test12();
	void test13();
	void test14();
	void test15();
	void test16();
}