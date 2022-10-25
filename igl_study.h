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
#include <deque>
#include <forward_list>
#include <algorithm>
#include <functional>

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
#include <igl/circulation.h>
#include <igl/remove_unreferenced.h>
#include <igl/connect_boundary_to_infinity.h>
#include <igl/is_edge_manifold.h>

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

#include <igl/decimate.h>
#include <igl/collapse_edge.h>
#include <igl/collapse_small_triangles.h>
#include <igl/qslim.h>
#include <igl/per_vertex_point_to_plane_quadrics.h>
#include <igl/quadric_binary_plus_operator.h>

#include <igl/winding_number.h>

#include <igl/copyleft/cgal/mesh_boolean.h>

#include <igl/topological_hole_fill.h>

#include "tutorial_shared_path.h"


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
	void test2();
	void test3();
	void test33();
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