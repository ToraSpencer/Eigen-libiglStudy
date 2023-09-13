#pragma once

#include "myEigen.h" 

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

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
#include <igl/readDMAT.h>

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

// 符号距离场相关
#include <igl/AABB.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/voxel_grid.h>
#include <igl/signed_distance.h>								
#include <igl/marching_cubes.h>
#include <igl/march_cube.h>
#include "MC_tables.h"

// 网格精简相关：
#include <igl/decimate.h>
#include <igl/collapse_edge.h>
#include <igl/qslim.h>
#include <igl/per_vertex_point_to_plane_quadrics.h>
#include <igl/quadric_binary_plus_operator.h>
#include <igl/decimate_trivial_callbacks.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/qslim_optimal_collapse_edge_callbacks.h>
#include <igl/hausdorff.h>

// mesh deformation相关
#include <igl/colon.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/arap.h>                               // ARAP变形

#include <igl/winding_number.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/topological_hole_fill.h>


// 静态库：



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
	void test555();
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
	void test55();
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


// IGL实现的生成数字模型相关算法；
namespace IGL_MODELLING
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


// IGL中网格变形相关；
namespace IGL_DEFORMATION
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