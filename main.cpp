#include "dense_mat.h"
#include "sparse_mat.h"
#include "scientific_calc.h"
#include "igl_study.h"

#define DATA_PATH "./data/"

// 当前问题-easy
/*
	
*/


// 当前问题-hard
/*
	1. 给定流形网格上两个点，找出它们之间的最短路径。
	
	2. 使用OPENGL实现上面的问题，可以用鼠标点击确定网格上的任意两个点作为输入。 


*/



// 项目信息
/*
	编译环境： x64 Relase
 
	使用的第三方库:
		eigen
		libigl				 
		glad
		glfw
*/


// 项目中几个预定义宏：
/*
	IGL_STATIC_LIBRARY
	
	NOMINMAX
			解决windows.h和<algorithm>中同时定义std::max()，std::min()的冲突。		
	
	TUTORIAL_SHARED_PATH="G:/gitRepositories/ligIgl/libigl_CGAL_openGL/cmake/../external/../tutorial/data"
	
	CMAKE_INTDIR="Release"
*/

 
int main()
{
	// DENSEMAT::test14();
 
	// SPARSEMAT::test0();
	
	// IGL_BASIC::test6();

	IGL_DIF_GEO::test0();

	// IGL_GRAPH::test2();

	// IGL_SPACE_PARTITION::test0();

	// SCIENTIFICCALC::test6();

	// IGL_BASIC_PMP::test3();
 

	std::cout << "main() finished." << std::endl;
}
