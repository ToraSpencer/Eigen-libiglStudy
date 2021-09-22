#include "dense_mat.h"
#include "sparse_mat.h"
#include "scientific_calc.h"
#include "igl_study.h"

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
	使用的第三方库:
		eigen
		libigl				x64
		glad
		glfw
*/


// 项目中几个预定义宏：
/*
	IGL_STATIC_LIBRARY
	
	NOMINMAX
			解决windows.h和<algorithm>中同时定义std::max()，std::min()的冲突。		
	
	TUTORIAL_SHARED_PATH="G:/gitRepositories/ligIgl/libigl-main/cmake/../external/../tutorial/data"
	
	CMAKE_INTDIR="Release"
*/


 

int main()
{
	DENSEMAT::test1();
 
	//SPARSEMAT::test1();

	//IGLSTUDY::test1();

	//SCIENTIFICCALC::test1();

}
