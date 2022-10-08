#include "igl_study.h"

#define DATA_PATH "./data/"


// libigl基本功能
namespace IGL_BASIC
{
	Eigen::MatrixXd vers, newVers, normals;
	Eigen::MatrixXi tris;
	Eigen::SparseMatrix<double> L;
	igl::opengl::glfw::Viewer viewer;		// libigl中的基于glfw的显示窗口；

	// 文件IO
	void test0()
	{
		//// readOBJ(), writeObj()——OBJ文件的IO，有多个重载，至少有路径、点云、三角片这三个数据。
		igl::readOBJ("./data/bunny.obj", vers, tris);
		igl::writeOBJ("./data/bunny_export.obj", vers, tris);
		igl::writeOBJ("./data/bunnyVers.obj", vers, Eigen::MatrixXi{});			// 只要点云不需要三角片的话，传入空矩阵；
	
		//vers.resize(0, 0);
		//tris.resize(0, 0);
		//igl::readOBJ("E:/fatTeeth1_预处理后.obj", vers, tris);
		//igl::writeOFF("E:/fatTeeth1_预处理后.off", vers, tris);

		// 读取stl文件；
		std::string fileName{"E:/材料/jawMeshSimplified"};
		vers.resize(0, 0);
		tris.resize(0, 0);
		std::ifstream fileIn((fileName + std::string{ ".stl" }).c_str(), std::ios::binary);			// stl文件是二进制文件；
		igl::readSTL(fileIn, vers, tris, normals);
		fileIn.close();
		igl::writeOBJ((fileName + std::string{".obj"}).c_str(), vers, tris);
 

		std::cout << "finished." << std::endl;
	}


	// libigl中的显示窗口类Viewer
	void test1() 
	{
		igl::readOBJ("bunny.obj", vers, tris);

		// Viewer::data()——返回viewer中数据对象的引用；

		// ViewerData::set_mesh()——输入顶点和三角片矩阵生成网格数据，写入到ViewerData的成员变量中；
		viewer.data().set_mesh(vers, tris);		// 窗口中装载数据；


		 // 默认是两轴旋转，set_rotation_type()方法可以指定其他旋转风格
		viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);				// 三轴旋转

		viewer.data().show_lines = 0;			// show_lines指定是否画出网格线；

		viewer.launch();
	}

 
	// 计算梯度
	void test2()
	{
		using namespace Eigen;
		using namespace std;

		igl::readOBJ("./data/rootTooth1.obj", vers, tris);

		VectorXd U;
		igl::readDMAT("./data/cheburashka-scalar.dmat", U);		// 一系列的函数值

		SparseMatrix<double> G;			// 梯度算子
		igl::grad(vers, tris, G);

		// 计算函数值U的梯度
		MatrixXd GU = Map<const MatrixXd>((G * U).eval().data(), tris.rows(), 3);

		// 函数值U的梯度GU的模长
		const VectorXd GU_mag = GU.rowwise().norm();

		viewer.data().set_mesh(vers, tris);
		viewer.data().set_data(U);

		// Average edge length divided by average gradient (for scaling)
		const double max_size = igl::avg_edge_length(vers, tris) / GU_mag.mean();

		// 每个三角片重心上画一根指示线，方向为梯度方向。 
		MatrixXd BC;
		igl::barycenter(vers, tris, BC);
		const RowVector3d black(0, 0, 0);
		viewer.data().add_edges(BC, BC + max_size * GU, black);
		viewer.data().show_lines = false;	  // 隐藏网格线

		viewer.launch();
	}
 

	// lambda——键盘事件：使用laplacian光顺网格
	bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)
	{
		switch (key)
		{
		case 'r':

		case 'R':			// 复位程序
			newVers = vers;
			break;

		case ' ':				// 空格键，执行一次laplace光顺
		{
			// 重新计算质量矩阵
			Eigen::SparseMatrix<double> mass;
			igl::massmatrix(newVers, tris, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

			// 解线性方程组 (mass - delta*L) * newVers = mass * newVers
			float delta = 0.001;
			const auto& S = (mass - delta * L);
			Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
			assert(solver.info() == Eigen::Success);
			newVers = solver.solve(mass * newVers).eval();
 
			break;
		}

		default:
			return false;
		}

		viewer.data().set_vertices(newVers);
		viewer.data().compute_normals();
		viewer.core().align_camera_center(newVers, tris);
		return true;
	};


	// 使用Laplacian光顺网格
	void test3() 
	{
		 igl::readOBJ( "./data/bunny.obj", vers, tris);
		newVers = vers;

		// 1.a 直接构造laplacian——Compute Laplace-Beltrami operator: 
		igl::cotmatrix(vers, tris, L);

		// 1.b 分步构造laplacian
		{
			SparseMatrix<double> Gradient, L2;

			igl::grad(vers, tris, Gradient);      // 离散梯度

			// Diagonal per-triangle "mass matrix"
			VectorXd dblA;
			igl::doublearea(vers, tris, dblA);             // 每个三角片面积的两倍

			// Place areas along diagonal  #dim times
			const auto& T = 1. * (dblA.replicate(3, 1) * 0.5).asDiagonal();

			L2 = -Gradient.transpose() * T * Gradient;         // discrete Dirichelet energy Hessian 离散狄利克雷能量海塞矩阵
			std::cout << "两种方法得到的laplacian的差的范数：" << std::endl;
			cout << "(L2 - L).norm() == " << (L2 - L).norm() << endl;
		}

		// 2. 根据原始的法向量，使用伪色
		MatrixXd norms;
		igl::per_vertex_normals(vers, tris, norms);
		MatrixXd colors = norms.rowwise().normalized().array() * 0.5 + 0.5;

		// 3. viewr填充初始数据
		newVers = vers;
		viewer.data().set_mesh(newVers, tris);
		viewer.data().set_colors(colors);
		viewer.callback_key_down = key_down;

		// 4. 运行
		cout << "Press [space] to smooth." << endl;;
		cout << "Press [r] to reset." << endl;;
		viewer.launch();
	}


	// 生成三维空间中的栅格：
	void test4() 
	{
		Eigen::MatrixXd vers, gridCenters;
		Eigen::MatrixXi tris;	
		Eigen::RowVector3i gridCounts;			// 三个维度上栅格的数目；
		unsigned num = 10;				// 跨度最大的那个维度(xyz中的一个)的栅格数；

		igl::readOBJ("E:/材料/tooth.obj", vers, tris);
		igl::voxel_grid(vers, 0, num, 1, gridCenters, gridCounts);

		igl::writeOBJ("E:/gridCenters.obj", gridCenters, Eigen::MatrixXi{});

		std::cout << "finished." << std::endl;
	}


	// 计算符号距离场、marching cubes提取等值面网格；
	void test5()
	{
		MatrixXi tris;
		MatrixXd vers;
		igl::readOBJ("E:/材料/jawMesh.obj", vers, tris);

		tiktok& tt = tiktok::getInstance();
		double gridStep = 0.5;			
		double range[3];
		for (unsigned i = 0; i < 3; ++i)
		{
			Eigen::VectorXd coors = vers.col(i);
			range[i] = coors.maxCoeff() - coors.minCoeff();
		}

		// 1. 生成栅格：
		double largestRange = std::max({ range[0], range[1], range[2] });
		int num = std::ceil((largestRange/gridStep));              // 跨度最大的那个维度(xyz中的一个)的栅格数；
		double gridsOffset = (gridStep * num - largestRange) / 2.;
		
		Eigen::MatrixXd gridCenters;
		Eigen::RowVector3i gridCounts;
		tt.start();
		igl::voxel_grid(vers, gridsOffset, num, 1, gridCenters, gridCounts);
		tt.endCout("Elapsed time of igl::voxel_grid() is ");

		VectorXd SDF, signValues;

		{
			VectorXi I;                     // useless
			MatrixXd C, N;              // useless

			// 2. 计算符号距离场
			tt.start();
			igl::signed_distance(gridCenters, vers, tris, igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, SDF, I, C, N);
			tt.endCout("Elapsed time of igl::signed_distance() is ");

			// 3. 符号距离场改写为符号场——网格内为-1，网格面上为0，外面为1：
			signValues = SDF;
			for_each(signValues.data(), signValues.data() + signValues.size(), [](double& b)\
			{
				b = (b > 0 ? 1 : (b < 0 ? -1 : 0));
			});
		}

		// 4. marching cubes算法生成最终曲面：
		MatrixXd versResult_SDF, versResults_signs;
		MatrixXi trisResult_SDF, trisResults_signs;
		double selectedSDF = -1.;
		tt.start();
		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), selectedSDF, versResult_SDF, trisResult_SDF);
		tt.endCout("Elapsed time of igl::marching_cubes() is ");
		
		// igl::marching_cubes(signValues, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), 0, versResults_signs, trisResults_signs);
		igl::writeOBJ("E:/shrinkedMesh.obj", versResult_SDF, trisResult_SDF);


		std::cout << "finished." << std::endl;
	}


	// 读取SDFGen.exe生成的.sdf距离场数据，使用igl::marching_cubes()提取等值面网格：

	//		解析.sdf文本文件：
	double parseSDF(std::vector<int>& stepCounts, Eigen::RowVector3d& gridsOri, Eigen::VectorXd& SDF, const char* filePath)
	{
		double SDFstep = -1;
		std::string readStr(1024, '\0');
		std::ifstream sdfFile(filePath);
		if (!sdfFile)
			return SDFstep;

		// 第一行：步数
		{
			std::string tmpStr;
			sdfFile.getline(&readStr[0], 1024);

			unsigned index = 0;
			for (const auto& ch : readStr)
			{
				if (ch >= '0' && ch <= '9' || ch == '.')
					tmpStr.push_back(ch);
				else
				{
					if (tmpStr.size() > 0)
					{
						stepCounts[index] = std::stoi(tmpStr);
						index++;
						tmpStr.clear();
					}
				}
			}
		}

		// 第二行：栅格原点
		{
			std::string tmpStr;
			sdfFile.getline(&readStr[0], 1024);

			unsigned index = 0;
			for (const auto& ch : readStr)
			{
				if (ch >= '0' && ch <= '9' || ch == '.')
					tmpStr.push_back(ch);
				else
				{
					if (tmpStr.size() > 0)
					{
						gridsOri(index) = std::stod(tmpStr);
						index++;
						tmpStr.clear();
					}
				}
			}
		}

		// 第三行：距离场空间步长：
		sdfFile.getline(&readStr[0], 1024);
		SDFstep = std::stod(readStr);

		// 第三行之后：距离场数据：
		unsigned dataCount = stepCounts[0] * stepCounts[1] * stepCounts[2];
		SDF.resize(dataCount);
		for (unsigned i = 0; i < dataCount; ++i)
		{
			sdfFile.getline(&readStr[0], 1024);
			SDF(i) = std::stod(readStr);
		}
		sdfFile.close();

		return SDFstep;
	}


	void test55() 
	{
		// 0. 解析SDFGen.exe生成的.sdf距离场数据文件：
		std::vector<int> stepCounts(3);				// xyz三个维度上栅格数
		Eigen::RowVector3d gridsOri;					// 栅格原点：
		Eigen::VectorXd SDF;
		const char* sdfFilePath = "E:/inputMesh.sdf";
		double SDFstep = parseSDF(stepCounts, gridsOri, SDF, sdfFilePath);

		// 1. 生成栅格：
		Eigen::RowVector3i gridCounts;
		Eigen::MatrixXd gridCenters;
		Eigen::RowVector3d minp = gridsOri - SDFstep * Eigen::RowVector3d(stepCounts[0] / 2.0, stepCounts[1] / 2.0, stepCounts[2] / 2.0);
		Eigen::RowVector3d maxp = gridsOri + SDFstep * Eigen::RowVector3d(stepCounts[0] / 2.0, stepCounts[1] / 2.0, stepCounts[2] / 2.0);
		Eigen::AlignedBox<double, 3> box(minp, maxp);
		igl::voxel_grid(box, std::max({stepCounts[0], stepCounts[1], stepCounts[2]}), 0, gridCenters, gridCounts);
		
		// 2. marching cubes算法生成最终曲面：
		tiktok& tt = tiktok::getInstance();
		MatrixXd versResult_SDF, versResults_signs;
		MatrixXi trisResult_SDF, trisResults_signs;
		double selectedSDF = -1.;
		tt.start();
		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), selectedSDF, versResult_SDF, trisResult_SDF);
		tt.endCout("Elapsed time of igl::marching_cubes() is ");

		igl::writeOBJ("E:/shrinkedMesh.obj", versResult_SDF, trisResult_SDF);
 
		std::cout << "finished." << std::endl;
	}


	// 轴向包围盒类，方向包围盒类
	void test6() 
	{
		// 生成包围盒网格：
		Eigen::MatrixXd aabbVers, obbVers;
		Eigen::MatrixXi aabbTris, obbTris;
		Eigen::RowVector3d minp = Eigen::RowVector3d(-1, 2, -3);
		Eigen::RowVector3d maxp = Eigen::RowVector3d(4, 5, 6);
		Eigen::AlignedBox<double, 3> aabb(minp, maxp);
		genAABBmesh(aabb, aabbVers, aabbTris);
		igl::writeOBJ("E:/aabbMesh.obj", aabbVers, aabbTris);
 
		OBB<double> obb(aabb, Eigen::RowVector3d(3,4,5).normalized(), Eigen::RowVector3d(1, 2, 3) );
		genOBBmesh(obb, obbVers, obbTris);
		igl::writeOBJ("E:/obbMesh.obj", obbVers, obbTris);
 
		// contains()——判断顶点与包围盒的关系，顶点需要用列向量表示；在包围盒表面上也视为在内部；
		std::cout << "aabb.contains (3,4,5) ? " << aabb.contains(Eigen::RowVector3d(3, 4, 5).transpose()) << std::endl;
		std::cout << "aabb.contains (4,5,5) ? " << aabb.contains(Eigen::RowVector3d(4, 5, 5).transpose()) << std::endl;
		std::cout << "aabb.contains (9,9,9) ? " << aabb.contains(Eigen::RowVector3d(9, 9, 9).transpose()) << std::endl;

		// OBB类重写的contains()方法：
		std::cout << "obb.contains(0, 0, 0) ? " << obb.contains(Eigen::RowVector3d(0, 0, 0)) << std::endl;
		std::cout << "obb.contains(4, 2, 6) ? " << obb.contains(Eigen::RowVector3d(4, 2, 6)) << std::endl;
		std::cout << "obb.contains(-5, -5, -5) ? " << obb.contains(Eigen::RowVector3d(-5, -5, -5)) << std::endl;
 
		igl::writeOBJ("E:/426.obj", MatrixXd{ RowVector3d(4, 2, 6) }, MatrixXi{});

		std::cout << "finished." << std::endl;
	}


	// 网格精简：
	void test7() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
		Eigen::VectorXi newOldVersInfo;						 
		tiktok& tt = tiktok::getInstance();

		// igl::readOBJ("E:/材料/jawMeshDense.obj", vers, tris);
		igl::readOBJ("E:/材料/tooth.obj", vers, tris);
		unsigned trisCount = tris.rows();
		unsigned tarTrisCount = std::round(trisCount * 0.5);


		// 当前使用igl::decimate()简化简单的网格可以成功，太复杂的网格会失败；
		tt.start();
		std::cout << "succeeded? " << igl::decimate(vers, tris, tarTrisCount, versOut, trisOut, newOldTrisInfo, newOldVersInfo) << std::endl;
		tt.endCout("Elapsed time of mesh simplification is ");
		std::vector<int> newOldTrisInfoVec = vec2Vec(newOldTrisInfo);
		Eigen::MatrixXd vers1;
		subFromIdxVec(vers1, vers, newOldVersInfo);

		igl::writeOBJ("E:/meshIn.obj", vers, tris);
		igl::writeOBJ("E:/meshSimplified.obj", versOut, trisOut);
		igl::writeOBJ("E:/meshSimplifiedVers.obj", vers1, Eigen::MatrixXi{});

		std::cout << "finished." << std::endl;
	}
}


// libigl中的微分几何相关
namespace IGL_DIF_GEO 
{

	// 质量矩阵和LB算子
	void test0() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		Eigen::SparseMatrix<double> L, M;

		igl::readOBJ("E:/材料/tooth.obj", vers, tris);
		igl::cotmatrix(vers, tris, L);
		igl::massmatrix(vers, tris, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);

		dispSpMat(M, 0, M.rows() - 1, 10);

		std::cout << "finished." << std::endl;
	}


	// libigl中的网格布尔操作；！！！需要配置CGAL库，当前未完成；
#if 0
	// https://stackoverflow.com/questions/37898084/how-to-perform-boolean-operation-on-thin-surface-using-libigl
	Eigen::MatrixXd VA, VB, VC;
	Eigen::VectorXi J, I;
	Eigen::MatrixXi FA, FB, FC;
	igl::MeshBooleanType boolean_type(
		igl::MESH_BOOLEAN_TYPE_UNION);

	const char* MESH_BOOLEAN_TYPE_NAMES[] =
	{
	  "Union",
	  "Intersect",
	  "Minus",
	  "XOR",
	  "Resolve",
	};

	const auto& key_down = [](igl::opengl::glfw::Viewer& viewer, unsigned char key, int mods) -> bool
	{
		switch (key)
		{
		default:
			return false;
		case 'A':
			viewer.data().clear();
			std::cout << "Loading A" << std::endl;
			viewer.data().set_mesh(VA, FA);
			break;
		case 'B':
			viewer.data().clear();
			std::cout << "Loading B" << std::endl;
			viewer.data().set_mesh(VB, FB);
			break;
		case 'C':
			viewer.data().clear();
			std::cout << "Loading C" << std::endl;
			viewer.data().set_mesh(VC, FC);
			return true;
		}

		return true;
	};

	void test1()
	{
		using namespace Eigen;
		using namespace std;

		double prismSize = 150;
		double Heigh = 300;
		VA.resize(6, 3);
		VA << -prismSize, prismSize, 0,
			prismSize, prismSize, 0,
			0, 2 * prismSize, 0,
			-prismSize, prismSize, Heigh,
			prismSize, prismSize, Heigh,
			0, 2 * prismSize, Heigh;
		FA.resize(8, 3);
		FA << 1, 0, 2,
			5, 3, 4,
			4, 1, 2,
			2, 5, 4,
			3, 5, 2,
			2, 0, 3,
			0, 1, 4,
			4, 3, 0;

		double tetsize = 300;
		VB.resize(4, 3);
		VB << 0, 0, tetsize,
			-tetsize, 0, 0,
			tetsize, 0, 0,
			0, tetsize * 2, 0;
		FB.resize(4, 3);
		FB << 2, 1, 3,
			2, 0, 1,
			3, 0, 2,
			1, 0, 3;


		igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, igl::MESH_BOOLEAN_TYPE_INTERSECT, VC, FC);

		std::cout
			<< "VA:" << std::endl << VA << std::endl << "==============" << std::endl
			<< "FA:" << std::endl << FA << std::endl << "==============" << std::endl
			<< "VB:" << std::endl << VB << std::endl << "==============" << std::endl
			<< "FB:" << std::endl << FB << std::endl << "==============" << std::endl
			<< "VC:" << std::endl << VC << std::endl << "==============" << std::endl
			<< "FC:" << std::endl << FC << std::endl << "==============" << std::endl;

		// Plot the mesh with pseudocolors
		igl::opengl::glfw::Viewer viewer;

		viewer.data().set_mesh(VA, FA);
 

		viewer.data().show_lines = 1;
		viewer.callback_key_down = key_down;
		viewer.core().camera_dnear = 3.9;
		cout <<
			"Press '.' to switch to next boolean operation type." << endl <<
			"Press ',' to switch to previous boolean operation type." << endl <<
			"Press ']' to push near cutting plane away from camera." << endl <<
			"Press '[' to pull near cutting plane closer to camera." << endl <<
			"Hint: investigate _inside_ the model to see orientation changes." << endl;
		viewer.launch();
	}
#endif

}


// 图算法
namespace IGL_GRAPH 
{
	// 图数据结构的转换：
	void test0() 
	{
		MatrixXd vers;
		MatrixXi tris;
		igl::readOBJ("E:/材料/cube重复三角片.obj", vers, tris);

		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		
		// 由面片信息得到图的邻接表、邻接矩阵
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);		// 有向边存在则元素为1，否则为0；

		// 手动计算邻接矩阵：
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		unsigned edgesCount = 3 * trisCount;
		Eigen::SparseMatrix<int> adjSM_eCount;				// 有向边邻接矩阵，权重为该有向边的数量；

		MatrixXi edges = MatrixXi::Zero(edgesCount, 2);
		edges.block(0, 0, trisCount, 1) = tris.col(0);
		edges.block(0, 1, trisCount, 1) = tris.col(1);
		edges.block(trisCount, 0, trisCount, 1) = tris.col(1);
		edges.block(trisCount, 1, trisCount, 1) = tris.col(2);
		edges.block(2 * trisCount, 0, trisCount, 1) = tris.col(2);
		edges.block(2 * trisCount, 1, trisCount, 1) = tris.col(0);

		// FOR DEBUG:
		dispMatBlock<int>(edges, 0, 10, 0, 1);
		dispMatBlock<int>(edges, 3 * trisCount - 10, 3 * trisCount - 1, 0, 1);

		std::vector<Eigen::Triplet<int>> elems;
		elems.reserve(2 * (versCount + trisCount + 100));				// euler formula: E = V + S - 2 + 2*g;
		for (unsigned i = 0; i < edgesCount; ++i)
			elems.push_back(Eigen::Triplet<int>(edges(i, 0), edges(i, 1), 1));

		adjSM_eCount.resize(versCount, versCount);
		adjSM_eCount.reserve(elems.size());
		adjSM_eCount.setFromTriplets(elems.begin(), elems.end());
		adjSM_eCount.makeCompressed();

		std::cout << "索引为0的顶点的1邻域顶点：" << std::endl;
		traverseSTL(adjList[0], disp<int>);
		std::cout << std::endl << std::endl;

		std::cout << "索引为1的顶点的入边（即邻接矩阵中第1列中所有不为零的元素下标）：" << std::endl;
		for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM, 0); it; ++it)
		{
			std::cout << "value == " << it.value() << std::endl;
			std::cout << "row == " << it.row() << std::endl;			 // row index
			std::cout << "col == " << it.col() << std::endl;			 // col index (here it is equal to k)
			std::cout << std::endl;
		} 

		std::cout << "索引为2的顶点的出边（即邻接矩阵中第2行中所有不为零的元素下标）：" << std::endl;
		for (int i = 0; i < adjSM_eCount.outerSize(); ++i)			// 对列的遍历：
		{
			for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM_eCount, i); it; ++it)		// 当前列的列内迭代器；
			{
				if (2 == it.row())
				{
					std::cout << "value == " << it.value() << std::endl;
					std::cout << "row == " << it.row() << std::endl;			 // row index
					std::cout << "col == " << it.col() << std::endl;			 // col index (here it is equal to k)
				}
			}
		}

		// 判断是否有重复三角片——若有重复三角片，则存在重复有向边：
		bool hasDupTri = false;;
		for (int i = 0; i < adjSM_eCount.outerSize(); ++i)			// 对列的遍历：
		{
			if (hasDupTri)
				break;
			for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM_eCount, i); it; ++it)		// 当前列的列内迭代器；
			{
				if (it.value() > 1)
				{
					std::cout << "存在重复三角片" << std::endl;
					hasDupTri = true;
					break;
				}
			}
		}
		if (!hasDupTri)
			std::cout << "没有重复三角片" << std::endl;

		std::cout << "finished." << std::endl;
	}


	// 基本的图算法；
	void test1() 
	{
		MatrixXd vers;
		MatrixXi tris;
		igl::readOBJ("E:/材料/roundSurf.obj", vers, tris);
 
		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);

		// dfs:
		Eigen::VectorXi disCoveredIdx, bfsTreeVec, dfsTreeVec, closedIdx;
		size_t startIdx = 165;			// 接近圆心的顶点；
		igl::dfs(adjList, startIdx, disCoveredIdx, dfsTreeVec, closedIdx);
		objWriteTreePath("E:/dfsTree.obj", bfsTreeVec, vers);
		Eigen::VectorXi retCrev = closedIdx.reverse();
		auto flag = (retCrev == disCoveredIdx);

		std::vector<int> vec1, vec2, vec3;
		vec1 = vec2Vec(disCoveredIdx);
		vec2 = vec2Vec(closedIdx);
		vec3 = vec2Vec(retCrev);

		// bfs:
		igl::bfs(adjList, startIdx, disCoveredIdx, bfsTreeVec);
		objWriteTreePath("E:/bfsTree.obj", bfsTreeVec, vers);

		// myDfs:
		std::vector<bool> visited(vers.rows(), false);
		std::list<int> discoveredVersIdx;				// 已被访问的顶点：
		std::list<int> closedVersIdx;

		std::function<void(const int, const int)> myDfs = [&](const int index, const int parentIdx)
		{
			// 递归终止；
			if (visited[index])			
				return;

			visited[index] = true;
			discoveredVersIdx.push_back(index);
			const std::vector<int>& adjVersIdx = adjList[index];

			// 递归递推；
			for (const auto& adjIdx : adjVersIdx)			
				myDfs(adjIdx, index);

			closedVersIdx.push_back(index);
		};


 
		std::cout << "finished." << std::endl;
	}


	// prime, dijkstra
	void test2() 
	{
		// dijkstra
		MatrixXd vers;
		MatrixXi tris;
		igl::readOBJ("E:/材料/roundSurf.obj", vers, tris);

		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);

		Eigen::VectorXd min_distance;				// 图中所有顶点到指定顶点的最短路径长度；
		Eigen::VectorXi mst;						// 最小生成树；
		int verIdx0 = 165;						// 设起点为接近圆心的顶点；

		// 不指定重点时， igl::dijkstra()返回图中以verIdx0为起点的最小生成树；
		int retIdx = igl::dijkstra(vers, adjList, verIdx0, std::set<int>{}, min_distance, mst);
		std::cout << "retIdx == " << retIdx << std::endl;
		objWriteTreePath("E:/mst.obj", mst, vers);


		std::cout << "finished." << std::endl;
	}
 
}


// 空间划分
namespace IGL_SPACE_PARTITION
{
	// aabb树实现的BVH(层次包围体)
	void test0()
	{
		Eigen::MatrixXd vers, vers0, minDisVers;
		Eigen::MatrixXi tris;
		Eigen::VectorXd minSqrDis;
		Eigen::VectorXi minDisIdx;
		igl::readOBJ("E:/材料/tooth.obj", vers, tris);

		vers0.resize(2, 3);
		vers0.row(0) = vers.row(0);
		vers0.row(1) = vers.row(99);

		// igl::point_mesh_squared_distance()——使用BVH求顶点到网格最小距离；
		igl::point_mesh_squared_distance(vers0, vers, tris, minSqrDis, minDisIdx, minDisVers);
		igl::writeOBJ("E:/inputMesh.obj", vers, tris);
		igl::writeOBJ("E:/vers0.obj", vers0, Eigen::MatrixXi{});
		igl::writeOBJ("E:/minDisVers.obj", minDisVers, Eigen::MatrixXi{});

		//// 生成网格的BVH对象：
		//igl::AABB<double, 3> bvh;

		

		std::cout << "finished." << std::endl;
	}

}


// IGL实现的基础三角网格处理算法；
namespace IGL_BASIC_PMP 
{
	// 边std::pair<int, int>的自定义哈希函数；
	struct edgeHash
	{
		bool operator()(const std::pair<int, int>& edge) const
		{
			return (std::hash<int>()(edge.first) + std::hash<int>()(edge.second));
		}
	};

	// 边std::pair<int, int>的等价比较器；
	struct edgeComparator 
	{
		bool operator()(const std::pair<int, int>& edge1, const std::pair<int, int>& edge2) const
		{
			return (edge1.first == edge2.first && edge1.second == edge2.second);
		}
	};


	// 非流形网格中的边-三角片邻接关系：
	void test1() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris, edges;
		igl::readOBJ("E:/材料/meshArranged_hole.obj", vers, tris);
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		unsigned edgesCount = 3 * trisCount;		// 非流形网格中edgesCount比真实的有向边数要大，因为对于重复的同一有向边进行了重复计数；

		// edges == [ea; eb; ec] == [vbIdxes, vcIdxes; vcIdxes, vaIdxes; vaIdxes, vbIdxes];
		/*
			若网格是流形网格，则edges列表里每条边都是unique的；
			若存在非流形边，则非流形边在edges里会重复存储，实际边数量也会少于edges行数；
			可以说使用这种representation的话，非流形边有不止一个索引，取决包含该非流形边的三角片数量；
		*/
		edges = Eigen::MatrixXi::Zero(edgesCount, 2);
		Eigen::MatrixXi vaIdxes = tris.col(0);
		Eigen::MatrixXi vbIdxes = tris.col(1);
		Eigen::MatrixXi vcIdxes = tris.col(2);
		edges.block(0, 0, trisCount, 1) = vbIdxes;
		edges.block(trisCount, 0, trisCount, 1) = vcIdxes;
		edges.block(trisCount * 2, 0, trisCount, 1) = vaIdxes;
		edges.block(0, 1, trisCount, 1) = vcIdxes;
		edges.block(trisCount, 1, trisCount, 1) = vaIdxes;
		edges.block(trisCount * 2, 1, trisCount, 1) = vbIdxes;

		// 1. 邻接矩阵：
		std::vector<Eigen::Triplet<int>> smElems;
		smElems.reserve(edgesCount);
		for (unsigned i = 0; i < edgesCount; ++i)
			smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});

		Eigen::SparseMatrix<int> adjSM_eCount;
		adjSM_eCount.resize(versCount, versCount);
		adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());		// 权重为该有向边重复的次数；

		// 2. 边索引-三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
		std::vector<int> etInfo(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

#if __cplusplus < 202002L
		// 3. 确定非流形边；
		std::unordered_set<std::pair<int, int>, edgeHash, edgeComparator> edgesNM;
#else
		// lambda——std::pair<int, int>表示的边的自定义哈希函数；
		auto edgeHashLamb = [](const std::pair<int, int>& edge)->std::size_t
		{
			return (std::hash<int>()(edge.first) + std::hash<int>()(edge.second));
		};

		// lambda——std::pair<int, int>表示的边的等价比较器；
		auto edgeComLamb = [](const std::pair<int, int>& edge1, const std::pair<int, int>& edge2)->bool
		{
			return (edge1.first == edge2.first && edge1.second == edge2.second);
		};

		// 3. 确定非流形边；
		std::unordered_set<std::pair<int, int>, decltype(edgeHashLamb), decltype(edgeComLamb)> edgesNM;				// 需要C++20;非流形边；使用unordered_set去重；
#endif
		 for (unsigned i = 0; i < adjSM_eCount.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_eCount, i); iter; ++iter)	// 第i列的列内迭代器；
			{
				if (iter.value() > 1)
				{
					int num1 = iter.row();
					int num2 = iter.col();
					edgesNM.insert({num1, num2});
				}
			}
		}

		
#if __cplusplus < 202002L
		 //		非流形有向边-其索引的键值对，一个非流形边对应着多个边索引；
		 std::unordered_map<std::pair<int, int>, std::vector<int>, edgeHash, edgeComparator> edgesNMmap;

		 //		非流形有向边-其所在的三角片的索引的键值对，一个非流形边对应着多个三角片索引；
		 std::unordered_map<std::pair<int, int>, std::vector<int>, edgeHash, edgeComparator> etNMmap;
#else
		 //		非流形有向边-其索引的键值对，一个非流形边对应着多个边索引；
		 std::unordered_map<std::pair<int, int>, std::vector<int>, decltype(edgeHashLamb), decltype(edgeComLamb)> edgesNMmap;
		 //		非流形有向边-其所在的三角片的索引的键值对，一个非流形边对应着多个三角片索引；
		 std::unordered_map<std::pair<int, int>, std::vector<int>, decltype(edgeHashLamb), decltype(edgeComLamb)> etNMmap;
#endif

		for (auto& nmEdge : edgesNM)
		{
			for (int i = 0; i < edgesCount; ++i)
			{
				if (edges(i, 0) == nmEdge.first && edges(i, 1) == nmEdge.second)
				{
					auto retPair = edgesNMmap.insert({nmEdge, std::vector<int>{i} });
					if (!retPair.second)			// 若插入失败，则说明已有此键；
					{
						auto iter = edgesNMmap.find(nmEdge);
						iter->second.push_back(i);
					}
				}
			}
		}
		for (const auto& pair : edgesNMmap)
		{
			auto copyPair = pair;
			for (auto& index : copyPair.second)
				index = etInfo[index];
			etNMmap.insert(copyPair);
		}

		// 4. 找出所有非流形边所在三角片：
		std::vector<int> trisIdxNM;
		for (auto& pair : etNMmap)
		{
			for (auto& index : pair.second)
				trisIdxNM.push_back(index);
		}
		Eigen::MatrixXi trisNM(Eigen::MatrixXi::Zero(trisIdxNM.size(), 3));
		for (int i = 0; i < trisIdxNM.size(); ++i)
			trisNM.row(i) = tris.row(trisIdxNM[i]);

		// for debug: 打印trisNM:
		objWriteMeshMat("E:/trisNM.obj", vers, trisNM);
		std::cout << "finished." << std::endl;

		// 5. 确定流形三角片、非流形三角片的邻接关系：


	}

	using ttTuple = std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>;

	bool buildAdjacency(const Eigen::MatrixXi& tris, Eigen::MatrixXi& ttAdj_nmEdge, \
		std::vector<ttTuple>& ttAdj_nmnEdge, std::vector<ttTuple>& ttAdj_nmnOppEdge)
	{
		/*
			bool buildAdjacency(
						const Eigen::MatrixXi& tris,												输入的三角片数据
						Eigen::MatrixXi& ttAdj_nmEdge,										三角片的非边缘流形有向边邻接的三角片索引；
						std::vector<ttTuple>& ttAdj_nmnEdge,							三角片的非流形有向边所在的三角片索引；
						std::vector<ttTuple>& ttAdj_nmnOppEdge					三角片的非流形有向边邻接的三角片索引（即对边所在的三角片的索引）；
						)

		*/
		const unsigned trisCount = tris.rows();
		const unsigned edgesCount = 3 * trisCount;
		const unsigned versCount = tris.maxCoeff() + 1;

		// 1. 求顶点邻接关系：

		// edges == [ea; eb; ec] == [vbIdxes, vcIdxes; vcIdxes, vaIdxes; vaIdxes, vbIdxes];
		/*
			若网格是流形网格，则edges列表里每条边都是unique的；
			若存在非流形边，则非流形边在edges里会重复存储，实际边数量也会少于edges行数；
			可以说使用这种representation的话，非流形边有不止一个索引，取决包含该非流形边的三角片数量；
		*/

		// 三角片三条边的索引：teIdx == [eaIdx, ebIdx, ecIdx] == [(0: trisCount-1)', (trisCount: 2*trisCount-1)', (2*trisCount, 3*trisCount-1)'];
		/*
			teIdx(i, j) = trisCount *j + i;
		*/
		Eigen::MatrixXi edges = Eigen::MatrixXi::Zero(edgesCount, 2);
		Eigen::MatrixXi vaIdxes = tris.col(0);
		Eigen::MatrixXi vbIdxes = tris.col(1);
		Eigen::MatrixXi vcIdxes = tris.col(2);
		edges.block(0, 0, trisCount, 1) = vbIdxes;
		edges.block(trisCount, 0, trisCount, 1) = vcIdxes;
		edges.block(trisCount * 2, 0, trisCount, 1) = vaIdxes;
		edges.block(0, 1, trisCount, 1) = vcIdxes;
		edges.block(trisCount, 1, trisCount, 1) = vaIdxes;
		edges.block(trisCount * 2, 1, trisCount, 1) = vbIdxes;

		std::vector<Eigen::Triplet<int>> smElems, smElems_weighted;
		smElems.reserve(edgesCount);
		smElems_weighted.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
		{
			smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});
			smElems_weighted.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), i});
		}

		Eigen::SparseMatrix<int> adjSM_eCount, adjSM_weighted;
		adjSM_eCount.resize(versCount, versCount);
		adjSM_weighted.resize(versCount, versCount);
		adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());		// 权重为该有向边重复的次数；
		adjSM_weighted.setFromTriplets(smElems_weighted.begin(), smElems_weighted.end());		// 权重为该有向边重复的次数(非流形边的值无效)
		Eigen::SparseMatrix<int> adjSM_weighted_opp = adjSM_weighted.transpose();

		Eigen::SparseMatrix<int> adjSM = adjSM_eCount;		// 有向边邻接矩阵；
		for (unsigned i = 0; i < adjSM.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM, i); iter; ++iter)
			{
				if (iter.value() > 0)
					adjSM.coeffRef(iter.row(), iter.col()) = 1;
			}
		}


		// 2. 确定所有非边缘的流形有向边：
		Eigen::SparseMatrix<int> adjSM_eCount_ND = adjSM_eCount + Eigen::SparseMatrix<int>(adjSM_eCount.transpose());		// 权重为无向边ij关联的三角片数量；

		//		非边缘流形无向边：
		Eigen::SparseMatrix<int> adjSM_MNnonBdry_ND = adjSM_eCount_ND;
		for (unsigned i = 0; i < adjSM_MNnonBdry_ND.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_MNnonBdry_ND, i); iter; ++iter)
			{
				if (iter.value() > 0)
				{
					if (iter.value() == 2)
						adjSM_MNnonBdry_ND.coeffRef(iter.row(), iter.col()) = 1;
					else
						adjSM_MNnonBdry_ND.coeffRef(iter.row(), iter.col()) = 0;
				}
			}
		}

		//		非边缘流形有向边邻接矩阵
		Eigen::SparseMatrix<int> adjSM_MNnonBdry = adjSM + adjSM_MNnonBdry_ND;		// adjSM & adjSM_MNnonBdry_ND
		adjSM_MNnonBdry.prune([](const Eigen::Index& row, const Eigen::Index& col, const float& value)->bool
			{
				if (2 == value)
					return true;
				else
					return false;
			});
		adjSM_MNnonBdry /= 2;

		unsigned edgesCount_MNnonBdry = adjSM_MNnonBdry.sum();
		Eigen::MatrixXi edges_MNnonBdry(edgesCount_MNnonBdry, 2);		// 非边缘流形有向边
		std::vector<int> edgesIdx_MNnonBdry;												// 非边缘流形有向边的索引；
		edgesIdx_MNnonBdry.reserve(edgesCount_MNnonBdry);
		unsigned index = 0;
		for (unsigned i = 0; i < adjSM_MNnonBdry.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_MNnonBdry, i); iter; ++iter)
			{
				edges_MNnonBdry(index, 0) = iter.row();
				edges_MNnonBdry(index, 1) = iter.col();
				edgesIdx_MNnonBdry.push_back(adjSM_weighted.coeffRef(iter.row(), iter.col()));
				index++;
			}
		}

		Eigen::SparseMatrix<int> adjSM_MNnonBdry_opp = adjSM_MNnonBdry.transpose();
		Eigen::MatrixXi edges_MNnonBdry_opp(edgesCount_MNnonBdry, 2);		// 非边缘流形有向边的对边；
		std::vector<int> edgesIdx_MNnonBdry_opp;								// 非边缘流形有向边的对边的索引；
		edgesIdx_MNnonBdry_opp.reserve(edgesCount_MNnonBdry);
		index = 0;
		for (unsigned i = 0; i < adjSM_MNnonBdry_opp.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_MNnonBdry_opp, i); iter; ++iter)
			{
				edges_MNnonBdry_opp(index, 0) = iter.row();
				edges_MNnonBdry_opp(index, 1) = iter.col();
				edgesIdx_MNnonBdry_opp.push_back(adjSM_weighted_opp.coeffRef(iter.row(), iter.col()));
				index++;
			}
		}

		// 3. 流形边-三角片邻接关系，三角片邻接关系：
		std::vector<int> etInfo(edgesCount);				// 边索引 - 三角片索引映射表；etInfo(i)是索引为i的边所在的三角片的索引；
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;

		Eigen::VectorXi etAdj_mnEdge(-Eigen::VectorXi::Ones(edgesCount));	// 边所在三角片的索引，非流形边或边缘边写-1；
		for (unsigned i = 0; i < edgesIdx_MNnonBdry.size(); ++i)
		{
			const int& edgeIdx = edgesIdx_MNnonBdry[i];
			const int& edgesIdxOpp = edgesIdx_MNnonBdry_opp[i];
			etAdj_mnEdge(edgeIdx) = etInfo[edgesIdxOpp];
		}

		//			三角片邻接矩阵，ttAdj_mnEdge(i, :)是索引为i的三角片三条边邻接的三个三角片，若其中有非流形边或边缘边则写为-1；
		ttAdj_nmEdge = Eigen::Map<Eigen::MatrixXi>(etAdj_mnEdge.data(), trisCount, 3);

		// 4. 非流形边（关联两个三角片的有向边）信息
		index = 0;
		std::vector<int> edgeIdx_nmn;
		edgeIdx_nmn.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			if (2 == adjSM_eCount.coeffRef(edges(i, 0), edges(i, 1)))
				edgeIdx_nmn.push_back(i);
		edgeIdx_nmn.shrink_to_fit();


		//			找出同一条非流形有向边对应的多个边索引：map<pair表示的非流形边数据， 边索引vector>
		std::unordered_map<std::pair<int, int>, std::vector<int>, edgeHash, edgeComparator> edges_nmn_map;
		for (const auto& eIdx : edgeIdx_nmn)
		{
			std::pair<int, int> edge{ edges(eIdx, 0), edges(eIdx, 1) };
			auto retPair = edges_nmn_map.insert({ edge , std::vector<int>{eIdx} });
			if (!retPair.second)		// 若插入失败，则说明已有此键；
			{
				auto iter = edges_nmn_map.find(edge);
				iter->second.push_back(eIdx);
			}
		}

		//			map<非流形边Idx, 该边的对边对应的多个边索引的vector>
		std::unordered_map<int, std::vector<int>> edgeIdx_nmn_map;
		for (const auto& eIdx : edgeIdx_nmn)
		{
			std::pair<int, int> edge{ edges(eIdx, 0), edges(eIdx, 1) };
			const std::vector<int>& commonEidxes = edges_nmn_map.find(edge)->second;
			edgeIdx_nmn_map.insert({ eIdx, commonEidxes });
		}
		std::unordered_map<int, std::vector<int>> edgeIdx_nmn_opp_map;
		for (const auto& pair : edgeIdx_nmn_map)
		{
			int eIdx = pair.first;
			int vaIdx = edges(eIdx, 0);
			int vbIdx = edges(eIdx, 1);
			int eOppIdx = -1;
			for (int i = 0; i < edgesCount; ++i)
				if (edges(i, 0) == vbIdx && edges(i, 1) == vaIdx)
				{
					eOppIdx = i;
					break;
				}
			auto iter = edgeIdx_nmn_map.find(eOppIdx);
			edgeIdx_nmn_opp_map.insert({eIdx, iter->second});
		}


		// 5. 含有非流形边的三角片的邻接关系：
		ttAdj_nmnEdge.resize(trisCount);
		for (const auto& pair : edgeIdx_nmn_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;

			std::vector<int> trisIdx_nmn = pair.second;
			for (auto& index : trisIdx_nmn)
				index = etInfo[index];

			switch (col)
			{
			case 0:	std::get<0>(ttAdj_nmnEdge[row]) = trisIdx_nmn;
				break;
			case 1:	std::get<1>(ttAdj_nmnEdge[row]) = trisIdx_nmn;
				break;
			case 2:	std::get<2>(ttAdj_nmnEdge[row]) = trisIdx_nmn;
				break;
			default:
				return false;
			}
		}

		ttAdj_nmnOppEdge.resize(trisCount);
		for (const auto& pair : edgeIdx_nmn_opp_map)
		{
			int eIdx = pair.first;
			int row = eIdx % trisCount;
			int col = eIdx / trisCount;

			std::vector<int> trisIdx_nmn_opp = pair.second;
			for (auto& index : trisIdx_nmn_opp)
				index = etInfo[index];

			switch (col)
			{
			case 0:	std::get<0>(ttAdj_nmnOppEdge[row]) = trisIdx_nmn_opp;
				break;
			case 1:	std::get<1>(ttAdj_nmnOppEdge[row]) = trisIdx_nmn_opp;
				break;
			case 2:	std::get<2>(ttAdj_nmnOppEdge[row]) = trisIdx_nmn_opp;
				break;
			default:
				return false;
			}
		}


		// for debug
		if(0)
		{
			Eigen::MatrixXd vers;
			Eigen::MatrixXi tris_useless;
			igl::readOBJ("E:/材料/meshArranged_hole.obj", vers, tris_useless);

			for (int i = 0; i < ttAdj_nmnEdge.size(); ++i)
			{
				const auto& tuple = ttAdj_nmnEdge[i];
				const auto& tupleOpp = ttAdj_nmnOppEdge[i];
				if (!std::get<0>(tuple).empty() || !std::get<1>(tuple).empty() || !std::get<2>(tuple).empty())
				{
					const std::vector<int>& tmpVec = (std::get<0>(tuple).empty() && std::get<1>(tuple).empty()) ? std::get<2>(tuple)\
								: (std::get<0>(tuple).empty() ? std::get<1>(tuple): std::get<0>(tuple));
					Eigen::MatrixXi tris0 = tris.row(i);
					Eigen::MatrixXi tris00(tmpVec.size(), 3);
					for (int k = 0; k < tmpVec.size(); ++k)
						tris00.row(k) = tris.row(tmpVec[k]);
					igl::writeOBJ("E:/tris0.obj", vers, tris0);
					igl::writeOBJ("E:/tris00.obj", vers, tris00);

					const std::vector<int>& tmpVecOpp = (std::get<0>(tupleOpp).empty() && std::get<1>(tupleOpp).empty()) ? std::get<2>(tupleOpp)\
						: (std::get<0>(tupleOpp).empty() ? std::get<1>(tupleOpp) : std::get<0>(tupleOpp));
					Eigen::MatrixXi tris000(tmpVecOpp.size(), 3);
					for (int k = 0; k < tmpVecOpp.size(); ++k)
						tris000.row(k) = tris.row(tmpVecOpp[k]);
					igl::writeOBJ("E:/tris000.obj", vers, tris000);

					break;
				}
			}
		}
 
		return true;
	}


	// test buildAdjacency();
	void test2() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris, edges;
		igl::readOBJ("E:/材料/meshArranged_hole.obj", vers, tris);

		Eigen::MatrixXi ttAdj_nmEdge;
		std::vector<ttTuple> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		bool retFlag = buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);
	}


	// boolean——select tris:
	enum BOOLEAN_TYPE
	{
		UNION = 0,
		DIFF = 1,
	};


	bool booleanSelectTris(const Eigen::MatrixXd& vers, const Eigen::MatrixXi& tris, const Eigen::VectorXi& trisLabel, \
		const int startIdx, const BOOLEAN_TYPE type, Eigen::MatrixXd& versOut, Eigen::MatrixXi& trisOut)
	{
		/*
			bool booleanSelectTris(
					const Eigen::MatrixXd& vers,
					const Eigen::MatrixXi& tris,
					const Eigen::VectorXi& trisLabel, 
					const int startIdx, 
					const BOOLEAN_TYPE type, 
					Eigen::MatrixXd& versOut, 
					Eigen::MatrixXi& trisOut
					)


		
		*/

		// 1. 计算三角片邻接关系
		Eigen::MatrixXi ttAdj_nmEdge;
		std::vector<ttTuple> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);
		unsigned trisCount = tris.rows();
		std::vector<bool> visited(trisCount, false);

		// 2. 使用区域生长的形式遍历三角片：
		std::deque<int> triIdxDeq;
		triIdxDeq.push_back(startIdx);
		visited[startIdx] = true;
		const int startLabel = trisLabel[startIdx];

		while (!triIdxDeq.empty())					// 区域生长的循环：
		{
			int currentTriIdx = triIdxDeq.front();
			int currentLabel = trisLabel[currentTriIdx];
			triIdxDeq.pop_front();				
			
			// 对当前三角片邻接的三个三角片的遍历：
			for (int i = 0; i < 3; ++i)
			{
				int nbrTriIdx = ttAdj_nmEdge(currentTriIdx, i);
				
				// wf1. 若当前边为非流形边：
				if (-1 == nbrTriIdx)
				{
					std::vector<std::vector<int>*> vecPtr(3), vecOppPtr(3);
					vecPtr[0] = &std::get<0>(ttAdj_nmnEdge[currentTriIdx]);
					vecPtr[1] = &std::get<1>(ttAdj_nmnEdge[currentTriIdx]);
					vecPtr[2] = &std::get<2>(ttAdj_nmnEdge[currentTriIdx]);
					vecOppPtr[0] = &std::get<0>(ttAdj_nmnOppEdge[currentTriIdx]);
					vecOppPtr[1] = &std::get<1>(ttAdj_nmnOppEdge[currentTriIdx]);
					vecOppPtr[2] = &std::get<2>(ttAdj_nmnOppEdge[currentTriIdx]);

					const std::vector<int>& relaTrisIdx = *vecPtr[i];				// 当前非流形边所在的所有三角片（正常的话应该为两个）；
					const std::vector<int>& relaOppTrisIdx = *vecOppPtr[i];	// 当前非流形边的所有对面三角片（正常的话应该为两个）；
					switch (type)
					{
					case BOOLEAN_TYPE::DIFF:
						{
							assert(2 == relaTrisIdx.size() && "Exceptional non-manifold edge detected!");
							assert(2 == relaOppTrisIdx.size() && "Exceptional non-manifold edge detected!");

							// 若nbrTriIdx两个三角片中至少有一个和startIdx三角片标签相同；往标签跳变的三角片上扩散；
							if (trisLabel[relaTrisIdx[0]] == startLabel || trisLabel[relaTrisIdx[1]] == startLabel)
							{
								nbrTriIdx = (trisLabel[relaTrisIdx[0]] == currentLabel ? relaTrisIdx[1] : relaTrisIdx[0]);
								assert(trisLabel[relaTrisIdx[0]] != trisLabel[relaTrisIdx[1]] && "Exceptional non-manifold edge detected!");
							}
							else    // 若nbrTriIdx两个三角片都和startIdx三角片标签不同；往标签跳变的对面三角片上扩散；
							{
								nbrTriIdx = (trisLabel[relaOppTrisIdx[0]] == currentLabel ? relaOppTrisIdx[1] : relaOppTrisIdx[0]);
								assert(trisLabel[relaOppTrisIdx[0]] != trisLabel[relaOppTrisIdx[1]] && "Exceptional non-manifold edge detected!");
							}
							break;
						}

					case BOOLEAN_TYPE::UNION:
						{
							// 若nbrTriIdx两个三角片都和startIdx三角片标签不同；往标签跳变的对面三角片上扩散；
							nbrTriIdx = (trisLabel[relaOppTrisIdx[0]] == currentLabel ? relaOppTrisIdx[1] : relaOppTrisIdx[0]);
							assert(trisLabel[relaOppTrisIdx[0]] != trisLabel[relaOppTrisIdx[1]] && "Exceptional non-manifold edge detected!");
							break;
						}

					default:
						assert("invalid BOOLEAN_TYPE!");
					}
				}

				// wf2. 标记访问之后的三角片，队尾加入新的三角片
				if (!visited[nbrTriIdx])
				{
					visited[nbrTriIdx] = true;
					triIdxDeq.push_back(nbrTriIdx);
				}
			}
		}

		// 3. 找出区域生长过程中被标记的三角片，区分是否和起始三角片标记相同；
		std::vector<int> selectedTrisIdx, selectedTrisIdx1, selectedTrisIdx2;
		Eigen::MatrixXi selectedTris1, selectedTris2, selectedTris;
		selectedTrisIdx.reserve(trisCount);
		selectedTrisIdx1.reserve(trisCount);
		selectedTrisIdx2.reserve(trisCount);
		for (int i = 0; i < trisCount; ++i)
		{
			if (visited[i])
			{
				selectedTrisIdx.push_back(i);
				if ((trisLabel[i] == startLabel))
					selectedTrisIdx1.push_back(i);
				else
					selectedTrisIdx2.push_back(i);
			}
		}
		selectedTrisIdx.shrink_to_fit();
		selectedTrisIdx1.shrink_to_fit();
		selectedTrisIdx2.shrink_to_fit();

		// 4. 选取三角片生成结果网格：
		switch (type)
		{
		case BOOLEAN_TYPE::DIFF:
			{
				// d1. 提取被标记的三角片，区分是否和起始三角片标记相同；
				subFromIdxVec(selectedTris1, tris, selectedTrisIdx1);
				subFromIdxVec(selectedTris2, tris, selectedTrisIdx2);

				// d2. 与起始三角片标记不同，则flip该三角片
				Eigen::VectorXi tmpVec = selectedTris2.col(2);
				selectedTris2.col(2) = selectedTris2.col(1);
				selectedTris2.col(1) = tmpVec;
				trisOut = selectedTris1;
				matInsertRows(trisOut, selectedTris2);

				// d3. 去除其他顶点，修正三角片中的索引；
				std::set<int> newTrisIdxSet;
				std::vector<int> newOldIdxInfo;
				std::vector<int> oldNewIdxInfo(trisCount, -1);
				int index = 0;
				int* intPtr = trisOut.data();
				for (unsigned i = 0; i < trisOut.size(); ++i)
					newTrisIdxSet.insert(*(intPtr+ i));
				newOldIdxInfo.insert(newOldIdxInfo.end(), newTrisIdxSet.begin(), newTrisIdxSet.end());

				for (const auto& oldIdx : newOldIdxInfo)
					oldNewIdxInfo[oldIdx] = index++;
				
				subFromIdxVec(versOut, vers, newOldIdxInfo);
				intPtr = trisOut.data();
				for (unsigned i = 0; i < trisOut.size(); ++i)
				{
					int& currentIdx = *(intPtr + i);
					currentIdx = oldNewIdxInfo[currentIdx];
				}

				break;
			}

		case BOOLEAN_TYPE::UNION:
			{
				// u1. 提取被标记的三角片：
				subFromIdxVec(selectedTris, tris, selectedTrisIdx);
				trisOut = selectedTris;

				// u2. 去除其他顶点，修正三角片中的索引；
				std::set<int> newTrisIdxSet;
				std::vector<int> newOldIdxInfo;
				std::vector<int> oldNewIdxInfo(trisCount, -1);
				int index = 0;
				int* intPtr = trisOut.data();
				for (unsigned i = 0; i < trisOut.size(); ++i)
					newTrisIdxSet.insert(*(intPtr + i));
				newOldIdxInfo.insert(newOldIdxInfo.end(), newTrisIdxSet.begin(), newTrisIdxSet.end());

				for (const auto& oldIdx : newOldIdxInfo)
					oldNewIdxInfo[oldIdx] = index++;

				subFromIdxVec(versOut, vers, newOldIdxInfo);
				intPtr = trisOut.data();
				for (unsigned i = 0; i < trisOut.size(); ++i)
				{
					int& currentIdx = *(intPtr + i);
					currentIdx = oldNewIdxInfo[currentIdx];
				}

				break;
			}

		default:
			break;
		}

		return true;
	}


	// test boolean:
	void test3() 
	{
		Eigen::MatrixXd vers, versDiff, versUnion;
		Eigen::MatrixXi tris, trisDiff, trisUnion;
		Eigen::VectorXi trisLabel;

		igl::readOBJ("E:/材料/meshArranged.obj", vers, tris);
		int trisCount = tris.rows();

		// 读取trisLabel
		std::string str(1024, '\0');
		trisLabel.resize(trisCount);
		std::ifstream file("E:/outputlabel.txt");
		for(int i = 0; i<trisCount; ++i)
		{
			file.getline(&str[0], 1024);
			trisLabel(i) = std::stoi(str);
		}

		int startIdx = 0;
		booleanSelectTris(vers, tris, trisLabel, startIdx, BOOLEAN_TYPE::UNION, versUnion, trisUnion);
		booleanSelectTris(vers, tris, trisLabel, startIdx, BOOLEAN_TYPE::DIFF, versDiff, trisDiff);

		igl::writeOBJ("E:/diffResult.obj", versDiff, trisDiff);
		igl::writeOBJ("E:/unionResult.obj", versUnion, trisUnion);

		std::cout << "finished." << std::endl;
	}


	// marching cubes:

	template <
		typename Derivedres,
		typename DerivedV>
		IGL_INLINE void grid(
			const Eigen::MatrixBase<Derivedres>& res,
			Eigen::PlainObjectBase<DerivedV>& GV)
	{
		using namespace Eigen;
		typedef typename DerivedV::Scalar Scalar;
		GV.resize(res.array().prod(), res.size());
		const auto lerp = [&res](const Scalar di, const int d)->Scalar {return di / (Scalar)(res(d) - 1); };

		int gi = 0;
		Derivedres sub;
		sub.resizeLike(res);
		sub.setConstant(0);

		for (int gi = 0; gi < GV.rows(); gi++)
		{
			// omg, I'm implementing addition...
			for (int c = 0; c < res.size() - 1; c++)
			{
				if (sub(c) >= res(c))
				{
					sub(c) = 0;
					// roll over
					sub(c + 1)++;
				}
			}

			for (int c = 0; c < res.size(); c++)
				GV(gi, c) = lerp(sub(c), c);

			sub(0)++;
		}
	}

	//		输入AABB，生成栅格：
	template <typename Scalar, typename DerivedV,	typename DerivedI>
	bool genGrids(const Eigen::AlignedBox<Scalar, 3>& box, const int largestCount,	const int pad_count,\
			Eigen::PlainObjectBase<DerivedV>& gridCenters, Eigen::PlainObjectBase<DerivedI>& gridCounts)
	{
		/*bool genGrids(																		成功返回true
				const Eigen::AlignedBox<Scalar, 3>&box,						输入的AABB对象；
				const int largestCount,													跨度最大的那个维度(xyz中的一个)的栅格数；		
				const int pad_count,														超出包围盒边界的栅格数；
				Eigen::PlainObjectBase<DerivedV>&gridCenters,			栅格中心；
				Eigen::PlainObjectBase<DerivedI>&gridCounts				xyz三个维度上栅格的数量；
				)
		*/
		using namespace Eigen;
		using namespace std;

		// 1. 计算包围盒对角线向量中的最大分量；
		gridCounts.resize(1, 3);
		typename DerivedV::Index maxCompIdx = -1;            // 包围盒box的对角线向量中最大的分量的索引；0,1,2分别对应xyz分量；
		box.diagonal().maxCoeff(&maxCompIdx);
		const Scalar maxComp = box.diagonal()(maxCompIdx);          // 包围盒box的对角线向量中最大的分量；
		assert(largestCount > (pad_count * 2 + 1) && "largestCount should be > 2*pad_count+1");


		// 2. 计算xyz三个维度上栅格的个数gridCounts
		const Scalar largestCount0 = largestCount - 2 * pad_count;
		gridCounts(maxCompIdx) = largestCount0;
		for (int i = 0; i < 3; i++)
		{
			if (i != maxCompIdx)
				gridCounts(i) = std::ceil(largestCount0 * (box.diagonal()(i)) / maxComp);
		}
		gridCounts.array() += 2 * pad_count;

		// 3. 计算gridCenters;
		grid(gridCounts, gridCenters);            // 计算中心在原点的gridCenters;

		/*
			 A *    p/largestCount  + B = min
			 A * (1-p/largestCount) + B = max
			 B = min - A * p/largestCount
			 A * (1-p/largestCount) + min - A * p/largestCount = max
			 A * (1-p/largestCount) - A * p/largestCount = max-min
			 A * (1-2p/largestCount) = max-min
			 A  = (max-min)/(1-2p/largestCount)
		*/

		auto tmp = gridCounts.transpose().template cast<Scalar>().array() - 1.;
		const Array<Scalar, 3, 1> ps = (Scalar)(pad_count) / tmp;
		const Array<Scalar, 3, 1> A = box.diagonal().array() / (1.0 - 2. * ps);

		/*
			// This would result in an "anamorphic",  but perfectly fit grid:
			const Array<Scalar, 3, 1> B = box.min().array() - A.array()*ps;
			gridCenters.array().rowwise() *= A.transpose();
			gridCenters.array().rowwise() += B.transpose();
			 Instead scale by largest factor and move to match center
		*/
		typename Array<Scalar, 3, 1>::Index ai = -1;
		Scalar a = A.maxCoeff(&ai);
		const Array<Scalar, 1, 3> ratio = a * (gridCounts.template cast<Scalar>().array() - 1.0) / (Scalar)(gridCounts(ai) - 1.0);
		gridCenters.array().rowwise() *= ratio;
		const Eigen::Matrix<Scalar, 1, 3> offset = (box.center().transpose() - gridCenters.colwise().mean()).eval();
		gridCenters.rowwise() += offset;

		return true;
	}


	bool marchingCubes(const Eigen::VectorXd& SDF, const Eigen::MatrixXd& gridCenters, const Eigen::RowVector3i& gridCounts, \
				const double selectedSDF, Eigen::MatrixXd& versResult_SDF, Eigen::MatrixXi& trisResult_SDF)
	{


		return true;
	}


	void test4()
	{
		// 0. 解析SDFGen.exe生成的.sdf距离场数据文件：
		std::vector<int> stepCounts(3);				// xyz三个维度上栅格数
		Eigen::RowVector3d gridsOri;					// 栅格原点：
		Eigen::VectorXd SDF;
		const char* sdfFilePath = "E:/inputMesh.sdf";
		double SDFstep = IGL_BASIC::parseSDF(stepCounts, gridsOri, SDF, sdfFilePath);

		// 1. 生成栅格：
		Eigen::RowVector3i gridCounts;
		Eigen::MatrixXd gridCenters;
		Eigen::RowVector3d minp = gridsOri - SDFstep * Eigen::RowVector3d(stepCounts[0] / 2.0, stepCounts[1] / 2.0, stepCounts[2] / 2.0);
		Eigen::RowVector3d maxp = gridsOri + SDFstep * Eigen::RowVector3d(stepCounts[0] / 2.0, stepCounts[1] / 2.0, stepCounts[2] / 2.0);
		Eigen::AlignedBox<double, 3> box(minp, maxp);
		genGrids(box, std::max({ stepCounts[0], stepCounts[1], stepCounts[2] }), 0, gridCenters, gridCounts);

		// gridCenters是所有栅格中点坐标的矩阵，每行都是一个中点坐标；存储优先级是x, y, z
		Eigen::MatrixXd gridCenters0;
		Eigen::RowVector3i gridCounts0{stepCounts[0], stepCounts[1], stepCounts[2]};
		Eigen::RowVector3d arrow = box.diagonal();
		gridCenters0.resize(gridCounts0(0) * gridCounts(1) * gridCounts(2), 3);
		Eigen::VectorXd xPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(0), minp(0), maxp(0));
		Eigen::VectorXd yPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(1), minp(1), maxp(1));
		Eigen::VectorXd zPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(2), minp(2), maxp(2));

		Eigen::MatrixXd tmpVec0, tmpVec1, tmpVec2;
		kron(tmpVec0, VectorXd::Ones(gridCounts(1) * gridCounts(2)), xPeriod);
 

		// for debug
		dispVec<double>(xPeriod);
		dispVecSeg<double>(Eigen::VectorXd{ tmpVec0 }, 0, 150);
		dispMatBlock<double>(gridCenters, 0, 144, 0, 2);

		// 2. marching cubes算法生成最终曲面：
		tiktok& tt = tiktok::getInstance();
		MatrixXd versResult_SDF, versResults_signs;
		MatrixXi trisResult_SDF, trisResults_signs;
		double selectedSDF = -1.;
		tt.start();
		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), selectedSDF, versResult_SDF, trisResult_SDF);
		tt.endCout("Elapsed time of igl::marching_cubes() is ");

		igl::writeOBJ("E:/shrinkedMesh.obj", versResult_SDF, trisResult_SDF);

	}
}
