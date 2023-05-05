#include "igl_study.h"
#define DATA_PATH "./data/"

static igl::opengl::glfw::Viewer viewer;				// libigl中的基于glfw的显示窗口；
static std::mutex g_mutex;


/// /////////////////////////////////////////////////////////////////////////////////////////// DEBUG 接口

static std::string g_debugPath = "E:/";

static void debugDisp()			// 递归终止
{						//		递归终止设为无参或者一个参数的情形都可以。
	std::cout << std::endl;
	return;
}

template <typename T, typename... Types>
static void debugDisp(const T& firstArg, const Types&... args)
{
	std::cout << firstArg << " ";
	debugDisp(args...);
}


template <typename T, int M, int N>
static void dispData(const Eigen::Matrix<T, M, N>& m)
{
	auto dataPtr = m.data();
	unsigned elemsCount = m.size();

	for (unsigned i = 0; i < elemsCount; ++i)
		std::cout << dataPtr[i] << ", ";

	std::cout << std::endl;
}


template <typename Derived>
static void dispData(const Eigen::PlainObjectBase<Derived>& m)
{
	int m0 = m.RowsAtCompileTime;
	int n0 = m.ColsAtCompileTime;

	auto dataPtr = m.data();
	unsigned elemsCount = m.size();

	for (unsigned i = 0; i < elemsCount; ++i)
		std::cout << dataPtr[i] << ", ";

	std::cout << std::endl;
}


template <typename Derived>
static void dispElem(const Eigen::MatrixBase<Derived>& m)
{
	const Derived& mm = m.derived();
	std::cout << mm(1, 1) << std::endl;
}


template<typename DerivedV>
static void debugWriteVers(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteVerticesMat(path, vers);
}


template<typename T>
static void debugWriteMesh(const char* name, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteMeshMat(path, vers, tris);
}
 

template<typename DerivedV>
static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteEdgesMat(path, edges, vers);
}



////////////////////////////////////////////////////////////////////////////// libigl基本功能
namespace IGL_BASIC
{
	Eigen::MatrixXd vers, newVers, normals;
	Eigen::MatrixXi tris;
	Eigen::SparseMatrix<double> L;


	// igl中基础的矩阵算子：
	void test00() 
	{
		Eigen::MatrixXi m1(3, 3);
		Eigen::MatrixXi m11;
		Eigen::VectorXi vec1(2);

		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		dispMat(m1);

		// slice()——使用索引矩阵、索引向量提取矩阵中的元素；
		/*
			  template <typename MatX,  typename DerivedR,  typename MatY>
			  IGL_INLINE void slice(
						const MatX& X,																	输入矩阵
						const Eigen::DenseBase<DerivedR> & R,							索引向量
						const int dim,															维度——1-行，2-列；
						MatY& Y																				输出矩阵；
					)
		*/

		vec1 << 0, 2;
		igl::slice(m1, vec1, 1, m11);				// 提取第0, 2行；
		dispMat(m11);

		igl::slice(m1, vec1, 2, m11);				// 提取第0, 2列；
		dispMat(m11);

		Eigen::VectorXi vec2(2);
		vec2 << 1, 2;
		igl::slice(m1, vec1, vec2, m11);		// 等价于：m11 == m1(0: 2, 1: 2);
		dispMat(m11);


		std::cout << "finished." << std::endl;
	}


	// igl中计算基础几何属性的接口：
	void test000()
	{
		bool retFlag = igl::readOBJ("E:/材料/tooth.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		Eigen::MatrixXd verNorms, triNorms, conerNorms;			// 各种法向；每行是一个法向量；
		Eigen::MatrixXi edges;									// 半边
		Eigen::MatrixXi uEdges;								// 无向边，既是undirected的又是unique的；
		Eigen::VectorXi edgeUeInfo;						// 半边在无向边中对应的索引；edgeUeInfo(i)是索引为i的半边对应的无向边的索引；
		Eigen::MatrixXi UeTrisInfo;							// 无向边-关联三角片映射表——第i行的元素：索引为i的无向边关联的三角片的索引；
		Eigen::MatrixXi UeCornersInfo;					// 无向边-对顶点映射表——第i行的元素：索引为i的无向边所对的两个顶点的索引；若为非边缘流形边则有两个对顶点；
		std::vector<std::vector<int>> UeEdgeInfo;
		Eigen::VectorXi uEC, uEE;
		Eigen::VectorXd dbArea;								// 网格中每个三角片的两倍面积；
		std::vector<int> nbrVersIdx, nbrTrisIdx;

		// per_vertex_normals
		igl::per_vertex_normals(vers, tris, verNorms);
		igl::per_face_normals(vers, tris, triNorms);
		igl::per_corner_normals(vers, tris, conerNorms);

		// unique_edge_map()——计算网格的半边、无向边、及其对应关系；
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, UeEdgeInfo);
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, uEC, uEE);

		// edge_flaps()——计算无向边关联的三角片；
		igl::edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
		objWriteEdgesMat("E:/uEdge0.obj", Eigen::RowVector2i(uEdges.row(0)), vers);
		Eigen::RowVector3i tri00 = tris.row(UeTrisInfo(0, 0));
		Eigen::RowVector3i tri01 = tris.row(UeTrisInfo(0, 1));
		Eigen::MatrixXi tmpTris;
		matInsertRows<int, 3>(tmpTris, tri00);
		matInsertRows<int, 3>(tmpTris, tri01);
		objWriteMeshMat("E:/uEdge0relaTris.obj", vers, tmpTris);

		Eigen::MatrixXd oppVers(2, 3);
		oppVers.row(0) = vers.row(tri00(UeCornersInfo(0, 0)));
		oppVers.row(1) = vers.row(tri01(UeCornersInfo(0, 1)));
		objWriteVerticesMat("E:/oppVers.obj", oppVers);


		// circulation()——寻找输入边的端点关联的三角片，及其1领域顶点；
		igl::circulation(0, true, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx, nbrTrisIdx);

		Eigen::MatrixXd vers0;
		Eigen::MatrixXi tris0;
		objWriteVerticesMat("E:/edgeHead0.obj", Eigen::MatrixXd{ vers.row(uEdges(0, 0)) });
		subFromIdxVec(vers0, vers, nbrVersIdx);
		objWriteVerticesMat("E:/vers0.obj", vers0);
		subFromIdxVec(tris0, tris, nbrTrisIdx);
		igl::writeOBJ("E:/tris0.obj", vers, tris0);

		igl::circulation(0, false, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx, nbrTrisIdx);

		Eigen::MatrixXd vers1;
		Eigen::MatrixXi tris1;
		objWriteVerticesMat("E:/edgeTail0.obj", Eigen::MatrixXd{ vers.row(uEdges(0, 1)) });
		subFromIdxVec(vers1, vers, nbrVersIdx);
		objWriteVerticesMat("E:/vers1.obj", vers1);
		subFromIdxVec(tris1, tris, nbrTrisIdx);
		igl::writeOBJ("E:/tris1.obj", vers, tris1);


		// igl::doubleArea()——计算输入网格每个三角片面积的两倍：
		igl::doublearea(vers, tris, dbArea);
		std::vector<double> dbAreaVec = vec2Vec(dbArea);
		std::sort(dbAreaVec.begin(), dbAreaVec.end());
		std::cout << "finished." << std::endl;
	}


	// 文件IO
	void test0()
	{
		//// readOBJ(), writeObj()——OBJ文件的IO，有多个重载，至少有路径、点云、三角片这三个数据。
		igl::readOBJ("./data/bunny.obj", vers, tris);
		igl::writeOBJ("./data/bunny_export.obj", vers, tris);
		igl::writeOBJ("./data/bunnyVers.obj", vers, Eigen::MatrixXi{});			// 只要点云不需要三角片的话，传入空矩阵；
	
		vers.resize(0, 0);
		tris.resize(0, 0);
		igl::readOBJ("E:/材料/jawMeshDense.obj", vers, tris);
		igl::writeOFF("E:/材料/jawMeshDense.off", vers, tris);

		// 读取stl文件；
		std::string fileName{"E:/材料/jawMeshSimplified"};
		vers.resize(0, 0);
		tris.resize(0, 0);
		std::ifstream fileIn((fileName + std::string{ ".stl" }).c_str(), std::ios::binary);			// stl文件是二进制文件；
		igl::readSTL(fileIn, vers, tris, normals);
		fileIn.close();
		igl::writeOBJ((fileName + std::string{".obj"}).c_str(), vers, tris);

		// 读.mesh文件
		fileName = "E:/材料/bunny";
		vers.resize(0, 0);
		tris.resize(0, 0);
		Eigen::MatrixXi tets;
		igl::readMESH((fileName + std::string{ ".mesh" }).c_str(), vers, tets, tris);
		igl::writeOBJ((fileName + std::string{ ".obj" }).c_str(), vers, tris);
		igl::writeOBJ((fileName + std::string{ "_tets.obj" }).c_str(), vers, tets);
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

		Eigen::VectorXd U;
		igl::readDMAT("./data/cheburashka-scalar.dmat", U);		// 一系列的函数值

		Eigen::SparseMatrix<double> G;			// 梯度算子
		igl::grad(vers, tris, G);

		// 计算函数值U的梯度
		Eigen::MatrixXd GU = Map<const Eigen::MatrixXd>((G * U).eval().data(), tris.rows(), 3);

		// 函数值U的梯度GU的模长
		const Eigen::VectorXd GU_mag = GU.rowwise().norm();

		viewer.data().set_mesh(vers, tris);
		viewer.data().set_data(U);

		// Average edge length divided by average gradient (for scaling)
		const double max_size = igl::avg_edge_length(vers, tris) / GU_mag.mean();

		// 每个三角片重心上画一根指示线，方向为梯度方向。 
		Eigen::MatrixXd BC;
		igl::barycenter(vers, tris, BC);
		const Eigen::RowVector3d black(0, 0, 0);
		viewer.data().add_edges(BC, BC + max_size * GU, black);
		viewer.data().show_lines = false;	  // 隐藏网格线

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


	// test5. 计算符号距离场、marching cubes提取等值面网格；
	void test5()
	{
		Eigen::MatrixXi tris;
		Eigen::MatrixXd vers;
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


		// ！！！注：libigl中的计算SDF接口igl::signed_distance生成的距离场貌似没有SDFgen中的好用，有时候会出错；
		Eigen::VectorXd SDF, signValues;
		{
			Eigen::VectorXi I;                     // useless
			Eigen::MatrixXd C, N;              // useless

			// 2. 计算符号距离场
			tt.start();
			igl::signed_distance(gridCenters, vers, tris, igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, SDF, I, C, N);
			tt.endCout("Elapsed time of igl::signed_distance() is ");

			// 3. 符号距离场改写为符号场——网格内为-1，网格面上为0，外面为1：
			signValues = SDF;
			std::for_each(signValues.data(), signValues.data() + signValues.size(), [](double& b)\
			{
				b = (b > 0 ? 1 : (b < 0 ? -1 : 0));
			});
		}

		// 4. marching cubes算法生成最终曲面：
		Eigen::MatrixXd versResult_SDF, versResults_signs;
		Eigen::MatrixXi trisResult_SDF, trisResults_signs;
		double selectedSDF = -1.;
		tt.start();
		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), selectedSDF, versResult_SDF, trisResult_SDF);
		tt.endCout("Elapsed time of igl::marching_cubes() is ");
		
		// igl::marching_cubes(signValues, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), 0, versResults_signs, trisResults_signs);
		igl::writeOBJ("E:/shrinkedMesh.obj", versResult_SDF, trisResult_SDF);


		std::cout << "finished." << std::endl;
	}


	// 读取SDFGen.exe生成的.sdf距离场数据，使用igl::marching_cubes()提取等值面网格：

	//解析.sdf文本文件：
	double parseSDF(std::vector<int>& stepCounts, Eigen::RowVector3d& gridsOri, Eigen::VectorXd& SDF, const char* filePath)
	{
		/*
			double parseSDF(												返回距离场的采样步长；
					std::vector<int>& stepCounts,					xyz三个方向上的步数
					Eigen::RowVector3d& gridsOri,					距离场空间的原点
					Eigen::VectorXd& SDF,								距离场数据
					const char* filePath										SDF文件目录
					)
		
		*/
		double SDFstep = -1;
		stepCounts.resize(3);
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
				if (ch >= '0' && ch <= '9' || ch == '.' || ch == '+' || ch == '-')
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
				if (ch >= '0' && ch <= '9' || ch == '.' || ch == '+' || ch == '-')
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


	// test55. 读取SDFgen生成的符号距离场数据，使用marching cubes生成网格：
	void test55() 
	{
		// 0. 解析SDFGen.exe生成的.sdf距离场数据文件：
		std::vector<int> stepCounts(3);				// xyz三个维度上栅格数
		Eigen::RowVector3d gridsOri;					// 栅格原点：
		Eigen::VectorXd SDF;
		const char* sdfFilePath = "E:/jawMeshUnionRepair1.sdf";
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
		Eigen::MatrixXd versResult_SDF, versResults_signs;
		Eigen::MatrixXi trisResult_SDF, trisResults_signs;
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
 
		igl::writeOBJ("E:/426.obj", Eigen::MatrixXd{ Eigen::RowVector3d(4, 2, 6) }, Eigen::MatrixXi{});

		std::cout << "finished." << std::endl;
	}

}


////////////////////////////////////////////////////////////////////////////// libigl中的微分几何相关
namespace IGL_DIF_GEO 
{
	Eigen::MatrixXd vers_g, newVers_g, normals_g;
	Eigen::MatrixXi tris_g;
	Eigen::SparseMatrix<double> L_g;
	float deltaLB_g;
	unsigned smoothLoopCount_g = 0;


	// test0: 质量矩阵和LB算子
	void test0() 
	{
		Eigen::MatrixXd vers_g;
		Eigen::MatrixXi tris_g;
		Eigen::SparseMatrix<double> L_g, M;

		igl::readOBJ("E:/材料/tooth.obj", vers_g, tris_g);
		igl::cotmatrix(vers_g, tris_g, L_g);
		igl::massmatrix(vers_g, tris_g, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);

		dispSpMat(M, 0, M.rows() - 1, 10);

		std::cout << "finished." << std::endl;
	}
 

	// test1: 使用Laplacian光顺网格
	
	// laplace光顺键盘事件：使用laplacian光顺网格
	bool key_down_laplacian(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)
	{
		std::string outputName;
		outputName.resize(512);
		switch (key)
		{
		case 'r':

		case 'R':							// 复位程序
			newVers_g = vers_g;
			smoothLoopCount_g = 0;
			break;

		case 's':

		case 'S':			
			sprintf_s(&outputName[0], 512, "E:/meshSmoothLB_step%d.obj", smoothLoopCount_g);
			igl::writeOBJ(outputName, newVers_g, tris_g);
			std::cout << "step " << smoothLoopCount_g << " result saved." << std::endl;
			break;

		case ' ':								// 空格键，执行一次laplace光顺
		{
			// 1. 重新计算质量矩阵
			Eigen::SparseMatrix<double> mass;
			igl::massmatrix(newVers_g, tris_g, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

			// 2. 解线性方程组 (mass - delta*L_g) * newVers_g = mass * newVers_g
			const auto& S = (mass - deltaLB_g * L_g);
			Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
			assert(solver.info() == Eigen::Success);
			newVers_g = solver.solve(mass * newVers_g).eval();

			std::cout << "smoothLoopCount_g == " << (++smoothLoopCount_g) << std::endl;

			break;
		}

		default:
			return false;
		}

		viewer.data().set_vertices(newVers_g);
		viewer.data().compute_normals();
		viewer.core().align_camera_center(newVers_g, tris_g);
		return true;
	};


	void test1()
	{
		Eigen::MatrixXd norms;
		Eigen::MatrixXd colors;

		igl::readOBJ("E:/gum.obj", vers_g, tris_g);
		deltaLB_g = 0.5;

		newVers_g = vers_g;

		const unsigned versCount = vers_g.rows();
		const unsigned trisCount = tris_g.rows();

		// 1.a 直接构造laplacian——Compute Laplace-Beltrami operator: 
		igl::cotmatrix(vers_g, tris_g, L_g);

		// 1.b 测试分步构造laplacian
		{
			Eigen::SparseMatrix<double> Gradient, L2;
			Eigen::VectorXd dblA;									  // 每个三角片面积的两倍
			igl::grad(vers_g, tris_g, Gradient);      // 离散梯度

			// Diagonal per-triangle "mass matrix"			
			igl::doublearea(vers_g, tris_g, dblA);           

			// Place areas along diagonal  #dim times
			const auto& T = 1. * (dblA.replicate(3, 1) * 0.5).asDiagonal();

			L2 = -Gradient.transpose() * T * Gradient;         // discrete Dirichelet energy Hessian 离散狄利克雷能量海塞矩阵
			std::cout << "两种方法得到的laplacian的差的范数：" << std::endl;
			std::cout << "(L2 - L_g).norm() == " << (L2 - L_g).norm() << std::endl;
		}

		// 2. 根据原始的法向量，使用伪色
		igl::per_vertex_normals(vers_g, tris_g, norms);
		colors = norms.rowwise().normalized().array() * 0.5 + 0.5;

		// 3. viewr填充初始数据
		newVers_g = vers_g;
		viewer.data().set_mesh(newVers_g, tris_g);
		viewer.data().set_colors(colors);
		viewer.callback_key_down = key_down_laplacian;

		// 4. 运行
		std::cout << "Press [space] to smooth." << std::endl;
		std::cout << "Press [r] to reset." << std::endl;

		viewer.launch();
	}

}


/////////////////////////////////////////////////////////////////////////////// 图算法
namespace IGL_GRAPH 
{
	// 图数据结构的转换：
	void test0() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("E:/材料/cube重复三角片.obj", vers, tris);

		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		
		// 由面片信息得到图的邻接表、邻接矩阵
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);				// 有向边存在则元素为1，否则为0；

		// 手动计算邻接矩阵：
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		unsigned edgesCount = 3 * trisCount;
		Eigen::SparseMatrix<int> adjSM_eCount;				// 有向边邻接矩阵，权重为该有向边的数量；

		Eigen::MatrixXi edges = Eigen::MatrixXi::Zero(edgesCount, 2);
		edges.block(0, 0, trisCount, 1) = tris.col(0);
		edges.block(0, 1, trisCount, 1) = tris.col(1);
		edges.block(trisCount, 0, trisCount, 1) = tris.col(1);
		edges.block(trisCount, 1, trisCount, 1) = tris.col(2);
		edges.block(2 * trisCount, 0, trisCount, 1) = tris.col(2);
		edges.block(2 * trisCount, 1, trisCount, 1) = tris.col(0);

		// FOR DEBUG:
		dispMatBlock(edges, 0, 10, 0, 1);
		dispMatBlock(edges, 3 * trisCount - 10, 3 * trisCount - 1, 0, 1);

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
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("E:/材料/roundSurf.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);
 
		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);

		size_t startIdx = 165;			// 接近圆心的顶点；
		Eigen::RowVector3d startVer = vers.row(startIdx);
		objWriteVerticesMat("E:/startVer.obj", startVer);

		// dfs:
		Eigen::VectorXi disCoveredIdx, bfsTreeVec, dfsTreeVec, closedIdx;
		igl::dfs(adjList, startIdx, disCoveredIdx, dfsTreeVec, closedIdx);
		objWriteTreePath("E:/dfsTree.obj", dfsTreeVec, vers);
		Eigen::VectorXi retCrev = closedIdx.reverse();
		auto cornerState = (retCrev == disCoveredIdx);

		std::vector<int> vec1, vec2, vec3;
		vec1 = vec2Vec(disCoveredIdx);
		vec2 = vec2Vec(closedIdx);
		vec3 = vec2Vec(retCrev);

		// bfs:
		igl::bfs(adjList, startIdx, disCoveredIdx, bfsTreeVec);
		objWriteTreePath("E:/bfsTree.obj", bfsTreeVec, vers);


		// 自己写的DFS算法——目前有问题
		std::vector<bool> visited(vers.rows(), false);
		std::list<int> discoveredVersIdx;						// 已被访问的顶点：
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


		myDfs(startIdx, startIdx);
		std::vector<int> tmpVec;
		tmpVec.insert(tmpVec.end(), closedVersIdx.begin(), closedVersIdx.end());
		Eigen::VectorXi myDfsTreeVec = vec2Vec(tmpVec);
		objWriteTreePath("E:/myDfsTree.obj", myDfsTreeVec, vers);

 
		std::cout << "finished." << std::endl;
	}


	// prime, dijkstra
	void test2() 
	{
		// dijkstra
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
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


/////////////////////////////////////////////////////////////////////////////// 空间划分
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


/////////////////////////////////////////////////////////////////////////////// IGL实现的基础三角网格处理算法；
namespace IGL_BASIC_PMP 
{
	// 非流形网格中的边-三角片邻接关系：
	void test1() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris, edges;
		igl::readOBJ("E:/材料/meshArranged.obj", vers, tris);
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


		// 3. 确定非流形边；
		std::unordered_set<std::int64_t> edgesNMN;					// 边用边编码表示
		 for (unsigned i = 0; i < adjSM_eCount.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_eCount, i); iter; ++iter)	// 第i列的列内迭代器；
			{
				if (iter.value() > 1)
				{
					int vaIdx = iter.row();
					int vbIdx = iter.col();
					edgesNMN.insert(encodeEdge(vaIdx, vbIdx));
				}
			}
		}

		 //		非流形有向边-其索引的键值对，一个非流形边对应着多个边索引；
		 std::unordered_map<std::int64_t, std::vector<int>> edgesNMMmap;

		 //		非流形有向边-其所在的三角片的索引的键值对，一个非流形边对应着多个三角片索引；
		 std::unordered_map<std::int64_t, std::vector<int>> etNMNmap;


		for (auto& eCode : edgesNMN)
		{
			for (int i = 0; i < edgesCount; ++i)
			{
				std::pair<int, int> retPair = decodeEdge(eCode);
				int vaIdx = retPair.first;
				int vbIdx = retPair.second;
				if (edges(i, 0) == vaIdx && edges(i, 1) == vbIdx)
				{
					auto retPair = edgesNMMmap.insert({eCode, std::vector<int>{i} });
					if (!retPair.second)			// 若插入失败，则说明已有此键；
					{
						auto iter = edgesNMMmap.find(eCode);
						iter->second.push_back(i);
					}
				}
			}
		}
		for (const auto& pair : edgesNMMmap)
		{
			auto copyPair = pair;
			for (auto& index : copyPair.second)
				index = etInfo[index];
			etNMNmap.insert(copyPair);
		}
	}


	// test buildAdjacency();
	void test2() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris, edges;
		igl::readOBJ("E:/材料/jawMeshArranged.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();

		Eigen::MatrixXi ttAdj_mnEdge;
		std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		bool retFlag = buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris);

		// 输出包含非流形边的三角片：
		std::vector<int> nmnTrisIdx;
		for (const auto& vec : ttAdj_nmnEdge)
		{
			nmnTrisIdx.insert(nmnTrisIdx.end(), vec[0].begin(), vec[0].end());
			nmnTrisIdx.insert(nmnTrisIdx.end(), vec[1].begin(), vec[1].end());
			nmnTrisIdx.insert(nmnTrisIdx.end(), vec[2].begin(), vec[2].end());
		}
		std::unique(nmnTrisIdx.begin(), nmnTrisIdx.end());
		Eigen::MatrixXi nmnTris;
		subFromIdxVec(nmnTris, tris, nmnTrisIdx);
		objWriteMeshMat("E:/nmnTris.obj", vers, nmnTris);

		std::vector<int> mnTrisIdx;
		std::sort(nmnTrisIdx.begin(), nmnTrisIdx.end());
		std::unordered_set<int> tmpSet;
		for (int i = 0; i < trisCount; ++i)
			tmpSet.insert(i);
		for (const auto& index : nmnTrisIdx)
			tmpSet.erase(index);

		mnTrisIdx.insert(mnTrisIdx.end(), tmpSet.begin(), tmpSet.end());
		Eigen::MatrixXi mnTris;
		subFromIdxVec(mnTris, tris, mnTrisIdx);
		objWriteMeshMat("E:/mnTris.obj", vers, mnTris);

		Eigen::MatrixXi nmnEdges;
		nonManifoldUEs(nmnEdges, tris);
		objWriteEdgesMat("E:/nmnEdges.obj", nmnEdges, vers); 

		std::cout << "finished." << std::endl;
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
		Eigen::MatrixXi ttAdj_mnEdge;
		std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris);
		unsigned trisCount = tris.rows();
		std::vector<bool> visited(trisCount, false);			// 三角片是否被访问的标记；

		// 2. 使用区域生长的形式遍历三角片：
		std::deque<int> triIdxDeq;						// 工作队列；
		triIdxDeq.push_back(startIdx);
		visited[startIdx] = true;
		const int startLabel = trisLabel[startIdx];
		while (!triIdxDeq.empty())					// 区域生长的循环：
		{
			int currentTriIdx = triIdxDeq.front();
			int currentLabel = trisLabel[currentTriIdx];
			triIdxDeq.pop_front();				
			
			// 对当前三角片三条边的遍历：
			for (int i = 0; i < 3; ++i)
			{
				// wf0. 若当前边为流形边，则在ttAdj_mnEdge中寻找其对面；
				int nbrTriIdx = ttAdj_mnEdge(currentTriIdx, i);
				
				// wf1. 若当前边为非流形边，根据布尔操作不同类型选取不同的对面；
				if (-1 == nbrTriIdx)
				{
					std::vector<std::vector<int>*> vecPtr(3), tmpPtrVec(3);
					vecPtr[0] = &ttAdj_nmnEdge[currentTriIdx][0];
					vecPtr[1] = &ttAdj_nmnEdge[currentTriIdx][1];
					vecPtr[2] = &ttAdj_nmnEdge[currentTriIdx][2];
					tmpPtrVec[0] = &ttAdj_nmnOppEdge[currentTriIdx][0];
					tmpPtrVec[1] = &ttAdj_nmnOppEdge[currentTriIdx][1];
					tmpPtrVec[2] = &ttAdj_nmnOppEdge[currentTriIdx][2];

					const std::vector<int>& relaTrisIdx = *vecPtr[i];				// 当前非流形边所在的所有三角片（正常的话应该为两个）；
					const std::vector<int>& relaOppTrisIdx = *tmpPtrVec[i];	// 当前非流形边的所有对面三角片（正常的话应该为两个）；
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


	// test solve self-intersection:
	enum class VISIT_FLAG
	{
		NOT_VISITED = 0,
		VISITED = 1,
		COLLECTED = 2,
	};


	// 33. 测试buildAdjacency
	void test33()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		igl::readOBJ("E:/材料/meshArrangeResult.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();

		Eigen::MatrixXi ttAdj_mnEdge;
		std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		bool retFlag = buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris);

		// 输出包含非流形边的三角片：
		std::vector<int> nmnTrisIdx;
		for (const auto& vec : ttAdj_nmnEdge)
		{
			nmnTrisIdx.insert(nmnTrisIdx.end(), vec[0].begin(), vec[0].end());
			nmnTrisIdx.insert(nmnTrisIdx.end(), vec[1].begin(), vec[1].end());
			nmnTrisIdx.insert(nmnTrisIdx.end(), vec[2].begin(), vec[2].end()); 
		}
		std::unique(nmnTrisIdx.begin(), nmnTrisIdx.end());
		Eigen::MatrixXi nmnTris;
		subFromIdxVec(nmnTris, tris, nmnTrisIdx);
		objWriteMeshMat("E:/nmnTris.obj", vers, nmnTris);

		Eigen::MatrixXi mnTris;				// 不包含非流形边的三角片：
		std::set<int> tmpSet;
		traverseMatrix(ttAdj_mnEdge, [&](const int num)
			{
				if (num > 0)
					tmpSet.insert(num);
			});

 
		std::cout << "finished." << std::endl;
	}


	// marching cubes:
	template <typename Derivedres, typename DerivedV>
		IGL_INLINE void grid(const Eigen::MatrixBase<Derivedres>& res,
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
			for (int i = 0; i < res.size() - 1; i++)
			{
				if (sub(i) >= res(i))
				{
					sub(i) = 0;
					// roll over
					sub(i + 1)++;
				}
			}

			for (int i = 0; i < res.size(); i++)
				GV(gi, i) = lerp(sub(i), i);

			sub(0)++;
		}
	}


	//		输入AABB，生成栅格；——from libigl
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
			if (i != maxCompIdx)
				gridCounts(i) = std::ceil(largestCount0 * (box.diagonal()(i)) / maxComp);
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
 
 
	//		marchingCubes算法中生成一个cube的三角片：
	template <typename DerivedGV, typename Scalar, typename Index, typename ScalarV, typename IndexF>
	void handleCube(const DerivedGV& gridCenters, const Eigen::Matrix<Scalar, 8, 1>& cornerSDF, \
		const Eigen::Matrix<Index, 8, 1>& cornerIdx, const Scalar& isovalue, \
		Eigen::Matrix<ScalarV, Eigen::Dynamic, Eigen::Dynamic>& versResult, Index& curVersCount, \
		Eigen::Matrix<IndexF, Eigen::Dynamic, Eigen::Dynamic>& trisResult, Index& curTrisCount, \
		std::unordered_map<int64_t, int>& edgeIsctMap)
	{
		/*
			const DerivedGV& gridCenters,															栅格数据
			const Eigen::Matrix<Scalar, 8, 1>& cornerSDF,									当前立方体八个顶点的SDF值
			const Eigen::Matrix<Index, 8, 1>& cornerIdx,									当前立方体八个顶点在栅格中的索引；
			const Scalar& isovalue,																		需要提取的等值面的SDF值
			Eigen::PlainObjectBase<DerivedV>& versResult,								输出网格的顶点
			Index& curVersCount,																			当前累计生成的输出网格顶点数
			Eigen::PlainObjectBase<DerivedF>& trisResult,									输出网格的三角片
			Index& curTrisCount,																			当前累计生成的输出网格三角片数
			std::unordered_map<int64_t, int>& edgeIsctMap								边编码-边交点索引的哈希表；

		*/

		Eigen::Matrix<Index, 12, 1> isctVerIdxes;		// 立方体边上的交点的绝对索引——是在最终输出网格中的索引；
		int cornerState = 0;											// 立方体顶点状态编码；256种情形；

		// 生成无向边编码
		const auto genMCedgeCode = [](int32_t vaIdx, int32_t vbIdx)
		{
			if (vaIdx > vbIdx)
				std::swap(vaIdx, vbIdx);
			std::int64_t edgeCode = 0;
			edgeCode |= vaIdx;
			edgeCode |= static_cast<std::int64_t>(vbIdx) << 32;
			return edgeCode;
		};

		// 1. 计算当前立方体的顶点状态编码，即8个顶点在等值面的内外状态1；
		for (int i = 0; i < 8; i++)
			if (cornerSDF(i) > isovalue)
				cornerState |= 1 << i;

		// 2. 确定当前立方体中和等值面相交的边；
		int edgeState = MC_TABLES::edgeStateCodes[cornerState];		// 立方体顶点状态编码映射为相交边编码；
		if (edgeState == 0)
			return;															// 表示当前立方体整体都在等值面外部或内部，没有交点；

		// 3. 确定等值面和当前立方体的边的交点； Find the point of intersection of the surface with each edge. Then find the normal to the surface at those points
		for (int i = 0; i < 12; i++)						// 对立方体所有边的遍历
		{
			if (edgeState & (1 << i))					// 若等值面和当前边相交：
			{
				int vaIdxRela = MC_TABLES::cubeEdges[i][0];			// 当前边两端点的相对索引；
				int vbIdxRela = MC_TABLES::cubeEdges[i][1];

				// 生成边上的顶点：
				int vaIdx = cornerIdx(vaIdxRela);				// 当前边两端点的绝对索引，是栅格中的顶点索引；
				int vbIdx = cornerIdx(vbIdxRela);
				std::int64_t edgeCode = genMCedgeCode(vaIdx, vbIdx);
				const auto iter = edgeIsctMap.find(edgeCode);

				if (iter == edgeIsctMap.end())								// 若当前边交点未生成；
				{
					if (curVersCount == versResult.rows())
						versResult.conservativeResize(versResult.rows() * 2 + 1, versResult.cols());

					// 插值生成新的顶点：find crossing point assuming linear interpolation along edges
					const Scalar& SDFa = cornerSDF(vaIdxRela);			// 当前边两端点的SDF值；
					const Scalar& SDFb = cornerSDF(vbIdxRela);
					const Scalar delta = SDFb - SDFa;
					Scalar t = (isovalue - SDFa) / delta;
					versResult.row(curVersCount) = (gridCenters.row(vaIdx) + t * (gridCenters.row(vbIdx) - gridCenters.row(vaIdx))).array().cast<ScalarV>();

					isctVerIdxes[i] = curVersCount;
					edgeIsctMap[edgeCode] = isctVerIdxes[i];
					curVersCount++;
				}
				else                                                                             // 若当前边交点已生成；
					isctVerIdxes[i] = iter->second;

				assert(isctVerIdxes[i] >= 0);
				assert(isctVerIdxes[i] < curVersCount);
			}
		}


		// 4. 生成当前立方体中的三角片，一个立方体中最多生成5个三角片；
		for (int i = 0; i < 5; i++)
		{
			if (MC_TABLES::cubeTriangles[cornerState][3 * i] < 0)
				break;

			if (curTrisCount == trisResult.rows())
				trisResult.conservativeResize(trisResult.rows() * 2 + 1, trisResult.cols());

			// 新增三角片数据中的顶点索引，是相对索引；
			int vaIdxRela = MC_TABLES::cubeTriangles[cornerState][3 * i + 0];
			int vbIdxRela = MC_TABLES::cubeTriangles[cornerState][3 * i + 1];
			int vcIdxRela = MC_TABLES::cubeTriangles[cornerState][3 * i + 2];

			assert(isctVerIdxes[vaIdxRela] >= 0);
			assert(isctVerIdxes[vbIdxRela] >= 0);
			assert(isctVerIdxes[vcIdxRela] >= 0);

			// 相对索引转换为绝对索引，插入新增三角片
			trisResult.row(curTrisCount) << isctVerIdxes[vaIdxRela], isctVerIdxes[vbIdxRela], isctVerIdxes[vcIdxRela];
			curTrisCount++;
		}

	}


	//		marchingCubes()——！！！还需要在最后去除重复顶点，太短的边。
	template <typename DerivedS, typename DerivedGV, typename ScalarV, typename IndexF>
	bool marchingCubes(Eigen::Matrix<ScalarV, Eigen::Dynamic, Eigen::Dynamic>& versResult, \
		Eigen::Matrix<IndexF, Eigen::Dynamic, Eigen::Dynamic>& trisResult, \
		const Eigen::MatrixBase<DerivedS>& scalarFied, const Eigen::MatrixBase<DerivedGV>& gridCenters, \
		const unsigned nx, const unsigned ny, const unsigned nz,
		const typename DerivedS::Scalar isovalue)
	{
		/*
			DerivedS				距离场数据的类型
			DerivedGV			栅格数据的类型
			ScalarV					顶点坐标的数据类型；
			IndexF					面片中的顶点索引的数据类型；
		

			const Eigen::MatrixBase<DerivedS>& scalarFied,							符号距离场数据
			const Eigen::MatrixBase<DerivedGV>& gridCenters,					栅格数据
			const unsigned nx,																			x方向上栅格个数
			const unsigned ny,
			const unsigned nz,
			const typename DerivedS::Scalar isovalue,										需要提取的水平集的SDF值；
			Eigen::PlainObjectBase<DerivedV>& versResult,							输出网格顶点
			Eigen::PlainObjectBase<DerivedF>& trisResult								输出网格三角片
		*/

		typedef typename DerivedS::Scalar Scalar;

		// lambda——栅格的三维索引映射到一维索引：
		const auto getGridIdx = [&nx, &ny, &nz](const int& x, const int& y, const int& z)->unsigned
		{
			return x + nx * (y + ny * (z));
		};

		const unsigned cornerIdxOffset[8] = { 0, 1, 1 + nx, nx, nx * ny, 1 + nx * ny, 1 + nx + nx * ny, nx + nx * ny };	// 立方体八个顶点的索引偏移量；
		std::unordered_map<int64_t, int> edgeIsctMap;							// 边编码-边交点索引的哈希表；

		unsigned curVersCount = 0;
		unsigned curTrisCount = 0;

		// 1. march over all cubes (loop order chosen to match memory)
		/*
			 Should be possible to parallelize safely if threads are "well separated".
			 Like red-black Gauss Seidel. Probably each thread need's their own edgeIsctMap, versResult, trisResult,
				   and then merge at the end.
			 Annoying part are the edges lying on the  interface between chunks.
		*/
		versResult.resize(std::pow(nx * ny * nz, 2. / 3.), 3);
		trisResult.resize(std::pow(nx * ny * nz, 2. / 3.), 3);
		for (int z = 0; z < nz - 1; z++)
		{
			for (int y = 0; y < ny - 1; y++)
			{
				for (int x = 0; x < nx - 1; x++)
				{
					// 1.1 计算当前栅格的索引：
					const unsigned gridIdx = getGridIdx(x, y, z);

					// 1.2 计算当前栅格对应的立方体的八个顶点的数据；
					static Eigen::Matrix<Scalar, 8, 1> cornerSDF;						// 方块的八个顶点的SDF值
					static Eigen::Matrix<unsigned, 8, 1> cornerIdx;               // 方块的八个顶点在栅格中的索引
					for (int i = 0; i < 8; i++)
					{
						const unsigned originIdx = gridIdx + cornerIdxOffset[i];
						cornerIdx(i) = originIdx;
						cornerSDF(i) = scalarFied(originIdx);
					}

					// 1.3 生成当前立方体内的三角片
					handleCube(gridCenters, cornerSDF, cornerIdx, isovalue, versResult, curVersCount, trisResult, curTrisCount, edgeIsctMap);
				}
			}
		}

		// 2. shrink_to_fit();
		versResult.conservativeResize(curVersCount, 3);
		trisResult.conservativeResize(curTrisCount, 3);


		// 3. 可能会有duplicated vertices:




		return true;
	}


	//		距离场数据的高斯滤波：
	bool gaussFilterSDFdata(Eigen::VectorXd& SDF, const int xCount, const int yCount, const int zCount)
	{
		Eigen::MatrixXd maskGauss(3, 3);
		maskGauss << 0.0113, 0.0838, 0.0113, 0.0838, 0.6193, 0.0838, 0.0113, 0.0838, 0.0113;
		int sliceSize = xCount * yCount;
		for (int i = 0; i < zCount; ++i)
		{
			Eigen::VectorXd slice = SDF.segment(sliceSize * i, sliceSize);
			Eigen::MatrixXd sliceMat = Eigen::Map<Eigen::MatrixXd>(slice.data(), xCount, yCount);
			Eigen::MatrixXd filteredMat;
			linearSpatialFilter(filteredMat, sliceMat, maskGauss);
			slice = Eigen::Map<Eigen::VectorXd>(filteredMat.data(), sliceSize, 1);
			SDF.segment(sliceSize * i, sliceSize) = slice;
		}
		return true;
	}


	// test4——测试生成栅格、marchingCubes算法：
	void test4()
	{
		tiktok& tt = tiktok::getInstance();
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		std::vector<int> stepCounts(3);				// xyz三个维度上栅格数
		Eigen::RowVector3d gridsOri;					// 栅格原点：
		Eigen::VectorXd SDF;
		Eigen::RowVector3i gridCounts;			// xyz三个维度上的栅格数量；
		Eigen::MatrixXd gridCenters;				//	 所有栅格中点坐标的矩阵，每行都是一个中点坐标；存储优先级是x, y, z
		Eigen::MatrixXd	boxVers;				// 栅格对应的包围盒的顶点；
		Eigen::MatrixXi boxTris;
		double selectedSDF = 0.75;						// 提取的水平集对应的距离场值；

		// 0. 解析SDFGen.exe生成的.sdf距离场数据文件：
		std::string fileName = "E:/材料/jawMesh6";
		std::string sdfFilePath = fileName + std::string{ ".sdf" };
		std::string objFilePath = fileName + std::string{ ".obj" };
		double SDFstep = IGL_BASIC::parseSDF(stepCounts, gridsOri, SDF, sdfFilePath.c_str());
		int xCount = stepCounts[0];
		int yCount = stepCounts[1];
		int zCount = stepCounts[2];
		igl::readOBJ(objFilePath, vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		// 0. 距离场高斯滤波	：
		tt.start();
		Eigen::MatrixXd maskGauss(3, 3);
		maskGauss << 0.0113, 0.0838, 0.0113, 0.0838, 0.6193, 0.0838, 0.0113, 0.0838, 0.0113;
		int sliceSize = xCount * yCount;
		/*PARALLEL_FOR(0, zCount, [&](int i)*/
		for(int i = 0; i<zCount; ++i)
			{
				Eigen::VectorXd slice = SDF.segment(sliceSize * i, sliceSize);
				Eigen::MatrixXd sliceMat = Eigen::Map<Eigen::MatrixXd>(slice.data(), xCount, yCount);
				Eigen::MatrixXd filteredMat;
				linearSpatialFilter(filteredMat, sliceMat, maskGauss);
				slice = Eigen::Map<Eigen::VectorXd>(filteredMat.data(), sliceSize, 1);
				SDF.segment(sliceSize * i, sliceSize) = slice;
			}
		tt.endCout("elapsed time of gaussian filtering is : ");

		// 1. 生成栅格：
		Eigen::RowVector3d minp = gridsOri;
		Eigen::RowVector3d maxp = gridsOri + SDFstep * Eigen::RowVector3d(stepCounts[0] - 1, stepCounts[1] - 1, stepCounts[2] - 1);
		Eigen::AlignedBox<double, 3> box(minp, maxp);		// 栅格对应的包围盒；
		genGrids(box, std::max({ stepCounts[0], stepCounts[1], stepCounts[2] }), 0, gridCenters, gridCounts);
		genAABBmesh(box, boxVers, boxTris);
		objWriteMeshMat("E:/AABB.obj", boxVers, boxTris);

		{
			//			栅格数据的分布：
			/*
				按索引增大排列的栅格中心点为：
				gc(000), gc(100), gc(200), gc(300),...... gc(010), gc(110), gc(210), gc(310),...... gc(001), gc(101), gc(201).....


				x坐标：
				x0, x1, x2, x3......x0, x1, x2, x3......x0, x1, x2, x3......
				周期为xCount;
				重复次数为(yCount * zCount)

				y坐标：
				y0, y0, y0...y1, y1, y1...y2, y2, y2.........y0, y0, y0...y1, y1, y1...
				周期为(xCount * yCount);
				重复次数为zCount;
				单个元素重复次数为xCount

				z坐标：
				z0, z0, z0......z1, z1, z1......z2, z2, z2......
				单个元素重复次数为(xCount * yCount)
			*/
			Eigen::MatrixXd gridCenters0;				// for try——尝试自己生成栅格数据：
			Eigen::RowVector3i gridCounts0{ stepCounts[0], stepCounts[1], stepCounts[2] };
			Eigen::VectorXd xPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(0), minp(0), maxp(0));
			Eigen::VectorXd yPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(1), minp(1), maxp(1));
			Eigen::VectorXd zPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(2), minp(2), maxp(2));

			Eigen::MatrixXd tmpVec0, tmpVec1, tmpVec2;
			kron(tmpVec0, Eigen::VectorXd::Ones(gridCounts(1) * gridCounts(2)), xPeriod);
			Eigen::VectorXd tmpVec11 = kron(yPeriod, Eigen::VectorXi::Ones(gridCounts(0)));
			kron(tmpVec1, Eigen::VectorXi::Ones(gridCounts(2)), tmpVec11);
			kron(tmpVec2, zPeriod, Eigen::VectorXd::Ones(gridCounts(0) * gridCounts(1)));
			gridCenters0.resize(stepCounts[0] * stepCounts[1] * stepCounts[2], 3);
			gridCenters0.col(0) = tmpVec0;
			gridCenters0.col(1) = tmpVec1;
			gridCenters0.col(2) = tmpVec2;

			//				提取栅格中SDF值小于0的顶点：
			Eigen::MatrixXd tmpVers(SDF.rows(), 3);
			int index = 0;
			for (int i = 0; i < SDF.rows(); ++i)
				if (SDF(i) <= 0)
					tmpVers.row(index++) = gridCenters0.row(i);
			tmpVers.conservativeResize(index, 3);
			igl::writeOBJ("E:/tmpVers.obj", tmpVers, Eigen::MatrixXi{});
		}

		// 2. marching cubes算法生成最终曲面：
		Eigen::MatrixXd versResult_SDF, versResults_signs, versResult_origin;
		Eigen::MatrixXi trisResult_SDF, trisResults_signs, trisResult_origin;
		tt.start();
		marchingCubes(versResult_SDF, trisResult_SDF, SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), selectedSDF);
		tt.endCout("Elapsed time of my marching_cubes() is ");
		igl::writeOBJ("E:/shrinkedMesh.obj", versResult_SDF, trisResult_SDF);

		// 原始的marching cubes
		versResult_SDF.resize(0, 0);
		trisResult_SDF.resize(0, 0);
		tt.start();
		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), selectedSDF, versResult_SDF, trisResult_SDF);
		tt.endCout("Elapsed time of igl::marching_cubes() is ");
		igl::writeOBJ("E:/shrinkedMeshOri.obj", versResult_SDF, trisResult_SDF);
 
		std::cout << "finished." << std::endl;
	}


	// test44
	void test44()
	{
		tiktok& tt = tiktok::getInstance();

		Eigen::VectorXd SDF;
		Eigen::RowVector3i gridCounts;			// xyz三个维度上的栅格数量；
		Eigen::MatrixXd gridCenters;				//	 所有栅格中点坐标的矩阵，每行都是一个中点坐标；存储优先级是x, y, z
		Eigen::MatrixXd	boxVers;						// 栅格对应的包围盒的顶点；
		Eigen::MatrixXi boxTris;
		double selectedSDF = 0.75;						// 提取的水平集对应的距离场值；

		std::vector<int> stepCounts{ 126, 132, 39 };														// xyz三个维度上栅格数
		Eigen::RowVector3d gridsOri{ -23.4242516, -40.1728897, -1.50000501 };					// 栅格原点：
		double SDFstep = 0.5;
		int xCount = stepCounts[0];
		int yCount = stepCounts[1];
		int zCount = stepCounts[2];
 
		// 0. 生成栅格：
		Eigen::RowVector3d minp = gridsOri;
		Eigen::RowVector3d maxp = gridsOri + SDFstep * Eigen::RowVector3d(stepCounts[0] - 1, stepCounts[1] - 1, stepCounts[2] - 1);
		Eigen::AlignedBox<double, 3> box(minp, maxp);		// 栅格对应的包围盒；
		genGrids(box, std::max({ stepCounts[0], stepCounts[1], stepCounts[2] }), 0, gridCenters, gridCounts);
		genAABBmesh(box, boxVers, boxTris);
		objWriteMeshMat("E:/AABB.obj", boxVers, boxTris);
		{
			//			栅格数据的分布：
			/*
				按索引增大排列的栅格中心点为：
				gc(000), gc(100), gc(200), gc(300),...... gc(010), gc(110), gc(210), gc(310),...... gc(001), gc(101), gc(201).....


				x坐标：
				x0, x1, x2, x3......x0, x1, x2, x3......x0, x1, x2, x3......
				周期为xCount;
				重复次数为(yCount * zCount)

				y坐标：
				y0, y0, y0...y1, y1, y1...y2, y2, y2.........y0, y0, y0...y1, y1, y1...
				周期为(xCount * yCount);
				重复次数为zCount;
				单个元素重复次数为xCount

				z坐标：
				z0, z0, z0......z1, z1, z1......z2, z2, z2......
				单个元素重复次数为(xCount * yCount)
			*/
			Eigen::MatrixXd gridCenters0;				// for try——尝试自己生成栅格数据：
			Eigen::RowVector3i gridCounts0{ stepCounts[0], stepCounts[1], stepCounts[2] };
			Eigen::VectorXd xPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(0), minp(0), maxp(0));
			Eigen::VectorXd yPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(1), minp(1), maxp(1));
			Eigen::VectorXd zPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(2), minp(2), maxp(2));

			Eigen::MatrixXd tmpVec0, tmpVec1, tmpVec2;
			kron(tmpVec0, Eigen::VectorXd::Ones(gridCounts(1) * gridCounts(2)), xPeriod);
			Eigen::VectorXd tmpVec11 = kron(yPeriod, Eigen::VectorXi::Ones(gridCounts(0)));
			kron(tmpVec1, Eigen::VectorXi::Ones(gridCounts(2)), tmpVec11);
			kron(tmpVec2, zPeriod, Eigen::VectorXd::Ones(gridCounts(0) * gridCounts(1)));
			gridCenters0.resize(stepCounts[0] * stepCounts[1] * stepCounts[2], 3);
			gridCenters0.col(0) = tmpVec0;
			gridCenters0.col(1) = tmpVec1;
			gridCenters0.col(2) = tmpVec2;

			//				提取栅格中SDF值小于0的顶点：
			Eigen::MatrixXd tmpVers(SDF.rows(), 3);
			int index = 0;
			for (int i = 0; i < SDF.rows(); ++i)
				if (SDF(i) <= 0)
					tmpVers.row(index++) = gridCenters0.row(i);
			tmpVers.conservativeResize(index, 3);
			igl::writeOBJ("E:/tmpVers.obj", tmpVers, Eigen::MatrixXi{});
		}


		// 1. 生成SDF数据：
		unsigned spVersCount = gridCenters.rows();
		SDF.resize(spVersCount);
		for (unsigned i = 0; i < spVersCount; ++i)
			SDF(i) = gridCenters(i, 2);

		gaussFilterSDFdata(SDF, xCount, yCount, zCount);

		// 2原始的marching cubes
		Eigen::MatrixXd versResult;
		Eigen::MatrixXi trisResult;
		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), 0, versResult, trisResult);
		igl::writeOBJ("E:/meshXOY.obj", versResult, trisResult);

		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), 1.7, versResult, trisResult);
		igl::writeOBJ("E:/meshSDF1.7.obj", versResult, trisResult);


		std::cout << "finished." << std::endl;
	}


	// 补洞：opological_hole_fill()——当前有问题
	void test5()
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris, bdrys;
		igl::readOBJ("E:/材料/holeMesh.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);		

		const unsigned versCount = vers.rows();

		// 计算洞的有序顶点索引
		std::vector<std::vector<int>> holes, bdrySegs;
		findHolesBdrySegs(holes, bdrySegs, vers, tris);

		// 需要手动为原网格添加洞的中心点，而且貌似处理多个洞的情形会出错；
		Eigen::MatrixXd holeVers;
		subFromIdxVec(holeVers, vers, holes[0]);
		Eigen::RowVector3d holeCenter = holeVers.colwise().mean();
		Eigen::MatrixXd versOut(versCount + 1, 3);
		versOut.topRows(versCount) = vers;
		versOut.bottomRows(1) = holeCenter;

		Eigen::MatrixXi trisOut;
		igl::topological_hole_fill(tris, Eigen::VectorXi{}, holes, trisOut);
		igl::writeOBJ("E:/meshFilled.obj", versOut, trisOut);
		std::cout << "finished." << std::endl;
	}


	// 缠绕数winding number:
	void test6() 
	{
		Eigen::MatrixXd vers0, vers1;
		Eigen::MatrixXi tris0, tris1;
		igl::readOBJ("E:/材料/tooth.obj", vers0, tris0);
		igl::readOBJ("E:/材料/cylinder1.obj", vers1, tris1);
		Eigen::VectorXd wNums;

		igl::winding_number(vers0, tris0, vers1, wNums);			// 网格外部则缠绕数为0；

		std::vector<double> wNumsVec = vec2Vec(wNums);
		std::vector<int> vec;
		vec.reserve(wNums.rows());
		const double eps = 1e-6;
		for (unsigned i = 0; i< wNums.rows(); ++i) 
		{
			if (std::abs(wNums(i)) < eps)
				vec.push_back(i);
		}
		vec.shrink_to_fit();

		Eigen::MatrixXd versOut;
		subFromIdxVec(versOut, vers1, vec);
		objWriteVerticesMat("E:/versOut.obj", versOut);

		std::cout << "finished." << std::endl;
	}


	// 测试网格单连通区域提取connected_components()
	void test7() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/材料/原型数据/originalMesh.obj");
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();


		// 1. 生成邻接矩阵：
		Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris);

		Eigen::SparseMatrix<int> adjSM = adjSM_eCount;
		traverseSparseMatrix(adjSM, [&](auto& iter)
			{
				iter.valueRef() = 1;
			});

		// 2. 确定单连通区域——connected_components()
		Eigen::VectorXi connectedLabels, connectedCount;
		int conCount = igl::connected_components(adjSM, connectedLabels, connectedCount);

		// 3. 提取最大的单连通区域中的顶点：
		std::vector<int> retVec1, retVec2;
		retVec1 = vec2Vec(connectedLabels);
		retVec2 = vec2Vec(connectedCount);

		int mainLabel = 0;						// 包含最大联通顶点数的标签；
		int mainLabelCount = 0;
		for (int i = 0; i < conCount; ++i)
		{
			if (connectedCount(i) > mainLabelCount)
			{
				mainLabel = i;
				mainLabelCount = connectedCount(i);
			}
		}

		std::unordered_set<int> indexSet;
		for (int i = 0; i < versCount; ++i)
		{
			if (mainLabel == connectedLabels(i))
				indexSet.insert(i);
		}

		std::vector<int> indexVec;
		indexVec.insert(indexVec.end(), indexSet.begin(), indexSet.end());
		Eigen::MatrixXd versOut;
		subFromIdxVec(versOut, vers, indexVec);

		std::vector<int> oldNewIdxInfo(versCount, -1);
		for (int i = 0; i < indexVec.size(); ++i)
		{
			int oldIdx = indexVec[i];
			oldNewIdxInfo[oldIdx] = i;
		}

		// 4. 提取最大单连通区域中的三角片：
		Eigen::MatrixXi trisCopy = tris;
		int* intPtr = trisCopy.data();
		for (int i = 0; i < trisCopy.size(); ++i)
		{
			int oldIdx = *intPtr;
			*intPtr = oldNewIdxInfo[oldIdx];
			intPtr++;
		}

		Eigen::MatrixXi trisOut = Eigen::MatrixXi::Zero(trisCount, 3);
		int trisCountNew = 0;
		for (int i = 0; i < trisCount; ++i)
		{
			if (trisCopy(i, 0) >= 0 && trisCopy(i, 1) >= 0 && trisCopy(i, 2) >= 0)
			{
				trisOut.row(trisCountNew) = trisCopy.row(i);
				trisCountNew++;
			}
		}

		trisOut.conservativeResize(trisCountNew, 3);
		objWriteMeshMat("E:/mainConnectedMesh.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// 自己实现的connected_conponents()
	void test77()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/材料/原型数据/originalMesh.obj");
		objWriteMeshMat("E:/meshInput.obj", vers, tris);

		bool retFlag = simplyConnectedLargest(versOut, trisOut, vers, tris);
		if (!retFlag)
			std::cout << "function failed!!!" << std::endl;

		objWriteMeshMat("E:/meshOut.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// 测量两个网格之间的hausdorff distance
	void test8() 
	{
		Eigen::MatrixXd vers0, vers1, vers2, vers3, vers4;
		Eigen::MatrixXi tris0, tris1, tris2, tris3, tris4;
		double hd0, hd1, hd2, hd3, hd4;
		hd0 = hd1 = hd2 = hd3 = hd4 = 0;

		tiktok& tt = tiktok::getInstance();
		tt.start();
		objReadMeshMat(vers0, tris0, "E:/材料/jawMeshDense.obj");
		objReadMeshMat(vers1, tris1, "E:/材料/jawMeshDenseSimplifyOutIsctCleaned.obj");
		objReadMeshMat(vers2, tris2, "E:/材料/jawMeshDenseQEMoutputIsctCleaned.obj");
		objReadMeshMat(vers3, tris3, "E:/材料/jawMeshDenseQslimOutputIsctCleaned.obj");
		objReadMeshMat(vers4, tris4, "E:/材料/jawMeshDensePMPsimplifiedIsctCleaned.obj");
		tt.endCout("meshes loading finished:");

		igl::hausdorff(vers0, tris0, vers0, tris0, hd0);
		igl::hausdorff(vers0, tris0, vers1, tris1, hd1);
		igl::hausdorff(vers0, tris0, vers2, tris2, hd2);
		igl::hausdorff(vers0, tris0, vers3, tris3, hd3);
		igl::hausdorff(vers0, tris0, vers4, tris4, hd4);
		debugDisp("Hausdorff distance0 == ", hd0);
		debugDisp("Hausdorff distance1 == ", hd1);
		debugDisp("Hausdorff distance2 == ", hd2);
		debugDisp("Hausdorff distance3 == ", hd3);
		debugDisp("Hausdorff distance4 == ", hd4);

		std::cout << "finished." << std::endl;
	}
 
}
