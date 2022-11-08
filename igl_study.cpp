#include "igl_study.h"
#define DATA_PATH "./data/"

igl::opengl::glfw::Viewer viewer;				// libigl中的基于glfw的显示窗口；
static std::mutex g_mutex;


// libigl基本功能
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

		// per_vertex_normals
		Eigen::MatrixXd verNorms, triNorms, conerNorms;			// 每行是一个法向量；
		igl::per_vertex_normals(vers, tris, verNorms);
		igl::per_face_normals(vers, tris, triNorms);
		igl::per_corner_normals(vers, tris, conerNorms);

		// unique_edge_map()——计算网格的半边、无向边、及其对应关系；
		Eigen::MatrixXi edges, uEdges;
		Eigen::VectorXi edgeUeInfo;
		std::vector<std::vector<int>> UeEdgeInfo;
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, UeEdgeInfo);
		//std::cout << "tris:" << std::endl;
		//dispMat(tris);
		//std::cout << "edges:" << std::endl;
		//dispMat(edges);				// 半边
		//std::cout << "uEdges:" << std::endl;
		//dispMat(uEdges);			// 无向边，既是undirected的又是unique的；
		//std::cout << "edgeUeInfo:" << std::endl;
		//dispMat(edgeUeInfo);		// 半边在无向边中对应的索引；edgeUeInfo(i)是索引为i的半边对应的无向边的索引；

		Eigen::VectorXi uEC, uEE;
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, uEC, uEE);
		//std::cout << "uEC:" << std::endl;
		//dispMat(uEC);		 
		//std::cout << "uEE:" << std::endl;
		//dispMat(uEE);

		// edge_flaps()——计算无向边关联的三角片；
		Eigen::MatrixXi UeTrisInfo, UeCornersInfo;
		igl::edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
		//std::cout << "UeTrisInfo:" << std::endl;
		//dispMat(UeTrisInfo);								// 第i行的元素——索引为i的无向边关联的三角片的索引；
		//std::cout << "UeCornersInfo:" << std::endl;
		//dispMat(UeCornersInfo);

		objWriteEdgesMat("E:/edge0.obj", Eigen::RowVector2i(uEdges.row(0)), vers);
		Eigen::RowVector3i tri00 = tris.row(UeTrisInfo(0, 0));
		Eigen::RowVector3i tri01 = tris.row(UeTrisInfo(0, 1));
		Eigen::MatrixXd oppVers(2, 3);
		oppVers.row(0) = vers.row(tri00(UeCornersInfo(0, 0)));
		oppVers.row(1) = vers.row(tri01(UeCornersInfo(0, 1)));
		objWriteVerticesMat("E:/oppVers.obj", oppVers);

		// circulation()——寻找输入边的端点关联的三角片，及其1领域顶点；
		std::vector<int> nbrVersIdx, nbrTrisIdx;
		igl:: circulation(0, true, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx, nbrTrisIdx);
		objWriteVerticesMat("E:/edgeHead0.obj", Eigen::MatrixXd{ vers.row(uEdges(0, 0)) });

		Eigen::MatrixXd vers0;
		subFromIdxVec(vers0, vers, nbrVersIdx);
		objWriteVerticesMat("E:/vers0.obj", vers0);
		Eigen::MatrixXi tris0;
		subFromIdxVec(tris0, tris, nbrTrisIdx);
		igl::writeOBJ("E:/tris0.obj", vers, tris0);

		// igl::doubleArea()——计算输入网格每个三角片面积的两倍：
		Eigen::VectorXd dbArea;
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


	// decimate()，qslim() 直接实现网格精简：
	void test7() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
		Eigen::VectorXi newOldVersInfo;						 
		tiktok& tt = tiktok::getInstance();

		std::string fileName = "tmpMesh";

		tt.start();
		igl::readOBJ((std::string{ "E:/材料/" } + fileName + std::string{".obj"}).c_str(), vers, tris);
		tt.endCout("elapsed time of loading mesh is: ");

		unsigned trisCount = tris.rows();
		unsigned tarTrisCount = std::round(trisCount * 0.5);			// 输出网格的最大三角片数；
		// unsigned tarTrisCount = 7516;
		igl::writeOBJ("E:/meshIn.obj", vers, tris);

#if 1
		// igl::decimate()——边折叠算法精简网格； 网格不可以有非流形边；				
		tt.start();
		std::cout << "succeeded? " << igl::decimate(vers, tris, tarTrisCount, versOut, trisOut, \
			newOldTrisInfo, newOldVersInfo) << std::endl;						// 重载1.1
		tt.endCout("Elapsed time of mesh simplification by igl::decimate() is ");
		igl::writeOBJ((std::string{ "E:/" } + fileName + std::string{ "_simplified_edgeCollapse.obj" }).c_str(), versOut, trisOut);

		std::vector<int> newOldTrisInfoVec = vec2Vec(newOldTrisInfo);
		Eigen::MatrixXd vers1;
		subFromIdxVec(vers1, vers, newOldVersInfo);
#endif

		// igl::qslim()——貌似比igl::decimate()慢一些，目前不知道优势在哪里；
		tt.start();
		std::cout << "qslim succeeded? " << igl::qslim(vers, tris, tarTrisCount, versOut, trisOut, newOldTrisInfo, newOldVersInfo) << std::endl;
		tt.endCout("Elapsed time of mesh simplification by igl::qslim() is ");
		igl::writeOBJ((std::string{ "E:/" } + fileName + std::string{ "_simplified_qslim.obj" }).c_str(), versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// collapse_edge()——网格精简中的边折叠算法：
	void test77() 
	{
		// 折叠一条边：
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		igl::readOBJ("E:/材料/cylinder1.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		// 生成边数据：
		Eigen::MatrixXi edges, uEdges, UeTrisInfo, UeCornersInfo;
		Eigen::VectorXi edgeUeInfo ;
		Eigen::VectorXi uEC, uEE;
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, uEC, uEE);
		igl::edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);

		int uEdgeIdx0 = 50;
		Eigen::RowVector2i uEdge0 = uEdges.row(uEdgeIdx0);
		Eigen::RowVector3d ver00 = vers.row(uEdge0(0));
		Eigen::RowVector3d ver01 = vers.row(uEdge0(1));
		Eigen::RowVector3d collapsedVer = (ver00 + ver01) / 2.0;
		objWriteEdgesMat("E:/uEdge0.obj", uEdge0, vers);
		
		// 1. igl::collapse_edge(); 折叠索引为uEdgeIdx0的边；
		int e1, e2, f1, f2;
		igl::collapse_edge(uEdgeIdx0, collapsedVer, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, e1, e2, f1, f2);
		igl::writeOBJ("E:/meshOut_删除三角片前.obj", vers, tris);

		// 2. 删除所有含有标记为IGL_COLLAPSE_EDGE_NULL边的三角片：
		MatrixXi tris0(tris.rows(), 3);
		int m = 0;
		for (int i = 0; i < tris.rows(); i++)
		{
			if (tris(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
				tris(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
				tris(i, 2) != IGL_COLLAPSE_EDGE_NULL)
			{
				tris0.row(m) = tris.row(i);
				m++;
			}
		}
		tris0.conservativeResize(m, tris0.cols());					// 这里相当于shrink_to_fit();
		igl::writeOBJ("E:/meshOut_beforeRemove.obj", vers, tris0);

		// 3. igl::remove_unreferenced()——删除网格中的孤立顶点；
		Eigen::VectorXi I, newOldVersInfo;
		igl::remove_unreferenced(vers, tris0, versOut, trisOut, I, newOldVersInfo);
		igl::writeOBJ("E:/meshOut.obj", versOut, trisOut);


		std::cout << "finished." << std::endl;
	}


	// qslim中逐个折叠边：
	void test777()
	{
		using namespace igl;

		using Quadric = std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>;

		Eigen::MatrixXd vers, versOut, vers0, versOri;
		Eigen::MatrixXi tris, trisOut, tris0, trisOri;
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
		Eigen::VectorXi newOldVersInfo;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
		Eigen::VectorXi timeStamps;
		MatrixXd collapsedVers;
		Eigen::VectorXd costs;
		std::vector<Quadric> quadrics;								// 每个顶点的Q矩阵；
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// 优先队列；
		decimate_cost_and_placement_callback cost_and_placement;
		decimate_pre_collapse_callback pre_collapse;
		decimate_post_collapse_callback post_collapse;
		decimate_stopping_condition_callback stop_condition;
		Eigen::VectorXi _1, I2;
		bool clean_finish = false;
		int v1 = -1;				// State variables keeping track of edge we just collapsed
		int v2 = -1;

		// 0.
		tiktok& tt = tiktok::getInstance();
		tt.start();
		igl::readOBJ("E:/材料/tmpMesh.obj", versOri, trisOri);
		igl::writeOBJ("E:/meshIn.obj", versOri, trisOri);
		tt.endCout("elapsed time of loading mesh is: ");

		// 1. 
		unsigned trisCount = trisOri.rows();
		igl::connect_boundary_to_infinity(versOri, trisOri, vers, tris);
		if (!igl::is_edge_manifold(tris))
			return;
		igl::edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);


		// 2. 计算每个顶点的Q矩阵：
		igl::per_vertex_point_to_plane_quadrics(vers, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, quadrics);


		// 3.  获取qslim算法中的函数子：
		int trisCountNew = trisCount;
		int tarTrisCount = std::round(0.5 * trisCount);
		igl::qslim_optimal_collapse_edge_callbacks(uEdges, quadrics, v1, v2, cost_and_placement, pre_collapse, post_collapse);
		stop_condition = max_faces_stopping_condition(trisCountNew, trisCount, tarTrisCount);


		// 4. 计算每条无向边的cost值，以此为优先级存入优先队列
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());
		collapsedVers.resize(uEdges.rows(), vers.cols());
		costs.resize(uEdges.rows());
		igl::parallel_for(uEdges.rows(), [&](const int e)
			{
				double cost = e;
				RowVectorXd p(1, 3);
				cost_and_placement(e, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, cost, p);
				collapsedVers.row(e) = p;
				costs(e) = cost;
			},
			10000);

		for (int i = 0; i < uEdges.rows(); i++)
			pQueue.emplace(costs(i), i, 0);


		// 5. 逐步执行边折叠；
		for (unsigned i =0; i <= 127; ++i) 
		{
			Eigen::MatrixXd versBC = vers;
			Eigen::MatrixXi trisBC = tris;
			Eigen::MatrixXi uEdgesBC = uEdges;
			auto tuple0 = pQueue.top();
			int uEdgeIdx0 = std::get<1>(tuple0);			// 本轮折叠的边的索引；

			// 5.1
			int e, e1, e2, f1, f2;
			if (!collapse_edge(cost_and_placement, pre_collapse, post_collapse, \
				vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, \
				pQueue, timeStamps, collapsedVers, e, e1, e2, f1, f2))          // collapse_edge()重载2；pre_collapse()操作，折叠边，执行post_collapse()，再更新相关的边的cost值；
				assert("edge collapse failed.");

			// 5.2. 删除所有含有标记为IGL_COLLAPSE_EDGE_NULL边的三角片：
			tris0.resize(tris.rows(), 3);
			newOldTrisInfo.resize(tris.rows());
			int index = 0;
			for (int i = 0; i < tris.rows(); i++)
			{
				if (tris(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
					tris(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
					tris(i, 2) != IGL_COLLAPSE_EDGE_NULL)
				{
					tris0.row(index) = tris.row(i);
					newOldTrisInfo(index) = i;
					index++;
				}
			}
			tris0.conservativeResize(index, tris0.cols());              // 这里相当于shrink_to_fit();
			newOldTrisInfo.conservativeResize(index);

			// 5.3. 删除网格中的孤立顶点：
			igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);

			//		打印当前步骤的结果：
			if (i == 126 || i == 127) 
			{
				char str[256];
				sprintf_s(str, 256, "E:/qslim折叠单条边loop%d.obj", i);
				igl::writeOBJ(str, versOut, trisOut);

				sprintf_s(str, 256, "E:/loop%d折叠的边.obj", i);
				Eigen::MatrixXi uEdge0 = uEdgesBC.row(uEdgeIdx0);
				objWriteEdgesMat(str, uEdge0, versBC);

				// for debug:
				Eigen::VectorXi uEdgeVec0 = uEdge0.col(0);
				Eigen::VectorXi uEdgeVec1 = uEdge0.col(1);
				std::vector<int> tmpVec0 = vec2Vec(uEdgeVec0);
				std::vector<int> tmpVec1 = vec2Vec(uEdgeVec0);

				std::cout << "loop " << i << " finished." << "折叠的边：" << uEdge0(0, 0) << ", " << uEdge0(0, 1) << std::endl;
			}
		}
 
		std::cout << "finished." << std::endl;
	}

 
	// 打印优先队列队首的若干元素：
	void dispCosts(const igl::min_heap<std::tuple<double, int ,int>>& pQueue, const unsigned num) 
	{
		igl::min_heap<std::tuple<double, int, int>> pCopy = pQueue;
		std::vector<double> costValues;
		std::vector<std::vector<int>> edges;
		costValues.reserve(num);
		edges.reserve(num);

		for (unsigned i = 0; i<num; ++i) 
		{
			auto& t = pCopy.top();
			costValues.push_back(std::get<0>(t));
			edges.push_back(std::vector<int>{std::get<1>(t), std::get<2>(t)});
			pCopy.pop();
		}

		for (unsigned i = 0; i < num; ++i)
			std::cout << costValues[i] << ", ";

		std::cout << std::endl << std::endl;
	}


	// 使用viewer对象的动画功能逐步查看精简过程中的边折叠：
	void test7777() 
	{
		using namespace igl;
		Eigen::MatrixXd vers;
		Eigen:MatrixXi tris;
		igl::readOBJ("E:/材料/tmpMesh.obj", vers, tris);

		bool keepWorkingFlag = true;
		bool simplest_decimate_flag = false;
		bool qslim_decimate_flag = ~simplest_decimate_flag;

		// 1. 准备数据
		Eigen::MatrixXd versCopy = vers;
		Eigen::MatrixXi trisCopy = tris;
		VectorXi edgeUeInfo;
		MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// 优先队列； 

		int v1 = -1;					// State variables keeping track of edge we just collapsed
		int v2 = -1;
		Eigen::VectorXi timeStamps;
		MatrixXd collapsedVers;
		typedef std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double> Quadric;
		std::vector<Quadric> quadrics;				// 每个顶点的Q矩阵；

		edge_flaps(trisCopy, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());						// Could reserve with https://stackoverflow.com/a/29236236/148668        
		collapsedVers.resize(uEdges.rows(), versCopy.cols());				// If an edge were collapsed, we'd collapse it to these points:

		// 检测是否有非流形边：
		{
			Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> BF;
			Eigen::Array<bool, Eigen::Dynamic, 1> BE;
			if (!is_edge_manifold(trisCopy, uEdges.rows(), edgeUeInfo, BF, BE))
				return;
		}
 
		decimate_cost_and_placement_callback cost_and_placement;
		decimate_pre_collapse_callback pre_collapse;
		decimate_post_collapse_callback post_collapse;

		if (simplest_decimate_flag)
		{
			// 获取最简边折叠精简算法的函数子
			cost_and_placement = shortest_edge_and_midpoint;
			decimate_trivial_callbacks(pre_collapse, post_collapse);    // 生成最简边折叠中的pre_collapse和post_collapse函数子——什么都不做；
		}
		else
		{
			//  计算每个顶点的Q矩阵：
			per_vertex_point_to_plane_quadrics(vers, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, quadrics);

			// 获取qslim算法中的函数子：
			igl::qslim_optimal_collapse_edge_callbacks(uEdges, quadrics, v1, v2, cost_and_placement, pre_collapse, post_collapse);
		}

		decimate_stopping_condition_callback stopping_condition;
 
		int loopCount = 0;
		int num_collapsed;                                // 边收缩循环的计数；

		const auto& exportCurrentMesh = [&]()
		{
			MatrixXi tris0(trisCopy.rows(), 3);
			VectorXi _1;
			Eigen::VectorXi newOldTrisInfo;
			Eigen::VectorXi newOldVersInfo;
			newOldTrisInfo.resize(trisCopy.rows());
			int m = 0;
			for (int i = 0; i < trisCopy.rows(); i++)
			{
				if (trisCopy(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
					trisCopy(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
					trisCopy(i, 2) != IGL_COLLAPSE_EDGE_NULL)
				{
					tris0.row(m) = trisCopy.row(i);
					newOldTrisInfo(m) = i;
					m++;
				}
			}
			tris0.conservativeResize(m, tris0.cols());              // 这里相当于shrink_to_fit();
			newOldTrisInfo.conservativeResize(m);

			// 3.2. 删除网格中的孤立顶点：
			Eigen::MatrixXd versOut;
			Eigen::MatrixXi trisOut;
			igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);
			char str[256];
			sprintf_s(str, 256, "E:/decimating_output_%d.obj", loopCount);
			igl::writeOBJ(str, versOut, trisOut);
		};
 
		// lambda——所有数据复位
		const auto& reset = [&]()
		{
			tris = trisCopy;
			vers = versCopy;
			edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
			collapsedVers.resize(uEdges.rows(), vers.cols());
			VectorXd costs(uEdges.rows());

			pQueue = {};            // https://stackoverflow.com/questions/2852140/priority-queue-clear-method
			timeStamps = Eigen::VectorXi::Zero(uEdges.rows());

			// r1. 计算每条边的折叠cost值，以及折叠之后的顶点坐标，存入优先队列pQueue；
			{
				Eigen::VectorXd costs(uEdges.rows());
				igl::parallel_for(uEdges.rows(), \
					[&](const int i)
					{
						double cost = i;
						RowVectorXd edgeCenter(1, 3);           // 取边的中点作为边折叠之后的顶点；
						cost_and_placement(i, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, cost, edgeCenter);
						collapsedVers.row(i) = edgeCenter;
						costs(i) = cost;
					}, \
					10000);

				for (int i = 0; i < uEdges.rows(); i++)
					pQueue.emplace(costs(i), i, 0);
			}

			num_collapsed = 0;
			viewer.data().clear();
			viewer.data().set_mesh(vers, tris);
			viewer.data().set_face_based(true);
		};


		// lambda——执行网格精简，准备渲染的数据；
		const auto& pre_draw = [&](igl::opengl::glfw::Viewer& viewer)->bool
		{
			// p1. 每一次动画循环中，收缩1%的边；
			if (viewer.core().is_animating && !pQueue.empty())
			{
				bool FlagCollapsed = false;           // 本次循环中边收缩是否执行成功；

				// p1.1 执行边收缩——collapse edge
				const int max_iter = std::ceil(0.01 * pQueue.size());
				for (int j = 0; j < max_iter; j++)
				{
					if (!collapse_edge(cost_and_placement, vers, tris, uEdges, \
						edgeUeInfo, UeTrisInfo, UeCornersInfo, pQueue, timeStamps, collapsedVers))              // collapse_edge()重载2.1
						break;
					FlagCollapsed = true;
					num_collapsed++;
				}

				// p1.2 
				if (FlagCollapsed)
				{
					viewer.data().clear();
					viewer.data().set_mesh(vers, tris);
					viewer.data().set_face_based(true);
				}
			}

			// p2. 计算折叠后的网格三角片数（三个点索引都没有被标记为IGL_COLLAPSE_EDGE_NULL的三角片数）
			unsigned currentTrisCount = 0;
			for (int i = 0; i < tris.rows(); i++)
			{
				if (tris(i, 0) != IGL_COLLAPSE_EDGE_NULL || tris(i, 1) != IGL_COLLAPSE_EDGE_NULL || tris(i, 2) != IGL_COLLAPSE_EDGE_NULL)
					currentTrisCount++;
			}

			if (loopCount < 10)
				exportCurrentMesh();
			std::cout << "loop : " << loopCount++ << ", current trisCount == " << currentTrisCount << std::endl;
			return false;
		};


		// lambda——键盘事件响应；
		const auto& key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)->bool
		{
			switch (key)
			{
			case ' ':
				viewer.core().is_animating ^= 1;
				break;
			case 'R':
			case 'r':
				reset();
				break;
			case '0':
				keepWorkingFlag = false;
				viewer.core().is_animating = 0;
				viewer.launch_shut();
				viewer.shutdown_plugins();
				break;
			default:
				return false;
			}
			return true;
		};


		// 1. 初始化（复位）；
		reset();

		// 2. 打开窗口，动画循环
		viewer.core().background_color.setConstant(1);
		viewer.core().is_animating = true;                                 // 动画
		viewer.core().animation_max_fps = 1.0;						// 指定最大动画帧率；
		viewer.callback_key_down = key_down;
		viewer.callback_pre_draw = pre_draw;
		viewer.launch(keepWorkingFlag);

		// 3. 动画终止后，输出网格：
		Eigen::VectorXi newOldTrisInfo;
		Eigen::VectorXi newOldVersInfo;
		
		// 3.1. 删除所有含有标记为IGL_COLLAPSE_EDGE_NULL边的三角片：
		MatrixXi tris0(trisCopy.rows(), 3);
		VectorXi _1;
		newOldTrisInfo.resize(trisCopy.rows());
		int m = 0;
		for (int i = 0; i < trisCopy.rows(); i++)
		{
			if (trisCopy(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
				trisCopy(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
				trisCopy(i, 2) != IGL_COLLAPSE_EDGE_NULL)
			{
				tris0.row(m) = trisCopy.row(i);
				newOldTrisInfo(m) = i;
				m++;
			}
		}
		tris0.conservativeResize(m, tris0.cols());              // 这里相当于shrink_to_fit();
		newOldTrisInfo.conservativeResize(m);

		// 3.2. 删除网格中的孤立顶点：
		Eigen::MatrixXd versOut;
		Eigen::MatrixXi trisOut;
		igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);
		igl::writeOBJ("E:/decimating_output.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// 解析qslim()
	void test77777() 
	{
		using namespace igl;

		using Quadric = std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>;

		Eigen::MatrixXd vers, versOut, vers0, versOri;
		Eigen::MatrixXi tris, trisOut, tris0, trisOri;
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]是精简后的网格中第i个三角片对应的原网格的三角片索引；
		Eigen::VectorXi newOldVersInfo;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
		Eigen::VectorXi timeStamps;
		MatrixXd collapsedVers;
		Eigen::VectorXd costs;
		std::vector<Quadric> quadrics;								// 每个顶点的Q矩阵；
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// 优先队列；
		decimate_cost_and_placement_callback cost_and_placement;
		decimate_pre_collapse_callback pre_collapse;
		decimate_post_collapse_callback post_collapse;
		decimate_stopping_condition_callback stop_condition;
		Eigen::VectorXi _1, I2;
		int prev_e = -1;
		bool clean_finish = false;
		int v1 = -1;				// State variables keeping track of edge we just collapsed
		int v2 = -1;

		// 0.
		tiktok& tt = tiktok::getInstance();
		tt.start();
		igl::readOBJ("E:/材料/tmpMesh.obj", versOri, trisOri);
		tt.endCout("elapsed time of loading mesh is: ");

		// 1. 
		unsigned trisCount = trisOri.rows();
		igl::connect_boundary_to_infinity(versOri, trisOri, vers, tris);
		if (!igl::is_edge_manifold(tris))
			return;
		igl::edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);


		// 2. 计算每个顶点的Q矩阵：
		igl::per_vertex_point_to_plane_quadrics(vers, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, quadrics);


		// 3.  获取qslim算法中的函数子：
		int trisCountNew = trisCount;
		int tarTrisCount = std::round(0.5 * trisCount);
		igl::qslim_optimal_collapse_edge_callbacks(uEdges, quadrics, v1, v2, cost_and_placement, pre_collapse, post_collapse);
		stop_condition = max_faces_stopping_condition(trisCountNew, trisCount, tarTrisCount);


		// 4. 计算每条无向边的cost值，以此为优先级存入优先队列
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());
		collapsedVers.resize(uEdges.rows(), vers.cols());
		costs.resize(uEdges.rows());
		igl::parallel_for(uEdges.rows(), [&](const int e)
			{
				double cost = e;
				RowVectorXd p(1, 3);
				cost_and_placement(e, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, cost, p);
				collapsedVers.row(e) = p;
				costs(e) = cost;
			},
			10000);

		for (int i = 0; i < uEdges.rows(); i++)
			pQueue.emplace(costs(i), i, 0);

		// 5. 边折叠的循环：
		while (true)
		{
			int e, e1, e2, f1, f2;

			// 执行pre_collapse()操作，折叠边，执行post_collapse()，再更新相关的边的cost值；
			if (collapse_edge(cost_and_placement, pre_collapse, post_collapse, \
				vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, \
				pQueue, timeStamps, collapsedVers, e, e1, e2, f1, f2))          // collapse_edge()重载2；
			{
				// 若满足终止条件——stopping_condition函数子返回true，则跳出折叠循环
				if (stop_condition(vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, pQueue, timeStamps, collapsedVers, e, e1, e2, f1, f2))
				{
					// LIBIGL中当前有两种终止条件函数子：max_faces_stopping_condition(), infinite_cost_stopping_condition();
					clean_finish = true;
					break;
				}
			}
			else  // 边折叠失败，退出循环：
			{
				if (e == -1)
					break;                // a candidate edge was not even found in pQueue.
				if (prev_e == e)
				{
					assert(false && "Edge collapse no progress... bad stopping condition?");
					break;
				}

				// Edge was not collapsed... must have been invalid. collapse_edge should have updated its cost to inf... continue
			}

			prev_e = e;
		}


		// 6. 删除所有含有标记为IGL_COLLAPSE_EDGE_NULL边的三角片：
		tris0.resize(tris.rows(), 3);
		newOldTrisInfo.resize(tris.rows());
		int index = 0;
		for (int i = 0; i < tris.rows(); i++)
		{
			if (tris(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
				tris(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
				tris(i, 2) != IGL_COLLAPSE_EDGE_NULL)
			{
				tris0.row(index) = tris.row(i);
				newOldTrisInfo(index) = i;
				index++;
			}
		}
		tris0.conservativeResize(index, tris0.cols());              // 这里相当于shrink_to_fit();
		newOldTrisInfo.conservativeResize(index);

		// 7. 删除网格中的孤立顶点：
		igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);

		//		for debug:
		igl::writeOBJ("E:/删除内部三角片前.obj", versOut, trisOut);

		// 8. ？？？删除内部三角片
		const Eigen::Array<bool, Eigen::Dynamic, 1> keep = (newOldTrisInfo.array() < trisCount);
		igl::slice_mask(Eigen::MatrixXi(trisOut), keep, 1, trisOut);
		igl::slice_mask(Eigen::VectorXi(newOldTrisInfo), keep, 1, newOldTrisInfo);
		igl::remove_unreferenced(Eigen::MatrixXd(versOut), Eigen::MatrixXi(trisOut), versOut, trisOut, _1, I2);
		igl::slice(Eigen::VectorXi(newOldVersInfo), I2, 1, newOldVersInfo);


		igl::writeOBJ("E:/qslimOutput.obj", versOut, trisOut);
		std::cout << "finished." << std::endl;
	}


	// 求精简后的网格和原始网格的近似误差（来自QEM的paper）：
	double calcSimpApproxError(const Eigen::MatrixXd& versSimp, const Eigen::MatrixXi& trisSimp, \
		const Eigen::MatrixXd& versOri, const Eigen::MatrixXi& trisOri)
	{
		// 注意不能出现太大的矩阵，否则容易bad alloc

		double appErr = 0;
		std::pair<unsigned, float> pair0;
		int versCount0 = versOri.rows();
		int versCount1 = versSimp.rows();
		const Eigen::MatrixXd& versMat0 = versOri;
		const Eigen::MatrixXd& versMat1 = versSimp;
		const Eigen::MatrixXi& trisMat0 = trisOri;
		const Eigen::MatrixXi& trisMat1 = trisSimp;
		Eigen::MatrixXd planeCoeff0, planeCoeff1;
		bool ret0 = trianglesPlane(planeCoeff0, versMat0, trisMat0);
		bool ret1 = trianglesPlane(planeCoeff1, versMat1, trisMat1);

		std::vector<long double> squaredDis0, squaredDis1;
		squaredDis0.reserve(versCount0);
		squaredDis1.reserve(versCount1);

		// 点到平面距离： dis == p.dot(v)，注意是有符号的距离； p == (a,b,c,d), v = (x0, y0, z0, 1), (a,b,c)是平面的归一化法向量；
		Eigen::MatrixXd versExt0{ Eigen::MatrixXd::Ones(versCount0, 4) };
		Eigen::MatrixXd versExt1{ Eigen::MatrixXd::Ones(versCount1, 4) };
		versExt0.leftCols(3) = versMat0.array().cast<double>();
		versExt1.leftCols(3) = versMat1.array().cast<double>();

		// for debug
		Eigen::Vector4d verColVec = versExt0.row(3).transpose();
		// 索引为verIdx的simp网格中的顶点，到ori网格中所有三角片平面的符号距离：
		Eigen::VectorXd signedDises = planeCoeff0 * verColVec;
		Eigen::VectorXd sqrDises = signedDises.array() * signedDises.array();

		// for new method:
		PARALLEL_FOR(0, versCount1, [&](int verIdx)
			{
				std::lock_guard<mutex> guard(g_mutex);
				Eigen::Vector4d verColVec = versExt1.row(verIdx).transpose();
				// 索引为verIdx的simp网格中的顶点，到ori网格中所有三角片平面的符号距离：
				Eigen::VectorXd signedDises = planeCoeff0 * verColVec;
				Eigen::VectorXd sqrDises = signedDises.array() * signedDises.array();
				squaredDis1.push_back(sqrDises.minCoeff());
			});

		PARALLEL_FOR(0, versCount0, [&](int verIdx)
			{
				std::lock_guard<mutex> guard(g_mutex);
				Eigen::Vector4d verColVec = versExt0.row(verIdx).transpose();
				// 索引为verIdx的simp网格中的顶点，到ori网格中所有三角片平面的符号距离：
				Eigen::VectorXd signedDises = planeCoeff1 * verColVec;
				Eigen::VectorXd sqrDises = signedDises.array() * signedDises.array();
				squaredDis0.push_back(sqrDises.minCoeff());
			});

		long double s0 = 0;
		long double s1 = 0;
		for (const auto& num : squaredDis0)
			s0 += num;
		for (const auto& num : squaredDis1)
			s1 += num;
		appErr = (s0 + s1) / (static_cast<long>(versCount0) + static_cast<long>(versCount1));

		return appErr;
	}


	// 测量精简网格的近似误差：
	void test777777() 
	{
		Eigen::MatrixXd vers0, vers1;
		Eigen::MatrixXi tris0, tris1;
		objReadMeshMat(vers0, tris0, "E:/材料/jawMeshDense4.obj");
		objReadMeshMat(vers1, tris1, "E:/jawMeshSimplified4.obj");

		double appErr = calcSimpApproxError(vers1, tris1, vers0, tris0);

		std::cout << "appErr == " << appErr << std::endl;
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
		igl::adjacency_matrix(tris, adjSM);				// 有向边存在则元素为1，否则为0；

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


		// 3. 确定非流形边；
#if __cplusplus < 202002L
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

		//	 nonManifoldEdges()——输出非流形半边：
		Eigen::MatrixXi  nmnEdges;
		nonManifoldEdges(tris, nmnEdges);
		dispMatBlock(nmnEdges, 0, 10, 0, 1);
		objWriteEdgesMat("E:/nmnEdges.obj", nmnEdges, vers);


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


	// test buildAdjacency();
	void test2() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris, edges;
		igl::readOBJ("E:/材料/jawMeshArranged.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();

		Eigen::MatrixXi ttAdj_nmEdge;
		std::vector<ttTuple> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		bool retFlag = buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);

		// 输出包含非流形边的三角片：
		std::vector<int> nmnTrisIdx;
		for (const auto& tuple : ttAdj_nmnEdge)
		{
			nmnTrisIdx.insert(nmnTrisIdx.end(), std::get<0>(tuple).begin(), std::get<0>(tuple).end());
			nmnTrisIdx.insert(nmnTrisIdx.end(), std::get<1>(tuple).begin(), std::get<1>(tuple).end());
			nmnTrisIdx.insert(nmnTrisIdx.end(), std::get<2>(tuple).begin(), std::get<2>(tuple).end());
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
		nonManifoldEdges(tris, nmnEdges);
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
		Eigen::MatrixXi ttAdj_nmEdge;
		std::vector<ttTuple> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);
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
				// wf0. 若当前边为流形边，则在ttAdj_nmEdge中寻找其对面；
				int nbrTriIdx = ttAdj_nmEdge(currentTriIdx, i);
				
				// wf1. 若当前边为非流形边，根据布尔操作不同类型选取不同的对面；
				if (-1 == nbrTriIdx)
				{
					std::vector<std::vector<int>*> vecPtr(3), tmpPtrVec(3);
					vecPtr[0] = &std::get<0>(ttAdj_nmnEdge[currentTriIdx]);
					vecPtr[1] = &std::get<1>(ttAdj_nmnEdge[currentTriIdx]);
					vecPtr[2] = &std::get<2>(ttAdj_nmnEdge[currentTriIdx]);
					tmpPtrVec[0] = &std::get<0>(ttAdj_nmnOppEdge[currentTriIdx]);
					tmpPtrVec[1] = &std::get<1>(ttAdj_nmnOppEdge[currentTriIdx]);
					tmpPtrVec[2] = &std::get<2>(ttAdj_nmnOppEdge[currentTriIdx]);

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


	void test33()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		igl::readOBJ("E:/材料/jawMeshArranged.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);
		igl::writeOBJ("E:/triNMN0.obj", vers, Eigen::MatrixXi{ tris.row(4431) });

		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();

		Eigen::MatrixXi ttAdj_nmEdge;
		std::vector<ttTuple> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		bool retFlag = buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);

		// 输出包含非流形边的三角片：
		std::vector<int> nmnTrisIdx;
		for (const auto& tuple : ttAdj_nmnEdge)
		{
			nmnTrisIdx.insert(nmnTrisIdx.end(), std::get<0>(tuple).begin(), std::get<0>(tuple).end());
			nmnTrisIdx.insert(nmnTrisIdx.end(), std::get<1>(tuple).begin(), std::get<1>(tuple).end());
			nmnTrisIdx.insert(nmnTrisIdx.end(), std::get<2>(tuple).begin(), std::get<2>(tuple).end());
		}
		std::unique(nmnTrisIdx.begin(), nmnTrisIdx.end());
		Eigen::MatrixXi nmnTris;
		subFromIdxVec(nmnTris, tris, nmnTrisIdx);
		objWriteMeshMat("E:/nmnTris.obj", vers, nmnTris);

		Eigen::MatrixXi mnTris;				// 不包含非流形边的三角片：
		std::set<int> tmpSet;
		traverseMatrix(ttAdj_nmEdge, [&](const int num)
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


	// test4——测试生成栅格、marchingCubes算法：
	void test4()
	{
		std::vector<int> stepCounts(3);				// xyz三个维度上栅格数
		Eigen::RowVector3d gridsOri;					// 栅格原点：
		Eigen::VectorXd SDF;
		Eigen::RowVector3i gridCounts;			// xyz三个维度上的栅格数量；
		Eigen::MatrixXd gridCenters;				//	 所有栅格中点坐标的矩阵，每行都是一个中点坐标；存储优先级是x, y, z
		Eigen::MatrixXd	boxVers;				// 栅格对应的包围盒的顶点；
		Eigen::MatrixXi boxTris;

		// 0. 解析SDFGen.exe生成的.sdf距离场数据文件：
		const char* sdfFilePath = "E:/tooth.sdf";
		double SDFstep = IGL_BASIC::parseSDF(stepCounts, gridsOri, SDF, sdfFilePath);

		// 1. 生成栅格：
		Eigen::RowVector3d minp = gridsOri;
		Eigen::RowVector3d maxp = gridsOri + SDFstep * Eigen::RowVector3d(stepCounts[0], stepCounts[1], stepCounts[2]);
		Eigen::AlignedBox<double, 3> box(minp, maxp);		// 栅格对应的包围盒；
		genGrids(box, std::max({ stepCounts[0], stepCounts[1], stepCounts[2] }), 0, gridCenters, gridCounts);
		genAABBmesh(box, boxVers, boxTris);
		objWriteMeshMat("E:/AABB.obj", boxVers, boxTris);

		// 栅格数据的分布：
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
		Eigen::RowVector3i gridCounts0{stepCounts[0], stepCounts[1], stepCounts[2]};
		Eigen::VectorXd xPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(0), minp(0), maxp(0));
		Eigen::VectorXd yPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(1), minp(1), maxp(1));
		Eigen::VectorXd zPeriod = Eigen::VectorXd::LinSpaced(gridCounts0(2), minp(2), maxp(2));

		Eigen::MatrixXd tmpVec0, tmpVec1, tmpVec2;
		kron(tmpVec0, VectorXd::Ones(gridCounts(1) * gridCounts(2)), xPeriod);
		Eigen::VectorXd tmpVec11 = kron(yPeriod, VectorXi::Ones(gridCounts(0)));
		kron(tmpVec1, VectorXi::Ones(gridCounts(2)), tmpVec11);
		kron(tmpVec2, zPeriod, VectorXd::Ones(gridCounts(0) * gridCounts(1)));
		gridCenters0.resize(stepCounts[0]* stepCounts[1] * stepCounts[2], 3);
		gridCenters0.col(0) = tmpVec0;
		gridCenters0.col(1) = tmpVec1;
		gridCenters0.col(2) = tmpVec2;

		// 提取栅格中SDF值小于0的顶点：
		Eigen::MatrixXd tmpVers(SDF.rows(), 3);
		int index = 0;
		for (int i = 0; i<SDF.rows(); ++i) 
			if (SDF(i) <= 0)
				tmpVers.row(index++) = gridCenters0.row(i);
		tmpVers.conservativeResize(index, 3);
		igl::writeOBJ("E:/tmpVers.obj", tmpVers, Eigen::MatrixXi{});

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


	// 补洞：opological_hole_fill()
	void test5()
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("E:/材料/holeMesh.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		// 确定洞的边——只关联一个三角片的边：
		const unsigned trisCount = tris.rows();
		const unsigned edgesCount = 3 * trisCount;
		const unsigned versCount = tris.maxCoeff() + 1;

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

		// 构建邻接矩阵：
		std::vector<Eigen::Triplet<int>> smElems;
		smElems.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});
		Eigen::SparseMatrix<int> adjSM_eCount;
		adjSM_eCount.resize(versCount, versCount);
		adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());		// 有向边邻接矩阵，权重为该有向边重复的次数；
		Eigen::SparseMatrix<int> tmpSp = adjSM_eCount.transpose();
		Eigen::SparseMatrix<int> adjSM_ueCount = adjSM_eCount + tmpSp;		// 无向边邻接矩阵，权重为该边重复的次数；

		// 寻找边缘无向边：
		std::unordered_set<std::pair<int, int>, edgeHash> ueSet;
		for (unsigned i = 0; i < adjSM_ueCount.outerSize(); ++i)
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_ueCount, i); iter; ++iter)
				if (1 == iter.value())
					ueSet.insert(std::make_pair(iter.row(), iter.col()));

		// 寻找边缘半边：
		std::unordered_set<std::pair<int, int>, edgeHash> bdryEdgeSet;
		for (const auto& pair : ueSet)
		{
			int vaIdx = pair.first;
			int vbIdx = pair.second;

			if (adjSM_eCount.coeff(vaIdx, vbIdx) > 0)
				bdryEdgeSet.insert(std::make_pair(vaIdx, vbIdx));
			if (adjSM_eCount.coeff(vbIdx, vaIdx) > 0)
				bdryEdgeSet.insert(std::make_pair(vbIdx, vaIdx));
		}
		Eigen::MatrixXi bdryEdges(bdryEdgeSet.size(), 2);
		unsigned index = 0;
		for (const auto& pair : bdryEdgeSet)
		{
			bdryEdges(index, 0) = pair.first;
			bdryEdges(index, 1) = pair.second;
			index++;
		}
		objWriteEdgesMat("E:/bdryEdges.obj", bdryEdges, vers);

		// 计算洞的有序顶点索引
		std::vector<std::vector<int>> hole(1);
		std::unordered_map<int, int> tmpMap;
		for (unsigned i = 0; i < bdryEdges.rows(); ++i)
			tmpMap.insert(std::make_pair(bdryEdges(i, 0), bdryEdges(i, 1)));

		unsigned headIdx = tmpMap.begin()->first;
		unsigned tailIdx = tmpMap.begin()->second;
		tmpMap.erase(tmpMap.begin());
		while (!tmpMap.empty())
		{
			hole[0].push_back(headIdx);
			headIdx = tailIdx;
			tailIdx = tmpMap[headIdx];
			tmpMap.erase(headIdx);
		}
		hole[0].push_back(headIdx);

		// 需要手动为原网格添加洞的中心点，而且貌似处理多个洞的情形会出错；
		Eigen::MatrixXd holeVers;
		subFromIdxVec(holeVers, vers, hole[0]);
		Eigen::RowVector3d holeCenter = holeVers.colwise().mean();
		Eigen::MatrixXd versOut(versCount + 1, 3);
		versOut.topRows(versCount) = vers;
		versOut.bottomRows(1) = holeCenter;

		Eigen::MatrixXi trisOut;
		igl::topological_hole_fill(tris, Eigen::VectorXi{}, hole, trisOut);
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
}
