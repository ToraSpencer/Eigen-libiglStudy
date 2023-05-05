#include "igl_study.h"
#define DATA_PATH "./data/"

static igl::opengl::glfw::Viewer viewer;				// libigl�еĻ���glfw����ʾ���ڣ�
static std::mutex g_mutex;


/// /////////////////////////////////////////////////////////////////////////////////////////// DEBUG �ӿ�

static std::string g_debugPath = "E:/";

static void debugDisp()			// �ݹ���ֹ
{						//		�ݹ���ֹ��Ϊ�޲λ���һ�����������ζ����ԡ�
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



////////////////////////////////////////////////////////////////////////////// libigl��������
namespace IGL_BASIC
{
	Eigen::MatrixXd vers, newVers, normals;
	Eigen::MatrixXi tris;
	Eigen::SparseMatrix<double> L;


	// igl�л����ľ������ӣ�
	void test00() 
	{
		Eigen::MatrixXi m1(3, 3);
		Eigen::MatrixXi m11;
		Eigen::VectorXi vec1(2);

		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		dispMat(m1);

		// slice()����ʹ��������������������ȡ�����е�Ԫ�أ�
		/*
			  template <typename MatX,  typename DerivedR,  typename MatY>
			  IGL_INLINE void slice(
						const MatX& X,																	�������
						const Eigen::DenseBase<DerivedR> & R,							��������
						const int dim,															ά�ȡ���1-�У�2-�У�
						MatY& Y																				�������
					)
		*/

		vec1 << 0, 2;
		igl::slice(m1, vec1, 1, m11);				// ��ȡ��0, 2�У�
		dispMat(m11);

		igl::slice(m1, vec1, 2, m11);				// ��ȡ��0, 2�У�
		dispMat(m11);

		Eigen::VectorXi vec2(2);
		vec2 << 1, 2;
		igl::slice(m1, vec1, vec2, m11);		// �ȼ��ڣ�m11 == m1(0: 2, 1: 2);
		dispMat(m11);


		std::cout << "finished." << std::endl;
	}


	// igl�м�������������ԵĽӿڣ�
	void test000()
	{
		bool retFlag = igl::readOBJ("E:/����/tooth.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		Eigen::MatrixXd verNorms, triNorms, conerNorms;			// ���ַ���ÿ����һ����������
		Eigen::MatrixXi edges;									// ���
		Eigen::MatrixXi uEdges;								// ����ߣ�����undirected������unique�ģ�
		Eigen::VectorXi edgeUeInfo;						// �����������ж�Ӧ��������edgeUeInfo(i)������Ϊi�İ�߶�Ӧ������ߵ�������
		Eigen::MatrixXi UeTrisInfo;							// �����-��������Ƭӳ�������i�е�Ԫ�أ�����Ϊi������߹���������Ƭ��������
		Eigen::MatrixXi UeCornersInfo;					// �����-�Զ���ӳ�������i�е�Ԫ�أ�����Ϊi����������Ե������������������Ϊ�Ǳ�Ե���α����������Զ��㣻
		std::vector<std::vector<int>> UeEdgeInfo;
		Eigen::VectorXi uEC, uEE;
		Eigen::VectorXd dbArea;								// ������ÿ������Ƭ�����������
		std::vector<int> nbrVersIdx, nbrTrisIdx;

		// per_vertex_normals
		igl::per_vertex_normals(vers, tris, verNorms);
		igl::per_face_normals(vers, tris, triNorms);
		igl::per_corner_normals(vers, tris, conerNorms);

		// unique_edge_map()������������İ�ߡ�����ߡ������Ӧ��ϵ��
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, UeEdgeInfo);
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, uEC, uEE);

		// edge_flaps()������������߹���������Ƭ��
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


		// circulation()����Ѱ������ߵĶ˵����������Ƭ������1���򶥵㣻
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


		// igl::doubleArea()����������������ÿ������Ƭ�����������
		igl::doublearea(vers, tris, dbArea);
		std::vector<double> dbAreaVec = vec2Vec(dbArea);
		std::sort(dbAreaVec.begin(), dbAreaVec.end());
		std::cout << "finished." << std::endl;
	}


	// �ļ�IO
	void test0()
	{
		//// readOBJ(), writeObj()����OBJ�ļ���IO���ж�����أ�������·�������ơ�����Ƭ���������ݡ�
		igl::readOBJ("./data/bunny.obj", vers, tris);
		igl::writeOBJ("./data/bunny_export.obj", vers, tris);
		igl::writeOBJ("./data/bunnyVers.obj", vers, Eigen::MatrixXi{});			// ֻҪ���Ʋ���Ҫ����Ƭ�Ļ�������վ���
	
		vers.resize(0, 0);
		tris.resize(0, 0);
		igl::readOBJ("E:/����/jawMeshDense.obj", vers, tris);
		igl::writeOFF("E:/����/jawMeshDense.off", vers, tris);

		// ��ȡstl�ļ���
		std::string fileName{"E:/����/jawMeshSimplified"};
		vers.resize(0, 0);
		tris.resize(0, 0);
		std::ifstream fileIn((fileName + std::string{ ".stl" }).c_str(), std::ios::binary);			// stl�ļ��Ƕ������ļ���
		igl::readSTL(fileIn, vers, tris, normals);
		fileIn.close();
		igl::writeOBJ((fileName + std::string{".obj"}).c_str(), vers, tris);

		// ��.mesh�ļ�
		fileName = "E:/����/bunny";
		vers.resize(0, 0);
		tris.resize(0, 0);
		Eigen::MatrixXi tets;
		igl::readMESH((fileName + std::string{ ".mesh" }).c_str(), vers, tets, tris);
		igl::writeOBJ((fileName + std::string{ ".obj" }).c_str(), vers, tris);
		igl::writeOBJ((fileName + std::string{ "_tets.obj" }).c_str(), vers, tets);
		std::cout << "finished." << std::endl;
	}


	// libigl�е���ʾ������Viewer
	void test1() 
	{
		igl::readOBJ("bunny.obj", vers, tris);

		// Viewer::data()��������viewer�����ݶ�������ã�

		// ViewerData::set_mesh()�������붥�������Ƭ���������������ݣ�д�뵽ViewerData�ĳ�Ա�����У�
		viewer.data().set_mesh(vers, tris);		// ������װ�����ݣ�


		 // Ĭ����������ת��set_rotation_type()��������ָ��������ת���
		viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);				// ������ת

		viewer.data().show_lines = 0;			// show_linesָ���Ƿ񻭳������ߣ�

		viewer.launch();
	}

 
	// �����ݶ�
	void test2()
	{
		using namespace Eigen;
		using namespace std;

		igl::readOBJ("./data/rootTooth1.obj", vers, tris);

		Eigen::VectorXd U;
		igl::readDMAT("./data/cheburashka-scalar.dmat", U);		// һϵ�еĺ���ֵ

		Eigen::SparseMatrix<double> G;			// �ݶ�����
		igl::grad(vers, tris, G);

		// ���㺯��ֵU���ݶ�
		Eigen::MatrixXd GU = Map<const Eigen::MatrixXd>((G * U).eval().data(), tris.rows(), 3);

		// ����ֵU���ݶ�GU��ģ��
		const Eigen::VectorXd GU_mag = GU.rowwise().norm();

		viewer.data().set_mesh(vers, tris);
		viewer.data().set_data(U);

		// Average edge length divided by average gradient (for scaling)
		const double max_size = igl::avg_edge_length(vers, tris) / GU_mag.mean();

		// ÿ������Ƭ�����ϻ�һ��ָʾ�ߣ�����Ϊ�ݶȷ��� 
		Eigen::MatrixXd BC;
		igl::barycenter(vers, tris, BC);
		const Eigen::RowVector3d black(0, 0, 0);
		viewer.data().add_edges(BC, BC + max_size * GU, black);
		viewer.data().show_lines = false;	  // ����������

		viewer.launch();
	}
 

	// ������ά�ռ��е�դ��
	void test4() 
	{
		Eigen::MatrixXd vers, gridCenters;
		Eigen::MatrixXi tris;	
		Eigen::RowVector3i gridCounts;			// ����ά����դ�����Ŀ��
		unsigned num = 10;				// ��������Ǹ�ά��(xyz�е�һ��)��դ������

		igl::readOBJ("E:/����/tooth.obj", vers, tris);
		igl::voxel_grid(vers, 0, num, 1, gridCenters, gridCounts);

		igl::writeOBJ("E:/gridCenters.obj", gridCenters, Eigen::MatrixXi{});

		std::cout << "finished." << std::endl;
	}


	// test5. ������ž��볡��marching cubes��ȡ��ֵ������
	void test5()
	{
		Eigen::MatrixXi tris;
		Eigen::MatrixXd vers;
		igl::readOBJ("E:/����/jawMesh.obj", vers, tris);

		tiktok& tt = tiktok::getInstance();
		double gridStep = 0.5;			
		double range[3];
		for (unsigned i = 0; i < 3; ++i)
		{
			Eigen::VectorXd coors = vers.col(i);
			range[i] = coors.maxCoeff() - coors.minCoeff();
		}

		// 1. ����դ��
		double largestRange = std::max({ range[0], range[1], range[2] });
		int num = std::ceil((largestRange/gridStep));              // ��������Ǹ�ά��(xyz�е�һ��)��դ������
		double gridsOffset = (gridStep * num - largestRange) / 2.;
		
		Eigen::MatrixXd gridCenters;
		Eigen::RowVector3i gridCounts;
		tt.start();
		igl::voxel_grid(vers, gridsOffset, num, 1, gridCenters, gridCounts);
		tt.endCout("Elapsed time of igl::voxel_grid() is ");


		// ������ע��libigl�еļ���SDF�ӿ�igl::signed_distance���ɵľ��볡ò��û��SDFgen�еĺ��ã���ʱ������
		Eigen::VectorXd SDF, signValues;
		{
			Eigen::VectorXi I;                     // useless
			Eigen::MatrixXd C, N;              // useless

			// 2. ������ž��볡
			tt.start();
			igl::signed_distance(gridCenters, vers, tris, igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, SDF, I, C, N);
			tt.endCout("Elapsed time of igl::signed_distance() is ");

			// 3. ���ž��볡��дΪ���ų�����������Ϊ-1����������Ϊ0������Ϊ1��
			signValues = SDF;
			std::for_each(signValues.data(), signValues.data() + signValues.size(), [](double& b)\
			{
				b = (b > 0 ? 1 : (b < 0 ? -1 : 0));
			});
		}

		// 4. marching cubes�㷨�����������棺
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


	// ��ȡSDFGen.exe���ɵ�.sdf���볡���ݣ�ʹ��igl::marching_cubes()��ȡ��ֵ������

	//����.sdf�ı��ļ���
	double parseSDF(std::vector<int>& stepCounts, Eigen::RowVector3d& gridsOri, Eigen::VectorXd& SDF, const char* filePath)
	{
		/*
			double parseSDF(												���ؾ��볡�Ĳ���������
					std::vector<int>& stepCounts,					xyz���������ϵĲ���
					Eigen::RowVector3d& gridsOri,					���볡�ռ��ԭ��
					Eigen::VectorXd& SDF,								���볡����
					const char* filePath										SDF�ļ�Ŀ¼
					)
		
		*/
		double SDFstep = -1;
		stepCounts.resize(3);
		std::string readStr(1024, '\0');
		std::ifstream sdfFile(filePath);
		if (!sdfFile)
			return SDFstep;

		// ��һ�У�����
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

		// �ڶ��У�դ��ԭ��
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

		// �����У����볡�ռ䲽����
		sdfFile.getline(&readStr[0], 1024);
		SDFstep = std::stod(readStr);

		// ������֮�󣺾��볡���ݣ�
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


	// test55. ��ȡSDFgen���ɵķ��ž��볡���ݣ�ʹ��marching cubes��������
	void test55() 
	{
		// 0. ����SDFGen.exe���ɵ�.sdf���볡�����ļ���
		std::vector<int> stepCounts(3);				// xyz����ά����դ����
		Eigen::RowVector3d gridsOri;					// դ��ԭ�㣺
		Eigen::VectorXd SDF;
		const char* sdfFilePath = "E:/jawMeshUnionRepair1.sdf";
		double SDFstep = parseSDF(stepCounts, gridsOri, SDF, sdfFilePath);

		// 1. ����դ��
		Eigen::RowVector3i gridCounts;
		Eigen::MatrixXd gridCenters;
		Eigen::RowVector3d minp = gridsOri - SDFstep * Eigen::RowVector3d(stepCounts[0] / 2.0, stepCounts[1] / 2.0, stepCounts[2] / 2.0);
		Eigen::RowVector3d maxp = gridsOri + SDFstep * Eigen::RowVector3d(stepCounts[0] / 2.0, stepCounts[1] / 2.0, stepCounts[2] / 2.0);
		Eigen::AlignedBox<double, 3> box(minp, maxp);
		igl::voxel_grid(box, std::max({stepCounts[0], stepCounts[1], stepCounts[2]}), 0, gridCenters, gridCounts);
		
		// 2. marching cubes�㷨�����������棺
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


	// �����Χ���࣬�����Χ����
	void test6() 
	{
		// ���ɰ�Χ������
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
 
		// contains()�����ж϶������Χ�еĹ�ϵ��������Ҫ����������ʾ���ڰ�Χ�б�����Ҳ��Ϊ���ڲ���
		std::cout << "aabb.contains (3,4,5) ? " << aabb.contains(Eigen::RowVector3d(3, 4, 5).transpose()) << std::endl;
		std::cout << "aabb.contains (4,5,5) ? " << aabb.contains(Eigen::RowVector3d(4, 5, 5).transpose()) << std::endl;
		std::cout << "aabb.contains (9,9,9) ? " << aabb.contains(Eigen::RowVector3d(9, 9, 9).transpose()) << std::endl;

		// OBB����д��contains()������
		std::cout << "obb.contains(0, 0, 0) ? " << obb.contains(Eigen::RowVector3d(0, 0, 0)) << std::endl;
		std::cout << "obb.contains(4, 2, 6) ? " << obb.contains(Eigen::RowVector3d(4, 2, 6)) << std::endl;
		std::cout << "obb.contains(-5, -5, -5) ? " << obb.contains(Eigen::RowVector3d(-5, -5, -5)) << std::endl;
 
		igl::writeOBJ("E:/426.obj", Eigen::MatrixXd{ Eigen::RowVector3d(4, 2, 6) }, Eigen::MatrixXi{});

		std::cout << "finished." << std::endl;
	}

}


////////////////////////////////////////////////////////////////////////////// libigl�е�΢�ּ������
namespace IGL_DIF_GEO 
{
	Eigen::MatrixXd vers_g, newVers_g, normals_g;
	Eigen::MatrixXi tris_g;
	Eigen::SparseMatrix<double> L_g;
	float deltaLB_g;
	unsigned smoothLoopCount_g = 0;


	// test0: ���������LB����
	void test0() 
	{
		Eigen::MatrixXd vers_g;
		Eigen::MatrixXi tris_g;
		Eigen::SparseMatrix<double> L_g, M;

		igl::readOBJ("E:/����/tooth.obj", vers_g, tris_g);
		igl::cotmatrix(vers_g, tris_g, L_g);
		igl::massmatrix(vers_g, tris_g, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);

		dispSpMat(M, 0, M.rows() - 1, 10);

		std::cout << "finished." << std::endl;
	}
 

	// test1: ʹ��Laplacian��˳����
	
	// laplace��˳�����¼���ʹ��laplacian��˳����
	bool key_down_laplacian(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)
	{
		std::string outputName;
		outputName.resize(512);
		switch (key)
		{
		case 'r':

		case 'R':							// ��λ����
			newVers_g = vers_g;
			smoothLoopCount_g = 0;
			break;

		case 's':

		case 'S':			
			sprintf_s(&outputName[0], 512, "E:/meshSmoothLB_step%d.obj", smoothLoopCount_g);
			igl::writeOBJ(outputName, newVers_g, tris_g);
			std::cout << "step " << smoothLoopCount_g << " result saved." << std::endl;
			break;

		case ' ':								// �ո����ִ��һ��laplace��˳
		{
			// 1. ���¼�����������
			Eigen::SparseMatrix<double> mass;
			igl::massmatrix(newVers_g, tris_g, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

			// 2. �����Է����� (mass - delta*L_g) * newVers_g = mass * newVers_g
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

		// 1.a ֱ�ӹ���laplacian����Compute Laplace-Beltrami operator: 
		igl::cotmatrix(vers_g, tris_g, L_g);

		// 1.b ���Էֲ�����laplacian
		{
			Eigen::SparseMatrix<double> Gradient, L2;
			Eigen::VectorXd dblA;									  // ÿ������Ƭ���������
			igl::grad(vers_g, tris_g, Gradient);      // ��ɢ�ݶ�

			// Diagonal per-triangle "mass matrix"			
			igl::doublearea(vers_g, tris_g, dblA);           

			// Place areas along diagonal  #dim times
			const auto& T = 1. * (dblA.replicate(3, 1) * 0.5).asDiagonal();

			L2 = -Gradient.transpose() * T * Gradient;         // discrete Dirichelet energy Hessian ��ɢ��������������������
			std::cout << "���ַ����õ���laplacian�Ĳ�ķ�����" << std::endl;
			std::cout << "(L2 - L_g).norm() == " << (L2 - L_g).norm() << std::endl;
		}

		// 2. ����ԭʼ�ķ�������ʹ��αɫ
		igl::per_vertex_normals(vers_g, tris_g, norms);
		colors = norms.rowwise().normalized().array() * 0.5 + 0.5;

		// 3. viewr����ʼ����
		newVers_g = vers_g;
		viewer.data().set_mesh(newVers_g, tris_g);
		viewer.data().set_colors(colors);
		viewer.callback_key_down = key_down_laplacian;

		// 4. ����
		std::cout << "Press [space] to smooth." << std::endl;
		std::cout << "Press [r] to reset." << std::endl;

		viewer.launch();
	}

}


/////////////////////////////////////////////////////////////////////////////// ͼ�㷨
namespace IGL_GRAPH 
{
	// ͼ���ݽṹ��ת����
	void test0() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("E:/����/cube�ظ�����Ƭ.obj", vers, tris);

		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		
		// ����Ƭ��Ϣ�õ�ͼ���ڽӱ��ڽӾ���
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);				// ����ߴ�����Ԫ��Ϊ1������Ϊ0��

		// �ֶ������ڽӾ���
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		unsigned edgesCount = 3 * trisCount;
		Eigen::SparseMatrix<int> adjSM_eCount;				// ������ڽӾ���Ȩ��Ϊ������ߵ�������

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

		std::cout << "����Ϊ0�Ķ����1���򶥵㣺" << std::endl;
		traverseSTL(adjList[0], disp<int>);
		std::cout << std::endl << std::endl;

		std::cout << "����Ϊ1�Ķ������ߣ����ڽӾ����е�1�������в�Ϊ���Ԫ���±꣩��" << std::endl;
		for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM, 0); it; ++it)
		{
			std::cout << "value == " << it.value() << std::endl;
			std::cout << "row == " << it.row() << std::endl;			 // row index
			std::cout << "col == " << it.col() << std::endl;			 // col index (here it is equal to k)
			std::cout << std::endl;
		} 

		std::cout << "����Ϊ2�Ķ���ĳ��ߣ����ڽӾ����е�2�������в�Ϊ���Ԫ���±꣩��" << std::endl;
		for (int i = 0; i < adjSM_eCount.outerSize(); ++i)			// ���еı�����
		{
			for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM_eCount, i); it; ++it)		// ��ǰ�е����ڵ�������
			{
				if (2 == it.row())
				{
					std::cout << "value == " << it.value() << std::endl;
					std::cout << "row == " << it.row() << std::endl;			 // row index
					std::cout << "col == " << it.col() << std::endl;			 // col index (here it is equal to k)
				}
			}
		}

		// �ж��Ƿ����ظ�����Ƭ���������ظ�����Ƭ��������ظ�����ߣ�
		bool hasDupTri = false;;
		for (int i = 0; i < adjSM_eCount.outerSize(); ++i)			// ���еı�����
		{
			if (hasDupTri)
				break;
			for (Eigen::SparseMatrix<int>::InnerIterator it(adjSM_eCount, i); it; ++it)		// ��ǰ�е����ڵ�������
			{
				if (it.value() > 1)
				{
					std::cout << "�����ظ�����Ƭ" << std::endl;
					hasDupTri = true;
					break;
				}
			}
		}
		if (!hasDupTri)
			std::cout << "û���ظ�����Ƭ" << std::endl;

		std::cout << "finished." << std::endl;
	}


	// ������ͼ�㷨��
	void test1() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("E:/����/roundSurf.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);
 
		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);

		size_t startIdx = 165;			// �ӽ�Բ�ĵĶ��㣻
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


		// �Լ�д��DFS�㷨����Ŀǰ������
		std::vector<bool> visited(vers.rows(), false);
		std::list<int> discoveredVersIdx;						// �ѱ����ʵĶ��㣺
		std::list<int> closedVersIdx;

		std::function<void(const int, const int)> myDfs = [&](const int index, const int parentIdx)
		{
			// �ݹ���ֹ��
			if (visited[index])			
				return;

			visited[index] = true;
			discoveredVersIdx.push_back(index);
			const std::vector<int>& adjVersIdx = adjList[index];

			// �ݹ���ƣ�
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
		igl::readOBJ("E:/����/roundSurf.obj", vers, tris);

		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);

		Eigen::VectorXd min_distance;				// ͼ�����ж��㵽ָ����������·�����ȣ�
		Eigen::VectorXi mst;						// ��С��������
		int verIdx0 = 165;						// �����Ϊ�ӽ�Բ�ĵĶ��㣻

		// ��ָ���ص�ʱ�� igl::dijkstra()����ͼ����verIdx0Ϊ������С��������
		int retIdx = igl::dijkstra(vers, adjList, verIdx0, std::set<int>{}, min_distance, mst);
		std::cout << "retIdx == " << retIdx << std::endl;
		objWriteTreePath("E:/mst.obj", mst, vers);

		std::cout << "finished." << std::endl;
	}
 
}


/////////////////////////////////////////////////////////////////////////////// �ռ仮��
namespace IGL_SPACE_PARTITION
{
	// aabb��ʵ�ֵ�BVH(��ΰ�Χ��)
	void test0()
	{
		Eigen::MatrixXd vers, vers0, minDisVers;
		Eigen::MatrixXi tris;
		Eigen::VectorXd minSqrDis;
		Eigen::VectorXi minDisIdx;
		igl::readOBJ("E:/����/tooth.obj", vers, tris);

		vers0.resize(2, 3);
		vers0.row(0) = vers.row(0);
		vers0.row(1) = vers.row(99);

		// igl::point_mesh_squared_distance()����ʹ��BVH�󶥵㵽������С���룻
		igl::point_mesh_squared_distance(vers0, vers, tris, minSqrDis, minDisIdx, minDisVers);
		igl::writeOBJ("E:/inputMesh.obj", vers, tris);
		igl::writeOBJ("E:/vers0.obj", vers0, Eigen::MatrixXi{});
		igl::writeOBJ("E:/minDisVers.obj", minDisVers, Eigen::MatrixXi{});

		//// ���������BVH����
		//igl::AABB<double, 3> bvh;

		

		std::cout << "finished." << std::endl;
	}

}


/////////////////////////////////////////////////////////////////////////////// IGLʵ�ֵĻ��������������㷨��
namespace IGL_BASIC_PMP 
{
	// �����������еı�-����Ƭ�ڽӹ�ϵ��
	void test1() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris, edges;
		igl::readOBJ("E:/����/meshArranged.obj", vers, tris);
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		unsigned edgesCount = 3 * trisCount;		// ������������edgesCount����ʵ���������Ҫ����Ϊ�����ظ���ͬһ����߽������ظ�������

		// edges == [ea; eb; ec] == [vbIdxes, vcIdxes; vcIdxes, vaIdxes; vaIdxes, vbIdxes];
		/*
			������������������edges�б���ÿ���߶���unique�ģ�
			�����ڷ����αߣ�������α���edges����ظ��洢��ʵ�ʱ�����Ҳ������edges������
			����˵ʹ������representation�Ļ��������α��в�ֹһ��������ȡ�������÷����αߵ�����Ƭ������
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

		// 1. �ڽӾ���
		std::vector<Eigen::Triplet<int>> smElems;
		smElems.reserve(edgesCount);
		for (unsigned i = 0; i < edgesCount; ++i)
			smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});

		Eigen::SparseMatrix<int> adjSM_eCount;
		adjSM_eCount.resize(versCount, versCount);
		adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());		// Ȩ��Ϊ��������ظ��Ĵ�����

		// 2. ������-����Ƭ����ӳ���etInfo(i)������Ϊi�ı����ڵ�����Ƭ��������
		std::vector<int> etInfo(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			etInfo[i] = i % trisCount;


		// 3. ȷ�������αߣ�
		std::unordered_set<std::int64_t> edgesNMN;					// ���ñ߱����ʾ
		 for (unsigned i = 0; i < adjSM_eCount.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_eCount, i); iter; ++iter)	// ��i�е����ڵ�������
			{
				if (iter.value() > 1)
				{
					int vaIdx = iter.row();
					int vbIdx = iter.col();
					edgesNMN.insert(encodeEdge(vaIdx, vbIdx));
				}
			}
		}

		 //		�����������-�������ļ�ֵ�ԣ�һ�������α߶�Ӧ�Ŷ����������
		 std::unordered_map<std::int64_t, std::vector<int>> edgesNMMmap;

		 //		�����������-�����ڵ�����Ƭ�������ļ�ֵ�ԣ�һ�������α߶�Ӧ�Ŷ������Ƭ������
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
					if (!retPair.second)			// ������ʧ�ܣ���˵�����д˼���
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
		igl::readOBJ("E:/����/jawMeshArranged.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();

		Eigen::MatrixXi ttAdj_mnEdge;
		std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		bool retFlag = buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris);

		// ������������αߵ�����Ƭ��
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


	// boolean����select tris:
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

		// 1. ��������Ƭ�ڽӹ�ϵ
		Eigen::MatrixXi ttAdj_mnEdge;
		std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris);
		unsigned trisCount = tris.rows();
		std::vector<bool> visited(trisCount, false);			// ����Ƭ�Ƿ񱻷��ʵı�ǣ�

		// 2. ʹ��������������ʽ��������Ƭ��
		std::deque<int> triIdxDeq;						// �������У�
		triIdxDeq.push_back(startIdx);
		visited[startIdx] = true;
		const int startLabel = trisLabel[startIdx];
		while (!triIdxDeq.empty())					// ����������ѭ����
		{
			int currentTriIdx = triIdxDeq.front();
			int currentLabel = trisLabel[currentTriIdx];
			triIdxDeq.pop_front();				
			
			// �Ե�ǰ����Ƭ�����ߵı�����
			for (int i = 0; i < 3; ++i)
			{
				// wf0. ����ǰ��Ϊ���αߣ�����ttAdj_mnEdge��Ѱ������棻
				int nbrTriIdx = ttAdj_mnEdge(currentTriIdx, i);
				
				// wf1. ����ǰ��Ϊ�����αߣ����ݲ���������ͬ����ѡȡ��ͬ�Ķ��棻
				if (-1 == nbrTriIdx)
				{
					std::vector<std::vector<int>*> vecPtr(3), tmpPtrVec(3);
					vecPtr[0] = &ttAdj_nmnEdge[currentTriIdx][0];
					vecPtr[1] = &ttAdj_nmnEdge[currentTriIdx][1];
					vecPtr[2] = &ttAdj_nmnEdge[currentTriIdx][2];
					tmpPtrVec[0] = &ttAdj_nmnOppEdge[currentTriIdx][0];
					tmpPtrVec[1] = &ttAdj_nmnOppEdge[currentTriIdx][1];
					tmpPtrVec[2] = &ttAdj_nmnOppEdge[currentTriIdx][2];

					const std::vector<int>& relaTrisIdx = *vecPtr[i];				// ��ǰ�����α����ڵ���������Ƭ�������Ļ�Ӧ��Ϊ��������
					const std::vector<int>& relaOppTrisIdx = *tmpPtrVec[i];	// ��ǰ�����αߵ����ж�������Ƭ�������Ļ�Ӧ��Ϊ��������
					switch (type)
					{
					case BOOLEAN_TYPE::DIFF:
						{
							assert(2 == relaTrisIdx.size() && "Exceptional non-manifold edge detected!");
							assert(2 == relaOppTrisIdx.size() && "Exceptional non-manifold edge detected!");

							// ��nbrTriIdx��������Ƭ��������һ����startIdx����Ƭ��ǩ��ͬ������ǩ���������Ƭ����ɢ��
							if (trisLabel[relaTrisIdx[0]] == startLabel || trisLabel[relaTrisIdx[1]] == startLabel)
							{
								nbrTriIdx = (trisLabel[relaTrisIdx[0]] == currentLabel ? relaTrisIdx[1] : relaTrisIdx[0]);
								assert(trisLabel[relaTrisIdx[0]] != trisLabel[relaTrisIdx[1]] && "Exceptional non-manifold edge detected!");
							}
							else    // ��nbrTriIdx��������Ƭ����startIdx����Ƭ��ǩ��ͬ������ǩ����Ķ�������Ƭ����ɢ��
							{
								nbrTriIdx = (trisLabel[relaOppTrisIdx[0]] == currentLabel ? relaOppTrisIdx[1] : relaOppTrisIdx[0]);
								assert(trisLabel[relaOppTrisIdx[0]] != trisLabel[relaOppTrisIdx[1]] && "Exceptional non-manifold edge detected!");
							}
							break;
						}

					case BOOLEAN_TYPE::UNION:
						{
							// ��nbrTriIdx��������Ƭ����startIdx����Ƭ��ǩ��ͬ������ǩ����Ķ�������Ƭ����ɢ��
							nbrTriIdx = (trisLabel[relaOppTrisIdx[0]] == currentLabel ? relaOppTrisIdx[1] : relaOppTrisIdx[0]);
							assert(trisLabel[relaOppTrisIdx[0]] != trisLabel[relaOppTrisIdx[1]] && "Exceptional non-manifold edge detected!");
							break;
						}

					default:
						assert("invalid BOOLEAN_TYPE!");
					}
				}

				// wf2. ��Ƿ���֮�������Ƭ����β�����µ�����Ƭ
				if (!visited[nbrTriIdx])
				{
					visited[nbrTriIdx] = true;
					triIdxDeq.push_back(nbrTriIdx);
				}
			}
		}

		// 3. �ҳ��������������б���ǵ�����Ƭ�������Ƿ����ʼ����Ƭ�����ͬ��
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

		// 4. ѡȡ����Ƭ���ɽ������
		switch (type)
		{
		case BOOLEAN_TYPE::DIFF:
			{
				// d1. ��ȡ����ǵ�����Ƭ�������Ƿ����ʼ����Ƭ�����ͬ��
				subFromIdxVec(selectedTris1, tris, selectedTrisIdx1);
				subFromIdxVec(selectedTris2, tris, selectedTrisIdx2);

				// d2. ����ʼ����Ƭ��ǲ�ͬ����flip������Ƭ
				Eigen::VectorXi tmpVec = selectedTris2.col(2);
				selectedTris2.col(2) = selectedTris2.col(1);
				selectedTris2.col(1) = tmpVec;
				trisOut = selectedTris1;
				matInsertRows(trisOut, selectedTris2);

				// d3. ȥ���������㣬��������Ƭ�е�������
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
				// u1. ��ȡ����ǵ�����Ƭ��
				subFromIdxVec(selectedTris, tris, selectedTrisIdx);
				trisOut = selectedTris;

				// u2. ȥ���������㣬��������Ƭ�е�������
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

		igl::readOBJ("E:/����/meshArranged.obj", vers, tris);
		int trisCount = tris.rows();

		// ��ȡtrisLabel
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


	// 33. ����buildAdjacency
	void test33()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		igl::readOBJ("E:/����/meshArrangeResult.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();

		Eigen::MatrixXi ttAdj_mnEdge;
		std::vector<tVec> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		bool retFlag = buildAdjacency(ttAdj_mnEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge, tris);

		// ������������αߵ�����Ƭ��
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

		Eigen::MatrixXi mnTris;				// �����������αߵ�����Ƭ��
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


	//		����AABB������դ�񣻡���from libigl
	template <typename Scalar, typename DerivedV,	typename DerivedI>
	bool genGrids(const Eigen::AlignedBox<Scalar, 3>& box, const int largestCount,	const int pad_count,\
			Eigen::PlainObjectBase<DerivedV>& gridCenters, Eigen::PlainObjectBase<DerivedI>& gridCounts)
	{
		/*bool genGrids(																		�ɹ�����true
				const Eigen::AlignedBox<Scalar, 3>&box,						�����AABB����
				const int largestCount,													��������Ǹ�ά��(xyz�е�һ��)��դ������		
				const int pad_count,														������Χ�б߽��դ������
				Eigen::PlainObjectBase<DerivedV>&gridCenters,			դ�����ģ�
				Eigen::PlainObjectBase<DerivedI>&gridCounts				xyz����ά����դ���������
				)
		*/
		using namespace Eigen;
		using namespace std;

		// 1. �����Χ�жԽ��������е���������
		gridCounts.resize(1, 3);
		typename DerivedV::Index maxCompIdx = -1;            // ��Χ��box�ĶԽ������������ķ�����������0,1,2�ֱ��Ӧxyz������
		box.diagonal().maxCoeff(&maxCompIdx);
		const Scalar maxComp = box.diagonal()(maxCompIdx);          // ��Χ��box�ĶԽ������������ķ�����
		assert(largestCount > (pad_count * 2 + 1) && "largestCount should be > 2*pad_count+1");


		// 2. ����xyz����ά����դ��ĸ���gridCounts
		const Scalar largestCount0 = largestCount - 2 * pad_count;
		gridCounts(maxCompIdx) = largestCount0;
		for (int i = 0; i < 3; i++)
			if (i != maxCompIdx)
				gridCounts(i) = std::ceil(largestCount0 * (box.diagonal()(i)) / maxComp);
		gridCounts.array() += 2 * pad_count;

		// 3. ����gridCenters;
		grid(gridCounts, gridCenters);            // ����������ԭ���gridCenters;

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
 
 
	//		marchingCubes�㷨������һ��cube������Ƭ��
	template <typename DerivedGV, typename Scalar, typename Index, typename ScalarV, typename IndexF>
	void handleCube(const DerivedGV& gridCenters, const Eigen::Matrix<Scalar, 8, 1>& cornerSDF, \
		const Eigen::Matrix<Index, 8, 1>& cornerIdx, const Scalar& isovalue, \
		Eigen::Matrix<ScalarV, Eigen::Dynamic, Eigen::Dynamic>& versResult, Index& curVersCount, \
		Eigen::Matrix<IndexF, Eigen::Dynamic, Eigen::Dynamic>& trisResult, Index& curTrisCount, \
		std::unordered_map<int64_t, int>& edgeIsctMap)
	{
		/*
			const DerivedGV& gridCenters,															դ������
			const Eigen::Matrix<Scalar, 8, 1>& cornerSDF,									��ǰ������˸������SDFֵ
			const Eigen::Matrix<Index, 8, 1>& cornerIdx,									��ǰ������˸�������դ���е�������
			const Scalar& isovalue,																		��Ҫ��ȡ�ĵ�ֵ���SDFֵ
			Eigen::PlainObjectBase<DerivedV>& versResult,								�������Ķ���
			Index& curVersCount,																			��ǰ�ۼ����ɵ�������񶥵���
			Eigen::PlainObjectBase<DerivedF>& trisResult,									������������Ƭ
			Index& curTrisCount,																			��ǰ�ۼ����ɵ������������Ƭ��
			std::unordered_map<int64_t, int>& edgeIsctMap								�߱���-�߽��������Ĺ�ϣ��

		*/

		Eigen::Matrix<Index, 12, 1> isctVerIdxes;		// ��������ϵĽ���ľ���������������������������е�������
		int cornerState = 0;											// �����嶥��״̬���룻256�����Σ�

		// ��������߱���
		const auto genMCedgeCode = [](int32_t vaIdx, int32_t vbIdx)
		{
			if (vaIdx > vbIdx)
				std::swap(vaIdx, vbIdx);
			std::int64_t edgeCode = 0;
			edgeCode |= vaIdx;
			edgeCode |= static_cast<std::int64_t>(vbIdx) << 32;
			return edgeCode;
		};

		// 1. ���㵱ǰ������Ķ���״̬���룬��8�������ڵ�ֵ�������״̬1��
		for (int i = 0; i < 8; i++)
			if (cornerSDF(i) > isovalue)
				cornerState |= 1 << i;

		// 2. ȷ����ǰ�������к͵�ֵ���ཻ�ıߣ�
		int edgeState = MC_TABLES::edgeStateCodes[cornerState];		// �����嶥��״̬����ӳ��Ϊ�ཻ�߱��룻
		if (edgeState == 0)
			return;															// ��ʾ��ǰ���������嶼�ڵ�ֵ���ⲿ���ڲ���û�н��㣻

		// 3. ȷ����ֵ��͵�ǰ������ıߵĽ��㣻 Find the point of intersection of the surface with each edge. Then find the normal to the surface at those points
		for (int i = 0; i < 12; i++)						// �����������бߵı���
		{
			if (edgeState & (1 << i))					// ����ֵ��͵�ǰ���ཻ��
			{
				int vaIdxRela = MC_TABLES::cubeEdges[i][0];			// ��ǰ�����˵�����������
				int vbIdxRela = MC_TABLES::cubeEdges[i][1];

				// ���ɱ��ϵĶ��㣺
				int vaIdx = cornerIdx(vaIdxRela);				// ��ǰ�����˵�ľ�����������դ���еĶ���������
				int vbIdx = cornerIdx(vbIdxRela);
				std::int64_t edgeCode = genMCedgeCode(vaIdx, vbIdx);
				const auto iter = edgeIsctMap.find(edgeCode);

				if (iter == edgeIsctMap.end())								// ����ǰ�߽���δ���ɣ�
				{
					if (curVersCount == versResult.rows())
						versResult.conservativeResize(versResult.rows() * 2 + 1, versResult.cols());

					// ��ֵ�����µĶ��㣺find crossing point assuming linear interpolation along edges
					const Scalar& SDFa = cornerSDF(vaIdxRela);			// ��ǰ�����˵��SDFֵ��
					const Scalar& SDFb = cornerSDF(vbIdxRela);
					const Scalar delta = SDFb - SDFa;
					Scalar t = (isovalue - SDFa) / delta;
					versResult.row(curVersCount) = (gridCenters.row(vaIdx) + t * (gridCenters.row(vbIdx) - gridCenters.row(vaIdx))).array().cast<ScalarV>();

					isctVerIdxes[i] = curVersCount;
					edgeIsctMap[edgeCode] = isctVerIdxes[i];
					curVersCount++;
				}
				else                                                                             // ����ǰ�߽��������ɣ�
					isctVerIdxes[i] = iter->second;

				assert(isctVerIdxes[i] >= 0);
				assert(isctVerIdxes[i] < curVersCount);
			}
		}


		// 4. ���ɵ�ǰ�������е�����Ƭ��һ�����������������5������Ƭ��
		for (int i = 0; i < 5; i++)
		{
			if (MC_TABLES::cubeTriangles[cornerState][3 * i] < 0)
				break;

			if (curTrisCount == trisResult.rows())
				trisResult.conservativeResize(trisResult.rows() * 2 + 1, trisResult.cols());

			// ��������Ƭ�����еĶ��������������������
			int vaIdxRela = MC_TABLES::cubeTriangles[cornerState][3 * i + 0];
			int vbIdxRela = MC_TABLES::cubeTriangles[cornerState][3 * i + 1];
			int vcIdxRela = MC_TABLES::cubeTriangles[cornerState][3 * i + 2];

			assert(isctVerIdxes[vaIdxRela] >= 0);
			assert(isctVerIdxes[vbIdxRela] >= 0);
			assert(isctVerIdxes[vcIdxRela] >= 0);

			// �������ת��Ϊ����������������������Ƭ
			trisResult.row(curTrisCount) << isctVerIdxes[vaIdxRela], isctVerIdxes[vbIdxRela], isctVerIdxes[vcIdxRela];
			curTrisCount++;
		}

	}


	//		marchingCubes()��������������Ҫ�����ȥ���ظ����㣬̫�̵ıߡ�
	template <typename DerivedS, typename DerivedGV, typename ScalarV, typename IndexF>
	bool marchingCubes(Eigen::Matrix<ScalarV, Eigen::Dynamic, Eigen::Dynamic>& versResult, \
		Eigen::Matrix<IndexF, Eigen::Dynamic, Eigen::Dynamic>& trisResult, \
		const Eigen::MatrixBase<DerivedS>& scalarFied, const Eigen::MatrixBase<DerivedGV>& gridCenters, \
		const unsigned nx, const unsigned ny, const unsigned nz,
		const typename DerivedS::Scalar isovalue)
	{
		/*
			DerivedS				���볡���ݵ�����
			DerivedGV			դ�����ݵ�����
			ScalarV					����������������ͣ�
			IndexF					��Ƭ�еĶ����������������ͣ�
		

			const Eigen::MatrixBase<DerivedS>& scalarFied,							���ž��볡����
			const Eigen::MatrixBase<DerivedGV>& gridCenters,					դ������
			const unsigned nx,																			x������դ�����
			const unsigned ny,
			const unsigned nz,
			const typename DerivedS::Scalar isovalue,										��Ҫ��ȡ��ˮƽ����SDFֵ��
			Eigen::PlainObjectBase<DerivedV>& versResult,							������񶥵�
			Eigen::PlainObjectBase<DerivedF>& trisResult								�����������Ƭ
		*/

		typedef typename DerivedS::Scalar Scalar;

		// lambda����դ�����ά����ӳ�䵽һά������
		const auto getGridIdx = [&nx, &ny, &nz](const int& x, const int& y, const int& z)->unsigned
		{
			return x + nx * (y + ny * (z));
		};

		const unsigned cornerIdxOffset[8] = { 0, 1, 1 + nx, nx, nx * ny, 1 + nx * ny, 1 + nx + nx * ny, nx + nx * ny };	// ������˸����������ƫ������
		std::unordered_map<int64_t, int> edgeIsctMap;							// �߱���-�߽��������Ĺ�ϣ��

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
					// 1.1 ���㵱ǰդ���������
					const unsigned gridIdx = getGridIdx(x, y, z);

					// 1.2 ���㵱ǰդ���Ӧ��������İ˸���������ݣ�
					static Eigen::Matrix<Scalar, 8, 1> cornerSDF;						// ����İ˸������SDFֵ
					static Eigen::Matrix<unsigned, 8, 1> cornerIdx;               // ����İ˸�������դ���е�����
					for (int i = 0; i < 8; i++)
					{
						const unsigned originIdx = gridIdx + cornerIdxOffset[i];
						cornerIdx(i) = originIdx;
						cornerSDF(i) = scalarFied(originIdx);
					}

					// 1.3 ���ɵ�ǰ�������ڵ�����Ƭ
					handleCube(gridCenters, cornerSDF, cornerIdx, isovalue, versResult, curVersCount, trisResult, curTrisCount, edgeIsctMap);
				}
			}
		}

		// 2. shrink_to_fit();
		versResult.conservativeResize(curVersCount, 3);
		trisResult.conservativeResize(curTrisCount, 3);


		// 3. ���ܻ���duplicated vertices:




		return true;
	}


	//		���볡���ݵĸ�˹�˲���
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


	// test4������������դ��marchingCubes�㷨��
	void test4()
	{
		tiktok& tt = tiktok::getInstance();
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		std::vector<int> stepCounts(3);				// xyz����ά����դ����
		Eigen::RowVector3d gridsOri;					// դ��ԭ�㣺
		Eigen::VectorXd SDF;
		Eigen::RowVector3i gridCounts;			// xyz����ά���ϵ�դ��������
		Eigen::MatrixXd gridCenters;				//	 ����դ���е�����ľ���ÿ�ж���һ���е����ꣻ�洢���ȼ���x, y, z
		Eigen::MatrixXd	boxVers;				// դ���Ӧ�İ�Χ�еĶ��㣻
		Eigen::MatrixXi boxTris;
		double selectedSDF = 0.75;						// ��ȡ��ˮƽ����Ӧ�ľ��볡ֵ��

		// 0. ����SDFGen.exe���ɵ�.sdf���볡�����ļ���
		std::string fileName = "E:/����/jawMesh6";
		std::string sdfFilePath = fileName + std::string{ ".sdf" };
		std::string objFilePath = fileName + std::string{ ".obj" };
		double SDFstep = IGL_BASIC::parseSDF(stepCounts, gridsOri, SDF, sdfFilePath.c_str());
		int xCount = stepCounts[0];
		int yCount = stepCounts[1];
		int zCount = stepCounts[2];
		igl::readOBJ(objFilePath, vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		// 0. ���볡��˹�˲�	��
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

		// 1. ����դ��
		Eigen::RowVector3d minp = gridsOri;
		Eigen::RowVector3d maxp = gridsOri + SDFstep * Eigen::RowVector3d(stepCounts[0] - 1, stepCounts[1] - 1, stepCounts[2] - 1);
		Eigen::AlignedBox<double, 3> box(minp, maxp);		// դ���Ӧ�İ�Χ�У�
		genGrids(box, std::max({ stepCounts[0], stepCounts[1], stepCounts[2] }), 0, gridCenters, gridCounts);
		genAABBmesh(box, boxVers, boxTris);
		objWriteMeshMat("E:/AABB.obj", boxVers, boxTris);

		{
			//			դ�����ݵķֲ���
			/*
				�������������е�դ�����ĵ�Ϊ��
				gc(000), gc(100), gc(200), gc(300),...... gc(010), gc(110), gc(210), gc(310),...... gc(001), gc(101), gc(201).....


				x���꣺
				x0, x1, x2, x3......x0, x1, x2, x3......x0, x1, x2, x3......
				����ΪxCount;
				�ظ�����Ϊ(yCount * zCount)

				y���꣺
				y0, y0, y0...y1, y1, y1...y2, y2, y2.........y0, y0, y0...y1, y1, y1...
				����Ϊ(xCount * yCount);
				�ظ�����ΪzCount;
				����Ԫ���ظ�����ΪxCount

				z���꣺
				z0, z0, z0......z1, z1, z1......z2, z2, z2......
				����Ԫ���ظ�����Ϊ(xCount * yCount)
			*/
			Eigen::MatrixXd gridCenters0;				// for try���������Լ�����դ�����ݣ�
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

			//				��ȡդ����SDFֵС��0�Ķ��㣺
			Eigen::MatrixXd tmpVers(SDF.rows(), 3);
			int index = 0;
			for (int i = 0; i < SDF.rows(); ++i)
				if (SDF(i) <= 0)
					tmpVers.row(index++) = gridCenters0.row(i);
			tmpVers.conservativeResize(index, 3);
			igl::writeOBJ("E:/tmpVers.obj", tmpVers, Eigen::MatrixXi{});
		}

		// 2. marching cubes�㷨�����������棺
		Eigen::MatrixXd versResult_SDF, versResults_signs, versResult_origin;
		Eigen::MatrixXi trisResult_SDF, trisResults_signs, trisResult_origin;
		tt.start();
		marchingCubes(versResult_SDF, trisResult_SDF, SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), selectedSDF);
		tt.endCout("Elapsed time of my marching_cubes() is ");
		igl::writeOBJ("E:/shrinkedMesh.obj", versResult_SDF, trisResult_SDF);

		// ԭʼ��marching cubes
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
		Eigen::RowVector3i gridCounts;			// xyz����ά���ϵ�դ��������
		Eigen::MatrixXd gridCenters;				//	 ����դ���е�����ľ���ÿ�ж���һ���е����ꣻ�洢���ȼ���x, y, z
		Eigen::MatrixXd	boxVers;						// դ���Ӧ�İ�Χ�еĶ��㣻
		Eigen::MatrixXi boxTris;
		double selectedSDF = 0.75;						// ��ȡ��ˮƽ����Ӧ�ľ��볡ֵ��

		std::vector<int> stepCounts{ 126, 132, 39 };														// xyz����ά����դ����
		Eigen::RowVector3d gridsOri{ -23.4242516, -40.1728897, -1.50000501 };					// դ��ԭ�㣺
		double SDFstep = 0.5;
		int xCount = stepCounts[0];
		int yCount = stepCounts[1];
		int zCount = stepCounts[2];
 
		// 0. ����դ��
		Eigen::RowVector3d minp = gridsOri;
		Eigen::RowVector3d maxp = gridsOri + SDFstep * Eigen::RowVector3d(stepCounts[0] - 1, stepCounts[1] - 1, stepCounts[2] - 1);
		Eigen::AlignedBox<double, 3> box(minp, maxp);		// դ���Ӧ�İ�Χ�У�
		genGrids(box, std::max({ stepCounts[0], stepCounts[1], stepCounts[2] }), 0, gridCenters, gridCounts);
		genAABBmesh(box, boxVers, boxTris);
		objWriteMeshMat("E:/AABB.obj", boxVers, boxTris);
		{
			//			դ�����ݵķֲ���
			/*
				�������������е�դ�����ĵ�Ϊ��
				gc(000), gc(100), gc(200), gc(300),...... gc(010), gc(110), gc(210), gc(310),...... gc(001), gc(101), gc(201).....


				x���꣺
				x0, x1, x2, x3......x0, x1, x2, x3......x0, x1, x2, x3......
				����ΪxCount;
				�ظ�����Ϊ(yCount * zCount)

				y���꣺
				y0, y0, y0...y1, y1, y1...y2, y2, y2.........y0, y0, y0...y1, y1, y1...
				����Ϊ(xCount * yCount);
				�ظ�����ΪzCount;
				����Ԫ���ظ�����ΪxCount

				z���꣺
				z0, z0, z0......z1, z1, z1......z2, z2, z2......
				����Ԫ���ظ�����Ϊ(xCount * yCount)
			*/
			Eigen::MatrixXd gridCenters0;				// for try���������Լ�����դ�����ݣ�
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

			//				��ȡդ����SDFֵС��0�Ķ��㣺
			Eigen::MatrixXd tmpVers(SDF.rows(), 3);
			int index = 0;
			for (int i = 0; i < SDF.rows(); ++i)
				if (SDF(i) <= 0)
					tmpVers.row(index++) = gridCenters0.row(i);
			tmpVers.conservativeResize(index, 3);
			igl::writeOBJ("E:/tmpVers.obj", tmpVers, Eigen::MatrixXi{});
		}


		// 1. ����SDF���ݣ�
		unsigned spVersCount = gridCenters.rows();
		SDF.resize(spVersCount);
		for (unsigned i = 0; i < spVersCount; ++i)
			SDF(i) = gridCenters(i, 2);

		gaussFilterSDFdata(SDF, xCount, yCount, zCount);

		// 2ԭʼ��marching cubes
		Eigen::MatrixXd versResult;
		Eigen::MatrixXi trisResult;
		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), 0, versResult, trisResult);
		igl::writeOBJ("E:/meshXOY.obj", versResult, trisResult);

		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), 1.7, versResult, trisResult);
		igl::writeOBJ("E:/meshSDF1.7.obj", versResult, trisResult);


		std::cout << "finished." << std::endl;
	}


	// ������opological_hole_fill()������ǰ������
	void test5()
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris, bdrys;
		igl::readOBJ("E:/����/holeMesh.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);		

		const unsigned versCount = vers.rows();

		// ���㶴�����򶥵�����
		std::vector<std::vector<int>> holes, bdrySegs;
		findHolesBdrySegs(holes, bdrySegs, vers, tris);

		// ��Ҫ�ֶ�Ϊԭ������Ӷ������ĵ㣬����ò�ƴ������������λ����
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


	// ������winding number:
	void test6() 
	{
		Eigen::MatrixXd vers0, vers1;
		Eigen::MatrixXi tris0, tris1;
		igl::readOBJ("E:/����/tooth.obj", vers0, tris0);
		igl::readOBJ("E:/����/cylinder1.obj", vers1, tris1);
		Eigen::VectorXd wNums;

		igl::winding_number(vers0, tris0, vers1, wNums);			// �����ⲿ�������Ϊ0��

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


	// ����������ͨ������ȡconnected_components()
	void test7() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		objReadMeshMat(vers, tris, "E:/����/ԭ������/originalMesh.obj");
		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();


		// 1. �����ڽӾ���
		Eigen::SparseMatrix<int> adjSM_eCount, adjSM_eIdx;
		adjMatrix(adjSM_eCount, adjSM_eIdx, tris);

		Eigen::SparseMatrix<int> adjSM = adjSM_eCount;
		traverseSparseMatrix(adjSM, [&](auto& iter)
			{
				iter.valueRef() = 1;
			});

		// 2. ȷ������ͨ���򡪡�connected_components()
		Eigen::VectorXi connectedLabels, connectedCount;
		int conCount = igl::connected_components(adjSM, connectedLabels, connectedCount);

		// 3. ��ȡ���ĵ���ͨ�����еĶ��㣺
		std::vector<int> retVec1, retVec2;
		retVec1 = vec2Vec(connectedLabels);
		retVec2 = vec2Vec(connectedCount);

		int mainLabel = 0;						// ���������ͨ�������ı�ǩ��
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

		// 4. ��ȡ�����ͨ�����е�����Ƭ��
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


	// �Լ�ʵ�ֵ�connected_conponents()
	void test77()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		objReadMeshMat(vers, tris, "E:/����/ԭ������/originalMesh.obj");
		objWriteMeshMat("E:/meshInput.obj", vers, tris);

		bool retFlag = simplyConnectedLargest(versOut, trisOut, vers, tris);
		if (!retFlag)
			std::cout << "function failed!!!" << std::endl;

		objWriteMeshMat("E:/meshOut.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// ������������֮���hausdorff distance
	void test8() 
	{
		Eigen::MatrixXd vers0, vers1, vers2, vers3, vers4;
		Eigen::MatrixXi tris0, tris1, tris2, tris3, tris4;
		double hd0, hd1, hd2, hd3, hd4;
		hd0 = hd1 = hd2 = hd3 = hd4 = 0;

		tiktok& tt = tiktok::getInstance();
		tt.start();
		objReadMeshMat(vers0, tris0, "E:/����/jawMeshDense.obj");
		objReadMeshMat(vers1, tris1, "E:/����/jawMeshDenseSimplifyOutIsctCleaned.obj");
		objReadMeshMat(vers2, tris2, "E:/����/jawMeshDenseQEMoutputIsctCleaned.obj");
		objReadMeshMat(vers3, tris3, "E:/����/jawMeshDenseQslimOutputIsctCleaned.obj");
		objReadMeshMat(vers4, tris4, "E:/����/jawMeshDensePMPsimplifiedIsctCleaned.obj");
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
