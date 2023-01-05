#include "igl_study.h"
#define DATA_PATH "./data/"

igl::opengl::glfw::Viewer viewer;				// libigl�еĻ���glfw����ʾ���ڣ�
static std::mutex g_mutex;


// libigl��������
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

		VectorXd U;
		igl::readDMAT("./data/cheburashka-scalar.dmat", U);		// һϵ�еĺ���ֵ

		SparseMatrix<double> G;			// �ݶ�����
		igl::grad(vers, tris, G);

		// ���㺯��ֵU���ݶ�
		MatrixXd GU = Map<const MatrixXd>((G * U).eval().data(), tris.rows(), 3);

		// ����ֵU���ݶ�GU��ģ��
		const VectorXd GU_mag = GU.rowwise().norm();

		viewer.data().set_mesh(vers, tris);
		viewer.data().set_data(U);

		// Average edge length divided by average gradient (for scaling)
		const double max_size = igl::avg_edge_length(vers, tris) / GU_mag.mean();

		// ÿ������Ƭ�����ϻ�һ��ָʾ�ߣ�����Ϊ�ݶȷ��� 
		MatrixXd BC;
		igl::barycenter(vers, tris, BC);
		const RowVector3d black(0, 0, 0);
		viewer.data().add_edges(BC, BC + max_size * GU, black);
		viewer.data().show_lines = false;	  // ����������

		viewer.launch();
	}
 

	// lambda���������¼���ʹ��laplacian��˳����
	bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)
	{
		switch (key)
		{
		case 'r':

		case 'R':			// ��λ����
			newVers = vers;
			break;

		case ' ':				// �ո����ִ��һ��laplace��˳
		{
			// ���¼�����������
			Eigen::SparseMatrix<double> mass;
			igl::massmatrix(newVers, tris, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

			// �����Է����� (mass - delta*L) * newVers = mass * newVers
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


	// ʹ��Laplacian��˳����
	void test3() 
	{
		 igl::readOBJ( "./data/bunny.obj", vers, tris);
		newVers = vers;

		// 1.a ֱ�ӹ���laplacian����Compute Laplace-Beltrami operator: 
		igl::cotmatrix(vers, tris, L);

		// 1.b �ֲ�����laplacian
		{
			SparseMatrix<double> Gradient, L2;

			igl::grad(vers, tris, Gradient);      // ��ɢ�ݶ�

			// Diagonal per-triangle "mass matrix"
			VectorXd dblA;
			igl::doublearea(vers, tris, dblA);             // ÿ������Ƭ���������

			// Place areas along diagonal  #dim times
			const auto& T = 1. * (dblA.replicate(3, 1) * 0.5).asDiagonal();

			L2 = -Gradient.transpose() * T * Gradient;         // discrete Dirichelet energy Hessian ��ɢ��������������������
			std::cout << "���ַ����õ���laplacian�Ĳ�ķ�����" << std::endl;
			cout << "(L2 - L).norm() == " << (L2 - L).norm() << endl;
		}

		// 2. ����ԭʼ�ķ�������ʹ��αɫ
		MatrixXd norms;
		igl::per_vertex_normals(vers, tris, norms);
		MatrixXd colors = norms.rowwise().normalized().array() * 0.5 + 0.5;

		// 3. viewr����ʼ����
		newVers = vers;
		viewer.data().set_mesh(newVers, tris);
		viewer.data().set_colors(colors);
		viewer.callback_key_down = key_down;

		// 4. ����
		cout << "Press [space] to smooth." << endl;;
		cout << "Press [r] to reset." << endl;;
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
		MatrixXi tris;
		MatrixXd vers;
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
		VectorXd SDF, signValues;
		{
			VectorXi I;                     // useless
			MatrixXd C, N;              // useless

			// 2. ������ž��볡
			tt.start();
			igl::signed_distance(gridCenters, vers, tris, igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, SDF, I, C, N);
			tt.endCout("Elapsed time of igl::signed_distance() is ");

			// 3. ���ž��볡��дΪ���ų�����������Ϊ-1����������Ϊ0������Ϊ1��
			signValues = SDF;
			for_each(signValues.data(), signValues.data() + signValues.size(), [](double& b)\
			{
				b = (b > 0 ? 1 : (b < 0 ? -1 : 0));
			});
		}

		// 4. marching cubes�㷨�����������棺
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


	// ��ȡSDFGen.exe���ɵ�.sdf���볡���ݣ�ʹ��igl::marching_cubes()��ȡ��ֵ������

	//		����.sdf�ı��ļ���
	double parseSDF(std::vector<int>& stepCounts, Eigen::RowVector3d& gridsOri, Eigen::VectorXd& SDF, const char* filePath)
	{
		double SDFstep = -1;
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
		MatrixXd versResult_SDF, versResults_signs;
		MatrixXi trisResult_SDF, trisResults_signs;
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
 
		igl::writeOBJ("E:/426.obj", MatrixXd{ RowVector3d(4, 2, 6) }, MatrixXi{});

		std::cout << "finished." << std::endl;
	}


	// decimate()��qslim() ֱ��ʵ�����񾫼�
	void test7() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]�Ǿ����������е�i������Ƭ��Ӧ��ԭ���������Ƭ������
		Eigen::VectorXi newOldVersInfo;						 
		tiktok& tt = tiktok::getInstance();

		std::string fileName = "tmpMesh";

		tt.start();
		igl::readOBJ((std::string{ "E:/����/" } + fileName + std::string{".obj"}).c_str(), vers, tris);
		tt.endCout("elapsed time of loading mesh is: ");

		unsigned trisCount = tris.rows();
		unsigned tarTrisCount = std::round(trisCount * 0.5);			// ���������������Ƭ����
		// unsigned tarTrisCount = 7516;
		igl::writeOBJ("E:/meshIn.obj", vers, tris);

#if 1
		// igl::decimate()�������۵��㷨�������� ���񲻿����з����αߣ�				
		tt.start();
		std::cout << "succeeded? " << igl::decimate(vers, tris, tarTrisCount, versOut, trisOut, \
			newOldTrisInfo, newOldVersInfo) << std::endl;						// ����1.1
		tt.endCout("Elapsed time of mesh simplification by igl::decimate() is ");
		igl::writeOBJ((std::string{ "E:/" } + fileName + std::string{ "_simplified_edgeCollapse.obj" }).c_str(), versOut, trisOut);

		std::vector<int> newOldTrisInfoVec = vec2Vec(newOldTrisInfo);
		Eigen::MatrixXd vers1;
		subFromIdxVec(vers1, vers, newOldVersInfo);
#endif

		// igl::qslim()����ò�Ʊ�igl::decimate()��һЩ��Ŀǰ��֪�����������
		tt.start();
		std::cout << "qslim succeeded? " << igl::qslim(vers, tris, tarTrisCount, versOut, trisOut, newOldTrisInfo, newOldVersInfo) << std::endl;
		tt.endCout("Elapsed time of mesh simplification by igl::qslim() is ");
		igl::writeOBJ((std::string{ "E:/" } + fileName + std::string{ "_simplified_qslim.obj" }).c_str(), versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// collapse_edge()�������񾫼��еı��۵��㷨��
	void test77() 
	{
		// �۵�һ���ߣ�
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		igl::readOBJ("E:/����/cylinder1.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		// ���ɱ����ݣ�
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
		
		// 1. igl::collapse_edge(); �۵�����ΪuEdgeIdx0�ıߣ�
		int e1, e2, f1, f2;
		igl::collapse_edge(uEdgeIdx0, collapsedVer, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, e1, e2, f1, f2);
		igl::writeOBJ("E:/meshOut_ɾ������Ƭǰ.obj", vers, tris);

		// 2. ɾ�����к��б��ΪIGL_COLLAPSE_EDGE_NULL�ߵ�����Ƭ��
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
		tris0.conservativeResize(m, tris0.cols());					// �����൱��shrink_to_fit();
		igl::writeOBJ("E:/meshOut_beforeRemove.obj", vers, tris0);

		// 3. igl::remove_unreferenced()����ɾ�������еĹ������㣻
		Eigen::VectorXi I, newOldVersInfo;
		igl::remove_unreferenced(vers, tris0, versOut, trisOut, I, newOldVersInfo);
		igl::writeOBJ("E:/meshOut.obj", versOut, trisOut);


		std::cout << "finished." << std::endl;
	}


	bool edge_collapse_is_valid(std::vector<int>& srcNbrIdx, 	std::vector<int>& desNbrIdx)
	{
		// Do we really need to check if edge is IGL_COLLAPSE_EDGE_NULL ?
		if (srcNbrIdx.size() < 2 || desNbrIdx.size() < 2)
		{
			// Bogus data
			assert(false);
			return false;
		}

		// determine if the first two vertices are the same before reordering.
		// If they are and there are 3 each, then (I claim) this is an edge on a single tet.
		const bool first_two_same = (srcNbrIdx[0] == desNbrIdx[0]) && (srcNbrIdx[1] == desNbrIdx[1]);
		if (srcNbrIdx.size() == 3 && desNbrIdx.size() == 3 && first_two_same)
			return false;           // single tet


		  // https://stackoverflow.com/a/19483741/148668
		std::sort(srcNbrIdx.begin(), srcNbrIdx.end());
		std::sort(desNbrIdx.begin(), desNbrIdx.end());
		std::vector<int> Nint;
		std::set_intersection(srcNbrIdx.begin(), srcNbrIdx.end(), desNbrIdx.begin(), desNbrIdx.end(), std::back_inserter(Nint));

		// check if edge collapse is valid: intersection of vertex neighbors of s and d should be exactly 2+(s,d) = 4

		// http://stackoverflow.com/a/27049418/148668
		if (Nint.size() != 2)
			return false;


		return true;
	}


	bool collapseSingleEdge(
		const int uEdgeIdx,
		const Eigen::RowVectorXd& collapsedVer,
		std::vector<int>& nbrVersIdx_src,
		const std::vector<int>& nbrTrisIdx_src,
		std::vector<int>& nbrVersIdx_des,
		const std::vector<int>& nbrTrisIdx_des,
		Eigen::MatrixXd& vers,
		Eigen::MatrixXi& tris,
		Eigen::MatrixXi& uEdges,
		Eigen::VectorXi& edgeUeInfo,
		Eigen::MatrixXi& UeTrisInfo,
		Eigen::MatrixXi& UeCornersInfo,
		int& a_e1,
		int& a_e2,
		int& a_f1,
		int& a_f2)
	{
		/*
		   Assign this to 0 rather than,  say,  -1 so that deleted elements will get draw as degenerate elements at vertex 0
				  (which should always exist and never get collapsed to anything else since it is the smallest index)
	   */
		using namespace Eigen;
		using namespace std;
		const int eFlipFlag = uEdges(uEdgeIdx, 0) > uEdges(uEdgeIdx, 1);

		// lambda����ĳһ���ߵ�������Ϣ����Ϊ��Ч��Ϣ��
		const auto& kill_edge = [&uEdges, &UeCornersInfo, &UeTrisInfo](const int uEdgeIdx)
		{
			uEdges(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			uEdges(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
			UeTrisInfo(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			UeTrisInfo(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
			UeCornersInfo(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			UeCornersInfo(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
		};

		// ��������ߵ������յ㡪������С��Ϊ��㣬�������Ϊ�յ㣻
		const int srcIdx = eFlipFlag ? uEdges(uEdgeIdx, 1) : uEdges(uEdgeIdx, 0);
		const int desIdx = eFlipFlag ? uEdges(uEdgeIdx, 0) : uEdges(uEdgeIdx, 1);

		// 0. �жϸñ��Ƿ���۵�
		if (!edge_collapse_is_valid(nbrVersIdx_src, nbrVersIdx_des))
			return false;

		// Important to grab neighbors of desIdx before monkeying with uEdges
		const std::vector<int>& nV2Fd = (!eFlipFlag ? nbrTrisIdx_src : nbrTrisIdx_des);

		assert(srcIdx < desIdx && "srcIdx should be less than desIdx");			 // The following implementation strongly relies on srcIdx<desIdx

		// 1. �����˵㶼�滻ΪcollapsedVer
		vers.row(srcIdx) = collapsedVer;
		vers.row(desIdx) = collapsedVer;

		// update edge info for each flap
		const int trisCount = tris.rows();
		for (int side = 0; side < 2; side++)
		{
			const int f = UeTrisInfo(uEdgeIdx, side);
			const int v = UeCornersInfo(uEdgeIdx, side);
			const int sign = (eFlipFlag == 0 ? 1 : -1) * (1 - 2 * side);

			// next edge emanating from desIdx
			const int e1 = edgeUeInfo(f + trisCount * ((v + sign * 1 + 3) % 3));

			// prev edge pointing to srcIdx
			const int e2 = edgeUeInfo(f + trisCount * ((v + sign * 2 + 3) % 3));
			assert(uEdges(e1, 0) == desIdx || uEdges(e1, 1) == desIdx);
			assert(uEdges(e2, 0) == srcIdx || uEdges(e2, 1) == srcIdx);

			// face adjacent to f on e1,  also incident on desIdx
			const bool flip1 = UeTrisInfo(e1, 1) == f;
			const int f1 = flip1 ? UeTrisInfo(e1, 0) : UeTrisInfo(e1, 1);
			assert(f1 != f);
			assert(tris(f1, 0) == desIdx || tris(f1, 1) == desIdx || tris(f1, 2) == desIdx);

			// across from which vertex of f1 does e1 appear?
			const int v1 = flip1 ? UeCornersInfo(e1, 0) : UeCornersInfo(e1, 1);

			// Kill e1
			kill_edge(e1);

			// Kill f
			tris(f, 0) = IGL_COLLAPSE_EDGE_NULL;
			tris(f, 1) = IGL_COLLAPSE_EDGE_NULL;
			tris(f, 2) = IGL_COLLAPSE_EDGE_NULL;

			// map f1'srcIdx edge on e1 to e2
			assert(edgeUeInfo(f1 + trisCount * v1) == e1);
			edgeUeInfo(f1 + trisCount * v1) = e2;

			// side opposite f2,  the face adjacent to f on e2,  also incident on srcIdx
			const int opp2 = (UeTrisInfo(e2, 0) == f ? 0 : 1);
			assert(UeTrisInfo(e2, opp2) == f);
			UeTrisInfo(e2, opp2) = f1;
			UeCornersInfo(e2, opp2) = v1;

			// remap e2 from desIdx to srcIdx
			uEdges(e2, 0) = uEdges(e2, 0) == desIdx ? srcIdx : uEdges(e2, 0);
			uEdges(e2, 1) = uEdges(e2, 1) == desIdx ? srcIdx : uEdges(e2, 1);
			if (side == 0)
			{
				a_e1 = e1;
				a_f1 = f;
			}
			else
			{
				a_e2 = e1;
				a_f2 = f;
			}
		}

		/*
			 finally,  reindex faces and uEdges incident on desIdx. Do this last so asserts make sense.
			 Could actually skip first and last,  since those are always the two collpased faces.
			 Nah,  this is handled by (tris(f, v) == desIdx)
			 Don't attempt to use Nde, Nse here because edgeUeInfo has changed
		*/
		{
			int p1 = -1;
			for (auto f : nV2Fd)
			{
				for (int v = 0; v < 3; v++)
				{
					if (tris(f, v) == desIdx)
					{
						const int e1 = edgeUeInfo(f + trisCount * ((v + 1) % 3));
						const int flip1 = (UeTrisInfo(e1, 0) == f) ? 1 : 0;
						assert(uEdges(e1, flip1) == desIdx || uEdges(e1, flip1) == srcIdx);
						uEdges(e1, flip1) = srcIdx;
						const int e2 = edgeUeInfo(f + trisCount * ((v + 2) % 3));

						// Skip if we just handled this edge (claim: this will be all except for the first non-trivial face)
						if (e2 != p1)
						{
							const int flip2 = (UeTrisInfo(e2, 0) == f) ? 0 : 1;
							assert(uEdges(e2, flip2) == desIdx || uEdges(e2, flip2) == srcIdx);
							uEdges(e2, flip2) = srcIdx;
						}

						tris(f, v) = srcIdx;
						p1 = e1;
						break;
					}
				}
			}
		}

		// Finally,  "remove" this edge and its information
		kill_edge(uEdgeIdx);
		return true;
	}


	// qslim������۵��ߣ�
	void test777()
	{
		using namespace igl;

		using Quadric = std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>;

		Eigen::MatrixXd vers, versOut, vers0, versOri;
		Eigen::MatrixXi tris, trisOut, tris0, trisOri;
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]�Ǿ����������е�i������Ƭ��Ӧ��ԭ���������Ƭ������
		Eigen::VectorXi newOldVersInfo;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
		Eigen::VectorXi timeStamps;
		MatrixXd collapsedVers;
		Eigen::VectorXd costs;
		std::vector<Quadric> quadrics;								// ÿ�������Q����
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// ���ȶ��У�
		Eigen::VectorXi _1, I2;
		bool clean_finish = false;
		int v1 = -1;				// State variables keeping track of edge we just collapsed
		int v2 = -1;

		// 0.
		tiktok& tt = tiktok::getInstance();
		tt.start();
		igl::readOBJ("E:/����/tmpMesh.obj", versOri, trisOri);
		igl::writeOBJ("E:/meshIn.obj", versOri, trisOri);
		tt.endCout("elapsed time of loading mesh is: ");

		// 1. 
		unsigned trisCount = trisOri.rows();
		igl::connect_boundary_to_infinity(versOri, trisOri, vers, tris);
		if (!igl::is_edge_manifold(tris))
			return;
		igl::edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);


		// 2. ����ÿ�������Q����
		igl::per_vertex_point_to_plane_quadrics(vers, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, quadrics);


		// 3.  ��ȡqslim�㷨�еĺ����ӣ�
		int trisCountNew = trisCount;
		int tarTrisCount = std::round(0.5 * trisCount);

		// 4. ����ÿ������ߵ�costֵ���Դ�Ϊ���ȼ��������ȶ���
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());
		collapsedVers.resize(uEdges.rows(), vers.cols());
		costs.resize(uEdges.rows());
		igl::parallel_for(uEdges.rows(), [&](const int ueIdx)
			{
				// ������cost_and_placement()�����ݣ�
				double cost = ueIdx;
				RowVectorXd p(1, 3);

				// Combined quadric
				Quadric quadric_p;
				quadric_p = quadrics[uEdges(ueIdx, 0)] + quadrics[uEdges(ueIdx, 1)];

				// Quadric: p'Ap + 2b'p + c,  optimal point: Ap = -b, or rather because we have row vectors: pA=-b
				const auto& A = std::get<0>(quadric_p);
				const auto& b = std::get<1>(quadric_p);
				const auto& c = std::get<2>(quadric_p);
				p = -b * A.inverse();
				cost = p.dot(p * A) + 2 * p.dot(b) + c;

				// Force infs and nans to infinity
				if (std::isinf(cost) || cost != cost)
				{
					cost = std::numeric_limits<double>::infinity();
					// Prevent NaNs. Actually NaNs might be useful for debugging.
					p.setConstant(0);
				}

				collapsedVers.row(ueIdx) = p;
				costs(ueIdx) = cost;
			},
			10000);

		for (int i = 0; i < uEdges.rows(); i++)
			pQueue.emplace(costs(i), i, 0);

#if 0
		// 5. ��ִ�б��۵���
		for (unsigned i =0; i <= 127; ++i) 
		{
			Eigen::MatrixXd versBC = vers;
			Eigen::MatrixXi trisBC = tris;
			Eigen::MatrixXi uEdgesBC = uEdges;
			auto tuple0 = pQueue.top();
			int uEdgeIdx0 = std::get<1>(tuple0);			// �����۵��ıߵ�������

			// 5.1
			int e, e1, e2, f1, f2;
			if (!collapse_edge(cost_and_placement, pre_collapse, post_collapse, \
				vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, \
				pQueue, timeStamps, collapsedVers, e, e1, e2, f1, f2))          // collapse_edge()����2��pre_collapse()�������۵��ߣ�ִ��post_collapse()���ٸ�����صıߵ�costֵ��
				assert("edge collapse failed.");

			// 5.2. ɾ�����к��б��ΪIGL_COLLAPSE_EDGE_NULL�ߵ�����Ƭ��
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
			tris0.conservativeResize(index, tris0.cols());              // �����൱��shrink_to_fit();
			newOldTrisInfo.conservativeResize(index);

			// 5.3. ɾ�������еĹ������㣺
			igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);

			//		��ӡ��ǰ����Ľ����
			if (i == 125 || i == 126 || i == 127) 
			{
				char str[256];
				sprintf_s(str, 256, "E:/qslim�۵�������loop%d.obj", i);
				igl::writeOBJ(str, versOut, trisOut);

				sprintf_s(str, 256, "E:/loop%d�۵��ı�.obj", i);
				Eigen::MatrixXi uEdge0 = uEdgesBC.row(uEdgeIdx0);
				objWriteEdgesMat(str, uEdge0, versBC);

				// for debug:
				Eigen::VectorXi uEdgeVec0 = uEdge0.col(0);
				Eigen::VectorXi uEdgeVec1 = uEdge0.col(1);
				std::vector<int> tmpVec0 = vec2Vec(uEdgeVec0);
				std::vector<int> tmpVec1 = vec2Vec(uEdgeVec0);

				std::cout << "loop " << i << " finished." << "�۵��ıߣ�" << uEdge0(0, 0) << ", " << uEdge0(0, 1) << std::endl;
			}
		}
#endif

		// 5. ���۵���ѭ����
		int uEdgeIdx0, e1, e2, f1, f2;
		for (unsigned i = 0; i <= 130; ++i)
		{
			// ��һ�ֵ����ݣ�
			Eigen::MatrixXd versBC = vers;
			Eigen::MatrixXi trisBC = tris;
			Eigen::MatrixXi uEdgesBC = uEdges;

			// 5.1. ȡ����Ԫ�أ�����Ԫ�س��ӣ�
			std::tuple<double, int, int> edgeTuple;
			while (true)
			{
				// 1.1 ������Ϊ�գ��˳�ѭ����
				if (pQueue.empty())   // no uEdges to collapse
					assert("���۵�����");

				// 1.2ȡ����Ԫ�أ�����Ԫ�س��ӣ�
				edgeTuple = pQueue.top();
				if (std::get<0>(edgeTuple) == std::numeric_limits<double>::infinity())
					assert("���ױߵ�cost�������");

				pQueue.pop();
				uEdgeIdx0 = std::get<1>(edgeTuple);          // ���׵����������;

				// 1.3 Check if matches timestamp
				if (std::get<2>(edgeTuple) == timeStamps(uEdgeIdx0))
					break;

				// 1.4 ������
				assert(std::get<2>(edgeTuple) < timeStamps(uEdgeIdx0) || timeStamps(uEdgeIdx0) == -1);          // must be stale or dead.
			}

			// 5.2. ���㵱ǰ�����˵�1����Ķ��㡢����Ƭ��
			std::vector<int> nbrTrisIdx_src, nbrVersIdx_src;
			igl::circulation(uEdgeIdx0, true, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_src, nbrTrisIdx_src);
			std::vector<int>  nbrTrisIdx_des, nbrVersIdx_des;
			igl::circulation(uEdgeIdx0, false, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_des, nbrTrisIdx_des);

			//		pre_collapse:
			v1 = uEdges(uEdgeIdx0, 0);
			v2 = uEdges(uEdgeIdx0, 1);

			// 5.3. �۵����׵ıߣ�
			bool collapsed = true;
			collapsed = collapse_edge(uEdgeIdx0, collapsedVers.row(uEdgeIdx0), nbrVersIdx_src, nbrTrisIdx_src, \
				nbrVersIdx_des, nbrTrisIdx_des, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, e1, e2, f1, f2);

			//		post_collapses:
			if (collapsed)
				quadrics[v1 < v2 ? v1 : v2] = quadrics[v1] + quadrics[v2];

			// 5.4. �۵�����֮�󣬸������timeStamp�� ������رߵ�costֵ
			if (collapsed)
			{
				// 4.1 Erase the two,  other collapsed uEdges by marking their timestamps as -1
				timeStamps(e1) = -1;
				timeStamps(e2) = -1;

				// 4.2
				std::vector<int> nbrTrisIdx;
				nbrTrisIdx.reserve(nbrTrisIdx_src.size() + nbrTrisIdx_des.size());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_src.begin(), nbrTrisIdx_src.end());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_des.begin(), nbrTrisIdx_des.end());
				std::sort(nbrTrisIdx.begin(), nbrTrisIdx.end());
				nbrTrisIdx.erase(std::unique(nbrTrisIdx.begin(), nbrTrisIdx.end()), nbrTrisIdx.end());

				// 4.3 Collect all uEdges that must be updated
				std::vector<int> Ne;
				Ne.reserve(3 * nbrTrisIdx.size());
				for (auto& triIdx : nbrTrisIdx)
				{
					if (tris(triIdx, 0) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 1) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 2) != IGL_COLLAPSE_EDGE_NULL)
					{
						for (int i = 0; i < 3; i++)
						{
							const int ueIdx = edgeUeInfo(i * tris.rows() + triIdx);
							Ne.push_back(ueIdx);
						}
					}
				}

				// Only process edge once
				std::sort(Ne.begin(), Ne.end());
				Ne.erase(std::unique(Ne.begin(), Ne.end()), Ne.end());             // ȥ���ظ�Ԫ�أ�
				for (auto& ueIdx : Ne)
				{
					// ������۵���costֵ�����۵���Ķ������꣺
					double cost;
					RowVectorXd place;

					// ������cost_and_placement

					 // Combined quadric
					Quadric quadric_p;
					quadric_p = quadrics[uEdges(ueIdx, 0)] + quadrics[uEdges(ueIdx, 1)];

					// Quadric: place'Ap + 2b'place + c,  optimal point: Ap = -b, or rather because we have row vectors: pA=-b
					const auto& A = std::get<0>(quadric_p);
					const auto& b = std::get<1>(quadric_p);
					const auto& c = std::get<2>(quadric_p);
					place = -b * A.inverse();
					cost = place.dot(place * A) + 2 * place.dot(b) + c;

					// Force infs and nans to infinity
					if (std::isinf(cost) || cost != cost)
					{
						cost = std::numeric_limits<double>::infinity();
						// Prevent NaNs. Actually NaNs might be useful for debugging.
						place.setConstant(0);
					}

					// Increment timestamp
					timeStamps(ueIdx)++;

					// Replace in queue
					pQueue.emplace(cost, ueIdx, timeStamps(ueIdx));
					collapsedVers.row(ueIdx) = place;
				}
			}
			else
				assert("edge collapse failed.");

			// 5.5 ִ��pre_collapse()�������۵��ߣ�ִ��post_collapse()���ٸ�����صıߵ�costֵ��
			if (collapsed)
			{
				// ��stopping_condition�����ӷ���true����������ֹ�����������۵�ѭ��
				if (f1 < trisCount)
					trisCountNew -= 1;
				if (f2 < trisCount)
					trisCountNew -= 1;

				bool stopConditionFlag = (trisCountNew <= (int)tarTrisCount);
				if (stopConditionFlag)
				{
					clean_finish = true;
					break;
				}
			}
			else			 // ���۵�ʧ�ܣ��˳�ѭ����
				assert("edge collapse failed.");

			// for debug����
			
			//		ɾ�����к��б��ΪIGL_COLLAPSE_EDGE_NULL�ߵ�����Ƭ��
			tris0.resize(tris.rows(), 3);
			newOldTrisInfo.resize(tris.rows());
			int index = 0;
			for (int k = 0; k < tris.rows(); k++)
			{
				if (tris(k, 0) != IGL_COLLAPSE_EDGE_NULL ||
					tris(k, 1) != IGL_COLLAPSE_EDGE_NULL ||
					tris(k, 2) != IGL_COLLAPSE_EDGE_NULL)
				{
					tris0.row(index) = tris.row(k);
					newOldTrisInfo(index) = k;
					index++;
				}
			}
			tris0.conservativeResize(index, tris0.cols());              // �����൱��shrink_to_fit();
			newOldTrisInfo.conservativeResize(index);

			//		 ɾ�������еĹ������㣺
			igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);

			//		��ӡ��ǰ����Ľ������ò����127��ʼ�����Խ���
			if (i >=125) 
			{
				char str[256];
				sprintf_s(str, 256, "E:/qslim�۵�������loop%d.obj", i);
				igl::writeOBJ(str, versOut, trisOut);

				sprintf_s(str, 256, "E:/loop%d�۵��ı�.obj", i);
				Eigen::MatrixXi uEdge0 = uEdgesBC.row(uEdgeIdx0);
				objWriteEdgesMat(str, uEdge0, versBC);

				// for debug:
				Eigen::VectorXi uEdgeVec0 = uEdge0.col(0);
				Eigen::VectorXi uEdgeVec1 = uEdge0.col(1);
				std::vector<int> tmpVec0 = vec2Vec(uEdgeVec0);
				std::vector<int> tmpVec1 = vec2Vec(uEdgeVec0);

				std::cout << "loop " << i << " finished." << "�۵��ıߣ�" << uEdge0(0, 0) << ", " << uEdge0(0, 1) << std::endl;
			}
		}
 
		std::cout << "finished." << std::endl;
	}

 
	// ��ӡ���ȶ��ж��׵�����Ԫ�أ�
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


	// ʹ��viewer����Ķ��������𲽲鿴��������еı��۵���
	void test7777() 
	{
		using namespace igl;
		Eigen::MatrixXd vers;
		Eigen:MatrixXi tris;
		igl::readOBJ("E:/����/tmpMesh.obj", vers, tris);

		bool keepWorkingFlag = true;
		bool simplest_decimate_flag = false;
		bool qslim_decimate_flag = ~simplest_decimate_flag;

		// 1. ׼������
		Eigen::MatrixXd versCopy = vers;
		Eigen::MatrixXi trisCopy = tris;
		VectorXi edgeUeInfo;
		MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// ���ȶ��У� 

		int v1 = -1;					// State variables keeping track of edge we just collapsed
		int v2 = -1;
		Eigen::VectorXi timeStamps;
		MatrixXd collapsedVers;
		typedef std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double> Quadric;
		std::vector<Quadric> quadrics;				// ÿ�������Q����

		edge_flaps(trisCopy, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());						// Could reserve with https://stackoverflow.com/a/29236236/148668        
		collapsedVers.resize(uEdges.rows(), versCopy.cols());				// If an edge were collapsed, we'd collapse it to these points:

		// ����Ƿ��з����αߣ�
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
			// ��ȡ�����۵������㷨�ĺ�����
			cost_and_placement = shortest_edge_and_midpoint;
			decimate_trivial_callbacks(pre_collapse, post_collapse);    // ���������۵��е�pre_collapse��post_collapse�����ӡ���ʲô��������
		}
		else
		{
			//  ����ÿ�������Q����
			per_vertex_point_to_plane_quadrics(vers, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, quadrics);

			// ��ȡqslim�㷨�еĺ����ӣ�
			igl::qslim_optimal_collapse_edge_callbacks(uEdges, quadrics, v1, v2, cost_and_placement, pre_collapse, post_collapse);
		}

		decimate_stopping_condition_callback stopping_condition;
 
		int loopCount = 0;
		int num_collapsed;                                // ������ѭ���ļ�����

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
			tris0.conservativeResize(m, tris0.cols());              // �����൱��shrink_to_fit();
			newOldTrisInfo.conservativeResize(m);

			// 3.2. ɾ�������еĹ������㣺
			Eigen::MatrixXd versOut;
			Eigen::MatrixXi trisOut;
			igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);
			char str[256];
			sprintf_s(str, 256, "E:/decimating_output_%d.obj", loopCount);
			igl::writeOBJ(str, versOut, trisOut);
		};
 
		// lambda�����������ݸ�λ
		const auto& reset = [&]()
		{
			tris = trisCopy;
			vers = versCopy;
			edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
			collapsedVers.resize(uEdges.rows(), vers.cols());
			VectorXd costs(uEdges.rows());

			pQueue = {};            // https://stackoverflow.com/questions/2852140/priority-queue-clear-method
			timeStamps = Eigen::VectorXi::Zero(uEdges.rows());

			// r1. ����ÿ���ߵ��۵�costֵ���Լ��۵�֮��Ķ������꣬�������ȶ���pQueue��
			{
				Eigen::VectorXd costs(uEdges.rows());
				igl::parallel_for(uEdges.rows(), \
					[&](const int i)
					{
						double cost = i;
						RowVectorXd edgeCenter(1, 3);           // ȡ�ߵ��е���Ϊ���۵�֮��Ķ��㣻
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


		// lambda����ִ�����񾫼�׼����Ⱦ�����ݣ�
		const auto& pre_draw = [&](igl::opengl::glfw::Viewer& viewer)->bool
		{
			// p1. ÿһ�ζ���ѭ���У�����1%�ıߣ�
			if (viewer.core().is_animating && !pQueue.empty())
			{
				bool FlagCollapsed = false;           // ����ѭ���б������Ƿ�ִ�гɹ���

				// p1.1 ִ�б���������collapse edge
				const int max_iter = std::ceil(0.01 * pQueue.size());
				for (int j = 0; j < max_iter; j++)
				{
					if (!collapse_edge(cost_and_placement, vers, tris, uEdges, \
						edgeUeInfo, UeTrisInfo, UeCornersInfo, pQueue, timeStamps, collapsedVers))              // collapse_edge()����2.1
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

			// p2. �����۵������������Ƭ����������������û�б����ΪIGL_COLLAPSE_EDGE_NULL������Ƭ����
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


		// lambda���������¼���Ӧ��
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


		// 1. ��ʼ������λ����
		reset();

		// 2. �򿪴��ڣ�����ѭ��
		viewer.core().background_color.setConstant(1);
		viewer.core().is_animating = true;                                 // ����
		viewer.core().animation_max_fps = 1.0;						// ָ����󶯻�֡�ʣ�
		viewer.callback_key_down = key_down;
		viewer.callback_pre_draw = pre_draw;
		viewer.launch(keepWorkingFlag);

		// 3. ������ֹ���������
		Eigen::VectorXi newOldTrisInfo;
		Eigen::VectorXi newOldVersInfo;
		
		// 3.1. ɾ�����к��б��ΪIGL_COLLAPSE_EDGE_NULL�ߵ�����Ƭ��
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
		tris0.conservativeResize(m, tris0.cols());              // �����൱��shrink_to_fit();
		newOldTrisInfo.conservativeResize(m);

		// 3.2. ɾ�������еĹ������㣺
		Eigen::MatrixXd versOut;
		Eigen::MatrixXi trisOut;
		igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);
		igl::writeOBJ("E:/decimating_output.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// ����qslim()
	void test77777() 
	{
		using namespace igl;

		using Quadric = std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>;

		Eigen::MatrixXd vers, versOut, vers0, versOri;
		Eigen::MatrixXi tris, trisOut, tris0, trisOri;
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]�Ǿ����������е�i������Ƭ��Ӧ��ԭ���������Ƭ������
		Eigen::VectorXi newOldVersInfo;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
		Eigen::VectorXi timeStamps;
		MatrixXd collapsedVers;
		Eigen::VectorXd costs;
		std::vector<Quadric> quadrics;													// ÿ�������Q����
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// ���ȶ��У�
		Eigen::VectorXi _1, I2;										// ���ɾ����������ڲ�����Ƭʱ�õ���
		bool clean_finish = false;


		// State variables keeping track of edge we just collapsed
		int v1 = -1;								
		int v2 = -1;

		// 0.
		tiktok& tt = tiktok::getInstance();
		tt.start();
		igl::readOBJ("E:/����/tmpMesh.obj", versOri, trisOri);
		tt.endCout("elapsed time of loading mesh is: ");
		igl::writeOBJ("E:/meshIn.obj", versOri, trisOri);

		// 1. 
		unsigned trisCount = trisOri.rows();
		igl::connect_boundary_to_infinity(versOri, trisOri, vers, tris);
		if (!igl::is_edge_manifold(tris))
			return;
		igl::edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);

		// 2. ����ÿ�������Q����
		igl::per_vertex_point_to_plane_quadrics(vers, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, quadrics);

		// 3.  ��ȡqslim�㷨�еĺ����ӣ�
		int trisCountNew = trisCount;
		int tarTrisCount = std::round(0.5 * trisCount);

		// 4. ����ÿ������ߵ�costֵ���Դ�Ϊ���ȼ��������ȶ���
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());
		collapsedVers.resize(uEdges.rows(), vers.cols());
		costs.resize(uEdges.rows());

		igl::parallel_for(uEdges.rows(), [&](const int ueIdx)
			{
				// ������cost_and_placement()�����ݣ�
				double cost = ueIdx;
				RowVectorXd p(1, 3);
 
				// Combined quadric
				Quadric quadric_p;
				quadric_p = quadrics[uEdges(ueIdx, 0)] + quadrics[uEdges(ueIdx, 1)];

				// Quadric: p'Ap + 2b'p + c,  optimal point: Ap = -b, or rather because we have row vectors: pA=-b
				const auto& A = std::get<0>(quadric_p);
				const auto& b = std::get<1>(quadric_p);
				const auto& c = std::get<2>(quadric_p);
				p = -b * A.inverse();
				cost = p.dot(p * A) + 2 * p.dot(b) + c;

				// Force infs and nans to infinity
				if (std::isinf(cost) || cost != cost)
				{
					cost = std::numeric_limits<double>::infinity();
					// Prevent NaNs. Actually NaNs might be useful for debugging.
					p.setConstant(0);
				}
 
				collapsedVers.row(ueIdx) = p;
				costs(ueIdx) = cost;
			},
			10000);

		for (int i = 0; i < uEdges.rows(); i++)
			pQueue.emplace(costs(i), i, 0);
	
		// 5. ���۵���ѭ����
		int uEdgeIdx0;							// ���ȶ����ײ��ıߣ�
		int e1, e2, f1, f2;
		while (true)
		{
			bool collapsed = true;
			std::tuple<double, int, int> edgeTuple;							// ���ȶ����ײ���Ԫ�أ�
			std::vector<int> nbrTrisIdx_src, nbrVersIdx_src;
			std::vector<int>  nbrTrisIdx_des, nbrVersIdx_des;

			// 5.1. ȡ����Ԫ�أ�����Ԫ�س��ӣ�
			while (true)
			{
				// 1.1 ������Ϊ�գ��˳�ѭ����
				if (pQueue.empty())   // no uEdges to collapse
					assert("���۵�����");

				// 1.2ȡ����Ԫ�أ�����Ԫ�س��ӣ�
				edgeTuple = pQueue.top();
				if (std::get<0>(edgeTuple) == std::numeric_limits<double>::infinity())
					assert("���ױߵ�cost�������");

				pQueue.pop();
				uEdgeIdx0 = std::get<1>(edgeTuple);          // ���׵����������;

				// 1.3 Check if matches timestamp
				if (std::get<2>(edgeTuple) == timeStamps(uEdgeIdx0))
					break;

				// 1.4 ������
				assert(std::get<2>(edgeTuple) < timeStamps(uEdgeIdx0) || timeStamps(uEdgeIdx0) == -1);          // must be stale or dead.
			}

			// 5.2. ���㵱ǰ�����˵�1����Ķ��㡢����Ƭ��
			igl::circulation(uEdgeIdx0, true, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_src, nbrTrisIdx_src);
			igl::circulation(uEdgeIdx0, false, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_des, nbrTrisIdx_des);

			//		pre_collapse:
			v1 = uEdges(uEdgeIdx0, 0);
			v2 = uEdges(uEdgeIdx0, 1);

			// 5.3. �۵����׵ıߣ� collapse_edge()����1
			collapsed = collapseSingleEdge(uEdgeIdx0, collapsedVers.row(uEdgeIdx0), nbrVersIdx_src, nbrTrisIdx_src, \
				nbrVersIdx_des, nbrTrisIdx_des, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, e1, e2, f1, f2);

			//		post_collapses:
			if (collapsed)
				quadrics[v1 < v2 ? v1 : v2] = quadrics[v1] + quadrics[v2];


			// �������ظ�����һ�κ�û��֮ǰ������Ƭ�Խ�����֪����ʲôԵ�ʣ�����
#if 0
			if (collapsed)
				quadrics[v1 < v2 ? v1 : v2] = quadrics[v1] + quadrics[v2];
#endif
 
			// 5.4. �۵�����֮�󣬸������timeStamp�� ������رߵ�costֵ
			if (collapsed)
			{
				// 4.1 Erase the two,  other collapsed uEdges by marking their timestamps as -1
				timeStamps(e1) = -1;
				timeStamps(e2) = -1;

				// 4.2
				std::vector<int> nbrTrisIdx;
				nbrTrisIdx.reserve(nbrTrisIdx_src.size() + nbrTrisIdx_des.size());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_src.begin(), nbrTrisIdx_src.end());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_des.begin(), nbrTrisIdx_des.end());
				std::sort(nbrTrisIdx.begin(), nbrTrisIdx.end());
				nbrTrisIdx.erase(std::unique(nbrTrisIdx.begin(), nbrTrisIdx.end()), nbrTrisIdx.end());

				// 4.3 Collect all uEdges that must be updated
				std::vector<int> Ne;
				Ne.reserve(3 * nbrTrisIdx.size());
				for (auto& triIdx : nbrTrisIdx)
				{
					if (tris(triIdx, 0) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 1) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 2) != IGL_COLLAPSE_EDGE_NULL)
					{
						for (int i = 0; i < 3; i++)
						{
							const int ueIdx = edgeUeInfo(i * tris.rows() + triIdx);
							Ne.push_back(ueIdx);
						}
					}
				}

				// Only process edge once
				std::sort(Ne.begin(), Ne.end());
				Ne.erase(std::unique(Ne.begin(), Ne.end()), Ne.end());             // ȥ���ظ�Ԫ�أ�
				for (auto& ueIdx : Ne)
				{
					// ������۵���costֵ�����۵���Ķ������꣺
					double cost;
					RowVectorXd place;
					
					// ������cost_and_placement

					 // Combined quadric
					Quadric quadric_p;
					quadric_p = quadrics[uEdges(ueIdx, 0)] + quadrics[uEdges(ueIdx, 1)];

					// Quadric: place'Ap + 2b'place + c,  optimal point: Ap = -b, or rather because we have row vectors: pA=-b
					const auto& A = std::get<0>(quadric_p);
					const auto& b = std::get<1>(quadric_p);
					const auto& c = std::get<2>(quadric_p);
					place = -b * A.inverse();
					cost = place.dot(place * A) + 2 * place.dot(b) + c;

					// Force infs and nans to infinity
					if (std::isinf(cost) || cost != cost)
					{
						cost = std::numeric_limits<double>::infinity();
						// Prevent NaNs. Actually NaNs might be useful for debugging.
						place.setConstant(0);
					}

					// Increment timestamp
					timeStamps(ueIdx)++;

					// Replace in queue
					pQueue.emplace(cost, ueIdx, timeStamps(ueIdx));
					collapsedVers.row(ueIdx) = place;
				}
			}
			else
				assert("edge collapse failed.");

			// 5.5 ִ��pre_collapse()�������۵��ߣ�ִ��post_collapse()���ٸ�����صıߵ�costֵ��
			if (collapsed)          
			{
				// ��stopping_condition�����ӷ���true����������ֹ�����������۵�ѭ��
				if (f1 < trisCount)
					trisCountNew -= 1;
				if (f2 < trisCount)
					trisCountNew -= 1;

				bool stopConditionFlag = (trisCountNew <= (int)tarTrisCount);
				if (stopConditionFlag)
				{
					clean_finish = true;
					break;
				}
			}
			else			 // ���۵�ʧ�ܣ��˳�ѭ����
				assert("edge collapse failed.");
		}

		// 6. ɾ�����к��б��ΪIGL_COLLAPSE_EDGE_NULL�ߵ�����Ƭ��
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
		tris0.conservativeResize(index, tris0.cols());              // �����൱��shrink_to_fit();
		newOldTrisInfo.conservativeResize(index);

		// 7. ɾ�������еĹ������㣺
		igl::remove_unreferenced(vers, tris0, versOut, trisOut, _1, newOldVersInfo);

		//		for debug:
		igl::writeOBJ("E:/ɾ���ڲ�����Ƭǰ.obj", versOut, trisOut);

		// 8. ������ɾ���ڲ�����Ƭ
		const Eigen::Array<bool, Eigen::Dynamic, 1> keep = (newOldTrisInfo.array() < trisCount);
		igl::slice_mask(Eigen::MatrixXi(trisOut), keep, 1, trisOut);
		igl::slice_mask(Eigen::VectorXi(newOldTrisInfo), keep, 1, newOldTrisInfo);
		igl::remove_unreferenced(Eigen::MatrixXd(versOut), Eigen::MatrixXi(trisOut), versOut, trisOut, _1, I2);
		igl::slice(Eigen::VectorXi(newOldVersInfo), I2, 1, newOldVersInfo);

		igl::writeOBJ("E:/qslimOutput.obj", versOut, trisOut);
		std::cout << "finished." << std::endl;
	}


	// �󾫼��������ԭʼ����Ľ���������QEM��paper����
	double calcSimpApproxError(const Eigen::MatrixXd& versSimp, const Eigen::MatrixXi& trisSimp, \
		const Eigen::MatrixXd& versOri, const Eigen::MatrixXi& trisOri)
	{
		// ע�ⲻ�ܳ���̫��ľ��󣬷�������bad alloc

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

		// �㵽ƽ����룺 dis == p.dot(v)��ע�����з��ŵľ��룻 p == (a,b,c,d), v = (x0, y0, z0, 1), (a,b,c)��ƽ��Ĺ�һ����������
		Eigen::MatrixXd versExt0{ Eigen::MatrixXd::Ones(versCount0, 4) };
		Eigen::MatrixXd versExt1{ Eigen::MatrixXd::Ones(versCount1, 4) };
		versExt0.leftCols(3) = versMat0.array().cast<double>();
		versExt1.leftCols(3) = versMat1.array().cast<double>();

		// for debug
		Eigen::Vector4d verColVec = versExt0.row(3).transpose();
		// ����ΪverIdx��simp�����еĶ��㣬��ori��������������Ƭƽ��ķ��ž��룺
		Eigen::VectorXd signedDises = planeCoeff0 * verColVec;
		Eigen::VectorXd sqrDises = signedDises.array() * signedDises.array();

		// for new method:
		PARALLEL_FOR(0, versCount1, [&](int verIdx)
			{
				std::lock_guard<mutex> guard(g_mutex);
				Eigen::Vector4d verColVec = versExt1.row(verIdx).transpose();
				// ����ΪverIdx��simp�����еĶ��㣬��ori��������������Ƭƽ��ķ��ž��룺
				Eigen::VectorXd signedDises = planeCoeff0 * verColVec;
				Eigen::VectorXd sqrDises = signedDises.array() * signedDises.array();
				squaredDis1.push_back(sqrDises.minCoeff());
			});

		PARALLEL_FOR(0, versCount0, [&](int verIdx)
			{
				std::lock_guard<mutex> guard(g_mutex);
				Eigen::Vector4d verColVec = versExt0.row(verIdx).transpose();
				// ����ΪverIdx��simp�����еĶ��㣬��ori��������������Ƭƽ��ķ��ž��룺
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


	// ������������Ľ�����
	void test777777() 
	{
		Eigen::MatrixXd vers0, vers1;
		Eigen::MatrixXi tris0, tris1;
		objReadMeshMat(vers0, tris0, "E:/����/jawMeshDense4.obj");
		objReadMeshMat(vers1, tris1, "E:/jawMeshSimplified4.obj");

		double appErr = calcSimpApproxError(vers1, tris1, vers0, tris0);

		std::cout << "appErr == " << appErr << std::endl;
		std::cout << "finished." << std::endl;
	}


	// ��������
	void test8() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut, trisCopy;
		Eigen::VectorXi selectedIdxes, oldNewIdxInfo;

		igl::readOBJ("E:/arrangeResultMesh.obj", vers, tris);

		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		// ��ӡ��ʼ��Ϣ��
		Eigen::MatrixXi bdrys, bdryTris;
		std::vector<int> bdryTriIdxes;
		bdryEdges(bdrys, bdryTriIdxes, tris);
		subFromIdxVec(bdryTris, tris, bdryTriIdxes);
		igl::writeOBJ("E:/bdryTris.obj", vers, bdryTris);
		std::cout << "versCount == " << versCount << std::endl;
		std::cout << "trisCount == " << trisCount << std::endl;
		std::cout << "bdrysCount == " << bdrys.rows() << std::endl;

		// for debug: ������ӡ��Ե�ߵĶ��㣺
		char str[512];
		for (unsigned i = 0; i<bdrys.rows(); ++i) 
		{
			Eigen::MatrixXd tmpVers;
			subFromIdxVec(tmpVers, vers, std::vector<int>{bdrys(i, 0), bdrys(i, 1)});
			sprintf_s(str, "E:/bdryVers%d.obj", i);
			objWriteVerticesMat(str, tmpVers);
		}
		
		// 1. ȥ��duplicated vertices:
		igl::remove_duplicate_vertices(vers, 0.001, versOut, selectedIdxes, oldNewIdxInfo);
		objWriteVerticesMat("E:/versCleaned.obj", versOut);

		trisCopy = tris;
		int* ptr = trisCopy.data();
		for (unsigned i = 0; i < 3*trisCount; ++i) 
		{
			int oldIdx = *ptr;
			*ptr = oldNewIdxInfo(oldIdx);
			ptr++;
		}

		// 2. ȥ���Ƿ�����Ƭ��
		std::vector<unsigned> sickTriIdxes;
		checkSickTris(sickTriIdxes, trisCopy);

		std::vector<int> tmpVec;
		tmpVec.resize(trisCount);
		for (unsigned i = 0; i < trisCount; ++i)
			tmpVec[i] = i;

		for (const auto& index : sickTriIdxes)
			tmpVec[index] = -1;

		std::vector<int> selectedTriIdx;
		selectedTriIdx.reserve(trisCount);
		for (const auto& index : tmpVec)
			if (index >= 0)
				selectedTriIdx.push_back(index);

		subFromIdxVec(trisOut, trisCopy, selectedTriIdx);
		bdrys.resize(0, 0);
		bdryEdges(bdrys, trisOut);

		igl::writeOBJ("E:/meshOut.obj", versOut, trisOut);

		// ��ӡ������Ϣ��
		versCount = versOut.rows();
		trisCount = trisOut.rows();
		std::cout << "final versCount == " << versCount << std::endl;
		std::cout << "final trisCount == " << trisCount << std::endl;
		std::cout << "final bdrysCount == " << bdrys.rows() << std::endl;

		std::cout << "finished." << std::endl;
	}
}


// libigl�е�΢�ּ������
namespace IGL_DIF_GEO 
{

	// ���������LB����
	void test0() 
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		Eigen::SparseMatrix<double> L, M;

		igl::readOBJ("E:/����/tooth.obj", vers, tris);
		igl::cotmatrix(vers, tris, L);
		igl::massmatrix(vers, tris, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);

		dispSpMat(M, 0, M.rows() - 1, 10);

		std::cout << "finished." << std::endl;
	}

}


// ͼ�㷨
namespace IGL_GRAPH 
{
	// ͼ���ݽṹ��ת����
	void test0() 
	{
		MatrixXd vers;
		MatrixXi tris;
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
		MatrixXd vers;
		MatrixXi tris;
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
		MatrixXd vers;
		MatrixXi tris;
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


// �ռ仮��
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


// IGLʵ�ֵĻ��������������㷨��
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

		//	 nonManifoldEdges()������������ΰ�ߣ�
		Eigen::MatrixXi  nmnEdges;
		nonManifoldEdges(tris, nmnEdges);
		dispMatBlock(nmnEdges, 0, 10, 0, 1);
		objWriteEdgesMat("E:/nmnEdges.obj", nmnEdges, vers);

		// 4. �ҳ����з����α���������Ƭ��
		std::vector<int> trisIdxNM;
		for (auto& pair : etNMNmap)
		{
			for (auto& index : pair.second)
				trisIdxNM.push_back(index);
		}
		Eigen::MatrixXi trisNM(Eigen::MatrixXi::Zero(trisIdxNM.size(), 3));
		for (int i = 0; i < trisIdxNM.size(); ++i)
			trisNM.row(i) = tris.row(trisIdxNM[i]);

		// for debug: ��ӡtrisNM:
		objWriteMeshMat("E:/trisNM.obj", vers, trisNM);
		std::cout << "finished." << std::endl;

		// 5. ȷ����������Ƭ������������Ƭ���ڽӹ�ϵ��

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

		Eigen::MatrixXi ttAdj_nmEdge;
		std::vector<ttTuple> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		bool retFlag = buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);

		// ������������αߵ�����Ƭ��
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
		Eigen::MatrixXi ttAdj_nmEdge;
		std::vector<ttTuple> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);
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
				// wf0. ����ǰ��Ϊ���αߣ�����ttAdj_nmEdge��Ѱ������棻
				int nbrTriIdx = ttAdj_nmEdge(currentTriIdx, i);
				
				// wf1. ����ǰ��Ϊ�����αߣ����ݲ���������ͬ����ѡȡ��ͬ�Ķ��棻
				if (-1 == nbrTriIdx)
				{
					std::vector<std::vector<int>*> vecPtr(3), tmpPtrVec(3);
					vecPtr[0] = &std::get<0>(ttAdj_nmnEdge[currentTriIdx]);
					vecPtr[1] = &std::get<1>(ttAdj_nmnEdge[currentTriIdx]);
					vecPtr[2] = &std::get<2>(ttAdj_nmnEdge[currentTriIdx]);
					tmpPtrVec[0] = &std::get<0>(ttAdj_nmnOppEdge[currentTriIdx]);
					tmpPtrVec[1] = &std::get<1>(ttAdj_nmnOppEdge[currentTriIdx]);
					tmpPtrVec[2] = &std::get<2>(ttAdj_nmnOppEdge[currentTriIdx]);

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
		igl::readOBJ("E:/����/jawMeshArranged.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);
		igl::writeOBJ("E:/triNMN0.obj", vers, Eigen::MatrixXi{ tris.row(4431) });

		unsigned versCount = vers.rows();
		unsigned trisCount = tris.rows();

		Eigen::MatrixXi ttAdj_nmEdge;
		std::vector<ttTuple> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		bool retFlag = buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);

		// ������������αߵ�����Ƭ��
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

		Eigen::MatrixXi mnTris;				// �����������αߵ�����Ƭ��
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


	//		����AABB������դ��
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
		{
			if (i != maxCompIdx)
				gridCounts(i) = std::ceil(largestCount0 * (box.diagonal()(i)) / maxComp);
		}
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
				std::int64_t edgeCode = encodeEdge(vaIdx, vbIdx);
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


	// marchingCubes()��������������Ҫ�����ȥ���ظ����㣬̫�̵ıߡ�
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
 

	// test4������������դ��marchingCubes�㷨��
	void test4()
	{
		std::vector<int> stepCounts(3);				// xyz����ά����դ����
		Eigen::RowVector3d gridsOri;					// դ��ԭ�㣺
		Eigen::VectorXd SDF;
		Eigen::RowVector3i gridCounts;			// xyz����ά���ϵ�դ��������
		Eigen::MatrixXd gridCenters;				//	 ����դ���е�����ľ���ÿ�ж���һ���е����ꣻ�洢���ȼ���x, y, z
		Eigen::MatrixXd	boxVers;				// դ���Ӧ�İ�Χ�еĶ��㣻
		Eigen::MatrixXi boxTris;

		// 0. ����SDFGen.exe���ɵ�.sdf���볡�����ļ���
		const char* sdfFilePath = "E:/jawMeshUnionRepair1.sdf";
		double SDFstep = IGL_BASIC::parseSDF(stepCounts, gridsOri, SDF, sdfFilePath);

		// 1. ����դ��
		Eigen::RowVector3d minp = gridsOri;
		Eigen::RowVector3d maxp = gridsOri + SDFstep * Eigen::RowVector3d(stepCounts[0], stepCounts[1], stepCounts[2]);
		Eigen::AlignedBox<double, 3> box(minp, maxp);		// դ���Ӧ�İ�Χ�У�
		genGrids(box, std::max({ stepCounts[0], stepCounts[1], stepCounts[2] }), 0, gridCenters, gridCounts);
		genAABBmesh(box, boxVers, boxTris);
		objWriteMeshMat("E:/AABB.obj", boxVers, boxTris);

		// դ�����ݵķֲ���
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

		// ��ȡդ����SDFֵС��0�Ķ��㣺
		Eigen::MatrixXd tmpVers(SDF.rows(), 3);
		int index = 0;
		for (int i = 0; i<SDF.rows(); ++i) 
			if (SDF(i) <= 0)
				tmpVers.row(index++) = gridCenters0.row(i);
		tmpVers.conservativeResize(index, 3);
		igl::writeOBJ("E:/tmpVers.obj", tmpVers, Eigen::MatrixXi{});

		// 2. marching cubes�㷨�����������棺
		tiktok& tt = tiktok::getInstance();
		MatrixXd versResult_SDF, versResults_signs, versResult_origin;
		MatrixXi trisResult_SDF, trisResults_signs, trisResult_origin;
		double selectedSDF = -1.5;
		tt.start();
		marchingCubes(versResult_SDF, trisResult_SDF, SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), selectedSDF);
		tt.endCout("Elapsed time of igl::marching_cubes() is ");
		igl::writeOBJ("E:/shrinkedMesh.obj", versResult_SDF, trisResult_SDF);

		// ԭʼ��marching cubes
		versResult_SDF.resize(0, 0);
		trisResult_SDF.resize(0, 0);
		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), selectedSDF, versResult_SDF, trisResult_SDF);
		igl::writeOBJ("E:/shrinkedMeshOri.obj", versResult_SDF, trisResult_SDF);
 

		std::cout << "finished." << std::endl;
	}


	// ������opological_hole_fill()������ǰ������
	void test5()
	{
		Eigen::MatrixXd vers;
		Eigen::MatrixXi tris;
		igl::readOBJ("E:/����/holeMesh.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		// ȷ�����ıߡ���ֻ����һ������Ƭ�ıߣ�
		const unsigned trisCount = tris.rows();
		const unsigned edgesCount = 3 * trisCount;
		const unsigned versCount = tris.maxCoeff() + 1;

		Eigen::MatrixXi edges;
		Eigen::SparseMatrix<int> adjSM_eCount;				// ������ڽӾ���Ȩ��Ϊ��������ظ��Ĵ�����
		Eigen::SparseMatrix<int> adjSM_ueCount;				// ������ڽӾ���Ȩ��Ϊ�ñ��ظ��Ĵ�����

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

		// �����ڽӾ���
		std::vector<Eigen::Triplet<int>> smElems;
		smElems.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});

		adjSM_eCount.resize(versCount, versCount);
		adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());		// ������ڽӾ���Ȩ��Ϊ��������ظ��Ĵ�����
		Eigen::SparseMatrix<int> tmpSp = adjSM_eCount.transpose();
		adjSM_ueCount = adjSM_eCount + tmpSp;		// ������ڽӾ���Ȩ��Ϊ�ñ��ظ��Ĵ�����

		// Ѱ�ұ�Ե����ߣ�
		std::unordered_set<std::int64_t> ueSet;
		for (unsigned i = 0; i < adjSM_ueCount.outerSize(); ++i)
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_ueCount, i); iter; ++iter)
				if (1 == iter.value())
					ueSet.insert(encodeEdge(iter.row(), iter.col()));

		// Ѱ�ұ�Ե��ߣ�
		std::unordered_set<std::int64_t> bdryEdgeSet;
		for (const auto& eCode : ueSet)
		{
			std::pair<int, int> retPair = decodeEdge(eCode);
			int vaIdx = retPair.first;
			int vbIdx = retPair.second;
			std::int64_t eCodeOpp = encodeEdge(vbIdx, vaIdx);

			if (adjSM_eCount.coeff(vaIdx, vbIdx) > 0)
				bdryEdgeSet.insert(eCode);
			if (adjSM_eCount.coeff(vbIdx, vaIdx) > 0)
				bdryEdgeSet.insert(eCodeOpp);
		}

		Eigen::MatrixXi bdryEdges(bdryEdgeSet.size(), 2);
		unsigned index = 0;
		for (const auto& eCode : bdryEdgeSet)
		{
			auto retPair = decodeEdge(eCode);
			bdryEdges(index, 0) = retPair.first;
			bdryEdges(index, 1) = retPair.second;
			index++;
		}
		objWriteEdgesMat("E:/bdryEdges.obj", bdryEdges, vers);

		// ���㶴�����򶥵�����
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

		// ��Ҫ�ֶ�Ϊԭ������Ӷ������ĵ㣬����ò�ƴ������������λ����
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
		adjMatrix(tris, adjSM_eCount, adjSM_eIdx);

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

		bool retFlag = simplyConnectedLargest(vers, tris, versOut, trisOut);
		if (!retFlag)
			std::cout << "function failed!!!" << std::endl;

		objWriteMeshMat("E:/meshOut.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}
}
