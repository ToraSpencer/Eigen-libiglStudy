#include "igl_study.h"

#define DATA_PATH "./data/"



// libigl��������
namespace IGL_BASIC
{
	Eigen::MatrixXd vers, newVers, normals;
	Eigen::MatrixXi tris;
	Eigen::SparseMatrix<double> L;
	igl::opengl::glfw::Viewer viewer;		// libigl�еĻ���glfw����ʾ���ڣ�

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
		// unique_edge_map()������������İ�ߡ�����ߡ������Ӧ��ϵ��
		bool retFlag = igl::readOBJ("E:/����/tooth.obj", vers, tris);
		Eigen::MatrixXi edges, uEdges;
		Eigen::VectorXi edgeUeInfo;
		std::vector<std::vector<int>> UeEdgeInfo;
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, UeEdgeInfo);
		//std::cout << "tris:" << std::endl;
		//dispMat(tris);
		//std::cout << "edges:" << std::endl;
		//dispMat(edges);				// ���
		//std::cout << "uEdges:" << std::endl;
		//dispMat(uEdges);			// ����ߣ�����undirected������unique�ģ�
		//std::cout << "edgeUeInfo:" << std::endl;
		//dispMat(edgeUeInfo);		// �����������ж�Ӧ��������edgeUeInfo(i)������Ϊi�İ�߶�Ӧ������ߵ�������

		Eigen::VectorXi uEC, uEE;
		igl::unique_edge_map(tris, edges, uEdges, edgeUeInfo, uEC, uEE);
		//std::cout << "uEC:" << std::endl;
		//dispMat(uEC);		 
		//std::cout << "uEE:" << std::endl;
		//dispMat(uEE);

		// edge_flaps()������������߹���������Ƭ��
		Eigen::MatrixXi UeTrisInfo, UeCornersInfo;
		igl::edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
		//std::cout << "UeTrisInfo:" << std::endl;
		//dispMat(UeTrisInfo);								// ��i�е�Ԫ�ء�������Ϊi������߹���������Ƭ��������
		//std::cout << "UeCornersInfo:" << std::endl;
		//dispMat(UeCornersInfo);

		objWriteEdgesMat("E:/edge0.obj", Eigen::RowVector2i(uEdges.row(0)), vers);
		Eigen::RowVector3i tri00 = tris.row(UeTrisInfo(0, 0));
		Eigen::RowVector3i tri01 = tris.row(UeTrisInfo(0, 1));
		Eigen::MatrixXd oppVers(2, 3);
		oppVers.row(0) = vers.row(tri00(UeCornersInfo(0, 0)));
		oppVers.row(1) = vers.row(tri01(UeCornersInfo(0, 1)));
		objWriteVerticesMat("E:/oppVers.obj", oppVers);

		// circulation()����Ѱ������ߵĶ˵����������Ƭ������1���򶥵㣻
		std::vector<int> nbrVersIdx, nbrTrisIdx;
		igl:: circulation(0, true, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx, nbrTrisIdx);
		objWriteVerticesMat("E:/edgeHead0.obj", Eigen::MatrixXd{ vers.row(uEdges(0, 0)) });

		Eigen::MatrixXd vers0;
		subFromIdxVec(vers0, vers, nbrVersIdx);
		objWriteVerticesMat("E:/vers0.obj", vers0);
		Eigen::MatrixXi tris0;
		subFromIdxVec(tris0, tris, nbrTrisIdx);
		igl::writeOBJ("E:/tris0.obj", vers, tris0);

		// igl::doubleArea()����������������ÿ������Ƭ�����������
		Eigen::VectorXd dbArea;
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
	
		//vers.resize(0, 0);
		//tris.resize(0, 0);
		//igl::readOBJ("E:/fatTeeth1_Ԥ�����.obj", vers, tris);
		//igl::writeOFF("E:/fatTeeth1_Ԥ�����.off", vers, tris);

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


	// ������ž��볡��marching cubes��ȡ��ֵ������
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

		// �ڶ��У�դ��ԭ��
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


	void test55() 
	{
		// 0. ����SDFGen.exe���ɵ�.sdf���볡�����ļ���
		std::vector<int> stepCounts(3);				// xyz����ά����դ����
		Eigen::RowVector3d gridsOri;					// դ��ԭ�㣺
		Eigen::VectorXd SDF;
		const char* sdfFilePath = "E:/inputMesh.sdf";
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

		std::string fileName = "jawMeshDenseRepaired";

		// igl::readOBJ("E:/����/jawMeshDense.obj", vers, tris);
		igl::readOBJ((std::string{ "E:/����/" } + fileName + std::string{".obj"}).c_str(), vers, tris);
		unsigned trisCount = tris.rows();
		// unsigned tarTrisCount = std::round(trisCount * 0.1);			// ���������������Ƭ����
		unsigned tarTrisCount = 10000;
		igl::writeOBJ("E:/meshIn.obj", vers, tris);

		// igl::decimate()�������۵��㷨�������� ���񲻿����з����αߣ�				
		tt.start();
		std::cout << "succeeded? " << igl::decimate(vers, tris, tarTrisCount, versOut, trisOut, \
			newOldTrisInfo, newOldVersInfo) << std::endl;						// ����1.1
		tt.endCout("Elapsed time of mesh simplification by igl::decimate() is ");
		igl::writeOBJ((std::string{ "E:/" } + fileName + std::string{ "_simplified_edgeCollapse.obj" }).c_str(), versOut, trisOut);

		std::vector<int> newOldTrisInfoVec = vec2Vec(newOldTrisInfo);
		Eigen::MatrixXd vers1;
		subFromIdxVec(vers1, vers, newOldVersInfo);

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


	// collapse_small_triangles()
	void test777() 
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		tiktok& tt = tiktok::getInstance();
		igl::readOBJ("E:/����/tooth.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);
	
		Eigen::VectorXd dbArea;
		igl::doublearea(vers, tris, dbArea);
		std::vector<double> areaVec = vec2Vec(dbArea);
		std::sort(areaVec.begin(), areaVec.end());
		std::for_each(areaVec.begin(), areaVec.end(), [&](double& num)
			{
				num /= 2.0;
			});

		// �۵�1/2����������Ƭ��
		unsigned trisCount = tris.rows();
		double areaThreshold = areaVec[std::ceil(trisCount / 2.0)];
		double boxArrow = igl::bounding_box_diagonal(vers);            // ��Χ�еĶԽ��߳��ȣ�

		//			const double min_dblarea = 2.0 * eps * boxArrow * boxArrow;
		double eps = areaThreshold / (boxArrow * boxArrow);

		for (unsigned i = 0; i < trisCount - 1; ++i)
		{
			if (areaVec[i] < eps && areaVec[i + 1] >= eps)
			{
				std::cout << "pos == " << i + 1 << std::endl;
				break;
			}
		}

		Eigen::MatrixXi trisNew;
		tt.start();
		igl::collapse_small_triangles(vers, tris, eps, trisNew);
		tt.endCout("Elapsed time of mesh simplification by igl::collapse_small_triangles()");

		igl::writeOBJ("E:/meshSimplified_trisCollapse.obj", vers, trisNew);
		std::cout << "trisNew.rows() == " << trisNew.rows() << std::endl;

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


	// libigl�е����񲼶���������������Ҫ����CGAL�⣬��ǰδ��ɣ�
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
		igl::adjacency_matrix(tris, adjSM);		// ����ߴ�����Ԫ��Ϊ1������Ϊ0��

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
 
		std::vector<std::vector<int>> adjList;
		Eigen::SparseMatrix<int> adjSM;
		igl::adjacency_list(tris, adjList);
		igl::adjacency_matrix(tris, adjSM);

		// dfs:
		Eigen::VectorXi disCoveredIdx, bfsTreeVec, dfsTreeVec, closedIdx;
		size_t startIdx = 165;			// �ӽ�Բ�ĵĶ��㣻
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
		std::list<int> discoveredVersIdx;				// �ѱ����ʵĶ��㣺
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
		igl::readOBJ("E:/����/meshArranged_hole.obj", vers, tris);
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

#if __cplusplus < 202002L
		// 3. ȷ�������αߣ�
		std::unordered_set<std::pair<int, int>, edgeHash, edgeComparator> edgesNM;
#else
		// lambda����std::pair<int, int>��ʾ�ıߵ��Զ����ϣ������
		auto edgeHashLamb = [](const std::pair<int, int>& edge)->std::size_t
		{
			return (std::hash<int>()(edge.first) + std::hash<int>()(edge.second));
		};

		// lambda����std::pair<int, int>��ʾ�ıߵĵȼ۱Ƚ�����
		auto edgeComLamb = [](const std::pair<int, int>& edge1, const std::pair<int, int>& edge2)->bool
		{
			return (edge1.first == edge2.first && edge1.second == edge2.second);
		};

		// 3. ȷ�������αߣ�
		std::unordered_set<std::pair<int, int>, decltype(edgeHashLamb), decltype(edgeComLamb)> edgesNM;				// ��ҪC++20;�����αߣ�ʹ��unordered_setȥ�أ�
#endif
		 for (unsigned i = 0; i < adjSM_eCount.outerSize(); ++i)
		{
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_eCount, i); iter; ++iter)	// ��i�е����ڵ�������
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
		 //		�����������-�������ļ�ֵ�ԣ�һ�������α߶�Ӧ�Ŷ����������
		 std::unordered_map<std::pair<int, int>, std::vector<int>, edgeHash, edgeComparator> edgesNMmap;

		 //		�����������-�����ڵ�����Ƭ�������ļ�ֵ�ԣ�һ�������α߶�Ӧ�Ŷ������Ƭ������
		 std::unordered_map<std::pair<int, int>, std::vector<int>, edgeHash, edgeComparator> etNMmap;
#else
		 //		�����������-�������ļ�ֵ�ԣ�һ�������α߶�Ӧ�Ŷ����������
		 std::unordered_map<std::pair<int, int>, std::vector<int>, decltype(edgeHashLamb), decltype(edgeComLamb)> edgesNMmap;
		 //		�����������-�����ڵ�����Ƭ�������ļ�ֵ�ԣ�һ�������α߶�Ӧ�Ŷ������Ƭ������
		 std::unordered_map<std::pair<int, int>, std::vector<int>, decltype(edgeHashLamb), decltype(edgeComLamb)> etNMmap;
#endif

		for (auto& nmEdge : edgesNM)
		{
			for (int i = 0; i < edgesCount; ++i)
			{
				if (edges(i, 0) == nmEdge.first && edges(i, 1) == nmEdge.second)
				{
					auto retPair = edgesNMmap.insert({nmEdge, std::vector<int>{i} });
					if (!retPair.second)			// ������ʧ�ܣ���˵�����д˼���
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

		// 4. �ҳ����з����α���������Ƭ��
		std::vector<int> trisIdxNM;
		for (auto& pair : etNMmap)
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
		igl::readOBJ("E:/����/meshArranged.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		Eigen::MatrixXi ttAdj_nmEdge;
		std::vector<ttTuple> ttAdj_nmnEdge, ttAdj_nmnOppEdge;
		bool retFlag = buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);

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
	void test33()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut, tris0;
		igl::readOBJ("E:/����/jawMeshArranged.obj", vers, tris);
		igl::writeOBJ("E:/meshInput.obj", vers, tris);

		Eigen::MatrixXd triNorms;
		igl::per_face_normals(vers, tris, triNorms);

		// 1. ȡz�������Ķ������ڵ�����ƬΪ������ɢ��ʼ����Ƭ��ȷ�������ⲿ����Ƭ
		Eigen::Index verIdx0 = 0;
		Eigen::MatrixXi trisTrans = tris.transpose();
		vers.col(2).maxCoeff(&verIdx0);
		int triIdx0 = 0;
		for (int i = 0; i < trisTrans.size(); ++i)
		{
			if (verIdx0 == *(trisTrans.data() + i))
			{
				triIdx0 = i / 3;
				break;
			}
		}
		objWriteMeshMat("E:/tri0.obj", vers, tris.row(triIdx0));

		// 2. �����ػ��ֺ�����������Ƭ�ڽӹ�ϵ��
		Eigen::MatrixXi ttAdj_nmEdge;
		std::vector<ttTuple> ttAdj_nmnEdge;
		std::vector<ttTuple> ttAdj_nmnOppEdge;
		buildAdjacency(tris, ttAdj_nmEdge, ttAdj_nmnEdge, ttAdj_nmnOppEdge);
#if 1
		{
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
		}
#endif

		// 3. ��������ѭ����
		std::deque<int> triIdxQueue;
		std::vector<bool> visited(tris.rows(), false);
		int currentTriIdx = -1;
		visited[triIdx0] = 1;
		triIdxQueue.push_back(triIdx0);
		while (!triIdxQueue.empty()) 
		{
			currentTriIdx = triIdxQueue.front();
			triIdxQueue.pop_front();

			// �Ե�ǰ����Ƭ�����ߵı�����
			for (int i = 0; i < 3; ++i)
			{
				// wf0: ����ǰ�������αߣ����ttAdj_nmEdge�л�ȡ�����������
				int nbrTriIdx = ttAdj_nmEdge(currentTriIdx, i);

				// wf1. ����ǰ��Ϊ�����α�;
				if (-1 == nbrTriIdx)
				{
					// �����α����������棬�Ƚ�����������ķ���͵�ǰ����Ƭ����ļнǣ�ȡ�нǸ�����Ǹ����棺
					std::vector<std::vector<int>*> tmpPtrVec(3);
					tmpPtrVec[0] = &std::get<0>(ttAdj_nmnOppEdge[currentTriIdx]);
					tmpPtrVec[1] = &std::get<1>(ttAdj_nmnOppEdge[currentTriIdx]);
					tmpPtrVec[2] = &std::get<2>(ttAdj_nmnOppEdge[currentTriIdx]);
					const std::vector<int>& relaOppTrisIdx = *tmpPtrVec[i];			// ��ǰ�����αߵ����ж�������Ƭ�������Ļ�Ӧ��Ϊ��������
					Eigen::RowVector3d currentNorm = triNorms.row(currentTriIdx);

					double minCosValue = 2.0;
					for (const auto& triIdx : relaOppTrisIdx)
					{
						Eigen::RowVector3d oppNorm = triNorms.row(triIdx);
						double cosValue = currentNorm.dot(oppNorm);				// ���ԽС���н�Խ��
						if (cosValue < minCosValue)
						{
							nbrTriIdx = triIdx;
							minCosValue = cosValue;
						}
					}
				}

				// wf2. ��Ƿ���֮�������Ƭ����β�����µ�����Ƭ
				if (!visited[nbrTriIdx])
				{
					visited[nbrTriIdx] = true;
					triIdxQueue.push_back(nbrTriIdx);
				}
			}
		}

		// 4. ��ȡ��������ѭ���з��ʵ�����Ƭ�����������
		std::vector<int> selectedTrisIdx;
		selectedTrisIdx.reserve(tris.rows());
		for (int i = 0; i < tris.rows(); ++i)
		{
			if (visited[i] > 0)
				selectedTrisIdx.push_back(i);
		}
		selectedTrisIdx.shrink_to_fit();
		subFromIdxVec(tris0, tris, selectedTrisIdx);
		
		// 5. ȥ���������㣺
		Eigen::VectorXi I, newOldVersInfo;
		igl::remove_unreferenced(vers, tris0, versOut, trisOut, I, newOldVersInfo);
		igl::writeOBJ("E:/meshOut.obj", versOut, trisOut);

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


	bool marchingCubes(const Eigen::VectorXd& SDF, const Eigen::MatrixXd& gridCenters, const Eigen::RowVector3i& gridCounts, \
				const double selectedSDF, Eigen::MatrixXd& versResult_SDF, Eigen::MatrixXi& trisResult_SDF)
	{


		return true;
	}


	void test4()
	{
		// 0. ����SDFGen.exe���ɵ�.sdf���볡�����ļ���
		std::vector<int> stepCounts(3);				// xyz����ά����դ����
		Eigen::RowVector3d gridsOri;					// դ��ԭ�㣺
		Eigen::VectorXd SDF;
		const char* sdfFilePath = "E:/inputMesh.sdf";
		double SDFstep = IGL_BASIC::parseSDF(stepCounts, gridsOri, SDF, sdfFilePath);

		// 1. ����դ��
		Eigen::RowVector3i gridCounts;
		Eigen::MatrixXd gridCenters;
		Eigen::RowVector3d minp = gridsOri - SDFstep * Eigen::RowVector3d(stepCounts[0] / 2.0, stepCounts[1] / 2.0, stepCounts[2] / 2.0);
		Eigen::RowVector3d maxp = gridsOri + SDFstep * Eigen::RowVector3d(stepCounts[0] / 2.0, stepCounts[1] / 2.0, stepCounts[2] / 2.0);
		Eigen::AlignedBox<double, 3> box(minp, maxp);
		genGrids(box, std::max({ stepCounts[0], stepCounts[1], stepCounts[2] }), 0, gridCenters, gridCounts);

		// gridCenters������դ���е�����ľ���ÿ�ж���һ���е����ꣻ�洢���ȼ���x, y, z
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
		dispMatBlock(gridCenters, 0, 144, 0, 2);

		// 2. marching cubes�㷨�����������棺
		tiktok& tt = tiktok::getInstance();
		MatrixXd versResult_SDF, versResults_signs;
		MatrixXi trisResult_SDF, trisResults_signs;
		double selectedSDF = -1.;
		tt.start();
		igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), selectedSDF, versResult_SDF, trisResult_SDF);
		tt.endCout("Elapsed time of igl::marching_cubes() is ");

		igl::writeOBJ("E:/shrinkedMesh.obj", versResult_SDF, trisResult_SDF);

	}


	// ������opological_hole_fill()
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

		// �����ڽӾ���
		std::vector<Eigen::Triplet<int>> smElems;
		smElems.reserve(edgesCount);
		for (int i = 0; i < edgesCount; ++i)
			smElems.push_back(Eigen::Triplet<int>{edges(i, 0), edges(i, 1), 1});
		Eigen::SparseMatrix<int> adjSM_eCount;
		adjSM_eCount.resize(versCount, versCount);
		adjSM_eCount.setFromTriplets(smElems.begin(), smElems.end());		// ������ڽӾ���Ȩ��Ϊ��������ظ��Ĵ�����
		Eigen::SparseMatrix<int> tmpSp = adjSM_eCount.transpose();
		Eigen::SparseMatrix<int> adjSM_ueCount = adjSM_eCount + tmpSp;		// ������ڽӾ���Ȩ��Ϊ�ñ��ظ��Ĵ�����

		// Ѱ�ұ�Ե����ߣ�
		std::unordered_set<std::pair<int, int>, edgeHash> ueSet;
		for (unsigned i = 0; i < adjSM_ueCount.outerSize(); ++i)
			for (auto iter = Eigen::SparseMatrix<int>::InnerIterator(adjSM_ueCount, i); iter; ++iter)
				if (1 == iter.value())
					ueSet.insert(std::make_pair(iter.row(), iter.col()));

		// Ѱ�ұ�Ե��ߣ�
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
}
