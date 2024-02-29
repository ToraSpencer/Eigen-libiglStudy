#include "test_sparse_mat.h"


////////////////////////////////////////////////////////////////////////////////////////////// DEBUG �ӿ�
namespace MY_DEBUG
{
	static std::string g_debugPath = "E:/";

	// lambda������ӡstd::cout֧�ֵ����ͱ�����
	template <typename T>
	static auto disp = [](const T& arg)
	{
		std::cout << arg << ", ";
	};

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


	template <typename Derived>
	static void dispData(const Eigen::MatrixBase<Derived>& m)
	{
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
	static void debugWriteVers(const char* name, const Eigen::MatrixBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat(path, vers);
	}


	template<typename DerivedV>
	static void debugWriteVers2D(const char* name, const Eigen::MatrixBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat2D(path, vers);
	}


	template<typename DerivedV>
	static void debugWriteMesh(const char* name, \
		const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteMeshMat(path, vers, tris);
	}


	template<typename DerivedV>
	static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, \
		const Eigen::MatrixBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteEdgesMat(path, edges, vers);
	}

}
using namespace MY_DEBUG;


namespace TEST_SPARSE_MAT
{
	// �����Լ�д��ϡ�������ص����ӣ�
	void test00() 
	{
		// 1. traverseSparseMatrix()�������뺯���ӱ���ϡ������е�Ԫ�أ�
		Eigen::MatrixXi m1(Eigen::MatrixXi::Zero(4, 5));
		m1(1, 1) = 3;
		m1(1, 2) = 2;
		Eigen::SparseMatrix<int> sm1 = m1.sparseView();
		traverseSparseMatrix(sm1, [&](auto& iter)
			{
				std::cout << "value == ";
				disp<int>(iter.value());
				std::cout << "row == ";
				disp<int>(iter.row());
				std::cout << "col == ";
				disp<int>(iter.col());
				std::cout << std::endl;
			});
		std::cout << std::endl;

		//		���з���Ԫ��+1��
		traverseSparseMatrix(sm1, [&](auto& iter)
			{
				iter.valueRef()++;
			});


		//	2. dispSpMat()������ӡϡ��������Ԫ��
		dispSpMat(sm1);
 
		Eigen::SparseMatrix<int, Eigen::RowMajor> sm11 = m1.sparseView();
		Eigen::SparseMatrix<int> sm2;
		bool retFlag = spMatTranspose(sm2, sm1);
		dispSpMat(sm2);

		std::cout << "finished." << std::endl;
	}


	// ϡ�����Ĺ��졢��������
	void test0()
	{
		std::vector<Eigen::Triplet<float>> tripVec;
		tripVec.push_back(Eigen::Triplet<float>(0, 0, 1));
		tripVec.push_back(Eigen::Triplet<float>(1, 1, 1.1));
		tripVec.push_back(Eigen::Triplet<float>(2, 2, 1.2));
		tripVec.push_back(Eigen::Triplet<float>(3, 3, 1.3));
		tripVec.push_back(Eigen::Triplet<float>(4, 4, 1.4));
		tripVec.push_back(Eigen::Triplet<float>(1, 1, 9.0));

		//   1. ����ϡ�����

		//					ʹ����Ԫ����������ϡ����󡪡�setFromTriplets()
		Eigen::SparseMatrix<float> sm1(6, 7);
		sm1.setFromTriplets(tripVec.begin(), tripVec.end());
		dispMat(sm1.toDense());

		//					����Ԫ�صķ�ʽ����ϡ����󡪡�insert()
		Eigen::SparseMatrix<float> sm2(9, 10);
		sm2.insert(1, 2) = 1;
		sm2.insert(2, 3) = 2;
		sm2.insert(8, 3) = 88;
		sm2.insert(3, 4) = 3;
		sm2.insert(4, 4) = 99;

		//					����һЩ�����ϡ�����
		Eigen::SparseMatrix<float> sm3(4, 4);

		//			setIdentity()
		sm3.setIdentity();
		std::cout << "sm3 == \n" << sm3 << std::endl;


		// 1+. ѹ����It is worth noting that most of our wrappers to external libraries requires compressed matrices as inputs

		//			isCompressed()���������鿴�Ƿ���ѹ��״̬������Ԫ��֮�������ֶ�ѹ�����Ƿ�ѹ��״̬��
		std::cout << "sm1.isCompressed() == " << sm1.isCompressed() << std::endl;
		std::cout << "sm2.isCompressed() == " << sm2.isCompressed() << std::endl;
 
		//			 makeCompressed()��������ѹ��ʣ��ռ䣺
		sm2.makeCompressed();				
		std::cout << "sm2.isCompressed() == " << sm2.isCompressed() << std::endl;


		//	 2. ϡ����������

		//			nonZeros()
		std::cout << "sm1����Ԫ������" << sm1.nonZeros() << std::endl;

		//			innerSize()
		std::cout << "sm1��ά�ȣ�" << sm1.innerSize() << std::endl;			// �����ȵ�ϡ�����������������ȵ�ϡ����������
		
		//			outerSize()
		std::cout << "sm1��ά�ȣ�" << sm1.outerSize() << std::endl;			// �����ȵ�ϡ�����������������ȵ�ϡ���������� 


		//		valuePtr()������������ϡ�����Ԫ��������׵�ַ��
		auto dataPtr = sm1.valuePtr();
		std::cout << "valuePtr()������������ϡ�����Ԫ��������׵�ַ\n data of sm1: " << std::endl;
		for(unsigned i = 0; i < sm1.nonZeros(); ++i)
			std::cout << dataPtr[i] << std::endl;
		std::cout << std::endl;

		// ע��ò��ֱ���޸�dataPtrָ���������Σ�յġ����Ҫ��ĳԪ������Ӧ��ʹ��prune()������
		Eigen::SparseMatrix<float> sm11 = sm1;
		dispMat(sm1.toDense());
		std::cout << "sm1.nonZeros() == " << sm1.nonZeros() << std::endl;
		std::cout << "sm1.isCompressed() == " << sm1.isCompressed() << std::endl << std::endl;;

		//		prune����()������������������Ԫ�����㣻
		sm1.prune([&](const Eigen::Index& row, const Eigen::Index& col, const float& value)->bool
			{
				if (value > 1.3)
					return true;					// ����true��Ԫ�ر�������
				else
					return false;
			});
		std::cout << "prune����()������������������Ԫ�����㣻" << std::endl;
		dispMat(sm1.toDense());
		std::cout << "sm1.nonZeros() == " << sm1.nonZeros() << std::endl;
		std::cout << "sm1.isCompressed() == " << sm1.isCompressed() << std::endl;

		sm11.prune([&](const Eigen::Index& row, const Eigen::Index& col, const float& value)->bool
			{
				if (0 == row && 0 == col)
					return false;
				else
					return true;
			});
		dispMat(sm11.toDense());
		std::cout << "sm11.nonZeros() == " << sm11.nonZeros() << std::endl;


		//	3. ����ϡ������Ԫ��

		//					coeffRef(row, col)��������ϡ�����ĳһ��Ԫ��
		std::cout << " sm2.coeffRef(3, 4) == " << sm2.coeffRef(3, 4) << std::endl;
		std::cout << " sm2.coeffRef(3, 1) == " << sm2.coeffRef(3, 1) << std::endl;

		//					ʹ��InnerIterator����ϡ�����ķ���Ԫ��
		dispMat(sm2.toDense());
		std::cout << "sm2����Ԫ������" << sm2.nonZeros() << std::endl;
		int outerSize2 = sm2.outerSize();				// Ĭ�ϵ������ȴ洢����outerSize��������
		for (int colIdx = 0; colIdx < outerSize2; ++colIdx)
		{
			for (Eigen::SparseMatrix<float>::InnerIterator it(sm2, colIdx); it; ++it)			// �����ȴ洢ʱ��InnerIterator�����ڵ�������
			{
				std::cout << "value == " << it.value() << std::endl;					// Ԫ�ص�ֵ��
				std::cout << "valueRef == " << it.valueRef() << std::endl;			// Ԫ�ص����ã�
				std::cout << "row == " << it.row() << std::endl;			 // row index
				std::cout << "col == " << it.col() << std::endl;				 // col index (here it is equal to k)
				std::cout << "index == " << it.index() << std::endl;
				std::cout << std::endl;
				//std::cout << "" << it.index() << std::endl;				// inner index, here it is equal to it.row()
			}
		}

		// 4.  ϡ�����ͳ��ܾ����ת����
		Eigen::MatrixXf m1(3, 4);
		m1 << 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0;

		//			sparseView()
		sm1 = m1.sparseView();
		std::cout << "sm1 == \n" << sm1 << std::endl;
		std::cout << "sm1�Ƿ���ѹ����ʽ��" << sm1.isCompressed() << std::endl;

		//			toDense();
		Eigen::MatrixXf m2 = sm1.toDense();
		std::cout << "m2 == \n" << m2 << std::endl;
	}


	// ϡ�����Ļ����任�Ͳ���
	void test1()
	{
		// 1. ��͡���û��dense matrix���.colwise().sum()��.rowwise().sum()����������ʹ����/�ҳ�������ʵ�֣�
		Eigen::MatrixXf m1 = Eigen::MatrixXf::Random(4, 5);
		Eigen::MatrixXf m2(4, 5);
		m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20;
		std::cout << "m1 == \n" << m1 << std::endl;

		// 2. sparseView()������
		Eigen::SparseMatrix<float> sm1 = m1.sparseView();
		Eigen::SparseMatrix<float> sm2 = m2.sparseView();
		Eigen::RowVectorXf colSum = Eigen::RowVectorXf::Ones(4) * sm1;			// ���һ��ȫ1���������õ�������͵Ľ����
		Eigen::VectorXf rowSum = sm1 * Eigen::VectorXf::Ones(5);							// �ҳ�һ��ȫ1���������õ�������͵Ľ����
		std::cout << "������ͽ����\n" << colSum << std::endl;
		std::cout << "������ͽ����\n" << rowSum << std::endl << std::endl;

		// 3. �����ȴ洢ϡ������row()�������ص�������ֻ���ģ������Ա���ֵ��
		std::cout << "������i�У�" << std::endl;
		for (Eigen::SparseMatrix<float>::InnerIterator it(sm1, 1); it; ++it)				// ��ά�Ⱦ��Ǵ洢���ȵ�ά�ȣ�������ʹ��Ĭ�ϵ������ȴ洢�������it����ĳһ�еĵ�������
			std::cout << it.value() << ", ";
		std::cout << std::endl;

		Eigen::RowVectorXf tempVec = sm1.row(1);
		std::cout << "tempVec == \n" << tempVec << std::endl;

		sm1.col(0) = Eigen::SparseVector<float>(4);									// ĳһ�п�����ϡ��������ֵ����������dense������ֵ��
		std::cout << "sm1 == \n" << sm1 << std::endl;
		sm1.makeCompressed();
		std::cout << "sm1.nonZeros()  == " << sm1.nonZeros() << std::endl;

		sm1.col(1) = sm2.col(1);
		std::cout << "sm1 == \n" << sm1 << std::endl;
		sm1.makeCompressed();
		std::cout << "sm1.nonZeros()  == " << sm1.nonZeros() << std::endl; 

		// ʹ��sp2��ĳһ�и�sp1��ĳһ�и�ֵʱ�����뽫��������ת��Ϊ�����������������ұ����������쳣��
		sm1.col(0) = sm2.row(0).transpose();
		dispMat(sm1.toDense());

		// ��Ҫ��ϡ������Դ���tranpose()������������
		auto tmpSm = sm1.transpose();					// ֻ�ܸ�ֵ�������ȴ洢��ϡ�����

		std::cout << "finished." << std::endl;
	}


	// ϡ������ƴ��
	void test2()
	{
		Eigen::MatrixXf m1(2, 3), m2(2, 1), m3(1, 4), m4(3, 4);
		m1 << 1, 2, 3, 4, 5, 6;
		m2 << 0.1, 0.2;
		m3 << 1, 2, 3, 4;
		m4 = Eigen::MatrixXf::Ones(3, 4);
		Eigen::SparseMatrix<float> sm1(2, 3), sm2(2, 1), sm3(1, 4), sm4(3, 4), sm(4, 4);
		sm1 = m1.sparseView();
		sm2 = m2.sparseView();
		sm3 = m3.sparseView();
		sm4 = m4.sparseView();
		sm.leftCols(3) = sm1;			// �����д洢��ϡ�������˵��leftCols()��rightCols()���ص�����ֻ���ģ�topRows()��bottomRows()���ص����ǿɶ�д�ġ�
		sm.rightCols(1) = sm2;		// �д洢��ϡ�������֮

		//sm.topRows(1) = sm3;
		//sm.bottomRows(3) = sm4;

		std::cout << sm << std::endl;
	}


	// ϡ�����ķֽ�
	void test3() 
	{
		using spMatD = Eigen::SparseMatrix<double>;
		Eigen::MatrixXd m1(4, 4);
		spMatD sm1;
		// m1 << 1, 2, -3, 4, 5, 6, 7, 8, 9, 10, -11, 12, 13, 14, 15, 0;		// �������ǶԳ�
		m1 << 1, 2, 99, 4, 5, 6, 7, 8, 9, 10, 88, 77, 13, 14, 15, 66; 
		m1 = (m1 * m1.transpose()/1000).eval();					// �Գ�����
		debugDisp("m1 == \n", m1, "\n");

		sm1 = m1.sparseView();
		debugDisp("m1.determinant() == ", m1.determinant());

		// 1. ϡ������LU�ֽ�����
		Eigen::SparseLU<spMatD> solverLU;
		solverLU.compute(sm1);
		if (solverLU.info() != Eigen::Success)
		{
			debugDisp("LU decomposition failed");
			return;
		}
		debugDisp("solverLU.determinant() == ", solverLU.determinant());		// ϡ������಻����������ʽ�ķ�������Ҫͨ��LU��LDLT�ֽ�����á� 
 
		// 2. Eigen::SparseQR< MatrixType_, OrderingType_ >�ࡪ��ϡ������QR�ֽ���
		/*
			OrderingType_��ϡ�������������㷨����
					��ϡ��������������������Ż�����ֽ⡢�����������Է���������ܡ�
					�����������ʹ�þ��������ӷֽ���������ɵ��������С�����Ӷ����ټ��������ڴ�ռ�á�
					һЩ����㷨��������������ǳ����С�ͨ��ѡ����ʵ�������ԣ����ǿ��Լ�����ⲽ�������
					ĳЩ�������㷨�����������㷨����ֵ�ȶ��ԣ�������ֵ����Ӱ�졣

			OrderingType_�����ֿ��õ������㷨��
					1. Eigen::NaturalOrdering< StorageIndex >
								��Ȼ���򣬲���������ʹ��ԭ����
					2. Eigen::AMDOrdering< StorageIndex >
								Ӧ���ڶԳƾ���
								������С������Approximate Minimum Degree Ordering����ͨ������ÿ��ѡȡ��С�ȵ��н��������Լ�����䡣
								�еĶȡ������еķ���Ԫ�ص�������
					3. Eigen::COLAMDOrdering< StorageIndex >
								��Ӧ���ڷǶԳƾ���Ľ�����С������
								���ܱ�AMDOrdering����ɱ����ߡ�
				
		*/
		Eigen::SparseQR<spMatD, Eigen::COLAMDOrdering<int>> solverQR;
		solverQR.compute(sm1);
		if (solverQR.info() != Eigen::Success)
		{
			debugDisp("QR decomposition failed");
			return;
		}
		debugDisp("solverQR.rank() == ", solverQR.rank());			// ϡ������಻�������ȵķ�������Ҫͨ��QR��á� 

		// 3. ϡ������LDLT�ֽ�����LLT�ֽ����ܴ�����ӣ����Բ����ˣ�
		Eigen::SimplicialLDLT<spMatD> solverLDLT;
		solverLDLT.compute(sm1);
		if (solverLDLT.info() != Eigen::Success)		
		{
			// ������������ֻ�������Գƾ���ſ�����LDLT�ֽ⣬�������������Ҫ��ʱ����Ҳ����ֽ�ʧ�ܣ��������������ʽ�Ǵ���ġ�
			debugDisp("LDLT decomposition failed");
			return;
		}
		debugDisp("solverLDLT.determinant() == ", solverLDLT.determinant());	// ϡ������಻����������ʽ�ķ�������Ҫͨ��LU��LDLT�ֽ�����á� 

		debugDisp("test3 finished.");
	}


	// ���ϡ�����Է����顪��ֱ�ӷ�
	void test4()
	{
		// 1. ֱ�ӷ����ϡ�����Է�����

		//				ϡ������ʾ�����Է����飺Ax = b
		Eigen::SparseMatrix<float> A;
		Eigen::VectorXf b, x;
		// Eigen::SimplicialLLT<Eigen::SparseMatrix<float>> solver;					// ����LLT�ֽ���������һ�㲻��LDLT�ֽ�
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;					// ����LDLT�ֽ���������ʹ���ڳ��ܶȺܵ��ҹ�ģ���Ǻܴ�ľ��󣻣�����ò����Ҫ�Գƾ���
		// Eigen::SparseLU<Eigen::SparseMatrix<float> > solver;						// ����LU�ֽ���������
		// Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int> > solver;		// ����QR�ֽ����������Ƽ���С��������ʹ�á�

		//				д�뷽�������ݣ�
		Eigen::MatrixXf Adense(2, 2);
		b.resize(2);
		Adense << 1, 2, 1, -1;
		b << 0, 3;
		A = Adense.sparseView();
		debugDisp("A == ");
		dispSpMat(A);
		debugDisp("b == \n", b, "\n");
		debugDisp("det(A) == ", Adense.determinant()); 

		solver.compute(A);
		if (solver.info() != Eigen::Success)
		{
			std::cout << " decomposition failed" << std::endl;
			return;
		}
		x = solver.solve(b);
		if (solver.info() != Eigen::Success)
		{
			std::cout << "solving failed" << std::endl;
			return;
		}
		debugDisp("x == \n",  x, "\n");

		debugDisp("test3() finished.");
	}


	// ���ϡ�����Է����顪��������
	void test5()
	{
		//	1. ����ϡ�����Է����飺Ax = b
		Eigen::SparseMatrix<float> A;
		Eigen::VectorXf b(2), x; 
		Eigen::MatrixXf Adense(2, 2);
		Adense << 1, 2, 1, -1;
		b << 0, 3;
		A = Adense.sparseView();
		debugDisp("A == \n", A, "\n");
		debugDisp("b == \n", b, "\n"); 
		debugDisp("det(A) == ", Adense.determinant());

		// 2. ���������ϡ�����Է����飺
		/*
			�������������
					ConjugateGradient				
							������������ݶ�������� 
							Ҫ�����ΪSPD����Symmetric positive definite matrices �Գơ�������

					LeastSquaresConjugateGradient
							��С���˹����ݶ������
							�����ڳ����ξ��󣨷Ƿ���

					BiCGSTAB						
							����˫�����ݶ��������Iterative stabilized bi-conjugate gradient��			
							Ҫ�����Ϊ����
		*/
		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float>> solverLSCG;		
		Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> solverBCG;

		//			2.1 ��С���˹����ݶ������
		solverLSCG.setTolerance(0.00001);				//	���õ�������
		solverLSCG.compute(A);
		if (solverLSCG.info() != Eigen::Success)
		{
			debugDisp("!!! Eigen::LeastSquaresConjugateGradient<> computing failed.");
			return;
		}
		x = solverLSCG.solve(b); 
		debugDisp("Eigen::LeastSquaresConjugateGradient<> result�� x == \n", x, "\n\n");
		debugDisp("solverLSCG.iterations() == ", solverLSCG.iterations());					// ��������
		debugDisp("solverLSCG.error() == ", solverLSCG.error());

		//			2.2 ����˫�����ݶ������
		solverBCG.compute(A);
		if(solverBCG.info() != Eigen::Success)
		{
			debugDisp("!!! Eigen::BiCGSTAB<> computing failed.");
			return;
		}
		x = solverBCG.solve(b);
		debugDisp("Eigen::BiCGSTAB<> result�� x == \n", x, "\n");
		debugDisp("solverBCG.iterations() == ", solverBCG.iterations());
		debugDisp("solverBCG.error() == ", solverBCG.error());

		debugDisp("\n\ntest4() finished.");
	}


	// �Ƚϲ�ͬ���ϵ�����Է����鷽�����ٶȣ�
	void test6() 
	{
		// ����һ��3n x 3n��A��������
		tiktok& tt = tiktok::getInstance();
		const double tolerance = 1e-3;			// �趨�������ľ��ȣ�
		double interval = 0;								// ��¼ʱ������

		// 1. ����ϵ�����Է������A�����b������
		const int n = 10000;
		const int rows = n * 3;
		const int cols = n * 3;
		int offset = 0;
		std::default_random_engine e;											// ��������������������
		std::uniform_real_distribution<double> URD_d(-100, 100);
		std::vector<Eigen::Triplet<double>> tripVec;
		Eigen::SparseMatrix<double> A(rows, cols);
		Eigen::VectorXd b(rows);
		b.setRandom();

		tt.start();
		tripVec.reserve(n * 9);
		for (int i = 0; i < n; ++i) 
		{
			offset = 3 * i;
			tripVec.push_back(Eigen::Triplet<double>(offset, offset, URD_d(e)));
			tripVec.push_back(Eigen::Triplet<double>(offset, offset + 1, URD_d(e)));
			tripVec.push_back(Eigen::Triplet<double>(offset, offset + 2, URD_d(e)));
			tripVec.push_back(Eigen::Triplet<double>(offset + 1, offset + 1, URD_d(e)));
			tripVec.push_back(Eigen::Triplet<double>(offset + 1, offset + 2, URD_d(e)));
			tripVec.push_back(Eigen::Triplet<double>(offset + 1, offset + 1, URD_d(e)));
			tripVec.push_back(Eigen::Triplet<double>(offset + 2, offset + 2, URD_d(e)));
			tripVec.push_back(Eigen::Triplet<double>(offset + 2, offset + 1, URD_d(e)));
			tripVec.push_back(Eigen::Triplet<double>(offset + 2, offset + 2, URD_d(e)));
		}
		A.setFromTriplets(tripVec.begin(), tripVec.end());
		tt.endCout("1. ����ϡ�����A��ʱ��");
		debugDisp("\n");

		// 2. QR�ֽⷨ
		debugDisp("QR�ֽⷨ��");
		{
			Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solverQR;
			interval = 0;

			//			2.1 QR�ֽ⣺
			tt.start();
			solverQR.compute(A);
			if (solverQR.info() != Eigen::Success)
			{
				std::cout << "QR decomposition failed" << std::endl;
				return;
			}
			tt.endCout("A����QR�ֽ��ʱ��");
			interval += tt.lastDur;

			debugDisp("solverQR.rank() == ", solverQR.rank());
			if (solverQR.rank() < 3 * n)
			{
				debugDisp("������A�������ȡ�");
				return;
			}

			//			2.2 �����Է����飺
			tt.start();
			Eigen::VectorXd x = solverQR.solve(b);
			tt.endCout("��ϡ�����Է������ʱ��");
			interval += tt.lastDur;

			//			2.3 ��֤���׼ȷ�ԣ�
			Eigen::VectorXd b_prime = A * x;
			Eigen::VectorXd residual = b - b_prime;
			double sqrErr = residual.squaredNorm();				// ƽ����
			debugDisp("ƽ�����sqrErr == ", sqrErr);

			debugDisp("QR�ֽⷨ�ܺ�ʱ��", interval, "s\n\n");
		}

		// 3. LU�ֽⷨ��
		debugDisp("LU�ֽⷨ��");
		{
			Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solverLU;
			interval = 0;

			//			3.1 LU�ֽ⣺
			tt.start();
			solverLU.compute(A);
			if (solverLU.info() != Eigen::Success)
			{
				std::cout << "LU decomposition failed" << std::endl;
				return;
			}
			tt.endCout("A����LU�ֽ��ʱ��");
			interval += tt.lastDur; 

			//			3.2 �����Է����飺
			tt.start();
			Eigen::VectorXd x = solverLU.solve(b);
			tt.endCout("��ϡ�����Է������ʱ��");
			interval += tt.lastDur;

			//			3.3 ��֤���׼ȷ�ԣ�
			Eigen::VectorXd b_prime = A * x;
			Eigen::VectorXd residual = b - b_prime;
			double sqrErr = residual.squaredNorm();				// ƽ����
			debugDisp("ƽ�����sqrErr == ", sqrErr);

			debugDisp("LU�ֽⷨ�ܺ�ʱ��", interval, "s\n\n");
		}

		// 4. ��С���˹����ݶȵ�������
		debugDisp("��С���˹����ݶȵ�������");
		{
			Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solverLSCG;
			Eigen::VectorXd x;
			interval = 0;
			solverLSCG.setTolerance(tolerance);

			// 4.1 compute
			tt.start();
			solverLSCG.compute(A);
			tt.end();
			interval += tt.lastDur;
			if (solverLSCG.info() != Eigen::Success)
			{
				debugDisp("!!! Eigen::LeastSquaresConjugateGradient<> computing failed.");
				return;
			}
			
			// 4.2 solve
			tt.start();
			x = solverLSCG.solve(b);
			tt.end();
			interval += tt.lastDur;
			debugDisp("solverLSCG.iterations() == ", solverLSCG.iterations());							// ��������
			debugDisp("solverLSCG.error() == ", solverLSCG.error());			

			// 4.3 ��֤���׼ȷ�ԣ�
			Eigen::VectorXd b_prime = A * x;
			Eigen::VectorXd residual = b - b_prime;
			double sqrErr = residual.squaredNorm();				// ƽ����
			debugDisp("ƽ�����sqrErr == ", sqrErr);

			debugDisp("��С���˹����ݶȵ������ܺ�ʱ��", interval, "s\n\n");
		}

		// 5. ˫�����ݶȵ�������
		debugDisp("˫�����ݶȵ�������");
		{
			Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solverBCG;
			Eigen::VectorXd x;
			interval = 0;
			solverBCG.setTolerance(tolerance);

			// 5.1 compute
			tt.start();
			solverBCG.compute(A);
			tt.end();
			interval += tt.lastDur;
			if (solverBCG.info() != Eigen::Success)
			{
				debugDisp("!!! Eigen::BiCGSTAB<Eigen::SparseMatrix<> computing failed.");
				return;
			}

			// 5.2 solve
			tt.start();
			x = solverBCG.solve(b);
			tt.end();
			interval += tt.lastDur;
			debugDisp("solverBCG.iterations() == ", solverBCG.iterations());							// ��������
			debugDisp("solverBCG.error() == ", solverBCG.error());

			// 5.3 ��֤���׼ȷ�ԣ�
			Eigen::VectorXd b_prime = A * x;
			Eigen::VectorXd residual = b - b_prime;
			double sqrErr = residual.squaredNorm();				// ƽ����
			debugDisp("ƽ�����sqrErr == ", sqrErr);

			debugDisp("˫�����ݶȵ������ܺ�ʱ��", interval, "s\n\n");
		}


		debugDisp("\n\ntest6() finished.");
	}
}

