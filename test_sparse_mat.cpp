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


namespace SPARSEMAT
{
	// �����Լ�д��ϡ�������ص����ӣ�
	void test00() 
	{
		// traverseSparseMatrix()�������뺯���ӱ���ϡ������е�Ԫ�أ�
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


		//	dispSpMat()������ӡϡ��������Ԫ��
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
		std::vector<Eigen::Triplet<float>> tripList;
		tripList.push_back(Eigen::Triplet<float>(0, 0, 1));
		tripList.push_back(Eigen::Triplet<float>(1, 1, 1.1));
		tripList.push_back(Eigen::Triplet<float>(2, 2, 1.2));
		tripList.push_back(Eigen::Triplet<float>(3, 3, 1.3));
		tripList.push_back(Eigen::Triplet<float>(4, 4, 1.4));
		tripList.push_back(Eigen::Triplet<float>(1, 1, 9.0));

		//   1. ����ϡ�����

		//					ʹ����Ԫ����������ϡ����󡪡�setFromTriplets()
		Eigen::SparseMatrix<float> sm1(6, 7);
		sm1.setFromTriplets(tripList.begin(), tripList.end());
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
		// ��͡���û��dense matrix���.colwise().sum()��.rowwise().sum()����������ʹ����/�ҳ�������ʵ�֣�
		Eigen::MatrixXf m1 = Eigen::MatrixXf::Random(4, 5);
		Eigen::MatrixXf m2(4, 5);
		m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20;
		std::cout << "m1 == \n" << m1 << std::endl;

		// sparseView()������
		Eigen::SparseMatrix<float> sm1 = m1.sparseView();
		Eigen::SparseMatrix<float> sm2 = m2.sparseView();
		Eigen::RowVectorXf colSum = Eigen::RowVectorXf::Ones(4) * sm1;			// ���һ��ȫ1���������õ�������͵Ľ����
		Eigen::VectorXf rowSum = sm1 * Eigen::VectorXf::Ones(5);							// �ҳ�һ��ȫ1���������õ�������͵Ľ����
		std::cout << "������ͽ����\n" << colSum << std::endl;
		std::cout << "������ͽ����\n" << rowSum << std::endl << std::endl;

		// �����ȴ洢ϡ������row()�������ص�������ֻ���ģ������Ա���ֵ��
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


	// ���ϡ�����Է�����
	void test3()
	{
		// 1. ֱ�ӷ����ϡ�����Է�����

		//				ϡ������ʾ�����Է����飺Ax = b
		Eigen::SparseMatrix<float> A;
		Eigen::VectorXf b(2), x;
		///Eigen::SimplicialLLT<Eigen::SparseMatrix<float>> solver;									//����LLT�ֽ���������һ���Ƽ�ʹ�ô��������
		//Eigen::SparseLU<Eigen::SparseMatrix<float> > solver;											// ����LU�ֽ���������
		Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int> > solver;		// ����QR�ֽ����������Ƽ���С��������ʹ�á�

		//				д�뷽�������ݣ�
		Eigen::MatrixXf Adense(2, 2);
		Adense << 1, 2, 1, -1;
		b << 0, 3;
		A = Adense.sparseView();
		std::cout << "A == \n" << A << std::endl;
		std::cout << "b == \n" << b << std::endl;
		std::cout << "det(A) == " << Adense.determinant() << std::endl;

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

		std::cout << "x == \n" << x << std::endl;


		// 2. ���������ϡ�����Է����飺

		/*
			�������������
				ConjugateGradient				������������ݶ�������� Ҫ�����ΪSPD����Symmetric positive definite matrices �Գơ�������
				LeastSquaresConjugateGradient
				BiCGSTAB						����˫�����ݶ��������			Ҫ�����Ϊ����
		*/
		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float> > IterSolver;

		// ���õ�������
		IterSolver.setTolerance(0.00001);
		IterSolver.compute(A);
		if (IterSolver.info() != Eigen::Success)
			return;

		x = IterSolver.solve(b);
		if (IterSolver.info() != Eigen::Success)
			return;
		std::cout << "x == \n" << x << std::endl;
	}
}

