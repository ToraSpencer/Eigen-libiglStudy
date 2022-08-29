#include "sparse_mat.h"

namespace SPARSEMAT
{
	using spMatf = SparseMatrix<float, ColMajor>;
	using TripF = Eigen::Triplet<float>;

	// ϡ�����
	void test0()
	{
		std::vector<TripF> tripList;
		tripList.push_back(TripF(0, 0, 1));
		tripList.push_back(TripF(1, 1, 1.1));
		tripList.push_back(TripF(2, 2, 1.2));
		tripList.push_back(TripF(3, 3, 1.3));
		tripList.push_back(TripF(4, 4, 1.4));

		//   1. ����ϡ�����

		//					ʹ����Ԫ����������ϡ����󡪡�setFromTriplets()
		spMatf sm1(120, 100);
		sm1.setFromTriplets(tripList.begin(), tripList.end());

		//					����Ԫ�صķ�ʽ����ϡ����󡪡�insert()
		spMatf sm2(120, 100);
		sm2.insert(1, 2) = 1;
		sm2.insert(2, 3) = 2;
		sm2.insert(3, 4) = 3;
		sm2.makeCompressed();					// ѹ��ʣ��ռ䣺

		//					����һЩ�����ϡ�����
		spMatf sm3(4, 4);
		sm3.setIdentity();
		std::cout << "sm3 == \n" << sm3 << std::endl;

		//	 2. ϡ����������
		std::cout << "sm1����Ԫ������" << sm1.nonZeros() << std::endl;
		std::cout << "sm1��ά�ȣ�" << sm1.innerSize() << std::endl;			// �����ȵ�ϡ�����������������ȵ�ϡ����������
		std::cout << "sm1��ά�ȣ�" << sm1.outerSize() << std::endl;			// �����ȵ�ϡ�����������������ȵ�ϡ����������
		std::cout << "sm1�Ƿ���ѹ����ʽ��" << sm1.isCompressed() << std::endl;


		//	3. ����ϡ������Ԫ��

		//					coeffRef(row, col)��������ϡ�����ĳһ��Ԫ��
		std::cout << " sm2.coeffRef(3, 4) == " << sm2.coeffRef(3, 4) << std::endl;

		//					ʹ��InnerIterator����ϡ�����ķ���Ԫ��
		std::cout << "sm2����Ԫ������" << sm2.nonZeros() << std::endl;
		for (int k = 0; k < sm2.outerSize(); ++k)
		{
			for (spMatf::InnerIterator it(sm2, k); it; ++it)
			{
				std::cout << "value == " << it.value() << std::endl;
				std::cout << "row == " << it.row() << std::endl;			 // row index
				std::cout << "col == " << it.col() << std::endl;			 // col index (here it is equal to k)
				std::cout << std::endl;
				//std::cout << "" << it.index() << std::endl;			// inner index, here it is equal to it.row()
			}
		}


		// 4.  ϡ�����ͳ��ܾ����ת����
		MatrixXf m1(3, 4);
		m1 << 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0;

		//			sparseView()
		sm1 = m1.sparseView();
		std::cout << "sm1 == \n" << sm1 << std::endl;
		std::cout << "sm1�Ƿ���ѹ����ʽ��" << sm1.isCompressed() << std::endl;

		//			toDense();
		MatrixXf m2 = sm1.toDense();
		std::cout << "m2 == \n" << m2 << std::endl;

	}


	// ϡ�����Ļ����任�Ͳ���
	void test1()
	{
		// ��͡���û��dense matrix���.colwise().sum()��.rowwise().sum()����������ʹ����/�ҳ�������ʵ�֣�
		MatrixXf m1 = MatrixXf::Random(4, 4);
		MatrixXf m2(4, 4);
		m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
		std::cout << "m1 == \n" << m1 << std::endl;
		SparseMatrix<float> sm1 = m1.sparseView();
		SparseMatrix<float> sm2 = m2.sparseView();
		RowVectorXf colSum = RowVectorXf::Ones(4) * sm1;			// ���һ��ȫ1���������õ�������͵Ľ����
		VectorXf rowSum = sm1 * VectorXf::Ones(4);							// �ҳ�һ��ȫ1���������õ�������͵Ľ����
		std::cout << "������ͽ����\n" << colSum << std::endl;
		std::cout << "������ͽ����\n" << rowSum << std::endl << std::endl;

		// �����ȴ洢ϡ������row()�������ص�������ֻ���ģ������Ա���ֵ��
		std::cout << "������i�У�" << std::endl;
		for (spMatf::InnerIterator it(sm1, 1); it; ++it)				// ��ά�Ⱦ��Ǵ洢���ȵ�ά�ȣ�������ʹ��Ĭ�ϵ������ȴ洢�������it����ĳһ�еĵ�������
			std::cout << it.value() << ", ";
		std::cout << std::endl;

		RowVectorXf tempVec = sm1.row(1);
		std::cout << "tempVec == \n" << tempVec << std::endl;

		sm1.col(0) = SparseVector<float>(4);									// ĳһ�п�����ϡ��������ֵ����������dense������ֵ��
		std::cout << "sm1 == \n" << sm1 << std::endl;
		sm1.makeCompressed();
		std::cout << "sm1.nonZeros()  == " << sm1.nonZeros() << std::endl;

		sm1.col(1) = sm2.col(1);
		std::cout << "sm1 == \n" << sm1 << std::endl;
		sm1.makeCompressed();
		std::cout << "sm1.nonZeros()  == " << sm1.nonZeros() << std::endl;

		// ʹ��sp2��ĳһ�и�sp1��ĳһ�и�ֵʱ�����뽫��������ת��Ϊ�����������������ұ����������쳣��
		sm1.col(0) = sm2.row(0).transpose();
		dispMat<float>(sm1.toDense());
	}


	// ϡ������ƴ��
	void test2()
	{
		MatrixXf m1(2, 3), m2(2, 1), m3(1, 4), m4(3, 4);
		m1 << 1, 2, 3, 4, 5, 6;
		m2 << 0.1, 0.2;
		m3 << 1, 2, 3, 4;
		m4 = MatrixXf::Ones(3, 4);
		spMatf sm1(2, 3), sm2(2, 1), sm3(1, 4), sm4(3, 4), sm(4, 4);
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
		SparseMatrix<float> A;
		VectorXf b(2), x;
		///SimplicialLLT<SparseMatrix<float>> solver;									//����LLT�ֽ���������һ���Ƽ�ʹ�ô��������
		//SparseLU<SparseMatrix<float> > solver;											// ����LU�ֽ���������
		SparseQR<SparseMatrix<float>, COLAMDOrdering<int> > solver;		// ����QR�ֽ����������Ƽ���С��������ʹ�á�

		//				д�뷽�������ݣ�
		MatrixXf Adense(2, 2);
		Adense << 1, 2, 1, -1;
		b << 0, 3;
		A = Adense.sparseView();
		std::cout << "A == \n" << A << std::endl;
		std::cout << "b == \n" << b << std::endl;
		std::cout << "det(A) == " << Adense.determinant() << std::endl;

		solver.compute(A);
		if (solver.info() != Success)
		{
			std::cout << " decomposition failed" << std::endl;
			return;
		}
		x = solver.solve(b);
		if (solver.info() != Success)
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
		LeastSquaresConjugateGradient<Eigen::SparseMatrix<float> > IterSolver;

		// ���õ�������
		IterSolver.setTolerance(0.00001);
		IterSolver.compute(A);
		if (IterSolver.info() != Success)
		{
			return;
		}

		x = IterSolver.solve(b);
		if (IterSolver.info() != Success)
		{
			return;
		}
		std::cout << "x == \n" << x << std::endl;

	}

}

