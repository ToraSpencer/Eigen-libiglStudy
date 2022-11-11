#include "sparse_mat.h"

namespace SPARSEMAT
{
	// �����Լ�д��ϡ�������ص����ӣ�
	void test00() 
	{
		// traverseSparseMatrix()�������뺯���ӱ���ϡ������е�Ԫ�أ�
		Eigen::Matrix3i m1(Eigen::Matrix3i::Zero());
		m1(1, 1) = 3;
		m1(1, 2) = 2;
		Eigen::SparseMatrix<int> sm1 = m1.sparseView();
		traverseSparseMatrix(sm1, [&](auto& iter)
			{
				disp<int>(iter.value());
			});

		traverseSparseMatrix(sm1, [&](auto& iter)
			{
				iter.valueRef()++;
			});

		traverseSparseMatrix(sm1, [&](auto& iter)
			{
				disp<int>(iter.value());
			});

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
		SparseMatrix<float> sm1(6, 7);
		sm1.setFromTriplets(tripList.begin(), tripList.end());
		dispMat(sm1.toDense());

		//					����Ԫ�صķ�ʽ����ϡ����󡪡�insert()
		SparseMatrix<float> sm2(120, 100);
		sm2.insert(1, 2) = 1;
		sm2.insert(2, 3) = 2;
		sm2.insert(8, 3) = 88;
		sm2.insert(3, 4) = 3;
		sm2.insert(4, 4) = 99;

		//					����һЩ�����ϡ�����
		SparseMatrix<float> sm3(4, 4);

		//			setIdentity()
		sm3.setIdentity();
		std::cout << "sm3 == \n" << sm3 << std::endl;

		// 1+. ѹ����

		//			isCompressed()���������鿴�Ƿ���ѹ��״̬������Ԫ��֮�������ֶ�ѹ�����Ƿ�ѹ��״̬��
		std::cout << "sm1.isCompressed() == " << sm1.isCompressed() << std::endl;
		std::cout << "sm2.isCompressed() == " << sm2.isCompressed() << std::endl;
 
		//			 makeCompressed()��������ѹ��ʣ��ռ䣺
		sm2.makeCompressed();				
		std::cout << "sm2.isCompressed() == " << sm2.isCompressed() << std::endl;


		//	 2. ϡ����������
		std::cout << "sm1����Ԫ������" << sm1.nonZeros() << std::endl;
		std::cout << "sm1��ά�ȣ�" << sm1.innerSize() << std::endl;			// �����ȵ�ϡ�����������������ȵ�ϡ����������
		std::cout << "sm1��ά�ȣ�" << sm1.outerSize() << std::endl;			// �����ȵ�ϡ�����������������ȵ�ϡ����������
		std::cout << "sm1�Ƿ���ѹ����ʽ��" << sm1.isCompressed() << std::endl;

		//		valuePtr()������������ϡ�����Ԫ��������׵�ַ��
		auto dataPtr = sm1.valuePtr();
		std::cout << "data of sm1: " << std::endl;
		for(unsigned i = 0; i < sm1.nonZeros(); ++i)
			std::cout << dataPtr[i] << std::endl;
		std::cout << std::endl;

		// ע��ò��ֱ���޸�dataPtrָ���������Σ�յġ����Ҫ��ĳԪ������Ӧ��ʹ��prune()������
		SparseMatrix<float> sm11 = sm1;

		//		prune����()������������������Ԫ�����㣻
		sm1.prune([&](const Index& row, const Index& col, const float& value)->bool
			{
				if (value > 1.3)
					return true;			// ����true��Ԫ�ر�������
				else
					return false;
			});
		dispMat(sm1.toDense());
		std::cout << "sm1.nonZeros() == " << sm1.nonZeros() << std::endl;
		std::cout << "sm1.isCompressed() == " << sm1.isCompressed() << std::endl;

		sm11.prune([&](const Index& row, const Index& col, const float& value)->bool
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

		//					ʹ��InnerIterator����ϡ�����ķ���Ԫ��
		std::cout << "sm2����Ԫ������" << sm2.nonZeros() << std::endl;
		int outerSize2 = sm2.outerSize();				// Ĭ�ϵ������ȴ洢����outerSize��������
		for (int k = 0; k < outerSize2; ++k)
		{
			for (SparseMatrix<float>::InnerIterator it(sm2, k); it; ++it)			// �����ȴ洢ʱ��InnerIterator�����ڵ�������
			{
				std::cout << "value == " << it.value() << std::endl;					// Ԫ�ص�ֵ��
				std::cout << "valueRef == " << it.valueRef() << std::endl;			// Ԫ�ص����ã�
				std::cout << "row == " << it.row() << std::endl;			 // row index
				std::cout << "col == " << it.col() << std::endl;				 // col index (here it is equal to k)
				std::cout << std::endl;
				//std::cout << "" << it.index() << std::endl;				// inner index, here it is equal to it.row()
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

		// sparseView()������
		SparseMatrix<float> sm1 = m1.sparseView();
		SparseMatrix<float> sm2 = m2.sparseView();
		RowVectorXf colSum = RowVectorXf::Ones(4) * sm1;			// ���һ��ȫ1���������õ�������͵Ľ����
		VectorXf rowSum = sm1 * VectorXf::Ones(4);							// �ҳ�һ��ȫ1���������õ�������͵Ľ����
		std::cout << "������ͽ����\n" << colSum << std::endl;
		std::cout << "������ͽ����\n" << rowSum << std::endl << std::endl;

		// �����ȴ洢ϡ������row()�������ص�������ֻ���ģ������Ա���ֵ��
		std::cout << "������i�У�" << std::endl;
		for (SparseMatrix<float>::InnerIterator it(sm1, 1); it; ++it)				// ��ά�Ⱦ��Ǵ洢���ȵ�ά�ȣ�������ʹ��Ĭ�ϵ������ȴ洢�������it����ĳһ�еĵ�������
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
		dispMat(sm1.toDense());

		// ò�������ȵ�ϡ������transpose()�������صľ����������ȵģ�
		SparseMatrix<float> tmpSm = sm1.transpose();
		sm2 = sm1 + tmpSm;						// sm2 = sm1 + sm1.transpose() ����벻ͨ����
		std::cout << "sm2 = sm1 + tmpSm: " << std::endl;
		dispSpMat(sm2);
	}


	// ϡ������ƴ��
	void test2()
	{
		MatrixXf m1(2, 3), m2(2, 1), m3(1, 4), m4(3, 4);
		m1 << 1, 2, 3, 4, 5, 6;
		m2 << 0.1, 0.2;
		m3 << 1, 2, 3, 4;
		m4 = MatrixXf::Ones(3, 4);
		SparseMatrix<float> sm1(2, 3), sm2(2, 1), sm3(1, 4), sm4(3, 4), sm(4, 4);
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
			return;

		x = IterSolver.solve(b);
		if (IterSolver.info() != Success)
			return;
		std::cout << "x == \n" << x << std::endl;
	}
}

