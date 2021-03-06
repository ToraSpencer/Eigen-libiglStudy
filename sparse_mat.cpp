#include "sparse_mat.h"

namespace SPARSEMAT
{
	using spMatf = SparseMatrix<float, ColMajor>;
	using TripF = Eigen::Triplet<float>;

	// 稀疏矩阵
	void test0()
	{
		std::vector<TripF> tripList;
		tripList.push_back(TripF(0, 0, 1));
		tripList.push_back(TripF(1, 1, 1.1));
		tripList.push_back(TripF(2, 2, 1.2));
		tripList.push_back(TripF(3, 3, 1.3));
		tripList.push_back(TripF(4, 4, 1.4));

		//   1. 生成稀疏矩阵：

		//					使用三元数数组生成稀疏矩阵——setFromTriplets()
		spMatf sm1(120, 100);
		sm1.setFromTriplets(tripList.begin(), tripList.end());

		//					插入元素的方式生成稀疏矩阵——insert()
		spMatf sm2(120, 100);
		sm2.insert(1, 2) = 1;
		sm2.insert(2, 3) = 2;
		sm2.insert(3, 4) = 3;
		sm2.makeCompressed();					// 压缩剩余空间：

		//					生成一些特殊的稀疏矩阵：
		spMatf sm3(4, 4);
		sm3.setIdentity();
		std::cout << "sm3 == \n" << sm3 << std::endl;

		//	 2. 稀疏矩阵的属性
		std::cout << "sm1非零元素数：" << sm1.nonZeros() << std::endl;
		std::cout << "sm1内维度：" << sm1.innerSize() << std::endl;			// 列优先的稀疏矩阵的列数，行优先的稀疏矩阵的行数
		std::cout << "sm1外维度：" << sm1.outerSize() << std::endl;			// 列优先的稀疏矩阵的行数，行优先的稀疏矩阵的列数
		std::cout << "sm1是否是压缩格式：" << sm1.isCompressed() << std::endl;


		//	3. 访问稀疏矩阵的元素

		//					coeffRef(row, col)——访问稀疏矩阵某一个元素
		std::cout << " sm2.coeffRef(3, 4) == " << sm2.coeffRef(3, 4) << std::endl;

		//					使用InnerIterator访问稀疏矩阵的非零元素
		std::cout << "sm2非零元素数：" << sm2.nonZeros() << std::endl;
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


		// 4.  稀疏矩阵和稠密矩阵的转换：
		MatrixXf m1(3, 4);
		m1 << 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0;

		//			sparseView()
		sm1 = m1.sparseView();
		std::cout << "sm1 == \n" << sm1 << std::endl;
		std::cout << "sm1是否是压缩格式：" << sm1.isCompressed() << std::endl;

		//			toDense();
		MatrixXf m2 = sm1.toDense();
		std::cout << "m2 == \n" << m2 << std::endl;

	}


	// 稀疏矩阵的基本变换和操作
	void test1()
	{
		// 求和——没有dense matrix里的.colwise().sum()和.rowwise().sum()操作，可以使用左/右乘向量来实现：
		MatrixXf m1 = MatrixXf::Random(4, 4);
		MatrixXf m2(4, 4);
		m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
		std::cout << "m1 == \n" << m1 << std::endl;
		SparseMatrix<float> sm1 = m1.sparseView();
		SparseMatrix<float> sm2 = m2.sparseView();
		RowVectorXf colSum = RowVectorXf::Ones(4) * sm1;			// 左乘一个全1行向量，得到按列求和的结果。
		VectorXf rowSum = sm1 * VectorXf::Ones(4);							// 右乘一个全1列向量，得到按行求和的结果。
		std::cout << "按列求和结果：\n" << colSum << std::endl;
		std::cout << "按行求和结果：\n" << rowSum << std::endl << std::endl;

		// 列优先存储稀疏矩阵的row()方法返回的引用是只读的，不可以被赋值。
		std::cout << "遍历第i列：" << std::endl;
		for (spMatf::InnerIterator it(sm1, 1); it; ++it)				// 内维度就是存储优先的维度，如这里使用默认的列优先存储，这里的it就是某一列的迭代器；
		{
			std::cout << it.value() << ", ";
		}
		std::cout << std::endl;

		RowVectorXf tempVec = sm1.row(1);
		std::cout << "tempVec == \n" << tempVec << std::endl;

		sm1.col(0) = SparseVector<float>(4);									// 某一列可以用稀疏向量赋值，不可以用dense向量赋值。
		std::cout << "sm1 == \n" << sm1 << std::endl;
		sm1.makeCompressed();
		std::cout << "sm1.nonZeros()  == " << sm1.nonZeros() << std::endl;

		sm1.col(1) = sm2.col(1);
		std::cout << "sm1 == \n" << sm1 << std::endl;
		sm1.makeCompressed();
		std::cout << "sm1.nonZeros()  == " << sm1.nonZeros() << std::endl;

		// 使用sp2的某一行给sp1的某一列赋值时，必须将该行向量转置为列向量，否则会出错，且编译器不报异常；
		sm1.col(0) = sm2.row(0).transpose();
		dispMat<float>(sm1.toDense());


	}


	// 稀疏矩阵的拼接
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
		sm.leftCols(3) = sm1;			// 对于行存储的稀疏矩阵来说，leftCols()和rightCols()返回的列是只读的，topRows()和bottomRows()返回的行是可读写的。
		sm.rightCols(1) = sm2;		// 列存储的稀疏矩阵则反之
		//sm.topRows(1) = sm3;
		//sm.bottomRows(3) = sm4;

		std::cout << sm << std::endl;
	}


	// 求解稀疏线性方程组
	void test3()
	{
		// 1. 直接法求解稀疏线性方程组

		//				稀疏矩阵表示的线性方程组：Ax = b
		SparseMatrix<float> A;
		VectorXf b(2), x;
		///SimplicialLLT<SparseMatrix<float>> solver;									//基于LLT分解的求解器，一般推荐使用此求解器。
		//SparseLU<SparseMatrix<float> > solver;											// 基于LU分解的求解器。
		SparseQR<SparseMatrix<float>, COLAMDOrdering<int> > solver;		// 基于QR分解的求解器，推荐最小方差问题使用。

		//				写入方程组数据：
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


		// 2. 迭代法求解稀疏线性方程组：

		/*
			迭代器求解器：
				ConjugateGradient				经典迭代共轭梯度求解器。 要求矩阵为SPD矩阵（Symmetric positive definite matrices 对称、正定）
				LeastSquaresConjugateGradient
				BiCGSTAB						迭代双共轭梯度求解器。			要求矩阵为方阵。
		*/
		LeastSquaresConjugateGradient<Eigen::SparseMatrix<float> > IterSolver;

		// 设置迭代精度
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

