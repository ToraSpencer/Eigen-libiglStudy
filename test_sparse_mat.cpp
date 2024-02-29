#include "test_sparse_mat.h"


////////////////////////////////////////////////////////////////////////////////////////////// DEBUG 接口
namespace MY_DEBUG
{
	static std::string g_debugPath = "E:/";

	// lambda——打印std::cout支持的类型变量。
	template <typename T>
	static auto disp = [](const T& arg)
	{
		std::cout << arg << ", ";
	};

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
	// 测试自己写的稀疏矩阵相关的轮子：
	void test00() 
	{
		// 1. traverseSparseMatrix()——传入函数子遍历稀疏矩阵中的元素；
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

		//		所有非零元素+1；
		traverseSparseMatrix(sm1, [&](auto& iter)
			{
				iter.valueRef()++;
			});


		//	2. dispSpMat()——打印稀疏矩阵非零元素
		dispSpMat(sm1);
 
		Eigen::SparseMatrix<int, Eigen::RowMajor> sm11 = m1.sparseView();
		Eigen::SparseMatrix<int> sm2;
		bool retFlag = spMatTranspose(sm2, sm1);
		dispSpMat(sm2);

		std::cout << "finished." << std::endl;
	}


	// 稀疏矩阵的构造、基本属性
	void test0()
	{
		std::vector<Eigen::Triplet<float>> tripVec;
		tripVec.push_back(Eigen::Triplet<float>(0, 0, 1));
		tripVec.push_back(Eigen::Triplet<float>(1, 1, 1.1));
		tripVec.push_back(Eigen::Triplet<float>(2, 2, 1.2));
		tripVec.push_back(Eigen::Triplet<float>(3, 3, 1.3));
		tripVec.push_back(Eigen::Triplet<float>(4, 4, 1.4));
		tripVec.push_back(Eigen::Triplet<float>(1, 1, 9.0));

		//   1. 生成稀疏矩阵：

		//					使用三元数数组生成稀疏矩阵——setFromTriplets()
		Eigen::SparseMatrix<float> sm1(6, 7);
		sm1.setFromTriplets(tripVec.begin(), tripVec.end());
		dispMat(sm1.toDense());

		//					插入元素的方式生成稀疏矩阵——insert()
		Eigen::SparseMatrix<float> sm2(9, 10);
		sm2.insert(1, 2) = 1;
		sm2.insert(2, 3) = 2;
		sm2.insert(8, 3) = 88;
		sm2.insert(3, 4) = 3;
		sm2.insert(4, 4) = 99;

		//					生成一些特殊的稀疏矩阵：
		Eigen::SparseMatrix<float> sm3(4, 4);

		//			setIdentity()
		sm3.setIdentity();
		std::cout << "sm3 == \n" << sm3 << std::endl;


		// 1+. 压缩：It is worth noting that most of our wrappers to external libraries requires compressed matrices as inputs

		//			isCompressed()方法——查看是否是压缩状态；插入元素之后若不手动压缩则是非压缩状态；
		std::cout << "sm1.isCompressed() == " << sm1.isCompressed() << std::endl;
		std::cout << "sm2.isCompressed() == " << sm2.isCompressed() << std::endl;
 
		//			 makeCompressed()方法——压缩剩余空间：
		sm2.makeCompressed();				
		std::cout << "sm2.isCompressed() == " << sm2.isCompressed() << std::endl;


		//	 2. 稀疏矩阵的属性

		//			nonZeros()
		std::cout << "sm1非零元素数：" << sm1.nonZeros() << std::endl;

		//			innerSize()
		std::cout << "sm1内维度：" << sm1.innerSize() << std::endl;			// 列优先的稀疏矩阵的列数，行优先的稀疏矩阵的行数
		
		//			outerSize()
		std::cout << "sm1外维度：" << sm1.outerSize() << std::endl;			// 列优先的稀疏矩阵的行数，行优先的稀疏矩阵的列数 


		//		valuePtr()方法——返回稀疏矩阵元素数组的首地址；
		auto dataPtr = sm1.valuePtr();
		std::cout << "valuePtr()方法——返回稀疏矩阵元素数组的首地址\n data of sm1: " << std::endl;
		for(unsigned i = 0; i < sm1.nonZeros(); ++i)
			std::cout << dataPtr[i] << std::endl;
		std::cout << std::endl;

		// 注：貌似直接修改dataPtr指向的数据是危险的。如果要将某元素清零应该使用prune()方法；
		Eigen::SparseMatrix<float> sm11 = sm1;
		dispMat(sm1.toDense());
		std::cout << "sm1.nonZeros() == " << sm1.nonZeros() << std::endl;
		std::cout << "sm1.isCompressed() == " << sm1.isCompressed() << std::endl << std::endl;;

		//		prune方法()——将不满足条件的元素清零；
		sm1.prune([&](const Eigen::Index& row, const Eigen::Index& col, const float& value)->bool
			{
				if (value > 1.3)
					return true;					// 返回true的元素被保留；
				else
					return false;
			});
		std::cout << "prune方法()——将不满足条件的元素清零；" << std::endl;
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


		//	3. 访问稀疏矩阵的元素

		//					coeffRef(row, col)——访问稀疏矩阵某一个元素
		std::cout << " sm2.coeffRef(3, 4) == " << sm2.coeffRef(3, 4) << std::endl;
		std::cout << " sm2.coeffRef(3, 1) == " << sm2.coeffRef(3, 1) << std::endl;

		//					使用InnerIterator访问稀疏矩阵的非零元素
		dispMat(sm2.toDense());
		std::cout << "sm2非零元素数：" << sm2.nonZeros() << std::endl;
		int outerSize2 = sm2.outerSize();				// 默认的列优先存储矩阵，outerSize即列数；
		for (int colIdx = 0; colIdx < outerSize2; ++colIdx)
		{
			for (Eigen::SparseMatrix<float>::InnerIterator it(sm2, colIdx); it; ++it)			// 列优先存储时，InnerIterator即列内迭代器；
			{
				std::cout << "value == " << it.value() << std::endl;					// 元素的值；
				std::cout << "valueRef == " << it.valueRef() << std::endl;			// 元素的引用；
				std::cout << "row == " << it.row() << std::endl;			 // row index
				std::cout << "col == " << it.col() << std::endl;				 // col index (here it is equal to k)
				std::cout << "index == " << it.index() << std::endl;
				std::cout << std::endl;
				//std::cout << "" << it.index() << std::endl;				// inner index, here it is equal to it.row()
			}
		}

		// 4.  稀疏矩阵和稠密矩阵的转换：
		Eigen::MatrixXf m1(3, 4);
		m1 << 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0;

		//			sparseView()
		sm1 = m1.sparseView();
		std::cout << "sm1 == \n" << sm1 << std::endl;
		std::cout << "sm1是否是压缩格式：" << sm1.isCompressed() << std::endl;

		//			toDense();
		Eigen::MatrixXf m2 = sm1.toDense();
		std::cout << "m2 == \n" << m2 << std::endl;
	}


	// 稀疏矩阵的基本变换和操作
	void test1()
	{
		// 1. 求和——没有dense matrix里的.colwise().sum()和.rowwise().sum()操作，可以使用左/右乘向量来实现：
		Eigen::MatrixXf m1 = Eigen::MatrixXf::Random(4, 5);
		Eigen::MatrixXf m2(4, 5);
		m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20;
		std::cout << "m1 == \n" << m1 << std::endl;

		// 2. sparseView()方法；
		Eigen::SparseMatrix<float> sm1 = m1.sparseView();
		Eigen::SparseMatrix<float> sm2 = m2.sparseView();
		Eigen::RowVectorXf colSum = Eigen::RowVectorXf::Ones(4) * sm1;			// 左乘一个全1行向量，得到按列求和的结果。
		Eigen::VectorXf rowSum = sm1 * Eigen::VectorXf::Ones(5);							// 右乘一个全1列向量，得到按行求和的结果。
		std::cout << "按列求和结果：\n" << colSum << std::endl;
		std::cout << "按行求和结果：\n" << rowSum << std::endl << std::endl;

		// 3. 列优先存储稀疏矩阵的row()方法返回的引用是只读的，不可以被赋值。
		std::cout << "遍历第i列：" << std::endl;
		for (Eigen::SparseMatrix<float>::InnerIterator it(sm1, 1); it; ++it)				// 内维度就是存储优先的维度，如这里使用默认的列优先存储，这里的it就是某一列的迭代器；
			std::cout << it.value() << ", ";
		std::cout << std::endl;

		Eigen::RowVectorXf tempVec = sm1.row(1);
		std::cout << "tempVec == \n" << tempVec << std::endl;

		sm1.col(0) = Eigen::SparseVector<float>(4);									// 某一列可以用稀疏向量赋值，不可以用dense向量赋值。
		std::cout << "sm1 == \n" << sm1 << std::endl;
		sm1.makeCompressed();
		std::cout << "sm1.nonZeros()  == " << sm1.nonZeros() << std::endl;

		sm1.col(1) = sm2.col(1);
		std::cout << "sm1 == \n" << sm1 << std::endl;
		sm1.makeCompressed();
		std::cout << "sm1.nonZeros()  == " << sm1.nonZeros() << std::endl; 

		// 使用sp2的某一行给sp1的某一列赋值时，必须将该行向量转置为列向量，否则会出错，且编译器不报异常；
		sm1.col(0) = sm2.row(0).transpose();
		dispMat(sm1.toDense());

		// 不要用稀疏矩阵自带的tranpose()方法，很垃圾
		auto tmpSm = sm1.transpose();					// 只能赋值给行优先存储的稀疏矩阵；

		std::cout << "finished." << std::endl;
	}


	// 稀疏矩阵的拼接
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
		sm.leftCols(3) = sm1;			// 对于行存储的稀疏矩阵来说，leftCols()和rightCols()返回的列是只读的，topRows()和bottomRows()返回的行是可读写的。
		sm.rightCols(1) = sm2;		// 列存储的稀疏矩阵则反之

		//sm.topRows(1) = sm3;
		//sm.bottomRows(3) = sm4;

		std::cout << sm << std::endl;
	}


	// 稀疏矩阵的分解
	void test3() 
	{
		using spMatD = Eigen::SparseMatrix<double>;
		Eigen::MatrixXd m1(4, 4);
		spMatD sm1;
		// m1 << 1, 2, -3, 4, 5, 6, 7, 8, 9, 10, -11, 12, 13, 14, 15, 0;		// 非正定非对称
		m1 << 1, 2, 99, 4, 5, 6, 7, 8, 9, 10, 88, 77, 13, 14, 15, 66; 
		m1 = (m1 * m1.transpose()/1000).eval();					// 对称正定
		debugDisp("m1 == \n", m1, "\n");

		sm1 = m1.sparseView();
		debugDisp("m1.determinant() == ", m1.determinant());

		// 1. 稀疏矩阵的LU分解器：
		Eigen::SparseLU<spMatD> solverLU;
		solverLU.compute(sm1);
		if (solverLU.info() != Eigen::Success)
		{
			debugDisp("LU decomposition failed");
			return;
		}
		debugDisp("solverLU.determinant() == ", solverLU.determinant());		// 稀疏矩阵类不包含求行列式的方法，需要通过LU或LDLT分解器求得。 
 
		// 2. Eigen::SparseQR< MatrixType_, OrderingType_ >类——稀疏矩阵的QR分解器
		/*
			OrderingType_是稀疏矩阵的列排序算法对象
					对稀疏矩阵的列重新排序可以优化矩阵分解、迭代法求线性方程组的性能。
					重新排序可以使得矩阵在因子分解过程中生成的填充量最小化，从而减少计算量和内存占用。
					一些求解算法对列向量的排序非常敏感。通过选择合适的排序策略，我们可以减少求解步骤的数量
					某些列排序算法可以提高求解算法的数值稳定性，减少数值误差的影响。

			OrderingType_有三种可用的排序算法：
					1. Eigen::NaturalOrdering< StorageIndex >
								自然排序，不重新排序，使用原矩阵；
					2. Eigen::AMDOrdering< StorageIndex >
								应用于对称矩阵。
								近似最小度排序（Approximate Minimum Degree Ordering），通过考虑每次选取最小度的列进行排序，以减少填充。
								列的度——该列的非零元素的数量。
					3. Eigen::COLAMDOrdering< StorageIndex >
								可应用于非对称矩阵的近似最小度排序。
								可能比AMDOrdering计算成本更高。
				
		*/
		Eigen::SparseQR<spMatD, Eigen::COLAMDOrdering<int>> solverQR;
		solverQR.compute(sm1);
		if (solverQR.info() != Eigen::Success)
		{
			debugDisp("QR decomposition failed");
			return;
		}
		debugDisp("solverQR.rank() == ", solverQR.rank());			// 稀疏矩阵类不包含求秩的方法，需要通过QR求得。 

		// 3. 稀疏矩阵的LDLT分解器（LLT分解器很挫的样子，所以不测了）
		Eigen::SimplicialLDLT<spMatD> solverLDLT;
		solverLDLT.compute(sm1);
		if (solverLDLT.info() != Eigen::Success)		
		{
			// ！！！理论上只有正定对称矩阵才可以做LDLT分解，但如果矩阵不满足要求时这里也不会分解失败，但计算出的行列式是错误的。
			debugDisp("LDLT decomposition failed");
			return;
		}
		debugDisp("solverLDLT.determinant() == ", solverLDLT.determinant());	// 稀疏矩阵类不包含求行列式的方法，需要通过LU或LDLT分解器求得。 

		debugDisp("test3 finished.");
	}


	// 求解稀疏线性方程组——直接法
	void test4()
	{
		// 1. 直接法求解稀疏线性方程组

		//				稀疏矩阵表示的线性方程组：Ax = b
		Eigen::SparseMatrix<float> A;
		Eigen::VectorXf b, x;
		// Eigen::SimplicialLLT<Eigen::SparseMatrix<float>> solver;					// 基于LLT分解的求解器，一般不如LDLT分解
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;					// 基于LDLT分解的求解器，使用于稠密度很低且规模不是很大的矩阵；？？？貌似需要对称矩阵？
		// Eigen::SparseLU<Eigen::SparseMatrix<float> > solver;						// 基于LU分解的求解器。
		// Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int> > solver;		// 基于QR分解的求解器，推荐最小二乘问题使用。

		//				写入方程组数据：
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


	// 求解稀疏线性方程组——迭代法
	void test5()
	{
		//	1. 生成稀疏线性方程组：Ax = b
		Eigen::SparseMatrix<float> A;
		Eigen::VectorXf b(2), x; 
		Eigen::MatrixXf Adense(2, 2);
		Adense << 1, 2, 1, -1;
		b << 0, 3;
		A = Adense.sparseView();
		debugDisp("A == \n", A, "\n");
		debugDisp("b == \n", b, "\n"); 
		debugDisp("det(A) == ", Adense.determinant());

		// 2. 迭代法求解稀疏线性方程组：
		/*
			迭代器求解器：
					ConjugateGradient				
							经典迭代共轭梯度求解器。 
							要求矩阵为SPD矩阵（Symmetric positive definite matrices 对称、正定）

					LeastSquaresConjugateGradient
							最小二乘共轭梯度求解器
							适用于长方形矩阵（非方阵）

					BiCGSTAB						
							迭代双共轭梯度求解器（Iterative stabilized bi-conjugate gradient）			
							要求矩阵为方阵。
		*/
		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float>> solverLSCG;		
		Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> solverBCG;

		//			2.1 最小二乘共轭梯度求解器
		solverLSCG.setTolerance(0.00001);				//	设置迭代精度
		solverLSCG.compute(A);
		if (solverLSCG.info() != Eigen::Success)
		{
			debugDisp("!!! Eigen::LeastSquaresConjugateGradient<> computing failed.");
			return;
		}
		x = solverLSCG.solve(b); 
		debugDisp("Eigen::LeastSquaresConjugateGradient<> result： x == \n", x, "\n\n");
		debugDisp("solverLSCG.iterations() == ", solverLSCG.iterations());					// 迭代次数
		debugDisp("solverLSCG.error() == ", solverLSCG.error());

		//			2.2 迭代双共轭梯度求解器
		solverBCG.compute(A);
		if(solverBCG.info() != Eigen::Success)
		{
			debugDisp("!!! Eigen::BiCGSTAB<> computing failed.");
			return;
		}
		x = solverBCG.solve(b);
		debugDisp("Eigen::BiCGSTAB<> result： x == \n", x, "\n");
		debugDisp("solverBCG.iterations() == ", solverBCG.iterations());
		debugDisp("solverBCG.error() == ", solverBCG.error());

		debugDisp("\n\ntest4() finished.");
	}


	// 比较不同求解系数线性方程组方法的速度：
	void test6() 
	{
		// 构造一个3n x 3n的A矩阵，其中
		tiktok& tt = tiktok::getInstance();
		const double tolerance = 1e-3;			// 设定迭代法的精度；
		double interval = 0;								// 记录时间间隔；

		// 1. 构造系数线性方程组的A矩阵和b向量；
		const int n = 10000;
		const int rows = n * 3;
		const int cols = n * 3;
		int offset = 0;
		std::default_random_engine e;											// 随机数生成器的引擎对象
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
		tt.endCout("1. 生成稀疏矩阵A耗时：");
		debugDisp("\n");

		// 2. QR分解法
		debugDisp("QR分解法：");
		{
			Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solverQR;
			interval = 0;

			//			2.1 QR分解：
			tt.start();
			solverQR.compute(A);
			if (solverQR.info() != Eigen::Success)
			{
				std::cout << "QR decomposition failed" << std::endl;
				return;
			}
			tt.endCout("A矩阵QR分解耗时：");
			interval += tt.lastDur;

			debugDisp("solverQR.rank() == ", solverQR.rank());
			if (solverQR.rank() < 3 * n)
			{
				debugDisp("！！！A矩阵不满秩。");
				return;
			}

			//			2.2 解线性方程组：
			tt.start();
			Eigen::VectorXd x = solverQR.solve(b);
			tt.endCout("解稀疏线性方程组耗时：");
			interval += tt.lastDur;

			//			2.3 验证结果准确性：
			Eigen::VectorXd b_prime = A * x;
			Eigen::VectorXd residual = b - b_prime;
			double sqrErr = residual.squaredNorm();				// 平方误差；
			debugDisp("平方误差sqrErr == ", sqrErr);

			debugDisp("QR分解法总耗时：", interval, "s\n\n");
		}

		// 3. LU分解法：
		debugDisp("LU分解法：");
		{
			Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solverLU;
			interval = 0;

			//			3.1 LU分解：
			tt.start();
			solverLU.compute(A);
			if (solverLU.info() != Eigen::Success)
			{
				std::cout << "LU decomposition failed" << std::endl;
				return;
			}
			tt.endCout("A矩阵LU分解耗时：");
			interval += tt.lastDur; 

			//			3.2 解线性方程组：
			tt.start();
			Eigen::VectorXd x = solverLU.solve(b);
			tt.endCout("解稀疏线性方程组耗时：");
			interval += tt.lastDur;

			//			3.3 验证结果准确性：
			Eigen::VectorXd b_prime = A * x;
			Eigen::VectorXd residual = b - b_prime;
			double sqrErr = residual.squaredNorm();				// 平方误差；
			debugDisp("平方误差sqrErr == ", sqrErr);

			debugDisp("LU分解法总耗时：", interval, "s\n\n");
		}

		// 4. 最小二乘共轭梯度迭代法：
		debugDisp("最小二乘共轭梯度迭代法：");
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
			debugDisp("solverLSCG.iterations() == ", solverLSCG.iterations());							// 迭代次数
			debugDisp("solverLSCG.error() == ", solverLSCG.error());			

			// 4.3 验证结果准确性：
			Eigen::VectorXd b_prime = A * x;
			Eigen::VectorXd residual = b - b_prime;
			double sqrErr = residual.squaredNorm();				// 平方误差；
			debugDisp("平方误差sqrErr == ", sqrErr);

			debugDisp("最小二乘共轭梯度迭代法总耗时：", interval, "s\n\n");
		}

		// 5. 双共轭梯度迭代法：
		debugDisp("双共轭梯度迭代法：");
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
			debugDisp("solverBCG.iterations() == ", solverBCG.iterations());							// 迭代次数
			debugDisp("solverBCG.error() == ", solverBCG.error());

			// 5.3 验证结果准确性：
			Eigen::VectorXd b_prime = A * x;
			Eigen::VectorXd residual = b - b_prime;
			double sqrErr = residual.squaredNorm();				// 平方误差；
			debugDisp("平方误差sqrErr == ", sqrErr);

			debugDisp("双共轭梯度迭代法总耗时：", interval, "s\n\n");
		}


		debugDisp("\n\ntest6() finished.");
	}
}

