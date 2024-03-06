#include "test_dense_mat.h"


#define MAXLEN 1024

using namespace Eigen;
using spMatf = Eigen::SparseMatrix<float, ColMajor>;
using TripF = Eigen::Triplet<float>;


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


////////////////////////////////////////////////////////////////////////////////////////////// 测试eigen中的稠密矩阵：
namespace TEST_DENSE_MAT
{
	// 暂时无法分类：
	void test000()
	{
		// 1. 矩阵索引：
		std::vector<int> rowVec{ 1, 3, 4 };
		std::vector<int> colVec{ 5, 2, 0 };

		std::vector<int> ind{ 4,2,5,5,3 };
		Eigen::MatrixXi A = Eigen::MatrixXi::Random(4, 6);
		std::cout << "Initial matrix A:\n" << A << "\n\n";

		// 2. 查看eigen库的版本：
		std::cout << EIGEN_WORLD_VERSION << std::endl;
		std::cout << EIGEN_MAJOR_VERSION << std::endl;
		std::cout << EIGEN_MINOR_VERSION << std::endl;		// 版本为3.4.0

		// 3. 
		{
			Eigen::MatrixXi m1(3, 4);
			Eigen::Matrix2i m2;
			Eigen::VectorXi v1(5);
			debugDisp("m1.RowsAtCompileTime == ", m1.RowsAtCompileTime);
			debugDisp("m1.ColsAtCompileTime == ", m1.ColsAtCompileTime);
			debugDisp("m2.RowsAtCompileTime == ", m2.RowsAtCompileTime);
			debugDisp("m2.ColsAtCompileTime == ", m2.ColsAtCompileTime);
			debugDisp("v1.RowsAtCompileTime == ", v1.RowsAtCompileTime);
			debugDisp("v1.ColsAtCompileTime == ", v1.ColsAtCompileTime);
		}
	}


	// test0——eigen库的基本数据结构
	void test0()
	{
		// 堆矩阵、向量——确定了尺寸，但未初始化,数据存在堆上
		/*
			基本模板——Matrix<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
			堆矩阵——typedef Matrix<double, Dynamic, Dynamic> Eigen::MatrixXd;
			堆向量——typedef Matrix<int, Dynamic, 1> Eigen::VectorXi;

			继承关系： 
				Matrix<> → PlainObjectBase<>
				MatrixBase<> → DenseBase<>

			Options_:
					Eigen::RowMajor		行优先存储；
					ColMajor			列优先存储（默认值）；
					AutoAlign
					DontAlign. 
		*/

		Eigen::MatrixXd m1(2, 2);
		Eigen::MatrixXf mf1(1, 2);
		Eigen::VectorXd v1(3);			// 注意是列向量
		std::cout << v1 << std::endl;

		// 1. 堆Array
		ArrayXXd a1(2, 2), a2(2, 2);
		a1 << 1, 2, 3, 4;											// operator << 填充元素是行优先的顺序；
		a2 << 1, 2, 3, 4;
		std::cout << "a1 = \n" << a1 << std::endl;
		std::cout << "a1*a2 = \n" << a1 * a2 << std::endl;

		// 2. 生成特殊向量的接口——LinSpaced(元素数，起点，终点)
		int start = 0;
		int end = 10;
		Eigen::VectorXi vi1 = Eigen::VectorXi::LinSpaced(end - start + 1, start, end);
		std::cout << "vi1 == " << vi1 << std::endl;

		Eigen::VectorXf vf1 = Eigen::VectorXf::LinSpaced(5, 0, 10.0);
		std::cout << "vf1 == " << vf1 << std::endl;

		// 3. 生成特殊矩阵的接口——Random(), Constant(), Ones()....
		Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(3, 3);              // 矩阵类static方法——Random()——返回随机矩阵
		Eigen::MatrixXd m3 = Eigen::MatrixXd::Constant(3, 3, 1.2);		// 常数矩阵，前面是尺寸，后面是值；
		Eigen::MatrixXd m4 = Eigen::MatrixXd::Ones(1, 2);					// 全1矩阵
 
		// 4. 数据存在栈上的矩阵类型
		Matrix3d mm1 = Matrix3d::Random();
		Vector3d vv1(1, 2, 3);
		std::cout << "m2 = \n" << m2 << std::endl << std::endl;
		std::cout << "mm1 = \n" << mm1 << std::endl << std::endl;

		mm1 = m2;				//		堆矩阵和栈矩阵可以相互赋值。
		std::cout << "mm1 = \n" << mm1 << std::endl << std::endl;


		// 初始化：

		//			栈向量可以使用大括号初始化
		Eigen::Vector3f vvv1{ 1,2,3 };
		RowVector3f vvv2{ 4,5,6 };

		std::cout << vvv1 << std::endl;
		std::cout << vvv2 << std::endl;

		//	赋值

		//			列向量对象可以和矩阵对象相互构造
		Vector3d vv2(1, 2, 3);

		//			列向量对象可以和矩阵对象相互赋值
		Eigen::MatrixXd mm(v1);
		vv2 = m2.block<3, 1>(0, 0);
		std::cout << "vv2 = \n" << vv2 << std::endl << std::endl;
		mm = vv2;
		std::cout << "mm = \n" << mm << std::endl << std::endl;
	}


	// test1——矩阵性质、元素访问。
	void test1()
	{
		Eigen::MatrixXd m1(3, 4);
		Eigen::VectorXd v1(5);

		//1.  输出流运算符赋值——行优先顺序填充
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		v1 << 1, 2;
		std::cout << "m1 = \n" << m1 << std::endl << std::endl;
		std::cout << "v1 = \n" << v1 << std::endl << std::endl;

		v1 << 3, 4, 5;			// 不会接着上面的赋值，而是从第一个元素开始赋值
		std::cout << "v1 = \n" << v1 << std::endl << std::endl;

		// 2. 矩阵元素的索引：

		//		下标运算符[]只能获取向量元素，矩阵对象无法使用，因为[]只支持一个参数。
		std::cout << "v1[0] == " << v1[0] << std::endl << std::endl;

		//		括号运算符访问元素，注意索引从0开始
		std::cout << "m1(0, 1) ==" << m1(0, 1) << std::endl;
		std::cout << "v1(3) == " << v1(3) << std::endl;


		// 3. 求矩阵的性质的类内接口
		m1 = Eigen::MatrixXd::Random(3, 4);
		Eigen::Matrix<double, Dynamic, Dynamic, Eigen::RowMajor> m1_rm(5, 6);		// 行优先存储；
		std::cout << "m1 == \n" << m1 << std::endl;
		std::cout << "元素数：m1.size() == " << m1.size() << std::endl;
		std::cout << "行数：m1.rows() == " << m1.rows() << std::endl;
		std::cout << "列数：m1.cols() == " << m1.cols() << std::endl;
		std::cout << "求和：sum():       " << m1.sum() << std::endl;
		std::cout << "连乘：prod():      " << m1.prod() << std::endl;
		std::cout << "均值：mean():      " << m1.mean() << std::endl;
		std::cout << "最小元素：minCoeff():  " << m1.minCoeff() << std::endl;
		std::cout << "最大元素：maxCoeff():  " << m1.maxCoeff() << std::endl;
		std::cout << "矩阵的迹：trace():     " << m1.trace() << std::endl << std::endl;
		std::cout << "行列式：m1.determinant() == " << m1.determinant() << std::endl;

		//			outerSize(), innerSize()——默认的列优先存储的矩阵，outerSize为列数，行优先存储的outerSize为行数；
		std::cout << "m1.outerSize() == " << m1.outerSize() << std::endl;			
		std::cout << "m1.innerSize() == " << m1.innerSize() << std::endl;
		std::cout << "m1_rm.outerSize() == " << m1_rm.outerSize() << std::endl;
		std::cout << "m1_rm.innerSize() == " << m1_rm.innerSize() << std::endl;
 
		//			逆矩阵——使用lu分解得到
		std::cout << "逆矩阵：m1.inverse() ==  \n" << m1.inverse() << std::endl;

		//			转置
		std::cout << "矩阵的转置：transpose() \n" << m1.transpose() << std::endl << std::endl;		// 返回转置的矩阵，矩阵自身不转置
		std::cout << m1 << std::endl;

		//			bug——transpose()在对自身赋值的时候有时会有bug，要想要矩阵自身改变，使用~InPlace()
		Eigen::Matrix3f mm;
		mm << 1, 2, 3, 4, 5, 6, 0, 8, 9;
		std::cout << mm << std::endl << std::endl;
		//mm = mm.transpose();
		mm.transposeInPlace();
		std::cout << mm << std::endl << std::endl;
		Eigen::Matrix3f mmCopy = mm;
		Eigen::Matrix3f mm1 = mm.inverse();			// inverse()也一样，但是只能创建其他变量来赋值。
		std::cout << mm1 << std::endl << std::endl;
		std::cout << mm1 * mmCopy << std::endl << std::endl;

		//			norm()——欧几里得范数，也是p==2时的lp范数。既所有元素平方和的开方。
		debugDisp("v1 == \n", v1, "\n");
		debugDisp("v1.norm() == ", v1.norm()); 		// 向量的范数等于向量的模长。
		m1.resize(3, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		std::cout << "m1 == \n" << m1 << std::endl;
		std::cout << "m1.norm() == " << m1.norm() << std::endl;
		std::cout << "m1.rowwise().norm() == \n" << m1.rowwise().norm() << std::endl;		// 所有行向量求norm，返回一个列向量；
		std::cout << "m1.colwise().norm() == \n" << m1.colwise().norm() << std::endl;

		//		通过范数求行向量的模长：
		Eigen::Matrix3f arrows;
		arrows << 0, 1, 0, 0, 3, 4, 0, 5, 12;
		Eigen::Vector3f arrowsLen = arrows.rowwise().norm();
		dispVec<float>(arrowsLen);

		// data()——获得矩阵数据数组的首指针

		//				快速遍历矩阵中的元素：
		Eigen::MatrixXf m2(5, 10);
		float* elemPtr = nullptr;
		for (unsigned i = 0; i < m2.size(); ++i)
		{
			elemPtr = m2.data() + i;
			*elemPtr = i;
		}
		std::cout << "m2 == \n" << m2 << std::endl << std::endl;;

		//				拷贝堆矩阵中的数据：
		Eigen::MatrixXf m22(5, 10);
		std::memcpy(m22.data(), m2.data(), m2.size() * sizeof(float));
		std::cout << m22 << std::endl;

		//	 minCoeff() —— 搜索矩阵元素中的最值
		std::cout << "m1 == \n" << m1 << std::endl;
		Eigen::MatrixXd::Index maxRow, maxCol;
		Eigen::MatrixXd::Index minRow, minCol;
		double min = m1.minCoeff(&minRow, &minCol);
		double max = m1.maxCoeff(&maxRow, &maxCol);
		std::cout << "最小元素行号" << minRow << "，列号" << minCol << std::endl;
		std::cout << "最大元素行号" << maxRow << "，列号" << maxCol << std::endl;

		v1.resize(10);
		v1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
		Eigen::VectorXd::Index minIdx, maxIdx;
		min = v1.minCoeff(&minIdx);
		debugDisp("向量中最小元素：min == " , min , ", 索引minIdx == ：" , minIdx); 

		v1 <<  -INFINITY, 2, 3, 4, 5, 6, 7, 8, INFINITY, 0;
		min = v1.minCoeff(&minIdx);					// 元素中可以有INFINITY，但是不可以有NAN，否则运行时会抛出异常；
		max = v1.maxCoeff(&maxIdx);
		debugDisp("向量中最小元素：min == ", min, ", 索引minIdx == ：", minIdx);
		debugDisp("向量中最大元素：max == ", max, ", 索引maxIdx == ：", maxIdx);

		// all(), any(), count()——矩阵元素的搜索
		ArrayXXf a1(3, 3);
		ArrayXXi a2(2, 2);
		a1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		a2 << 1, 2, 3, 4;
		bool flag1 = (a1 > 0).all();						// all()——如果所有元素都为true，则返回true,
		bool flag2 = (a2 < 3).any();					// any()——如果至少有一个元素为true，则返回true
		Eigen::Index count = (a1 > 5).count();		// count()——返回满足条件的元素数

		Eigen::VectorXi vec(9);
		vec << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		ArrayXi aa = vec.array();
		flag1 = (aa > 2).any();
		std::cout << "flag1" << flag1 << std::endl;

		// 貌似只能对单个元素使用，不能对一组元素使用。
		m1.resize(3, 3);
		m1 << 1, 2, 3, 1, 2, 3, 4, 5, 6;
		Eigen::RowVector3d vd1(1, 2, 3);
		std::cout << "m1中是否有行向量(1,2,3)：" << (m1.array() == vd1.array()).any() << std::endl;
		m1 = Eigen::MatrixXd::Ones(3, 3);
		std::cout << "m1中是否有行向量(1,2,3)：" << (m1.array() == vd1.array()).any() << std::endl;
	}


	// test2——矩阵基本变换、运算
	void test2()
	{
		Eigen::MatrixXf m1(3, 4), m3(4, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

		std::cout << "m1 = \n" << std::endl;
		std::cout << m1 << std::endl << std::endl;

		// resize()——更改矩阵尺寸——只适用于堆矩阵、堆向量  ！！！不管矩阵尺寸更改后如何，所有元素全部变成随机值。
		m3.resize(3, 3);
		m3.row(2) = Eigen::VectorXf::Ones(3);
		std::cout << "m3.size() == " << m3.size() << std::endl;
		std::cout << "m3.rows() == " << m3.rows() << std::endl;
		std::cout << "m3.cols() == " << m3.cols() << std::endl;
		std::cout << "m3 = \n" << m3 << std::endl << std::endl;

		// conservativeResize() 更改矩阵尺寸，保留原有数据
		m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m3.conservativeResize(4, 3);
		m3.row(3) = 5 * Eigen::VectorXf::Ones(3);
		std::cout << "m3 = \n" << m3 << std::endl << std::endl;

		// colwise(), rowwise()——矩阵按行、按列操作
		std::cout << m1.colwise().sum() << std::endl << std::endl;				// 按列求和，压成一个行向量
		std::cout << m1.rowwise().sum() << std::endl << std::endl;				// 按行求和，压成一个列向量。

		// 矩阵扩张——通过矩阵乘法，或克罗内克积实现；类似于matlab中的repmat();
		std::cout << "矩阵扩张：" << std::endl;
		RowVector3f a(1,2,3);
		Eigen::Vector3f b(4, 5, 6);
		
		auto aa = Eigen::MatrixXf::Ones(3, 1) * a;			//		行向量扩张N行 == 左乘ones(N,1)
		auto bb = b * Eigen::MatrixXf::Ones(1, 4);			//		列向量扩张为N列 == 右乘ones(1,N)
		std::cout << aa << std::endl << std::endl;;
		std::cout << bb << std::endl << std::endl;

		Eigen::MatrixXd aaa, bbbb;
		aaa = kron(RowVector3f::Ones(), a); 
		bbbb = kron(Matrix2i::Ones(), b); 
		std::cout << aaa << std::endl << std::endl;;
		std::cout << bbbb << std::endl << std::endl;

		// reverse();
		m1.resize(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		std::cout << "m1 = \n" << m1 << std::endl << std::endl;
		Eigen::MatrixXf result1 = m1.reverse();								// 行列都倒装
		Eigen::MatrixXf result2 = m1.colwise().reverse();				// 矩阵视为列向量，然后每一行内部保持不变倒装；
		Eigen::MatrixXf result3 = m1.rowwise().reverse();				// 矩阵视为行向量....
		std::cout << "result1 = \n" << result1 << std::endl << std::endl;
		std::cout << "result2 = \n" << result2 << std::endl << std::endl;
		std::cout << "result3 = \n" << result3 << std::endl << std::endl;
		m1.colwise().reverseInPlace();
		std::cout << "m1 = \n" << m1 << std::endl << std::endl;

		// 向量点乘、叉乘
		Eigen::Vector3f v1(1, 2, 3);
		Eigen::Vector3f v2(3, 2, 1);
		std::cout << "v1.dot(v2) == \n" << v1.dot(v2) << std::endl;
		std::cout << "v1.cross(v2) == " << v1.cross(v2) << std::endl << std::endl;;

		// 向量归一化
		v1.normalize();					// normalize()会将向量本身归一化，normalized()只是返回归一化向量，不改变自身
		std::cout << "v1 = \n" << v1 << std::endl << std::endl;
		std::cout << "v2.normalized() = \n" << v2.normalized() << std::endl << std::endl;
		std::cout << "v2 = \n" << v2 << std::endl << std::endl;

		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		Eigen::MatrixXf m11 = m1.rowwise().normalized();			// 每条行向量归一化；
		Eigen::MatrixXf m12 = m1.colwise().normalized();				// 每条列向量归一化；
		std::cout << "m1 = \n" << m1 << std::endl << std::endl;
		std::cout << "m11 = \n" << m11 << std::endl << std::endl;
		std::cout << "m12 = \n" << m12 << std::endl << std::endl;
		
	}


	// test3——稠密矩阵的分解、求稠密线性方程组：
	void test3()
	{
		// 几种分解方法求线性方程组的性能比较：
		/*
			LU分解：		
					效率高；
					适用于稠密度高的矩阵；

			LDLT分解：	
					适用于堆成正定矩阵；

			QR分解：		
					在处理奇异矩阵时数值稳定性好；
					在求逆的时效率更高；
		
			奇异值分解（SVD）
					若矩阵是稀疏的、奇异值分解的低秩近似能够提供足够的精度，
							并且矩阵的特征不适合其他分解方法，那么奇异值分解可能是一个合适的选择
		*/
		Eigen::MatrixXd A(3, 3);
		Vector3d b{1, 2, 3};
		Vector3d x;
		A << 1, 22, 3, 4, 55, 6, 77, 8, 9;
		debugDisp("A == \n", A, "\n"); 

		// 1. 特征分解：
		EigenSolver<Matrix3d> solverEigen(A);
		{
			Matrix3d D = solverEigen.pseudoEigenvalueMatrix();			// 对角线元素是特征值
			Matrix3d V = solverEigen.pseudoEigenvectors();					// 每一个列向量都是特征向量。
			debugDisp("特征值矩阵D == \n", D, "\n");
			debugDisp("特征值矩阵V == \n", V, "\n");
			debugDisp("V * D * V.inverse() == \n", V * D * V.inverse(), "\n");
			debugDisp("A * V.block<3, 1>(0, 0) == \n", A * V.block<3, 1>(0, 0), "\n");
			debugDisp("D(0, 0) * V.block<3, 1>(0, 0) == \n", D(0, 0) * V.block<3, 1>(0, 0), "\n");
			debugDisp("\n\n");
		}

		// 2. 矩阵的LU分解——Eigen::FullPivLU
		Eigen::FullPivLU<Eigen::MatrixXd> solverLU;
		{
			Eigen::MatrixXd L, U, P_inv, Q_inv;
			solverLU.compute(A);
			debugDisp("执行LU分解：solverLU.matrixLU() == \n", solverLU.matrixLU());
			L = Eigen::MatrixXd::Identity(3, 3);
			L.triangularView<StrictlyLower>() = solverLU.matrixLU();
			U = solverLU.matrixLU().triangularView<Upper>();
			P_inv = solverLU.permutationP().inverse();
			Q_inv = solverLU.permutationQ().inverse();
			debugDisp("P_inv == \n", P_inv, "\n");
			debugDisp("L == \n", L, "\n");
			debugDisp("U == \n", U, "\n");
			debugDisp("Q_inv == \n", Q_inv, "\n");
			debugDisp("P_inv * L * U * Q_inv == \n", P_inv * L * U * Q_inv, "\n");

			// LU分解解线性方程组；
			x = solverLU.solve(b);
			debugDisp("线性方程组Ax == b的解：x == \n", x, "\n");
			debugDisp("x == \n", x, "\n");
			debugDisp("A * x == \n", A * x, "\n");

			debugDisp("\n\n");
		}

		// 3. 矩阵的奇异值分解（SVD）——  
		JacobiSVD<Eigen::MatrixXd> solverSVD(A, ComputeThinU | ComputeThinV);
		{
			debugDisp("奇异值：solverSVD.singularValues() == \n ", solverSVD.singularValues(), "\n");
			debugDisp("solverSVD.matrixU() == \n ", solverSVD.matrixU(), "\n");
			debugDisp("solverSVD.matrixV() == \n ", solverSVD.matrixV(), "\n");
			debugDisp("\n");

			//		3.1 通过SVD求解线性方程组： 
			x = solverSVD.solve(b);
			debugDisp("使用奇异值分解解线性方程组：b == \n", b, "\n");
			debugDisp("A * x == \n", A * x, "\n");
			debugDisp("\n");

			//		3.2 通过SVD求矩阵的条件数：
			double condNum = calcCondNum(A);
			debugDisp("矩阵条件数：condNum == ", condNum);
			debugDisp("\n\n");

			//		3.3 求秩——矩阵类中没有求秩的方法，只能通过SVD分解来求秩。
			debugDisp("solverSVD.rank() == ", solverSVD.rank());  
		}

		// 4. LLT分解、LDLT分解：
		Eigen::LLT<Eigen::MatrixXd> solverLLT;
		Eigen::LDLT<Eigen::MatrixXd> solverLDLT;
		{
			Eigen::MatrixXd A1(2, 2); 
			A1 << 2, -1, -1, 3;							// ？？？不知道为什么非正定的矩阵也可以做LLT分解；
			debugDisp("A1 == \n", A1, "\n");

			solverLLT.compute(A1);
			solverLDLT.compute(A1);
			if (solverLLT.info() != Eigen::Success)
			{
				debugDisp("LLT decomposition failed");
				return;
			}
			if (solverLDLT.info() != Eigen::Success)
			{
				debugDisp("LDLT decomposition failed");
				return;
			}
			Eigen::MatrixXd L_llt = solverLLT.matrixL();
			debugDisp("LLT分解：L_llt == \n", L_llt, "\n");
		}

		// 5. QR分解——基于householder变换： 
		Eigen::HouseholderQR<Eigen::MatrixXd> solverQR1;				// 没有rank()方法；稳定性依赖于矩阵条件数；
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solverQR2;		// 有rank()方法；稳定性不错；
		{
			//	5.1 测试自己封装的使用QR分解的求解恰定稠密线性方程组的轮子： 
			Eigen::MatrixXd B = Eigen::MatrixXd::Random(3, 3);
			Eigen::MatrixXd X;
			solveLinearEquation(x, A, b);
			solveLinearEquations<double>(X, A, B);
			debugDisp("线性方程组Ax == b的解：x == \n", x, "\n");
			debugDisp("x == \n", x, "\n");
			debugDisp("A * x == \n", A * x, "\n");
			debugDisp("B == \n", B, "\n");
			debugDisp("线性方程组AX == B的解：X == \n", X, "\n");
			debugDisp("A * X == \n", A * X, "\n");

			// 5. 2 Eigen::HouseholderQR<>
			A.resize(2, 3);
			A << 1, 2, 3, 4, 5, 6;
			solverQR1.compute(A);			// A = Q* R;
			solverQR2.compute(A);			// A * P = Q * R;	其中P是置换矩阵；
			if (solverQR2.info() != Eigen::Success)
			{
				debugDisp("ColPivHouseholderQR decomposition failed");
				return;
			}
			Eigen::MatrixXd Q = solverQR1.householderQ();
			Eigen::MatrixXd R = solverQR1.matrixQR().triangularView<Eigen::Upper>();
			debugDisp("\n\nQR分解：Eigen::HouseholderQR<>");
			debugDisp("A == \n", A, "\n");
			debugDisp("Q == \n", Q, "\n");
			debugDisp("R == \n", R, "\n");
			debugDisp("Q * R == \n", Q * R, "\n");

			// 5.3 Eigen::ColPivHouseholderQR<>
			Q = solverQR2.householderQ();
			R = solverQR2.matrixQR().triangularView<Eigen::Upper>();
			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P = solverQR2.colsPermutation();			// 置换矩阵
			Eigen::MatrixXd Pd = P;								// Eigen::PermutationMatrix通过赋值运算符来转换为普通矩阵类型，没有cast方法；
			Eigen::MatrixXd PdT = Pd.transpose();		// 置换矩阵是正交阵，其逆等于其转置；
			debugDisp("\n\nQR分解：Eigen::ColPivHouseholderQR<>");
			debugDisp("Q == \n", Q, "\n");
			debugDisp("R == \n", R, "\n");
			debugDisp("P.rows() == ", P.rows());
			debugDisp("P.cols() == ", P.cols());
			// debugDisp("A * P == \n", A * P, "\n");
			debugDisp("A * Pd == \n", A * Pd, "\n");
			debugDisp("Q * R == \n", Q * R, "\n");
			debugDisp("Q * R * PdT == \n", Q * R * PdT, "\n");
			debugDisp("solverQR2.rank() == ", solverQR2.rank()); 
		}

		 
		debugDisp("test3() finished.");
	}


	// test4——矩阵块操作：
	void test4()
	{
		Eigen::MatrixXf m1(5, 6);
		float* pdata = m1.data();
		for (unsigned i = 0; i < m1.size(); ++i)
			pdata[i] = static_cast<float>(i);
		std::cout << "m1 == \n" << std::endl;
		dispMat(m1);

		// 1. setConstant()方法——块赋值：
		m1.setConstant(1.0);
		m1.topRightCorner(2, 3).setConstant(2.0);
		dispMat(m1);

		//	2. VectorXT::segment<>()方法——提取向量的子向量片段，返回其左值引用。两个重载
		Eigen::VectorXi v1(9);
		v1 << 1, 2, 3, 4, 5, 6, 7, 8, 8;

		//			segment()重载1——segment<elemCount>( startIdx)，静态segment，片段长度elemCount必须在编译期已知，是constExpr
		std::cout << "v1.segment<3>(1)  == \n" << v1.segment<3>(1) << std::endl;

		//			segment()重载2——segment(startIdx, elemCount)，动态segment，片段长度elemCount可以在编译器内未知。
		std::cout << "v1.segment(1, v1(2)) == \n" << v1.segment(1, v1(2)) << std::endl;

		//			segment()返回的向量片段是左值引用：
		v1.segment<5>(2) = 999 * Eigen::VectorXi::Ones(5);
		std::cout << "v1 == \n" << v1 << std::endl << std::endl;

		// 2+. 向量的head(), tail()方法——行为上看是返回向量片段的引用；
		v1.head(3) = Eigen::Vector3i{-1, -1, -1};
		v1.tail(2) = Eigen::RowVector2i{88, 88};
		debugDisp("v1 == \n", v1);

		// 3. block()方法——提取子矩阵块，返回其左值引用。有两个重载
		std::cout << "block()方法" << std::endl;

		//			重载1：block<rows, cols>(startRow, startCol)——静态block，矩阵块的尺寸必须是编译期已知的，必须是constExpr
		std::cout << "m1.block<2, 2>(1, 1)  == \n" << m1.block<2, 2>(1, 1) << std::endl << std::endl;

		//			重载2: block(startRow, startCol, rows, cols)——动态block，矩阵块的尺寸可以是编译期未知的
		std::cout << "m1.block(1, 0, m1(2, 0), m1(3, 0)) == \n" << m1.block(1, 0, m1(2, 0), m1(3, 0)) << std::endl << std::endl;


		//			 返回的矩阵块是左值引用。可以被<<赋值；
		m1.block<1, 4>(1, 0) << 88, 99, 111, 222;
		std::cout << "m1 = \n" << m1 << std::endl << std::endl;
		m1.block<1, 4>(1, 0) << 5, 6, 7, 8;


		//			返回的矩阵块引用可以直接用' = ' 被另一个矩阵赋值。
		Eigen::MatrixXf m2 = Eigen::MatrixXf::Ones(1, 6);
		m1.block<1, 6>(0, 0) = m2;
		std::cout << "m1 == \n" << m1 << std::endl;


		// 4. row(), col()——提取某一行或者某一列
		std::cout << "row(), col()——提取某一行或者某一列\n m1.col(3) == \n" << m1.col(3) << std::endl;

		//			提取出的行、列是原数据的引用，不是新的拷贝，并且是左值。
		m1.col(3) << 1, 2, 3, 4, 5;
		std::cout << "m1.col(3) == \n" << m1.col(3) << std::endl;

		// 5. topRows(), bottomRows(), leftCols(), rightCols();
		std::cout << "m1.leftCols(2) ==   \n" << m1.leftCols(2) << std::endl;
		std::cout << "m1.topRows(2) == \n" << m1.topRows(2) << std::endl;
 
		//	6. minCoeff()求最值对矩阵分块也同样适用
		Eigen::MatrixXf::Index maxRow, maxCol;
		Eigen::MatrixXf::Index minRow, minCol;
		float min = m1.col(1).minCoeff(&minRow, &minCol);
		float max = m1.block<3, 3>(2, 2).maxCoeff(&maxRow, &maxCol);
		std::cout << "m1 == \n" << m1 << std::endl;

		std::cout << "m1.col(1) == \n" << m1.col(1) << std::endl;
		std::cout << "min == " << min << std::endl;
		std::cout << "minRow == " << minRow << std::endl;
		std::cout << "minCol == " << minCol << std::endl;
		std::cout << "max == " << max << std::endl;
		std::cout << "maxRow == " << maxRow << std::endl;
		std::cout << "maxCol == " << maxCol << std::endl;

		// 7. rowwise(), colwise()对矩阵逐行、列的操作：
		m1.resize(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		std::cout << m1 << std::endl << std::endl;
		m1.rowwise() -= RowVector3f(1, 2, 3);
		std::cout << m1 << std::endl;
	}


	// test5——矩阵的Array
	void test5()
	{
		Eigen::MatrixXf m1(5, 6);
		float* pdata = m1.data();
		for (unsigned i = 0; i < m1.size(); ++i) 
			pdata[i] = static_cast<float>(i); 
		debugDisp("m1 == \n", m1, "\n"); 

		// 1. 矩阵类的array()方法——返回矩阵对象的数组部分？？？
		Eigen::MatrixXi m111{ Eigen::MatrixXi::Ones(3, 4) };
		auto retA = m111.array();						// 返回类型为Eigen::ArrayWrapper<>
		debugDisp("typeid(aaa).name() == ", typeid(retA).name());
		debugDisp("retA == \n", retA);

		m111(0, 0) = 999; 
		debugDisp("retA == \n", retA);
		debugDisp("\n\n");

		// 按元素操作需要使用array()
		Eigen::VectorXf v1(6), v2(6);
		v1 << 1, 2, 3, 4, 5, 6;
		v2 << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
		Eigen::VectorXf result1 = v1.array() * v2.array();
		Eigen::VectorXf result2 = v1.array() / v2.array();
		debugDisp("按元素运算： \n result1 == \n", result1, "\n"); 
		debugDisp( "result2 == \n", result2, "\n");

		// array的sqrt(), pow(), isX()方法： 
		Eigen::MatrixXf result;
		result = m1.array().sqrt();
		debugDisp("按元素开方： \n", result, "\n"); 

		result = m1.array().pow(2);
		debugDisp("按元素平方： \n", result, "\n"); 

		m1.col(0).array() = -1;
		m1.col(1).array() = 0;
		m1.col(2).array() = 1;
		m1.rightCols(3).array() = 2;
		debugDisp(".array() = 常数来给矩阵元素赋值：\n", m1, "\n"); 

		Eigen::MatrixXf m2 = Eigen::MatrixXf::Ones(5, 6);
		auto result3 = (m1.array() < m2.array());
		debugDisp("类似于MATLAB中的逻辑矩阵： \n", result3, "\n"); 
		debugDisp("typeid(result3).name() == ", typeid(result3).name());

		debugDisp("finished.");
	}


	// test6—— 泛型矩阵
	void test6()
	{
		Eigen::MatrixXf m1 = Eigen::MatrixXf::Ones(5, 6);
		Matrix3i m2 = Matrix3i::Random();
		dispMat(m1);
		dispMat(m2);
	}


	// test7——flag矩阵
	void test7()
	{
		// eigen 3.3.7还没有支持flag矩阵，但是可以使用select()方法来实现类似的效果：
		Eigen::MatrixXi m1(3, 3);
		Eigen::MatrixXf m11(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m11 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

		//	1. (m1的array条件判断).select(矩阵a,矩阵b)——若m1下标为ij的元素条件判断为真，则返回矩阵a中ij位的元素，否则为矩阵b中ij位的元素。
		Eigen::MatrixXi flag1 = (m1.array() > 4).select(Matrix3i::Ones(), Matrix3i::Zero());
		debugDisp("m1 == \n", m1, "\n");
		debugDisp("flag1 == \n", flag1, "\n");

		m11 = (m11.array() < 6).select(INFINITY * Eigen::Matrix3f::Ones(), m11);
		debugDisp("m11 == \n", m11, "\n\n");

		m11(2, 2) = NAN;
		m11 = (m11.array().isNaN() || m11.array().isInf()).select(-Eigen::Matrix3f::Ones(), m11);
		debugDisp("m11 == \n", m11, "\n\n");

		// 2.
		Eigen::MatrixXf m111 = 10 * Eigen::MatrixXf::Random(4, 4);
		m11 = 10 * Eigen::MatrixXf::Random(4, 4);
		debugDisp("矩阵之间逐元素比较条件判断：");
		debugDisp("m11 == \n", m11, "\n");
		debugDisp("m111 == \n", m111, "\n");
		m11 = (m11.array() < m111.array()).select(-Eigen::Matrix4f::Ones(), m11);
		debugDisp("m11 == \n", m11, "\n\n");

		// 2. rowInMat()——检查矩阵内是否包含某一行向量——返回一个索引列向量。
		debugDisp("检查矩阵内是否包含某一行向量——返回一个索引列向量。" );
		Eigen::MatrixXi m2(5, 5);
		Eigen::Matrix4i m3;
		Eigen::MatrixXi m33;
		Eigen::RowVectorXi v(5);
		Eigen::RowVector4i vec;
		m2 << 1, 2, 3, 4, 5, 3, 1, 2, 9, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 1,1,1,1,1;
		m3 << 1, 2, 3, 4, 5, 3, 1, 2, 9, 0, 0, 0, 0, 0, 0, 4;
		m33 = m3;
		std::cout << "m2 == \n" << m2 << std::endl << std::endl;

		v << 1, 2, 3, 4, 5;
		std::cout << rowInMat(m2, v) << std::endl << std::endl;

		v << 1, 1, 1, 1, 1;
		std::cout << rowInMat(m2, v) << std::endl << std::endl;

		v << 2, 2, 8, 4, 4;
		std::cout << rowInMat(m2, v) << std::endl << std::endl;

		v.resize(4);
		v << 1, 2, 3, 4;
		vec = v;
		std::cout << rowInMat(m3, v) << std::endl << std::endl;
		std::cout << rowInMat(m33, v) << std::endl << std::endl;

		std::cout << "finished." << std::endl;
	}


	// test8——Eigen::Map类
	void test8()
	{
		Eigen::MatrixXf	m1(4, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
		m1.transposeInPlace();

		// 1. 矩阵的reshape——使用Eigen::Map类实现。
		Eigen::Map<Eigen::MatrixXf>  map1(m1.data(), 8, 2);				// reshape时元素是按存储顺序获取的，默认即按列获取。
		std::cout << "reshape之前：m1 == \n" << m1 << std::endl << std::endl;
		std::cout << "reshape之后：m11 == \n" << map1 << std::endl << std::endl;
		std::cout << "typeid(map1).name()  == " << typeid(map1).name() << std::endl;

		Eigen::MatrixXf m11{ Eigen::Map<Eigen::MatrixXf>(m1.data(), 1, 16) };
		std::cout << "m111 == \n" << m11 << std::endl << std::endl;
		std::cout << "typeid(m111).name()  == " << typeid(m11).name() << std::endl;

		// 2. map过程中没有发生数据拷贝，但是使用map对象给矩阵对象赋值，会发生数据拷贝；
		debugDisp("reinterpret_cast<unsigned>(m1.data()) == ", reinterpret_cast<unsigned>(m1.data()));
		debugDisp("reinterpret_cast<unsigned>(map1.data()) == ", reinterpret_cast<unsigned>(map1.data()));
		debugDisp("reinterpret_cast<unsigned>(m11.data()) == ", reinterpret_cast<unsigned>(m11.data()));

		// 3. map实现矩阵和向量的转换：
		Eigen::VectorXf v2 = Eigen::Map<Eigen::VectorXf>(m1.data(), 16, 1);
		debugDisp("v2 == \n", v2);
		Eigen::MatrixXf m2 = Eigen::Map<Eigen::MatrixXf>(v2.data(), 4, 4);
		debugDisp("m2 == \n", m2);

		// 4. 无法对const矩阵生成map对象
		const Eigen::MatrixXf m3 = m2;
		// Eigen::Map<Eigen::MatrixXf> map3 = Eigen::Map<Eigen::MatrixXf>(m3.data(), 16, 1);		// 无法编译；

		debugDisp("finished.");
	}


	// test9——空间变换基础（旋转、平移、缩放）
	void test9()
	{
		// 旋转(rotation)变换
		Eigen::Quaterniond q1;
		Vector3d pos(1, 2, 3);
		q1.setFromTwoVectors(Eigen::Vector3d(0, 1, 0), pos);

		Eigen::Matrix<double, 3, 3> rm1;

		// toRotationMatrix()——提取四元数中的旋转矩阵：
		rm1 = q1.toRotationMatrix();
		std::cout << "rm1 == \n" << rm1 << std::endl;

		// w(), vec()——提取四元数的分量。
		Eigen::Quaterniond q2(2, 0, 1, -3);
		std::cout << "This quaternion consists of a scalar " << q2.w()
			<< " and a vector " << std::endl << q2.vec() << std::endl;

		q2.normalize();

		std::cout << "To represent rotation, we need to normalize it such that its length is " << q2.norm() << std::endl;

		Eigen::Vector3d vec(1, 2, -1);
		Eigen::Quaterniond q3;
		q3.w() = 0;
		q3.vec() = vec;
		Eigen::Quaterniond q4 = q2 * q3 * q2.inverse();
		Eigen::Vector3d rotatedV = q4.vec();
		std::cout << "We can now use it to rotate a vector " << std::endl
			<< vec << " to " << std::endl << rotatedV << std::endl;

		// 转动到目标向量：
		Eigen::Matrix3f rotation;
		getRotationMat(rotation, RowVector3f(1, 0, 0), RowVector3f(0, 0, 1));
		RowVector3f v1 = (rotation * RowVector3f(1, 0, 0).transpose()).transpose();
		std::cout << "旋转后的v1 == " << v1 << std::endl;
		getRotationMat(rotation, RowVector3f(1, 1, 1), RowVector3f(3, 4, 5));
		v1 = (rotation * RowVector3f(1, 1, 1).transpose()).transpose();
		float scale = 5.0 / v1(2);
		v1 *= scale;
		std::cout << "旋转后的v1 == " << v1 << std::endl;

		// eigen中不支持四元数加法，只能自己逐个分量相加：
		q3.w() = q1.w() + q2.w();
		q3.vec() = q1.vec() + q2.vec();
		std::cout << std::endl << std::endl;
		dispQuat(q3);


		// Eigen::AngleAxisf类——一个类对象表示一次绕着某轴进行一定角度的旋转，可用于构造旋转矩阵。
		const float pi = 3.14159;
		Eigen::Matrix3f m;
		m = Eigen::AngleAxisf(0.25 * pi, Eigen::Vector3f::UnitX())
			* Eigen::AngleAxisf(0.5 * pi, Eigen::Vector3f::UnitY())
			* Eigen::AngleAxisf(0.33 * pi, Eigen::Vector3f::UnitZ());				// 表示先后绕xyz三个轴旋转一定角度。
		std::cout << m << std::endl << "is unitary: " << m.isUnitary() << std::endl;

		// 平移(translation)变换


		// 缩放(scaling)

	}


	// test10——空间变换应用：
	void test10()
	{
		// 已知向量在旋转变换前后的坐标，求这个旋转变换的旋转矩阵：
		Eigen::MatrixXf vers, gumline, axis;
		Eigen::MatrixXi tris;

		objReadMeshMat(vers, tris, "E:/材料/patientTooth11.obj");
		objReadVerticesMat(gumline, "E:/材料/gumline11.obj");
		objReadVerticesMat(axis, "E:/材料/axisPatient11.obj");
		
		// 以牙齿重心为原点：
		RowVector3f bary = vers.colwise().mean();
		vers.rowwise() -= bary;
		gumline.rowwise() -= bary;

		objWriteMeshMat("E:/tooth_ori.obj", vers, tris);
		objWriteVerticesMat("E:/gumline_ori.obj", gumline);

		objWriteCoorSys("E:/coordinate_ori.obj", RowVector3f::Zero(), axis.row(0), axis.row(1), axis.row(2));

		// 原世界坐标系对应着一个线性空间omega1，三轴方向向量组成的矩阵是一个单位矩阵：M1 == I;
		Eigen::Matrix3f M1 = Eigen::Matrix3f::Identity();

		// 目标坐标系对应的新线性空间omega2，三轴方向向量组成的矩阵为M2 == axis.transpose();
		Eigen::Matrix3f M2 = axis.transpose();

		// 旋转变换rotation可以将线性空间omega1变换到omega2：rotation * M1 == M2; → rotation == M2;
		Eigen::Matrix3f rotation = M2;

		// 线性空间omega2变回到omega1的变换矩阵为rotation的逆矩阵：rotationIev == rotation.inverse();
		Eigen::Matrix3f rotationInv = rotation.inverse();

		Eigen::MatrixXf versNew = rotationInv * vers.transpose();
		versNew.transposeInPlace();
		Eigen::MatrixXf gumlineNew = rotationInv * gumline.transpose();
		gumlineNew.transposeInPlace();
		Eigen::Matrix3f axisNew = rotationInv * axis.transpose();
		axisNew.transposeInPlace();
		objWriteMeshMat("E:/tooth_new.obj", versNew, tris);
		objWriteVerticesMat("E:/gumline_new.obj", gumlineNew);
		objWriteCoorSys("E:/coordinate_new.obj", RowVector3f::Zero(), axisNew.row(0), axisNew.row(1), axisNew.row(2));
	}


	// test11——使用eval()生成临时矩阵、向量
	void test11()
	{
		Eigen::MatrixXf m1(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m1 = m1.transpose();						// 实际上transpose()不产生返回值，这样做得到的结果无法预测。
		std::cout << "m1 == \n" << m1 << std::endl;

		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m1 = m1.transpose().eval();			// 使用eval()，转置的结果存储在临时矩阵中。
		std::cout << "m1 == \n" << m1 << std::endl;

		m1 = 2 * Eigen::MatrixXf::Ones(3, 3);
		Eigen::MatrixXf m2(3, 2);
		m2 << 1, 2, 3, 4, 5, 6;
		m2 = m1 * m2;							// 有时候会出现混淆的现象。
		std::cout << "m2 == \n" << m2 << std::endl;

		m1 = 2 * Eigen::MatrixXf::Ones(3, 3);
		m2 << 1, 2, 3, 4, 5, 6;
		m2 = (m1 * m2).eval();				// 这样可以确保不会出现混淆的现象。
		std::cout << "m2 == \n" << m2 << std::endl;


	}


	// test12——map类——复用数组中的数据，将其映射成向量或矩阵。
	void test12()
	{
		int intArr[] = { 1,2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
		Map<RowVectorXi> mv1(intArr, 4);
		Map<Eigen::VectorXi> mv2(intArr + 2, 5);
		std::cout << "mv1 == \n" << mv1 << std::endl;
		std::cout << "mv2 == \n" << mv2 << std::endl << std::endl;

		Map<Matrix<int, 3, 4>> mm1(intArr);
		Map<Matrix<int, 2, 6>> mm2(intArr);
		Map<Matrix<int, 3, 4, Eigen::RowMajor>> mm3(intArr);
		std::cout << "mm1 == \n" << mm1 << std::endl;
		std::cout << "mm2 == \n" << mm2 << std::endl;
		std::cout << "mm3 == \n" << mm3 << std::endl << std::endl;

		// map矩阵类对象可以赋值给矩阵类对象，发生数据的拷贝。
		Eigen::VectorXi v1 = mv1;
		Eigen::MatrixXi m1 = mm1;
		Eigen::MatrixXi m2 = Map<Matrix<int, 4, 3, Eigen::RowMajor>>(intArr);		// 这样可以快速拷贝数组中的数据，赋值到新创建的矩阵中
		dispVec<int>(v1);
		dispMat(m1);
		dispMat(m2);

		std::cout << "数组首元素地址：" << reinterpret_cast<size_t>(&intArr[0]) << std::endl;
		std::cout << "mm1首元素地址：" << reinterpret_cast<size_t>(&mm1(0, 0)) << std::endl;
		std::cout << "m1首元素地址：" << reinterpret_cast<size_t>(&m1(0, 0)) << std::endl << std::endl;

		// 映射矩阵的数据被改变，原数组的数据也改变
		mm1.setRandom();
		std::cout << "mm1 == \n" << mm1 << std::endl;
		for (const auto& num : intArr)
			std::cout << num << ", ";
		std::cout << std::endl << std::endl;

		// 不希望改变原数组中的数据的话，可以映射成常数矩阵：
		Map<const Matrix<int, 3, 4>> cmm1(intArr);				// 不是左值

	}


	// test13——类型转换：
	void test13() 
	{
		// cast()方法：
		Eigen::MatrixXf m1(5, 6);

		for (unsigned i = 0; i < m1.rows(); ++i)
			for (unsigned j = 0; j < m1.cols(); ++j)
				m1(i, j) = 10 * i + j;
		std::cout << "m1 == \n" << m1 << std::endl;

		Eigen::MatrixXi indices = (m1.array() < 33).cast<int>();
		std::cout << indices << std::endl;

		auto mu1 = m1.array().cast<unsigned>();
		std::cout << "mu1 == \n" << mu1 << std::endl;
		std::cout << typeid(mu1).name() << std::endl;

	}


	// test14——dense mat的类层次结构：
	/*
		Eigen::Matrix<>的继承链：
				Matrix<> → PlainObjectBase<> → MatrixBase<> → DenseBase<> → DenseCoeffBase<>三层 → EigenBase<>;

		Eigen::Block<>的继承链：
				Block<> → BlockImpl<> → internal::BlockImpl_dense<> → MapBase<>两层 → MatrixBase<> →...

		Eigen::CwiseNullaryOp<>的继承链——多继承：
				多继承：internal::no_assignment_operator + MatrixBase
				MatrixBase → ...

		Eigen::CwiseBinaryOp<>的继承链
				 多继承：internal::no_assignment_operator + CwiseBinaryOpImpl<>
				 CwiseBinaryOpImpl<> → MatrixBase<> → ...
	*/
	

	// 若传入参数中希望包含各种return type的对象，应使用Eigen::MatrixBase<>类型；
	template<typename Derived>										// Derived是具体的矩阵类型，如Eigen::Matrix<int,1,-1,1,1,-1>
	void foo(Eigen::MatrixBase<Derived>& base)				// Eigen::MatrixBase<>是任意稠密矩阵、向量的基类；
	{
		const Derived& originMat = base.derived();						// derived()返回原矩阵的引用，类型为原类型；
		debugDisp("typeid(base.derived()).name()  == ", typeid(originMat).name());
		debugDisp("typeid(base).name() == ", typeid(base).name());
	
		// base.resize(0, 0);

		debugDisp("foo() passed. \n\n"); 
	}


	template <typename DerivedV>
	void goo(Eigen::PlainObjectBase<DerivedV>& obj)
	{
		debugDisp("goo() passed. \n\n");
	}


	void test14() 
	{
		Eigen::RowVector3f v1;
		Eigen::VectorXi v2{ Eigen::RowVectorXi::Ones(6)};
		Eigen::Matrix3d m1(Matrix3d::Ones());
		Eigen::MatrixXi m2(5, 6);		 
		
		// 各种各样的return type:
		auto ret1 = m1.row(1);												// Eigen::Block<>
		auto ret2 = m2.col(3);													// Eigen::Block<>
		auto ret3 = m2.topRightCorner(2, 3);							// Eigen::Block<>
		auto ret4 = v1.normalized();										// Eigen::Matrix<>
		auto ret5 = Eigen::MatrixXf::Ones(2, 3);						// Eigen::CwiseNullaryOp<>
		auto ret6 = m2 * v2;													// Eigen::Product<>
		auto ret7 = m2 + Eigen::MatrixXi::Ones(5, 6);				// Eigen::CwiseBinaryOp<>
		auto ret8 = m2 - Eigen::MatrixXi::Ones(5, 6);				// Eigen::CwiseBinaryOp<>
		auto ret9 = 3 * m2;														// Eigen::CwiseBinaryOp<>

		debugDisp("typeid(ret1).name() == ", typeid(ret1).name());
		debugDisp("typeid(ret2).name() == ", typeid(ret2).name());
		debugDisp("typeid(ret3).name() == ", typeid(ret3).name());
		debugDisp("typeid(ret4).name() == ", typeid(ret4).name());
		debugDisp("typeid(ret5).name() == ", typeid(ret5).name());
		debugDisp("typeid(ret6).name() == ", typeid(ret6).name());
		debugDisp("typeid(ret7).name() == ", typeid(ret7).name());
		debugDisp("typeid(ret8).name() == ", typeid(ret8).name());
		debugDisp("typeid(ret9).name() == ", typeid(ret9).name());
		debugDisp("\n\n");

		foo(m2); 
		foo(v2); 
		foo(m1); 
		foo(v1);
		foo(ret1); foo(ret2); foo(ret3); foo(ret4); foo(ret5); foo(ret6); foo(ret7); foo(ret8); foo(ret9);  

		goo(v1); goo(v2); goo(m1); goo(m2);
		// goo(ret1); goo(ret2); goo(ret3); 
		goo(ret4); 
		// goo(ret5); goo(ret6); goo(ret7); goo(ret8); goo(ret9);


		std::cout << "finished." << std::endl;
	}

	 
}
