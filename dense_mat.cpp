#include "dense_mat.h"

namespace DENSEMAT
{
#define MAXLEN 1024
	using namespace Eigen;
	using spMatf = SparseMatrix<float, ColMajor>;
	using TripF = Eigen::Triplet<float>;


	// 写一个接口，将矩阵数据保存到.dat文件中，方便python读取然后画图
	void writeData2D(const VectorXd& x, const VectorXd& y, const char* filename)
	{
		// 顺序挨个写入x和y向量中的数据，先写x再写y，因为两条向量是对应的，所以肯定前一半是x坐标，后一半是y坐标。
		double darr1[MAXLEN];
		double darr2[MAXLEN];
		unsigned int size = x.rows();
		string str1 = filename;
		string str2 = str1;

		auto iter = find(str1.begin(), str1.end(), '.');
		if (iter == str1.end())
		{
			cout << "错误，输出的二进制文件必须有后缀名。" << endl;
			return;
		}


		auto dis = distance(str1.begin(), iter);
		str1.insert(dis, "_x");
		str2.insert(dis, "_y");


		for (unsigned int i = 0; i < size; i++)
		{
			darr1[i] = x(i);
		}
		for (unsigned int i = 0; i < size; i++)
		{
			darr2[i] = y(i);
		}

		ofstream file1(str1, ios::out | ios::binary);
		ofstream file2(str2, ios::out | ios::binary);

		file1.write(reinterpret_cast<char*>(&darr1[0]), size * sizeof(double));
		file2.write(reinterpret_cast<char*>(&darr2[0]), size * sizeof(double));
		file1.close();
		file2.close();
	}


	void readData(VectorXd& x, const char* filename)
	{
		ifstream file(filename, ios::in | ios::binary);
		file.seekg(0, file.end);					// 追溯到文件流的尾部
		unsigned int size = file.tellg();			// 获取文件流的长度。
		file.seekg(0, file.beg);					// 回到文件流的头部	

													// 这一块以后考虑用alloctor改写
		char* pc = (char*)malloc(size);
		file.read(pc, size);

		double* pd = reinterpret_cast<double*>(pc);
		for (unsigned int i = 0; i < size / sizeof(double); i++)
		{
			x[i] = *pd;
			pd++;
		}

	}


	// 一些泛型的矩阵接口：
	template <typename T, unsigned N >
	void dispVec(const Matrix<T, N, 1>& vec)		// 元素类型不定，行数不定的列向量，且数据可以在栈上也可以在堆上。
	{
		std::cout << vec << std::endl;
	}


	template<typename T, unsigned M, unsigned N>
	void dispMat(const Matrix<T, M, N>& mat)
	{
		cout << mat << endl;
	}


	// test0——eigen库的基本数据结构
	void test0()
	{
		// 堆矩阵、向量——确定了尺寸，但未初始化,数据存在堆上
		/*
			最基本模板——Matrix<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
			堆矩阵——typedef Matrix<double, Dynamic, Dynamic> MatrixXd;
			堆向量——typedef Matrix<int, Dynamic, 1> VectorXi;
		*/

		MatrixXd m1(2, 2);
		MatrixXf mf1(1, 2);
		VectorXd v1(3);			// 注意是列向量
		cout << v1 << endl;

		// 堆Array
		ArrayXXd a1(2, 2), a2(2, 2);
		a1 << 1, 2, 3, 4;
		a2 << 1, 2, 3, 4;
		cout << "a1 = \n" << a1 << endl;
		cout << "a1*a2 = \n" << a1 * a2 << endl;

		// 生成特殊向量的接口——LinSpaced()
		int start = 0;
		int end = 10;
		VectorXi vi1 = VectorXi::LinSpaced(end - start + 1, start, end);
		std::cout << "vi1 == " << vi1 << std::endl;

		VectorXf vf1 = VectorXf::LinSpaced(5, 0, 10.0);
		std::cout << "vf1 == " << vf1 << std::endl;

		// 生成特殊矩阵的接口——Random(), Constant(), Ones()....
		MatrixXd m2 = MatrixXd::Random(3, 3);              // 矩阵类static方法——Random()——返回随机矩阵
		MatrixXd m3 = MatrixXd::Constant(3, 3, 1.2);		// 常数矩阵，前面是尺寸，后面是值；
		MatrixXd m4 = MatrixXd::Ones(1, 2);					// 全1矩阵


		// 数据存在栈上的矩阵类型
		Matrix3d mm1 = Matrix3d::Random();
		Vector3d vv1(1, 2, 3);
		cout << "m2 = \n" << m2 << endl << endl;
		cout << "mm1 = \n" << mm1 << endl << endl;

		mm1 = m2;				//		堆矩阵和栈矩阵可以相互赋值。
		cout << "mm1 = \n" << mm1 << endl << endl;


		// 初始化：

		//			栈向量可以使用大括号初始化
		Vector3f vvv1{ 1,2,3 };
		RowVector3f vvv2{ 4,5,6 };

		std::cout << vvv1 << std::endl;
		std::cout << vvv2 << std::endl;

		//	赋值

		//			列向量对象可以和矩阵对象相互构造
		Vector3d vv2(1, 2, 3);


		//			列向量对象可以和矩阵对象相互赋值
		MatrixXd mm(v1);
		vv2 = m2.block<3, 1>(0, 0);
		cout << "vv2 = \n" << vv2 << endl << endl;
		mm = vv2;
		cout << "mm = \n" << mm << endl << endl;
	}


	// test1——矩阵性质、元素访问。
	void test1()
	{
		MatrixXd m1(3, 4);
		VectorXd v1(5);

		// 输出流运算符赋值——行优先顺序填充
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		v1 << 1, 2;
		cout << "m1 = \n" << m1 << endl << endl;
		cout << "v1 = \n" << v1 << endl << endl;

		v1 << 3, 4, 5;			// 不会接着上面的赋值，而是从第一个元素开始赋值
		cout << "v1 = \n" << v1 << endl << endl;


		// 下标运算符[]只能获取向量元素，矩阵对象无法使用，因为[]只支持一个参数。
		cout << "v1[0] == " << v1[0] << endl << endl;

		// 括号运算符访问元素，注意索引从0开始
		cout << "m1(0, 1) ==" << m1(0, 1) << endl;
		cout << "v1(3) == " << v1(3) << endl;

		// 求矩阵的性质的类内接口
		m1 = MatrixXd::Random(3, 4);
		cout << "m1 == \n" << m1 << endl;
		cout << "元素数：m1.size() == " << m1.size() << endl;
		cout << "行数：m1.rows() == " << m1.rows() << endl;
		cout << "列数：m1.cols() == " << m1.cols() << endl;
		cout << "求和：sum():       " << m1.sum() << endl;
		cout << "连乘：prod():      " << m1.prod() << endl;
		cout << "均值：mean():      " << m1.mean() << endl;
		cout << "最小元素：minCoeff():  " << m1.minCoeff() << endl;
		cout << "最大元素：maxCoeff():  " << m1.maxCoeff() << endl;
		cout << "矩阵的迹：trace():     " << m1.trace() << endl << endl;
		std::cout << "行列式：m1.determinant() == " << m1.determinant() << std::endl;

		// 逆矩阵——使用lu分解得到
		cout << "逆矩阵：m1.inverse() ==  \n" << m1.inverse() << endl;

		// 基本的矩阵变换
		cout << "矩阵的转置：transpose() \n" << m1.transpose() << endl << endl;		// 返回转置的矩阵，矩阵自身不转置
		cout << m1 << endl;


		// bug——transpose()在对自身赋值的时候有时会有bug，要想要矩阵自身改变，使用~InPlace()
		Matrix3f mm;
		mm << 1, 2, 3, 4, 5, 6, 0, 8, 9;
		std::cout << mm << std::endl << std::endl;
		//mm = mm.transpose();
		mm.transposeInPlace();
		std::cout << mm << std::endl << std::endl;
		Matrix3f mmCopy = mm;
		Matrix3f mm1 = mm.inverse();			// inverse()也一样，但是只能创建其他变量来赋值。
		std::cout << mm1 << std::endl << std::endl;
		std::cout << mm1 * mmCopy << std::endl << std::endl;

		// data()——获得矩阵数据数组的首指针

		//				快速遍历矩阵中的元素：
		MatrixXf m2(5, 10);
		float* elemPtr = nullptr;
		for (unsigned i = 0; i < m2.size(); ++i)
		{
			elemPtr = m2.data() + i;
			*elemPtr = i;
		}
		std::cout << m2 << std::endl << std::endl;;


		//				拷贝堆矩阵中的数据：
		MatrixXf m22(5, 10);
		std::memcpy(m22.data(), m2.data(), m2.size() * sizeof(float));
		std::cout << m22 << std::endl;


		//	 minCoeff() —— 搜索矩阵元素中的最值
		cout << "m1 == \n" << m1 << endl;
		MatrixXd::Index maxRow, maxCol;
		MatrixXd::Index minRow, minCol;
		double min = m1.minCoeff(&minRow, &minCol);
		double max = m1.maxCoeff(&maxRow, &maxCol);
		std::cout << "最小元素行号" << minRow << "，列号" << minCol << std::endl;
		std::cout << "最大元素行号" << maxRow << "，列号" << maxCol << std::endl;

		v1.resize(10);
		v1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
		VectorXd::Index minIdx;
		min = v1.minCoeff(&minIdx);
		std::cout << "向量中最小元素：" << min << ", 索引：" << minIdx << std::endl;

		// all(), any(), count()——矩阵元素的搜索
		ArrayXXf a1(3, 3), a2(2, 2);
		a1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		a2 << 1, 2, 3, 4;
		bool flag1 = (a1 > 0).all();						// all()——如果所有元素都为true，则返回true,
		bool flag2 = (a2 == 4).any();					// any()——如果至少有一个元素为true，则返回true;
		Eigen::Index count = (a1 > 5).count();		// count()——返回满足条件的元素数

		VectorXi vec(9);
		vec << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		ArrayXi aa = vec.array();
		flag1 = (aa == 2).any();
		std::cout << "flag1" << flag1 << std::endl;

		// 貌似只能对单个元素使用，不能对一组元素使用。
		m1.resize(3, 3);
		m1 << 1, 2, 3, 1, 2, 3, 4, 5, 6;
		RowVector3d vd1(1, 2, 3);
		std::cout << "m1中是否有行向量(1,2,3)：" << (m1.array() == vd1.array()).any() << std::endl;
		m1 = MatrixXd::Ones(3, 3);
		std::cout << "m1中是否有行向量(1,2,3)：" << (m1.array() == vd1.array()).any() << std::endl;
	}


	// test2——矩阵基本变换、运算
	void test2()
	{
		MatrixXf m1(3, 4), m3(4, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

		cout << "m1 = \n" << endl;
		cout << m1 << endl << endl;

		// resize()——更改矩阵尺寸——只适用于堆矩阵、堆向量  ！！！不管矩阵尺寸更改后如何，所有元素全部变成随机值。
		m3.resize(3, 3);
		m3.row(2) = VectorXf::Ones(3);
		cout << "m3.size() == " << m3.size() << endl;
		cout << "m3.rows() == " << m3.rows() << endl;
		cout << "m3.cols() == " << m3.cols() << endl;
		cout << "m3 = \n" << m3 << endl << endl;

		// conservativeResize() 更改矩阵尺寸，保留原有数据
		m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m3.conservativeResize(4, 3);
		m3.row(3) = 5 * VectorXf::Ones(3);
		cout << "m3 = \n" << m3 << endl << endl;

		// colwise(), rowwise()——矩阵按行、按列操作
		cout << m1.colwise().sum() << endl << endl;				// 按列求和，压成一个行向量
		cout << m1.rowwise().sum() << endl << endl;				// 按行求和，压成一个列向量。

		// 矩阵扩张——通过矩阵乘法实现
		MatrixXf  a(1, 3), b(3, 1);
		a << 1, 2, 3;
		b << 4, 5, 6;
		//		行向量扩张N行 == 左乘ones(N,1)
		auto aa = MatrixXf::Ones(3, 1) * a;
		//		列向量扩张为N列 == 右乘ones(1,N)
		auto bb = b * MatrixXf::Ones(1, 4);

		cout << aa << endl << endl;;
		cout << bb << endl << endl;

		// reverse();
		m1.resize(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		cout << "m1 = \n" << m1 << endl << endl;
		MatrixXf result1 = m1.reverse();								// 行列都倒装
		MatrixXf result2 = m1.colwise().reverse();				// 矩阵视为列向量，然后每一行内部保持不变倒装；
		MatrixXf result3 = m1.rowwise().reverse();				// 矩阵视为行向量....
		cout << "result1 = \n" << result1 << endl << endl;
		cout << "result2 = \n" << result2 << endl << endl;
		cout << "result3 = \n" << result3 << endl << endl;
		m1.colwise().reverseInPlace();
		cout << "m1 = \n" << m1 << endl << endl;

		// 向量点乘、叉乘
		Vector3f v1(1, 2, 3);
		Vector3f v2(3, 2, 1);
		cout << "v1.dot(v2) == \n" << v1.dot(v2) << endl;
		cout << "v1.cross(v2) == " << v1.cross(v2) << endl << endl;;

		// 向量归一化
		v1.normalize();					// normalize()会将向量本身归一化，normalized()只是返回归一化向量，不改变自身
		cout << "v1 = \n" << v1 << endl << endl;
		cout << "v2.normalized() = \n" << v2.normalized() << endl << endl;
		cout << "v2 = \n" << v2 << endl << endl;
	}


	// test3——矩阵相关科学计算的接口
	void test3()
	{
		MatrixXd A(3, 3);
		A << 1, 2, 3, 4, 5, 6, 7, 8, 9;

		// 1. 求矩阵的特征值、特征向量。
		EigenSolver<Matrix3d> es(A);
		Matrix3d D = es.pseudoEigenvalueMatrix();			// 对角线元素是特征值
		Matrix3d V = es.pseudoEigenvectors();				// 每一个列向量都是特征向量。
		cout << "特征值矩阵D：" << endl << D << endl;
		cout << "特征向量矩阵V: " << endl << V << endl;
		cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;

		cout << A * V.block<3, 1>(0, 0) << endl << endl;
		cout << D(0, 0) * V.block<3, 1>(0, 0) << endl << endl;

		// 2. 矩阵的LU分解——Eigen::FullPivLU
		MatrixXf m = MatrixXf::Random(5, 5);
		cout << "Here is the matrix m:" << endl << m << endl;
		Eigen::FullPivLU<MatrixXf> lu(m);
		cout << "执行LU分解：lu.matrixLU() == :"
			<< endl << lu.matrixLU() << endl;

		MatrixXf L, U, P_inv, Q_inv;
		L = MatrixXf::Identity(5, 5);
		L.triangularView<StrictlyLower>() = lu.matrixLU();
		U = lu.matrixLU().triangularView<Upper>();
		P_inv = lu.permutationP().inverse();
		Q_inv = lu.permutationQ().inverse();

		std::cout << "P_inv  == \n" << P_inv << std::endl << std::endl;
		std::cout << "L == \n" << L << std::endl << std::endl;
		std::cout << "U == \n" << U << std::endl << std::endl;
		std::cout << "Q_inv  == \n" << Q_inv << std::endl << std::endl;

		cout << "Let us now reconstruct the original matrix m:" << endl;
		cout << P_inv * L * U * Q_inv << endl;


		// 3. 矩阵的奇异值分解——
		m = MatrixXf::Random(3, 3);
		cout << "m == \n" << m << endl;
		JacobiSVD<MatrixXf> svd(m, ComputeThinU | ComputeThinV);
		cout << "奇异值 == " << endl << svd.singularValues() << endl;
		cout << "U == " << endl << svd.matrixU() << endl;
		cout << "V == " << endl << svd.matrixV() << endl;
		Vector3f rhs(1, 0, 0);
		cout << "使用奇异值分解解线性方程组：" << endl << svd.solve(rhs) << endl << endl;

		// 3+. 求秩：
		std::cout << "svd.rank == " << svd.rank() << std::endl << std::endl;

		// 4. 测试自己封装的使用QR分解的求解恰定稠密线性方程组的轮子：
		A = MatrixXd::Random(3, 3);
		MatrixXd B = MatrixXd::Random(3, 3);
		Vector3d b = Vector3d::Random();
		Vector3d x;
		MatrixXd X;
		solveLinearEquation<double>(x, A, b);
		solveLinearEquations<double>(X, A, B);
		std::cout << "A == \n" << A << std::endl;
		std::cout << "b == \n" << b << std::endl;
		std::cout << "B == \n" << B << std::endl;
		cout << "线性方程组Ax == b的解：" << endl << x << endl << endl;
		cout << "线性方程组AX == B的解：" << endl << X << endl << endl;

		// norm()——欧几里得范数，也是p==2时的lp范数。既所有元素平方和的开方。
		m.resize(2, 2);
		m << 1, 2, 3, 4;
		std::cout << "m.norm() == " << m.norm() << std::endl;
		Vector3f v1(1, 2, 3);
		std::cout << "v1.norm() == " << v1.norm() << std::endl;	// 向量的范数等于向量的模长。

		//			求行向量的模长：
		Matrix3f arrows;
		arrows << 0, 1, 0, 0, 3, 4, 0, 5, 12;
		Vector3f arrowsLen = arrows.rowwise().norm();
		dispVec<float>(arrowsLen);
	}


	// test4——矩阵块操作：
	void test4()
	{
		MatrixXf m1(5, 6);
		float* pdata = m1.data();
		for (unsigned i = 0; i < m1.size(); ++i)
		{
			pdata[i] = static_cast<float>(i);
		}
		std::cout << "m1 == \n" << m1 << std::endl << std::endl;

		// VectorXT::segment<>()方法——提取向量的子向量片段，返回其左值引用。两个重载
		VectorXi v1(9);
		v1 << 1, 2, 3, 4, 5, 6, 7, 8, 8;

		//			segment()重载1——segment<elemCount>( startIdx)，静态segment，片段长度elemCount必须在编译期已知，是constExpr
		std::cout << "v1.segment<3>(1)  == \n" << v1.segment<3>(1) << std::endl;

		//			segment()重载2——segment(startIdx, elemCount)，动态segment，片段长度elemCount可以在编译器内未知。
		std::cout << "v1.segment(1, v1(2)) == \n" << v1.segment(1, v1(2)) << std::endl;

		//			segment()返回的向量片段是左值引用：
		v1.segment<5>(2) = 999 * VectorXi::Ones(5);
		std::cout << "v1 == \n" << v1 << std::endl << std::endl;


		// block()方法——提取子矩阵块，返回其左值引用。有两个重载
		std::cout << "block()方法" << std::endl;

		//			重载1：block<rows, cols>(startRow, startCol)——静态block，矩阵块的尺寸必须是编译期已知的，必须是constExpr
		cout << "m1.block<2, 2>(1, 1)  == \n" << m1.block<2, 2>(1, 1) << endl << endl;

		//			重载2: block(startRow, startCol, rows, cols)——动态block，矩阵块的尺寸可以是编译期未知的
		cout << "m1.block(1, 0, m1(2, 0), m1(3, 0)) == \n" << m1.block(1, 0, m1(2, 0), m1(3, 0)) << endl << endl;


		//			 返回的矩阵块是左值引用。可以被<<赋值；
		m1.block<1, 4>(1, 0) << 88, 99, 111, 222;
		cout << "m1 = \n" << m1 << endl << endl;
		m1.block<1, 4>(1, 0) << 5, 6, 7, 8;


		//			返回的矩阵块引用可以直接用' = ' 被另一个矩阵赋值。
		MatrixXf m2 = MatrixXf::Ones(1, 6);
		m1.block<1, 6>(0, 0) = m2;
		std::cout << "m1 == \n" << m1 << std::endl;


		// row(), col()——提取某一行或者某一列
		std::cout << "m1.col(3) == \n" << m1.col(3) << std::endl;

		//			提取出的行、列是原数据的引用，不是新的拷贝，并且是左值。
		m1.col(3) << 1, 2, 3, 4, 5;
		std::cout << "m1.col(3) == \n" << m1.col(3) << std::endl;

		//	minCoeff()求最值对矩阵分块也同样适用
		MatrixXf::Index maxRow, maxCol;
		MatrixXf::Index minRow, minCol;
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


		// rowwise(), colwise()对矩阵逐行、列的操作：
		m1.resize(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		std::cout << m1 << std::endl << std::endl;
		m1.rowwise() -= RowVector3f(1, 2, 3);
		std::cout << m1 << std::endl;
	}


	// test5——矩阵的Array
	void test5()
	{
		MatrixXf m1(5, 6);
		float* pdata = m1.data();
		for (unsigned i = 0; i < m1.size(); ++i)
		{
			pdata[i] = static_cast<float>(i);
		}
		std::cout << "m1 == \n" << m1 << std::endl;

		// 按元素操作需要使用array()
		MatrixXf result = m1.array() * MatrixXf::Ones(5, 6).array();
		std::cout << "按元素相乘： \n" << result << std::endl << std::endl;

		result = m1.array().sqrt();
		std::cout << "按元素开方： \n" << result << std::endl << std::endl;

		result = m1.array().pow(2);
		std::cout << "按元素平方： \n" << result << std::endl << std::endl;

		m1.col(0).array() = -1;
		m1.col(1).array() = 0;
		m1.col(2).array() = 1;
		m1.rightCols(3).array() = 2;
		std::cout << ".array() = 常数来给矩阵元素赋值：\n" << m1 << std::endl << std::endl;


		MatrixXf m2 = MatrixXf::Ones(5, 6);
		auto result2 = (m1.array() < m2.array());
		std::cout << "类似于MATLAB中的逻辑矩阵： \n" << result2 << std::endl << std::endl;
		std::cout << typeid(result2).name() << std::endl;
	}


	// test6—— 泛型矩阵
	void test6()
	{
		MatrixXf m1 = MatrixXf::Ones(5, 6);
		Matrix3i m2 = Matrix3i::Random();
		dispMat(m1);
		dispMat(m2);
	}


	// test7——索引矩阵
	void test7()
	{
		// eigen 3.3.7还没有支持索引矩阵，但是可以使用select()方法来实现类似的效果：
		Matrix3i m1;
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

		//		(array条件判断).select(矩阵a,矩阵b)，array中下标为ij的元素条件判断为真，则返回矩阵ij位的元素为矩阵a中ij位的元素，否则为矩阵b中ij位的元素。
		Matrix3i flag1 = (m1.array() > 4).select(Matrix3i::Ones(), Matrix3i::Zero());
		std::cout << "m1 == \n" << m1 << std::endl;
		std::cout << "flag1 == \n" << flag1 << std::endl;
		std::cout << std::endl;

		// vecInMat()——检查矩阵内是否包含某一行向量——返回一个索引列向量。
		std::cout << "检查矩阵内是否包含某一行向量——返回一个索引列向量。" << std::endl;
		MatrixXi m2(3, 5);
		m2 << 1, 2, 3, 4, 5, 3, 1, 2, 9, 0, 0, 0, 0, 0, 0;
		RowVectorXi v(5);
		v << 1, 2, 3, 4, 5;

		std::cout << vecInMat<int>(m2, v) << std::endl << std::endl;;

		MatrixXi m3(1, 5);
		m3 << 1, 2, 3, 4, 5;
		std::cout << vecInMat<int>(m3, v) << std::endl << std::endl;;

		m3 << 1, 2, 3, 4, 4;
		std::cout << vecInMat<int>(m3, v) << std::endl << std::endl;
	}


	// test8——Eigen::Map类
	void test8()
	{
		// 矩阵的reshape——使用Eigen::Map类实现。
		MatrixXf	m1(4, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
		m1.transposeInPlace();
		Eigen::Map<MatrixXf>  m11(m1.data(), 8, 2);			// reshape时元素是按存储顺序获取的，默认即按列获取。
		std::cout << "reshape之前：m1 == \n" << m1 << std::endl << std::endl;
		std::cout << "reshape之后：m11 == \n" << m11 << std::endl << std::endl;


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
		Matrix3f rotation = getRotationMat(RowVector3f(1, 0, 0), RowVector3f(0, 0, 1));
		RowVector3f v1 = (rotation * RowVector3f(1, 0, 0).transpose()).transpose();
		std::cout << "旋转后的v1 == " << v1 << std::endl;
		rotation = getRotationMat(RowVector3f(1, 1, 1), RowVector3f(3, 4, 5));
		v1 = (rotation * RowVector3f(1, 1, 1).transpose()).transpose();
		float scale = 5.0 / v1(2);
		v1 *= scale;
		std::cout << "旋转后的v1 == " << v1 << std::endl;

		// eigen中不支持四元数加法，只能自己逐个分量相加：
		q3.w() = q1.w() + q2.w();
		q3.vec() = q1.vec() + q2.vec();
		std::cout << std::endl << std::endl;
		dispQuat(q3);


		// AngleAxisf类——一个类对象表示一次绕着某轴进行一定角度的旋转，可用于构造旋转矩阵。
		const float M_PI = 3.14159;
		Matrix3f m;
		m = AngleAxisf(0.25 * M_PI, Vector3f::UnitX())
			* AngleAxisf(0.5 * M_PI, Vector3f::UnitY())
			* AngleAxisf(0.33 * M_PI, Vector3f::UnitZ());				// 表示先后绕xyz三个轴旋转一定角度。
		cout << m << endl << "is unitary: " << m.isUnitary() << endl;

		// 平移(translation)变换


		// 缩放(scaling)

	}


	// test10——空间变换应用：
	void test10()
	{
		// 已知向量在旋转变换前后的坐标，求这个旋转变换的旋转矩阵：

		MatrixXf vers, gumline, axis;
		MatrixXi tris;

		objReadMeshMat(vers, tris, "E:/材料/patientTooth11.obj");
		objReadVerticesMat(gumline, "E:/材料/gumline11.obj");
		objReadVerticesMat(axis, "E:/材料/axisPatient11.obj");
		
		// 以牙齿重心为原点：
		RowVector3f bary = vers.colwise().mean();
		vers.rowwise() -= bary;
		gumline.rowwise() -= bary;

		objWriteMeshMat("E:/tooth_ori.obj", vers, tris);
		objWriteVerticesMat("E:/gumline_ori.obj", gumline);

		printCoordinateEigen("E:/coordinate_ori.obj", RowVector3f::Zero(), axis.row(0), axis.row(1), axis.row(2));

		// 原世界坐标系对应着一个线性空间omega1，三轴方向向量组成的矩阵是一个单位矩阵：M1 == I;
		Matrix3f M1 = Matrix3f::Identity();

		// 目标坐标系对应的新线性空间omega2，三轴方向向量组成的矩阵为M2 == axis.transpose();
		Matrix3f M2 = axis.transpose();

		// 旋转变换rotation可以将线性空间omega1变换到omega2：rotation * M1 == M2; → rotation == M2;
		Matrix3f rotation = M2;

		// 线性空间omega2变回到omega1的变换矩阵为rotation的逆矩阵：rotationIev == rotation.inverse();
		Matrix3f rotationInv = rotation.inverse();

		MatrixXf versNew = rotationInv * vers.transpose();
		versNew.transposeInPlace();
		MatrixXf gumlineNew = rotationInv * gumline.transpose();
		gumlineNew.transposeInPlace();
		Matrix3f axisNew = rotationInv * axis.transpose();
		axisNew.transposeInPlace();
		objWriteMeshMat("E:/tooth_new.obj", versNew, tris);
		objWriteVerticesMat("E:/gumline_new.obj", gumlineNew);
		printCoordinateEigen("E:/coordinate_new.obj", RowVector3f::Zero(), axisNew.row(0), axisNew.row(1), axisNew.row(2));
	}


	// test11——使用eval()生成临时矩阵、向量
	void test11()
	{
		MatrixXf m1(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m1 = m1.transpose();						// 实际上transpose()不产生返回值，这样做得到的结果无法预测。
		std::cout << "m1 == \n" << m1 << std::endl;

		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m1 = m1.transpose().eval();			// 使用eval()，转置的结果存储在临时矩阵中。
		std::cout << "m1 == \n" << m1 << std::endl;

		m1 = 2 * MatrixXf::Ones(3, 3);
		MatrixXf m2(3, 2);
		m2 << 1, 2, 3, 4, 5, 6;
		m2 = m1 * m2;							// 有时候会出现混淆的现象。
		std::cout << "m2 == \n" << m2 << std::endl;

		m1 = 2 * MatrixXf::Ones(3, 3);
		m2 << 1, 2, 3, 4, 5, 6;
		m2 = (m1 * m2).eval();				// 这样可以确保不会出现混淆的现象。
		std::cout << "m2 == \n" << m2 << std::endl;


	}


	// test12——map类——复用数组中的数据，将其映射成向量或矩阵。
	void test12()
	{
		int intArr[] = { 1,2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
		Map<RowVectorXi> mv1(intArr, 4);
		Map<VectorXi> mv2(intArr + 2, 5);
		std::cout << "mv1 == \n" << mv1 << std::endl;
		std::cout << "mv2 == \n" << mv2 << std::endl << std::endl;

		Map<Matrix<int, 3, 4>> mm1(intArr);
		Map<Matrix<int, 2, 6>> mm2(intArr);
		Map<Matrix<int, 3, 4, RowMajor>> mm3(intArr);
		std::cout << "mm1 == \n" << mm1 << std::endl;
		std::cout << "mm2 == \n" << mm2 << std::endl;
		std::cout << "mm3 == \n" << mm3 << std::endl << std::endl;

		// map矩阵类对象可以赋值给矩阵类对象，发生数据的拷贝。
		VectorXi v1 = mv1;
		MatrixXi m1 = mm1;
		MatrixXi m2 = Map<Matrix<int, 4, 3, RowMajor>>(intArr);		// 这样可以快速拷贝数组中的数据，赋值到新创建的矩阵中
		dispVec<int>(v1);
		dispMat<int>(m1);
		dispMat<int>(m2);

		std::cout << "数组首元素地址：" << reinterpret_cast<unsigned>(&intArr[0]) << std::endl;
		std::cout << "mm1首元素地址：" << reinterpret_cast<unsigned>(&mm1(0, 0)) << std::endl;
		std::cout << "m1首元素地址：" << reinterpret_cast<unsigned>(&m1(0, 0)) << std::endl << std::endl;

		// 映射矩阵的数据被改变，原数组的数据也改变
		mm1.setRandom();
		std::cout << "mm1 == \n" << mm1 << std::endl;
		for (const auto& num : intArr)
		{
			std::cout << num << ", ";
		}
		std::cout << std::endl << std::endl;

		// 不希望改变原数组中的数据的话，可以映射成常数矩阵：
		Map<const Matrix<int, 3, 4>> cmm1(intArr);				// 不是左值

	}


	// 暂时无法分类：
	void test000()
	{
		MatrixXf m1(5, 6);

		for (unsigned i = 0; i < m1.rows(); ++i)
		{
			for (unsigned j = 0; j < m1.cols(); ++j)
			{
				m1(i, j) = 10 * i + j;
			}
		}
		std::cout << "m1 == \n" << m1 << std::endl;


		MatrixXi indices = (m1.array() < 33).cast<int>();
		std::cout << indices << std::endl;

		// cast()
		auto mu1 = m1.array().cast<unsigned>();
		std::cout << "mu1 == \n" << mu1 << std::endl;
		std::cout << typeid(mu1).name() << std::endl;

		// 矩阵索引：
		std::vector<int> rowVec{ 1, 3, 4 };
		std::vector<int> colVec{ 5, 2, 0 };

		std::vector<int> ind{ 4,2,5,5,3 };
		MatrixXi A = MatrixXi::Random(4, 6);
		cout << "Initial matrix A:\n" << A << "\n\n";


		// 查看eigen库的版本：
		std::cout << EIGEN_WORLD_VERSION << std::endl;
		std::cout << EIGEN_MAJOR_VERSION << std::endl;
		std::cout << EIGEN_MINOR_VERSION << std::endl;		// 版本为3.3.7;
	}

}
