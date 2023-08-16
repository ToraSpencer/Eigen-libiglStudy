#include "dense_mat.h"

namespace DENSEMAT
{
#define MAXLEN 1024
	using namespace Eigen;
	using spMatf = Eigen::SparseMatrix<float, ColMajor>;
	using TripF = Eigen::Triplet<float>;

	////////////////////////////////////////////////////////////////////////////////////////////// DEBUG 接口

	namespace MY_DEBUG
	{
		static std::string g_debugPath = "E:/";

		// 写一个接口，将矩阵数据保存到.dat文件中，方便python读取然后画图
		void writeData2D(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const char* filename)
		{
			// 顺序挨个写入x和y向量中的数据，先写x再写y，因为两条向量是对应的，所以肯定前一半是x坐标，后一半是y坐标。
			double darr1[MAXLEN];
			double darr2[MAXLEN];
			unsigned int size = x.rows();
			std::string str1 = filename;
			std::string str2 = str1;

			auto iter = find(str1.begin(), str1.end(), '.');
			if (iter == str1.end())
			{
				std::cout << "错误，输出的二进制文件必须有后缀名。" << std::endl;
				return;
			}


			auto dis = distance(str1.begin(), iter);
			str1.insert(dis, "_x");
			str2.insert(dis, "_y");

			for (unsigned int i = 0; i < size; i++)
				darr1[i] = x(i);
			for (unsigned int i = 0; i < size; i++)
				darr2[i] = y(i);

			std::ofstream file1(str1, std::ios::out | std::ios::binary);
			std::ofstream file2(str2, std::ios::out | std::ios::binary);

			file1.write(reinterpret_cast<char*>(&darr1[0]), size * sizeof(double));
			file2.write(reinterpret_cast<char*>(&darr2[0]), size * sizeof(double));
			file1.close();
			file2.close();
		}


		void readData(Eigen::VectorXd& x, const char* filename)
		{
			std::ifstream file(filename, std::ios::in | std::ios::binary);
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


		template <typename T, int M, int N>
		static void dispData(const Eigen::Matrix<T, M, N>& m)
		{
			auto dataPtr = m.data();
			unsigned elemsCount = m.size();

			for (unsigned i = 0; i < elemsCount; ++i)
				std::cout << dataPtr[i] << ", ";

			std::cout << std::endl;
		}


		template <typename Derived>
		static void dispData(const Eigen::PlainObjectBase<Derived>& m)
		{
			int m0 = m.RowsAtCompileTime;
			int n0 = m.ColsAtCompileTime;

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
		static void debugWriteVers(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
		{
			char path[512] = { 0 };
			sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
			objWriteVerticesMat(path, vers);
		}


		template<typename T>
		static void debugWriteMesh(const char* name, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
		{
			char path[512] = { 0 };
			sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
			objWriteMeshMat(path, vers, tris);
		}



		template<typename DerivedV>
		static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers)
		{
			char path[512] = { 0 };
			sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
			objWriteEdgesMat(path, edges, vers);
		}
	}
	using namespace MY_DEBUG;
 

	// test0——eigen库的基本数据结构
	void test0()
	{
		// 堆矩阵、向量——确定了尺寸，但未初始化,数据存在堆上
		/*
			基本模板——Matrix<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
			堆矩阵——typedef Matrix<double, Dynamic, Dynamic> Eigen::MatrixXd;
			堆向量——typedef Matrix<int, Dynamic, 1> Eigen::VectorXi;

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

		// 堆Array
		ArrayXXd a1(2, 2), a2(2, 2);
		a1 << 1, 2, 3, 4;											// operator << 填充元素是行优先的顺序；
		a2 << 1, 2, 3, 4;
		std::cout << "a1 = \n" << a1 << std::endl;
		std::cout << "a1*a2 = \n" << a1 * a2 << std::endl;

		// 生成特殊向量的接口——LinSpaced(元素数，起点，终点)
		int start = 0;
		int end = 10;
		Eigen::VectorXi vi1 = Eigen::VectorXi::LinSpaced(end - start + 1, start, end);
		std::cout << "vi1 == " << vi1 << std::endl;

		Eigen::VectorXf vf1 = Eigen::VectorXf::LinSpaced(5, 0, 10.0);
		std::cout << "vf1 == " << vf1 << std::endl;

		// 生成特殊矩阵的接口——Random(), Constant(), Ones()....
		Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(3, 3);              // 矩阵类static方法——Random()——返回随机矩阵
		Eigen::MatrixXd m3 = Eigen::MatrixXd::Constant(3, 3, 1.2);		// 常数矩阵，前面是尺寸，后面是值；
		Eigen::MatrixXd m4 = Eigen::MatrixXd::Ones(1, 2);					// 全1矩阵
 
		// 数据存在栈上的矩阵类型
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

		// 输出流运算符赋值——行优先顺序填充
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		v1 << 1, 2;
		std::cout << "m1 = \n" << m1 << std::endl << std::endl;
		std::cout << "v1 = \n" << v1 << std::endl << std::endl;

		v1 << 3, 4, 5;			// 不会接着上面的赋值，而是从第一个元素开始赋值
		std::cout << "v1 = \n" << v1 << std::endl << std::endl;


		// 下标运算符[]只能获取向量元素，矩阵对象无法使用，因为[]只支持一个参数。
		std::cout << "v1[0] == " << v1[0] << std::endl << std::endl;

		// 括号运算符访问元素，注意索引从0开始
		std::cout << "m1(0, 1) ==" << m1(0, 1) << std::endl;
		std::cout << "v1(3) == " << v1(3) << std::endl;

		// 求矩阵的性质的类内接口
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

		// outerSize(), innerSize()——默认的列优先存储的矩阵，outerSize为列数，行优先存储的outerSize为行数；
		std::cout << "m1.outerSize() == " << m1.outerSize() << std::endl;			
		std::cout << "m1.innerSize() == " << m1.innerSize() << std::endl;
		std::cout << "m1_rm.outerSize() == " << m1_rm.outerSize() << std::endl;
		std::cout << "m1_rm.innerSize() == " << m1_rm.innerSize() << std::endl;
 
		// 逆矩阵——使用lu分解得到
		std::cout << "逆矩阵：m1.inverse() ==  \n" << m1.inverse() << std::endl;

		// 基本的矩阵变换
		std::cout << "矩阵的转置：transpose() \n" << m1.transpose() << std::endl << std::endl;		// 返回转置的矩阵，矩阵自身不转置
		std::cout << m1 << std::endl;


		// bug——transpose()在对自身赋值的时候有时会有bug，要想要矩阵自身改变，使用~InPlace()
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
		Eigen::VectorXd::Index minIdx;
		min = v1.minCoeff(&minIdx);
		std::cout << "向量中最小元素：" << min << ", 索引：" << minIdx << std::endl;

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

		Eigen::MatrixXf aaa, bbbb;
		kron(aaa, RowVector3f::Ones(), a);
		kron(bbbb, Matrix2i::Ones(), b);
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


	// test3——矩阵相关科学计算的接口
	void test3()
	{
		Eigen::MatrixXd A(3, 3);
		A << 1, 2, 3, 4, 5, 6, 7, 8, 9;

		// 1. 求矩阵的特征值、特征向量。
		EigenSolver<Matrix3d> es(A);
		Matrix3d D = es.pseudoEigenvalueMatrix();			// 对角线元素是特征值
		Matrix3d V = es.pseudoEigenvectors();					// 每一个列向量都是特征向量。
		std::cout << "特征值矩阵D：" << std::endl << D << std::endl;
		std::cout << "特征向量矩阵V: " << std::endl << V << std::endl;
		std::cout << "Finally, V * D * V^(-1) = " << std::endl << V * D * V.inverse() << std::endl;

		std::cout << A * V.block<3, 1>(0, 0) << std::endl << std::endl;
		std::cout << D(0, 0) * V.block<3, 1>(0, 0) << std::endl << std::endl;
		debugDisp("\n\n");


		// 2. 矩阵的LU分解——Eigen::FullPivLU
		Eigen::MatrixXf m = Eigen::MatrixXf::Random(5, 5);
		std::cout << "Here is the matrix m:" << std::endl << m << std::endl;
		Eigen::FullPivLU<Eigen::MatrixXf> lu(m);
		std::cout << "执行LU分解：lu.matrixLU() == :"
			<< std::endl << lu.matrixLU() << std::endl;

		Eigen::MatrixXf L, U, P_inv, Q_inv;
		L = Eigen::MatrixXf::Identity(5, 5);
		L.triangularView<StrictlyLower>() = lu.matrixLU();
		U = lu.matrixLU().triangularView<Upper>();
		P_inv = lu.permutationP().inverse();
		Q_inv = lu.permutationQ().inverse();

		std::cout << "P_inv  == \n" << P_inv << std::endl << std::endl;
		std::cout << "L == \n" << L << std::endl << std::endl;
		std::cout << "U == \n" << U << std::endl << std::endl;
		std::cout << "Q_inv  == \n" << Q_inv << std::endl << std::endl;

		std::cout << "Let us now reconstruct the original matrix m:" << std::endl;
		std::cout << P_inv * L * U * Q_inv << std::endl;
		debugDisp("\n\n");


		// 3. 矩阵的奇异值分解（SVD）——
		m = Eigen::MatrixXf::Random(3, 3);
		std::cout << "m == \n" << m << std::endl;
		JacobiSVD<Eigen::MatrixXf> svdSolver(m, ComputeThinU | ComputeThinV);
		std::cout << "奇异值 == " << std::endl << svdSolver.singularValues() << std::endl;
		std::cout << "U == " << std::endl << svdSolver.matrixU() << std::endl;
		std::cout << "V == " << std::endl << svdSolver.matrixV() << std::endl;
		debugDisp("\n");

		//		3.1 通过SVD求解线性方程组：
		Eigen::Vector3f rhs(1, 0, 0);
		std::cout << "使用奇异值分解解线性方程组：" << std::endl << svdSolver.solve(rhs) << std::endl << std::endl;
		debugDisp("\n");

		//		3.2 通过SVD求矩阵的条件数：
		double condNum = calcCondNum(m);
		debugDisp("矩阵条件数：condNum == ", condNum);
		debugDisp("\n\n");

		// 4. 求秩：
		std::cout << "svdSolver.rank == " << svdSolver.rank() << std::endl << std::endl;

		// 5. 测试自己封装的使用QR分解的求解恰定稠密线性方程组的轮子：
		A = Eigen::MatrixXd::Random(3, 3);
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(3, 3);
		Vector3d b = Vector3d::Random();
		Vector3d x;
		Eigen::MatrixXd X;
		solveLinearEquation<double>(x, A, b);
		solveLinearEquations<double>(X, A, B);
		std::cout << "A == \n" << A << std::endl;
		std::cout << "b == \n" << b << std::endl;
		std::cout << "B == \n" << B << std::endl;
		std::cout << "线性方程组Ax == b的解：" << std::endl << x << std::endl << std::endl;
		std::cout << "线性方程组AX == B的解：" << std::endl << X << std::endl << std::endl;
		debugDisp("\n\n");

		// 6. norm()——欧几里得范数，也是p==2时的lp范数。既所有元素平方和的开方。
		Eigen::Vector3f v1(1, 2, 3);
		std::cout << "v1.norm() == " << v1.norm() << std::endl;	// 向量的范数等于向量的模长。
		m.resize(3, 4);
		m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 , 11, 12;
		std::cout << "m == \n" << m << std::endl;
		std::cout << "m.norm() == " << m.norm() << std::endl;
		std::cout << "m.rowwise().norm() == \n" << m.rowwise().norm() << std::endl;		// 所有行向量求norm，返回一个列向量；
		std::cout << "m.colwise().norm() == \n" << m.colwise().norm() << std::endl;

		//	6.1. 通过范数求行向量的模长：
		Eigen::Matrix3f arrows;
		arrows << 0, 1, 0, 0, 3, 4, 0, 5, 12;
		Eigen::Vector3f arrowsLen = arrows.rowwise().norm();
		dispVec<float>(arrowsLen);
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


		// setConstant()方法——块赋值：
		m1.setConstant(1.0);
		m1.topRightCorner(2, 3).setConstant(2.0);
		dispMat(m1);


		// VectorXT::segment<>()方法——提取向量的子向量片段，返回其左值引用。两个重载
		Eigen::VectorXi v1(9);
		v1 << 1, 2, 3, 4, 5, 6, 7, 8, 8;

		//			segment()重载1——segment<elemCount>( startIdx)，静态segment，片段长度elemCount必须在编译期已知，是constExpr
		std::cout << "v1.segment<3>(1)  == \n" << v1.segment<3>(1) << std::endl;

		//			segment()重载2——segment(startIdx, elemCount)，动态segment，片段长度elemCount可以在编译器内未知。
		std::cout << "v1.segment(1, v1(2)) == \n" << v1.segment(1, v1(2)) << std::endl;

		//			segment()返回的向量片段是左值引用：
		v1.segment<5>(2) = 999 * Eigen::VectorXi::Ones(5);
		std::cout << "v1 == \n" << v1 << std::endl << std::endl;


		// block()方法——提取子矩阵块，返回其左值引用。有两个重载
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


		// row(), col()——提取某一行或者某一列
		std::cout << "row(), col()——提取某一行或者某一列\n m1.col(3) == \n" << m1.col(3) << std::endl;

		//			提取出的行、列是原数据的引用，不是新的拷贝，并且是左值。
		m1.col(3) << 1, 2, 3, 4, 5;
		std::cout << "m1.col(3) == \n" << m1.col(3) << std::endl;

		// topRows(), bottomRows(), leftCols(), rightCols();
		std::cout << "m1.leftCols(2) ==   \n" << m1.leftCols(2) << std::endl;
		std::cout << "m1.topRows(2) == \n" << m1.topRows(2) << std::endl;
 

		//	minCoeff()求最值对矩阵分块也同样适用
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
		Eigen::MatrixXf m1(5, 6);
		float* pdata = m1.data();
		for (unsigned i = 0; i < m1.size(); ++i)
		{
			pdata[i] = static_cast<float>(i);
		}
		std::cout << "m1 == \n" << m1 << std::endl;

		// 按元素操作需要使用array()
		Eigen::MatrixXf result = m1.array() * Eigen::MatrixXf::Ones(5, 6).array();
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


		Eigen::MatrixXf m2 = Eigen::MatrixXf::Ones(5, 6);
		auto result2 = (m1.array() < m2.array());
		std::cout << "类似于MATLAB中的逻辑矩阵： \n" << result2 << std::endl << std::endl;
		std::cout << typeid(result2).name() << std::endl;
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
		Matrix3i m1;
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

		//		(array条件判断).select(矩阵a,矩阵b)，array中下标为ij的元素条件判断为真，则返回矩阵ij位的元素为矩阵a中ij位的元素，否则为矩阵b中ij位的元素。
		Matrix3i flag1 = (m1.array() > 4).select(Matrix3i::Ones(), Matrix3i::Zero());
		std::cout << "m1 == \n" << m1 << std::endl << std::endl;
		std::cout << "flag1 == \n" << flag1 << std::endl << std::endl;
		std::cout << std::endl;

		// rowInMat()——检查矩阵内是否包含某一行向量——返回一个索引列向量。
		std::cout << "检查矩阵内是否包含某一行向量——返回一个索引列向量。" << std::endl;
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
		// 矩阵的reshape——使用Eigen::Map类实现。
		Eigen::MatrixXf	m1(4, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
		m1.transposeInPlace();
		Eigen::Map<Eigen::MatrixXf>  m11(m1.data(), 8, 2);			// reshape时元素是按存储顺序获取的，默认即按列获取。
		std::cout << "reshape之前：m1 == \n" << m1 << std::endl << std::endl;
		std::cout << "reshape之后：m11 == \n" << m11 << std::endl << std::endl;
		std::cout << "typeid(m11).name()  == " << typeid(m11).name() << std::endl;

		Eigen::MatrixXf m111{ Eigen::Map<Eigen::MatrixXf>(m1.data(), 1, 16) };
		std::cout << "m111 == \n" << m111 << std::endl << std::endl;
		std::cout << "typeid(m111).name()  == " << typeid(m111).name() << std::endl;

		Eigen::VectorXf v1 = Eigen::Map<Eigen::VectorXf>(m1.data(), 16, 1);
		dispVec(v1);

		std::cout << "finished." << std::endl;
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
		Eigen::Matrix3f rotation = getRotationMat(RowVector3f(1, 0, 0), RowVector3f(0, 0, 1));
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


	template<typename Derived>											// Derived是具体的矩阵类型，如Eigen::Matrix<int,1,-1,1,1,-1>
	void foo(const Eigen::MatrixBase<Derived>& base)			// Eigen::MatrixBase<>是任意稠密矩阵、向量的基类；
	{
		const Derived& originMat = base.derived();			// derived()返回原矩阵的引用，类型为原类型；
		std::cout << "typeid(base).name() == " << typeid(base).name() << std::endl;
		std::cout << "typeid(base.derived()).name()  == " << typeid(originMat).name() << std::endl;
		std::cout << "typeid(Derived).name() == " << typeid(Derived).name() << std::endl;

		std::cout << "address of base: " << reinterpret_cast<size_t>(&base) << std::endl;
		std::cout << "address of base.derived(): " << reinterpret_cast<size_t>(&originMat) << std::endl;
		std::cout << std::endl;
	}

	// test14——dense mat的类层次结构：
	void test14() 
	{
		Eigen::Vector3f v1;
		Eigen::RowVectorXi v2;
		Eigen::Matrix3d m1(Matrix3d::Ones());
		Eigen::MatrixXi m2(5, 6);

		foo(v1);
		foo(v2);
		foo(m1);
		foo(m2);
		std::cout << "finished." << std::endl;
	}


	// 暂时无法分类：
	void test000()
	{
		// 矩阵索引：
		std::vector<int> rowVec{ 1, 3, 4 };
		std::vector<int> colVec{ 5, 2, 0 };

		std::vector<int> ind{ 4,2,5,5,3 };
		Eigen::MatrixXi A = Eigen::MatrixXi::Random(4, 6);
		std::cout << "Initial matrix A:\n" << A << "\n\n";

		// 查看eigen库的版本：
		std::cout << EIGEN_WORLD_VERSION << std::endl;
		std::cout << EIGEN_MAJOR_VERSION << std::endl;
		std::cout << EIGEN_MINOR_VERSION << std::endl;		// 版本为3.3.7;
	}

}
