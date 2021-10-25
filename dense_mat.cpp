#include "dense_mat.h"

namespace DENSEMAT
{
#define MAXLEN 1024
	using namespace Eigen;
	using spMatf = SparseMatrix<float, ColMajor>;
	using TripF = Eigen::Triplet<float>;

	// дһ���ӿڣ����������ݱ��浽.dat�ļ��У�����python��ȡȻ��ͼ
	void writeData2D(const VectorXd& x, const VectorXd& y, const char* filename)
	{
		// ˳�򰤸�д��x��y�����е����ݣ���дx��дy����Ϊ���������Ƕ�Ӧ�ģ����Կ϶�ǰһ����x���꣬��һ����y���ꡣ
		double darr1[MAXLEN];
		double darr2[MAXLEN];
		unsigned int size = x.rows();
		string str1 = filename;
		string str2 = str1;

		auto iter = find(str1.begin(), str1.end(), '.');
		if (iter == str1.end())
		{
			cout << "��������Ķ������ļ������к�׺����" << endl;
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
		file.seekg(0, file.end);					// ׷�ݵ��ļ�����β��
		unsigned int size = file.tellg();			// ��ȡ�ļ����ĳ��ȡ�
		file.seekg(0, file.beg);					// �ص��ļ�����ͷ��	

													// ��һ���Ժ�����alloctor��д
		char* pc = (char*)malloc(size);
		file.read(pc, size);

		double* pd = reinterpret_cast<double*>(pc);
		for (unsigned int i = 0; i < size / sizeof(double); i++)
		{
			x[i] = *pd;
			pd++;
		}

	}


	// һЩ���͵ľ���ӿڣ�
	template <typename T, unsigned N >
	void dispVec(const Matrix<T, N, 1>& vec)		// Ԫ�����Ͳ����������������������������ݿ�����ջ��Ҳ�����ڶ��ϡ�
	{
		std::cout << vec << std::endl;
	}


	template<typename T, unsigned M, unsigned N>
	void dispMat(const Matrix<T, M, N>& mat)
	{
		cout << mat << endl;
	}


	// test0����eigen��Ļ������ݽṹ
	void test0()
	{
		// �Ѿ�����������ȷ���˳ߴ磬��δ��ʼ��,���ݴ��ڶ���
		//			�����ģ�塪��Matrix<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
		//			�Ѿ��󡪡�typedef Matrix<double, Dynamic, Dynamic> MatrixXd;
		//			����������typedef Matrix<int, Dynamic, 1> VectorXi;
		MatrixXd m1(2, 2);
		MatrixXf mf1(1, 2);
		VectorXd v1(3);			// ע����������
		cout << v1 << endl;


		// ��Array
		ArrayXXd a1(2, 2), a2(2, 2);
		a1 << 1, 2, 3, 4;
		a2 << 1, 2, 3, 4;
		cout << "a1 = \n" << a1 << endl;
		cout << "a1*a2 = \n" << a1 * a2 << endl;


		// �������������Ľӿڡ���LinSpaced()
		int start = 0;
		int end = 10;
		VectorXi vi1 = VectorXi::LinSpaced(end - start + 1, start, end);
		std::cout << "vi1 == " << vi1 << std::endl;

		VectorXf vf1 = VectorXf::LinSpaced(5, 0, 10.0);
		std::cout << "vf1 == " << vf1 << std::endl;

		// �����������Ľӿڡ���Random(), Constant(), Ones()....
		MatrixXd m2 = MatrixXd::Random(3, 3);              // ������static��������Random()���������������
		MatrixXd m3 = MatrixXd::Constant(3, 3, 1.2);		// ��������ǰ���ǳߴ磬������ֵ��
		MatrixXd m4 = MatrixXd::Ones(1, 2);					// ȫ1����


		// ���ݴ���ջ�ϵľ�������
		Matrix3d mm1 = Matrix3d::Random();
		Vector3d vv1(1, 2, 3);
		cout << "m2 = \n" << m2 << endl << endl;
		cout << "mm1 = \n" << mm1 << endl << endl;

		mm1 = m2;				//		�Ѿ����ջ��������໥��ֵ��
		cout << "mm1 = \n" << mm1 << endl << endl;


		// ��ʼ����

		//			ջ��������ʹ�ô����ų�ʼ��
		Vector3f vvv1{ 1,2,3 };
		RowVector3f vvv2{ 4,5,6 };

		std::cout << vvv1 << std::endl;
		std::cout << vvv2 << std::endl;

		//	��ֵ

		//			������������Ժ;�������໥����
		Vector3d vv2(1, 2, 3);


		//			������������Ժ;�������໥��ֵ
		MatrixXd mm(v1);
		vv2 = m2.block<3, 1>(0, 0);
		cout << "vv2 = \n" << vv2 << endl << endl;
		mm = vv2;
		cout << "mm = \n" << mm << endl << endl;
	}


	// test1�����������ʡ�Ԫ�ط��ʡ�
	void test1()
	{
		MatrixXd m1(3, 4);
		VectorXd v1(5);

		// ������������ֵ����������˳��
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		v1 << 1, 2;
		cout << "m1 = \n" << m1 << endl << endl;
		cout << "v1 = \n" << v1 << endl << endl;

		v1 << 3, 4, 5;			// �����������ĸ�ֵ�����Ǵӵ�һ��Ԫ�ؿ�ʼ��ֵ
		cout << "v1 = \n" << v1 << endl << endl;


		// �±������[]ֻ�ܻ�ȡ����Ԫ�أ���������޷�ʹ�ã���Ϊ[]ֻ֧��һ��������
		cout << "v1[0] == " << v1[0] << endl << endl;


		// �������������Ԫ�أ�ע��������0��ʼ
		cout << "m1(0, 1) ==" << m1(0, 1) << endl;
		cout << "v1(3) == " << v1(3) << endl;


		// ���������ʵ����ڽӿ�
		m1 = MatrixXd::Random(3, 4);
		cout << "m1 == \n" << m1 << endl;
		cout << "Ԫ������m1.size() == " << m1.size() << endl;
		cout << "������m1.rows() == " << m1.rows() << endl;
		cout << "������m1.cols() == " << m1.cols() << endl;
		cout << "��ͣ�sum():       " << m1.sum() << endl;
		cout << "���ˣ�prod():      " << m1.prod() << endl;
		cout << "��ֵ��mean():      " << m1.mean() << endl;
		cout << "��СԪ�أ�minCoeff():  " << m1.minCoeff() << endl;
		cout << "���Ԫ�أ�maxCoeff():  " << m1.maxCoeff() << endl;
		cout << "����ļ���trace():     " << m1.trace() << endl << endl;
		std::cout << "����ʽ��m1.determinant() == " << m1.determinant() << std::endl;

		// ����󡪡�ʹ��lu�ֽ�õ�
		cout << "�����m1.inverse() ==  \n" << m1.inverse() << endl;

		// �����ľ���任
		cout << "�����ת�ã�transpose() \n" << m1.transpose() << endl << endl;		// ����ת�õľ��󣬾�������ת��
		cout << m1 << endl;


		// bug����transpose()�ڶ�����ֵ��ʱ����ʱ����bug��Ҫ��Ҫ��������ı䣬ʹ��~InPlace()
		Matrix3f mm;
		mm << 1, 2, 3, 4, 5, 6, 0, 8, 9;
		std::cout << mm << std::endl << std::endl;
		//mm = mm.transpose();
		mm.transposeInPlace();
		std::cout << mm << std::endl << std::endl;
		Matrix3f mmCopy = mm;
		Matrix3f mm1 = mm.inverse();			// inverse()Ҳһ��������ֻ�ܴ���������������ֵ��
		std::cout << mm1 << std::endl << std::endl;
		std::cout << mm1 * mmCopy << std::endl << std::endl;

		// data()������þ��������������ָ��

		//				���ٱ��������е�Ԫ�أ�
		MatrixXf m2(5, 10);
		float* elemPtr = nullptr;
		for (unsigned i = 0; i < m2.size(); ++i)
		{
			elemPtr = m2.data() + i;
			*elemPtr = i;
		}
		std::cout << m2 << std::endl << std::endl;;


		//				�����Ѿ����е����ݣ�
		MatrixXf m22(5, 10);
		std::memcpy(m22.data(), m2.data(), m2.size() * sizeof(float));
		std::cout << m22 << std::endl;


		//	 minCoeff() ���� ��������Ԫ���е���ֵ
		cout << "m1 == \n" << m1 << endl;
		MatrixXd::Index maxRow, maxCol;
		MatrixXd::Index minRow, minCol;
		double min = m1.minCoeff(&minRow, &minCol);
		double max = m1.maxCoeff(&maxRow, &maxCol);
		std::cout << "��СԪ���к�" << minRow << "���к�" << minCol << std::endl;
		std::cout << "���Ԫ���к�" << maxRow << "���к�" << maxCol << std::endl;

		v1.resize(10);
		v1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
		VectorXd::Index minIdx;
		min = v1.minCoeff(&minIdx);
		std::cout << "��������СԪ�أ�" << min << ", ������" << minIdx << std::endl;

		// all(), any(), count()��������Ԫ�ص�����
		ArrayXXf a1(3, 3), a2(2, 2);
		a1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		a2 << 1, 2, 3, 4;
		bool flag1 = (a1 > 0).all();						// all()�����������Ԫ�ض�Ϊtrue���򷵻�true,
		bool flag2 = (a2 == 4).any();					// any()�������������һ��Ԫ��Ϊtrue���򷵻�true;
		Eigen::Index count = (a1 > 5).count();		// count()������������������Ԫ����

		VectorXi vec(9);
		vec << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		ArrayXi aa = vec.array();
		flag1 = (aa == 2).any();
		std::cout << "flag1" << flag1 << std::endl;

		// ò��ֻ�ܶԵ���Ԫ��ʹ�ã����ܶ�һ��Ԫ��ʹ�á�
		m1.resize(3, 3);
		m1 << 1, 2, 3, 1, 2, 3, 4, 5, 6;
		RowVector3d vd1(1, 2, 3);
		std::cout << "m1���Ƿ���������(1,2,3)��" << (m1.array() == vd1.array()).any() << std::endl;
		m1 = MatrixXd::Ones(3, 3);
		std::cout << "m1���Ƿ���������(1,2,3)��" << (m1.array() == vd1.array()).any() << std::endl;
	}


	// test2������������任������
	void test2()
	{
		MatrixXf m1(3, 4), m3(4, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

		cout << "m1 = \n" << endl;
		cout << m1 << endl << endl;

		// resize()�������ľ���ߴ硪��ֻ�����ڶѾ��󡢶�����  ���������ܾ���ߴ���ĺ���Σ�����Ԫ��ȫ��������ֵ��
		m3.resize(3, 3);
		m3.row(2) = VectorXf::Ones(3);
		cout << "m3.size() == " << m3.size() << endl;
		cout << "m3.rows() == " << m3.rows() << endl;
		cout << "m3.cols() == " << m3.cols() << endl;
		cout << "m3 = \n" << m3 << endl << endl;

		// conservativeResize() ���ľ���ߴ磬����ԭ������
		m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m3.conservativeResize(4, 3);
		m3.row(3) = 5 * VectorXf::Ones(3);
		cout << "m3 = \n" << m3 << endl << endl;

		// colwise(), rowwise()���������С����в���
		cout << m1.colwise().sum() << endl << endl;				// ������ͣ�ѹ��һ��������
		cout << m1.rowwise().sum() << endl << endl;				// ������ͣ�ѹ��һ����������

		// �������š���ͨ������˷�ʵ��
		MatrixXf  a(1, 3), b(3, 1);
		a << 1, 2, 3;
		b << 4, 5, 6;
		//		����������N�� == ���ones(N,1)
		auto aa = MatrixXf::Ones(3, 1) * a;
		//		����������ΪN�� == �ҳ�ones(1,N)
		auto bb = b * MatrixXf::Ones(1, 4);

		cout << aa << endl << endl;;
		cout << bb << endl << endl;


		// reverse();
		m1.resize(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		cout << "m1 = \n" << m1 << endl << endl;
		MatrixXf result1 = m1.reverse();								// ���ж���װ
		MatrixXf result2 = m1.colwise().reverse();				// ������Ϊ��������Ȼ��ÿһ���ڲ����ֲ��䵹װ��
		MatrixXf result3 = m1.rowwise().reverse();				// ������Ϊ������....
		cout << "result1 = \n" << result1 << endl << endl;
		cout << "result2 = \n" << result2 << endl << endl;
		cout << "result3 = \n" << result3 << endl << endl;
		m1.colwise().reverseInPlace();
		cout << "m1 = \n" << m1 << endl << endl;




		// ������ˡ����
		Vector3f v1(1, 2, 3);
		Vector3f v2(3, 2, 1);
		cout << "v1.dot(v2) == \n" << v1.dot(v2) << endl;
		cout << "v1.cross(v2) == " << v1.cross(v2) << endl << endl;;

		// ������һ��
		v1.normalize();					// normalize()�Ὣ���������һ����normalized()ֻ�Ƿ��ع�һ�����������ı�����
		cout << "v1 = \n" << v1 << endl << endl;
		cout << "v2.normalized() = \n" << v2.normalized() << endl << endl;
		cout << "v2 = \n" << v2 << endl << endl;



	}


	// test3����������ؿ�ѧ����Ľӿ�
	void test3()
	{
		MatrixXd A(3, 3);
		A << 1, 2, 3, 4, 5, 6, 7, 8, 9;

		// ����������ֵ������������
		EigenSolver<Matrix3d> es(A);
		Matrix3d D = es.pseudoEigenvalueMatrix();			// �Խ���Ԫ��������ֵ
		Matrix3d V = es.pseudoEigenvectors();				// ÿһ����������������������
		cout << "����ֵ����D��" << endl << D << endl;
		cout << "������������V: " << endl << V << endl;
		cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;

		cout << A * V.block<3, 1>(0, 0) << endl << endl;
		cout << D(0, 0) * V.block<3, 1>(0, 0) << endl << endl;

		// �����LU�ֽ⡪��Eigen::FullPivLU
		MatrixXf m = MatrixXf::Random(5, 5);
		cout << "Here is the matrix m:" << endl << m << endl;
		Eigen::FullPivLU<MatrixXf> lu(m);
		cout << "ִ��LU�ֽ⣺lu.matrixLU() == :"
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


		// ���������ֵ�ֽ⡪��
		m = MatrixXf::Random(3, 3);
		cout << "m == \n" << m << endl;
		JacobiSVD<MatrixXf> svd(m, ComputeThinU | ComputeThinV);
		cout << "����ֵ == " << endl << svd.singularValues() << endl;
		cout << "U == " << endl << svd.matrixU() << endl;
		cout << "V == " << endl << svd.matrixV() << endl;
		Vector3f rhs(1, 0, 0);
		cout << "ʹ������ֵ�ֽ�����Է����飺" << endl << svd.solve(rhs) << endl << endl;


		// �����Լ���װ��ʹ��QR�ֽ�����ǡ���������Է���������ӣ�
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
		cout << "���Է�����Ax == b�Ľ⣺" << endl << x << endl << endl;
		cout << "���Է�����AX == B�Ľ⣺" << endl << X << endl << endl;

		// norm()����ŷ����÷�����Ҳ��p==2ʱ��lp������������Ԫ��ƽ���͵Ŀ�����
		m.resize(2, 2);
		m << 1, 2, 3, 4;
		std::cout << "m.norm() == " << m.norm() << std::endl;
		Vector3f v1(1, 2, 3);
		std::cout << "v1.norm() == " << v1.norm() << std::endl;	// �����ķ�������������ģ����

	}


	// test4��������������
	void test4()
	{
		MatrixXf m1(5, 6);
		float* pdata = m1.data();
		for (unsigned i = 0; i < m1.size(); ++i)
		{
			pdata[i] = static_cast<float>(i);
		}
		std::cout << "m1 == \n" << m1 << std::endl << std::endl;

		// VectorXT::segment<>()����������ȡ������������Ƭ�Σ���������ֵ���á���������
		VectorXi v1(9);
		v1 << 1, 2, 3, 4, 5, 6, 7, 8, 8;

		//			segment()����1����segment<elemCount>( startIdx)����̬segment��Ƭ�γ���elemCount�����ڱ�������֪����constExpr
		std::cout << "v1.segment<3>(1)  == \n" << v1.segment<3>(1) << std::endl;

		//			segment()����2����segment(startIdx, elemCount)����̬segment��Ƭ�γ���elemCount�����ڱ�������δ֪��
		std::cout << "v1.segment(1, v1(2)) == \n" << v1.segment(1, v1(2)) << std::endl;

		//			segment()���ص�����Ƭ������ֵ���ã�
		v1.segment<5>(2) = 999 * VectorXi::Ones(5);
		std::cout << "v1 == \n" << v1 << std::endl << std::endl;


		// block()����������ȡ�Ӿ���飬��������ֵ���á�����������
		std::cout << "block()����" << std::endl;

		//			����1��block<rows, cols>(startRow, startCol)������̬block�������ĳߴ�����Ǳ�������֪�ģ�������constExpr
		cout << "m1.block<2, 2>(1, 1)  == \n" << m1.block<2, 2>(1, 1) << endl << endl;

		//			����2: block(startRow, startCol, rows, cols)������̬block�������ĳߴ�����Ǳ�����δ֪��
		cout << "m1.block(1, 0, m1(2, 0), m1(3, 0)) == \n" << m1.block(1, 0, m1(2, 0), m1(3, 0)) << endl << endl;


		//			 ���صľ��������ֵ���á����Ա�<<��ֵ��
		m1.block<1, 4>(1, 0) << 88, 99, 111, 222;
		cout << "m1 = \n" << m1 << endl << endl;
		m1.block<1, 4>(1, 0) << 5, 6, 7, 8;


		//			���صľ�������ÿ���ֱ����' = ' ����һ������ֵ��
		MatrixXf m2 = MatrixXf::Ones(1, 6);
		m1.block<1, 6>(0, 0) = m2;
		std::cout << "m1 == \n" << m1 << std::endl;


		// row(), col()������ȡĳһ�л���ĳһ��
		std::cout << "m1.col(3) == \n" << m1.col(3) << std::endl;

		//			��ȡ�����С�����ԭ���ݵ����ã������µĿ�������������ֵ��
		m1.col(3) << 1, 2, 3, 4, 5;
		std::cout << "m1.col(3) == \n" << m1.col(3) << std::endl;

		//	minCoeff()����ֵ�Ծ���ֿ�Ҳͬ������
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


		// rowwise(), colwise()�Ծ������С��еĲ�����
		m1.resize(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		std::cout << m1 << std::endl << std::endl;
		m1.rowwise() -= RowVector3f(1, 2, 3);
		std::cout << m1 << std::endl;
	}


	// test5���������Array
	void test5()
	{
		MatrixXf m1(5, 6);
		float* pdata = m1.data();
		for (unsigned i = 0; i < m1.size(); ++i)
		{
			pdata[i] = static_cast<float>(i);
		}
		std::cout << "m1 == \n" << m1 << std::endl;


		// ��Ԫ�ز�����Ҫʹ��array()
		MatrixXf result = m1.array() * MatrixXf::Ones(5, 6).array();
		std::cout << "��Ԫ����ˣ� \n" << result << std::endl << std::endl;

		result = m1.array().sqrt();
		std::cout << "��Ԫ�ؿ����� \n" << result << std::endl << std::endl;

		result = m1.array().pow(2);
		std::cout << "��Ԫ��ƽ���� \n" << result << std::endl << std::endl;

		m1.col(0).array() = -1;
		m1.col(1).array() = 0;
		m1.col(2).array() = 1;
		m1.rightCols(3).array() = 2;
		std::cout << ".array() = ������������Ԫ�ظ�ֵ��\n" << m1 << std::endl << std::endl;


		MatrixXf m2 = MatrixXf::Ones(5, 6);
		auto result2 = (m1.array() < m2.array());
		std::cout << "������MATLAB�е��߼����� \n" << result2 << std::endl << std::endl;
		std::cout << typeid(result2).name() << std::endl;
	}


	// test6���� ���;���
	void test6()
	{
		MatrixXf m1 = MatrixXf::Ones(5, 6);
		Matrix3i m2 = Matrix3i::Random();
		dispMat(m1);
		dispMat(m2);
	}


	// test7������������
	void test7()
	{
		// eigen 3.3.7��û��֧���������󣬵��ǿ���ʹ��select()������ʵ�����Ƶ�Ч����
		Matrix3i m1;
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

		//		(array�����ж�).select(����a,����b)��array���±�Ϊij��Ԫ�������ж�Ϊ�棬�򷵻ؾ���ijλ��Ԫ��Ϊ����a��ijλ��Ԫ�أ�����Ϊ����b��ijλ��Ԫ�ء�
		Matrix3i flag1 = (m1.array() > 4).select(Matrix3i::Ones(), Matrix3i::Zero());
		std::cout << "m1 == \n" << m1 << std::endl;
		std::cout << "flag1 == \n" << flag1 << std::endl;
		std::cout << std::endl;

		// vecInMat()�������������Ƿ����ĳһ��������������һ��������������
		std::cout << "���������Ƿ����ĳһ��������������һ��������������" << std::endl;
		MatrixXi m2(3, 5);
		m2 << 1, 2, 3, 4, 5, 3, 1, 2, 9, 0, 0, 0, 0, 0, 0;
		RowVectorXi v(5);
		v << 1, 2, 3, 4, 5;

		std::cout << vecInMat<int>(m2, v) << std::endl << std::endl;;

		MatrixXi m3(1, 5);
		m3 << 1, 2, 3, 4, 5;
		std::cout << vecInMat<int>(m3, v) << std::endl << std::endl;;

		m3 << 1, 2, 3, 4, 4;
		std::cout << vecInMat<int>(m3, v) << std::endl << std::endl;;


	}


	// test8����Eigen::Map��
	void test8()
	{
		// �����reshape����ʹ��Eigen::Map��ʵ�֡�
		MatrixXf	m1(4, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
		m1.transposeInPlace();
		Eigen::Map<MatrixXf>  m11(m1.data(), 8, 2);			// reshapeʱԪ���ǰ��洢˳���ȡ�ģ�Ĭ�ϼ����л�ȡ��
		std::cout << "reshape֮ǰ��m1 == \n" << m1 << std::endl << std::endl;
		std::cout << "reshape֮��m11 == \n" << m11 << std::endl << std::endl;


	}


	// test9�����ռ�任��������ת��ƽ�ơ����ţ�
	void test9()
	{
		// ��ת(rotation)�任
		Eigen::Quaterniond q1;
		Vector3d pos(1, 2, 3);
		q1.setFromTwoVectors(Eigen::Vector3d(0, 1, 0), pos);

		Eigen::Matrix<double, 3, 3> rm1;

		// toRotationMatrix()������ȡ��Ԫ���е���ת����
		rm1 = q1.toRotationMatrix();
		std::cout << "rm1 == \n" << rm1 << std::endl;

		// w(), vec()������ȡ��Ԫ���ķ�����
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


		// eigen�в�֧����Ԫ���ӷ���ֻ���Լ����������ӣ�
		q3.w() = q1.w() + q2.w();
		q3.vec() = q1.vec() + q2.vec();
		std::cout << std::endl << std::endl;
		dispQuat(q3);


		// AngleAxisf�ࡪ��һ��������ʾһ������ĳ�����һ���Ƕȵ���ת�������ڹ�����ת����
		const float M_PI = 3.14159;
		Matrix3f m;
		m = AngleAxisf(0.25 * M_PI, Vector3f::UnitX())
			* AngleAxisf(0.5 * M_PI, Vector3f::UnitY())
			* AngleAxisf(0.33 * M_PI, Vector3f::UnitZ());				// ��ʾ�Ⱥ���xyz��������תһ���Ƕȡ�
		cout << m << endl << "is unitary: " << m.isUnitary() << endl;

		// ƽ��(translation)�任


		// ����(scaling)

	}


	// test10�����ռ�任Ӧ�ã�
	void test10()
	{
		// ��֪��������ת�任ǰ������꣬�������ת�任����ת����

		MatrixXf vers, gumline, axis;
		MatrixXi tris;

		objReadMeshMat(vers, tris, "E:/����/patientTooth11.obj");
		objReadVerticesMat(gumline, "E:/����/gumline11.obj");
		objReadVerticesMat(axis, "E:/����/axisPatient11.obj");
		
		// ����������Ϊԭ�㣺
		RowVector3f bary = vers.colwise().mean();
		vers.rowwise() -= bary;
		gumline.rowwise() -= bary;

		objWriteMeshMat("E:/tooth_ori.obj", vers, tris);
		objWriteVerticesMat("E:/gumline_ori.obj", gumline);

		printCoordinateEigen("E:/coordinate_ori.obj", RowVector3f::Zero(), axis.row(0), axis.row(1), axis.row(2));

		// ԭ��������ϵ��Ӧ��һ�����Կռ�omega1�����᷽��������ɵľ�����һ����λ����M1 == I;
		Matrix3f M1 = Matrix3f::Identity();

		// Ŀ������ϵ��Ӧ�������Կռ�omega2�����᷽��������ɵľ���ΪM2 == axis.transpose();
		Matrix3f M2 = axis.transpose();

		// ��ת�任rotation���Խ����Կռ�omega1�任��omega2��rotation * M1 == M2; �� rotation == M2;
		Matrix3f rotation = M2;

		// ���Կռ�omega2��ص�omega1�ı任����Ϊrotation�������rotationIev == rotation.inverse();
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


	// test11����ʹ��eval()������ʱ��������
	void test11()
	{
		MatrixXf m1(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m1 = m1.transpose();						// ʵ����transpose()����������ֵ���������õ��Ľ���޷�Ԥ�⡣
		std::cout << "m1 == \n" << m1 << std::endl;

		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m1 = m1.transpose().eval();			// ʹ��eval()��ת�õĽ���洢����ʱ�����С�
		std::cout << "m1 == \n" << m1 << std::endl;

		m1 = 2 * MatrixXf::Ones(3, 3);
		MatrixXf m2(3, 2);
		m2 << 1, 2, 3, 4, 5, 6;
		m2 = m1 * m2;							// ��ʱ�����ֻ���������
		std::cout << "m2 == \n" << m2 << std::endl;

		m1 = 2 * MatrixXf::Ones(3, 3);
		m2 << 1, 2, 3, 4, 5, 6;
		m2 = (m1 * m2).eval();				// ��������ȷ��������ֻ���������
		std::cout << "m2 == \n" << m2 << std::endl;


	}


	// test12����map�ࡪ�����������е����ݣ�����ӳ������������
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

		// map�����������Ը�ֵ����������󣬷������ݵĿ�����
		VectorXi v1 = mv1;
		MatrixXi m1 = mm1;
		MatrixXi m2 = Map<Matrix<int, 4, 3, RowMajor>>(intArr);		// �������Կ��ٿ��������е����ݣ���ֵ���´����ľ�����
		dispVec<int>(v1);
		dispMat<int>(m1);
		dispMat<int>(m2);

		std::cout << "������Ԫ�ص�ַ��" << reinterpret_cast<unsigned>(&intArr[0]) << std::endl;
		std::cout << "mm1��Ԫ�ص�ַ��" << reinterpret_cast<unsigned>(&mm1(0, 0)) << std::endl;
		std::cout << "m1��Ԫ�ص�ַ��" << reinterpret_cast<unsigned>(&m1(0, 0)) << std::endl << std::endl;

		// ӳ���������ݱ��ı䣬ԭ���������Ҳ�ı�
		mm1.setRandom();
		std::cout << "mm1 == \n" << mm1 << std::endl;
		for (const auto& num : intArr)
		{
			std::cout << num << ", ";
		}
		std::cout << std::endl << std::endl;

		// ��ϣ���ı�ԭ�����е����ݵĻ�������ӳ��ɳ�������
		Map<const Matrix<int, 3, 4>> cmm1(intArr);				// ������ֵ

	}


	// ��ʱ�޷����ࣺ
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

		// ����������
		std::vector<int> rowVec{ 1, 3, 4 };
		std::vector<int> colVec{ 5, 2, 0 };

		std::vector<int> ind{ 4,2,5,5,3 };
		MatrixXi A = MatrixXi::Random(4, 6);
		cout << "Initial matrix A:\n" << A << "\n\n";


		// �鿴eigen��İ汾��
		std::cout << EIGEN_WORLD_VERSION << std::endl;
		std::cout << EIGEN_MAJOR_VERSION << std::endl;
		std::cout << EIGEN_MINOR_VERSION << std::endl;		// �汾Ϊ3.3.7;
	}

}
