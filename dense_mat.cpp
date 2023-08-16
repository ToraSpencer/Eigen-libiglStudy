#include "dense_mat.h"

namespace DENSEMAT
{
#define MAXLEN 1024
	using namespace Eigen;
	using spMatf = Eigen::SparseMatrix<float, ColMajor>;
	using TripF = Eigen::Triplet<float>;

	////////////////////////////////////////////////////////////////////////////////////////////// DEBUG �ӿ�

	namespace MY_DEBUG
	{
		static std::string g_debugPath = "E:/";

		// дһ���ӿڣ����������ݱ��浽.dat�ļ��У�����python��ȡȻ��ͼ
		void writeData2D(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const char* filename)
		{
			// ˳�򰤸�д��x��y�����е����ݣ���дx��дy����Ϊ���������Ƕ�Ӧ�ģ����Կ϶�ǰһ����x���꣬��һ����y���ꡣ
			double darr1[MAXLEN];
			double darr2[MAXLEN];
			unsigned int size = x.rows();
			std::string str1 = filename;
			std::string str2 = str1;

			auto iter = find(str1.begin(), str1.end(), '.');
			if (iter == str1.end())
			{
				std::cout << "��������Ķ������ļ������к�׺����" << std::endl;
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
 

	// test0����eigen��Ļ������ݽṹ
	void test0()
	{
		// �Ѿ�����������ȷ���˳ߴ磬��δ��ʼ��,���ݴ��ڶ���
		/*
			����ģ�塪��Matrix<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
			�Ѿ��󡪡�typedef Matrix<double, Dynamic, Dynamic> Eigen::MatrixXd;
			����������typedef Matrix<int, Dynamic, 1> Eigen::VectorXi;

			Options_:
					Eigen::RowMajor		�����ȴ洢��
					ColMajor			�����ȴ洢��Ĭ��ֵ����
					AutoAlign
					DontAlign. 
		*/

		Eigen::MatrixXd m1(2, 2);
		Eigen::MatrixXf mf1(1, 2);
		Eigen::VectorXd v1(3);			// ע����������
		std::cout << v1 << std::endl;

		// ��Array
		ArrayXXd a1(2, 2), a2(2, 2);
		a1 << 1, 2, 3, 4;											// operator << ���Ԫ���������ȵ�˳��
		a2 << 1, 2, 3, 4;
		std::cout << "a1 = \n" << a1 << std::endl;
		std::cout << "a1*a2 = \n" << a1 * a2 << std::endl;

		// �������������Ľӿڡ���LinSpaced(Ԫ��������㣬�յ�)
		int start = 0;
		int end = 10;
		Eigen::VectorXi vi1 = Eigen::VectorXi::LinSpaced(end - start + 1, start, end);
		std::cout << "vi1 == " << vi1 << std::endl;

		Eigen::VectorXf vf1 = Eigen::VectorXf::LinSpaced(5, 0, 10.0);
		std::cout << "vf1 == " << vf1 << std::endl;

		// �����������Ľӿڡ���Random(), Constant(), Ones()....
		Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(3, 3);              // ������static��������Random()���������������
		Eigen::MatrixXd m3 = Eigen::MatrixXd::Constant(3, 3, 1.2);		// ��������ǰ���ǳߴ磬������ֵ��
		Eigen::MatrixXd m4 = Eigen::MatrixXd::Ones(1, 2);					// ȫ1����
 
		// ���ݴ���ջ�ϵľ�������
		Matrix3d mm1 = Matrix3d::Random();
		Vector3d vv1(1, 2, 3);
		std::cout << "m2 = \n" << m2 << std::endl << std::endl;
		std::cout << "mm1 = \n" << mm1 << std::endl << std::endl;

		mm1 = m2;				//		�Ѿ����ջ��������໥��ֵ��
		std::cout << "mm1 = \n" << mm1 << std::endl << std::endl;


		// ��ʼ����

		//			ջ��������ʹ�ô����ų�ʼ��
		Eigen::Vector3f vvv1{ 1,2,3 };
		RowVector3f vvv2{ 4,5,6 };

		std::cout << vvv1 << std::endl;
		std::cout << vvv2 << std::endl;

		//	��ֵ

		//			������������Ժ;�������໥����
		Vector3d vv2(1, 2, 3);

		//			������������Ժ;�������໥��ֵ
		Eigen::MatrixXd mm(v1);
		vv2 = m2.block<3, 1>(0, 0);
		std::cout << "vv2 = \n" << vv2 << std::endl << std::endl;
		mm = vv2;
		std::cout << "mm = \n" << mm << std::endl << std::endl;
	}


	// test1�����������ʡ�Ԫ�ط��ʡ�
	void test1()
	{
		Eigen::MatrixXd m1(3, 4);
		Eigen::VectorXd v1(5);

		// ������������ֵ����������˳�����
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		v1 << 1, 2;
		std::cout << "m1 = \n" << m1 << std::endl << std::endl;
		std::cout << "v1 = \n" << v1 << std::endl << std::endl;

		v1 << 3, 4, 5;			// �����������ĸ�ֵ�����Ǵӵ�һ��Ԫ�ؿ�ʼ��ֵ
		std::cout << "v1 = \n" << v1 << std::endl << std::endl;


		// �±������[]ֻ�ܻ�ȡ����Ԫ�أ���������޷�ʹ�ã���Ϊ[]ֻ֧��һ��������
		std::cout << "v1[0] == " << v1[0] << std::endl << std::endl;

		// �������������Ԫ�أ�ע��������0��ʼ
		std::cout << "m1(0, 1) ==" << m1(0, 1) << std::endl;
		std::cout << "v1(3) == " << v1(3) << std::endl;

		// ���������ʵ����ڽӿ�
		m1 = Eigen::MatrixXd::Random(3, 4);
		Eigen::Matrix<double, Dynamic, Dynamic, Eigen::RowMajor> m1_rm(5, 6);		// �����ȴ洢��
		std::cout << "m1 == \n" << m1 << std::endl;
		std::cout << "Ԫ������m1.size() == " << m1.size() << std::endl;
		std::cout << "������m1.rows() == " << m1.rows() << std::endl;
		std::cout << "������m1.cols() == " << m1.cols() << std::endl;
		std::cout << "��ͣ�sum():       " << m1.sum() << std::endl;
		std::cout << "���ˣ�prod():      " << m1.prod() << std::endl;
		std::cout << "��ֵ��mean():      " << m1.mean() << std::endl;
		std::cout << "��СԪ�أ�minCoeff():  " << m1.minCoeff() << std::endl;
		std::cout << "���Ԫ�أ�maxCoeff():  " << m1.maxCoeff() << std::endl;
		std::cout << "����ļ���trace():     " << m1.trace() << std::endl << std::endl;
		std::cout << "����ʽ��m1.determinant() == " << m1.determinant() << std::endl;

		// outerSize(), innerSize()����Ĭ�ϵ������ȴ洢�ľ���outerSizeΪ�����������ȴ洢��outerSizeΪ������
		std::cout << "m1.outerSize() == " << m1.outerSize() << std::endl;			
		std::cout << "m1.innerSize() == " << m1.innerSize() << std::endl;
		std::cout << "m1_rm.outerSize() == " << m1_rm.outerSize() << std::endl;
		std::cout << "m1_rm.innerSize() == " << m1_rm.innerSize() << std::endl;
 
		// ����󡪡�ʹ��lu�ֽ�õ�
		std::cout << "�����m1.inverse() ==  \n" << m1.inverse() << std::endl;

		// �����ľ���任
		std::cout << "�����ת�ã�transpose() \n" << m1.transpose() << std::endl << std::endl;		// ����ת�õľ��󣬾�������ת��
		std::cout << m1 << std::endl;


		// bug����transpose()�ڶ�����ֵ��ʱ����ʱ����bug��Ҫ��Ҫ��������ı䣬ʹ��~InPlace()
		Eigen::Matrix3f mm;
		mm << 1, 2, 3, 4, 5, 6, 0, 8, 9;
		std::cout << mm << std::endl << std::endl;
		//mm = mm.transpose();
		mm.transposeInPlace();
		std::cout << mm << std::endl << std::endl;
		Eigen::Matrix3f mmCopy = mm;
		Eigen::Matrix3f mm1 = mm.inverse();			// inverse()Ҳһ��������ֻ�ܴ���������������ֵ��
		std::cout << mm1 << std::endl << std::endl;
		std::cout << mm1 * mmCopy << std::endl << std::endl;

		// data()������þ��������������ָ��

		//				���ٱ��������е�Ԫ�أ�
		Eigen::MatrixXf m2(5, 10);
		float* elemPtr = nullptr;
		for (unsigned i = 0; i < m2.size(); ++i)
		{
			elemPtr = m2.data() + i;
			*elemPtr = i;
		}
		std::cout << "m2 == \n" << m2 << std::endl << std::endl;;


		//				�����Ѿ����е����ݣ�
		Eigen::MatrixXf m22(5, 10);
		std::memcpy(m22.data(), m2.data(), m2.size() * sizeof(float));
		std::cout << m22 << std::endl;


		//	 minCoeff() ���� ��������Ԫ���е���ֵ
		std::cout << "m1 == \n" << m1 << std::endl;
		Eigen::MatrixXd::Index maxRow, maxCol;
		Eigen::MatrixXd::Index minRow, minCol;
		double min = m1.minCoeff(&minRow, &minCol);
		double max = m1.maxCoeff(&maxRow, &maxCol);
		std::cout << "��СԪ���к�" << minRow << "���к�" << minCol << std::endl;
		std::cout << "���Ԫ���к�" << maxRow << "���к�" << maxCol << std::endl;

		v1.resize(10);
		v1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
		Eigen::VectorXd::Index minIdx;
		min = v1.minCoeff(&minIdx);
		std::cout << "��������СԪ�أ�" << min << ", ������" << minIdx << std::endl;

		// all(), any(), count()��������Ԫ�ص�����
		ArrayXXf a1(3, 3);
		ArrayXXi a2(2, 2);
		a1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		a2 << 1, 2, 3, 4;
		bool flag1 = (a1 > 0).all();						// all()�����������Ԫ�ض�Ϊtrue���򷵻�true,
		bool flag2 = (a2 < 3).any();					// any()�������������һ��Ԫ��Ϊtrue���򷵻�true
		Eigen::Index count = (a1 > 5).count();		// count()������������������Ԫ����

		Eigen::VectorXi vec(9);
		vec << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		ArrayXi aa = vec.array();
		flag1 = (aa > 2).any();
		std::cout << "flag1" << flag1 << std::endl;

		// ò��ֻ�ܶԵ���Ԫ��ʹ�ã����ܶ�һ��Ԫ��ʹ�á�
		m1.resize(3, 3);
		m1 << 1, 2, 3, 1, 2, 3, 4, 5, 6;
		Eigen::RowVector3d vd1(1, 2, 3);
		std::cout << "m1���Ƿ���������(1,2,3)��" << (m1.array() == vd1.array()).any() << std::endl;
		m1 = Eigen::MatrixXd::Ones(3, 3);
		std::cout << "m1���Ƿ���������(1,2,3)��" << (m1.array() == vd1.array()).any() << std::endl;
	}


	// test2������������任������
	void test2()
	{
		Eigen::MatrixXf m1(3, 4), m3(4, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

		std::cout << "m1 = \n" << std::endl;
		std::cout << m1 << std::endl << std::endl;

		// resize()�������ľ���ߴ硪��ֻ�����ڶѾ��󡢶�����  ���������ܾ���ߴ���ĺ���Σ�����Ԫ��ȫ��������ֵ��
		m3.resize(3, 3);
		m3.row(2) = Eigen::VectorXf::Ones(3);
		std::cout << "m3.size() == " << m3.size() << std::endl;
		std::cout << "m3.rows() == " << m3.rows() << std::endl;
		std::cout << "m3.cols() == " << m3.cols() << std::endl;
		std::cout << "m3 = \n" << m3 << std::endl << std::endl;

		// conservativeResize() ���ľ���ߴ磬����ԭ������
		m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m3.conservativeResize(4, 3);
		m3.row(3) = 5 * Eigen::VectorXf::Ones(3);
		std::cout << "m3 = \n" << m3 << std::endl << std::endl;

		// colwise(), rowwise()���������С����в���
		std::cout << m1.colwise().sum() << std::endl << std::endl;				// ������ͣ�ѹ��һ��������
		std::cout << m1.rowwise().sum() << std::endl << std::endl;				// ������ͣ�ѹ��һ����������

		// �������š���ͨ������˷���������ڿ˻�ʵ�֣�������matlab�е�repmat();
		std::cout << "�������ţ�" << std::endl;
		RowVector3f a(1,2,3);
		Eigen::Vector3f b(4, 5, 6);
		
		auto aa = Eigen::MatrixXf::Ones(3, 1) * a;			//		����������N�� == ���ones(N,1)
		auto bb = b * Eigen::MatrixXf::Ones(1, 4);			//		����������ΪN�� == �ҳ�ones(1,N)
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
		Eigen::MatrixXf result1 = m1.reverse();								// ���ж���װ
		Eigen::MatrixXf result2 = m1.colwise().reverse();				// ������Ϊ��������Ȼ��ÿһ���ڲ����ֲ��䵹װ��
		Eigen::MatrixXf result3 = m1.rowwise().reverse();				// ������Ϊ������....
		std::cout << "result1 = \n" << result1 << std::endl << std::endl;
		std::cout << "result2 = \n" << result2 << std::endl << std::endl;
		std::cout << "result3 = \n" << result3 << std::endl << std::endl;
		m1.colwise().reverseInPlace();
		std::cout << "m1 = \n" << m1 << std::endl << std::endl;

		// ������ˡ����
		Eigen::Vector3f v1(1, 2, 3);
		Eigen::Vector3f v2(3, 2, 1);
		std::cout << "v1.dot(v2) == \n" << v1.dot(v2) << std::endl;
		std::cout << "v1.cross(v2) == " << v1.cross(v2) << std::endl << std::endl;;

		// ������һ��
		v1.normalize();					// normalize()�Ὣ���������һ����normalized()ֻ�Ƿ��ع�һ�����������ı�����
		std::cout << "v1 = \n" << v1 << std::endl << std::endl;
		std::cout << "v2.normalized() = \n" << v2.normalized() << std::endl << std::endl;
		std::cout << "v2 = \n" << v2 << std::endl << std::endl;

		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		Eigen::MatrixXf m11 = m1.rowwise().normalized();			// ÿ����������һ����
		Eigen::MatrixXf m12 = m1.colwise().normalized();				// ÿ����������һ����
		std::cout << "m1 = \n" << m1 << std::endl << std::endl;
		std::cout << "m11 = \n" << m11 << std::endl << std::endl;
		std::cout << "m12 = \n" << m12 << std::endl << std::endl;
		
	}


	// test3����������ؿ�ѧ����Ľӿ�
	void test3()
	{
		Eigen::MatrixXd A(3, 3);
		A << 1, 2, 3, 4, 5, 6, 7, 8, 9;

		// 1. ����������ֵ������������
		EigenSolver<Matrix3d> es(A);
		Matrix3d D = es.pseudoEigenvalueMatrix();			// �Խ���Ԫ��������ֵ
		Matrix3d V = es.pseudoEigenvectors();					// ÿһ����������������������
		std::cout << "����ֵ����D��" << std::endl << D << std::endl;
		std::cout << "������������V: " << std::endl << V << std::endl;
		std::cout << "Finally, V * D * V^(-1) = " << std::endl << V * D * V.inverse() << std::endl;

		std::cout << A * V.block<3, 1>(0, 0) << std::endl << std::endl;
		std::cout << D(0, 0) * V.block<3, 1>(0, 0) << std::endl << std::endl;
		debugDisp("\n\n");


		// 2. �����LU�ֽ⡪��Eigen::FullPivLU
		Eigen::MatrixXf m = Eigen::MatrixXf::Random(5, 5);
		std::cout << "Here is the matrix m:" << std::endl << m << std::endl;
		Eigen::FullPivLU<Eigen::MatrixXf> lu(m);
		std::cout << "ִ��LU�ֽ⣺lu.matrixLU() == :"
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


		// 3. ���������ֵ�ֽ⣨SVD������
		m = Eigen::MatrixXf::Random(3, 3);
		std::cout << "m == \n" << m << std::endl;
		JacobiSVD<Eigen::MatrixXf> svdSolver(m, ComputeThinU | ComputeThinV);
		std::cout << "����ֵ == " << std::endl << svdSolver.singularValues() << std::endl;
		std::cout << "U == " << std::endl << svdSolver.matrixU() << std::endl;
		std::cout << "V == " << std::endl << svdSolver.matrixV() << std::endl;
		debugDisp("\n");

		//		3.1 ͨ��SVD������Է����飺
		Eigen::Vector3f rhs(1, 0, 0);
		std::cout << "ʹ������ֵ�ֽ�����Է����飺" << std::endl << svdSolver.solve(rhs) << std::endl << std::endl;
		debugDisp("\n");

		//		3.2 ͨ��SVD��������������
		double condNum = calcCondNum(m);
		debugDisp("������������condNum == ", condNum);
		debugDisp("\n\n");

		// 4. ���ȣ�
		std::cout << "svdSolver.rank == " << svdSolver.rank() << std::endl << std::endl;

		// 5. �����Լ���װ��ʹ��QR�ֽ�����ǡ���������Է���������ӣ�
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
		std::cout << "���Է�����Ax == b�Ľ⣺" << std::endl << x << std::endl << std::endl;
		std::cout << "���Է�����AX == B�Ľ⣺" << std::endl << X << std::endl << std::endl;
		debugDisp("\n\n");

		// 6. norm()����ŷ����÷�����Ҳ��p==2ʱ��lp������������Ԫ��ƽ���͵Ŀ�����
		Eigen::Vector3f v1(1, 2, 3);
		std::cout << "v1.norm() == " << v1.norm() << std::endl;	// �����ķ�������������ģ����
		m.resize(3, 4);
		m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 , 11, 12;
		std::cout << "m == \n" << m << std::endl;
		std::cout << "m.norm() == " << m.norm() << std::endl;
		std::cout << "m.rowwise().norm() == \n" << m.rowwise().norm() << std::endl;		// ������������norm������һ����������
		std::cout << "m.colwise().norm() == \n" << m.colwise().norm() << std::endl;

		//	6.1. ͨ����������������ģ����
		Eigen::Matrix3f arrows;
		arrows << 0, 1, 0, 0, 3, 4, 0, 5, 12;
		Eigen::Vector3f arrowsLen = arrows.rowwise().norm();
		dispVec<float>(arrowsLen);
	}


	// test4��������������
	void test4()
	{
		Eigen::MatrixXf m1(5, 6);
		float* pdata = m1.data();
		for (unsigned i = 0; i < m1.size(); ++i)
			pdata[i] = static_cast<float>(i);
		std::cout << "m1 == \n" << std::endl;
		dispMat(m1);


		// setConstant()���������鸳ֵ��
		m1.setConstant(1.0);
		m1.topRightCorner(2, 3).setConstant(2.0);
		dispMat(m1);


		// VectorXT::segment<>()����������ȡ������������Ƭ�Σ���������ֵ���á���������
		Eigen::VectorXi v1(9);
		v1 << 1, 2, 3, 4, 5, 6, 7, 8, 8;

		//			segment()����1����segment<elemCount>( startIdx)����̬segment��Ƭ�γ���elemCount�����ڱ�������֪����constExpr
		std::cout << "v1.segment<3>(1)  == \n" << v1.segment<3>(1) << std::endl;

		//			segment()����2����segment(startIdx, elemCount)����̬segment��Ƭ�γ���elemCount�����ڱ�������δ֪��
		std::cout << "v1.segment(1, v1(2)) == \n" << v1.segment(1, v1(2)) << std::endl;

		//			segment()���ص�����Ƭ������ֵ���ã�
		v1.segment<5>(2) = 999 * Eigen::VectorXi::Ones(5);
		std::cout << "v1 == \n" << v1 << std::endl << std::endl;


		// block()����������ȡ�Ӿ���飬��������ֵ���á�����������
		std::cout << "block()����" << std::endl;

		//			����1��block<rows, cols>(startRow, startCol)������̬block�������ĳߴ�����Ǳ�������֪�ģ�������constExpr
		std::cout << "m1.block<2, 2>(1, 1)  == \n" << m1.block<2, 2>(1, 1) << std::endl << std::endl;

		//			����2: block(startRow, startCol, rows, cols)������̬block�������ĳߴ�����Ǳ�����δ֪��
		std::cout << "m1.block(1, 0, m1(2, 0), m1(3, 0)) == \n" << m1.block(1, 0, m1(2, 0), m1(3, 0)) << std::endl << std::endl;


		//			 ���صľ��������ֵ���á����Ա�<<��ֵ��
		m1.block<1, 4>(1, 0) << 88, 99, 111, 222;
		std::cout << "m1 = \n" << m1 << std::endl << std::endl;
		m1.block<1, 4>(1, 0) << 5, 6, 7, 8;


		//			���صľ�������ÿ���ֱ����' = ' ����һ������ֵ��
		Eigen::MatrixXf m2 = Eigen::MatrixXf::Ones(1, 6);
		m1.block<1, 6>(0, 0) = m2;
		std::cout << "m1 == \n" << m1 << std::endl;


		// row(), col()������ȡĳһ�л���ĳһ��
		std::cout << "row(), col()������ȡĳһ�л���ĳһ��\n m1.col(3) == \n" << m1.col(3) << std::endl;

		//			��ȡ�����С�����ԭ���ݵ����ã������µĿ�������������ֵ��
		m1.col(3) << 1, 2, 3, 4, 5;
		std::cout << "m1.col(3) == \n" << m1.col(3) << std::endl;

		// topRows(), bottomRows(), leftCols(), rightCols();
		std::cout << "m1.leftCols(2) ==   \n" << m1.leftCols(2) << std::endl;
		std::cout << "m1.topRows(2) == \n" << m1.topRows(2) << std::endl;
 

		//	minCoeff()����ֵ�Ծ���ֿ�Ҳͬ������
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
		Eigen::MatrixXf m1(5, 6);
		float* pdata = m1.data();
		for (unsigned i = 0; i < m1.size(); ++i)
		{
			pdata[i] = static_cast<float>(i);
		}
		std::cout << "m1 == \n" << m1 << std::endl;

		// ��Ԫ�ز�����Ҫʹ��array()
		Eigen::MatrixXf result = m1.array() * Eigen::MatrixXf::Ones(5, 6).array();
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


		Eigen::MatrixXf m2 = Eigen::MatrixXf::Ones(5, 6);
		auto result2 = (m1.array() < m2.array());
		std::cout << "������MATLAB�е��߼����� \n" << result2 << std::endl << std::endl;
		std::cout << typeid(result2).name() << std::endl;
	}


	// test6���� ���;���
	void test6()
	{
		Eigen::MatrixXf m1 = Eigen::MatrixXf::Ones(5, 6);
		Matrix3i m2 = Matrix3i::Random();
		dispMat(m1);
		dispMat(m2);
	}


	// test7����flag����
	void test7()
	{
		// eigen 3.3.7��û��֧��flag���󣬵��ǿ���ʹ��select()������ʵ�����Ƶ�Ч����
		Matrix3i m1;
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

		//		(array�����ж�).select(����a,����b)��array���±�Ϊij��Ԫ�������ж�Ϊ�棬�򷵻ؾ���ijλ��Ԫ��Ϊ����a��ijλ��Ԫ�أ�����Ϊ����b��ijλ��Ԫ�ء�
		Matrix3i flag1 = (m1.array() > 4).select(Matrix3i::Ones(), Matrix3i::Zero());
		std::cout << "m1 == \n" << m1 << std::endl << std::endl;
		std::cout << "flag1 == \n" << flag1 << std::endl << std::endl;
		std::cout << std::endl;

		// rowInMat()�������������Ƿ����ĳһ��������������һ��������������
		std::cout << "���������Ƿ����ĳһ��������������һ��������������" << std::endl;
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


	// test8����Eigen::Map��
	void test8()
	{
		// �����reshape����ʹ��Eigen::Map��ʵ�֡�
		Eigen::MatrixXf	m1(4, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
		m1.transposeInPlace();
		Eigen::Map<Eigen::MatrixXf>  m11(m1.data(), 8, 2);			// reshapeʱԪ���ǰ��洢˳���ȡ�ģ�Ĭ�ϼ����л�ȡ��
		std::cout << "reshape֮ǰ��m1 == \n" << m1 << std::endl << std::endl;
		std::cout << "reshape֮��m11 == \n" << m11 << std::endl << std::endl;
		std::cout << "typeid(m11).name()  == " << typeid(m11).name() << std::endl;

		Eigen::MatrixXf m111{ Eigen::Map<Eigen::MatrixXf>(m1.data(), 1, 16) };
		std::cout << "m111 == \n" << m111 << std::endl << std::endl;
		std::cout << "typeid(m111).name()  == " << typeid(m111).name() << std::endl;

		Eigen::VectorXf v1 = Eigen::Map<Eigen::VectorXf>(m1.data(), 16, 1);
		dispVec(v1);

		std::cout << "finished." << std::endl;
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

		// ת����Ŀ��������
		Eigen::Matrix3f rotation = getRotationMat(RowVector3f(1, 0, 0), RowVector3f(0, 0, 1));
		RowVector3f v1 = (rotation * RowVector3f(1, 0, 0).transpose()).transpose();
		std::cout << "��ת���v1 == " << v1 << std::endl;
		rotation = getRotationMat(RowVector3f(1, 1, 1), RowVector3f(3, 4, 5));
		v1 = (rotation * RowVector3f(1, 1, 1).transpose()).transpose();
		float scale = 5.0 / v1(2);
		v1 *= scale;
		std::cout << "��ת���v1 == " << v1 << std::endl;

		// eigen�в�֧����Ԫ���ӷ���ֻ���Լ����������ӣ�
		q3.w() = q1.w() + q2.w();
		q3.vec() = q1.vec() + q2.vec();
		std::cout << std::endl << std::endl;
		dispQuat(q3);


		// Eigen::AngleAxisf�ࡪ��һ��������ʾһ������ĳ�����һ���Ƕȵ���ת�������ڹ�����ת����
		const float pi = 3.14159;
		Eigen::Matrix3f m;
		m = Eigen::AngleAxisf(0.25 * pi, Eigen::Vector3f::UnitX())
			* Eigen::AngleAxisf(0.5 * pi, Eigen::Vector3f::UnitY())
			* Eigen::AngleAxisf(0.33 * pi, Eigen::Vector3f::UnitZ());				// ��ʾ�Ⱥ���xyz��������תһ���Ƕȡ�
		std::cout << m << std::endl << "is unitary: " << m.isUnitary() << std::endl;

		// ƽ��(translation)�任


		// ����(scaling)

	}


	// test10�����ռ�任Ӧ�ã�
	void test10()
	{
		// ��֪��������ת�任ǰ������꣬�������ת�任����ת����
		Eigen::MatrixXf vers, gumline, axis;
		Eigen::MatrixXi tris;

		objReadMeshMat(vers, tris, "E:/����/patientTooth11.obj");
		objReadVerticesMat(gumline, "E:/����/gumline11.obj");
		objReadVerticesMat(axis, "E:/����/axisPatient11.obj");
		
		// ����������Ϊԭ�㣺
		RowVector3f bary = vers.colwise().mean();
		vers.rowwise() -= bary;
		gumline.rowwise() -= bary;

		objWriteMeshMat("E:/tooth_ori.obj", vers, tris);
		objWriteVerticesMat("E:/gumline_ori.obj", gumline);

		objWriteCoorSys("E:/coordinate_ori.obj", RowVector3f::Zero(), axis.row(0), axis.row(1), axis.row(2));

		// ԭ��������ϵ��Ӧ��һ�����Կռ�omega1�����᷽��������ɵľ�����һ����λ����M1 == I;
		Eigen::Matrix3f M1 = Eigen::Matrix3f::Identity();

		// Ŀ������ϵ��Ӧ�������Կռ�omega2�����᷽��������ɵľ���ΪM2 == axis.transpose();
		Eigen::Matrix3f M2 = axis.transpose();

		// ��ת�任rotation���Խ����Կռ�omega1�任��omega2��rotation * M1 == M2; �� rotation == M2;
		Eigen::Matrix3f rotation = M2;

		// ���Կռ�omega2��ص�omega1�ı任����Ϊrotation�������rotationIev == rotation.inverse();
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


	// test11����ʹ��eval()������ʱ��������
	void test11()
	{
		Eigen::MatrixXf m1(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m1 = m1.transpose();						// ʵ����transpose()����������ֵ���������õ��Ľ���޷�Ԥ�⡣
		std::cout << "m1 == \n" << m1 << std::endl;

		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m1 = m1.transpose().eval();			// ʹ��eval()��ת�õĽ���洢����ʱ�����С�
		std::cout << "m1 == \n" << m1 << std::endl;

		m1 = 2 * Eigen::MatrixXf::Ones(3, 3);
		Eigen::MatrixXf m2(3, 2);
		m2 << 1, 2, 3, 4, 5, 6;
		m2 = m1 * m2;							// ��ʱ�����ֻ���������
		std::cout << "m2 == \n" << m2 << std::endl;

		m1 = 2 * Eigen::MatrixXf::Ones(3, 3);
		m2 << 1, 2, 3, 4, 5, 6;
		m2 = (m1 * m2).eval();				// ��������ȷ��������ֻ���������
		std::cout << "m2 == \n" << m2 << std::endl;


	}


	// test12����map�ࡪ�����������е����ݣ�����ӳ������������
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

		// map�����������Ը�ֵ����������󣬷������ݵĿ�����
		Eigen::VectorXi v1 = mv1;
		Eigen::MatrixXi m1 = mm1;
		Eigen::MatrixXi m2 = Map<Matrix<int, 4, 3, Eigen::RowMajor>>(intArr);		// �������Կ��ٿ��������е����ݣ���ֵ���´����ľ�����
		dispVec<int>(v1);
		dispMat(m1);
		dispMat(m2);

		std::cout << "������Ԫ�ص�ַ��" << reinterpret_cast<size_t>(&intArr[0]) << std::endl;
		std::cout << "mm1��Ԫ�ص�ַ��" << reinterpret_cast<size_t>(&mm1(0, 0)) << std::endl;
		std::cout << "m1��Ԫ�ص�ַ��" << reinterpret_cast<size_t>(&m1(0, 0)) << std::endl << std::endl;

		// ӳ���������ݱ��ı䣬ԭ���������Ҳ�ı�
		mm1.setRandom();
		std::cout << "mm1 == \n" << mm1 << std::endl;
		for (const auto& num : intArr)
			std::cout << num << ", ";
		std::cout << std::endl << std::endl;

		// ��ϣ���ı�ԭ�����е����ݵĻ�������ӳ��ɳ�������
		Map<const Matrix<int, 3, 4>> cmm1(intArr);				// ������ֵ

	}


	// test13��������ת����
	void test13() 
	{
		// cast()������
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


	template<typename Derived>											// Derived�Ǿ���ľ������ͣ���Eigen::Matrix<int,1,-1,1,1,-1>
	void foo(const Eigen::MatrixBase<Derived>& base)			// Eigen::MatrixBase<>��������ܾ��������Ļ��ࣻ
	{
		const Derived& originMat = base.derived();			// derived()����ԭ��������ã�����Ϊԭ���ͣ�
		std::cout << "typeid(base).name() == " << typeid(base).name() << std::endl;
		std::cout << "typeid(base.derived()).name()  == " << typeid(originMat).name() << std::endl;
		std::cout << "typeid(Derived).name() == " << typeid(Derived).name() << std::endl;

		std::cout << "address of base: " << reinterpret_cast<size_t>(&base) << std::endl;
		std::cout << "address of base.derived(): " << reinterpret_cast<size_t>(&originMat) << std::endl;
		std::cout << std::endl;
	}

	// test14����dense mat�����νṹ��
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


	// ��ʱ�޷����ࣺ
	void test000()
	{
		// ����������
		std::vector<int> rowVec{ 1, 3, 4 };
		std::vector<int> colVec{ 5, 2, 0 };

		std::vector<int> ind{ 4,2,5,5,3 };
		Eigen::MatrixXi A = Eigen::MatrixXi::Random(4, 6);
		std::cout << "Initial matrix A:\n" << A << "\n\n";

		// �鿴eigen��İ汾��
		std::cout << EIGEN_WORLD_VERSION << std::endl;
		std::cout << EIGEN_MAJOR_VERSION << std::endl;
		std::cout << EIGEN_MINOR_VERSION << std::endl;		// �汾Ϊ3.3.7;
	}

}
