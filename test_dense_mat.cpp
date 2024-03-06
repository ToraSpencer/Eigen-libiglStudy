#include "test_dense_mat.h"


#define MAXLEN 1024

using namespace Eigen;
using spMatf = Eigen::SparseMatrix<float, ColMajor>;
using TripF = Eigen::Triplet<float>;


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


////////////////////////////////////////////////////////////////////////////////////////////// ����eigen�еĳ��ܾ���
namespace TEST_DENSE_MAT
{
	// ��ʱ�޷����ࣺ
	void test000()
	{
		// 1. ����������
		std::vector<int> rowVec{ 1, 3, 4 };
		std::vector<int> colVec{ 5, 2, 0 };

		std::vector<int> ind{ 4,2,5,5,3 };
		Eigen::MatrixXi A = Eigen::MatrixXi::Random(4, 6);
		std::cout << "Initial matrix A:\n" << A << "\n\n";

		// 2. �鿴eigen��İ汾��
		std::cout << EIGEN_WORLD_VERSION << std::endl;
		std::cout << EIGEN_MAJOR_VERSION << std::endl;
		std::cout << EIGEN_MINOR_VERSION << std::endl;		// �汾Ϊ3.4.0

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


	// test0����eigen��Ļ������ݽṹ
	void test0()
	{
		// �Ѿ�����������ȷ���˳ߴ磬��δ��ʼ��,���ݴ��ڶ���
		/*
			����ģ�塪��Matrix<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
			�Ѿ��󡪡�typedef Matrix<double, Dynamic, Dynamic> Eigen::MatrixXd;
			����������typedef Matrix<int, Dynamic, 1> Eigen::VectorXi;

			�̳й�ϵ�� 
				Matrix<> �� PlainObjectBase<>
				MatrixBase<> �� DenseBase<>

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

		// 1. ��Array
		ArrayXXd a1(2, 2), a2(2, 2);
		a1 << 1, 2, 3, 4;											// operator << ���Ԫ���������ȵ�˳��
		a2 << 1, 2, 3, 4;
		std::cout << "a1 = \n" << a1 << std::endl;
		std::cout << "a1*a2 = \n" << a1 * a2 << std::endl;

		// 2. �������������Ľӿڡ���LinSpaced(Ԫ��������㣬�յ�)
		int start = 0;
		int end = 10;
		Eigen::VectorXi vi1 = Eigen::VectorXi::LinSpaced(end - start + 1, start, end);
		std::cout << "vi1 == " << vi1 << std::endl;

		Eigen::VectorXf vf1 = Eigen::VectorXf::LinSpaced(5, 0, 10.0);
		std::cout << "vf1 == " << vf1 << std::endl;

		// 3. �����������Ľӿڡ���Random(), Constant(), Ones()....
		Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(3, 3);              // ������static��������Random()���������������
		Eigen::MatrixXd m3 = Eigen::MatrixXd::Constant(3, 3, 1.2);		// ��������ǰ���ǳߴ磬������ֵ��
		Eigen::MatrixXd m4 = Eigen::MatrixXd::Ones(1, 2);					// ȫ1����
 
		// 4. ���ݴ���ջ�ϵľ�������
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

		//1.  ������������ֵ����������˳�����
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		v1 << 1, 2;
		std::cout << "m1 = \n" << m1 << std::endl << std::endl;
		std::cout << "v1 = \n" << v1 << std::endl << std::endl;

		v1 << 3, 4, 5;			// �����������ĸ�ֵ�����Ǵӵ�һ��Ԫ�ؿ�ʼ��ֵ
		std::cout << "v1 = \n" << v1 << std::endl << std::endl;

		// 2. ����Ԫ�ص�������

		//		�±������[]ֻ�ܻ�ȡ����Ԫ�أ���������޷�ʹ�ã���Ϊ[]ֻ֧��һ��������
		std::cout << "v1[0] == " << v1[0] << std::endl << std::endl;

		//		�������������Ԫ�أ�ע��������0��ʼ
		std::cout << "m1(0, 1) ==" << m1(0, 1) << std::endl;
		std::cout << "v1(3) == " << v1(3) << std::endl;


		// 3. ���������ʵ����ڽӿ�
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

		//			outerSize(), innerSize()����Ĭ�ϵ������ȴ洢�ľ���outerSizeΪ�����������ȴ洢��outerSizeΪ������
		std::cout << "m1.outerSize() == " << m1.outerSize() << std::endl;			
		std::cout << "m1.innerSize() == " << m1.innerSize() << std::endl;
		std::cout << "m1_rm.outerSize() == " << m1_rm.outerSize() << std::endl;
		std::cout << "m1_rm.innerSize() == " << m1_rm.innerSize() << std::endl;
 
		//			����󡪡�ʹ��lu�ֽ�õ�
		std::cout << "�����m1.inverse() ==  \n" << m1.inverse() << std::endl;

		//			ת��
		std::cout << "�����ת�ã�transpose() \n" << m1.transpose() << std::endl << std::endl;		// ����ת�õľ��󣬾�������ת��
		std::cout << m1 << std::endl;

		//			bug����transpose()�ڶ�����ֵ��ʱ����ʱ����bug��Ҫ��Ҫ��������ı䣬ʹ��~InPlace()
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

		//			norm()����ŷ����÷�����Ҳ��p==2ʱ��lp������������Ԫ��ƽ���͵Ŀ�����
		debugDisp("v1 == \n", v1, "\n");
		debugDisp("v1.norm() == ", v1.norm()); 		// �����ķ�������������ģ����
		m1.resize(3, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
		std::cout << "m1 == \n" << m1 << std::endl;
		std::cout << "m1.norm() == " << m1.norm() << std::endl;
		std::cout << "m1.rowwise().norm() == \n" << m1.rowwise().norm() << std::endl;		// ������������norm������һ����������
		std::cout << "m1.colwise().norm() == \n" << m1.colwise().norm() << std::endl;

		//		ͨ����������������ģ����
		Eigen::Matrix3f arrows;
		arrows << 0, 1, 0, 0, 3, 4, 0, 5, 12;
		Eigen::Vector3f arrowsLen = arrows.rowwise().norm();
		dispVec<float>(arrowsLen);

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
		Eigen::VectorXd::Index minIdx, maxIdx;
		min = v1.minCoeff(&minIdx);
		debugDisp("��������СԪ�أ�min == " , min , ", ����minIdx == ��" , minIdx); 

		v1 <<  -INFINITY, 2, 3, 4, 5, 6, 7, 8, INFINITY, 0;
		min = v1.minCoeff(&minIdx);					// Ԫ���п�����INFINITY�����ǲ�������NAN����������ʱ���׳��쳣��
		max = v1.maxCoeff(&maxIdx);
		debugDisp("��������СԪ�أ�min == ", min, ", ����minIdx == ��", minIdx);
		debugDisp("���������Ԫ�أ�max == ", max, ", ����maxIdx == ��", maxIdx);

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

		Eigen::MatrixXd aaa, bbbb;
		aaa = kron(RowVector3f::Ones(), a); 
		bbbb = kron(Matrix2i::Ones(), b); 
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


	// test3�������ܾ���ķֽ⡢��������Է����飺
	void test3()
	{
		// ���ַֽⷽ�������Է���������ܱȽϣ�
		/*
			LU�ֽ⣺		
					Ч�ʸߣ�
					�����ڳ��ܶȸߵľ���

			LDLT�ֽ⣺	
					�����ڶѳ���������

			QR�ֽ⣺		
					�ڴ����������ʱ��ֵ�ȶ��Ժã�
					�������ʱЧ�ʸ��ߣ�
		
			����ֵ�ֽ⣨SVD��
					��������ϡ��ġ�����ֵ�ֽ�ĵ��Ƚ����ܹ��ṩ�㹻�ľ��ȣ�
							���Ҿ�����������ʺ������ֽⷽ������ô����ֵ�ֽ������һ�����ʵ�ѡ��
		*/
		Eigen::MatrixXd A(3, 3);
		Vector3d b{1, 2, 3};
		Vector3d x;
		A << 1, 22, 3, 4, 55, 6, 77, 8, 9;
		debugDisp("A == \n", A, "\n"); 

		// 1. �����ֽ⣺
		EigenSolver<Matrix3d> solverEigen(A);
		{
			Matrix3d D = solverEigen.pseudoEigenvalueMatrix();			// �Խ���Ԫ��������ֵ
			Matrix3d V = solverEigen.pseudoEigenvectors();					// ÿһ����������������������
			debugDisp("����ֵ����D == \n", D, "\n");
			debugDisp("����ֵ����V == \n", V, "\n");
			debugDisp("V * D * V.inverse() == \n", V * D * V.inverse(), "\n");
			debugDisp("A * V.block<3, 1>(0, 0) == \n", A * V.block<3, 1>(0, 0), "\n");
			debugDisp("D(0, 0) * V.block<3, 1>(0, 0) == \n", D(0, 0) * V.block<3, 1>(0, 0), "\n");
			debugDisp("\n\n");
		}

		// 2. �����LU�ֽ⡪��Eigen::FullPivLU
		Eigen::FullPivLU<Eigen::MatrixXd> solverLU;
		{
			Eigen::MatrixXd L, U, P_inv, Q_inv;
			solverLU.compute(A);
			debugDisp("ִ��LU�ֽ⣺solverLU.matrixLU() == \n", solverLU.matrixLU());
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

			// LU�ֽ�����Է����飻
			x = solverLU.solve(b);
			debugDisp("���Է�����Ax == b�Ľ⣺x == \n", x, "\n");
			debugDisp("x == \n", x, "\n");
			debugDisp("A * x == \n", A * x, "\n");

			debugDisp("\n\n");
		}

		// 3. ���������ֵ�ֽ⣨SVD������  
		JacobiSVD<Eigen::MatrixXd> solverSVD(A, ComputeThinU | ComputeThinV);
		{
			debugDisp("����ֵ��solverSVD.singularValues() == \n ", solverSVD.singularValues(), "\n");
			debugDisp("solverSVD.matrixU() == \n ", solverSVD.matrixU(), "\n");
			debugDisp("solverSVD.matrixV() == \n ", solverSVD.matrixV(), "\n");
			debugDisp("\n");

			//		3.1 ͨ��SVD������Է����飺 
			x = solverSVD.solve(b);
			debugDisp("ʹ������ֵ�ֽ�����Է����飺b == \n", b, "\n");
			debugDisp("A * x == \n", A * x, "\n");
			debugDisp("\n");

			//		3.2 ͨ��SVD��������������
			double condNum = calcCondNum(A);
			debugDisp("������������condNum == ", condNum);
			debugDisp("\n\n");

			//		3.3 ���ȡ�����������û�����ȵķ�����ֻ��ͨ��SVD�ֽ������ȡ�
			debugDisp("solverSVD.rank() == ", solverSVD.rank());  
		}

		// 4. LLT�ֽ⡢LDLT�ֽ⣺
		Eigen::LLT<Eigen::MatrixXd> solverLLT;
		Eigen::LDLT<Eigen::MatrixXd> solverLDLT;
		{
			Eigen::MatrixXd A1(2, 2); 
			A1 << 2, -1, -1, 3;							// ��������֪��Ϊʲô�������ľ���Ҳ������LLT�ֽ⣻
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
			debugDisp("LLT�ֽ⣺L_llt == \n", L_llt, "\n");
		}

		// 5. QR�ֽ⡪������householder�任�� 
		Eigen::HouseholderQR<Eigen::MatrixXd> solverQR1;				// û��rank()�������ȶ��������ھ�����������
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solverQR2;		// ��rank()�������ȶ��Բ���
		{
			//	5.1 �����Լ���װ��ʹ��QR�ֽ�����ǡ���������Է���������ӣ� 
			Eigen::MatrixXd B = Eigen::MatrixXd::Random(3, 3);
			Eigen::MatrixXd X;
			solveLinearEquation(x, A, b);
			solveLinearEquations<double>(X, A, B);
			debugDisp("���Է�����Ax == b�Ľ⣺x == \n", x, "\n");
			debugDisp("x == \n", x, "\n");
			debugDisp("A * x == \n", A * x, "\n");
			debugDisp("B == \n", B, "\n");
			debugDisp("���Է�����AX == B�Ľ⣺X == \n", X, "\n");
			debugDisp("A * X == \n", A * X, "\n");

			// 5. 2 Eigen::HouseholderQR<>
			A.resize(2, 3);
			A << 1, 2, 3, 4, 5, 6;
			solverQR1.compute(A);			// A = Q* R;
			solverQR2.compute(A);			// A * P = Q * R;	����P���û�����
			if (solverQR2.info() != Eigen::Success)
			{
				debugDisp("ColPivHouseholderQR decomposition failed");
				return;
			}
			Eigen::MatrixXd Q = solverQR1.householderQ();
			Eigen::MatrixXd R = solverQR1.matrixQR().triangularView<Eigen::Upper>();
			debugDisp("\n\nQR�ֽ⣺Eigen::HouseholderQR<>");
			debugDisp("A == \n", A, "\n");
			debugDisp("Q == \n", Q, "\n");
			debugDisp("R == \n", R, "\n");
			debugDisp("Q * R == \n", Q * R, "\n");

			// 5.3 Eigen::ColPivHouseholderQR<>
			Q = solverQR2.householderQ();
			R = solverQR2.matrixQR().triangularView<Eigen::Upper>();
			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P = solverQR2.colsPermutation();			// �û�����
			Eigen::MatrixXd Pd = P;								// Eigen::PermutationMatrixͨ����ֵ�������ת��Ϊ��ͨ�������ͣ�û��cast������
			Eigen::MatrixXd PdT = Pd.transpose();		// �û����������������������ת�ã�
			debugDisp("\n\nQR�ֽ⣺Eigen::ColPivHouseholderQR<>");
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


	// test4��������������
	void test4()
	{
		Eigen::MatrixXf m1(5, 6);
		float* pdata = m1.data();
		for (unsigned i = 0; i < m1.size(); ++i)
			pdata[i] = static_cast<float>(i);
		std::cout << "m1 == \n" << std::endl;
		dispMat(m1);

		// 1. setConstant()���������鸳ֵ��
		m1.setConstant(1.0);
		m1.topRightCorner(2, 3).setConstant(2.0);
		dispMat(m1);

		//	2. VectorXT::segment<>()����������ȡ������������Ƭ�Σ���������ֵ���á���������
		Eigen::VectorXi v1(9);
		v1 << 1, 2, 3, 4, 5, 6, 7, 8, 8;

		//			segment()����1����segment<elemCount>( startIdx)����̬segment��Ƭ�γ���elemCount�����ڱ�������֪����constExpr
		std::cout << "v1.segment<3>(1)  == \n" << v1.segment<3>(1) << std::endl;

		//			segment()����2����segment(startIdx, elemCount)����̬segment��Ƭ�γ���elemCount�����ڱ�������δ֪��
		std::cout << "v1.segment(1, v1(2)) == \n" << v1.segment(1, v1(2)) << std::endl;

		//			segment()���ص�����Ƭ������ֵ���ã�
		v1.segment<5>(2) = 999 * Eigen::VectorXi::Ones(5);
		std::cout << "v1 == \n" << v1 << std::endl << std::endl;

		// 2+. ������head(), tail()����������Ϊ�Ͽ��Ƿ�������Ƭ�ε����ã�
		v1.head(3) = Eigen::Vector3i{-1, -1, -1};
		v1.tail(2) = Eigen::RowVector2i{88, 88};
		debugDisp("v1 == \n", v1);

		// 3. block()����������ȡ�Ӿ���飬��������ֵ���á�����������
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


		// 4. row(), col()������ȡĳһ�л���ĳһ��
		std::cout << "row(), col()������ȡĳһ�л���ĳһ��\n m1.col(3) == \n" << m1.col(3) << std::endl;

		//			��ȡ�����С�����ԭ���ݵ����ã������µĿ�������������ֵ��
		m1.col(3) << 1, 2, 3, 4, 5;
		std::cout << "m1.col(3) == \n" << m1.col(3) << std::endl;

		// 5. topRows(), bottomRows(), leftCols(), rightCols();
		std::cout << "m1.leftCols(2) ==   \n" << m1.leftCols(2) << std::endl;
		std::cout << "m1.topRows(2) == \n" << m1.topRows(2) << std::endl;
 
		//	6. minCoeff()����ֵ�Ծ���ֿ�Ҳͬ������
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

		// 7. rowwise(), colwise()�Ծ������С��еĲ�����
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
			pdata[i] = static_cast<float>(i); 
		debugDisp("m1 == \n", m1, "\n"); 

		// 1. �������array()�����������ؾ����������鲿�֣�����
		Eigen::MatrixXi m111{ Eigen::MatrixXi::Ones(3, 4) };
		auto retA = m111.array();						// ��������ΪEigen::ArrayWrapper<>
		debugDisp("typeid(aaa).name() == ", typeid(retA).name());
		debugDisp("retA == \n", retA);

		m111(0, 0) = 999; 
		debugDisp("retA == \n", retA);
		debugDisp("\n\n");

		// ��Ԫ�ز�����Ҫʹ��array()
		Eigen::VectorXf v1(6), v2(6);
		v1 << 1, 2, 3, 4, 5, 6;
		v2 << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
		Eigen::VectorXf result1 = v1.array() * v2.array();
		Eigen::VectorXf result2 = v1.array() / v2.array();
		debugDisp("��Ԫ�����㣺 \n result1 == \n", result1, "\n"); 
		debugDisp( "result2 == \n", result2, "\n");

		// array��sqrt(), pow(), isX()������ 
		Eigen::MatrixXf result;
		result = m1.array().sqrt();
		debugDisp("��Ԫ�ؿ����� \n", result, "\n"); 

		result = m1.array().pow(2);
		debugDisp("��Ԫ��ƽ���� \n", result, "\n"); 

		m1.col(0).array() = -1;
		m1.col(1).array() = 0;
		m1.col(2).array() = 1;
		m1.rightCols(3).array() = 2;
		debugDisp(".array() = ������������Ԫ�ظ�ֵ��\n", m1, "\n"); 

		Eigen::MatrixXf m2 = Eigen::MatrixXf::Ones(5, 6);
		auto result3 = (m1.array() < m2.array());
		debugDisp("������MATLAB�е��߼����� \n", result3, "\n"); 
		debugDisp("typeid(result3).name() == ", typeid(result3).name());

		debugDisp("finished.");
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
		Eigen::MatrixXi m1(3, 3);
		Eigen::MatrixXf m11(3, 3);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m11 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

		//	1. (m1��array�����ж�).select(����a,����b)������m1�±�Ϊij��Ԫ�������ж�Ϊ�棬�򷵻ؾ���a��ijλ��Ԫ�أ�����Ϊ����b��ijλ��Ԫ�ء�
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
		debugDisp("����֮����Ԫ�رȽ������жϣ�");
		debugDisp("m11 == \n", m11, "\n");
		debugDisp("m111 == \n", m111, "\n");
		m11 = (m11.array() < m111.array()).select(-Eigen::Matrix4f::Ones(), m11);
		debugDisp("m11 == \n", m11, "\n\n");

		// 2. rowInMat()�������������Ƿ����ĳһ��������������һ��������������
		debugDisp("���������Ƿ����ĳһ��������������һ��������������" );
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
		Eigen::MatrixXf	m1(4, 4);
		m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
		m1.transposeInPlace();

		// 1. �����reshape����ʹ��Eigen::Map��ʵ�֡�
		Eigen::Map<Eigen::MatrixXf>  map1(m1.data(), 8, 2);				// reshapeʱԪ���ǰ��洢˳���ȡ�ģ�Ĭ�ϼ����л�ȡ��
		std::cout << "reshape֮ǰ��m1 == \n" << m1 << std::endl << std::endl;
		std::cout << "reshape֮��m11 == \n" << map1 << std::endl << std::endl;
		std::cout << "typeid(map1).name()  == " << typeid(map1).name() << std::endl;

		Eigen::MatrixXf m11{ Eigen::Map<Eigen::MatrixXf>(m1.data(), 1, 16) };
		std::cout << "m111 == \n" << m11 << std::endl << std::endl;
		std::cout << "typeid(m111).name()  == " << typeid(m11).name() << std::endl;

		// 2. map������û�з������ݿ���������ʹ��map������������ֵ���ᷢ�����ݿ�����
		debugDisp("reinterpret_cast<unsigned>(m1.data()) == ", reinterpret_cast<unsigned>(m1.data()));
		debugDisp("reinterpret_cast<unsigned>(map1.data()) == ", reinterpret_cast<unsigned>(map1.data()));
		debugDisp("reinterpret_cast<unsigned>(m11.data()) == ", reinterpret_cast<unsigned>(m11.data()));

		// 3. mapʵ�־����������ת����
		Eigen::VectorXf v2 = Eigen::Map<Eigen::VectorXf>(m1.data(), 16, 1);
		debugDisp("v2 == \n", v2);
		Eigen::MatrixXf m2 = Eigen::Map<Eigen::MatrixXf>(v2.data(), 4, 4);
		debugDisp("m2 == \n", m2);

		// 4. �޷���const��������map����
		const Eigen::MatrixXf m3 = m2;
		// Eigen::Map<Eigen::MatrixXf> map3 = Eigen::Map<Eigen::MatrixXf>(m3.data(), 16, 1);		// �޷����룻

		debugDisp("finished.");
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
		Eigen::Matrix3f rotation;
		getRotationMat(rotation, RowVector3f(1, 0, 0), RowVector3f(0, 0, 1));
		RowVector3f v1 = (rotation * RowVector3f(1, 0, 0).transpose()).transpose();
		std::cout << "��ת���v1 == " << v1 << std::endl;
		getRotationMat(rotation, RowVector3f(1, 1, 1), RowVector3f(3, 4, 5));
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


	// test14����dense mat�����νṹ��
	/*
		Eigen::Matrix<>�ļ̳�����
				Matrix<> �� PlainObjectBase<> �� MatrixBase<> �� DenseBase<> �� DenseCoeffBase<>���� �� EigenBase<>;

		Eigen::Block<>�ļ̳�����
				Block<> �� BlockImpl<> �� internal::BlockImpl_dense<> �� MapBase<>���� �� MatrixBase<> ��...

		Eigen::CwiseNullaryOp<>�ļ̳���������̳У�
				��̳У�internal::no_assignment_operator + MatrixBase
				MatrixBase �� ...

		Eigen::CwiseBinaryOp<>�ļ̳���
				 ��̳У�internal::no_assignment_operator + CwiseBinaryOpImpl<>
				 CwiseBinaryOpImpl<> �� MatrixBase<> �� ...
	*/
	

	// �����������ϣ����������return type�Ķ���Ӧʹ��Eigen::MatrixBase<>���ͣ�
	template<typename Derived>										// Derived�Ǿ���ľ������ͣ���Eigen::Matrix<int,1,-1,1,1,-1>
	void foo(Eigen::MatrixBase<Derived>& base)				// Eigen::MatrixBase<>��������ܾ��������Ļ��ࣻ
	{
		const Derived& originMat = base.derived();						// derived()����ԭ��������ã�����Ϊԭ���ͣ�
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
		
		// ���ָ�����return type:
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
