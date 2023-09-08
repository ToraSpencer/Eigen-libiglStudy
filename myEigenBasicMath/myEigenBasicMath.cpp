#include "myEigenBasicMath.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////// supportive

// ���뺯���ӱ���ϡ������еķ���Ԫ�أ������ӽ��ܵĲ�����Eigen::SparseMatrix<T>::InnerIterator&
template<typename spMat, typename F>
void traverseSparseMatrix(spMat& sm, F f)
{
	for (unsigned i = 0; i < sm.outerSize(); ++i)
		for (auto iter = spMat::InnerIterator(sm, i); iter; ++iter)
			f(iter);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////// ��ģ�庯����ʵ��

// ��������ת��Ϊ����������
Eigen::VectorXi flagVec2IdxVec(const Eigen::VectorXi& flag)
{
	Eigen::VectorXi idxVec;

	for (unsigned i = 0; i < flag.rows(); ++i)
		if (0 != flag(i) && 1 != flag(i))
			return idxVec;

	std::vector<int> tempVec;
	tempVec.reserve(flag.rows());
	for (unsigned i = 0; i < flag.rows(); ++i)
		if (flag(i) > 0)
			tempVec.push_back(i);
	idxVec = vec2EigenVec(tempVec);

	return idxVec;
}


// ��������ת��Ϊ����������
Eigen::VectorXi IdxVec2FlagVec(const Eigen::VectorXi& idxVec, const unsigned size)
{
	Eigen::VectorXi flag;

	if (idxVec.rows() > size)
		return flag;

	for (unsigned i = 0; i < idxVec.rows(); ++i)
		if (idxVec(i) >= size)
			return flag;

	flag = Eigen::VectorXi::Zero(size);
	for (unsigned i = 0; i < idxVec.rows(); ++i)
		flag(idxVec(i)) = 1;

	return flag;
}


// eigen��������std::vector<T>�໥ת��������1: Eigen����ת��Ϊstd::vector
template<typename Derived>
std::vector<typename Derived::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase<Derived>& vIn)
{
	using Scalar = typename Derived::Scalar;
	assert(1 == vIn.rows() || 1 == vIn.cols(), "Error!!! Input arg vIn should be a vector.");

	unsigned elemCount = vIn.rows() * vIn.cols();
	std::vector<Scalar> vOut(elemCount, 1);
	std::memcpy(&vOut[0], vIn.data(), sizeof(Scalar) * elemCount);

	return vOut;
}


// eigen��������std::vector<T>�໥ת��������2��std::vectorת��ΪEigen��������
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<T>& vIn)
{
	unsigned elemCount = vIn.size();
	Eigen::Matrix<T, Eigen::Dynamic, 1> vOut;
	vOut.resize(elemCount, 1);
	std::memcpy(vOut.data(), &vIn[0], sizeof(T) * elemCount);

	return vOut;
}


// ����ʽ��ֵ
void polyInterpolation() {}


// ��˹��ֵ
void gaussInterpolation() {}


// ��С���˶���ʽ������ߣ�
void leastSquarePolyFitting() {}


/////////////////////////////////////////////////////////////////////////////////////////////////////////// ģ�庯����ʵ��

//		��ά�ѿ�������ϵת��Ϊ������ϵ
template <typename T>
std::pair<double, double> cart2polar(const T x, const T y)
{
	// theta���귶ΧΪ[0, 2*pi]��
	double x0 = static_cast<double>(x);
	double y0 = static_cast<double>(y);
	const double pi = 3.14159;
	double eps = 10e-9;

	double radius = std::sqrt(x0 * x0 + y0 * y0);

	if (std::abs(x0) < eps)
	{
		if (std::abs(y0) < eps)
			return { radius, NAN };
		else if (y > 0)
			return { radius, pi / 2 };
		else
			return { radius, 3 * pi / 2 };
	}

	if (std::abs(y0) < eps)
	{
		if (x > 0)
			return { radius, 0 };
		else
			return { radius, -pi };
	}

	// ��radius�ӽ�0�� �����Ӧ��theta��0��
	if (radius < eps)
		return { 0, 0 };

	double theta = std::acos(x / radius);
	if (y < 0)
		theta = 2 * pi - theta;

	return { radius, theta };
}


template <typename T>
std::pair<double, double> polar2cart(const T radius, const T theta)
{
	return { radius * cos(theta), radius * sin(theta) };
}


// ע���������������ת����ȫ��Ĭ��Ϊ����������������v,uΪ��������r,sΪ��������R * v == u; �ȼ��� r * R.transpose() == s;

// getRotationMat() ����1����������ת����������ת�Ƕȣ�������ת����
template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& axisArrow, const float theta)
{
	Eigen::Matrix<T, 3, 1> axis = axisArrow.transpose().normalized();
	return Eigen::AngleAxis<T>(theta, axis).toRotationMatrix();
}


// getRotationMat() ����2�����õ���originArrow��ת��targetArrow����ת����
template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& originArrow, const Eigen::Matrix<T, 1, 3>& targetArrow)
{
	Eigen::Matrix<T, 3, 3> rotation = Eigen::Matrix<T, 3, 3>::Zero();
	if (0 == originArrow.norm() || 0 == targetArrow.norm())
		return rotation;

	Eigen::Matrix<T, 1, 3> axisArrow = originArrow.cross(targetArrow);		// ��ת�᣻

	if (0 == axisArrow.norm())
		return Eigen::Matrix<T, 3, 3>::Identity();

	axisArrow.normalize();
	T x0 = axisArrow(0);
	T y0 = axisArrow(1);
	T z0 = axisArrow(2);
	T cosTheta = originArrow.dot(targetArrow) / (originArrow.norm() * targetArrow.norm());
	T sinTheta = std::sqrt(1 - cosTheta * cosTheta);

	// �ȼ���Eigen::AngleAxis<T>(theta, axis).toRotationMatrix()��������������������תtheta�Ƕȣ�
	rotation << cosTheta + (1 - cosTheta) * x0 * x0, (1 - cosTheta)* x0* y0 - sinTheta * z0, (1 - cosTheta)* x0* z0 + sinTheta * y0, \
		(1 - cosTheta)* y0* x0 + sinTheta * z0, cosTheta + (1 - cosTheta) * y0 * y0, (1 - cosTheta)* y0* z0 - sinTheta * x0, \
		(1 - cosTheta)* z0* x0 - sinTheta * y0, (1 - cosTheta)* z0* y0 + sinTheta * x0, cosTheta + (1 - cosTheta) * z0 * z0;
	return rotation;
}


// ��������������Դ��������ȡԪ�������������
template <typename Derived>
bool subFromIdxVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec)
{
	Derived& matOut = matBaseOut.derived();
	matOut.resize(vec.size(), matBaseIn.cols());
	for (unsigned i = 0; i < vec.rows(); ++i)
	{
		const int& index = vec(i);
		matOut.row(i) = matBaseIn.row(index);
	}

	return true;
}


template <typename Derived, typename Index>
bool subFromIdxVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const std::vector<Index>& vec)
{
	Derived& matOut = matBaseOut.derived();
	matOut.resize(vec.size(), matBaseIn.cols());
	for (unsigned i = 0; i < vec.size(); ++i)
	{
		const Eigen::Index& index = vec[i];
		matOut.row(i) = matBaseIn.row(index);
	}

	return true;
}


// ���������е�������Դ��������ȡԪ�������������
template <typename Derived, typename IndexContainer>
bool subFromIdxCon(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const IndexContainer& con)
{
	Derived& matOut = matBaseOut.derived();
	matOut.resize(con.size(), matBaseIn.cols());

	auto iter = con.begin();
	for (unsigned i = 0; iter != con.end(); ++i)
	{
		const Eigen::Index& index = *iter++;
		matOut.row(i) = matBaseIn.row(index);
	}

	return true;
}


// ����flag������Դ��������ȡԪ�������������
template <typename Derived>
bool subFromFlagVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec)
{
	Derived& matOut = matBaseOut.derived();
	matOut.resize(vec.sum(), matBaseIn.cols());

	int count = 0;
	for (unsigned i = 0; i < vec.rows(); ++i)
		if (vec(i) > 0)
			matOut.row(count++) = matBaseIn.row(i);

	return true;
}


template <typename Derived>
bool subFromFlagVec(Eigen::MatrixBase<Derived>& matBaseOut, std::vector<int>& oldNewIdxInfo, std::vector<int>& newOldIdxInfo, \
	const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& flag)
{
	if (flag.rows() != matBaseIn.rows())
		return false;

	// 1. ��ȡflag���Ϊ1�Ĳ�����
	const int N = flag.rows();
	const int M = flag.sum();
	int count = 0;
	Derived& matOut = matBaseOut.derived();
	matOut.resize(0, 0);
	oldNewIdxInfo.clear();
	newOldIdxInfo.clear();
	matOut.resize(M, matBaseIn.cols());
	for (unsigned i = 0; i < flag.rows(); ++i)
		if (flag(i) > 0)
			matOut.row(count++) = matBaseIn.row(i);

	// 2. ������������ӳ���
	oldNewIdxInfo.resize(N, -1);
	newOldIdxInfo.reserve(M);
	int index = 0;
	for (int k = 0; k < N; ++k)
	{
		int oldIdx = k;
		if (flag(oldIdx) > 0)
		{
			int newIdx = index++;
			oldNewIdxInfo[oldIdx] = newIdx;
			newOldIdxInfo.push_back(oldIdx);
		}
	}

	return true;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////// �������ɾ���

// ������������
template<typename T>
bool vecInsertNum(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const T num)
{
	vec.conservativeResize(vec.rows() + 1, 1);
	vec(vec.rows() - 1) = num;
	return true;
}


// ������������
template<typename T>
bool vecInsertVec(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec2)
{
	unsigned currentRows = vec1.rows();
	unsigned addRows = vec2.rows();
	unsigned finalRows = currentRows + addRows;
	vec1.conservativeResize(finalRows);

	// �������ݣ�
	T* dataPtr = vec1.data();
	dataPtr += currentRows;
	std::memcpy(dataPtr, vec2.data(), sizeof(T) * addRows);

	return true;
}


// ����ĩβ�������/��������
template <typename Derived1, typename Derived2>
bool matInsertRows(Eigen::PlainObjectBase<Derived1>& mat, const Eigen::PlainObjectBase<Derived2>& mat1)
{
	assert(mat.cols() == mat1.cols(), "Error!!! Matrix size not match.");
	unsigned cols = mat1.cols();
	unsigned currentRows = mat.rows();
	unsigned addRows = mat1.rows();
	mat.conservativeResize(currentRows + addRows, cols);
	for (unsigned i = 0; i < addRows; ++i)
		mat.row(currentRows + i) = mat1.row(i);

	return true;
}


// ����һ��flag������retVec����mat�ĵ�i�к�������vec��ȣ���retVec(i)==1���������0��
template <typename Derived1, typename Derived2>
Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Derived1>& mat, const Eigen::PlainObjectBase<Derived2>& rowVec)
{
	int rows = mat.rows();
	int cols = mat.cols();
	Eigen::VectorXi retVec(rows);
	assert(rowVec.cols() == cols, "Error!!! Tow mats size do not match.");

	// ���бȽϣ�
	Eigen::MatrixXi tempMat(rows, cols);
	for (int i = 0; i < cols; ++i)
		tempMat.col(i) = (mat.col(i).array() == rowVec(i)).select(Eigen::VectorXi::Ones(rows), Eigen::VectorXi::Zero(rows));

	retVec = tempMat.col(0);

	if (cols > 1)
		for (int i = 1; i < cols; ++i)
			retVec = retVec.array() * tempMat.col(i).array();			// ������ˣ�

	return retVec;
}


// ����ϡ����󣬺�matlab�е�sparse()���ƣ�
template <class ValueVector, typename T>
void sparse(Eigen::SparseMatrix<T>& SM, const Eigen::VectorXi& I, const Eigen::VectorXi& J,
	const ValueVector& values, const size_t m, const size_t n)
{
	assert((int)I.maxCoeff() < (int)m);
	assert((int)I.minCoeff() >= 0);
	assert((int)J.maxCoeff() < (int)n);
	assert((int)J.minCoeff() >= 0);
	assert(I.size() == J.size());
	assert(J.size() == values.size());

	// Really we just need .size() to be the same, but this is safer
	assert(I.rows() == J.rows());
	assert(J.rows() == values.rows());
	assert(J.cols() == values.cols());

	std::vector<Eigen::Triplet<T>> IJV;
	IJV.reserve(I.size());
	for (int x = 0; x < I.size(); x++)
		IJV.push_back(Eigen::Triplet<T>(I(x), J(x), values(x)));
	SM.resize(m, n);
	SM.setFromTriplets(IJV.begin(), IJV.end());
}


// ϡ�����ת�ã�Eigen::SparseMatrix�Դ���transpose()����̫������
template<typename T>
bool spMatTranspose(Eigen::SparseMatrix<T>& smOut, const Eigen::SparseMatrix<T>& smIn)
{
	smOut.resize(0, 0);
	smOut.resize(smIn.cols(), smIn.rows());
	std::vector<Eigen::Triplet<T>> trips;
	trips.reserve(smIn.nonZeros());
	traverseSparseMatrix(smIn, [&smIn, &trips](auto& iter)
		{
			trips.push_back(Eigen::Triplet<T>{static_cast<int>(iter.col()), static_cast<int>(iter.row()), iter.value()});
		});
	smOut.setFromTriplets(trips.begin(), trips.end());

	return true;
}


// ��ǡ���ĳ������Է�����Ax == b;
template<typename T, int N>
bool solveLinearEquation(Eigen::Matrix<T, N, 1>& x, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<T, N, 1>& b)
{
	// �����Է�����Ax == b;
	Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svdSolver(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	x = svdSolver.solve(b);

	return true;
}


// ��һϵ��ǡ���ĳ������Է�����AX == B;
template <typename T>
bool solveLinearEquations(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& X, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& B)
{
	if (A.rows() != B.rows())
		return false;

	Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svdSolver(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	X.resize(A.cols(), B.cols());
	for (int i = 0; i < B.cols(); ++i)
	{
		Eigen::Matrix < T, Eigen::Dynamic, 1> x = svdSolver.solve(B.col(i));
		X.col(i) = x;
	}

	return true;
}


#include "templateSpecialization.cpp"
