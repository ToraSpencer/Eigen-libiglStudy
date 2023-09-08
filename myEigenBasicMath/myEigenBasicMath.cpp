#include "myEigenBasicMath.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////// supportive

// 传入函数子遍历稀疏矩阵中的非零元素，函数子接受的参数是Eigen::SparseMatrix<T>::InnerIterator&
template<typename spMat, typename F>
void traverseSparseMatrix(spMat& sm, F f)
{
	for (unsigned i = 0; i < sm.outerSize(); ++i)
		for (auto iter = spMat::InnerIterator(sm, i); iter; ++iter)
			f(iter);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////// 非模板函数的实现

// 布尔向量转化为索引向量；
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


// 索引向量转化为布尔向量；
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


// eigen的向量和std::vector<T>相互转换，重载1: Eigen向量转换为std::vector
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


// eigen的向量和std::vector<T>相互转换，重载2：std::vector转换为Eigen列向量；
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<T>& vIn)
{
	unsigned elemCount = vIn.size();
	Eigen::Matrix<T, Eigen::Dynamic, 1> vOut;
	vOut.resize(elemCount, 1);
	std::memcpy(vOut.data(), &vIn[0], sizeof(T) * elemCount);

	return vOut;
}


// 多项式插值
void polyInterpolation() {}


// 高斯插值
void gaussInterpolation() {}


// 最小二乘多项式拟合曲线：
void leastSquarePolyFitting() {}


/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板函数的实现

//		二维笛卡尔坐标系转换为极坐标系
template <typename T>
std::pair<double, double> cart2polar(const T x, const T y)
{
	// theta坐标范围为[0, 2*pi]；
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

	// 若radius接近0， 将其对应的theta置0；
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


// 注！！！计算出的旋转矩阵全部默认为作用于列向量；若v,u为列向量，r,s为行向量：R * v == u; 等价于 r * R.transpose() == s;

// getRotationMat() 重载1——输入旋转轴向量，旋转角度，返回旋转矩阵：
template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& axisArrow, const float theta)
{
	Eigen::Matrix<T, 3, 1> axis = axisArrow.transpose().normalized();
	return Eigen::AngleAxis<T>(theta, axis).toRotationMatrix();
}


// getRotationMat() 重载2——得到将originArrow旋转到targetArrow的旋转矩阵
template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& originArrow, const Eigen::Matrix<T, 1, 3>& targetArrow)
{
	Eigen::Matrix<T, 3, 3> rotation = Eigen::Matrix<T, 3, 3>::Zero();
	if (0 == originArrow.norm() || 0 == targetArrow.norm())
		return rotation;

	Eigen::Matrix<T, 1, 3> axisArrow = originArrow.cross(targetArrow);		// 旋转轴；

	if (0 == axisArrow.norm())
		return Eigen::Matrix<T, 3, 3>::Identity();

	axisArrow.normalize();
	T x0 = axisArrow(0);
	T y0 = axisArrow(1);
	T z0 = axisArrow(2);
	T cosTheta = originArrow.dot(targetArrow) / (originArrow.norm() * targetArrow.norm());
	T sinTheta = std::sqrt(1 - cosTheta * cosTheta);

	// 等价于Eigen::AngleAxis<T>(theta, axis).toRotationMatrix()，计算绕任意轴向量旋转theta角度；
	rotation << cosTheta + (1 - cosTheta) * x0 * x0, (1 - cosTheta)* x0* y0 - sinTheta * z0, (1 - cosTheta)* x0* z0 + sinTheta * y0, \
		(1 - cosTheta)* y0* x0 + sinTheta * z0, cosTheta + (1 - cosTheta) * y0 * y0, (1 - cosTheta)* y0* z0 - sinTheta * x0, \
		(1 - cosTheta)* z0* x0 - sinTheta * y0, (1 - cosTheta)* z0* y0 + sinTheta * x0, cosTheta + (1 - cosTheta) * z0 * z0;
	return rotation;
}


// 根据索引向量从源矩阵中提取元素生成输出矩阵。
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


// 根据容器中的索引从源矩阵中提取元素生成输出矩阵。
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


// 根据flag向量从源矩阵中提取元素生成输出矩阵。
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

	// 1. 抽取flag标记为1的层数：
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

	// 2. 生成新老索引映射表：
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



///////////////////////////////////////////////////////////////////////////////////////////////////////// 矩阵的增删查改

// 向量插入数据
template<typename T>
bool vecInsertNum(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const T num)
{
	vec.conservativeResize(vec.rows() + 1, 1);
	vec(vec.rows() - 1) = num;
	return true;
}


// 向量插入向量
template<typename T>
bool vecInsertVec(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec2)
{
	unsigned currentRows = vec1.rows();
	unsigned addRows = vec2.rows();
	unsigned finalRows = currentRows + addRows;
	vec1.conservativeResize(finalRows);

	// 拷贝数据：
	T* dataPtr = vec1.data();
	dataPtr += currentRows;
	std::memcpy(dataPtr, vec2.data(), sizeof(T) * addRows);

	return true;
}


// 矩阵末尾插入矩阵/行向量：
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


// 返回一个flag列向量retVec，若mat的第i行和行向量vec相等，则retVec(i)==1，否则等于0；
template <typename Derived1, typename Derived2>
Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Derived1>& mat, const Eigen::PlainObjectBase<Derived2>& rowVec)
{
	int rows = mat.rows();
	int cols = mat.cols();
	Eigen::VectorXi retVec(rows);
	assert(rowVec.cols() == cols, "Error!!! Tow mats size do not match.");

	// 逐列比较：
	Eigen::MatrixXi tempMat(rows, cols);
	for (int i = 0; i < cols; ++i)
		tempMat.col(i) = (mat.col(i).array() == rowVec(i)).select(Eigen::VectorXi::Ones(rows), Eigen::VectorXi::Zero(rows));

	retVec = tempMat.col(0);

	if (cols > 1)
		for (int i = 1; i < cols; ++i)
			retVec = retVec.array() * tempMat.col(i).array();			// 逐列相乘：

	return retVec;
}


// 生成稀疏矩阵，和matlab中的sparse()类似；
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


// 稀疏矩阵转置；Eigen::SparseMatrix自带的transpose()方法太垃圾了
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


// 解恰定的稠密线性方程组Ax == b;
template<typename T, int N>
bool solveLinearEquation(Eigen::Matrix<T, N, 1>& x, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<T, N, 1>& b)
{
	// 解线性方程组Ax == b;
	Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svdSolver(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	x = svdSolver.solve(b);

	return true;
}


// 解一系列恰定的稠密线性方程组AX == B;
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
