

//////////////////////////////////////////////////////////////////////////////////////////// 科学计算相关





// 通过SVD求稠密矩阵的条件数：
template<typename T>
double calcCondNum(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m)
{
	Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svdSolver(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::Matrix<T, Eigen::Dynamic, 1> sValues = svdSolver.singularValues();
	T sMax = sValues.maxCoeff();
	T sMin = sValues.minCoeff();
	return static_cast<double>(sMax / sMin);
}


// 霍纳方法（秦九昭算法）求多项式的值
template<typename T, int N>
double hornersPoly(const Eigen::Matrix<T, N, 1>& coeffs, const double x)
{
	// coeffs是{a0, a1, a2, ..., an}组成的(n+1)列向量，多项式为p == a0 + a1*x + a2*x^2 + ... + an* x^n; 
	int n = coeffs.rows() - 1;
	if (n < 0)
		return NAN;

	double result = coeffs(n);
	for (int i = n; i - 1 >= 0; i--)
	{
		result *= x;
		result += coeffs(i - 1);
	}

	return result;
}


// 求多项式的一阶微分
template<typename T, int N>
float polyDiff(const Eigen::Matrix<T, N, 1>& coeffs, const float x)
{
	// 多项式p == a0 + a1*x + a2*x^2 + ... + an* x^n 一阶微分为：p' == a1 + 2*a2*x + 3*a3*x^2 ... n*an*x^(n-1)
	int coeffsCount = coeffs.rows() * coeffs.cols();
	if (coeffsCount <= 0)
		return NAN;

	if (coeffsCount == 1)
		return 1.0f;

	Eigen::VectorXf diffCoeffs(coeffsCount - 1);
	for (int i = 0; i < diffCoeffs.rows(); ++i)
		diffCoeffs(i) = (i + 1) * coeffs(i + 1);

	float result = hornersPoly(diffCoeffs, x);

	return result;
}


// 计算稠密矩阵的克罗内克积（Kronecker product）
template <typename T, typename Derived1, typename Derived2>
void kron(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& result, \
	const Eigen::MatrixBase<Derived1>& mm1, \
	const Eigen::MatrixBase<Derived2>& mm2)
{
	// A(m×n) tensor B(p×q) == C(m*p × n*q);
	/*
		其中C.block(i * m, j* n, p, q) == Aij * B;
		如：
		1 2 3
		4 5 6

		tensor

		100 10
		1		0.1

		==
		100,	10,	200,		20,	 300,		 30,
		1,		0.1,		2,		0.2,		 3,		0.3,
		400,	 40,	500,		 50,	 600, 	60,
		4,		 0.4,		5,		0.5,		 6,		0.6,

	*/
	const Derived1& m1 = mm1.derived();
	const Derived2& m2 = mm2.derived();

	unsigned m = m1.rows();
	unsigned n = m1.cols();
	unsigned p = m2.rows();
	unsigned q = m2.cols();

	result.resize(0, 0);						// 重新分配内存；
	result.resize(m * p, n * q);
	Eigen::VectorXi rowIdxVec = Eigen::VectorXi::LinSpaced(0, m - 1, m - 1);
	auto calcKronCol = [&](const unsigned colIdx, const std::tuple<unsigned>& pt)
	{
		unsigned rowIdx0 = std::get<0>(pt);
		result.block(rowIdx0 * p, colIdx * q, p, q) = static_cast<T>(m1(rowIdx0, colIdx)) * m2.array().cast<T>();		// 转换成输出矩阵中的数据类型；
	};

	auto calcKronRow = [&](const unsigned rowIdx)
	{
		std::tuple<unsigned> pt0 = std::make_tuple(rowIdx);
		PARALLEL_FOR(0, n, calcKronCol, pt0);
	};

	PARALLEL_FOR(0, m, calcKronRow);
}


template <typename Derived1, typename Derived2>
Eigen::MatrixXd kron(const Eigen::MatrixBase<Derived1>& mm1, \
	const Eigen::MatrixBase<Derived2>& mm2)
{
	Eigen::MatrixXd result;
	kron(result, mm1, mm2);
	return result;
}



// 岭回归多项式拟合曲线
template <typename T>
void ridgeRegressionPolyFitting(Eigen::VectorXd& theta, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	unsigned m = 4)
{
	/*
		void ridgeRegressionPolyFitting(
					Eigen::VectorXd & theta,							拟合的多项式函数
					const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers			离散样本点
					const unsigned m						拟合的多项式项数；
					)
	*/

	const double lambda = 0.1;

	int versCount = vers.rows();
	if (versCount == 0)
		return;
	if (m >= versCount)
		m = versCount - 1;

	Eigen::MatrixXd X(versCount, m);
	for (int i = 0; i < versCount; ++i)
		for (int j = 0; j < m; ++j)
			X(i, j) = std::powf(static_cast<double>(vers(i, 0)), j);

	Eigen::VectorXd Y = vers.col(1).cast<double>();
	Eigen::MatrixXd Id(m, m);
	Id.setIdentity();
	theta = (X.transpose() * X + Id * lambda).inverse() * X.transpose() * Y;
}
 





