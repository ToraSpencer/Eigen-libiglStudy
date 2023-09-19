
///////////////////////////////////////////////////////////////////////////////////////////////////// disp interface

template <typename Derived>
void dispMat(const Eigen::PlainObjectBase<Derived>& mat)
{
	std::cout << ": rows == " << mat.rows() << ",  cols == " << mat.cols() << std::endl;
	for (int i = 0; i < mat.rows(); ++i)
	{
		std::cout << i << " ---\t";
		for (int j = 0; j < mat.cols(); ++j)
			std::cout << mat.coeff(i, j) << ", ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


template <typename Derived>
void dispMatBlock(const Eigen::PlainObjectBase<Derived>& mat, const int rowStart, \
	const int rowEnd, const int colStart, const int colEnd)
{
	if (rowEnd > mat.rows() - 1 || rowStart < 0 || rowStart >= rowEnd)
		return;
	if (colEnd > mat.cols() - 1 || colStart < 0 || colStart >= colEnd)
		return;

	std::cout << ": rows == " << mat.rows() << ",  cols == " << mat.cols() << std::endl;
	for (int i = rowStart; i <= rowEnd; ++i)
	{
		std::cout << i << " ---\t";
		for (int j = colStart; j <= colEnd; ++j)
			std::cout << mat(i, j) << ", ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


// 打印稀疏矩阵中的所有非零元素：
template<typename spMat>				// 如果模板参数为T, 函数参数为Eigen::SparseMatrix<T>，编译会报错，不知道什么缘故；
void dispSpMat(const spMat& sm, const unsigned startCol, const unsigned endCol, const unsigned showElemsCount)
{
	unsigned count = 0;
	std::cout << "rows == " << sm.rows() << ", cols == " << sm.cols() << std::endl;

	if (showElemsCount <= 0)
	{
		for (unsigned i = startCol; i <= endCol; ++i)
			for (auto iter = spMat::InnerIterator(sm, i); iter; ++iter)
				std::cout << "(" << iter.row() << ", " << iter.col() << ") " << iter.value() << std::endl;
	}
	else
	{
		for (unsigned i = startCol; i <= endCol; ++i)
			for (auto iter = spMat::InnerIterator(sm, i); iter; ++iter)
			{
				std::cout << "(" << iter.row() << ", " << iter.col() << ") " << iter.value() << std::endl;
				count++;
				if (count >= showElemsCount)
					break;
			}
	}

	std::cout << std::endl;
}


template<typename spMat>
void dispSpMat(const spMat& sm)
{
	dispSpMat(sm, 0, sm.rows() - 1);
}


template<typename T, int N>
void dispVec(const Eigen::Matrix<T, N, 1>& vec)
{
	std::cout << ": rows == " << vec.rows() << std::endl;
	for (int i = 0; i < vec.rows(); ++i)
		std::cout << i << "---\t" << vec(i) << std::endl;
	std::cout << std::endl;
}


template<typename T, int N>
void dispVec(const Eigen::Matrix<T, 1, N>& vec)
{
	std::cout << ": cols == " << vec.cols() << std::endl;
	for (int i = 0; i < vec.cols(); ++i)
		std::cout << i << "---\t" << vec(i) << std::endl;
	std::cout << std::endl;
}


template<typename T, int N>
void dispVecSeg(const Eigen::Matrix<T, N, 1>& vec, const int start, const int end)
{
	if (start < 0 || start >= end)
	{
		std::cout << "input error." << std::endl;
		return;
	}

	int endIdx = (vec.rows() - 1 < end) ? (vec.rows() - 1) : end;
	std::cout << ": rows == " << vec.rows() << std::endl;
	for (int i = start; i <= endIdx; ++i)
		std::cout << i << "---\t" << vec(i) << std::endl;
	std::cout << std::endl;
}


template<typename T, int N>
void dispVecSeg(const Eigen::Matrix<T, 1, N>& vec, const int start, const int end)
{
	if (start < 0 || start >= end)
	{
		std::cout << "input error." << std::endl;
		return;
	}
	int endIdx = (vec.cols() - 1 < end) ? (vec.cols() - 1) : end;
	std::cout << ": cols == " << vec.cols() << std::endl;
	for (int i = start; i <= endIdx; ++i)
		std::cout << i << "---\t" << vec(i) << std::endl;
	std::cout << std::endl;
}

