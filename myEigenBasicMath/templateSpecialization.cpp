#include "myEigenBasicMath.h"

// 模板函数需要特化之后才能在静态库中输出： 

/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化：
template std::pair<double, double> cart2polar(const float x, const float y);
template std::pair<double, double> cart2polar(const double x, const double y);

template std::pair<double, double> polar2cart(const float radius, const float theta);
template std::pair<double, double> polar2cart(const double radius, const double theta); 
 
template bool vecInsertNum(Eigen::Matrix<int, Eigen::Dynamic, 1>& vec, const int num);
template bool vecInsertNum(Eigen::Matrix<float, Eigen::Dynamic, 1>& vec, const float num);
template bool vecInsertNum(Eigen::Matrix<double, Eigen::Dynamic, 1>& vec, const double num);

template bool vecInsertVec(Eigen::Matrix<int, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<int, Eigen::Dynamic, 1>& vec2);
template bool vecInsertVec(Eigen::Matrix<float, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<float, Eigen::Dynamic, 1>& vec2);
template bool vecInsertVec(Eigen::Matrix<double, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<double, Eigen::Dynamic, 1>& vec2);


template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXi>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXi>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXi>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXi>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4i>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4i>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXf>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXf>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXf>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXf>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4f>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4f>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXd>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2d>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3d>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4d>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXd>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2d>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3d>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix2d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4d>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXd>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2d>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3d>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix3d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4d>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVectorXd>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector2d>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector3d>& rowVec);
template Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Eigen::Matrix4d>& mat, \
	const Eigen::PlainObjectBase<Eigen::RowVector4d>& rowVec);

template void sparse(Eigen::SparseMatrix<float>& SM, const Eigen::VectorXi& I, const Eigen::VectorXi& J,
	const Eigen::VectorXf& values, const size_t m, const size_t n);
template void sparse(Eigen::SparseMatrix<float>& SM, const Eigen::VectorXi& I, const Eigen::VectorXi& J,
	const Eigen::VectorXd& values, const size_t m, const size_t n);
template void sparse(Eigen::SparseMatrix<double>& SM, const Eigen::VectorXi& I, const Eigen::VectorXi& J,
	const Eigen::VectorXf& values, const size_t m, const size_t n);
template void sparse(Eigen::SparseMatrix<double>& SM, const Eigen::VectorXi& I, const Eigen::VectorXi& J,
	const Eigen::VectorXd& values, const size_t m, const size_t n);

template bool spMatTranspose(Eigen::SparseMatrix<int>& smOut, const Eigen::SparseMatrix<int>& smIn);
template bool spMatTranspose(Eigen::SparseMatrix<float>& smOut, const Eigen::SparseMatrix<float>& smIn);
template bool spMatTranspose(Eigen::SparseMatrix<double>& smOut, const Eigen::SparseMatrix<double>& smIn);

template double calcCondNum(const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& m);
template double calcCondNum(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& m);
 
template double hornersPoly(const Eigen::Matrix<float, 2, 1>& coeffs, const double x);
template double hornersPoly(const Eigen::Matrix<float, 3, 1>& coeffs, const double x);
template double hornersPoly(const Eigen::Matrix<float, 4, 1>& coeffs, const double x); 
template double hornersPoly(const Eigen::Matrix<float, -1, 1>& coeffs, const double x);
template double hornersPoly(const Eigen::Matrix<double, 2, 1>& coeffs, const double x);
template double hornersPoly(const Eigen::Matrix<double, 3, 1>& coeffs, const double x);
template double hornersPoly(const Eigen::Matrix<double, 4, 1>& coeffs, const double x);
template double hornersPoly(const Eigen::Matrix<double, -1, 1>& coeffs, const double x);

template float polyDiff(const Eigen::Matrix<float, 2, 1>& coeffs, const float x);
template float polyDiff(const Eigen::Matrix<float, 3, 1>& coeffs, const float x);
template float polyDiff(const Eigen::Matrix<float, 4, 1>& coeffs, const float x);
template float polyDiff(const Eigen::Matrix<float, -1, 1>& coeffs, const float x);
template float polyDiff(const Eigen::Matrix<double, 2, 1>& coeffs, const float x);
template float polyDiff(const Eigen::Matrix<double, 3, 1>& coeffs, const float x);
template float polyDiff(const Eigen::Matrix<double, 4, 1>& coeffs, const float x);
template float polyDiff(const Eigen::Matrix<double, -1, 1>& coeffs, const float x);

template void ridgeRegressionPolyFitting(Eigen::VectorXd& theta, \
	const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, unsigned m);
template void ridgeRegressionPolyFitting(Eigen::VectorXd& theta, \
	const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, unsigned m);
