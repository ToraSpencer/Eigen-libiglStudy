#include "myEigenBasicMath.h"

// 模板函数需要特化之后才能在静态库中输出： 

/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化：
template std::pair<double, double> cart2polar(const float x, const float y);
template std::pair<double, double> cart2polar(const double x, const double y);

template std::pair<double, double> polar2cart(const float radius, const float theta);
template std::pair<double, double> polar2cart(const double radius, const double theta);

template void eigenVec2Vec(std::vector<int>&, const Eigen::PlainObjectBase<Eigen::RowVectorXi>&);
template void eigenVec2Vec(std::vector<int>&, const Eigen::PlainObjectBase<Eigen::RowVector2i>&);
template void eigenVec2Vec(std::vector<int>&, const Eigen::PlainObjectBase<Eigen::RowVector3i>&);
template void eigenVec2Vec(std::vector<int>&, const Eigen::PlainObjectBase<Eigen::RowVector4i>&);
template void eigenVec2Vec(std::vector<int>&, const Eigen::PlainObjectBase<Eigen::VectorXi>&);
template void eigenVec2Vec(std::vector<int>&, const Eigen::PlainObjectBase<Eigen::Vector2i>&);
template void eigenVec2Vec(std::vector<int>&, const Eigen::PlainObjectBase<Eigen::Vector3i>&);
template void eigenVec2Vec(std::vector<int>&, const Eigen::PlainObjectBase<Eigen::Vector4i>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVectorXf>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVector2f>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVector3f>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVector4f>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::VectorXf>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::Vector2f>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::Vector3f>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::Vector4f>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVectorXi>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVector2i>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVector3i>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVector4i>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::VectorXi>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::Vector2i>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::Vector3i>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::Vector4i>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVectorXd>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVector2d>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVector3d>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::RowVector4d>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::VectorXd>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::Vector2d>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::Vector3d>&);
template void eigenVec2Vec(std::vector<float>&, const Eigen::PlainObjectBase<Eigen::Vector4d>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVectorXd>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVector2d>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVector3d>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVector4d>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::VectorXd>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::Vector2d>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::Vector3d>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::Vector4d>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVectorXi>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVector2i>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVector3i>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVector4i>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::VectorXi>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::Vector2i>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::Vector3i>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::Vector4i>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVectorXf>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVector2f>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVector3f>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::RowVector4f>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::VectorXf>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::Vector2f>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::Vector3f>&);
template void eigenVec2Vec(std::vector<double>&, const Eigen::PlainObjectBase<Eigen::Vector4f>&);

template std::vector<Eigen::RowVectorXi::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVectorXi>& vIn);
template std::vector<Eigen::RowVector2i::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVector2i>& vIn);
template std::vector<Eigen::RowVector3i::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVector3i>& vIn);
template std::vector<Eigen::RowVector4i::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVector4i>& vIn);
template std::vector<Eigen::VectorXi::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::VectorXi>& vIn);
template std::vector<Eigen::Vector2i::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::Vector2i>& vIn);
template std::vector<Eigen::Vector3i::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::Vector3i>& vIn);
template std::vector<Eigen::Vector4i::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::Vector4i>& vIn);
template std::vector<Eigen::RowVectorXf::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVectorXf>& vIn);
template std::vector<Eigen::RowVector2f::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVector2f>& vIn);
template std::vector<Eigen::RowVector3f::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVector3f>& vIn);
template std::vector<Eigen::RowVector4f::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVector4f>& vIn);
template std::vector<Eigen::VectorXf::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::VectorXf>& vIn);
template std::vector<Eigen::Vector2f::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::Vector2f>& vIn);
template std::vector<Eigen::Vector3f::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::Vector3f>& vIn);
template std::vector<Eigen::Vector4f::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::Vector4f>& vIn);
template std::vector<Eigen::RowVectorXd::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVectorXd>& vIn);
template std::vector<Eigen::RowVector2d::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVector2d>& vIn);
template std::vector<Eigen::RowVector3d::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVector3d>& vIn);
template std::vector<Eigen::RowVector4d::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::RowVector4d>& vIn);
template std::vector<Eigen::VectorXd::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::VectorXd>& vIn);
template std::vector<Eigen::Vector2d::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::Vector2d>& vIn);
template std::vector<Eigen::Vector3d::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::Vector3d>& vIn);
template std::vector<Eigen::Vector4d::Scalar>  eigenVec2Vec(const Eigen::PlainObjectBase < Eigen::Vector4d>& vIn);

template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::RowVectorXi>&, const std::vector<int>&);
template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::VectorXi>&, const std::vector<int>&);
template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::RowVectorXf>&, const std::vector<float>&);
template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::VectorXf>&, const std::vector<float>&);
template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::RowVectorXf>&, const std::vector<int>&);
template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::VectorXf>&, const std::vector<int>&);
template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::RowVectorXf>&, const std::vector<double>&);
template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::VectorXf>&, const std::vector<double>&);
template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::RowVectorXd>&, const std::vector<double>&);
template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::VectorXd>&, const std::vector<double>&);
template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::RowVectorXd>&, const std::vector<float>&);
template void vec2EigenVec(Eigen::PlainObjectBase<Eigen::VectorXd>&, const std::vector<float>&);

template Eigen::Matrix<int, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<int>& vIn);
template Eigen::Matrix<float, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<float>& vIn);
template Eigen::Matrix<float, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<double>& vIn); 
 
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

template bool solveLinearEquation(Eigen::Matrix<float, 2, 1>& x, const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<float, 2, 1>& b);
template bool solveLinearEquation(Eigen::Matrix<float, 3, 1>& x, const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<float, 3, 1>& b);
template bool solveLinearEquation(Eigen::Matrix<float, 4, 1>& x, const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<float, 4, 1>& b);
template bool solveLinearEquation(Eigen::Matrix<double, 2, 1>& x, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<double, 2, 1>& b);
template bool solveLinearEquation(Eigen::Matrix<double, 3, 1>& x, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<double, 3, 1>& b);
template bool solveLinearEquation(Eigen::Matrix<double, 4, 1>& x, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<double, 4, 1>& b);

template bool solveLinearEquations(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& X, \
	const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& B);
template bool solveLinearEquations(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& X, \
	const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& B);

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
 
template Eigen::VectorXd fittingStandardEllipse(const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& sampleVers);
template Eigen::VectorXd fittingStandardEllipse(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& sampleVers);