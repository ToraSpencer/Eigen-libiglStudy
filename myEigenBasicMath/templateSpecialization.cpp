#include "myEigenBasicMath.h"

// 模板函数需要特化之后才能在静态库中输出： 

/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化：
template std::pair<double, double> cart2polar(const float x, const float y);
template std::pair<double, double> cart2polar(const double x, const double y);

template std::pair<double, double> polar2cart(const float radius, const float theta);
template std::pair<double, double> polar2cart(const double radius, const double theta);

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

template Eigen::Matrix<int, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<int>& vIn);
template Eigen::Matrix<float, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<float>& vIn);
template Eigen::Matrix<float, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<double>& vIn);

template Eigen::Matrix<float, 3, 3> getRotationMat(const Eigen::Matrix<float, 1, 3>& axisArrow, const float theta);
template Eigen::Matrix<double, 3, 3> getRotationMat(const Eigen::Matrix<double, 1, 3>& axisArrow, const float theta);

template Eigen::Matrix<float, 3, 3> getRotationMat(const Eigen::Matrix<float, 1, 3>& originArrow, \
	const Eigen::Matrix<float, 1, 3>& targetArrow);
template Eigen::Matrix<double, 3, 3> getRotationMat(const Eigen::Matrix<double, 1, 3>& originArrow, \
	const Eigen::Matrix<double, 1, 3>& targetArrow);

template bool subFromIdxVec(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const Eigen::VectorXi& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const Eigen::VectorXi& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const Eigen::VectorXi& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const Eigen::VectorXi& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const Eigen::VectorXi& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const Eigen::VectorXi& vec);

template bool subFromIdxVec(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::vector<int>& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::vector<unsigned>& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::vector<int>& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::vector<unsigned>& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::vector<int>& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::vector<unsigned>& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::vector<int>& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::vector<unsigned>& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::vector<int>& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::vector<unsigned>& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::vector<int>& vec);
template bool subFromIdxVec(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::vector<unsigned>& vec);

template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::vector<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::vector<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::list<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::list<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::deque<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::deque<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::set<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::unordered_set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const std::unordered_set<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::vector<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::vector<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::list<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::list<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::deque<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::deque<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::set<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::unordered_set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const std::unordered_set<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::vector<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::vector<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::list<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::list<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::deque<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::deque<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::set<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::unordered_set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const std::unordered_set<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::vector<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::vector<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::list<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::list<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::deque<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::deque<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::set<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::unordered_set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const std::unordered_set<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::vector<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::vector<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::list<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::list<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::deque<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::deque<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::set<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::unordered_set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const std::unordered_set<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::vector<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::vector<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::list<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::list<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::deque<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::deque<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::set<unsigned int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::unordered_set<int>& con);
template bool subFromIdxCon(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const std::unordered_set<unsigned int>& con);

template bool subFromFlagVec(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, \
	const Eigen::VectorXi& vec);
template bool subFromFlagVec(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, \
	const Eigen::VectorXi& vec);
template bool subFromFlagVec(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, \
	const Eigen::VectorXi& vec);
template bool subFromFlagVec(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, \
	const Eigen::VectorXi& vec);
template bool subFromFlagVec(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, \
	const Eigen::VectorXi& vec);
template bool subFromFlagVec(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, \
	const Eigen::VectorXi& vec);

template bool subFromFlagVec(Eigen::MatrixBase<Eigen::MatrixXi>& matBaseOut, std::vector<int>& oldNewIdxInfo, \
	std::vector<int>& newOldIdxInfo, const Eigen::MatrixBase<Eigen::MatrixXi>& matBaseIn, const Eigen::VectorXi& vec);
template bool subFromFlagVec(Eigen::MatrixBase<Eigen::MatrixXf>& matBaseOut, std::vector<int>& oldNewIdxInfo, \
	std::vector<int>& newOldIdxInfo, const Eigen::MatrixBase<Eigen::MatrixXf>& matBaseIn, const Eigen::VectorXi& vec);
template bool subFromFlagVec(Eigen::MatrixBase<Eigen::MatrixXd>& matBaseOut, std::vector<int>& oldNewIdxInfo, \
	std::vector<int>& newOldIdxInfo, const Eigen::MatrixBase<Eigen::MatrixXd>& matBaseIn, const Eigen::VectorXi& vec);
template bool subFromFlagVec(Eigen::MatrixBase<Eigen::VectorXi>& matBaseOut, std::vector<int>& oldNewIdxInfo, \
	std::vector<int>& newOldIdxInfo, const Eigen::MatrixBase<Eigen::VectorXi>& matBaseIn, const Eigen::VectorXi& vec);
template bool subFromFlagVec(Eigen::MatrixBase<Eigen::VectorXf>& matBaseOut, std::vector<int>& oldNewIdxInfo, \
	std::vector<int>& newOldIdxInfo, const Eigen::MatrixBase<Eigen::VectorXf>& matBaseIn, const Eigen::VectorXi& vec);
template bool subFromFlagVec(Eigen::MatrixBase<Eigen::VectorXd>& matBaseOut, std::vector<int>& oldNewIdxInfo, \
	std::vector<int>& newOldIdxInfo, const Eigen::MatrixBase<Eigen::VectorXd>& matBaseIn, const Eigen::VectorXi& vec);
 
template bool vecInsertNum(Eigen::Matrix<int, Eigen::Dynamic, 1>& vec, const int num);
template bool vecInsertNum(Eigen::Matrix<float, Eigen::Dynamic, 1>& vec, const float num);
template bool vecInsertNum(Eigen::Matrix<double, Eigen::Dynamic, 1>& vec, const double num);

template bool vecInsertVec(Eigen::Matrix<int, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<int, Eigen::Dynamic, 1>& vec2);
template bool vecInsertVec(Eigen::Matrix<float, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<float, Eigen::Dynamic, 1>& vec2);
template bool vecInsertVec(Eigen::Matrix<double, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<double, Eigen::Dynamic, 1>& vec2);
  
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, const Eigen::PlainObjectBase<Eigen::MatrixXi>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, const Eigen::PlainObjectBase<Eigen::Matrix2i>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, const Eigen::PlainObjectBase<Eigen::Matrix3i>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, const Eigen::PlainObjectBase<Eigen::Matrix4i>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, const Eigen::PlainObjectBase<Eigen::RowVectorXi>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, const Eigen::PlainObjectBase<Eigen::RowVector2i>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, const Eigen::PlainObjectBase<Eigen::RowVector3i>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXi>& mat, const Eigen::PlainObjectBase<Eigen::RowVector4i>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, const Eigen::PlainObjectBase<Eigen::MatrixXf>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, const Eigen::PlainObjectBase<Eigen::Matrix2f>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, const Eigen::PlainObjectBase<Eigen::Matrix3f>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, const Eigen::PlainObjectBase<Eigen::Matrix4f>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, const Eigen::PlainObjectBase<Eigen::RowVectorXf>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, const Eigen::PlainObjectBase<Eigen::RowVector2f>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, const Eigen::PlainObjectBase<Eigen::RowVector3f>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXf>& mat, const Eigen::PlainObjectBase<Eigen::RowVector4f>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, const Eigen::PlainObjectBase<Eigen::MatrixXd>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, const Eigen::PlainObjectBase<Eigen::Matrix2d>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, const Eigen::PlainObjectBase<Eigen::Matrix3d>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, const Eigen::PlainObjectBase<Eigen::Matrix4d>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, const Eigen::PlainObjectBase<Eigen::RowVectorXd>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, const Eigen::PlainObjectBase<Eigen::RowVector2d>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, const Eigen::PlainObjectBase<Eigen::RowVector3d>& mat1);
template bool matInsertRows(Eigen::PlainObjectBase<Eigen::MatrixXd>& mat, const Eigen::PlainObjectBase<Eigen::RowVector4d>& mat1);

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
