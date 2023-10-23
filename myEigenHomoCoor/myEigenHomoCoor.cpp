#include "myEigenHomoCoor.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////// IO 

///////////////////////////////////////////////////////////////////////////////////////////////////////// 表象变换

// 笛卡尔坐标→齐次坐标系
template <typename DerivedVo, typename DerivedVi>
bool vers2HomoVers(Eigen::PlainObjectBase<DerivedVo>& versOut, const Eigen::PlainObjectBase<DerivedVi>& versIn)
{
	using ScalarI = typename DerivedVi::Scalar;
	using ScalarO = typename DerivedVo::Scalar;

	const int dim = versIn.cols();
	const int versCount = versIn.rows();
	assert(2 == dim || 3 == dim, "assert!!! input vertices dimension should be 2 or 3.");

	versOut.resize(0, 0);
	if (0 == versCount)
		return true;

	versOut.resize(dim + 1, versCount);
	versOut.setOnes();
	versOut.topRows(dim) = versIn.transpose().array().cast<ScalarO>();

	return true;
}


// 齐次坐标系→笛卡尔坐标系
template <typename DerivedVo, typename DerivedVi>
bool homoVers2Vers(Eigen::PlainObjectBase<DerivedVo>& versOut, const Eigen::PlainObjectBase<DerivedVi>& versIn)
{
	using ScalarI = typename DerivedVi::Scalar;
	using ScalarO = typename DerivedVo::Scalar;

	const int dim = versIn.rows() - 1;
	const int versCount = versIn.cols();
	assert(2 == dim || 3 == dim, "assert!!! input vertices dimension should be 2 or 3.");

	versOut.resize(0, 0);
	if (0 == versCount)
		return true;

	//versOut.resize(versCount, dim);
	//versOut.setOnes();
	//versOut.topRows(dim) = versIn.transpose().array().cast<ScalarO>();

	versOut = versIn.transpose().leftCols(dim).array().cast<ScalarO>();

	return true; 
}







#include "templateSpecialization.cpp"