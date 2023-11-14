#include "myEigenModeling.h"

// 模板函数需要特化之后才能在静态库中输出： 


/////////////////////////////////////////////////////////////////////////////////////////////////////////// 模板特化：

template void genCubeMesh(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);
template void genCubeMesh(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris);
 
template void genAABBmesh(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<float, 3>& aabb);
template void genAABBmesh(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<double, 3>& aabb);


#ifdef USE_TRIANGLE_H
template bool circuit2mesh<Eigen::MatrixXd, Eigen::MatrixXd>(Eigen::PlainObjectBase<Eigen::MatrixXd>&,\
	Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXd>&);
template bool circuit2mesh<Eigen::MatrixXf, Eigen::MatrixXf>(Eigen::PlainObjectBase<Eigen::MatrixXf>&, \
	Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXf>&);
template bool circuit2mesh<Eigen::MatrixXf, Eigen::MatrixXd>(Eigen::PlainObjectBase<Eigen::MatrixXf>&, \
	Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXd>&);
template bool circuit2mesh<Eigen::MatrixXd, Eigen::MatrixXf>(Eigen::PlainObjectBase<Eigen::MatrixXd>&, \
	Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXf>&);

template bool circuit2mesh<Eigen::MatrixXf, Eigen::MatrixXf, Eigen::RowVector3f>(\
	Eigen::PlainObjectBase<Eigen::MatrixXf>&, Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXf>&, \
	const Eigen::PlainObjectBase<Eigen::RowVector3f>&, const float);
template bool circuit2mesh<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::RowVector3d>(\
	Eigen::PlainObjectBase<Eigen::MatrixXd>&, Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXd>&, \
	const Eigen::PlainObjectBase<Eigen::RowVector3d>&, const float);
template bool circuit2mesh<Eigen::MatrixXf, Eigen::MatrixXd, Eigen::RowVector3d>(\
	Eigen::PlainObjectBase<Eigen::MatrixXf>&, Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXd>&, \
	const Eigen::PlainObjectBase<Eigen::RowVector3d>&, const float);
template bool circuit2mesh<Eigen::MatrixXd, Eigen::MatrixXf, Eigen::RowVector3d>(\
	Eigen::PlainObjectBase<Eigen::MatrixXd>&, Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXf>&, \
	const Eigen::PlainObjectBase<Eigen::RowVector3d>&, const float);
template bool circuit2mesh<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::RowVector3f>(\
	Eigen::PlainObjectBase<Eigen::MatrixXd>&, Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXd>&, \
	const Eigen::PlainObjectBase<Eigen::RowVector3f>&, const float);
template bool circuit2mesh<Eigen::MatrixXd, Eigen::MatrixXf, Eigen::RowVector3f>(\
	Eigen::PlainObjectBase<Eigen::MatrixXd>&, Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXf>&, \
	const Eigen::PlainObjectBase<Eigen::RowVector3f>&, const float);
template bool circuit2mesh<Eigen::MatrixXf, Eigen::MatrixXd, Eigen::RowVector3f>(\
	Eigen::PlainObjectBase<Eigen::MatrixXf>&, Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXd>&, \
	const Eigen::PlainObjectBase<Eigen::RowVector3f>&, const float);
template bool circuit2mesh<Eigen::MatrixXf, Eigen::MatrixXf, Eigen::RowVector3d>(\
	Eigen::PlainObjectBase<Eigen::MatrixXf>&, Eigen::MatrixXi&, const Eigen::PlainObjectBase<Eigen::MatrixXf>&, \
	const Eigen::PlainObjectBase<Eigen::RowVector3d>&, const float);

#endif