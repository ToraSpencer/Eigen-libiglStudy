#include "myEigenHomoCoor.h"


// ģ�庯����Ҫ�ػ�֮������ھ�̬��������� 

/////////////////////////////////////////////////////////////////////////////////////////////////////////// ģ���ػ��� 
template bool vers2HomoVers<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>>(\
	Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>>&, \
	const Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>>&);
template bool vers2HomoVers<Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1>>(\
	Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&, \
	const Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&);
template bool vers2HomoVers<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>>(\
	Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&, \
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&);
template bool vers2HomoVers<Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>>(\
	Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&, \
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&);
template bool vers2HomoVers<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1>>(\
	Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&, \
	const Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&);

template bool homoVers2Vers<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>>(\
	Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>>&, \
	const Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>>&);
template bool homoVers2Vers<Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1>>(\
	Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&, \
	const Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&);
template bool homoVers2Vers<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>>(\
	Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&, \
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&);
template bool homoVers2Vers<Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>>(\
	Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&, \
	const Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&);
template bool homoVers2Vers<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1>>(\
	Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&, \
	const Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&);