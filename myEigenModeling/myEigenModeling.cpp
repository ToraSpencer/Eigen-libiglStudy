#include "myEigenModeling.h"
 
///////////////////////////////////////////////////////////////////////////////////////////////////// modeling接口：

// 生成中心在原点，边长为1，三角片数为12的正方体网格； 
template	<typename T>
void genCubeMesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris)
{
	vers.resize(0, 0);
	tris.resize(0, 0);
	vers.resize(8, 3);
	vers << -0.5000000, -0.5000000, -0.5000000, -0.5000000, 0.5000000, -0.5000000, \
		0.5000000, -0.5000000, -0.5000000, 0.5000000, 0.5000000, -0.5000000, \
		0.5000000, -0.5000000, 0.5000000, 0.5000000, 0.5000000, 0.5000000, \
		- 0.5000000, -0.5000000, 0.5000000, -0.5000000, 0.5000000, 0.5000000;

	tris.resize(12, 3);
	tris << 1, 2, 0, 1, 3, 2, 3, 4, 2, 3, 5, 4, 0, 4, 6, 0, 2, 4, 7, 3, 1, 7, 5, 3, 7, 0, 6, 7, 1, 0, 5, 6, 4, 5, 7, 6;
}


// 生成轴向包围盒的三角网格；
template <typename T>
void genAABBmesh(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, \
	const Eigen::AlignedBox<T, 3>& aabb)
{
	vers.resize(0, 0);
	tris.resize(0, 0);
	Eigen::Matrix<T, 3, 1> minp = aabb.min();
	Eigen::Matrix<T, 3, 1> maxp = aabb.max();
	Eigen::Matrix<T, 3, 1> newOri = (minp + maxp) / 2.0;
	Eigen::Matrix<T, 3, 1> sizeVec = maxp - minp;

	genCubeMesh(vers, tris);
	vers.col(0) *= sizeVec(0);
	vers.col(1) *= sizeVec(1);
	vers.col(2) *= sizeVec(2);
	vers.rowwise() += newOri.transpose();
}
 

bool SDFvec2mat(std::vector<Eigen::MatrixXf>& matLayers, \
	const std::vector<float>& SDFvec, const std::vector<unsigned>& stepsCount)
{
	const int xCount = static_cast<int>(stepsCount[0]);
	const int yCount = static_cast<int>(stepsCount[1]);
	const int zCount = static_cast<int>(stepsCount[2]);
	matLayers.clear();
	matLayers.resize(zCount, Eigen::MatrixXf(yCount, xCount));		// ！！！注意：x是列，y是行   

	Eigen::VectorXf SDFvecCopy;
	vec2EigenVec(SDFvecCopy, SDFvec);
	Eigen::MatrixXf tmpMat{ Eigen::Map<Eigen::MatrixXf>(SDFvecCopy.data(), xCount * yCount, zCount)};
	for (int i = 0; i < zCount; ++i)
	{
		Eigen::VectorXf tmpVec = tmpMat.col(i);
		matLayers[i] = Eigen::Map<Eigen::MatrixXf>(tmpVec.data(), yCount, xCount);
	}
 
	return true;
}


bool SDFmat2vec(std::vector<float>& SDFvec, \
	const std::vector<Eigen::MatrixXf>& matLayers, const std::vector<unsigned>& stepsCount)
{
	const int xCount = static_cast<int>(stepsCount[0]);
	const int yCount = static_cast<int>(stepsCount[1]);
	const int zCount = static_cast<int>(stepsCount[2]);
	const int elemCount = xCount * yCount * zCount;
	const int sliceElemCount = xCount * yCount;
	SDFvec.clear();
	SDFvec.resize(elemCount);
	int pos = 0;
	const float* dataPtr = nullptr;
	for (int i = 0; i < zCount; ++i)
	{
		dataPtr = matLayers[i].data();
		std::memcpy(&SDFvec[pos], dataPtr, sizeof(float) * sliceElemCount);
		pos += sliceElemCount;
	}

	return true;
}


// 模板特化输出：
#include "templateSpecialization.cpp"