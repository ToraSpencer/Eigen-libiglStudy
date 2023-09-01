#pragma once
#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <fstream>
#include <queue>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
#include <initializer_list>
#include <memory>
#include <thread>
#include <limits>
#include <windows.h> 

#include "Eigen/Dense"
#include "Eigen/Sparse"

/////////////////////////////////////////////////////////////////////////////////////////////////// basic math tools:
template <typename T>
std::pair<double, double> cart2polar(const T x, const T y);
template <typename T>
std::pair<double, double> polar2cart(const T radius, const T theta);

template<typename Derived>
std::vector<typename Derived::Scalar> eigenVec2Vec(const Eigen::PlainObjectBase<Derived>& vIn);

template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> vec2EigenVec(const std::vector<T>& vIn);

template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& axisArrow, const float theta);
template <typename T>
Eigen::Matrix<T, 3, 3> getRotationMat(const Eigen::Matrix<T, 1, 3>& originArrow, const Eigen::Matrix<T, 1, 3>& targetArrow);

template <typename Derived>
bool subFromIdxVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec);
template <typename Derived, typename Index>
bool subFromIdxVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const std::vector<Index>& vec);
template <typename Derived, typename IndexContainer>
bool subFromIdxCon(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const IndexContainer& con);
template <typename Derived>
bool subFromFlagVec(Eigen::MatrixBase<Derived>& matBaseOut, const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec);
template <typename Derived>
bool subFromFlagVec(Eigen::MatrixBase<Derived>& matBaseOut, std::vector<int>& oldNewIdxInfo, std::vector<int>& newOldIdxInfo, \
	const Eigen::MatrixBase<Derived>& matBaseIn, const Eigen::VectorXi& vec);
 

template<typename T>
bool vecInsertNum(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, const T num);
template<typename T>
bool vecInsertVec(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec2);
 
template <typename Derived1, typename Derived2>
bool matInsertRows(Eigen::PlainObjectBase<Derived1>& mat, const Eigen::PlainObjectBase<Derived2>& mat1);
 
template <typename Derived1, typename Derived2>
Eigen::VectorXi rowInMat(const Eigen::PlainObjectBase<Derived1>& mat, const Eigen::PlainObjectBase<Derived2>& rowVec);

template <class ValueVector, typename T>
void sparse(Eigen::SparseMatrix<T>& SM, const Eigen::VectorXi& I, const Eigen::VectorXi& J,
	const ValueVector& values, const size_t m, const size_t n);

template<typename T>
bool spMatTranspose(Eigen::SparseMatrix<T>& smOut, const Eigen::SparseMatrix<T>& smIn);
template<typename T, int N>
bool solveLinearEquation(Eigen::Matrix<T, N, 1>& x, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, const Eigen::Matrix<T, N, 1>& b);
template <typename T>
bool solveLinearEquations(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& X, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, \
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& B);
template<typename T>
double calcCondNum(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m);
template<typename T, int N>
double hornersPoly(const Eigen::Matrix<T, N, 1>& coeffs, const double x);
template<typename T, int N>
float polyDiff(const Eigen::Matrix<T, N, 1>& coeffs, const float x);
template <typename T, typename Derived1, typename Derived2>
void kron(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& result, \
	const Eigen::MatrixBase<Derived1>& mm1, \
	const Eigen::MatrixBase<Derived2>& mm2);
template <typename Derived1, typename Derived2>
Eigen::MatrixXd kron(const Eigen::MatrixBase<Derived1>& mm1, \
	const Eigen::MatrixBase<Derived2>& mm2);
template <typename T>
void ridgeRegressionPolyFitting(Eigen::VectorXd& theta, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, unsigned);


Eigen::VectorXi flagVec2IdxVec(const Eigen::VectorXi& flag);
Eigen::VectorXi IdxVec2FlagVec(const Eigen::VectorXi& idxVec, const unsigned size);
void polyInterpolation();
void gaussInterpolation();
void leastSquarePolyFitting();



#include "myEigenBasicMath.tpp"