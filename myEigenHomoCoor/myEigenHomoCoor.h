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

//#define USE_TRIANGLE_H

#ifdef USE_TRIANGLE_H
// 和algorithm工程一样，使用单精度；libigl库中封装的三角剖分使用的是双精度；
#define ANSI_DECLARATORS
#define REAL DOUBLE
#define VOID int
#include "triangulate.h"
#endif

#include "Eigen/Dense"
#include "Eigen/Sparse"
#define WIN32_LEAN_AND_MEAN             // 从 Windows 头文件中排除极少使用的内容
 

///////////////////////////////////////////////////////////////////////////////////////////////////////// IO 

///////////////////////////////////////////////////////////////////////////////////////////////////////// 表象变换

// 笛卡尔坐标→齐次坐标系
template <typename DerivedVo, typename DerivedVi>
bool vers2HomoVers(Eigen::PlainObjectBase<DerivedVo>& versOut, const Eigen::PlainObjectBase<DerivedVi>& versIn);

// 齐次坐标系→笛卡尔坐标系
template <typename DerivedVo, typename DerivedVi>
bool homoVers2Vers(Eigen::PlainObjectBase<DerivedVo>& versOut, const Eigen::PlainObjectBase<DerivedVi>& versIn);

template <typename DerivedV>
Eigen::MatrixXd homoVers2Vers(const Eigen::PlainObjectBase<DerivedV>& versIn);
 
 
