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



///////////////////////////////////////////////////////////////////////////////////////////////////////// 齐次坐标系相关接口
void objReadVerticesHomoMat(Eigen::MatrixXf& vers, const char* fileName);

void objWriteVerticesHomoMat(const char* fileName, const Eigen::MatrixXf& vers);

void vers2homoVers(Eigen::MatrixXf& homoVers, const Eigen::MatrixXf& vers);

Eigen::MatrixXf vers2homoVers(const Eigen::MatrixXf& vers);

void homoVers2vers(Eigen::MatrixXf& vers, const Eigen::MatrixXf& homoVers);

Eigen::MatrixXf homoVers2vers(const Eigen::MatrixXf& homoVers);

