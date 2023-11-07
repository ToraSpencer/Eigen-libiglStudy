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
// ��algorithm����һ����ʹ�õ����ȣ�libigl���з�װ�������ʷ�ʹ�õ���˫���ȣ�
#define ANSI_DECLARATORS
#define REAL DOUBLE
#define VOID int
#include "triangulate.h"
#endif

#include "Eigen/Dense"
#include "Eigen/Sparse"
#define WIN32_LEAN_AND_MEAN             // �� Windows ͷ�ļ����ų�����ʹ�õ�����
 

///////////////////////////////////////////////////////////////////////////////////////////////////////// IO 

///////////////////////////////////////////////////////////////////////////////////////////////////////// ����任

// �ѿ���������������ϵ
template <typename DerivedVo, typename DerivedVi>
bool vers2HomoVers(Eigen::PlainObjectBase<DerivedVo>& versOut, const Eigen::PlainObjectBase<DerivedVi>& versIn);

// �������ϵ���ѿ�������ϵ
template <typename DerivedVo, typename DerivedVi>
bool homoVers2Vers(Eigen::PlainObjectBase<DerivedVo>& versOut, const Eigen::PlainObjectBase<DerivedVi>& versIn);

template <typename DerivedV>
Eigen::MatrixXd homoVers2Vers(const Eigen::PlainObjectBase<DerivedV>& versIn);
 
 
