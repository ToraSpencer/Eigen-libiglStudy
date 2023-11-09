#include "test_libigl.h" 

#include<stdio.h>
#include<assert.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <typeinfo>

#include <windows.h>
#include <atlstr.h>					// ����CString�ࡣ����microsoft ATL(�ģ���avtive template library)
#include <atlconv.h>
#include <io.h>
#include <winuser.h>


#define DATA_PATH "./data/"

static igl::opengl::glfw::Viewer viewer;				// ȫ�ֵ�����鿴������
static std::mutex g_mutex;						// ȫ�ֵĻ�������


// ��ǰ����-easy
/*

*/


// ��ǰ����-hard
/*
	1. �������������������㣬�ҳ�����֮������·����

	2. ʹ��OPENGLʵ����������⣬�����������ȷ�������ϵ�������������Ϊ���롣


*/


// ��Ŀ��Ϣ
/*
	���뻷���� x64 Relase

	ʹ�õĵ�������:
		eigen
		libigl
		glad
		glfw
*/


// ��Ŀ�м���Ԥ����꣺
/*
	IGL_STATIC_LIBRARY

	NOMINMAX
			���windows.h��<algorithm>��ͬʱ����std::max()��std::min()�ĳ�ͻ��

	TUTORIAL_SHARED_PATH="G:/gitRepositories/ligIgl/libigl_CGAL_openGL/cmake/../external/../tutorial/data"

	CMAKE_INTDIR="Release"
*/


////////////////////////////////////////////////////////////////////////////// ����WINAPI��һЩ�ӿڣ�
namespace MY_WIN_API
{
	// ��ȡĳ��Ŀ¼�������ļ�����Ŀ¼����
	void getFileNames(std::string path, std::vector<std::string>& files, bool blRecur = true)
	{
		std::string str;
		struct _finddata_t fileinfo;			// �ļ���Ϣ����_finddata_t����ͷ�ļ�<io.h>
		intptr_t hFile = _findfirst(str.assign(path).append("/*").c_str(), &fileinfo);							// �ļ����	
		bool blFileValid = (hFile != -1);

		if (blFileValid)
		{
			do
			{
				bool isSubDir = (fileinfo.attrib & _A_SUBDIR);
				//�����Ŀ¼,�ݹ���ң��������,���ļ�����·������vector��
				if (isSubDir & blRecur)
				{
					if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
						getFileNames(str.assign(path).append("/").append(fileinfo.name), files);
				}
				else
					if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
						files.push_back(str.assign(path).append("/").append(fileinfo.name));
			} while (_findnext(hFile, &fileinfo) == 0);
			_findclose(hFile);
		}
	}


	// config�ļ�·���ַ��������� �� std::wstring�ַ�����
	static void GetConfigDir(std::wstring& confDir, const TCHAR* pszConfFile)
	{
#ifdef WIN32
		std::wstring strRet(pszConfFile);
		if ((std::string::npos != strRet.find(L"./")) ||
			(std::string::npos != strRet.find(L".\\")))
		{
			TCHAR szCurDir[256] = { 0 };
			GetCurrentDirectory(256, szCurDir);
			confDir = szCurDir;
			confDir += strRet.substr(1, strRet.length());

			return;
		}
		confDir = strRet;
#endif
		return;
	}


	// ��ȡconfig�ļ��е�ĳһ������Ϊfloat�ļ�ֵ
	float INIGetFloat(const TCHAR* pszKey, const TCHAR* pszConfFile)
	{
#ifdef WIN32
		std::wstring strFileName;
		GetConfigDir(strFileName, pszConfFile);
		TCHAR szValue[256] = { 0 };
		GetPrivateProfileString(L"ZHENGYADENTALCONFIG", pszKey, L"", szValue, 256, strFileName.c_str());

#ifdef UNICODE
		return _wtof(szValue);
#else
		return atof(szValue);
#endif
#else
		return 0.0f;
#endif

	}


	// ��ȡconfig�ļ��е�ĳһ������Ϊint�ļ�ֵ
	DWORD INIGetInt(const TCHAR* pszKey, const TCHAR* pszConfFile)
	{
#ifdef WIN32
		std::wstring strFileName;
		GetConfigDir(strFileName, pszConfFile);
		return GetPrivateProfileInt(L"ZHENGYADENTALCONFIG", pszKey, 0, strFileName.c_str());
#else
		return 0;
#endif
	}


}
using namespace MY_WIN_API;



////////////////////////////////////////////////////////////////////////////////////////////// DEBUG �ӿ�
namespace MY_DEBUG
{
	static std::string g_debugPath = "E:/";


	static void debugDisp()			// �ݹ���ֹ
	{						//		�ݹ���ֹ��Ϊ�޲λ���һ�����������ζ����ԡ�
		std::cout << std::endl;
		return;
	}

	template <typename T, typename... Types>
	static void debugDisp(const T& firstArg, const Types&... args)
	{
		std::cout << firstArg << " ";
		debugDisp(args...);
	}


	template <typename T, int M, int N>
	static void dispData(const Eigen::Matrix<T, M, N>& m)
	{
		auto dataPtr = m.data();
		unsigned elemsCount = m.size();

		for (unsigned i = 0; i < elemsCount; ++i)
			std::cout << dataPtr[i] << ", ";

		std::cout << std::endl;
	}


	template <typename Derived>
	static void dispData(const Eigen::PlainObjectBase<Derived>& m)
	{
		int m0 = m.RowsAtCompileTime;
		int n0 = m.ColsAtCompileTime;

		auto dataPtr = m.data();
		unsigned elemsCount = m.size();

		for (unsigned i = 0; i < elemsCount; ++i)
			std::cout << dataPtr[i] << ", ";

		std::cout << std::endl;
	}


	template <typename Derived>
	static void dispElem(const Eigen::MatrixBase<Derived>& m)
	{
		const Derived& mm = m.derived();
		std::cout << mm(1, 1) << std::endl;
	}


	template<typename DerivedV>
	static void debugWriteVers(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat(path, vers);
	}

	template<typename DerivedV>
	static void debugWriteVers2D(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat2D(path, vers);
	}


	template<typename DerivedV>
	static void debugWriteMesh(const char* name, const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteMeshMat(path, vers, tris);
	}


	template<typename DerivedV>
	static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteEdgesMat(path, edges, vers);
	}
}
using namespace MY_DEBUG;




////////////////////////////////////////////////////////////////////////////// TEST: ���񾫼�
namespace DECIMATION
{
#define IGL_COLLAPSE_EDGE_NULL -1

	// 
	/*
		std::get<0>(q) == A;
		std::get<0>(q) == d0 * n;
		std::get<0>(q) == d0 ^ 2;
	*/
	using Quadric = std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>;


	// q���������������
	Quadric operator+(const Quadric& q1, const Quadric& q2)
	{
		Quadric result;
		std::get<0>(result) = (std::get<0>(q1) + std::get<0>(q2)).eval();
		std::get<1>(result) = (std::get<1>(q1) + std::get<1>(q2)).eval();
		std::get<2>(result) = (std::get<2>(q1) + std::get<2>(q2));
		return result;
	}

	// lambda�����õ�ĳһ����Ƭ��quadric
	/*
		 Inputs:
		   va		1 by n row point on the subspace
		   S		m by n matrix where rows coorespond to orthonormal spanning
					vectors of the subspace to which we're measuring distance (usually a plane,  m=2)
		   weight  scalar weight
		 Returns quadric triple {A, b, c} so that A-2*b+c measures the quadric


		 Weight face's quadric (v'*A*v + 2*b'*v + c) by area
	*/
	Quadric subspace_quadric(const Eigen::RowVector3d& va, const Eigen::RowVector3d& dir_ab, const Eigen::RowVector3d& dir_h, const double weight)
	{
		Eigen::MatrixXd A = Eigen::MatrixXd::Identity(3, 3) - dir_ab.transpose() * dir_ab - dir_h.transpose() * dir_h;
		Eigen::RowVector3d b = -va + va.dot(dir_ab) * dir_ab + va.dot(dir_h) * dir_h;
		double c = va.dot(va) - pow(va.dot(dir_ab), 2) - pow(va.dot(dir_h), 2);

		return Quadric{ weight * A,  weight * b,  weight * c };
	};


	// ���㶥���Q����
	void getQuadrics(const Eigen::MatrixXd& vers, const Eigen::MatrixXi& tris,
		const Eigen::MatrixXi& edgeUeInfo, const Eigen::MatrixXi& UeTrisInfo,
		const Eigen::MatrixXi& EI, std::vector<std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>>& quadrics)
	{
		using namespace std;
		using Quadric = std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>;
		const int versCount = vers.rows();
		Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);

		// ����Quadrics��ʼ��Ϊ0��
		quadrics.resize(versCount, Quadric{ Eigen::MatrixXd::Zero(3, 3), Eigen::RowVectorXd::Zero(3), 0 });

		// Rather initial with zeros,  initial with a small amount of energy pull toward original vertex position
		const double w = 1e-10;
		for (int i = 0; i < vers.rows(); i++)
		{
			std::get<0>(quadrics[i]) = w * I;
			Eigen::RowVectorXd ver = vers.row(i);
			std::get<1>(quadrics[i]) = -w * ver;
			std::get<2>(quadrics[i]) = w * ver.dot(ver);
		}

		// Generic nD qslim from "Simplifying Surfaces with Color and Texture using Quadric Error Metric" (follow up to original QSlim)
		for (int i = 0; i < tris.rows(); i++)
		{
			// f1. don't know what the fuck is...
			int infinite_corner = -1;
			for (int k = 0; k < 3; k++)
			{
				if (std::isinf(vers(tris(i, k), 0)) || std::isinf(vers(tris(i, k), 1)) || std::isinf(vers(tris(i, k), 2)))
				{
					assert(infinite_corner == -1 && "Should only be one infinite corner");
					infinite_corner = k;
				}
			}

			// f2. 
			int vaIdx = tris(i, 0);
			int vbIdx = tris(i, 1);
			int vcIdx = tris(i, 2);
			if (infinite_corner == -1)
			{
				// Finite (non-boundary) face
				Eigen::RowVector3d va = vers.row(vaIdx);
				Eigen::RowVector3d vb = vers.row(vbIdx);
				Eigen::RowVector3d vc = vers.row(vcIdx);
				Eigen::RowVector3d ab = vb - va;
				Eigen::RowVector3d ac = vc - va;

				double area = (ab.cross(ac)).norm();										// ����Ƭ���������;

				Eigen::RowVector3d dir_ab = ab.normalized();
				Eigen::RowVector3d dir_h = (ac - dir_ab.dot(ac) * dir_ab).normalized();
				Quadric face_quadric = subspace_quadric(va, dir_ab, dir_h, area);

				//// for debug:
				//if (37640 == vaIdx || 37640 == vbIdx || 37640 == vcIdx)
				//{ 
				//	Eigen::MatrixXd A = std::get<0>(face_quadric) / area;
				//	Eigen::RowVectorXd b = std::get<1>(face_quadric) / area;
				//	double c = std::get<2>(face_quadric) / area;

				//	double d0 = std::sqrt(c);
				//	Eigen::RowVector4d p{ b(0), b(1), b(2), c };
				//	p = (p / d0).eval();
				//	debugDisp("p == ", p);
				//	debugDisp("area == ", area);
				//}

				// for debug:
				if (36394 == vaIdx || 36394 == vbIdx || 36394 == vcIdx)
				{
					Eigen::MatrixXd A = std::get<0>(face_quadric) / area;
					Eigen::RowVectorXd b = std::get<1>(face_quadric) / area;
					double c = std::get<2>(face_quadric) / area;

					double d0 = std::sqrt(c);
					Eigen::RowVector4d p{ b(0), b(1), b(2), c };
					p = (p / d0).eval();
					debugDisp("p == ", p);
					debugDisp("area == ", area);
				}

				// Throw at each corner 
				quadrics[vaIdx] = quadrics[vaIdx] + face_quadric;
				quadrics[vbIdx] = quadrics[vbIdx] + face_quadric;
				quadrics[vcIdx] = quadrics[vcIdx] + face_quadric;
			}
			else
			{
#if 0
				// cth corner is infinite --> edge opposite cth corner is boundary				
				const Eigen::RowVectorXd va = vers.row(tris(i, (infinite_corner + 1) % 3));		// Boundary edge vector
				Eigen::RowVectorXd ev = vers.row(tris(i, (infinite_corner + 2) % 3)) - va;
				const double length = ev.norm();
				ev /= length;

				// Face neighbor across boundary edge
				int e = edgeUeInfo(i + tris.rows() * infinite_corner);
				int opp = UeTrisInfo(e, 0) == i ? 1 : 0;
				int n = UeTrisInfo(e, opp);
				int nc = EI(e, opp);
				assert(
					((tris(i, (infinite_corner + 1) % 3) == tris(n, (nc + 1) % 3) &&
						tris(i, (infinite_corner + 2) % 3) == tris(n, (nc + 2) % 3)) ||
						(tris(i, (infinite_corner + 1) % 3) == tris(n, (nc + 2) % 3)
							&& tris(i, (infinite_corner + 2) % 3) == tris(n, (nc + 1) % 3))) &&
					"Edge flaps not agreeing on shared edge");

				// Edge vector on opposite face
				const Eigen::RowVectorXd eu = vers.row(tris(n, nc)) - va;
				assert(!std::isinf(eu(0)));

				// Matrix with vectors spanning plane as columns
				Eigen::MatrixXd A(ev.size(), 2);
				A << ev.transpose(), eu.transpose();

				// Use QR decomposition to find basis for orthogonal space
				Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
				const Eigen::MatrixXd Q = qr.householderQ();
				const Eigen::MatrixXd N = Q.topRightCorner(ev.size(), ev.size() - 2).transpose();
				assert(N.cols() == ev.size());
				assert(N.rows() == ev.size() - 2);

				Eigen::MatrixXd S(N.rows() + 1, ev.size());
				S << ev, N;
				Quadric boundary_edge_quadric = subspace_quadric(va, S, length * length);
				for (int k = 0; k < 3; k++)
					if (k != infinite_corner)
						quadrics[tris(i, k)] = quadrics[tris(i, k)] + boundary_edge_quadric;
#endif
				assert("should not happen.");
			}
		}
	}


	void circulation(const int uEdgeIdx, const bool ccw, const Eigen::MatrixXi& tris, const Eigen::VectorXi& edgeUeInfo,
		const Eigen::MatrixXi& UeTrisInfo, const Eigen::MatrixXi& UeCornersInfo, std::vector<int>& nbrVersIdx, std::vector<int>& nbrTrisIdx)
	{
		nbrVersIdx.clear(); nbrVersIdx.reserve(10);
		nbrTrisIdx.clear(); nbrTrisIdx.reserve(10);
		const int m = edgeUeInfo.size() / 3;
		assert(m * 3 == edgeUeInfo.size());

		const auto& step = [&](const int uEdgeIdx, const int ff, int& ne, int& rv, int& nf)
		{
			assert((UeTrisInfo(uEdgeIdx, 1) == ff || UeTrisInfo(uEdgeIdx, 0) == ff) && "uEdgeIdx should touch ff");
			const int nside = UeTrisInfo(uEdgeIdx, 0) == ff ? 1 : 0;
			const int nv = UeCornersInfo(uEdgeIdx, nside);
			nf = UeTrisInfo(uEdgeIdx, nside);
			const int dir = ccw ? -1 : 1;
			rv = tris(nf, nv);
			ne = edgeUeInfo(nf + m * ((nv + dir + 3) % 3));
		};

		const int triIdx = UeTrisInfo(uEdgeIdx, 0);
		int cTriIdx = triIdx;
		int cUeIdx = uEdgeIdx;

		while (true)
		{
			int re, rv;
			step(cUeIdx, cTriIdx, cUeIdx, rv, cTriIdx);
			nbrTrisIdx.push_back(cTriIdx);
			nbrVersIdx.push_back(rv);

			if (cTriIdx == triIdx)
			{
				assert(cUeIdx == uEdgeIdx);
				break;
			}
		}
	}


	void getAllTheShit(const int uEdgeIdx, const Eigen::MatrixXi& uEdges, const Eigen::MatrixXd& vers, const Eigen::MatrixXi& tris, \
		std::vector<int>& nbrVersIdx_src, std::vector<int>& nbrTrisIdx_src, std::vector<int>& nbrVersIdx_des, std::vector<int>& nbrTrisIdx_des)
	{
		auto blEdgeInTri = [&](const int vaIdx0, const int vbIdx0, const int triIdx0) -> bool
		{
			if (triIdx0 < 0)
				return false;
			const int A = tris(triIdx0, 0);
			const int B = tris(triIdx0, 1);
			const int C = tris(triIdx0, 2);
			if (A == vaIdx0 && B == vbIdx0)
				return true;
			if (B == vaIdx0 && C == vbIdx0)
				return true;
			if (C == vaIdx0 && A == vbIdx0)
				return true;

			return false;
		};

		auto getEdgeArrow = [&](const int triIdx, const int centerIdx)->std::pair<int, int>
		{
			if (triIdx < 0 || centerIdx < 0)
				return std::make_pair(-1, -1);

			const int A = tris(triIdx, 0);
			const int B = tris(triIdx, 1);
			const int C = tris(triIdx, 2);
			if (centerIdx == A)
				return std::make_pair(B, C);
			if (centerIdx == B)
				return std::make_pair(C, A);
			if (centerIdx == C)
				return std::make_pair(A, B);

			return std::make_pair(-1, -1);
		};

		const int vaIdx = uEdges(uEdgeIdx, 0);
		const int vbIdx = uEdges(uEdgeIdx, 1);
		std::unordered_set<int> nbrTriIdxesA = oneRingTriIdxes(vaIdx, vers, tris);
		std::unordered_set<int> nbrTriIdxesB = oneRingTriIdxes(vbIdx, vers, tris);

		int arrowStart;
		std::unordered_map<int, std::pair<int, int>> srcTrisMap;				// ��BΪ���ĵ��һȦ����Ƭ��
		std::unordered_map<int, std::pair<int, int>> desTrisMap;				// ��AΪ���ĵ��һȦ����Ƭ��

		arrowStart = vaIdx;
		for (const auto& triIdx : nbrTriIdxesB)
		{
			std::pair<int, int> edgeArrow = getEdgeArrow(triIdx, vbIdx);
			if (edgeArrow.first < 0)
				continue;
			srcTrisMap.insert(std::make_pair(edgeArrow.first, std::make_pair(edgeArrow.second, triIdx)));
		}
		while (true)
		{
			auto retIter = srcTrisMap.find(arrowStart);
			if (retIter == srcTrisMap.end())
				break;
			std::pair<int, int> triInfo = retIter->second;
			nbrTrisIdx_src.push_back(triInfo.second);
			arrowStart = triInfo.first;
			nbrVersIdx_src.push_back(arrowStart);
			srcTrisMap.erase(retIter);
		}

		arrowStart = vbIdx;
		for (const auto& triIdx : nbrTriIdxesA)
		{
			std::pair<int, int> edgeArrow = getEdgeArrow(triIdx, vaIdx);
			if (edgeArrow.first < 0)
				continue;
			desTrisMap.insert(std::make_pair(edgeArrow.first, std::make_pair(edgeArrow.second, triIdx)));
		}
		while (true)
		{
			auto retIter = desTrisMap.find(arrowStart);
			if (retIter == desTrisMap.end())
				break;
			nbrVersIdx_des.push_back(arrowStart);
			std::pair<int, int> triInfo = retIter->second;
			nbrTrisIdx_des.push_back(triInfo.second);
			arrowStart = triInfo.first;
			desTrisMap.erase(retIter);
		}
		std::reverse(nbrTrisIdx_des.begin(), nbrTrisIdx_des.end());
		std::reverse(nbrVersIdx_des.begin(), nbrVersIdx_des.end());
	}


	// ����ߵ��۵�cost���Լ��۵�����¶��㡪����򵥵Ĳ��ԣ�costΪ�߳����۵�����¶���Ϊ�е㣻
	template <typename T>
	double costAndNewPos_simplest(Eigen::Matrix<T, 1, 3>& newPos, const int ueIdx, const Eigen::MatrixXi& uEdges, \
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
	{
		int vaIdx = uEdges(ueIdx, 0);
		int vbIdx = uEdges(ueIdx, 1);
		double cost = (vers.row(vaIdx) - vers.row(vbIdx)).norm();
		newPos = (vers.row(vaIdx) + vers.row(vbIdx)) / 2;

		// Force infs and nans to infinity
		if (std::isinf(cost) || cost != cost)
		{
			cost = std::numeric_limits<double>::infinity();
			newPos.setConstant(0);
		}

		return cost;
	}


	// ����ߵ��۵�cost���Լ��۵�����¶��㡪��qslim���ԣ�
	template <typename T>
	double costAndNewPos_qslim(Eigen::Matrix<T, 1, 3>& newPos, const int ueIdx, const Eigen::MatrixXi& uEdges, \
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris, const std::vector<Quadric>& quadrics)
	{
		// Quadric: newPos'Ap + 2b'newPos + c,  optimal point: Ap = -b, or rather because we have row vectors: pA=-b
		Quadric quadric_p = quadrics[uEdges(ueIdx, 0)] + quadrics[uEdges(ueIdx, 1)];		// Combined quadric
		const auto& A = std::get<0>(quadric_p);
		const auto& b = std::get<1>(quadric_p);
		const auto& c = std::get<2>(quadric_p);

		// for debug
#if 1
		newPos = -b * A.inverse();
		double cost = newPos.dot(newPos * A) + 2 * newPos.dot(b) + c;
#else
		Eigen::Matrix4d Qmat;
		Qmat.block<3, 3>(0, 0) = A;
		double d0 = std::sqrt(c);
		Qmat(3, 3) = c;
		Qmat(0, 3) = b(0);
		Qmat(1, 3) = b(1);
		Qmat(2, 3) = b(2);
		Qmat(3, 0) = b(0);
		Qmat(3, 1) = b(1);
		Qmat(3, 2) = b(2);

		// �����۵��󶥵㣺
		Eigen::MatrixX4d Q_solve = Qmat;
		Eigen::Vector4d new_vec;
		Q_solve(3, 0) = 0.0, Q_solve(3, 1) = 0.0, Q_solve(3, 2) = 0.0, Q_solve(3, 3) = 1.0;
		double det = Q_solve.determinant();					// Q_solve������ʽ��
		if (std::abs(det) < 1e-6)
		{
			// ��Q�������ȣ�ʹ�ñߵ��е���Ϊ�¶���
			int vaIdx = uEdges(ueIdx, 0);
			int vbIdx = uEdges(ueIdx, 1);
			newPos = (vers.row(vaIdx) + vers.row(vbIdx)) / 2.0;
		}
		else
		{
			// ��Q��������: new_vec = Q_solve.inverse() * {0, 0, 0, 1};
			Eigen::Vector4d tmpArrow = { 0.0, 0.0, 0.0, 1.0 };
			new_vec = Q_solve.inverse() * tmpArrow;
			newPos = Eigen::Matrix<T, 1, 3>{ new_vec[0], new_vec[1], new_vec[2] };
		}

		// ����costֵ��
		double cost = new_vec.transpose() * Qmat * new_vec;
#endif

		// Force infs and nans to infinity
		if (std::isinf(cost) || cost != cost)
		{
			cost = std::numeric_limits<double>::infinity();
			newPos.setConstant(0);
		}

		return cost;
	}


	// �жϱ��Ƿ�����۵���
	bool edge_collapse_is_valid(std::vector<int>& srcNbrIdx, std::vector<int>& desNbrIdx)
	{
		// Do we really need to check if edge is IGL_COLLAPSE_EDGE_NULL ?
		if (srcNbrIdx.size() < 2 || desNbrIdx.size() < 2)
			return false;

		// determine if the first two vertices are the same before reordering.
		// If they are and there are 3 each, then (I claim) this is an edge on a single tet.
		const bool first_two_same = (srcNbrIdx[0] == desNbrIdx[0]) && (srcNbrIdx[1] == desNbrIdx[1]);
		if (srcNbrIdx.size() == 3 && desNbrIdx.size() == 3 && first_two_same)
			return false;           // single tet


		  // https://stackoverflow.com/a/19483741/148668
		std::sort(srcNbrIdx.begin(), srcNbrIdx.end());
		std::sort(desNbrIdx.begin(), desNbrIdx.end());
		std::vector<int> Nint;
		std::set_intersection(srcNbrIdx.begin(), srcNbrIdx.end(), desNbrIdx.begin(), desNbrIdx.end(), std::back_inserter(Nint));

		// check if edge collapse is valid: intersection of vertex neighbors of s and d should be exactly 2+(s,d) = 4

		// http://stackoverflow.com/a/27049418/148668
		if (Nint.size() != 2)
			return false;


		return true;
	}


	// ���������һ��������񡪡������ڱ�Ե����ߣ������һ������Զ����������������Ƭ
	template <typename T>
	bool connect_boundary_to_infinity(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& versOut, Eigen::MatrixXi& trisOut,
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
	{
		// 1. ����ĩβ����һ������Զ���㣺
		versOut.resize(vers.rows() + 1, vers.cols());
		versOut.topLeftCorner(vers.rows(), vers.cols()) = vers;
		auto inf = std::numeric_limits<T>::infinity();
		versOut.row(vers.rows()).setConstant(inf);

		// 2. �����б�Ե����ߺ�����Զ���������γ��µ�����Ƭ��
		Eigen::MatrixXi bdrys, newTris;
		bdryEdges(bdrys, tris);
		if (0 == bdrys.size())
		{
			versOut = vers;
			trisOut = tris;
			return true;
		}

		newTris.resize(bdrys.rows(), 3);
		newTris.col(0) = bdrys.col(1);
		newTris.col(1) = bdrys.col(0);
		newTris.col(2).array() = vers.rows();

		trisOut = tris;
		matInsertRows(trisOut, newTris);

		return true;
	}


	// ��ӡ���ȶ��ж��׵�����Ԫ�أ�
	void dispCosts(const igl::min_heap<std::tuple<double, int, int>>& pQueue, const unsigned num)
	{
		igl::min_heap<std::tuple<double, int, int>> pCopy = pQueue;
		std::vector<double> costValues;
		std::vector<std::vector<int>> edges;
		costValues.reserve(num);
		edges.reserve(num);

		for (unsigned i = 0; i < num; ++i)
		{
			auto& t = pCopy.top();
			costValues.push_back(std::get<0>(t));
			edges.push_back(std::vector<int>{std::get<1>(t), std::get<2>(t)});
			pCopy.pop();
		}

		for (unsigned i = 0; i < num; ++i)
			std::cout << costValues[i] << ", ";

		std::cout << std::endl << std::endl;
	}


	// �۵�����ΪuEdgeIdx������ߣ�
	bool collapseSingleEdge(const int uEdgeIdx, const Eigen::RowVectorXd& collapsedVer,
		std::vector<int>& nbrVersIdx_src, const std::vector<int>& nbrTrisIdx_src,
		std::vector<int>& nbrVersIdx_des, const std::vector<int>& nbrTrisIdx_des,
		Eigen::MatrixXd& vers, Eigen::MatrixXi& tris, Eigen::MatrixXi& uEdges, Eigen::MatrixXi& edges, \
		Eigen::VectorXi& edgeUeInfo, Eigen::MatrixXi& UeTrisInfo, Eigen::MatrixXi& UeCornersInfo,
		int& a_e1, int& a_e2, int& a_f1, int& a_f2)
	{
		// lambda����ĳһ������ߵ�������Ϣ����Ϊ��Ч��Ϣ��
		const auto& kill_edge = [&uEdges, &UeCornersInfo, &UeTrisInfo](const int uEdgeIdx)
		{
			uEdges(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			uEdges(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
			UeTrisInfo(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			UeTrisInfo(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
			UeCornersInfo(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			UeCornersInfo(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
		};

		// ע�⣡����ǰ��������߶�������ģ���vaIdx < vbIdx; ����������߶���Ӧ����������ߣ���Ϊ��Ե���ѱ�������

	   // ��������ߵ������յ㡪������С��Ϊ��㣬�������Ϊ�յ㣻
		const int trisCount = tris.rows();
		const int versCount = vers.rows();
		const int vaIdx = uEdges(uEdgeIdx, 0);
		const int vbIdx = uEdges(uEdgeIdx, 1);
		const std::vector<int>& nbrTriIdxes = nbrTrisIdx_src;

		// 0. �жϸñ��Ƿ���۵�
		if (!edge_collapse_is_valid(nbrVersIdx_src, nbrVersIdx_des))
			return false;

		// 1. Ҫ�۵��ı����˵㶼�滻ΪcollapsedVer
		vers.row(vaIdx) = collapsedVer;
		vers.row(vbIdx) = collapsedVer;

		int corner = UeCornersInfo(uEdgeIdx, 0);
		int triIdx = UeTrisInfo(uEdgeIdx, 0);
		int eIdx = corner * trisCount + triIdx;
		int corner1 = (corner + 1) % 3;
		int corner2 = (corner + 2) % 3;
		int eIdx1 = triIdx + trisCount * corner1;
		int eIdx2 = triIdx + trisCount * corner2;
		const int ueIdx1 = edgeUeInfo(eIdx1);
		const int ueIdx2 = edgeUeInfo(eIdx2);

		{
			const bool flip1 = UeTrisInfo(ueIdx1, 1) == triIdx;
			const int triOppIdx1 = flip1 ? UeTrisInfo(ueIdx1, 0) : UeTrisInfo(ueIdx1, 1);
			// const int triOppIdx1 = eOppIdx1 % trisCount;
			const int cornerOpp1 = flip1 ? UeCornersInfo(ueIdx1, 0) : UeCornersInfo(ueIdx1, 1);
			const int eOppIdx1 = triOppIdx1 + trisCount * cornerOpp1;
			//const int eOppIdx1 = eOppIdx1;

			// kill edge, kill tri:
			kill_edge(ueIdx1);
			tris(triIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
			tris(triIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
			tris(triIdx, 2) = IGL_COLLAPSE_EDGE_NULL;

			// �������������ߵ�����ӳ�䣻
			edgeUeInfo(eOppIdx1) = ueIdx2;

			// ����������
			const int opp2 = (UeTrisInfo(ueIdx2, 0) == triIdx ? 0 : 1);
			UeTrisInfo(ueIdx2, opp2) = triOppIdx1;
			UeCornersInfo(ueIdx2, opp2) = cornerOpp1;

			a_e1 = ueIdx1;
			a_f1 = triIdx;
		}

		int CORNER = UeCornersInfo(uEdgeIdx, 1);
		int TRIIDX = UeTrisInfo(uEdgeIdx, 1);
		int EIDX = CORNER * trisCount + TRIIDX;
		int CORNER1 = (CORNER + 2) % 3;
		int CORNER2 = (CORNER + 1) % 3;
		int EIDX1 = TRIIDX + trisCount * CORNER1;
		int EIDX2 = TRIIDX + trisCount * CORNER2;
		const int UEIDX1 = edgeUeInfo(EIDX1);
		const int UEIDX2 = edgeUeInfo(EIDX2);

		{
			const bool flip1 = UeTrisInfo(UEIDX1, 1) == TRIIDX;
			const int TRIOPPIDX1 = flip1 ? UeTrisInfo(UEIDX1, 0) : UeTrisInfo(UEIDX1, 1);
			const int CORNEROPP1 = flip1 ? UeCornersInfo(UEIDX1, 0) : UeCornersInfo(UEIDX1, 1);
			const int EOPPIDX1 = TRIOPPIDX1 + trisCount * CORNEROPP1;

			kill_edge(UEIDX1);
			tris(TRIIDX, 0) = IGL_COLLAPSE_EDGE_NULL;
			tris(TRIIDX, 1) = IGL_COLLAPSE_EDGE_NULL;
			tris(TRIIDX, 2) = IGL_COLLAPSE_EDGE_NULL;

			edgeUeInfo(EOPPIDX1) = UEIDX2;

			const int opp2 = (UeTrisInfo(UEIDX2, 0) == TRIIDX ? 0 : 1);
			UeTrisInfo(UEIDX2, opp2) = TRIOPPIDX1;
			UeCornersInfo(UEIDX2, opp2) = CORNEROPP1;

			a_e2 = UEIDX1;
			a_f2 = TRIIDX;
		}

		// 3. ���۵��������uEdgeIdx��, vaIdx��������vbIdx��Ϊ�����㣻
		for (auto i : nbrTriIdxes)
		{
			for (int k = 0; k < 3; k++)
			{
				if (tris(i, k) == vbIdx)
				{
					const int eIdx1 = i + trisCount * ((k + 1) % 3);
					const int eIdx2 = i + trisCount * ((k + 2) % 3);
					const int ueIdx2 = edgeUeInfo(eIdx2);
					const int ueIdx1 = edgeUeInfo(eIdx1);
					if (vbIdx == uEdges(ueIdx1, 0))
						uEdges(ueIdx1, 0) = vaIdx;
					if (vbIdx == uEdges(ueIdx1, 1))
						uEdges(ueIdx1, 1) = vaIdx;

					// Skip if we just handled this edge (claim: this will be all except for the first non-trivial face)
					if (ueIdx2 != -1)
					{
						if (vbIdx == uEdges(ueIdx2, 0))
							uEdges(ueIdx2, 0) = vaIdx;
						if (vbIdx == uEdges(ueIdx2, 1))
							uEdges(ueIdx2, 1) = vaIdx;
					}

					tris(i, k) = vaIdx;
					break;
				}
			}
		}

		// 4. Finally,  "remove" this edge and its information
		kill_edge(uEdgeIdx);

		return true;
	}


	// ��򵥵ı��۵������㷨����ʹ�ñ߳��������۵���ģ�
	void test1()
	{
		Eigen::MatrixXd vers, versOut, versOri;
		Eigen::MatrixXi tris, trisOut, trisTmp, trisOri;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi edges, uEdges, UeTrisInfo, UeCornersInfo, ueEdgeInfo;
		Eigen::VectorXi timeStamps;
		Eigen::MatrixXd collapsedVers;
		Eigen::VectorXd costs;
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// ���ȶ��У�(����ߵ�cost, ���������, ʱ���)
		bool clean_finish = false;
		tiktok& tt = tiktok::getInstance();

		igl::readOBJ("E:/����/tooth.obj", versOri, trisOri);
		igl::writeOBJ("E:/meshIn.obj", versOri, trisOri);

		tt.start();
		const unsigned trisCountOri = trisOri.rows();
		int trisCountNew = trisCountOri;
		int tarTrisCount = std::round(trisCountOri / 3);

		// 0. ����Ƿ��з���������ߣ�����ֱ���˳���
		Eigen::MatrixXi nmnUedges;
		nonManifoldUEs(nmnUedges, trisOri);
		if (nmnUedges.size() > 0)
			return;

		// 1. ���������һ��������񡪡������ڱ�Ե����ߣ������һ������Զ����������������Ƭ
		if (!connect_boundary_to_infinity(vers, tris, versOri, trisOri))
			return;

		// 2. �����������Ϣ��
		if (!getEdges(edges, tris))
			return;
		if (!getUedges(uEdges, edgeUeInfo, ueEdgeInfo, edges))
			return;
		if (!getUeInfos(UeTrisInfo, UeCornersInfo, edgeUeInfo, uEdges, tris))
			return;

		// 3. ����ÿ������ߵ�costֵ���߳������Դ�Ϊ���ȼ��������ȶ���
		const unsigned ueCount = uEdges.rows();
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());
		collapsedVers.resize(uEdges.rows(), vers.cols());
		costs.resize(uEdges.rows());
		igl::parallel_for(ueCount, [&](const int ueIdx)
			{
				Eigen::RowVector3d newPos;
				double cost = costAndNewPos_simplest(newPos, ueIdx, uEdges, vers, tris);
				collapsedVers.row(ueIdx) = newPos;
				costs(ueIdx) = cost;
			},
			10000);
		for (int i = 0; i < ueCount; i++)
			pQueue.emplace(costs(i), i, 0);

		// 4. ���۵���ѭ����
		int uEdgeIdx0;							// ���ȶ����ײ��ıߣ�
		int ueIdx1, ueIdx02, triIdx1, triIdx02;
		unsigned loopCount = 0;
		while (true)
		{
			bool collapsed = true;
			std::tuple<double, int, int> edgeTuple;							// ���ȶ����ײ���Ԫ�أ�
			std::vector<int> nbrTrisIdx_src, nbrVersIdx_src;				// vaIdx��1���������Ƭ��������������
			std::vector<int>  nbrTrisIdx_des, nbrVersIdx_des;			// vbIdx��1���������Ƭ��������������

			// w1. ȡ����Ԫ�أ�����Ԫ�س��ӣ�
			while (true)
			{
				// ww1. ������Ϊ�գ��˳�ѭ����
				if (pQueue.empty())
					assert("���۵�����");

				// ww2. ȡ����Ԫ�أ�����Ԫ�س��ӣ�
				edgeTuple = pQueue.top();
				if (std::get<0>(edgeTuple) == std::numeric_limits<double>::infinity())
					assert("���ױߵ�cost�������");
				pQueue.pop();
				uEdgeIdx0 = std::get<1>(edgeTuple);          // ���׵����������;

				// ww3. Check if matches timestamp
				if (std::get<2>(edgeTuple) == timeStamps(uEdgeIdx0))
					break;
			}

			// w2. ���㵱ǰ�����˵�1����Ķ��㡢����Ƭ��
			getAllTheShit(uEdgeIdx0, uEdges, vers, tris, nbrVersIdx_src, nbrTrisIdx_src, nbrVersIdx_des, nbrTrisIdx_des);

			//	w3. �۵����׵ıߣ�
			collapsed = collapseSingleEdge(uEdgeIdx0, collapsedVers.row(uEdgeIdx0), nbrVersIdx_src, nbrTrisIdx_src, \
				nbrVersIdx_des, nbrTrisIdx_des, vers, tris, uEdges, edges, edgeUeInfo, UeTrisInfo, UeCornersInfo, \
				ueIdx1, ueIdx02, triIdx1, triIdx02);
			loopCount++;

			//	w4. post_collapses; Ĭ��������ʲô��������

			// w5. �۵�����֮�󣬸������timeStamp�� ������رߵ�costֵ
			if (collapsed)
			{
				// w5.1. Erase the two,  other collapsed uEdges by marking their timestamps as -1
				timeStamps(ueIdx1) = -1;
				timeStamps(ueIdx02) = -1;

				// w5.2.ȷ�����۵��ı�1�������������Ƭ
				std::vector<int> nbrTrisIdx;							// ���۵��ı�1�������������Ƭ������
				nbrTrisIdx.reserve(nbrTrisIdx_src.size() + nbrTrisIdx_des.size());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_src.begin(), nbrTrisIdx_src.end());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_des.begin(), nbrTrisIdx_des.end());
				std::sort(nbrTrisIdx.begin(), nbrTrisIdx.end());
				nbrTrisIdx.erase(std::unique(nbrTrisIdx.begin(), nbrTrisIdx.end()), nbrTrisIdx.end());

				// w5.3. ȷ�����۵��ı�1�������������� 
				std::vector<int> nbrUeIdxes;						// ���۵��ı�1��������������������
				nbrUeIdxes.reserve(3 * nbrTrisIdx.size());
				for (auto& triIdx : nbrTrisIdx)
				{
					if (tris(triIdx, 0) == IGL_COLLAPSE_EDGE_NULL || tris(triIdx, 1) == IGL_COLLAPSE_EDGE_NULL || tris(triIdx, 2) == IGL_COLLAPSE_EDGE_NULL)
						continue;

					for (int i = 0; i < 3; i++)
					{
						const int eIdx = i * tris.rows() + triIdx;
						const int ueIdx = edgeUeInfo(eIdx);
						nbrUeIdxes.push_back(ueIdx);
					}
				}
				std::sort(nbrUeIdxes.begin(), nbrUeIdxes.end());
				nbrUeIdxes.erase(std::unique(nbrUeIdxes.begin(), nbrUeIdxes.end()), nbrUeIdxes.end());             // ȥ���ظ�Ԫ�أ�

				// w5.4. ���±��۵��ı�1�������������ߵ���Ϣ����costֵ���۵��󶥵㡢ʱ�����
				for (auto& ueIdx : nbrUeIdxes)
				{
					// wf1. ������۵���costֵ�����۵���Ķ������꣺
					Eigen::RowVector3d newPos;
					double cost = costAndNewPos_simplest(newPos, ueIdx, uEdges, vers, tris);

					// wf2. ���µ�����ߵ�ʱ���+1
					timeStamps(ueIdx)++;

					// wf3. ���º������߱����²�����У�֮ǰ���������ʱ�����ƥ�������Ѿ���Ч��
					pQueue.emplace(cost, ueIdx, timeStamps(ueIdx));
					collapsedVers.row(ueIdx) = newPos;
				}
			}
			else
				assert("edge collapse failed.");

			// w6. �ж��Ƿ�������ֹ����stopping_condition�����������������۵�ѭ����
			if (collapsed)
			{
				// w6.1. ��������Ƭ��
				if (triIdx1 < trisCountOri)
					trisCountNew -= 1;
				if (triIdx02 < trisCountOri)
					trisCountNew -= 1;

				// w6.2. ������Ƭ���ﵽĿ�꣬���˳����۵�ѭ����
				bool stopConditionFlag = (trisCountNew <= (int)tarTrisCount);
				if (stopConditionFlag)
				{
					clean_finish = true;
					break;
				}
			}
			else
				assert("edge collapse failed.");		 // ���۵�ʧ�ܣ��˳�ѭ����
		}


		// 5. ɾ�����к���IGL_COLLAPSE_EDGE_NULL������Ƭ���Լ���������Զ�������Ƭ
		trisTmp.resize(tris.rows(), 3);
		int index = 0;
		for (int i = 0; i < tris.rows(); i++)
		{
			int vaIdx = tris(i, 0);
			int vbIdx = tris(i, 1);
			int vcIdx = tris(i, 2);
			bool continueFlag = false;
			Eigen::Matrix3d triVers;					// ��ǰ����Ƭ���������㣻
			double* coorPtr = nullptr;
			if (vaIdx == vbIdx || vbIdx == vcIdx || vcIdx == vaIdx)
				continue;
			triVers.row(0) = vers.row(vaIdx);
			triVers.row(1) = vers.row(vbIdx);
			triVers.row(2) = vers.row(vcIdx);
			coorPtr = triVers.data();
			for (int k = 0; k < 9; ++k)
			{
				if (std::isinf(*coorPtr++))
				{
					continueFlag = true;
					break;
				}
			}
			if (continueFlag)
				continue;
			trisTmp.row(index) = tris.row(i);
			index++;
		}
		trisTmp.conservativeResize(index, trisTmp.cols());					  // �����൱��shrink_to_fit(); 


		// 6. ɾ���������㣬����������Ƭ��
		std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, trisTmp);
		if (!removeIsoVers(versOut, trisOut, vers, trisTmp, isoVerIdxes))
			return;

		debugDisp("loopCount == ", loopCount);
		tt.endCout("�����ʱ��");
		igl::writeOBJ("E:/decimateOutput.obj", versOut, trisOut);
		std::cout << "finished." << std::endl;
	}


	// igl::qslim()
	void test0()
	{
		Eigen::MatrixXd vers, versOut;
		Eigen::MatrixXi tris, trisOut;
		std::string filePath = "E:/����/";
		std::string fileName = "jawBracket1";
		igl::readOBJ(filePath + fileName + std::string{ ".obj" }, vers, tris);

		int trisCount = tris.rows();
		int tarTrisCount = 37734;						// ����Ŀ������Ƭ 
		Eigen::VectorXi newOldTrisInfo;				// newOldTrisInfo[i]�Ǿ����������е�i������Ƭ��Ӧ��ԭ���������Ƭ������
		Eigen::VectorXi newOldVersInfo;
		if (!igl::qslim(vers, tris, tarTrisCount, versOut, trisOut, newOldTrisInfo, newOldVersInfo))
			std::cout << "igl::qslim() failed!!!" << std::endl;

		std::stringstream ss;
		ss << "E:/" << fileName << "_qslimOrigin_" << tarTrisCount << ".obj";
		igl::writeOBJ(ss.str(), versOut, trisOut);


		std::cout << "finished." << std::endl;
	}


	// ʹ��viewer����Ķ��������𲽲鿴��������еı��۵���
	void test000()
	{
		using namespace igl;
		Eigen::MatrixXd vers;
		Eigen:Eigen::MatrixXi tris;
		igl::readOBJ("E:/����/tmpMesh.obj", vers, tris);

		bool keepWorkingFlag = true;
		bool simplest_decimate_flag = false;
		bool qslim_decimate_flag = ~simplest_decimate_flag;

		// 1. ׼������
		Eigen::MatrixXd versCopy = vers;
		Eigen::MatrixXi trisCopy = tris;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// ���ȶ��У� 

		int vaIdx = -1;					// State variables keeping track of edge we just collapsed
		int vbIdx = -1;
		Eigen::VectorXi timeStamps;
		Eigen::MatrixXd collapsedVers;
		typedef std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double> Quadric;
		std::vector<Quadric> quadrics;				// ÿ�������Q����

		igl::edge_flaps(trisCopy, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());						// Could reserve with https://stackoverflow.com/a/29236236/148668        
		collapsedVers.resize(uEdges.rows(), versCopy.cols());				// If an edge were collapsed, we'd collapse it to these points:

		// ����Ƿ��з����αߣ�
		{
			Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> BF;
			Eigen::Array<bool, Eigen::Dynamic, 1> BE;
			if (!igl::is_edge_manifold(trisCopy, uEdges.rows(), edgeUeInfo, BF, BE))
				return;
		}

		decimate_cost_and_placement_callback cost_and_placement;
		decimate_pre_collapse_callback pre_collapse;
		decimate_post_collapse_callback post_collapse;

		if (simplest_decimate_flag)
		{
			// ��ȡ�����۵������㷨�ĺ�����
			cost_and_placement = shortest_edge_and_midpoint;
			decimate_trivial_callbacks(pre_collapse, post_collapse);    // ���������۵��е�pre_collapse��post_collapse�����ӡ���ʲô��������
		}
		else
		{
			//  ����ÿ�������Q����
			per_vertex_point_to_plane_quadrics(vers, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, quadrics);

			// ��ȡqslim�㷨�еĺ����ӣ�
			igl::qslim_optimal_collapse_edge_callbacks(uEdges, quadrics, vaIdx, vbIdx, cost_and_placement, pre_collapse, post_collapse);
		}

		decimate_stopping_condition_callback stopping_condition;

		int loopCount = 0;
		int num_collapsed;                                // ������ѭ���ļ�����

		const auto& exportCurrentMesh = [&]()
		{
			Eigen::MatrixXi trisTmp(trisCopy.rows(), 3);
			Eigen::VectorXi _1;
			Eigen::VectorXi newOldTrisInfo;
			Eigen::VectorXi newOldVersInfo;
			newOldTrisInfo.resize(trisCopy.rows());
			int m = 0;
			for (int i = 0; i < trisCopy.rows(); i++)
			{
				if (trisCopy(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
					trisCopy(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
					trisCopy(i, 2) != IGL_COLLAPSE_EDGE_NULL)
				{
					trisTmp.row(m) = trisCopy.row(i);
					newOldTrisInfo(m) = i;
					m++;
				}
			}
			trisTmp.conservativeResize(m, trisTmp.cols());              // �����൱��shrink_to_fit();
			newOldTrisInfo.conservativeResize(m);

			// 3.2. ɾ�������еĹ������㣺
			Eigen::MatrixXd versOut;
			Eigen::MatrixXi trisOut;
			igl::remove_unreferenced(vers, trisTmp, versOut, trisOut, _1, newOldVersInfo);
			char str[256];
			sprintf_s(str, 256, "E:/decimating_output_%d.obj", loopCount);
			igl::writeOBJ(str, versOut, trisOut);
		};

		// lambda�����������ݸ�λ
		const auto& reset = [&]()
		{
			tris = trisCopy;
			vers = versCopy;
			edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
			collapsedVers.resize(uEdges.rows(), vers.cols());
			Eigen::VectorXd costs(uEdges.rows());

			pQueue = {};            // https://stackoverflow.com/questions/2852140/priority-queue-clear-method
			timeStamps = Eigen::VectorXi::Zero(uEdges.rows());

			// r1. ����ÿ���ߵ��۵�costֵ���Լ��۵�֮��Ķ������꣬�������ȶ���pQueue��
			{
				Eigen::VectorXd costs(uEdges.rows());
				igl::parallel_for(uEdges.rows(), \
					[&](const int i)
					{
						double cost = i;
						Eigen::RowVectorXd edgeCenter(1, 3);           // ȡ�ߵ��е���Ϊ���۵�֮��Ķ��㣻
						cost_and_placement(i, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, cost, edgeCenter);
						collapsedVers.row(i) = edgeCenter;
						costs(i) = cost;
					}, \
					10000);

				for (int i = 0; i < uEdges.rows(); i++)
					pQueue.emplace(costs(i), i, 0);
			}

			num_collapsed = 0;
			viewer.data().clear();
			viewer.data().set_mesh(vers, tris);
			viewer.data().set_face_based(true);
		};


		// lambda����ִ�����񾫼�׼����Ⱦ�����ݣ�
		const auto& pre_draw = [&](igl::opengl::glfw::Viewer& viewer)->bool
		{
			// p1. ÿһ�ζ���ѭ���У�����1%�ıߣ�
			if (viewer.core().is_animating && !pQueue.empty())
			{
				bool FlagCollapsed = false;           // ����ѭ���б������Ƿ�ִ�гɹ���

				// p1.1 ִ�б���������collapse edge
				const int max_iter = std::ceil(0.01 * pQueue.size());
				for (int j = 0; j < max_iter; j++)
				{
					if (!collapse_edge(cost_and_placement, vers, tris, uEdges, \
						edgeUeInfo, UeTrisInfo, UeCornersInfo, pQueue, timeStamps, collapsedVers))              // collapse_edge()����2.1
						break;
					FlagCollapsed = true;
					num_collapsed++;
				}

				// p1.2 
				if (FlagCollapsed)
				{
					viewer.data().clear();
					viewer.data().set_mesh(vers, tris);
					viewer.data().set_face_based(true);
				}
			}

			// p2. �����۵������������Ƭ����������������û�б����ΪIGL_COLLAPSE_EDGE_NULL������Ƭ����
			unsigned currentTrisCount = 0;
			for (int i = 0; i < tris.rows(); i++)
			{
				if (tris(i, 0) != IGL_COLLAPSE_EDGE_NULL || tris(i, 1) != IGL_COLLAPSE_EDGE_NULL || tris(i, 2) != IGL_COLLAPSE_EDGE_NULL)
					currentTrisCount++;
			}

			if (loopCount < 10)
				exportCurrentMesh();
			std::cout << "loop : " << loopCount++ << ", current trisCount == " << currentTrisCount << std::endl;
			return false;
		};


		// lambda���������¼���Ӧ��
		const auto& key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)->bool
		{
			switch (key)
			{
			case ' ':
				viewer.core().is_animating ^= 1;
				break;
			case 'R':
			case 'r':
				reset();
				break;
			case '0':
				keepWorkingFlag = false;
				viewer.core().is_animating = 0;
				viewer.launch_shut();
				viewer.shutdown_plugins();
				break;
			default:
				return false;
			}
			return true;
		};


		// 1. ��ʼ������λ����
		reset();

		// 2. �򿪴��ڣ�����ѭ��
		viewer.core().background_color.setConstant(1);
		viewer.core().is_animating = true;                                 // ����
		viewer.core().animation_max_fps = 1.0;						// ָ����󶯻�֡�ʣ�
		viewer.callback_key_down = key_down;
		viewer.callback_pre_draw = pre_draw;
		viewer.launch(keepWorkingFlag);

		// 3. ������ֹ���������
		Eigen::VectorXi newOldTrisInfo;
		Eigen::VectorXi newOldVersInfo;

		// 3.1. ɾ�����к��б��ΪIGL_COLLAPSE_EDGE_NULL�ߵ�����Ƭ��
		Eigen::MatrixXi trisTmp(trisCopy.rows(), 3);
		Eigen::VectorXi _1;
		newOldTrisInfo.resize(trisCopy.rows());
		int m = 0;
		for (int i = 0; i < trisCopy.rows(); i++)
		{
			if (trisCopy(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
				trisCopy(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
				trisCopy(i, 2) != IGL_COLLAPSE_EDGE_NULL)
			{
				trisTmp.row(m) = trisCopy.row(i);
				newOldTrisInfo(m) = i;
				m++;
			}
		}
		trisTmp.conservativeResize(m, trisTmp.cols());              // �����൱��shrink_to_fit();
		newOldTrisInfo.conservativeResize(m);

		// 3.2. ɾ�������еĹ������㣺
		Eigen::MatrixXd versOut;
		Eigen::MatrixXi trisOut;
		igl::remove_unreferenced(vers, trisTmp, versOut, trisOut, _1, newOldVersInfo);
		igl::writeOBJ("E:/decimating_output.obj", versOut, trisOut);

		std::cout << "finished." << std::endl;
	}


	// ����qslim()
	void test0000()
	{
		Eigen::MatrixXd vers, versOut, vers0, versOri;
		Eigen::MatrixXi tris, trisOut, trisTmp, trisOri;
		Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]�Ǿ����������е�i������Ƭ��Ӧ��ԭ���������Ƭ������
		Eigen::VectorXi newOldVersInfo;
		Eigen::VectorXi edgeUeInfo;
		Eigen::MatrixXi edges, uEdges, ueEdgeInfo, UeTrisInfo, UeCornersInfo;
		Eigen::VectorXi timeStamps;
		Eigen::MatrixXd collapsedVers;
		Eigen::VectorXd costs;
		std::vector<Quadric> quadrics;													// ÿ�������Q����
		igl::min_heap<std::tuple<double, int, int> > pQueue;				// ���ȶ��У�
		Eigen::VectorXi _1, I2;											// ���ɾ����������ڲ�����Ƭʱ�õ���
		bool clean_finish = false;
		tiktok& tt = tiktok::getInstance();

		// for debug:
		std::vector<triplet<double>> versMoni;
		std::vector<triplet<int>>	trisMoni;
		std::vector<doublet<int>> edgesMoni;
		std::vector<int>			vecMoni;
		int retIdx = -1;

		tt.start();
		igl::readOBJ("E:/����/jawBracket1.obj", versOri, trisOri);
		tt.endCout("elapsed time of loading mesh is: ");
		igl::writeOBJ("E:/meshIn.obj", versOri, trisOri);

		unsigned versCount = versOri.rows();
		unsigned trisCount = trisOri.rows();
		int trisCountNew = trisCount;
		int tarTrisCount = 37734;

		// 0. ����Ƿ��з���������ߣ�����ֱ���˳���
		Eigen::MatrixXi nmnUedges;
		nonManifoldUEs(nmnUedges, trisOri);
		if (nmnUedges.size() > 0)
			return;

		// 1. ���������һ��������񡪡������ڱ�Ե����ߣ������һ������Զ����������������Ƭ
		connect_boundary_to_infinity(vers, tris, versOri, trisOri);

		// 1.1 �����������Ϣ��
		getEdges(edges, tris);
		getUedges(uEdges, edgeUeInfo, ueEdgeInfo, edges);
		getUeInfos(UeTrisInfo, UeCornersInfo, edgeUeInfo, uEdges, tris);

		// 2. ����ÿ�������Q����
		getQuadrics(vers, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, quadrics);

		// 3. ����ÿ������ߵ�costֵ���Դ�Ϊ���ȼ��������ȶ���
		timeStamps = Eigen::VectorXi::Zero(uEdges.rows());
		collapsedVers.resize(uEdges.rows(), vers.cols());
		costs.resize(uEdges.rows());
		igl::parallel_for(uEdges.rows(), [&](const int ueIdx)
			{
				// ������cost_and_placement()�����ݣ�
				double cost = ueIdx;
				Eigen::RowVectorXd p(1, 3);

				// Combined quadric
				Quadric quadric_p;
				quadric_p = quadrics[uEdges(ueIdx, 0)] + quadrics[uEdges(ueIdx, 1)];

				// Quadric: p'Ap + 2b'p + c,  optimal point: Ap = -b, or rather because we have row vectors: pA=-b
				const auto& A = std::get<0>(quadric_p);
				const auto& b = std::get<1>(quadric_p);
				const auto& c = std::get<2>(quadric_p);
				p = -b * A.inverse();
				cost = p.dot(p * A) + 2 * p.dot(b) + c;

				// Force infs and nans to infinity
				if (std::isinf(cost) || cost != cost)
				{
					cost = std::numeric_limits<double>::infinity();
					p.setConstant(0);
				}

				collapsedVers.row(ueIdx) = p;
				costs(ueIdx) = cost;
			},
			10000);

		// FOR DEBUG������ӡsick edge��
		double theCost = 0;
		int a, b;
		int theUeIdx = 0;
		for (int i = 0; i < uEdges.rows(); ++i)
		{
			a = uEdges(i, 0);
			b = uEdges(i, 1);

			if ((a == 37640 && b == 36394) || (a == 36394 && b == 37640))
			{
				theUeIdx = i;
				Eigen::RowVector3d verA = vers.row(a);
				Eigen::RowVector3d verB = vers.row(b);
				Eigen::RowVector3d newPos = collapsedVers.row(i);
				theCost = costs(i);
				//debugDisp("verA == ", verA);
				//debugDisp("verB == ", verB);
				//debugDisp("newPos == ", newPos);
				debugDisp("the cost == ", theCost);
				break;
			}
		}


		//// for debug������ӡsick edge���˵��quadric:
		//for (int i = 0; i < vers.rows(); ++i)
		//{
		   // if (i == 37640 || i == 36394)
		   // {
		   //	 Quadric q = quadrics[i];
		   //	 debugDisp("index == ", i);
		   //	 debugDisp("A == ", std::get<0>(q));
		   //	 debugDisp("vec == ", std::get<1>(q));
		   //	 debugDisp("\n\n\n");
		   // }
		//}


		for (int i = 0; i < uEdges.rows(); i++)
			pQueue.emplace(costs(i), i, 0);


		//// for debug������ӡǰʮ��cost
		//auto pQueueCopy = pQueue;
		//int order = 0; 
		//auto ttt = pQueueCopy.top();
		//while (order < 10) 
		//{
		//	debugDisp("cost == ", std::get<0>(ttt));
		//	pQueueCopy.pop();
		//	ttt = pQueueCopy.top();
		//	order++;
		//}
		//debugDisp("order == ", order);


		// 4. ���۵���ѭ����
		int uEdgeIdx0;							// ���ȶ����ײ��ıߣ�
		int ueIdx1, ueIdx02, triIdx1, triIdx02;
		while (true)
		{
			bool collapsed = true;
			std::tuple<double, int, int> edgeTuple;							// ���ȶ����ײ���Ԫ�أ�
			std::vector<int> nbrTrisIdx_src, nbrVersIdx_src;				// vaIdx��1���������Ƭ��������������
			std::vector<int>  nbrTrisIdx_des, nbrVersIdx_des;			// vbIdx��1���������Ƭ��������������

			// w1. ȡ����Ԫ�أ�����Ԫ�س��ӣ�
			while (true)
			{
				// ww1. ������Ϊ�գ��˳�ѭ����
				if (pQueue.empty())
					assert("���۵�����");

				// ww2. ȡ����Ԫ�أ�����Ԫ�س��ӣ�
				edgeTuple = pQueue.top();
				if (std::get<0>(edgeTuple) == std::numeric_limits<double>::infinity())
					assert("���ױߵ�cost�������");
				pQueue.pop();
				uEdgeIdx0 = std::get<1>(edgeTuple);          // ���׵����������;

				// ww3. Check if matches timestamp
				if (std::get<2>(edgeTuple) == timeStamps(uEdgeIdx0))
					break;
			}

			// w2. ���㵱ǰ�����˵�1����Ķ��㡢����Ƭ��
			int vaIdx = uEdges(uEdgeIdx0, 0);
			int vbIdx = uEdges(uEdgeIdx0, 1);
			int smallIdx = (vaIdx < vbIdx ? vaIdx : vbIdx);
			igl::circulation(uEdgeIdx0, true, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_src, nbrTrisIdx_src);
			igl::circulation(uEdgeIdx0, false, tris, edgeUeInfo, UeTrisInfo, UeCornersInfo, nbrVersIdx_des, nbrTrisIdx_des);

			//	w3. �۵����׵ıߣ�
			collapsed = collapseSingleEdge(uEdgeIdx0, collapsedVers.row(uEdgeIdx0), nbrVersIdx_src, nbrTrisIdx_src, \
				nbrVersIdx_des, nbrTrisIdx_des, vers, tris, uEdges, edges, edgeUeInfo, UeTrisInfo, UeCornersInfo, ueIdx1, ueIdx02, triIdx1, triIdx02);

			//	w4. post_collapses��������q����
			if (collapsed)
				quadrics[smallIdx] = quadrics[vaIdx] + quadrics[vbIdx];

			// w5. �۵�����֮�󣬸������timeStamp�� ������رߵ�costֵ
			if (collapsed)
			{
				// w5.1. Erase the two,  other collapsed uEdges by marking their timestamps as -1
				timeStamps(ueIdx1) = -1;
				timeStamps(ueIdx02) = -1;

				// w5.2.ȷ�����۵��ı�1�������������Ƭ
				std::vector<int> nbrTrisIdx;							// ���۵��ı�1�������������Ƭ������
				nbrTrisIdx.reserve(nbrTrisIdx_src.size() + nbrTrisIdx_des.size());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_src.begin(), nbrTrisIdx_src.end());
				nbrTrisIdx.insert(nbrTrisIdx.end(), nbrTrisIdx_des.begin(), nbrTrisIdx_des.end());
				std::sort(nbrTrisIdx.begin(), nbrTrisIdx.end());
				nbrTrisIdx.erase(std::unique(nbrTrisIdx.begin(), nbrTrisIdx.end()), nbrTrisIdx.end());

				// w5.3. ȷ�����۵��ı�1�������������� 
				std::vector<int> nbrUeIdxes;						// ���۵��ı�1��������������������
				nbrUeIdxes.reserve(3 * nbrTrisIdx.size());
				for (auto& triIdx : nbrTrisIdx)
				{
					if (tris(triIdx, 0) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 1) != IGL_COLLAPSE_EDGE_NULL ||
						tris(triIdx, 2) != IGL_COLLAPSE_EDGE_NULL)
					{
						for (int i = 0; i < 3; i++)
						{
							const int eIdx = i * tris.rows() + triIdx;
							const int ueIdx = edgeUeInfo(eIdx);
							nbrUeIdxes.push_back(ueIdx);
						}
					}
				}
				std::sort(nbrUeIdxes.begin(), nbrUeIdxes.end());
				nbrUeIdxes.erase(std::unique(nbrUeIdxes.begin(), nbrUeIdxes.end()), nbrUeIdxes.end());             // ȥ���ظ�Ԫ�أ�

				// w5.4. ���±��۵��ı�1�������������ߵ���Ϣ����costֵ���۵��󶥵㡢ʱ�����
				for (auto& ueIdx : nbrUeIdxes)
				{
					// wf1. ������۵���costֵ�����۵���Ķ������꣺
					Eigen::RowVector3d newPos;
					double cost = costAndNewPos_qslim(newPos, ueIdx, uEdges, vers, tris, quadrics);

					// wf2. ���µ�����ߵ�ʱ���+1
					timeStamps(ueIdx)++;

					// wf3. ���º������߱����²�����У�֮ǰ���������ʱ�����ƥ�������Ѿ���Ч��
					pQueue.emplace(cost, ueIdx, timeStamps(ueIdx));
					collapsedVers.row(ueIdx) = newPos;
				}
			}
			else
				assert("edge collapse failed.");

			// w6. �ж��Ƿ�������ֹ����stopping_condition�����������������۵�ѭ����
			if (collapsed)
			{
				// w6.1. ��������Ƭ��
				if (triIdx1 < trisCount)
					trisCountNew -= 1;
				if (triIdx02 < trisCount)
					trisCountNew -= 1;

				// w6.2. ������Ƭ���ﵽĿ�꣬���˳����۵�ѭ����
				bool stopConditionFlag = (trisCountNew <= (int)tarTrisCount);
				if (stopConditionFlag)
				{
					clean_finish = true;
					break;
				}
			}
			else
				assert("edge collapse failed.");		 // ���۵�ʧ�ܣ��˳�ѭ����
		}

		// 5. ɾ�����к���IGL_COLLAPSE_EDGE_NULL������Ƭ���Լ���������Զ�������Ƭ
		trisTmp.resize(tris.rows(), 3);
		int index = 0;
		for (int i = 0; i < tris.rows(); i++)
		{
			int vaIdx = tris(i, 0);
			int vbIdx = tris(i, 1);
			int vcIdx = tris(i, 2);
			bool continueFlag = false;
			Eigen::Matrix3d triVers;					// ��ǰ����Ƭ���������㣻
			double* coorPtr = nullptr;
			if (vaIdx == IGL_COLLAPSE_EDGE_NULL || vbIdx == IGL_COLLAPSE_EDGE_NULL || vcIdx == IGL_COLLAPSE_EDGE_NULL)
				continue;
			triVers.row(0) = vers.row(vaIdx);
			triVers.row(1) = vers.row(vbIdx);
			triVers.row(2) = vers.row(vcIdx);
			coorPtr = triVers.data();
			for (int k = 0; k < 9; ++k)
			{
				if (std::isinf(*coorPtr++))
				{
					continueFlag = true;
					break;
				}
			}
			if (continueFlag)
				continue;
			trisTmp.row(index) = tris.row(i);
			index++;
		}
		trisTmp.conservativeResize(index, trisTmp.cols());					  // �����൱��shrink_to_fit(); 

		// 6. ɾ���������㣬����������Ƭ��
		std::vector<unsigned> isoVerIdxes = checkIsoVers(vers, trisTmp);
		if (!removeIsoVers(versOut, trisOut, vers, trisTmp, isoVerIdxes))
			return;

		igl::writeOBJ("E:/qslimOutput.obj", versOut, trisOut);
		std::cout << "finished." << std::endl;
	}


	// �󾫼��������ԭʼ����Ľ���������QEM��paper�������ٶ�̫��
	double calcSimpApproxError(const Eigen::MatrixXd& versSimp, const Eigen::MatrixXi& trisSimp, \
		const Eigen::MatrixXd& versOri, const Eigen::MatrixXi& trisOri)
	{
		// ע�ⲻ�ܳ���̫��ľ��󣬷�������bad alloc

		double appErr = 0;
		std::pair<unsigned, float> pair0;
		int versCount0 = versOri.rows();
		int versCount1 = versSimp.rows();
		const Eigen::MatrixXd& versMat0 = versOri;
		const Eigen::MatrixXd& versMat1 = versSimp;
		const Eigen::MatrixXi& trisMat0 = trisOri;
		const Eigen::MatrixXi& trisMat1 = trisSimp;
		Eigen::MatrixXd planeCoeff0, planeCoeff1;
		bool ret0 = trianglesPlane(planeCoeff0, versMat0, trisMat0);
		bool ret1 = trianglesPlane(planeCoeff1, versMat1, trisMat1);

		std::vector<long double> squaredDis0, squaredDis1;
		squaredDis0.reserve(versCount0);
		squaredDis1.reserve(versCount1);

		// �㵽ƽ����룺 dis == p.dot(v)��ע�����з��ŵľ��룻 p == (a,b,c,d), v = (x0, y0, z0, 1), (a,b,c)��ƽ��Ĺ�һ����������
		Eigen::MatrixXd versExt0{ Eigen::MatrixXd::Ones(versCount0, 4) };
		Eigen::MatrixXd versExt1{ Eigen::MatrixXd::Ones(versCount1, 4) };
		versExt0.leftCols(3) = versMat0.array().cast<double>();
		versExt1.leftCols(3) = versMat1.array().cast<double>();

		// for new method:
		PARALLEL_FOR(0, versCount1, [&](int verIdx)
			{
				std::lock_guard<std::mutex> guard(g_mutex);
				Eigen::Vector4d verColVec = versExt1.row(verIdx).transpose();
				// ����ΪverIdx��simp�����еĶ��㣬��ori��������������Ƭƽ��ķ��ž��룺
				Eigen::VectorXd signedDises = planeCoeff0 * verColVec;
				Eigen::VectorXd sqrDises = signedDises.array() * signedDises.array();
				squaredDis1.push_back(sqrDises.minCoeff());
			});

		PARALLEL_FOR(0, versCount0, [&](int verIdx)
			{
				std::lock_guard<std::mutex> guard(g_mutex);
				Eigen::Vector4d verColVec = versExt0.row(verIdx).transpose();
				// ����ΪverIdx��simp�����еĶ��㣬��ori��������������Ƭƽ��ķ��ž��룺
				Eigen::VectorXd signedDises = planeCoeff1 * verColVec;
				Eigen::VectorXd sqrDises = signedDises.array() * signedDises.array();
				squaredDis0.push_back(sqrDises.minCoeff());
			});

		long double s0 = 0;
		long double s1 = 0;
		for (const auto& num : squaredDis0)
			s0 += num;
		for (const auto& num : squaredDis1)
			s1 += num;
		appErr = (s0 + s1) / (static_cast<long>(versCount0) + static_cast<long>(versCount1));

		return appErr;
	}


	// ������������Ľ�����
	void test2()
	{
		Eigen::MatrixXd vers1, vers2;
		Eigen::MatrixXi tris1, tris2;
		double appErr = 0;

		igl::readOBJ("E:/����.obj", vers1, tris1);
		igl::readOBJ("E:/����_remeshed.obj", vers2, tris2);
		appErr = calcSimpApproxError(vers2, tris2, vers1, tris1);
		std::cout << "appErr == " << appErr << std::endl;

		vers1.resize(0, 0);
		tris1.resize(0, 0);
		vers2.resize(0, 0);
		tris2.resize(0, 0);
		igl::readOBJ("E:/rootTooth1.obj", vers1, tris1);
		igl::readOBJ("E:/rootTooth1_remeshed.obj", vers2, tris2);
		appErr = calcSimpApproxError(vers2, tris2, vers1, tris1);
		std::cout << "appErr == " << appErr << std::endl;

		vers1.resize(0, 0);
		tris1.resize(0, 0);
		vers2.resize(0, 0);
		tris2.resize(0, 0);
		igl::readOBJ("E:/fatTeeth.obj", vers1, tris1);
		igl::readOBJ("E:/fatTeeth_remeshed.obj", vers2, tris2);
		appErr = calcSimpApproxError(vers2, tris2, vers1, tris1);
		std::cout << "appErr == " << appErr << std::endl;

		std::cout << "finished." << std::endl;
	}


	// ������ȡ��������ִ��qslim���۵�����
	int testCmd_qslimDecimation(int argc, char** argv)
	{
		tiktok& tt = tiktok::getInstance();
		CString   cPath, fileConfig;
		std::string path, pathOBJ, pathOutput;
		std::vector<std::string> fileNames, tmpStrVec, OBJfileNames;
		std::stringstream ss;
		bool debugFlag = false;
		bool blTarVersOrRatio = false;			// 0 - ���򵽰ٷֱȶ�����, 1 - ����ȷ��������
		int meshesCount = 0;
		int tarVersCount = 0;
		float deciRatio = 0;

		// 00. ��ȡ·����������
		{
			GetModuleFileName(NULL, cPath.GetBufferSetLength(MAX_PATH + 1), MAX_PATH);		// ��ȡ��ǰ���̼��ص�ģ���·����
			int nPos = cPath.ReverseFind('\\');
			cPath = cPath.Left(nPos);
			path = CT2CA{ cPath };
			pathOBJ = path + "\\inputOBJ";
			pathOutput = path + "\\outputData";
			fileConfig = cPath + "\\config.ini";

			// ��ȡ�����ļ��еĲ�����
			tarVersCount = INIGetInt(TEXT("decimationTarVers"), fileConfig);			// ����Ŀ�궥������
			deciRatio = INIGetFloat(TEXT("decimationRatio"), fileConfig);					// ����ٷֱȣ�
			unsigned tmpInt = INIGetInt(TEXT("deciFlagTarVersOrRatio"), fileConfig);
			blTarVersOrRatio = tmpInt > 0 ? true : false;
			unsigned debugInt = INIGetInt(TEXT("debugFlag"), fileConfig);
			debugFlag = debugInt > 0 ? true : false;
			if (debugFlag)
				debugDisp("Debug mode: ");

			getFileNames(pathOBJ.c_str(), tmpStrVec, false);
			meshesCount = tmpStrVec.size();
			OBJfileNames.reserve(meshesCount);
			for (const auto& str : tmpStrVec)
			{
				std::string tailStr = str.substr(str.size() - 4, 4);				//	".obj"
				if (".obj" == tailStr)
				{
					fileNames.push_back(str);
					unsigned index = str.find_last_of("/");
					std::string OBJfileName = str.substr(index, str.size() - index - 4);			// "/" + obj�ļ���������·����.obj��׺��
					OBJfileNames.push_back(OBJfileName);
				}
			}
		}

		// 0. ��ȡ��������
		tt.start();
		std::cout << "��ȡ��������..." << std::endl;
		std::vector<Eigen::MatrixXd> meshesVers;
		std::vector<Eigen::MatrixXi> meshesTris;
		{
			meshesVers.resize(meshesCount);
			meshesTris.resize(meshesCount);
			for (unsigned i = 0; i < meshesCount; ++i)
				objReadMeshMat(meshesVers[i], meshesTris[i], fileNames[i].c_str());
			tt.endCout("��ȡ���������ʱ��");
		}

		// 2. ִ��qslim����
		tt.start();
		debugDisp("ִ��qslim����");
		for (int i = 0; i < meshesCount; ++i)
		{
			Eigen::MatrixXd outVers;
			Eigen::MatrixXi outTris;
			Eigen::VectorXi newOldTrisInfo;						// newOldTrisInfo[i]�Ǿ����������е�i������Ƭ��Ӧ��ԭ���������Ƭ������
			Eigen::VectorXi newOldVersInfo;
			double tvRatio = static_cast<double>(meshesTris[i].rows()) / static_cast<double>(meshesVers[i].rows());
			int tarTrisCount = 0;
			if (blTarVersOrRatio)
				tarTrisCount = static_cast<int>(tarVersCount * tvRatio);
			else
				tarTrisCount = static_cast<int>(meshesTris[i].rows() * deciRatio);
			igl::qslim(meshesVers[i], meshesTris[i], tarTrisCount, outVers, outTris, newOldTrisInfo, newOldVersInfo);

			debugDisp(OBJfileNames[i], ".obj���񾫼���ɣ���������", meshesVers[i].rows(), "��", outVers.rows());
			if (debugFlag)
				debugWriteMesh((OBJfileNames[i] + std::string{ "_qslimDeci" }).c_str(), outVers, outTris);
			ss.str("");
			ss << pathOutput << OBJfileNames[i] << "_qslimDeci.obj";
			objWriteMeshMat(ss.str().c_str(), outVers, outTris);
		}
		tt.endCout("qslim�����ʱ��");

		debugDisp("finished.");
		getchar();

		return 0;
	}
}


 

////////////////////////////////////////////////////////////////////////////// ���ɿ���̨���򹤾ߣ�

// ������ȡ��������ִ��laplace��˳��
int testCmd_laplaceFaring(int argc, char** argv)
{
	tiktok& tt = tiktok::getInstance();
	float deciRatio = 0.5;							// ������

	// ����·����
	int   nPos;
	CString   cPath;
	GetModuleFileName(NULL, cPath.GetBufferSetLength(MAX_PATH + 1), MAX_PATH);		// ��ȡ��ǰ���̼��ص�ģ���·����
	nPos = cPath.ReverseFind('\\');
	cPath = cPath.Left(nPos);
	std::string path{ CT2CA{cPath} };
	std::string pathOBJ = path + "\\inputOBJ";
	std::string pathOutput = path + "\\outputData";
	CString fileConfig = cPath + "\\config.ini";

	// 1. ��ȡ�������� 
	tt.start();
	std::cout << "��ȡ��������..." << std::endl;
	std::vector<std::string> fileNames, tmpStrVec, OBJfileNames;
	std::vector<Eigen::MatrixXd> meshesVers, outVers;
	std::vector<Eigen::MatrixXi> meshesTris, outTris;

	getFileNames(pathOBJ.c_str(), tmpStrVec, false);
	const unsigned meshesCount = tmpStrVec.size();
	OBJfileNames.reserve(meshesCount);
	for (const auto& str : tmpStrVec)
	{
		std::string tailStr = str.substr(str.size() - 4, 4);				//	".obj"
		if (".obj" == tailStr)
		{
			fileNames.push_back(str);
			unsigned index = str.find_last_of("/");
			std::string OBJfileName = str.substr(index, str.size() - index - 4);			// "/" + obj�ļ���������·����.obj��׺��
			OBJfileNames.push_back(OBJfileName);
		}
	}
	meshesVers.resize(meshesCount);
	meshesTris.resize(meshesCount);
	outVers.resize(meshesCount);
	outTris.resize(meshesCount);
	for (unsigned i = 0; i < meshesCount; ++i)
		objReadMeshMat(meshesVers[i], meshesTris[i], fileNames[i].c_str());
	tt.endCout("��ȡ���������ʱ��");

	// 2. �������ִ��laplace��˳��
	tt.start();
	debugDisp("ִ��laplace��˳��");
	float deltaLB = INIGetFloat(TEXT("deltaLB"), fileConfig);
	unsigned loopCount = INIGetInt(TEXT("laplaceFaringLoopCount"), fileConfig);
	for (int i = 0; i < meshesCount; ++i)
		laplaceFaring(outVers[i], meshesVers[i], meshesTris[i], deltaLB, loopCount);
	outTris = std::move(meshesTris);
	tt.endCout("laplace��˳��");

	// 3. �����
	for (unsigned i = 0; i < meshesCount; ++i)
	{
		std::string str = pathOutput + OBJfileNames[i] + std::string{ ".obj" };
		objWriteMeshMat(str.c_str(), outVers[i], outTris[i]);
	}

	debugDisp("finished.");

	return 0;
}


int testCmd_hausdorffDistance(int argc, char** argv)
{
	// ����·����
	int   nPos;
	CString   cPath;
	GetModuleFileName(NULL, cPath.GetBufferSetLength(MAX_PATH + 1), MAX_PATH);		// ��ȡ��ǰ���̼��ص�ģ���·����
	nPos = cPath.ReverseFind('\\');
	cPath = cPath.Left(nPos);
	std::string path{ CT2CA{cPath} };
	std::string pathOBJ = path + "\\inputOBJ";
	std::string pathOutput = path + "\\outputData";
	CString fileConfig = cPath + "\\config.ini";

	// 1. ��ȡ��������  
	std::cout << "��ȡ��������..." << std::endl;
	std::vector<std::string> fileNames, tmpStrVec, OBJfileNames;
	std::vector<Eigen::MatrixXd> meshesVers, outVers;
	std::vector<Eigen::MatrixXi> meshesTris, outTris;

	getFileNames(pathOBJ.c_str(), tmpStrVec, false);
	const unsigned meshesCount = tmpStrVec.size();
	if (2 != meshesCount)
	{
		debugDisp("error!!! ��������������Ϊ2");
		return -1;
	}

	OBJfileNames.reserve(meshesCount);
	for (const auto& str : tmpStrVec)
	{
		std::string tailStr = str.substr(str.size() - 4, 4);				//	".obj"
		if (".obj" == tailStr)
		{
			fileNames.push_back(str);
			unsigned index = str.find_last_of("/");
			std::string OBJfileName = str.substr(index, str.size() - index - 4);			// "/" + obj�ļ���������·����.obj��׺��
			OBJfileNames.push_back(OBJfileName);
		}
	}
	meshesVers.resize(meshesCount);
	meshesTris.resize(meshesCount);
	outVers.resize(meshesCount);
	outTris.resize(meshesCount);
	for (unsigned i = 0; i < meshesCount; ++i)
		objReadMeshMat(meshesVers[i], meshesTris[i], fileNames[i].c_str());

	// 2. �������������hausdorff���룺
	double hd = 0;
	igl::hausdorff(meshesVers[0], meshesTris[0], meshesVers[1], meshesTris[1], hd);
	debugDisp("hausdorff distance == ", hd);

	debugDisp("finished.");
	getchar();

	return 0;
}



int main(int argc, char** argv)
{
	// IGL_DEFORMATION::test0();

	// IGL_MODELLING::test1();

	IGL_BASIC::test1();

	std::cout << "main() finished." << std::endl;
}

