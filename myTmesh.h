#pragma once

#include "myEigen.h"
#include "tmesh.h"								// ����������TMESH
#include "detectIntersections.h"			// TMESHS�����ཻ��⹦�ܣ�


static std::string g_debugPath = "E:/";



////////////////////////////////////////////////////////////////////////////////////////////////// ǰ������
template <typename F>
void traverseVersList(const T_MESH::List& list, F f);
template <typename F>
void traverseVersList(T_MESH::List& list, F f);
template <typename F>
void traverseEdgesList(const T_MESH::List& list, F f);
template <typename F>
void traverseEdgesList(T_MESH::List& list, F f);
template <typename F>
void traverseTrisList(const T_MESH::List& list, F f);
template <typename F>
void traverseTrisList(T_MESH::List& list, F f);
template <typename T>
void TMesh2MeshMat(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, T_MESH::Basic_TMesh& mesh);
template <typename T>
void meshMat2tMesh(T_MESH::Basic_TMesh& mesh, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const Eigen::MatrixXi& tris);
template <typename T>
double orient3D(const Eigen::Matrix<T, 1, 3>& v1, const Eigen::Matrix<T, 1, 3>& v2, const Eigen::Matrix<T, 1, 3>& v3, \
	const Eigen::Matrix<T, 1, 3>& v4);
template <typename T>
double orient2D(const Eigen::Matrix<T, 1, 2>& v1, const Eigen::Matrix<T, 1, 2>& v2, const Eigen::Matrix<T, 1, 2>& v3);
void genAABBmesh(const T_MESH::di_cell& cell, Eigen::MatrixXd& vers, Eigen::MatrixXi& tris);


/// /////////////////////////////////////////////////////////////////////////////////////////// DEBUG �ӿ�

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

template<typename T>
static void debugWriteMesh(const char* name, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteMeshMat(path, vers, tris);
}


static void debugWriteMesh(const char* name, T_MESH::Basic_TMesh& mesh)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);

	// tmesh�Լ���IO�ӿھ������ޣ�ת��ΪEigen�����ʾȻ�������
	Eigen::MatrixXd vers;
	Eigen::MatrixXi tris;
	TMesh2MeshMat(vers, tris, mesh);
	objWriteMeshMat(path, vers, tris);
}


template<typename DerivedV>
static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers)
{
	char path[512] = { 0 };
	sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
	objWriteEdgesMat(path, edges, vers);
}





//////////////////////////////////////////////////////////////////////////////////////////////////TMESH��ؽӿڣ�

// ��TMESH�����ֻ��������
template <typename F>
void traverseVersList(const T_MESH::List& list, F f)
{
	const T_MESH::Vertex* verPtr = nullptr;					// �ײ�constָ��
	const T_MESH::Node* nodePtr = nullptr;
	const T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), verPtr = nodePtr ? ((const T_MESH::Vertex*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), verPtr = (nodePtr) ? ((const T_MESH::Vertex*)nodePtr->data) : nullptr)
	{
		f(verPtr);
	}
}


// ��TMESH����ķ�ֻ��������
template <typename F>
void traverseVersList(T_MESH::List& list, F f)
{
	T_MESH::Vertex* verPtr = nullptr;
	T_MESH::Node* nodePtr = nullptr;
	T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), verPtr = nodePtr ? ((T_MESH::Vertex*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), verPtr = (nodePtr) ? ((T_MESH::Vertex*)nodePtr->data) : nullptr)
	{
		f(verPtr);
	}
}


template <typename F>
void traverseEdgesList(const T_MESH::List& list, F f)
{
	const T_MESH::Edge* ePtr = nullptr;
	const T_MESH::Node* nodePtr = nullptr;
	const T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), ePtr = nodePtr ? ((const T_MESH::Edge*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), ePtr = (nodePtr) ? ((const T_MESH::Edge*)nodePtr->data) : nullptr)
	{
		f(ePtr);
	}
}

template <typename F>
void traverseEdgesList(T_MESH::List& list, F f)
{
	T_MESH::Edge* ePtr = nullptr;
	T_MESH::Node* nodePtr = nullptr;
	T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), ePtr = nodePtr ? ((T_MESH::Edge*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), ePtr = (nodePtr) ? ((T_MESH::Edge*)nodePtr->data) : nullptr)
	{
		f(ePtr);
	}
}


template <typename F>
void traverseTrisList(const T_MESH::List& list, F f)
{
	const T_MESH::Triangle* triPtr = nullptr;
	const T_MESH::Node* nodePtr = nullptr;
	const T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), triPtr = nodePtr ? ((const T_MESH::Triangle*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), triPtr = (nodePtr) ? ((const T_MESH::Triangle*)nodePtr->data) : nullptr)
	{
		f(triPtr);
	}
}


template <typename F>
void traverseTrisList(T_MESH::List& list, F f)
{
	T_MESH::Triangle* triPtr = nullptr;
	T_MESH::Node* nodePtr = nullptr;
	T_MESH::List* listPtr = &list;
	for (nodePtr = listPtr->head(), triPtr = nodePtr ? ((T_MESH::Triangle*)nodePtr->data) : nullptr; \
		nodePtr != nullptr; nodePtr = nodePtr->next(), triPtr = (nodePtr) ? ((T_MESH::Triangle*)nodePtr->data) : nullptr)
	{
		f(triPtr);
	}
}


// TMESH����ת��Ϊ���� 
template <typename T>
void TMesh2MeshMat(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, Eigen::MatrixXi& tris, T_MESH::Basic_TMesh& mesh)
{
	const int versCount = mesh.V.numels();
	const int trisCount = mesh.T.numels();
	int index = 0;
	vers.resize(versCount, 3);
	tris.resize(trisCount, 3);

	// 1. ���ɶ������
	traverseVersList(mesh.V, [&](T_MESH::Vertex* verPtr)
		{
			vers(index, 0) = verPtr->x;
			vers(index, 1) = verPtr->y;
			vers(index, 2) = verPtr->z;
			index++;
		});

	// 2. ����ڵ��x�����ݴ���xValues�У�Ȼ�����дΪ����������
	std::vector<double> xValues;
	xValues.reserve(versCount);

	index = 0;
	traverseVersList(mesh.V, [&](T_MESH::Vertex* verPtr)
		{
			xValues.push_back(verPtr->x);
			verPtr->x = static_cast<double>(index++);
		});


	// 3. ��������Ƭ���ݣ��������
	index = 0;
	traverseTrisList(mesh.T, [&](T_MESH::Triangle* triPtr)
		{
			tris(index, 0) = static_cast<int>(triPtr->v1()->x);
			tris(index, 1) = static_cast<int>(triPtr->v2()->x);
			tris(index, 2) = static_cast<int>(triPtr->v3()->x);
			index++;
		});

	// 4. ����ڵ��е�x���ݻ�ԭ��
	index = 0;
	traverseVersList(mesh.V, [&](T_MESH::Vertex* verPtr)
		{
			verPtr->x = xValues[index++];
		});
}


template <typename T>
void meshMat2tMesh(T_MESH::Basic_TMesh& mesh, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, \
	const Eigen::MatrixXi& tris)
{
	const unsigned versCount = vers.rows();
	const unsigned trisCount = tris.rows();

	for (unsigned i = 0; i < versCount; ++i)
		mesh.V.appendTail(mesh.newVertex(vers(i, 0), vers(i, 1), vers(i, 2)));

	T_MESH::Vertex* vPtr = nullptr;
	T_MESH::Node* nodePtr = nullptr;
	T_MESH::ExtVertex** var = (T_MESH::ExtVertex**)malloc(sizeof(T_MESH::ExtVertex*) * versCount);

	unsigned index = 0;
	traverseVersList(mesh.V, [&](T_MESH::Vertex* vPtr)
		{
			var[index++] = new T_MESH::ExtVertex(vPtr);
		});

	for (unsigned i = 0; i < trisCount; ++i)
		mesh.CreateIndexedTriangle(var, tris(i, 0), tris(i, 1), tris(i, 2));

	mesh.fixConnectivity();
}


// �����ĸ�3D��Χ�ɵ�������ķ��������6������v2,v3,v4Χ�ɵ�������Ϊ����������ķ��������������ֵ��v1�����������ࣻ��ֵ��v1�������θ��ࣻ0���ĵ㹲�棻
/*
	Returns a positive value if the point d lies above the plane passing through a, b, and c, meaning that a, b,
			and c appear in counterclockwise order when viewed from d.

	Returns a negative value if d lies below the plane.

	Returns zero if the points are coplanar.

	The result is also an approximation of six times the signed volume of the tetrahedron defined by the four points.

*/
template <typename T>
double orient3D(const Eigen::Matrix<T, 1, 3>& v1, const Eigen::Matrix<T, 1, 3>& v2, const Eigen::Matrix<T, 1, 3>& v3, const Eigen::Matrix<T, 1, 3>& v4)
{
	std::vector<double> p1{ v1(0), v1(1), v1(2) };
	std::vector<double> p2{ v2(0), v2(1), v2(2) };
	std::vector<double> p3{ v3(0), v3(1), v3(2) };
	std::vector<double> p4{ v4(0), v4(1), v4(2) };
	return orient3d(&p1[0], &p2[0], &p3[0], &p4[0]);
}


// ��������2D��Χ�������εķ��������2�����������λ�ù�ϵ��������: ��ʱ�룻 ������˳ʱ�룻 0�����ߣ� 
/*
	Returns a positive value if the points a, b, and c occur in counterclockwise order
			(c lies to the left of the directed line defined by points a and b).

	Returns a negative value if they occur in clockwise order (c lies to the right of the directed line ab).

	Returns zero if they are collinear.

	The result is also an approximation of twice the signed area of the triangle defined by the three points.
*/
template <typename T>
double orient2D(const Eigen::Matrix<T, 1, 2>& v1, const Eigen::Matrix<T, 1, 2>& v2, const Eigen::Matrix<T, 1, 2>& v3)
{
	std::vector<double> p1{ v1(0), v1(1) };
	std::vector<double> p2{ v2(0), v2(1) };
	std::vector<double> p3{ v3(0), v3(1) };
	return orient2d(&p1[0], &p2[0], &p3[0]);
}


// ����cell��Ӧ�������Χ������
void genAABBmesh(const T_MESH::di_cell& cell, Eigen::MatrixXd& vers, Eigen::MatrixXi& tris);



// ��������������tmesh
namespace TEST_TMESH
{
	void test0();
	void test1();
	void test2();
	void test3();
	void test4();
	void test44();
	void test5();
	void test6();
	void test7();
	void test8();
	void test9();
	void test10();
	void test11();
	void test12();
	void test13();
	void test14();
	void test15();
	void test16();
}