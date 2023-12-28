#include "test_imgui.h"

#include <iostream>
#include <vector>
#include <algorithm>

#include "glad/glad.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imgui_internal.h"
#include "GLFW/glfw3.h"

#include "triMesh.h"
#include "representations.h"
#include "TEST_IMGUI/MathUtil.h"
#include "TEST_IMGUI/fitcurve.h"
#include "TEST_IMGUI/myapp.h" 

#include "myEigenIO/myEigenIO.h"
#pragma comment(lib,"myEigenIO.lib")	


////////////////////////////////////////////////////////////////////////////////////////////// DEBUG 接口
namespace MY_DEBUG
{
    static std::string g_debugPath = "E:/";

    // lambda——打印std::cout支持的类型变量。
    template <typename T>
    static auto disp = [](const T& arg)
    {
        std::cout << arg << ", ";
    };

    static void debugDisp()			// 递归终止
    {						//		递归终止设为无参或者一个参数的情形都可以。
        std::cout << std::endl;
        return;
    }


    template <typename T, typename... Types>
    static void debugDisp(const T& firstArg, const Types&... args)
    {
        std::cout << firstArg << " ";
        debugDisp(args...);
    }


    template <typename Derived>
    static void dispData(const Eigen::MatrixBase<Derived>& m)
    {
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
    static void debugWriteVers(const char* name, const Eigen::MatrixBase<DerivedV>& vers)
    {
        char path[512] = { 0 };
        sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
        objWriteVerticesMat(path, vers);
    }


    template<typename DerivedV>
    static void debugWriteVers2D(const char* name, const Eigen::MatrixBase<DerivedV>& vers)
    {
        char path[512] = { 0 };
        sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
        objWriteVerticesMat2D(path, vers);
    }


    template<typename DerivedV>
    static void debugWriteMesh(const char* name, \
        const Eigen::MatrixBase<DerivedV>& vers, const Eigen::MatrixXi& tris)
    {
        char path[512] = { 0 };
        sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
        objWriteMeshMat(path, vers, tris);
    }


    template<typename DerivedV>
    static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, \
        const Eigen::MatrixBase<DerivedV>& vers)
    {
        char path[512] = { 0 };
        sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
        objWriteEdgesMat(path, edges, vers);
    }

}
using namespace MY_DEBUG;


// 三次B样条插值
namespace BSPLINE_INTERP
{
    using namespace TRIANGLE_MESH;


    std::vector<verF> SmoothLoop(const std::vector<verF>& vertices, float weight)
    {
        auto versCount = vertices.size();
        std::vector<verF> result(versCount, verF{ 0, 0, 0 });
        Eigen::MatrixXf circuit(versCount, 3);
        for (size_t i = 0; i < versCount; i++)
            circuit.row(i) << vertices[i].x, vertices[i].y, vertices[i].z;
        Eigen::SparseMatrix<float> W(versCount, versCount);
        std::vector<Eigen::Triplet<float>> data;
        data.reserve(2 * versCount);
        for (int i = 0; i < versCount; ++i)
        {
            if (i == 0)
            {
                data.push_back(Eigen::Triplet<float>(i, versCount - 1, 1));
                data.push_back(Eigen::Triplet<float>(i, i + 1, 1));
            }
            else if (i == versCount - 1)
            {
                data.push_back(Eigen::Triplet<float>(i, i - 1, 1));
                data.push_back(Eigen::Triplet<float>(i, 0, 1));
            }
            else
            {
                data.push_back(Eigen::Triplet<float>(i, i - 1, 1));
                data.push_back(Eigen::Triplet<float>(i, i + 1, 1));
            }
        }
        W.setFromTriplets(data.begin(), data.end());
        Eigen::SparseMatrix<float> A1(versCount, versCount);
        A1.setIdentity();
        A1 -= 0.5 * W;
        Eigen::SparseMatrix<float> B1(versCount, 3);
        B1.setZero();
        Eigen::SparseMatrix<float> A2(versCount, versCount);
        A2.setIdentity();
        Eigen::SparseMatrix<float> B2 = circuit.sparseView();
        Eigen::SparseMatrix<float> tempMat(versCount, 2 * versCount);
        tempMat.leftCols(versCount) = A1.transpose();
        tempMat.rightCols(versCount) = weight * A2.transpose();
        Eigen::SparseMatrix<float> A = tempMat.transpose();
        tempMat.resize(3, 2 * versCount);
        tempMat.leftCols(versCount) = B1.transpose();
        tempMat.rightCols(versCount) = weight * B2.transpose();
        Eigen::SparseMatrix<float> B = tempMat.transpose();
        A.makeCompressed();
        B.makeCompressed();

        // 最小二乘法解超定线性方程组A*X = B;
        Eigen::MatrixXf ATA = (A.transpose() * A).toDense();
        Eigen::MatrixXf ATB = (A.transpose() * B).toDense();
        circuit = ATA.inverse() * ATB;

        for (size_t i = 0; i < versCount; i++)
        {
            result[i].x = circuit(i, 0);
            result[i].y = circuit(i, 1);
            result[i].z = circuit(i, 2);
        }
        return result;
    }


    void PointTriangleDistance(const std::vector<verF>& triangle, const verF& point, float& distance, verF& PP0)
    {
        constexpr float zero = 0;
        constexpr float one = 1;
        constexpr float two = 2;
        auto diff = triangle[0] - point;
        auto edge0 = triangle[1] - triangle[0];
        auto edge1 = triangle[2] - triangle[0];
        auto a00 = edge0.dot(edge0);
        auto a01 = edge0.dot(edge1);
        auto a11 = edge1.dot(edge1);
        auto b0 = diff.dot(edge0);
        auto b1 = diff.dot(edge1);
        auto det = (a00 * a11 - a01 * a01) > zero ? (a00 * a11 - a01 * a01) : zero;
        auto s = a01 * b1 - a11 * b0;
        auto t = a01 * b0 - a00 * b1;

        if (s + t <= det)
        {
            if (s < zero)
            {
                if (t < zero)  // region 4
                {
                    if (b0 < zero)
                    {
                        t = zero;
                        s = -b0 >= a00 ? one : -b0 / a00;
                    }
                    else
                    {
                        s = zero;
                        if (b1 >= zero)
                            t = zero;
                        else if (-b1 >= a11)
                            t = one;
                        else
                            t = -b1 / a11;
                    }
                }
                else  // region 3
                {
                    s = zero;
                    if (b1 >= zero)
                        t = zero;
                    else if (-b1 >= a11)
                        t = one;
                    else
                        t = -b1 / a11;
                }
            }
            else if (t < zero)  // region 5
            {
                t = zero;
                if (b0 >= zero)
                    s = zero;
                else if (-b0 >= a00)
                    s = one;
                else
                    s = -b0 / a00;
            }
            else  // region 0
            {
                // minimum at interior point
                s /= det;
                t /= det;
            }
        }
        else
        {
            float tmp0{}, tmp1{}, numer{}, denom{};

            if (s < zero)  // region 2
            {
                tmp0 = a01 + b0;
                tmp1 = a11 + b1;
                if (tmp1 > tmp0)
                {
                    numer = tmp1 - tmp0;
                    denom = a00 - two * a01 + a11;
                    if (numer >= denom)
                    {
                        s = one;
                        t = zero;
                    }
                    else
                    {
                        s = numer / denom;
                        t = one - s;
                    }
                }
                else
                {
                    s = zero;
                    if (tmp1 <= zero)
                        t = one;
                    else if (b1 >= zero)
                        t = zero;
                    else
                        t = -b1 / a11;
                }
            }
            else if (t < zero)  // region 6
            {
                tmp0 = a01 + b1;
                tmp1 = a00 + b0;
                if (tmp1 > tmp0)
                {
                    numer = tmp1 - tmp0;
                    denom = a00 - two * a01 + a11;
                    if (numer >= denom)
                    {
                        t = one;
                        s = zero;
                    }
                    else
                    {
                        t = numer / denom;
                        s = one - t;
                    }
                }
                else
                {
                    t = zero;
                    if (tmp1 <= zero)
                        s = one;
                    else if (b0 >= zero)
                        s = zero;
                    else
                        s = -b0 / a00;
                }
            }
            else  // region 1
            {
                numer = a11 + b1 - a01 - b0;
                if (numer <= zero)
                {
                    s = zero;
                    t = one;
                }
                else
                {
                    denom = a00 - two * a01 + a11;
                    if (numer >= denom)
                    {
                        s = one;
                        t = zero;
                    }
                    else
                    {
                        s = numer / denom;
                        t = one - s;
                    }
                }
            }
        }

        PP0 = triangle.at(0) + s * edge0 + t * edge1;
        distance = PP0.distance(point);
    }


    void SmoothFeaturePoints(std::vector<std::vector<std::vector<verF>>>& featurePoints)
    {
        int meshCnt = featurePoints[0].size();
        int totalStep = featurePoints.size();
        for (int i = 0; i < meshCnt; ++i)
        {
            for (int j = 0; j < 8; ++j)
            {
                std::vector<verF> V;
                for (int k = 0; k < totalStep; ++k)
                {
                    auto Fp = featurePoints.at(k).at(i);
                    V.push_back(Fp.at(j));
                }
                auto smoothV = SmoothLoop(V, 0.1);
                for (int i = 0; i < 3; ++i)
                {
                    smoothV.at(i) = V.at(i);
                    smoothV.at(totalStep - i - 1) = V.at(totalStep - i - 1);
                }
                for (int k = 0; k < totalStep; ++k)
                {
                    V = featurePoints.at(k).at(i);
                    V.at(j) = smoothV.at(k);
                    featurePoints.at(k).at(i) = V;
                }
            }
        }
    }


    // to be optimized...搜索点云中距离目标顶点targetVer最近的那个顶点；当前使用暴力算法，可以进一步提升性能；
    template <typename T1, typename T2>
    size_t nearestVerIdx(const std::vector<triplet<T1>>& vers, const triplet<T2>& targetVer)
    {
        Eigen::MatrixXf versMat;
        Eigen::RowVector3f tarVerVec;
        vers2mat(versMat, vers);
        ver2vec(tarVerVec, targetVer);
        Eigen::VectorXf dises = (versMat.rowwise() - tarVerVec).rowwise().norm();
        Eigen::VectorXf::Index minIdx;
        dises.minCoeff(&minIdx);
        return static_cast<size_t>(minIdx);
    }


    void ClosestPoint2Mesh(const triMeshF mesh, const std::vector<verF>& qPoint, std::vector<verF>& P)
    {
        std::pair<unsigned, float> pair_;
        std::vector<size_t> indices;
        std::vector<verF> CandidateP;
        for (int j = 0; j < qPoint.size(); ++j)
        {
            size_t idx = nearestVerIdx(mesh.vertices, qPoint.at(j));			// mesh.vertice[idx]是网格点云中距离qPoint.at(j)最近的点；
            indices.push_back(idx);
            CandidateP.push_back(mesh.vertices.at(idx));
        }

        for (int i = 0; i < indices.size(); ++i)
        {
            std::vector<unsigned> connectedFaces;
            for (int j = 0; j < mesh.triangles.size(); ++j)
                if (mesh.triangles.at(j).x == indices.at(i) ||
                    mesh.triangles.at(j).y == indices.at(i) ||
                    mesh.triangles.at(j).z == indices.at(i))
                    connectedFaces.push_back(j);
            float minDist = 0xFFFF;
            for (const auto& f : connectedFaces)
            {
                std::vector<verF> TRI(3);
                TRI = { mesh.vertices.at(mesh.triangles.at(f).x),
                    mesh.vertices.at(mesh.triangles.at(f).y),
                    mesh.vertices.at(mesh.triangles.at(f).z) };
                float dist;
                verF PP0;
                PointTriangleDistance(TRI, qPoint.at(i), dist, PP0);
                if (dist < minDist)
                {
                    minDist = dist;
                    CandidateP.at(i) = PP0;
                }
            }
            if (i == 0)
                P.push_back(CandidateP.at(i));
            else
                if (CandidateP.at(i).distance(CandidateP.back()) > 0.1)
                    P.push_back(CandidateP.at(i));
        }

    }


    std::vector<verF> TanVecDir(const std::vector<verF>& sampleVers)
    {
        auto sampleVersCount = sampleVers.size();
        auto qVecNum = sampleVersCount + 3;
        std::vector<verF> qVec(qVecNum, verF{ 0, 0, 0 }), tanVec(sampleVersCount, verF{ 0, 0, 0 });

        for (int i = 2; i < qVecNum - 2; ++i)
            qVec.at(i) = sampleVers.at(i - 1) - sampleVers.at(i - 2);

        if (sampleVers[0].distance(sampleVers[sampleVersCount - 1]) <= FLT_MIN)
        {
            qVec.at(1) = sampleVers.at(sampleVersCount - 1) - sampleVers.at(sampleVersCount - 2);
            qVec.at(0) = sampleVers.at(sampleVersCount - 2) - sampleVers.at(sampleVersCount - 3);
            qVec.at(qVecNum - 2) = sampleVers.at(1) - sampleVers.at(0);
            qVec.at(qVecNum - 1) = sampleVers.at(2) - sampleVers.at(1);
        }
        else
        {
            qVec.at(1) = 2 * qVec.at(2) - qVec.at(3);
            qVec.at(0) = 2 * qVec.at(1) - qVec.at(2);
            qVec.at(qVecNum - 2) = 2 * qVec.at(qVecNum - 3) - qVec.at(qVecNum - 4);
            qVec.at(qVecNum - 1) = 2 * qVec.at(qVecNum - 2) - qVec.at(qVecNum - 3);
        }
        for (int i = 0; i < sampleVersCount; ++i)
        {
            float interppar = 0.5;
            auto numerL = qVec.at(i).cross(qVec.at(i + 1)).length();
            auto numerR = qVec.at(i + 2).cross(qVec.at(i + 3)).length();
            if (auto numer = numerL + numerR; numer != 0)
                interppar = numerL / numer;
            tanVec.at(i) = (1 - interppar) * qVec.at(i + 1) + interppar * qVec.at(i + 2);
            if (auto normTanVec = tanVec.at(i).length(); normTanVec == 0)
                if (i == 0)
                    tanVec.at(i) = qVec.at(i + 2) / (qVec.at(i + 2).length());
                else
                    tanVec.at(i) = qVec.at(i + 1) / (qVec.at(i + 1).length());
            else
                tanVec.at(i) /= normTanVec;
        }

        return tanVec;
    }


    std::vector<verF> ContlPoint(const std::vector<verF>& sampleVers, const std::vector<verF>& tanVecDir)
    {
        auto sampleVersCount = sampleVers.size();
        auto curNum = sampleVersCount - 1;
        std::vector<verF> contlPt(2 * curNum + 2, verF{ 0, 0, 0 });

        contlPt.at(0) = sampleVers.at(0);
        contlPt.at(contlPt.size() - 1) = sampleVers.at(curNum);

        std::vector<verF> preTanVD(tanVecDir.begin(), tanVecDir.begin() + curNum);
        std::vector<verF> postTanVD(tanVecDir.begin() + 1, tanVecDir.end());
        std::vector<verF> TanVD(curNum), DataPt(curNum);
        for (int i = 0; i < curNum; ++i)
            TanVD.at(i) = preTanVD.at(i) + postTanVD.at(i);
        std::vector<verF> preDataPt(sampleVers.begin(), sampleVers.begin() + curNum);
        std::vector<verF> postDataPt(sampleVers.begin() + 1, sampleVers.end());
        for (int i = 0; i < curNum; ++i)
            DataPt.at(i) = postDataPt.at(i) - preDataPt.at(i);

        for (int i = 0; i < curNum; ++i)
        {
            auto a = 16 - pow(TanVD.at(i).length(), 2);
            auto b = 12 * DataPt.at(i).dot(TanVD.at(i));
            auto c = -36 * pow(DataPt.at(i).length(), 2);
            auto normTanV = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
            contlPt.at(2 * i + 1) = sampleVers.at(i) + normTanV * tanVecDir.at(i) / 3;
            contlPt.at(2 * i + 2) = sampleVers.at(i + 1) - normTanV * tanVecDir.at(i + 1) / 3;
        }

        return contlPt;
    }


    std::vector<float> NodeVec(const std::vector<verF>& sampleVers, const std::vector<verF>& contlPoint)
    {
        int sampleVersCount = sampleVers.size();
        std::vector parDataPt(sampleVersCount, 0.f), nodeV(2 * sampleVersCount, 0.f);

        for (int _ = 0; _ < 4; ++_)
            nodeV.push_back(1.f);
        for (int i = 0; i < sampleVersCount - 1; ++i)
            parDataPt.at(i + 1) = parDataPt.at(i) + 3 * contlPoint.at(2 * i + 1).distance(sampleVers.at(i));
        auto endParDP = parDataPt.at(sampleVersCount - 1);
        for (int i = 1; i < sampleVersCount - 1; ++i)
        {
            nodeV.at(2 * i + 2) = parDataPt.at(i) / endParDP;
            nodeV.at(2 * i + 3) = parDataPt.at(i) / endParDP;
        }

        return nodeV;
    }


    std::vector<float> BasisFuns(int span, float u, int p, const std::vector<float>& nodeVec)
    {
        std::vector N(p + 2, 0.f);
        N[1] = 1;
        for (int count = 1; count <= p; ++count)
        {
            std::vector preN = N;
            N.assign(p + 2, 0.f);
            for (int reg = 1; reg <= count + 1; ++reg)
            {
                auto denLeft = nodeVec.at(span + reg - 1) - nodeVec.at(span + reg - 1 - count);
                auto denRight = nodeVec.at(span + reg) - nodeVec.at(span + reg - count);
                auto numerL = u - nodeVec.at(span + reg - 1 - count);
                auto numerR = nodeVec.at(span + reg) - u;
                float itemL = 0, itemR = 0;
                if (denLeft != 0 && reg != 1)
                    itemL = numerL / denLeft * preN.at(reg - 1);
                if (denRight != 0 && reg != count + 1)
                    itemR = numerR / denRight * preN.at(reg);
                N.at(reg) = itemL + itemR;
            }
        }
        return N;
    }


    int FindSpan(int p, float u, const std::vector<float>& nodeVec)
    {
        auto nodeNum = nodeVec.size();
        if (u == nodeVec.at(nodeNum - p))
            return nodeNum - p - 1;
        else
        {
            auto low = p + 1;
            auto high = nodeNum - p;
            auto mid = floor((low + high) / 2);
            while (u < nodeVec.at(mid) || u >= nodeVec.at(mid + 1))
            {
                if (u < nodeVec.at(mid))
                    high = mid;
                else
                    low = mid;
                mid = floor((low + high) / 2);
            }
            return mid;
        }
    }


    verF SinCurPoint(int p, float u, const std::vector<float>& nodeVec, const std::vector<verF>& contlPoint)
    {
        std::vector nodeVec_ = { 0.f };
        nodeVec_.insert(nodeVec_.end(), nodeVec.begin(), nodeVec.end());
        std::vector<verF> contlPoint_ = { verF{0, 0, 0} };
        contlPoint_.insert(contlPoint_.end(), contlPoint.begin(), contlPoint.end());
        int span = FindSpan(p, u, nodeVec_);
        std::vector<float> func = BasisFuns(span, u, p, nodeVec_);
        verF result = verF{ 0, 0, 0 };
        for (int i = 1; i <= p + 1; ++i)
            result += func[i] * contlPoint_[span + i - 1 - p];
        return result;
    }
}


// 三次B样条插值拟合回路曲线	
template<typename DerivedVo, typename DerivedVi>
bool cubicBSplineInterpLoop3D(Eigen::PlainObjectBase<DerivedVo>& versOut, \
    const Eigen::PlainObjectBase<DerivedVi>& versIn, const unsigned versOutCount)
{
    using namespace BSPLINE_INTERP;
    using ScalarO = typename DerivedVo::Scalar;
    versOut.resize(0, 0);
    if (0 == versIn.rows())
        return true;
    if (versIn.rows() < 3)
    {
        versOut = versIn.cast<ScalarO>();
        return true;
    }

    std::vector<verF> versLoop = mat2versF(versIn);
    versLoop.push_back(versLoop[0]);

    std::vector<verF> tanVecDir = TanVecDir(versLoop);
    std::vector<verF> cntrlPts = ContlPoint(versLoop, tanVecDir);
    std::vector<float> nodeVector = NodeVec(versLoop, cntrlPts);

    int n = versOutCount + 1;
    std::vector<float> u;
    for (int i = 0; i < n; ++i)
        u.push_back(static_cast<float>(i) / (n - 1));

    std::vector<verF> curveFitted(versOutCount, verF{ 0, 0, 0 });
    for (int count = 0; count < versOutCount; ++count)
        curveFitted[count] = SinCurPoint(3, u.at(count), nodeVector, cntrlPts);

    vers2mat(versOut, curveFitted);

    return true;
}


// 三次B样条插值拟合非回路曲线	
template<typename DerivedVo, typename DerivedVi>
bool cubicBSplineInterpCurve3D(Eigen::PlainObjectBase<DerivedVo>& versOut, \
    const Eigen::PlainObjectBase<DerivedVi>& versIn, const unsigned versOutCount)
{
    using namespace BSPLINE_INTERP;
    using ScalarO = typename DerivedVo::Scalar;
    versOut.resize(0, 0);

    if (0 == versIn.rows())
        return true;
    if (1 == versIn.rows())
    {
        versOut = versIn.cast<ScalarO>();
        return true;
    }

    std::vector<verF> versSample = mat2versF(versIn); 
    std::vector<verF> tanVecDir = TanVecDir(versSample);
    std::vector<verF> cntrlPts = ContlPoint(versSample, tanVecDir);
    std::vector<float> nodeVector = NodeVec(versSample, cntrlPts);

    int n = versOutCount + 1;
    std::vector<float> u;
    for (int i = 0; i < n; ++i)
        u.push_back(static_cast<float>(i) / (n - 1));

    std::vector<verF> curveFitted(versOutCount, verF{ 0, 0, 0 });
    for (int count = 0; count < versOutCount; ++count)
        curveFitted[count] = SinCurPoint(3, u.at(count), nodeVector, cntrlPts);

    vers2mat(versOut, curveFitted);

    return true;
}


namespace MY_IMGUI 
{
    GLFWwindow* g_window = nullptr;                         // 全局的窗口对象
    int g_width = 1600;
    int g_height = 1200;

    ImVec2 g_canvas_pos_ul = { 0.0f, 0.0f };
    ImVec2 g_canvas_pos_br = { 0.0f, 0.0f };


    // 输入2D曲线点集，在画布上画出曲线
    void PlotLineSegments(const std::vector<Eigen::Vector2f>& poss, ImDrawList* draw_list, \
        ImU32 line_col, ImU32 point_col)
    {
        for (size_t i = 1; i < poss.size(); i++)
        {
            draw_list->AddLine({ g_canvas_pos_ul.x + poss[i - 1].x(), g_canvas_pos_br.y - poss[i - 1].y() },
                { g_canvas_pos_ul.x + poss[i].x(), g_canvas_pos_br.y - poss[i].y() }, line_col, 2.0f);
        }
        for (const auto& pos : poss)
        {
            draw_list->AddCircleFilled({ g_canvas_pos_ul.x + pos.x(), g_canvas_pos_br.y - pos.y() }, 5.0f, point_col);
            draw_list->AddCircle({ g_canvas_pos_ul.x + pos.x(), g_canvas_pos_br.y - pos.y() }, 5.0f, point_col);
        }
    }


    // 输入2D曲线点集，在画布上画出曲线
    void PlotLineSegments(const std::vector<float>& pos_x, const std::vector<float>& pos_y, \
        ImDrawList* draw_list, ImU32 line_col, ImU32 point_col, ImVec2 ul = g_canvas_pos_ul,\
        ImVec2 br = g_canvas_pos_br) 
    {
        const size_t n = pos_x.size();
        for (size_t i = 1; i < n; i++) 
        {
            draw_list->AddLine({ ul.x + pos_x[i - 1], br.y - pos_y[i - 1] },
                { ul.x + pos_x[i], br.y - pos_y[i] }, line_col, 2.0f);
        }
        for (size_t i = 0; i < n; i++) 
        {
            draw_list->AddCircleFilled({ ul.x + pos_x[i], br.y - pos_y[i] }, 5.0f, point_col);
            draw_list->AddCircle({ ul.x + pos_x[i], br.y - pos_y[i] }, 5.0f, point_col);
        }
    }


    bool Initialize()
    {
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
        g_window = glfwCreateWindow(g_width, g_height, "GAMES102 hw1", nullptr, nullptr);

        if (g_window == nullptr)
        {
            std::cerr << "Failed to create GLFW window" << std::endl;
            glfwTerminate();
            return false;
        }

        glfwMakeContextCurrent(g_window);
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
        {
            std::cout << "Failed to load glad" << std::endl;
            glfwTerminate();
            return false;
        }

        glViewport(0, 0, g_width, g_height);
        std::cout << "GL_VERSION: " << glGetString(GL_VERSION) << std::endl;
        std::cout << "GL_VENDOR: " << glGetString(GL_VENDOR) << std::endl;
        std::cout << "GL_RENDERER: " << glGetString(GL_RENDERER) << std::endl;

        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO();
        (void)io;

        ImGui::StyleColorsDark();

        ImGui_ImplGlfw_InitForOpenGL(g_window, true);
        ImGui_ImplOpenGL3_Init("#version 330 core");

        return true;
    }


    bool InitializeHW3() 
    {
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

        // 1. 创建窗口：
        g_window = glfwCreateWindow(g_width, g_height, "GAMES102 hw3", nullptr, nullptr);
        if (g_window == nullptr) 
        {
            std::cerr << "Failed to create GLFW window" << std::endl;
            glfwTerminate();
            return false;
        }

        // 2. 
        glfwMakeContextCurrent(g_window);
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) 
        {
            std::cout << "Failed to load glad" << std::endl;
            glfwTerminate();
            return false;
        }

        // 3. 
        glViewport(0, 0, g_width, g_height);
        std::cout << "GL_VERSION: " << glGetString(GL_VERSION) << std::endl;
        std::cout << "GL_VENDOR: " << glGetString(GL_VENDOR) << std::endl;
        std::cout << "GL_RENDERER: " << glGetString(GL_RENDERER) << std::endl;

        // 4. 
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO();
        (void)io;

        // 5. 
        ImGui::StyleColorsDark();

        ImGui_ImplGlfw_InitForOpenGL(g_window, true);
        ImGui_ImplOpenGL3_Init("#version 330 core");

        return true;
    }


    void BeginFrame()
    {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    }


    void EndFrame() {
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(g_window);
    }


    void Finalize() {
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();

        glfwDestroyWindow(g_window);
        glfwTerminate();
    }
}
using namespace MY_IMGUI;


// 曲线类；
struct CurveData 
{
    ImU32 line_color = 0;
    ImU32 point_color = 0;
    std::vector<float> pos_x;                   // 
    std::vector<float> pos_y;
    std::vector<float> in_pos_t;
    bool visible = false;

    void Clear() 
    {
        pos_x.clear();
        pos_y.clear();
        in_pos_t.clear();
    }


    // 在画布上画曲线：
    void Plot(ImDrawList* draw_list) const 
    {
        if (visible && !pos_x.empty())  
            PlotLineSegments(pos_x, pos_y, draw_list, line_color, point_color); 
    }


    void PlotXT(const std::vector<float>& pos_t, ImDrawList* draw_list, ImVec2 canvas_ul, ImVec2 canvas_br) const 
    {
        const size_t n = pos_x.size();
        if (!visible || n == 0)  
            return;

        std::vector<float> pos_x_t = pos_x;
        std::vector<float> pos_t_t = pos_t;
        for (size_t i = 0; i < n; i++) 
        {
            pos_x_t[i] = pos_x_t[i] / (g_canvas_pos_br.x - g_canvas_pos_ul.x) * (canvas_br.y - canvas_ul.y);
            pos_t_t[i] = pos_t_t[i] * (canvas_br.x - canvas_ul.x);
        }
        PlotLineSegments(pos_t_t, pos_x_t, draw_list, line_color, point_color, canvas_ul, canvas_br);
    }


    void PlotYT(const std::vector<float>& pos_t, ImDrawList* draw_list, ImVec2 canvas_ul, ImVec2 canvas_br) const
    {
        const size_t n = pos_y.size();
        if (!visible || n == 0) 
            return;
        
        std::vector<float> pos_y_t = pos_y;
        std::vector<float> pos_t_t = pos_t;
        for (size_t i = 0; i < n; i++) 
        {
            pos_y_t[i] = pos_y_t[i] / (g_canvas_pos_br.y - g_canvas_pos_ul.y) * (canvas_br.y - canvas_ul.y);
            pos_t_t[i] = pos_t_t[i] * (canvas_br.x - canvas_ul.x);
        }
        PlotLineSegments(pos_t_t, pos_y_t, draw_list, line_color, point_color, canvas_ul, canvas_br);
    }
};


namespace TEST_IMGUI
{
	void test0()
    {
        // 0. initialize and prepare data:
        if (!Initialize())
            return;

        std::vector<Eigen::Vector2f> in_pos;        // 样本点数组

        struct
        {
            std::vector<Eigen::Vector2f> pos;
            bool visible = false;
        } inter_poly;                         // 多项式插值拟合曲线

        struct
        {
            std::vector<Eigen::Vector2f> pos;
            int m = 0;
            int m_temp = 0;
            float sigma2 = 10000.0f;
            float sigma2_temp = 10000.0f;
            bool visible = false;
            bool update = false;
        } inter_gauss;                 // 高斯线性插值拟合曲线

        struct
        {
            std::vector<Eigen::Vector2f> pos;
            int m = 0;
            int m_temp = 0;
            bool visible = false;
            bool update = false;
        } approx_poly;                // 幂基函数最小二乘逼近拟合曲线

        struct
        {
            std::vector<Eigen::Vector2f> pos;
            int m = 0;
            int m_temp = 0;
            float lambda = 10.0f;
            float lambda_temp = 10.0f;
            bool visible = false;
            bool update = false;
        } approx_norm;              // 幂基函数岭回归逼近拟合曲线

        // 1. 准备窗口对象：
        ImGuiIO& io = ImGui::GetIO();           // 获取io使用权？？？
        ImFontConfig font_config;
        font_config.SizePixels = 24.0f;                             // 像素大小
        io.Fonts->AddFontDefault(&font_config);            // 设置字体；

        // 2. 窗口循环：
        while (!glfwWindowShouldClose(g_window))
        {
            BeginFrame();

            if (glfwGetKey(g_window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
                break;

            if (ImGui::Begin("Config"))
            {
                ImGui::Checkbox("Interpolation - Polygon (Red)", &inter_poly.visible);
                ImGui::Checkbox("Interpolation - Gauss (Green)", &inter_gauss.visible);

                ImGui::SameLine();
                ImGui::PushItemWidth(500.0f);
                ImGui::InputInt("m##1", &inter_gauss.m_temp);
                inter_gauss.m_temp = std::clamp(inter_gauss.m_temp, \
                    0, std::max<int>(0, inter_gauss.pos.size() - 1));
                if (inter_gauss.m_temp != inter_gauss.m)
                {
                    inter_gauss.m = inter_gauss.m_temp;
                    inter_gauss.update = true;
                }

                ImGui::SameLine();
                ImGui::InputFloat("sigma2", &inter_gauss.sigma2_temp);
                inter_gauss.sigma2_temp = std::max(inter_gauss.sigma2_temp, 1.0f);
                if (inter_gauss.sigma2_temp != inter_gauss.sigma2)
                {
                    inter_gauss.sigma2 = inter_gauss.sigma2_temp;
                    inter_gauss.update = true;
                }

                ImGui::Checkbox("Fitting - Polygon (Blue)", &approx_poly.visible);
                ImGui::SameLine();
                ImGui::InputInt("m##2", &approx_poly.m_temp);
                approx_poly.m_temp = std::clamp(approx_poly.m_temp, 0, std::max<int>(0, approx_poly.pos.size() - 1));
                if (approx_poly.m_temp != approx_poly.m)
                {
                    approx_poly.m = approx_poly.m_temp;
                    approx_poly.update = true;
                }

                ImGui::Checkbox("Fitting - Normalized (Yellow)", &approx_norm.visible);
                ImGui::SameLine();
                ImGui::InputInt("m##3", &approx_norm.m_temp);
                approx_norm.m_temp = std::clamp(approx_norm.m_temp, 0, std::max<int>(0, inter_gauss.pos.size() - 1));
                if (approx_norm.m_temp != approx_norm.m)
                {
                    approx_norm.m = approx_norm.m_temp;
                    approx_norm.update = true;
                }

                ImGui::SameLine();
                ImGui::InputFloat("lambda", &approx_norm.lambda_temp);
                approx_norm.lambda_temp = std::max(approx_norm.lambda_temp, 0.0f);
                if (approx_norm.lambda_temp != approx_norm.lambda)
                {
                    approx_norm.lambda = approx_norm.lambda_temp;
                    approx_norm.update = true;
                }

                ImGui::PopItemWidth();
                ImGui::End();
            }

            if (ImGui::Begin("Canvas"))
            {
                g_canvas_pos_ul = ImGui::GetCursorScreenPos();
                ImVec2 canvas_size = ImGui::GetContentRegionAvail();
                if (canvas_size.x < 50.0f)
                    canvas_size.x = 50.0f;

                if (canvas_size.y < 50.0f)
                    canvas_size.y = 50.0f;

                g_canvas_pos_br = ImVec2(g_canvas_pos_ul.x + canvas_size.x, g_canvas_pos_ul.y + canvas_size.y);

                ImGuiIO& io = ImGui::GetIO();
                ImDrawList* draw_list = ImGui::GetWindowDrawList();
                draw_list->AddRectFilled(g_canvas_pos_ul, g_canvas_pos_br, IM_COL32(50, 50, 50, 255));
                draw_list->AddRect(g_canvas_pos_ul, g_canvas_pos_br, IM_COL32(255, 255, 255, 255));

                float step = 20.0f;
                float lb = step;
                float rb = g_canvas_pos_br.x - step - g_canvas_pos_ul.x;
                ImGui::InvisibleButton("canvas", canvas_size);
                const bool is_hovered = ImGui::IsItemHovered();
                if (is_hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
                {
                    in_pos.emplace_back(io.MousePos.x - g_canvas_pos_ul.x, g_canvas_pos_br.y - io.MousePos.y);
                    std::sort(in_pos.begin(), in_pos.end(), [](const Eigen::Vector2f& a, const Eigen::Vector2f& b)
                        {
                            return a.x() < b.x();
                        });

                    inter_poly.pos = MathUtil::InterpolationPolygon(in_pos, lb, rb, step);
                    inter_gauss.update = true;
                    approx_poly.update = true;
                    approx_norm.update = true;
                }
                else if (is_hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Right))
                {
                    Eigen::Vector2f pos(io.MousePos.x - g_canvas_pos_ul.x, g_canvas_pos_br.y - io.MousePos.y);
                    size_t index = 0;
                    float min_dist = std::numeric_limits<float>::max();
                    for (size_t i = 0; i < in_pos.size(); i++)
                    {
                        float dist = (pos - in_pos[i]).squaredNorm();
                        if (dist < min_dist)
                        {
                            min_dist = dist;
                            index = i;
                        }
                    }
                    if (min_dist <= 100.0f)
                    {
                        in_pos.erase(in_pos.begin() + index);
                        inter_poly.pos = MathUtil::InterpolationPolygon(in_pos, lb, rb, step);
                        inter_gauss.update = true;
                        approx_poly.update = true;
                        approx_norm.update = true;
                    }
                }
                if (inter_gauss.update)
                {
                    inter_gauss.pos = MathUtil::InterpolationGauss(in_pos, inter_gauss.sigma2, inter_gauss.m, lb, rb, step);
                    inter_gauss.update = false;
                }
                if (approx_poly.update)
                {
                    approx_poly.pos = MathUtil::ApproximationPolygon(in_pos, approx_poly.m, lb, rb, step);
                    approx_poly.update = false;
                }
                if (approx_norm.update)
                {
                    approx_norm.pos = MathUtil::ApproximationNormalized(in_pos, approx_norm.m, approx_norm.lambda,
                        lb, rb, step);
                    approx_norm.update = false;
                }

                if (inter_poly.visible)
                {
                    PlotLineSegments(inter_poly.pos, draw_list, IM_COL32(255, 50, 50, 255), IM_COL32(255, 80, 80, 255));
                }
                if (inter_gauss.visible)
                {
                    PlotLineSegments(inter_gauss.pos, draw_list, IM_COL32(50, 255, 50, 255), IM_COL32(80, 255, 80, 255));
                }
                if (approx_poly.visible)
                {
                    PlotLineSegments(approx_poly.pos, draw_list, IM_COL32(50, 50, 255, 255), IM_COL32(80, 80, 255, 255));
                }
                if (approx_norm.visible)
                {
                    PlotLineSegments(approx_norm.pos, draw_list, IM_COL32(255, 255, 50, 255), IM_COL32(255, 255, 80, 255));
                }
                PlotLineSegments(in_pos, draw_list, IM_COL32(255, 255, 255, 0), IM_COL32(255, 255, 255, 255));

                ImGui::End();
            }

            EndFrame();
        }

        Finalize();
    }


    void test1()
    {
        // 0. initialize & prepare data:
        if (!InitializeHW3())
            return; 

        //      样本点数据；
        std::vector<float> in_pos_x;        
        std::vector<float> in_pos_y;

        //      插值拟合曲线数据：
        bool inter_update = false;              // 更新标识符，插入/删除顶点时会变为true;
        CurveData inter_uniform{ IM_COL32(255, 50, 50, 255), IM_COL32(255, 80, 80, 255) };
        CurveData inter_chordal{ IM_COL32(50, 255, 50, 255), IM_COL32(80, 255, 80, 255) };
        CurveData inter_centripetal{ IM_COL32(50, 50, 255, 255), IM_COL32(80, 80, 255, 255) };
        CurveData inter_foley{ IM_COL32(150, 150, 255, 255), IM_COL32(180, 180, 255, 255) };

        //      逼近拟合曲线数据：
        bool approx_update = false;
        int approx_m = 0;
        int approx_m_temp = 0;
        CurveData approx_uniform{ IM_COL32(50, 255, 255, 255), IM_COL32(80, 255, 255, 255) };
        CurveData approx_chordal{ IM_COL32(255, 50, 255, 255), IM_COL32(255, 80, 255, 255) };
        CurveData approx_centripetal{ IM_COL32(255, 255, 50, 255), IM_COL32(255, 255, 80, 255) };
        CurveData approx_foley{ IM_COL32(255, 150, 150, 255), IM_COL32(255, 180, 180, 255) };

        //      样条插值拟合曲线数据：
        bool spline_update = false;
        CurveData spline_uniform{ IM_COL32(255, 100, 50, 255), IM_COL32(255, 130, 80, 255) };
        CurveData spline_chordal{ IM_COL32(50, 255, 100, 255), IM_COL32(80, 255, 130, 255) };
        CurveData spline_centripetal{ IM_COL32(100, 50, 255, 255), IM_COL32(130, 80, 255, 255) };
        CurveData spline_foley{ IM_COL32(170, 255, 170, 255), IM_COL32(200, 255, 200, 255) };

        //      自己写的三次B样条插值：
        bool my_CBS_update = false;
        CurveData mySpline{ IM_COL32(100, 50, 255, 255), IM_COL32(130, 80, 255, 255) };
        CurveData mySpline_loop{ IM_COL32(170, 255, 170, 255), IM_COL32(200, 255, 200, 255) };

        // 设定采样空间（画布空间）   
        float lb = 0.0f;                    
        float rb = 1.0f;
        float step = (rb - lb) / 60.0f;
        std::vector<float> out_pos_t;
        for (float x = lb; x <= rb; x += step) 
            out_pos_t.push_back(x);       
        
        ImGuiIO& io = ImGui::GetIO();
        ImFontConfig font_config;
        font_config.SizePixels = 24.0f;
        io.Fonts->AddFontDefault(&font_config);

        // 1. 主循环：
        bool blIsLoop = true;
        bool blShowMCS = false;
        while (!glfwWindowShouldClose(g_window)) 
        {
            // w1.
            BeginFrame();

            // w2.
            if (glfwGetKey(g_window, GLFW_KEY_ESCAPE) == GLFW_PRESS) 
                break;
            
            // w3. 添加各种控件：
            if (ImGui::Begin("Config"))
            {
                // w3.1. 插值拟合控件栏
                ImGui::Text("Interpolation");
                ImGui::ColorButton("##1", ImColor(inter_uniform.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Uniform##1", &inter_uniform.visible);
                ImGui::ColorButton("##2", ImColor(inter_chordal.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Chordal##1", &inter_chordal.visible);
                ImGui::ColorButton("##3", ImColor(inter_centripetal.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Centripetal##1", &inter_centripetal.visible);
                ImGui::ColorButton("##4", ImColor(inter_foley.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Foley##1", &inter_foley.visible);

                // w3.2. 逼近拟合控件栏
                ImGui::Separator();
                ImGui::Text("Approx");
                ImGui::ColorButton("##5", ImColor(approx_uniform.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Uniform##2", &approx_uniform.visible);
                ImGui::ColorButton("##6", ImColor(approx_chordal.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Chordal##2", &approx_chordal.visible);
                ImGui::ColorButton("##7", ImColor(approx_centripetal.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Centripetal##2", &approx_centripetal.visible);
                ImGui::ColorButton("##8", ImColor(approx_foley.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Foley##2", &approx_foley.visible);
                ImGui::InputInt("m", &approx_m_temp);
                approx_m_temp = std::clamp(approx_m_temp, 0, std::max<int>(0, in_pos_x.size() - 1));
                if (approx_m_temp != approx_m) 
                {
                    approx_m = approx_m_temp;
                    approx_update = true;
                }

                // w3.3. 样条插值控件栏
                ImGui::Separator();
                ImGui::Text("Cubic Spline");
                ImGui::ColorButton("##9", ImColor(spline_uniform.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Uniform##3", &spline_uniform.visible);
                ImGui::ColorButton("##10", ImColor(spline_chordal.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Chordal##3", &spline_chordal.visible);
                ImGui::ColorButton("##11", ImColor(spline_centripetal.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Centripetal##3", &spline_centripetal.visible);
                ImGui::ColorButton("##12", ImColor(spline_foley.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("Foley##3", &spline_foley.visible);

                // w3.4. 自己写的三次B样条插值控件栏
                ImGui::Separator();
                ImGui::Text("my Cubic Spline");
                ImGui::ColorButton("##13", ImColor(mySpline.line_color).Value);
                ImGui::SameLine();
                ImGui::Checkbox("MyCBS##3", &blShowMCS);
                ImGui::SameLine();
                ImGui::Checkbox("isLoop##4", &blIsLoop);
                mySpline.visible = blShowMCS && (!blIsLoop);
                mySpline_loop.visible = blShowMCS && blIsLoop;

                // 重置按钮
                ImGui::Separator();
                bool reset = ImGui::Button("Reset");
                if (reset)
                {
                    in_pos_x.clear();
                    in_pos_y.clear();
                    inter_uniform.Clear();
                    inter_chordal.Clear();
                    inter_centripetal.Clear();
                    inter_foley.Clear();
                    approx_uniform.Clear();
                    approx_chordal.Clear();
                    approx_centripetal.Clear();
                    approx_foley.Clear();
                    spline_uniform.Clear();
                    spline_chordal.Clear();
                    spline_centripetal.Clear();
                    spline_foley.Clear();
                    mySpline.Clear();
                    mySpline_loop.Clear();
                }
                ImGui::End();
            }

            // w4. 主画布
            if (ImGui::Begin("Canvas")) 
            {
                // w4.1. 设置画布尺寸：
                g_canvas_pos_ul = ImGui::GetCursorScreenPos();
                ImVec2 canvas_size = ImGui::GetContentRegionAvail();
                if (canvas_size.x < 50.0f) 
                    canvas_size.x = 50.0f;                
                if (canvas_size.y < 50.0f) 
                    canvas_size.y = 50.0f;                
                g_canvas_pos_br = ImVec2(g_canvas_pos_ul.x + canvas_size.x, g_canvas_pos_ul.y + canvas_size.y);

                // w4.2
                ImGuiIO& io = ImGui::GetIO();
                ImDrawList* draw_list = ImGui::GetWindowDrawList();
                draw_list->AddRectFilled(g_canvas_pos_ul, g_canvas_pos_br, IM_COL32(50, 50, 50, 255));
                draw_list->AddRect(g_canvas_pos_ul, g_canvas_pos_br, IM_COL32(255, 255, 255, 255));

                // w4.3
                ImGui::InvisibleButton("canvas", canvas_size);
                const bool is_hovered = ImGui::IsItemHovered();
                if (is_hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
                {
                    in_pos_x.push_back(io.MousePos.x - g_canvas_pos_ul.x);
                    in_pos_y.push_back(g_canvas_pos_br.y - io.MousePos.y);
                    inter_update = true;
                    approx_update = true;
                    spline_update = true;
                    my_CBS_update = true;
                }
                else if (is_hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) 
                {
                    Eigen::Vector2f pos(io.MousePos.x - g_canvas_pos_ul.x, g_canvas_pos_br.y - io.MousePos.y);
                    size_t index = 0;
                    float min_dist = std::numeric_limits<float>::max();
                    for (size_t i = 0; i < in_pos_x.size(); i++)
                    {
                        float dist = (pos - Eigen::Vector2f(in_pos_x[i], in_pos_y[i])).squaredNorm();
                        if (dist < min_dist) 
                        {
                            min_dist = dist;
                            index = i;
                        }
                    }
                    if (min_dist <= 100.0f) 
                    {
                        in_pos_x.erase(in_pos_x.begin() + index);
                        in_pos_y.erase(in_pos_y.begin() + index);
                        inter_update = true;
                        approx_update = true;
                        spline_update = true;
                        my_CBS_update = true;
                    }
                }

                // w4.4. 更新插值拟合曲线
                if (inter_update) 
                {
                    inter_uniform.in_pos_t = MathUtil::ParameterizationUniform(in_pos_x, in_pos_y);
                    inter_uniform.pos_x = MathUtil::InterpolationPolygon(inter_uniform.in_pos_t, in_pos_x, lb, rb, step);
                    inter_uniform.pos_y = MathUtil::InterpolationPolygon(inter_uniform.in_pos_t, in_pos_y, lb, rb, step);

                    inter_chordal.in_pos_t = MathUtil::ParameterizationChoral(in_pos_x, in_pos_y);
                    inter_chordal.pos_x = MathUtil::InterpolationPolygon(inter_chordal.in_pos_t, in_pos_x, lb, rb, step);
                    inter_chordal.pos_y = MathUtil::InterpolationPolygon(inter_chordal.in_pos_t, in_pos_y, lb, rb, step);

                    inter_centripetal.in_pos_t = MathUtil::ParameterizationCentripetal(in_pos_x, in_pos_y);
                    inter_centripetal.pos_x = MathUtil::InterpolationPolygon(inter_centripetal.in_pos_t, in_pos_x, lb, rb, step);
                    inter_centripetal.pos_y = MathUtil::InterpolationPolygon(inter_centripetal.in_pos_t, in_pos_y, lb, rb, step);

                    inter_foley.in_pos_t = MathUtil::ParameterizationFoley(in_pos_x, in_pos_y);
                    inter_foley.pos_x = MathUtil::InterpolationPolygon(inter_foley.in_pos_t, in_pos_x, lb, rb, step);
                    inter_foley.pos_y = MathUtil::InterpolationPolygon(inter_foley.in_pos_t, in_pos_y, lb, rb, step);

                    inter_update = false;
                }

                // w4.5. 更新逼近拟合曲线：
                if (approx_update)
                {
                    approx_uniform.in_pos_t = inter_uniform.in_pos_t;
                    approx_uniform.pos_x = MathUtil::ApproximationPolygon(approx_uniform.in_pos_t, in_pos_x,
                        approx_m, lb, rb, step);
                    approx_uniform.pos_y = MathUtil::ApproximationPolygon(approx_uniform.in_pos_t, in_pos_y,
                        approx_m, lb, rb, step);

                    approx_chordal.in_pos_t = inter_chordal.in_pos_t;
                    approx_chordal.pos_x = MathUtil::ApproximationPolygon(approx_chordal.in_pos_t, in_pos_x,
                        approx_m, lb, rb, step);
                    approx_chordal.pos_y = MathUtil::ApproximationPolygon(approx_chordal.in_pos_t, in_pos_y,
                        approx_m, lb, rb, step);

                    approx_centripetal.in_pos_t = inter_centripetal.in_pos_t;
                    approx_centripetal.pos_x = MathUtil::ApproximationPolygon(approx_centripetal.in_pos_t, in_pos_x,
                        approx_m, lb, rb, step);
                    approx_centripetal.pos_y = MathUtil::ApproximationPolygon(approx_centripetal.in_pos_t, in_pos_y,
                        approx_m, lb, rb, step);

                    approx_foley.in_pos_t = inter_foley.in_pos_t;
                    approx_foley.pos_x = MathUtil::ApproximationPolygon(approx_foley.in_pos_t, in_pos_x,
                        approx_m, lb, rb, step);
                    approx_foley.pos_y = MathUtil::ApproximationPolygon(approx_foley.in_pos_t, in_pos_y,
                        approx_m, lb, rb, step);

                    approx_update = false;
                }

                // w4.6. 更新样条插值拟合曲线：
                if (spline_update)
                {
                    spline_uniform.in_pos_t = inter_uniform.in_pos_t;
                    spline_uniform.pos_x = MathUtil::CubicSpline(spline_uniform.in_pos_t, in_pos_x, lb, rb, step);
                    spline_uniform.pos_y = MathUtil::CubicSpline(spline_uniform.in_pos_t, in_pos_y, lb, rb, step);

                    spline_chordal.in_pos_t = inter_chordal.in_pos_t;
                    spline_chordal.pos_x = MathUtil::CubicSpline(spline_chordal.in_pos_t, in_pos_x, lb, rb, step);
                    spline_chordal.pos_y = MathUtil::CubicSpline(spline_chordal.in_pos_t, in_pos_y, lb, rb, step);

                    spline_centripetal.in_pos_t = inter_centripetal.in_pos_t;
                    spline_centripetal.pos_x = MathUtil::CubicSpline(spline_centripetal.in_pos_t, in_pos_x, lb, rb, step);
                    spline_centripetal.pos_y = MathUtil::CubicSpline(spline_centripetal.in_pos_t, in_pos_y, lb, rb, step);

                    spline_foley.in_pos_t = inter_foley.in_pos_t;
                    spline_foley.pos_x = MathUtil::CubicSpline(spline_foley.in_pos_t, in_pos_x, lb, rb, step);
                    spline_foley.pos_y = MathUtil::CubicSpline(spline_foley.in_pos_t, in_pos_y, lb, rb, step);

                    spline_update = false;
                }

                // w4.7 更新自己写的三次B样条插值曲线：
                if (my_CBS_update)
                {
                    Eigen::VectorXf tmpVec;
                    Eigen::MatrixXf versFitted;
                    Eigen::MatrixXf versSample(in_pos_x.size(), 3);  
                    vec2EigenVec(tmpVec, in_pos_x);
                    versSample.col(0) = tmpVec;
                    vec2EigenVec(tmpVec, in_pos_y);
                    versSample.col(1) = tmpVec;
                    versSample.col(2).setZero();

                    cubicBSplineInterpCurve3D(versFitted, versSample, 100); 
                    eigenVec2Vec(mySpline.pos_x, versFitted.col(0));
                    eigenVec2Vec(mySpline.pos_y, versFitted.col(1));
                     
                    cubicBSplineInterpLoop3D(versFitted, versSample, 100);
                    eigenVec2Vec(mySpline_loop.pos_x, versFitted.col(0));
                    eigenVec2Vec(mySpline_loop.pos_y, versFitted.col(1));
 
                    my_CBS_update = false;
                }

                // 生成绘图曲线：
                inter_uniform.Plot(draw_list);
                inter_chordal.Plot(draw_list);
                inter_centripetal.Plot(draw_list);
                inter_foley.Plot(draw_list);

                approx_uniform.Plot(draw_list);
                approx_chordal.Plot(draw_list);
                approx_centripetal.Plot(draw_list);
                approx_foley.Plot(draw_list);

                spline_uniform.Plot(draw_list);
                spline_chordal.Plot(draw_list);
                spline_centripetal.Plot(draw_list);
                spline_foley.Plot(draw_list);

                mySpline.Plot(draw_list);
                mySpline_loop.Plot(draw_list);

                PlotLineSegments(in_pos_x, in_pos_y, draw_list, IM_COL32(255, 255, 255, 0), IM_COL32(255, 255, 255, 255));
                ImGui::End();
            }

            // w5. X-T画布
            if (ImGui::Begin("X-T"))
            {
                ImVec2 canvas_pos_ul = ImGui::GetCursorScreenPos();
                ImVec2 canvas_size = ImGui::GetContentRegionAvail();
                if (canvas_size.x < 50.0f) 
                    canvas_size.x = 50.0f;
                
                if (canvas_size.y < 50.0f) 
                    canvas_size.y = 50.0f;
                
                ImVec2 canvas_pos_br = ImVec2(canvas_pos_ul.x + canvas_size.x, canvas_pos_ul.y + canvas_size.y);

                ImDrawList* draw_list = ImGui::GetWindowDrawList();
                draw_list->AddRectFilled(canvas_pos_ul, canvas_pos_br, IM_COL32(50, 50, 50, 255));
                draw_list->AddRect(canvas_pos_ul, canvas_pos_br, IM_COL32(255, 255, 255, 255));

                inter_uniform.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                inter_chordal.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                inter_centripetal.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                inter_foley.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);

                approx_uniform.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                approx_chordal.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                approx_centripetal.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                approx_foley.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);

                spline_uniform.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                spline_chordal.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                spline_centripetal.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                spline_foley.PlotXT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);

                ImGui::End();
            }

            // w6. Y-T画布
            if (ImGui::Begin("Y-T")) 
            {
                ImVec2 canvas_pos_ul = ImGui::GetCursorScreenPos();
                ImVec2 canvas_size = ImGui::GetContentRegionAvail();
                if (canvas_size.x < 50.0f) 
                    canvas_size.x = 50.0f;
                
                if (canvas_size.y < 50.0f) 
                    canvas_size.y = 50.0f;
                
                ImVec2 canvas_pos_br = ImVec2(canvas_pos_ul.x + canvas_size.x, canvas_pos_ul.y + canvas_size.y);

                ImDrawList* draw_list = ImGui::GetWindowDrawList();
                draw_list->AddRectFilled(canvas_pos_ul, canvas_pos_br, IM_COL32(50, 50, 50, 255));
                draw_list->AddRect(canvas_pos_ul, canvas_pos_br, IM_COL32(255, 255, 255, 255));

                inter_uniform.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                inter_chordal.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                inter_centripetal.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                inter_foley.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);

                approx_uniform.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                approx_chordal.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                approx_centripetal.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                approx_foley.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);

                spline_uniform.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                spline_chordal.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                spline_centripetal.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);
                spline_foley.PlotYT(out_pos_t, draw_list, canvas_pos_ul, canvas_pos_br);

                ImGui::End();
            }

            EndFrame();
        }

        Finalize(); 
    }


    void test11() 
    {
        Eigen::MatrixXf versSample(3, 3);
        Eigen::MatrixXf versFitted;
        versSample << 1, 2, 0, 2, 3, 0, 4, 5, 0;

        cubicBSplineInterpCurve3D(versFitted, versSample, 100);

        debugWriteVers("versIn", versSample);
        debugWriteVers("versOut", versFitted);
 

        debugDisp("test11 finished.");
    }

    void test2()
    {
        MyApp app;
        app.initAll();
        app.mainLoop();
        app.cleanUp();

    }
}
