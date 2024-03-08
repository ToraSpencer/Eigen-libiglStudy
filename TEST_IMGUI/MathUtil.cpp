#include "TEST_IMGUI/MathUtil.h"
#include <vector>
#include "Eigen/Dense"


namespace 
{
constexpr float Sqr(float x)
{
    return x * x;
}


Eigen::MatrixXf LeastSquares(const std::vector<Eigen::Vector2f>&in_pos,\
    int m, float lambda = 0.0f) 
{
    const int n = in_pos.size();
    Eigen::MatrixXf A(m + 1, m + 1);
    std::vector<float> pow_temp(n, 1.0f);
    for (int i = 0; i < 2 * m + 1; i++) 
    {
        float sum = 0;
        for (int j = 0; j < n; j++) 
        {
            sum += pow_temp[j];
            pow_temp[j] *= in_pos[j].x();
        }
        for (int j = 0; j <= i; j++) 
            if (j <= m && i - j <= m) 
                A(j, i - j) = sum;       
    }

    Eigen::MatrixXf norm = Eigen::MatrixXf::Identity(m + 1, m + 1);
    A += lambda * norm;

    Eigen::MatrixXf Y(m + 1, 1);
    std::fill(pow_temp.begin(), pow_temp.end(), 1.0f);
    for (int i = 0; i <= m; i++) 
    {
        Y(i, 0) = 0.0f;
        for (int j = 0; j < n; j++) 
        {
            Y(i, 0) += in_pos[j].y() * pow_temp[j];
            pow_temp[j] *= in_pos[j].x();
        }
    }

    Eigen::MatrixXf B = A.inverse() * Y;
    return B;
}


Eigen::VectorXf LeastSquares(const std::vector<float>& pos_x,\
    const std::vector<float>& pos_y, int m) 
{
    const int n = pos_x.size();
    Eigen::MatrixXf A(m + 1, m + 1);
    std::vector<float> pow_temp(n, 1.0f);
    for (int i = 0; i < 2 * m + 1; i++) {
        float sum = 0;
        for (int j = 0; j < n; j++) {
            sum += pow_temp[j];
            pow_temp[j] *= pos_x[j];
        }
        for (int j = 0; j <= i; j++) {
            if (j <= m && i - j <= m) {
                A(j, i - j) = sum;
            }
        }
    }

    Eigen::MatrixXf norm = Eigen::MatrixXf::Identity(m + 1, m + 1);

    Eigen::MatrixXf Y(m + 1, 1);
    std::fill(pow_temp.begin(), pow_temp.end(), 1.0f);
    for (int i = 0; i <= m; i++) {
        Y(i, 0) = 0.0f;
        for (int j = 0; j < n; j++) {
            Y(i, 0) += pos_y[j] * pow_temp[j];
            pow_temp[j] *= pos_x[j];
        }
    }

    Eigen::VectorXf B = A.colPivHouseholderQr().solve(Y);
    return B;
}

}


// 多项式插值
std::vector<Eigen::Vector2f> MathUtil::InterpolationPolygon(\
    const std::vector<Eigen::Vector2f> &in_pos, float lb, float rb, float step) 
{
    /*
    MathUtil::InterpolationPolygon(
            const std::vector<Eigen::Vector2f> &in_pos, 
            float lb, 
            float rb, 
            float step
            ) 
    
    */
    std::vector<Eigen::Vector2f> result;
    for (float x = lb; x <= rb; x += step) 
    {
        float y = 0;
        for (int i = 0; i < in_pos.size(); i++) 
        {
            float temp = in_pos[i].y();
            for (int j = 0; j < in_pos.size(); j++)
                if (i != j) 
                    temp = temp * (x - in_pos[j].x()) / (in_pos[i].x() - in_pos[j].x());
            y += temp;
        }
        result.emplace_back(x, y);
    }
    return result;
}


// 多项式插值
/*
    std::vector<float> MathUtil::InterpolationPolygon(                返回参数曲线的采样函数值
        const std::vector<float>& pos_t,            样本点对应的参数空间t坐标序列；
        const std::vector<float>& pos_u,           样本点的x坐标或y坐标序列；
        float lb,                   t区间左边界
        float rb,                   t区间右边界
        float step                采样步长；
        ) 
*/
std::vector<float> MathUtil::InterpolationPolygon(\
    const std::vector<float>& pos_t, const std::vector<float>& pos_u,
    float lb, float rb, float step) 
{
    const int m = pos_t.size();                           //  样本容量；
    const size_t n = std::ceil((rb - lb) / step) + 1;
    std::vector<float> result(n);       

    for (float t = lb; t <= rb; t += step) 
    {
        float u = 0;
        for (int i = 0; i < m; i++)
        {
            float temp = pos_u[i];
            for (int j = 0; j < m; j++) 
                if (i != j)                 
                    temp = temp * (t - pos_t[j]) / (pos_t[i] - pos_t[j]);                           
            u += temp;
        }
        result.push_back(u);
    } 

    return result;
}


// 高斯插值
std::vector<Eigen::Vector2f> MathUtil::InterpolationGauss(\
    const std::vector<Eigen::Vector2f> &in_pos, float sigma2, \
    int m, float lb, float rb, float step) {
    const int n = in_pos.size();
    m = std::min(m, std::max(n - 1, 0));

    Eigen::MatrixXf B_poly = LeastSquares(in_pos, m);
    std::vector<float> y_approx(n);
    for (int i = 0; i < n; i++) {
        float y = 0, x_temp = 1.0f;
        for (int j = 0; j <= m; j++) {
            y += B_poly(j, 0) * x_temp;
            x_temp *= in_pos[i].x();
        }
        y_approx[i] = y;
    }

    Eigen::MatrixXf A(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A(i, j) = std::exp(-Sqr(in_pos[i].x() - in_pos[j].x()) / (2 * sigma2));
        }
    }

    Eigen::MatrixXf Y(n, 1);
    for (int i = 0; i < n; i++) {
        Y(i, 0) = in_pos[i].y() - y_approx[i];
    }

    Eigen::MatrixXf B = A.inverse() * Y;
    std::vector<Eigen::Vector2f> result;
    for (float x = lb; x <= rb; x += step) {
        float y = 0, x_temp = 1.0f;
        for (int i = 0; i <= m; i++) {
            y += B_poly(i, 0) * x_temp;
            x_temp *= x;
        }
        for (int i = 0; i < n; i++) {
            y += B(i, 0) * std::exp(-Sqr(x - in_pos[i].x()) / (2 * sigma2));
        }
        result.emplace_back(x, y);
    }
    return result;
}


// 三次样条插值拟合：
std::vector<float> MathUtil::CubicSpline(\
    const std::vector<float>& pos_x, const std::vector<float>& pos_y,
    float lb, float rb, float step)
{
    // 详情见GAMES102的数学资料2：三次样条插值函数
    const size_t m = pos_x.size();
    const size_t n = std::ceil((rb - lb) / step) + 1;

    // 0. simple case:
    if (m == 1)
        return { 0.0f };

    // 1. 
    std::vector<float> diff_x(m - 1);
    std::vector<float> coe0(m);
    std::vector<float> coe1(m);
    std::vector<float> coe2(m);
    std::vector<float> coe3(m);
    for (int i = 0; i < m - 1; i++)
        diff_x[i] = pos_x[i + 1] - pos_x[i];
    for (int i = 0; i < m; i++)
        coe0[i] = pos_y[i];

    Eigen::MatrixXf A(m, m);
    Eigen::VectorXf B(m);
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < m; j++)
            A(i, j) = 0.0f;

        if (i == 0)
        {
            A(i, i) = 1.0f;
            B(i, 0) = 0.0f;
        }
        else if (i == m - 1)
        {
            A(i, i) = 1.0f;
            B(i, 0) = 0.0f;
        }
        else
        {
            A(i, i) = 2 * (diff_x[i] + diff_x[i - 1]);
            A(i, i + 1) = diff_x[i];
            A(i, i - 1) = diff_x[i - 1];
            B(i, 0) = 3.0 / diff_x[i] * \
                (coe0[i + 1] - coe0[i]) - 3.0 / diff_x[i - 1] * (coe0[i] - coe0[i - 1]);
        }
    }

    Eigen::VectorXf C = A.colPivHouseholderQr().solve(B);
    for (int i = 0; i < m; i++)
        coe2[i] = C(i, 0);

    for (int i = 0; i < m - 1; i++)
    {
        coe3[i] = (coe2[i + 1] - coe2[i]) / 3.0 / diff_x[i];
        coe1[i] = (coe0[i + 1] - coe0[i]) / diff_x[i] - coe2[i] * diff_x[i] - coe3[i] * diff_x[i] * diff_x[i];
    }

    size_t curr = 0;
    std::vector<float> result(n);
    for (size_t i = 0; i < n; ++i)
    {
        float x = lb + i * step;
        if (x > pos_x[curr + 1])
            ++curr;
        float X = x - pos_x[curr];
        float y = coe0[curr] + coe1[curr] * X + coe2[curr] * X * X + coe3[curr] * X * X * X;
        result[i] = y;
    }


    return result;
}


// 最小二乘多项式逼近：
std::vector<Eigen::Vector2f> MathUtil::ApproximationPolygon(\
    const std::vector<Eigen::Vector2f> &in_pos, int m, float lb, float rb, float step) 
{
    const int n = in_pos.size();
    m = std::min(m, std::max(n - 1, 0));
    Eigen::MatrixXf B = LeastSquares(in_pos, m);
    std::vector<Eigen::Vector2f> result;
    for (float x = lb; x <= rb; x += step) 
    {
        float y = 0, x_temp = 1.0f;
        for (int i = 0; i <= m; i++) 
        {
            y += B(i, 0) * x_temp;
            x_temp *= x;
        }
        result.emplace_back(x, y);
    }
    return result;
}


// 最小二乘多项式逼近：
std::vector<float> MathUtil::ApproximationPolygon(\
    const std::vector<float>& pos_x, const std::vector<float>& pos_y,
    int m, float lb, float rb, float step) 
{
    const int n = pos_x.size();
    m = std::min(m, std::max(n - 1, 0));
    Eigen::VectorXf B = LeastSquares(pos_x, pos_y, m);
    std::vector<float> result;
    for (float x = lb; x <= rb; x += step) 
    {
        float y = 0, x_temp = 1.0f;
        for (int i = 0; i <= m; i++) 
        {
            y += B(i, 0) * x_temp;
            x_temp *= x;
        }
        result.push_back(y);
    }
    return result;
}


// 幂基函数岭回归逼近
std::vector<Eigen::Vector2f> MathUtil::ApproximationNormalized(\
    const std::vector<Eigen::Vector2f> &in_pos, int m, float lambda, \
    float lb, float rb, float step) {
    const int n = in_pos.size();
    m = std::min(m, std::max(n - 1, 0));
    Eigen::MatrixXf B = LeastSquares(in_pos, m, lambda);
    std::vector<Eigen::Vector2f> result;
    for (float x = lb; x <= rb; x += step) {
        float y = 0, x_temp = 1.0f;
        for (int i = 0; i <= m; i++) {
            y += B(i, 0) * x_temp;
            x_temp *= x;
        }
        result.emplace_back(x, y);
    }
    return result;
}




// 均匀(Uniform)点列参数化：
/*
    std::vector<float> ParameterizationUniform(\        返回向量为参数t的采样点列
        const std::vector<float>& pos_x, 
        const std::vector<float>& pos_y
        ) 

    输入曲线点云，输出参数t的参数化点列；
    t取值在[0, 1]之间；
    参数点列点数和输入曲线点云的点数相同；
    等距点列参数化中，参数点列中的点间隔相等；
*/
std::vector<float> MathUtil::ParameterizationUniform(\
    const std::vector<float>& pos_x, const std::vector<float>& pos_y) 
{
    const size_t n = pos_x.size();          
    if (n == 1) 
        return { 0.0f };
    
    float inv = 1.0f / (n - 1);
    std::vector<float> result(n);
    for (size_t i = 0; i < n; i++) 
        result[i] = i * inv;
    
    return result;
}


// 等弦长点列参数化：
/*
    std::vector<float> ParameterizationChoral(\        返回向量为参数t的采样点列
        const std::vector<float>& pos_x,
        const std::vector<float>& pos_y
        )

    输入曲线点云，输出参数t的参数化点列；
    t取值在[0, 1]之间；
    参数点列点数和输入曲线点云的点数相同；
    
*/
std::vector<float> MathUtil::ParameterizationChoral(\
    const std::vector<float>& pos_x, const std::vector<float>& pos_y) 
{
    const size_t n = pos_x.size();              // 采样点数
    if (n == 1) 
        return { 0.0f };
    
    float sum = 0.0f;
    std::vector<float> result(n);
    std::vector<float> dist(n - 1);
    for (size_t i = 1; i < n; i++) 
    {
        dist[i - 1] = (Eigen::Vector2f(pos_x[i - 1], pos_y[i - 1]) \
            - Eigen::Vector2f(pos_x[i], pos_y[i])).norm();                      // 相邻采样点的欧氏距离；
        sum += dist[i - 1];                                 
    }
    result[0] = 0.0f;
    for (size_t i = 1; i < n - 1; i++) 
    {
        result[i] = dist[i - 1] / sum;
        result[i] += result[i - 1];
    }
    result[n - 1] = 1.0f;
    return result;
}


std::vector<float> MathUtil::ParameterizationCentripetal(\
    const std::vector<float>& pos_x, const std::vector<float>& pos_y) {
    const size_t n = pos_x.size();
    if (n == 1) {
        return { 0.0f };
    }
    float sum = 0.0f;
    std::vector<float> result(n);
    std::vector<float> dist_sqrt(n - 1);
    for (size_t i = 1; i < n; i++) {
        float dist = (Eigen::Vector2f(pos_x[i - 1], pos_y[i - 1]) - Eigen::Vector2f(pos_x[i], pos_y[i])).norm();
        dist_sqrt[i - 1] = std::sqrt(dist);
        sum += dist_sqrt[i - 1];
    }
    result[0] = 0.0f;
    for (size_t i = 1; i < n - 1; i++) {
        result[i] = dist_sqrt[i - 1] / sum;
        result[i] += result[i - 1];
    }
    result[n - 1] = 1.0f;
    return result;
}


std::vector<float> MathUtil::ParameterizationFoley(\
    const std::vector<float>& pos_x, const std::vector<float>& pos_y) 
{
    const size_t n = pos_x.size();
    if (n == 1) {
        return { 0.0f };
    }
    std::vector<float> dist(n + 1);
    for (size_t i = 1; i < n; i++) {
        dist[i] = (Eigen::Vector2f(pos_x[i - 1], pos_y[i - 1]) - Eigen::Vector2f(pos_x[i], pos_y[i])).norm();
    }
    dist[0] = dist[n] = 0.0f;
    std::vector<float> angle(n);
    for (size_t i = 1; i < n - 1; i++) {
        Eigen::Vector2f a = Eigen::Vector2f(pos_x[i - 1], pos_y[i - 1]) - Eigen::Vector2f(pos_x[i], pos_y[i]);
        Eigen::Vector2f b = Eigen::Vector2f(pos_x[i + 1], pos_y[i + 1]) - Eigen::Vector2f(pos_x[i], pos_y[i]);
        angle[i] = a.dot(b) / dist[i] / dist[i + 1];
        angle[i] = std::min(Pi - angle[i], Pi / 2.0f);
    }
    angle[0] = angle[n - 1] = 0.0f;
    float sum = 0.0f;
    std::vector<float> diff(n - 1);
    for (size_t i = 1; i < n; i++) {
        diff[i - 1] = dist[i] * (1.0f + 1.5f * (angle[i - 1] * dist[i - 1]) / (dist[i - 1] + dist[i]) +
            1.5f * (angle[i] * dist[i + 1]) / (dist[i] + dist[i + 1]));
        sum += diff[i - 1];
    }
    std::vector<float> result(n);
    result[0] = 0.0f;
    for (size_t i = 1; i < n - 1; i++) {
        result[i] = diff[i - 1] / sum;
        result[i] += result[i - 1];
    }
    result[n - 1] = 1.0f;
    return result;
}