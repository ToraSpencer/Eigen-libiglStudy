#pragma once

#include <vector>
#include "Eigen/Dense"


// 计算各种类型拟合的系数：
class MathUtil 
{
public:
    constexpr static const float Pi = 3.141592653589793f;


    /////////////////////////////////////////////////////////////////////////////////////// 插值拟合方法：

    // 多项式插值，重载1
    static std::vector<Eigen::Vector2f> InterpolationPolygon(\
        const std::vector<Eigen::Vector2f> &in_pos, float lb, float rb, float step);

    // 多项式插值，重载2
    static std::vector<float> InterpolationPolygon(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y,
        float lb, float rb, float step);

    // 高斯插值
    static std::vector<Eigen::Vector2f> InterpolationGauss(\
        const std::vector<Eigen::Vector2f> &in_pos, float sigma2, \
        int m, float lb, float rb, float step);

    // 三次B样条插值：
    static std::vector<float> CubicSpline(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y,
        float lb, float rb, float step);




    /////////////////////////////////////////////////////////////////////////////////////// 逼近拟合方法：
    
    // 多项式逼近
    static std::vector<float> ApproximationPolygon(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y,
        int m, float lb, float rb, float step);


    // 幂基函数最小二乘逼近
    static std::vector<Eigen::Vector2f> ApproximationPolygon(\
        const std::vector<Eigen::Vector2f>& in_pos, int m, float lb, float rb, float step);


    // 幂基函数岭回归逼近
    static std::vector<Eigen::Vector2f> ApproximationNormalized(\
        const std::vector<Eigen::Vector2f> &in_pos, int m, float lambda, \
        float lb, float rb, float step);


    /////////////////////////////////////////////////////////////////////////////////////// 点列参数化方法：

    // 点列均匀(Uniform)参数化——参数值均匀分配给控制点
    static std::vector<float> ParameterizationUniform(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y);

    // 点列Choral参数化：
    static std::vector<float> ParameterizationChoral(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y);

    // 点列Centripetal参数化：
    static std::vector<float> ParameterizationCentripetal(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y);

    // 点列Foley参数化：
    static std::vector<float> ParameterizationFoley(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y);
};


