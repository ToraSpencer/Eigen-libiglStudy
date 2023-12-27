#pragma once

#include <vector>
#include "Eigen/Dense"


// �������������ϵ�ϵ����
class MathUtil 
{
public:
    constexpr static const float Pi = 3.141592653589793f;

    // ����ʽ��ֵ
    static std::vector<Eigen::Vector2f> InterpolationPolygon(\
        const std::vector<Eigen::Vector2f> &in_pos, float lb, float rb, float step);

    static std::vector<float> InterpolationPolygon(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y,
        float lb, float rb, float step);

    // ��˹��ֵ
    static std::vector<Eigen::Vector2f> InterpolationGauss(\
        const std::vector<Eigen::Vector2f> &in_pos, float sigma2, \
        int m, float lb, float rb, float step);

    // �ݻ�������С���˱ƽ�
    static std::vector<Eigen::Vector2f> ApproximationPolygon(\
        const std::vector<Eigen::Vector2f> &in_pos, int m, float lb, float rb, float step);

    static std::vector<float> ApproximationPolygon(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y,
        int m, float lb, float rb, float step);

    // �ݻ�������ع�ƽ�
    static std::vector<Eigen::Vector2f> ApproximationNormalized(\
        const std::vector<Eigen::Vector2f> &in_pos, int m, float lambda, \
        float lb, float rb, float step);

    static std::vector<float> CubicSpline(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y,
        float lb, float rb, float step);


    static std::vector<float> ParameterizationUniform(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y);

    static std::vector<float> ParameterizationChoral(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y);

    static std::vector<float> ParameterizationCentripetal(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y);

    static std::vector<float> ParameterizationFoley(\
        const std::vector<float>& pos_x, const std::vector<float>& pos_y);
};

