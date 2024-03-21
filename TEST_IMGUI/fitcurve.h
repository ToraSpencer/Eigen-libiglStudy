#pragma once
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
 

struct NodeArr 
{
    NodeArr() :size(0){}
    NodeArr(int _size) : size(_size), xs(std::vector<double>(_size)), ys(std::vector<double>(_size)) {}
    NodeArr(std::vector<double> x, std::vector<double> y) : xs(x), ys(y), size(x.size()) {}
    unsigned int size;
    std::vector<double> xs;
    std::vector<double> ys;
};


struct controlPoint 
{
    double val = 0;
    double ldiff = 0;
    double rdiff = 0;
    bool fixed_diff = false;
};


std::vector<double> ThreeOrderSample(const std::vector<controlPoint>& ctrlarr, int begin, int end);


class ControlPointArray2D 
{
// 成员数据：
public:
    std::vector<controlPoint> xs;
    std::vector<controlPoint> ys;
    std::vector<bool> fixed;

    // 绘图的点
    std::vector<NodeArr> drawPoints;
    std::vector<std::vector<double>> param;
    uint64_t nodePerRange = 100;

    //
    NodeArr drawarr;

    // ctrlpoints
    std::vector<double> ctrlxs;
    std::vector<double> ctrlys;
    double maxChoosedist = 0.02;

    double showpara = 0.4;

    std::vector<double> ctrl_l_xs;
    std::vector<double> ctrl_l_ys;
    std::vector<double> ctrl_r_xs;
    std::vector<double> ctrl_r_ys;
    std::vector<std::vector<double>> controlbars;

public:
    void calculateParam_extend(int);


    int32_t size() 
    {
        return param.size();
    }


    uint32_t nodenum()
    {
        assert(xs.size()==ys.size() and ys.size()==fixed.size());
        return fixed.size();
    }


    void calculateRange(int p);


    void clear() 
    {
        while (nodenum() != 0)
        {
            xs.clear();
            ys.clear();
            fixed.clear();
            drawPoints.clear();
            param.clear();
            drawarr.xs.clear();
            drawarr.ys.clear();
            drawarr.size = 0;
            ctrlxs.clear();
            ctrlys.clear();
            ctrl_l_xs.clear();
            ctrl_l_ys.clear();
            ctrl_r_xs.clear();
            ctrl_r_ys.clear();
            controlbars.clear();
        }
    }


    void push_back(double x, double y);


    void delete_at(int pos);


    void setPoint(int p, double x, double y)
    {
        xs[p].val = x;
        ys[p].val = y;
        if (p > 0) {
            calculateParam_extend(p - 1);
        }
        else if (p == 0) {
            calculateParam_extend(0);
        }
        if (p < nodenum() - 1) {
            calculateParam_extend(p);
        }
    }


    void setLDiff(int p, double x, double y) {
        fixed[p] = true;
        xs[p].fixed_diff = true;
        ys[p].fixed_diff = true;
        xs[p].ldiff = (xs[p].val - x)/showpara;
        ys[p].ldiff = (ys[p].val - y)/showpara;
        if (p > 0) {
            calculateParam_extend(p - 1);
        }
    }


    void setRDiff(int p, double x, double y) {
        fixed[p] = true;
        xs[p].fixed_diff = true;
        ys[p].fixed_diff = true;
        xs[p].rdiff = (x - xs[p].val)/showpara;
        ys[p].rdiff = (y - ys[p].val)/showpara;
        if (p < fixed.size() - 1) {
            calculateParam_extend(p);
        }
    }


    NodeArr& getDrawPoints() {
        return drawarr;
    }


    std::vector<double>& getControlPointx() {
        return ctrlxs;
    }


    std::vector<double>& getControlPointy() {
        return ctrlys;
    }


    std::vector<std::vector<double>>& getControlBar() { return controlbars; }


    int findSuitableCtrlPoint(double x, double y);
     

    void setFixDiff(int pos) {
        fixed.at(pos) = true;
        xs.at(pos).fixed_diff = true;
        ys.at(pos).fixed_diff = true;
    }
     

    void update_drawarr()
    {
        drawarr.xs.clear();
        drawarr.ys.clear();
        for (int i = 0; i < size(); ++i) {
            drawarr.xs.insert(drawarr.xs.end(), drawPoints[i].xs.begin(), drawPoints[i].xs.end());
            drawarr.ys.insert(drawarr.ys.end(), drawPoints[i].ys.begin(), drawPoints[i].ys.end());
        }
        drawarr.size = drawarr.xs.size();
    }


    void updateCtrlPoints()
    {
        ctrlxs.resize(xs.size());
        ctrlys.resize(xs.size());
        ctrl_l_xs.resize(xs.size());
        ctrl_l_ys.resize(xs.size());
        ctrl_r_xs.resize(xs.size());
        ctrl_r_ys.resize(xs.size());
        controlbars.resize(xs.size());
        for (int i = 0; i < xs.size(); ++i) 
        {
            ctrlxs[i] = xs[i].val;
            ctrlys[i] = ys[i].val;
            ctrl_l_xs[i] = xs[i].val - xs[i].ldiff * showpara;
            ctrl_l_ys[i] = ys[i].val - ys[i].ldiff * showpara;
            ctrl_r_xs[i] = xs[i].val + xs[i].rdiff * showpara;
            ctrl_r_ys[i] = ys[i].val + ys[i].rdiff * showpara;
            controlbars[i] = { ctrl_l_xs[i], ctrl_l_ys[i], ctrlxs[i], ctrlys[i], ctrl_r_xs[i], ctrl_r_ys[i] };
        }
        
    }
 
};

