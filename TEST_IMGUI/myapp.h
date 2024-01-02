#pragma once
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "glad/glad.h"
#include <GLFW/glfw3.h>
#include "implot.h"
#include "implot_internal.h"

#include "TEST_IMGUI/fitcurve.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>

// hw5
/*
 
实现两种细分曲线的生成方法
    - 逼近型细分：Chaiukin 方法（二次 B 样条），三次 B 样条细分方法
    - 插值型细分：4 点细分方法


*/
 

class MyApp
{
public:

    void initAll();
    void mainLoop();
    void cleanUp();


private:
    GLFWwindow *window;

};