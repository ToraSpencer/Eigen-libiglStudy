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
 
ʵ������ϸ�����ߵ����ɷ���
    - �ƽ���ϸ�֣�Chaiukin ���������� B ������������ B ����ϸ�ַ���
    - ��ֵ��ϸ�֣�4 ��ϸ�ַ���


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