#include "test_imgui.h"
#include "TEST_IMGUI/MathUtil.h"

#include <iostream>
#include <vector>
#include <algorithm>

#include "glad/glad.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imgui_internal.h"
#include "GLFW/glfw3.h"

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


namespace MY_IMGUI 
{
    GLFWwindow* g_window = nullptr;                         // 全局的窗口对象
    int g_width = 1600;
    int g_height = 1200;

    ImVec2 g_canvas_pos_ul = { 0.0f, 0.0f };
    ImVec2 g_canvas_pos_br = { 0.0f, 0.0f };

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


    void PlotLineSegments(const std::vector<float>& pos_x, const std::vector<float>& pos_y, ImDrawList* draw_list,
        ImU32 line_col, ImU32 point_col, ImVec2 ul = g_canvas_pos_ul, ImVec2 br = g_canvas_pos_br) {
        const size_t n = pos_x.size();
        for (size_t i = 1; i < n; i++) {
            draw_list->AddLine({ ul.x + pos_x[i - 1], br.y - pos_y[i - 1] },
                { ul.x + pos_x[i], br.y - pos_y[i] }, line_col, 2.0f);
        }
        for (size_t i = 0; i < n; i++) {
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
    std::vector<float> pos_x;
    std::vector<float> pos_y;
    std::vector<float> in_pos_t;
    bool visible = false;

    void Clear() 
    {
        pos_x.clear();
        pos_y.clear();
        in_pos_t.clear();
    }

    void Plot(ImDrawList* draw_list) const 
    {
        if (visible && !pos_x.empty()) 
        {
            PlotLineSegments(pos_x, pos_y, draw_list, line_color, point_color);
        }
    }

    void PlotXT(const std::vector<float>& pos_t, ImDrawList* draw_list, ImVec2 canvas_ul, ImVec2 canvas_br) const {
        const size_t n = pos_x.size();
        if (!visible || n == 0) {
            return;
        }
        std::vector<float> pos_x_t = pos_x;
        std::vector<float> pos_t_t = pos_t;
        for (size_t i = 0; i < n; i++) {
            pos_x_t[i] = pos_x_t[i] / (g_canvas_pos_br.x - g_canvas_pos_ul.x) * (canvas_br.y - canvas_ul.y);
            pos_t_t[i] = pos_t_t[i] * (canvas_br.x - canvas_ul.x);
        }
        PlotLineSegments(pos_t_t, pos_x_t, draw_list, line_color, point_color, canvas_ul, canvas_br);
    }

    void PlotYT(const std::vector<float>& pos_t, ImDrawList* draw_list, ImVec2 canvas_ul, ImVec2 canvas_br) const {
        const size_t n = pos_y.size();
        if (!visible || n == 0) {
            return;
        }
        std::vector<float> pos_y_t = pos_y;
        std::vector<float> pos_t_t = pos_t;
        for (size_t i = 0; i < n; i++) {
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
        if (!Initialize())  
            return; 

        std::vector<float> in_pos_x;
        std::vector<float> in_pos_y;

        bool inter_update = false;
        CurveData inter_uniform{ IM_COL32(255, 50, 50, 255), IM_COL32(255, 80, 80, 255) };
        CurveData inter_chordal{ IM_COL32(50, 255, 50, 255), IM_COL32(80, 255, 80, 255) };
        CurveData inter_centripetal{ IM_COL32(50, 50, 255, 255), IM_COL32(80, 80, 255, 255) };
        CurveData inter_foley{ IM_COL32(150, 150, 255, 255), IM_COL32(180, 180, 255, 255) };

        bool approx_update = false;
        int approx_m = 0;
        int approx_m_temp = 0;
        CurveData approx_uniform{ IM_COL32(50, 255, 255, 255), IM_COL32(80, 255, 255, 255) };
        CurveData approx_chordal{ IM_COL32(255, 50, 255, 255), IM_COL32(255, 80, 255, 255) };
        CurveData approx_centripetal{ IM_COL32(255, 255, 50, 255), IM_COL32(255, 255, 80, 255) };
        CurveData approx_foley{ IM_COL32(255, 150, 150, 255), IM_COL32(255, 180, 180, 255) };

        bool spline_update = false;
        CurveData spline_uniform{ IM_COL32(255, 100, 50, 255), IM_COL32(255, 130, 80, 255) };
        CurveData spline_chordal{ IM_COL32(50, 255, 100, 255), IM_COL32(80, 255, 130, 255) };
        CurveData spline_centripetal{ IM_COL32(100, 50, 255, 255), IM_COL32(130, 80, 255, 255) };
        CurveData spline_foley{ IM_COL32(170, 255, 170, 255), IM_COL32(200, 255, 200, 255) };

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

        while (!glfwWindowShouldClose(g_window)) 
        {
            BeginFrame();
            if (glfwGetKey(g_window, GLFW_KEY_ESCAPE) == GLFW_PRESS) 
                break;
            

            if (ImGui::Begin("Config"))
            {
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

                ImGui::Separator();
                ImGui::Text("Fitting");
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
                }
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

                ImGui::InvisibleButton("canvas", canvas_size);
                const bool is_hovered = ImGui::IsItemHovered();
                if (is_hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
                {
                    in_pos_x.push_back(io.MousePos.x - g_canvas_pos_ul.x);
                    in_pos_y.push_back(g_canvas_pos_br.y - io.MousePos.y);
                    inter_update = true;
                    approx_update = true;
                    spline_update = true;
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
                    }
                }

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

                PlotLineSegments(in_pos_x, in_pos_y, draw_list, IM_COL32(255, 255, 255, 0), IM_COL32(255, 255, 255, 255));
                ImGui::End();
            }

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


    void test2()
    {
        MyApp app;
        app.initAll();
        app.mainLoop();
        app.cleanUp();

    }
}
