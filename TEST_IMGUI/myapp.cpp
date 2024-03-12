#include "TEST_IMGUI/myapp.h"

static void glfw_error_callback(int error, const char* description) 
{
    std::cerr << "Glfw error"+std::to_string(error)+
                 ": "+std::string(description) << std::endl;
}


void MyApp::initAll() 
{ 
    glfwInit();

    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only

    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); 

    window = glfwCreateWindow(1200, 1200, "Hw4 by dongmo", NULL, NULL);

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    gladLoadGL();

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO &io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);
}


void MyApp::cleanUp() 
{
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
}


void MyApp::mainLoop() 
{
    // 预定变量
    ControlPointArray2D arr;

    // control para
    bool changeTan = false;
    int moveNodeNum = -1;

    while (!glfwWindowShouldClose(window)) 
    {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        {            
            ImGui::SetNextWindowSize(ImVec2(1000,1000), ImGuiCond_Appearing);
            ImGui::Begin("main");
            bool cal = false;
            if (ImPlot::BeginPlot("", ImVec2(900, 900))) 
            {
                if (changeTan) 
                {
                    if (moveNodeNum==-1 && ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(0)) 
                    {
                        ImPlotPoint pt = ImPlot::GetPlotMousePos();
                        moveNodeNum = arr.findSuitableCtrlPoint(pt.x, pt.y);
                    }
                    else if (moveNodeNum!=-1 && ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(0))  
                        moveNodeNum = -1; 
                    if (moveNodeNum != -1) 
                    {
                        ImPlotPoint pt = ImPlot::GetPlotMousePos();
                        if (moveNodeNum < arr.nodenum())
                            arr.setPoint(moveNodeNum, pt.x, pt.y);
                        else if (moveNodeNum < 2 * arr.nodenum())
                            arr.setLDiff(moveNodeNum - arr.nodenum(), pt.x, pt.y);
                        else
                            arr.setRDiff(moveNodeNum - 2 * arr.nodenum(), pt.x, pt.y);
                    }
                }
                else 
                {
                    if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(0)) 
                    {
                        ImPlotPoint pt = ImPlot::GetPlotMousePos();
                        arr.push_back(pt.x, pt.y);
                    }
                }

                {
                    NodeArr& draw_arr = arr.getDrawPoints();
                    std::vector<double> ctrlxs = arr.getControlPointx();
                    std::vector<double> ctrlys = arr.getControlPointy();
                    auto ctrlbar = arr.getControlBar();
                    ImPlot::PlotLine("", draw_arr.xs.data(), draw_arr.ys.data(), draw_arr.size, 0, sizeof(double));
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, -1.0f);
                    ImPlot::PlotScatter("", ctrlxs.data(), ctrlys.data(), ctrlxs.size(), 0, sizeof(double));
                    //ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, -1.0f);
                    if (changeTan) 
                    {
                        //ImPlot::PlotScatter("", &ctrlxs[moveNodeNum], &ctrlys[moveNodeNum], 1, 0, sizeof(double));
                        for (int i = 0; i < ctrlbar.size(); ++i) 
                        {
                            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, -1.0f);
                            ImPlot::PlotLine("", &ctrlbar[i][0], &ctrlbar[i][1], 3, 0, 2 * sizeof(double));
                        }
                    } 
                }
                
                ImPlot::EndPlot();
            }
            ImGui::Checkbox("edit curve", &changeTan);
            ImGui::SameLine();
            if (ImGui::Button("clear")) 
                arr.clear();
             
            ImGui::End();
        }

        ImGui::Render();
        int _width, _height;
        glfwGetFramebufferSize(window, &_width, &_height);
        glViewport(0, 0, _width, _height);
        glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }
}