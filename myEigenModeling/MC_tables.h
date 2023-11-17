#pragma once

namespace MC_TABLES
{
    // 相交边编码——立方体中12条边和三角片相交，有256中情形；
    extern const int edgeStateCodes[256];

    // 立方体的12条边，用端点索引对表示：
    extern const int cubeEdges[12][2];

    // 立方体中新生成的三角片数据，最多5个三角片；
    extern const int cubeTriangles[256][16];
}