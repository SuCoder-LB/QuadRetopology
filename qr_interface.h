//
// Created by 苏立彪 on 2023/5/9.
//

#ifndef QUADRETOPOLOGY__QR_INTERFACE_H_
#define QUADRETOPOLOGY__QR_INTERFACE_H_

#include <vector>

namespace quad_retopology {

struct ChartData;
class PolyMesh;
class Projector;

// 输入分区信息，与投影对象，没有投影对象则仅根据输入分区边的信息的光顺结果进行网格化，
// 投影类是一个虚类，qr_projector中给出一个实现样例
void QuadRangulate(const ChartData &chart_data, PolyMesh &quad_mesh, Projector *projector = nullptr);

}

#endif //QUADRETOPOLOGY__QR_INTERFACE_H_
