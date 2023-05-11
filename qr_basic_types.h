//
// Created by 苏立彪 on 2023/5/9.
//

#ifndef QUADRETOPOLOGY__QR_BASIC_TYPES_H_
#define QUADRETOPOLOGY__QR_BASIC_TYPES_H_

#include <unordered_map>
#include <vector>
#include <array>


#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/clean.h>

namespace quad_retopology {

struct ChartNode {
  double x;
  double y;
  double z;
  size_t id;
  ChartNode() { x = y = z = 0, id = -1; }
};

struct ChartSubside {
  std::vector<ChartNode> vertices;  //分区线的离散点
  size_t start;                     //起始点编号
  size_t end;                       //终止点编号
  ChartSubside() { start = end = -1; }
};

struct ChartSide {
  std::vector<size_t> subsides;          //相关的子边编号
  std::vector<bool> reversed_subside;     //该子边在这条边上是否逆转
};

struct Chart {
  std::vector<ChartSide> sides;     //分区边结构
};

struct ChartData {
  std::vector<ChartNode> nodes;      //所有分区节点信息
  std::vector<ChartSubside> subsides;     //所有分区子边
  std::vector<Chart> charts;              //所有分区信息
};

class PolyVertex;
class PolyFace;
class PolyEdge;

struct MyPolyTypes : public vcg::UsedTypes<
    vcg::Use<PolyVertex>::AsVertexType,
    vcg::Use<PolyEdge>::AsEdgeType,
    vcg::Use<PolyFace>::AsFaceType> {
};

class PolyVertex : public vcg::Vertex<MyPolyTypes,
                                      vcg::vertex::Coord3d,
                                      vcg::vertex::Normal3f,
                                      vcg::vertex::Color4b,
                                      vcg::vertex::Qualityf,
                                      vcg::vertex::BitFlags,
                                      vcg::vertex::VFAdj,
                                      vcg::vertex::CurvatureDirf> {
};

class PolyFace : public vcg::Face<
    MyPolyTypes,
    vcg::face::PolyInfo,
    vcg::face::VertexRef,
    vcg::face::Normal3f,
    vcg::face::Color4b,
    vcg::face::Qualityf,
    vcg::face::BitFlags,
    vcg::face::PFVAdj,
    vcg::face::PFFAdj,
    vcg::face::PVFAdj,
    vcg::face::CurvatureDirf,
    vcg::face::Mark> {
};

class PolyEdge : public vcg::Edge<
    MyPolyTypes,
    vcg::edge::VertexRef,
    vcg::edge::BitFlags> {
};

class PolyMesh : public vcg::tri::TriMesh<
    std::vector<PolyVertex>,
    std::vector<PolyEdge>,
    std::vector<PolyFace>> {
};

}

#endif //QUADRETOPOLOGY__QR_BASIC_TYPES_H_
