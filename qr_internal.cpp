//
// Created by 苏立彪 on 2023/5/9.
//
#include "qr_interface.h"
#include "qr_basic_types.h"
#include "qr_projector.h"

#include "patterns/generate_patch.h"

namespace quad_retopology {

void ComputePattern(const std::vector<size_t> &patch_info, PolyMesh &poly_mesh, std::vector<size_t> &borders);

void QuadRangulate(const ChartData &chart_data, PolyMesh &quad_mesh, Projector *projector) {

  std::vector<std::vector<size_t>> subside_vertex_map(chart_data.subsides.size());  //分区子边与四边形点对应
  std::unordered_map<size_t, size_t> corners_map;                                //分区点的id与四边形点对应，考虑到分区点的乱序

  quad_mesh.Clear();

  //分区点对应关系，默认第i个分区点为第i个四边形点
  for (auto node : chart_data.nodes) {
    size_t new_vertex_id = quad_mesh.vert.size();
    vcg::tri::Allocator<PolyMesh>::AddVertex(
        quad_mesh, {node.x, node.y, node.z});
    corners_map[node.id] = new_vertex_id;
  }

  //查看这些子边数据是否正确
  std::vector<bool> valid_side(chart_data.subsides.size(), true);

  //分区边点对应，从分区点开始后推
  for (size_t subside_id = 0; subside_id < chart_data.subsides.size(); ++subside_id) {
    const ChartSubside &subside = chart_data.subsides[subside_id];
    if (corners_map.find(subside.start) == corners_map.end() || corners_map.find(subside.end) == corners_map.end()) {
      printf("%s:%i(%s) (error):输入分区子边的首尾点未找到，分区子边号%zu\n",
             __FILE__, __LINE__, __FUNCTION__, subside_id);
      valid_side[subside_id] = false;
      continue;
    }
    subside_vertex_map[subside_id].emplace_back(corners_map[subside.start]);
    for (size_t k = 1; k < subside.vertices.size() - 1; ++k) {
      size_t new_vertex_id = quad_mesh.vert.size();
      vcg::tri::Allocator<PolyMesh>::AddVertex(
          quad_mesh, {subside.vertices[k].x, subside.vertices[k].y, subside.vertices[k].z});
      subside_vertex_map[subside_id].emplace_back(new_vertex_id);
    }
    subside_vertex_map[subside_id].emplace_back(corners_map[subside.end]);
  }


  //开始计算每个分区的四边形网格
  for (size_t cId = 0; cId < chart_data.charts.size(); cId++) {
    const Chart &chart = chart_data.charts[cId];

    const std::vector<ChartSide> &chart_sides = chart.sides;

    //检查边数
    if (chart_sides.size() < 3 || chart_sides.size() > 6) {
      printf("%s:%i(%s) (error):输入分区边数量错误，分区号%zu\n", __FILE__, __LINE__, __FUNCTION__, cId);
      continue;
    }

    bool pattern_solve = true;

    //检查输入
    std::vector<size_t> patch_info(chart_sides.size());
    for (size_t side_id = 0; side_id < chart_sides.size() && pattern_solve; ++side_id) {
      const auto &chart_side = chart_sides[side_id];
      if (chart_side.subsides.size() != chart_side.reversed_subside.size()) {
        printf("%s:%i(%s) (error):分区边子边数量与相反判定数量不一致，分区号%zu,分区边号%zu\n",
               __FILE__, __LINE__, __FUNCTION__, cId, side_id);
        pattern_solve = false;
        break;
      }
      for (auto &subside_id : chart_side.subsides) {
        if (subside_id < 0 || subside_id >= chart_data.subsides.size()) {
          printf("%s:%i(%s) (error):分区子边号错误，分区号%zu,分区边号%zu,分区子边号%zu\n",
                 __FILE__, __LINE__, __FUNCTION__, cId, side_id, subside_id);
          pattern_solve = false;
          break;
        }
        if (!valid_side[subside_id]) {
          printf("%s:%i(%s) (error):含无效分区边(首尾点不在分区节点中)，分区号%zu,分区线号%zu,分区子边号%zu\n",
                 __FILE__, __LINE__, __FUNCTION__, cId, side_id, subside_id);
          pattern_solve = false;
          break;
        }
        if (chart_data.subsides[subside_id].vertices.size() < 2) {
          printf("%s:%i(%s) (error):分区子边点数量错误，分区号%zu,分区线号%zu,分区子边号%zu\n",
                 __FILE__, __LINE__, __FUNCTION__, cId, side_id, subside_id);
          pattern_solve = false;
          break;
        }
      }
    }
    //检查分区子边的首位是否相连,以及总段数是否符合要求
    if (pattern_solve) {
      size_t segment = 0;
      size_t pre_node = chart_sides.back().reversed_subside.back() ?
                        chart_data.subsides[chart_sides.back().subsides.back()].start :
                        chart_data.subsides[chart_sides.back().subsides.back()].end;
      for (size_t side_id = 0; side_id < chart_sides.size() && pattern_solve; ++side_id) {
        const auto &chart_side = chart_sides[side_id];
        for (size_t chart_subside_id = 0; chart_subside_id < chart_side.subsides.size(); ++chart_subside_id) {
          size_t cur_node = chart_side.reversed_subside[chart_subside_id] ?
                            chart_data.subsides[chart_side.subsides[chart_subside_id]].end :
                            chart_data.subsides[chart_side.subsides[chart_subside_id]].start;
          if (cur_node != pre_node) {
            printf("%s:%i(%s) (error):分区子边未连接上，分区号%zu,分区线号%zu,分区子边号%zu\n",
                   __FILE__, __LINE__, __FUNCTION__, cId, side_id, chart_subside_id);
            pattern_solve = false;
            break;
          }
          pre_node = chart_side.reversed_subside[chart_subside_id] ?
                     chart_data.subsides[chart_side.subsides[chart_subside_id]].start :
                     chart_data.subsides[chart_side.subsides[chart_subside_id]].end;
          segment += chart_data.subsides[chart_side.subsides[chart_subside_id]].vertices.size() - 1;
          patch_info[side_id] += chart_data.subsides[chart_side.subsides[chart_subside_id]].vertices.size() - 1;
        }
      }
      if (pattern_solve && (segment < 4 || (segment & 1))) {
        printf("%s:%i(%s) (error):分区子边总段数错误%zu，分区号%zu\n",
               __FILE__, __LINE__, __FUNCTION__, segment, cId);
        pattern_solve = false;
      }
    }

    if (!pattern_solve) {
      continue;
    }


    //计算参数域内的分区四边形网格存储在patch_mesh中
    std::vector<size_t> patch_borders;
    PolyMesh patch_mesh;
    ComputePattern(patch_info, patch_mesh, patch_borders);


    //记录patch网格与最终网格的点关系对应
    std::vector<size_t> current_vertex_map(patch_mesh.vert.size(), -1);

    //记录边界点与四边形网格点的对应关系
    size_t current_patch_side_vertex = 0;
    for (const auto &side : chart_sides) {

      for (size_t j = 0; j < side.subsides.size(); j++) {
        const size_t &subsideId = side.subsides[j];
        const bool &reversed = side.reversed_subside[j];
        const ChartSubside &subside = chart_data.subsides[subsideId];

        for (int k = 0; k < subside.vertices.size() - 1; k++) {
          size_t side_v_id = patch_borders[current_patch_side_vertex];
          size_t vertex_index = reversed ? subside.vertices.size() - 1 - k : k;
          current_vertex_map[side_v_id] = subside_vertex_map[subsideId][vertex_index];

//          printf("%s:%i(%s) (debug):边界点%zu,原坐标%lf,%lf,%lf,后坐标%lf,%lf,%lf\n",
//                 __FILE__, __LINE__, __FUNCTION__, currentPatchSideVertex,
//                 patch_mesh.vert[patchSideVId].P()[0], patch_mesh.vert[patchSideVId].P()[1],
//                 patch_mesh.vert[patchSideVId].P()[2],
//                 quad_mesh.vert[currentVertexMap[patchSideVId]].P()[0],
//                 quad_mesh.vert[currentVertexMap[patchSideVId]].P()[1],
//                 quad_mesh.vert[currentVertexMap[patchSideVId]].P()[2]);

          patch_mesh.vert[side_v_id].P() = quad_mesh.vert[current_vertex_map[side_v_id]].P();

          current_patch_side_vertex++;
        }
      }
    }

    // 固定边界进行一次laplace光顺，直到光顺不起作用，并将全局光顺后的点根据投影对象进行一次投影,没有投影对象就默认使用光顺后的作为背景
    //todo : 这里理应使用全局laplace光顺，vcg库里没找到，不想再额外实现
    // 先固定100次局部laplace光顺凑合着用（适当调整，10次不够），之后再进行5次带权重的局部光顺
    {
      vcg::tri::UpdateSelection<PolyMesh>::VertexClear(patch_mesh);
      for (size_t vId : patch_borders) patch_mesh.vert[vId].SetS();

      vcg::tri::UpdateTopology<PolyMesh>::FaceFace(patch_mesh);

      vcg::PolygonalAlgorithm<PolyMesh>::Laplacian(patch_mesh, true, 100, 0);
      if (projector) {
        for (size_t v_id = 0; v_id < patch_mesh.vn; ++v_id) {
          if (patch_mesh.vert[v_id].IsS() || patch_mesh.vert[v_id].IsD())continue;
          //todo:投影对象暂时不要求法向信息，可能会投影错位置
          auto new_pos = projector->Project({patch_mesh.vert[v_id].P()[0],
                                             patch_mesh.vert[v_id].P()[1],
                                             patch_mesh.vert[v_id].P()[2]}, cId);
          patch_mesh.vert[v_id].P() = {new_pos[0], new_pos[1], new_pos[2]};
        }
      }

      for (int i = 0; i < 5; ++i) {
        vcg::PolygonalAlgorithm<PolyMesh>::Laplacian(patch_mesh, true, 1, 0.5);
        if (projector) {
          for (size_t v_id = 0; v_id < patch_mesh.vn; ++v_id) {
            if (patch_mesh.vert[v_id].IsS() || patch_mesh.vert[v_id].IsD())continue;
            //todo:同上
            auto new_pos = projector->Project({patch_mesh.vert[v_id].P()[0],
                                               patch_mesh.vert[v_id].P()[1],
                                               patch_mesh.vert[v_id].P()[2]}, cId);
            patch_mesh.vert[v_id].P() = {new_pos[0], new_pos[1], new_pos[2]};
          }
        }
      }
    }

    //内部点增加与关系对应
    for (size_t i = 0; i < patch_mesh.vert.size(); i++) {
      if (current_vertex_map[i] == -1) {
        size_t new_id = quad_mesh.vert.size();
        const typename PolyMesh::CoordType &coord = patch_mesh.vert[i].P();
        vcg::tri::Allocator<PolyMesh>::AddVertex(quad_mesh, coord);
        current_vertex_map[i] = new_id;
      }
    }

    //增加单元面片
    for (size_t i = 0; i < patch_mesh.face.size(); i++) {
      size_t new_face_id = quad_mesh.face.size();
      vcg::tri::Allocator<PolyMesh>::AddFaces(quad_mesh, 1);
      quad_mesh.face[new_face_id].Alloc(4);
      for (int j = 0; j < 4; j++) {
        size_t vId = current_vertex_map[vcg::tri::Index(patch_mesh, patch_mesh.face[i].V(j))];
        quad_mesh.face[new_face_id].V(j) = &quad_mesh.vert[vId];
      }
    }
  }
}

void ComputePattern(const std::vector<size_t> &patch_info, PolyMesh &poly_mesh, std::vector<size_t> &borders) {

  patchgen::PatchParam param;
  patterns::Patch patch;
  Eigen::VectorXi l(patch_info.size());
  for (int i = 0; i < patch_info.size(); ++i) l(i) = static_cast<int>(patch_info[i]);
  patterns::generatePatch(l, param, patch);

  std::vector<size_t> corners(l.size());

  vcg::tri::Allocator<PolyMesh>::AddVertices(poly_mesh, patch.n_vertices());
  for (patterns::Patch::VertexIter v_it = patch.vertices_begin(); v_it != patch.vertices_end(); ++v_it) {
    const int &corner_index = patch.data(*v_it).patchgen.corner_index;
    if (corner_index >= 0) corners[corner_index] = v_it->idx();
    patterns::Patch::Point p = patch.point(*v_it);
    poly_mesh.vert[v_it->idx()].P() = {p[0], p[1], p[2]};
  }
  vcg::tri::Allocator<PolyMesh>::AddFaces(poly_mesh, patch.n_faces());
  int face_id = 0;
  for (patterns::Patch::FaceIter f_it = patch.faces_begin(); f_it != patch.faces_end(); ++f_it) {
    poly_mesh.face[face_id].Alloc(4);
    int j = 0;
    for (patterns::Patch::FaceVertexIter fv_it = patch.fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
      if (j >= 4)break;
      poly_mesh.face[face_id].V(j++) = &(poly_mesh.vert[fv_it->idx()]);
    }
    ++face_id;
  }

  vcg::tri::UpdateTopology<PolyMesh>::FaceFace(poly_mesh);

  //todo :找循环边界可以根据上面得到的patch的半边信息去找，可以省一些空间
  std::vector<size_t> next_map(poly_mesh.vn, poly_mesh.vn + 1);
  for (size_t i = 0; i < poly_mesh.face.size(); i++) {
    for (int j = 0; j < 4; j++) {
      if (vcg::face::IsBorder(poly_mesh.face[i], j)) {
        size_t start_id = vcg::tri::Index(poly_mesh, poly_mesh.face[i].V0(j));
        size_t end_id = vcg::tri::Index(poly_mesh, poly_mesh.face[i].V1(j));
        next_map[start_id] = end_id;
      }
    }
  }

  size_t current_id = corners[0];
  do {
    borders.push_back(current_id);
    current_id = next_map[current_id];
  } while (current_id != corners[0]);
}

}