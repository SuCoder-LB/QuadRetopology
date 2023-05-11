#include <string>
#include <array>
#include <map>

#include "qr_basic_types.h"
#include "qr_interface.h"
#include "qr_projector.h"

#include <wrap/io_trimesh/import.h>

//n_nodes n_subsides n_patches  //n_nodes个点，n_subsides个分区子边，n_patches个分区
//接下来n_nodes行，每行四个数
//node_id x y z        //点的编号，坐标x y z
//接下来n_subsides组数据，每组数据第1行3个数
//start end mid       //分区子边起始点编号，终止点编号（上面的分区节点编号），该条分区子边内部的点个数
//每组接下来mid行
//x y z                //中间节点依序的坐标x y z
//接下来n_patches组数据，每组数据第1行 n_sides+1 个数据
//n_sides n_subsides_1 ... n_subsides_n
//接下来n_subside1行以此类推
//subside_id reversed  //表示分区子边id，在该分区中，该条子边是否逆置

void LoadPatches(const std::string &filename, quad_retopology::ChartData &chart_data) {
  FILE *fp = fopen(filename.c_str(), "r");
  if (fp == nullptr)return;
  int n_nodes, n_subsides, n_patches;
  fscanf(fp, "%d%d%d\n", &n_nodes, &n_subsides, &n_patches);
  std::map<size_t, int> node_mp;
  for (int i = 0; i < n_nodes; ++i) {
    quad_retopology::ChartNode node;
    node_mp[node.id] = i;
    fscanf(fp, "%zd%lf%lf%lf", &node.id, &node.x, &node.y, &node.z);
    chart_data.nodes.emplace_back(node);
  }

  for (int i = 0; i < n_subsides; ++i) {
    quad_retopology::ChartSubside subside;
    int n_mid;
    fscanf(fp, "%zd%zd%d", &subside.start, &subside.end, &n_mid);
    subside.vertices.push_back(chart_data.nodes[node_mp[subside.start]]);
    for (int j = 0; j < n_mid; ++j) {
      quad_retopology::ChartNode node;
      fscanf(fp, "%lf%lf%lf", &node.x, &node.y, &node.z);
      subside.vertices.emplace_back(node);
    }
    subside.vertices.push_back(chart_data.nodes[node_mp[subside.end]]);
    chart_data.subsides.emplace_back(subside);
  }

  for (int i = 0; i < n_patches; ++i) {
    quad_retopology::Chart chart;
    int n_sides;
    fscanf(fp, "%d", &n_sides);
    chart.sides.resize(n_sides);
    for (int j = 0; j < n_sides; ++j) {
      int n_subsides_j;
      fscanf(fp, "%d", &n_subsides_j);
      chart.sides[j].subsides.resize(n_subsides_j);
      chart.sides[j].reversed_subside.resize(n_subsides_j);
    }
    for (int j = 0; j < n_sides; ++j) {
      for (int k = 0; k < chart.sides[j].subsides.size(); ++k) {
        int subside_id, reversed;
        fscanf(fp, "%d%d", &subside_id, &reversed);
        chart.sides[j].subsides[k] = subside_id;
        chart.sides[j].reversed_subside[k] = reversed;
      }
    }
    chart_data.charts.emplace_back(chart);
  }
  fclose(fp);
}

int main() {

  //读取分区信息
  quad_retopology::ChartData chart_data;
  LoadPatches("../patchInformation.patch", chart_data);

  //读取原背景网格初始化一个投影对象
  bool have_projector = true;
  quad_retopology::SurfaceProjector sp;
  quad_retopology::PolyMesh background_mesh;
  int mask;
  vcg::tri::io::ImporterOBJ<quad_retopology::PolyMesh>::LoadMask("../originTriMesh.obj", mask);
  int err = vcg::tri::io::ImporterOBJ<quad_retopology::PolyMesh>::Open(background_mesh, "../originTriMesh.obj", mask);
  if (err != 0)have_projector = false;
  if (have_projector) sp.Init(background_mesh);

  if (!have_projector) {
    printf("%s:%d(%s),无背景网格投影对象，仅根据分区生成网格\n", __FILE__, __LINE__, __FUNCTION__);
  }

  //输出的四边形网格
  quad_retopology::PolyMesh quad_mesh;
  quad_retopology::QuadRangulate(chart_data, quad_mesh, have_projector ? &sp : nullptr);

  vcg::tri::io::ExporterOBJ<quad_retopology::PolyMesh>::Save(quad_mesh, "../quad.obj", 0);

  return 0;
}
