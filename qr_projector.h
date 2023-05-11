//
// Created by 苏立彪 on 2023/5/9.
//

#ifndef QUADRETOPOLOGY__QR_PROJECTOR_H_
#define QUADRETOPOLOGY__QR_PROJECTOR_H_

#include <vector>
#include <array>


namespace quad_retopology {

class Projector {
 public:
  virtual std::array<double, 3> Project(std::array<double, 3> query, size_t patch_id) = 0;
};
}

#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>

#include "qr_basic_types.h"

namespace quad_retopology {
class SurfaceProjector : public Projector {

 private:

  typedef typename TempMesh::FaceType TriFaceType;
  typedef typename TempMesh::ScalarType TriScalarType;
  typedef typename TempMesh::CoordType TriCoordType;
  typedef vcg::GridStaticPtr<TriFaceType, TriScalarType> TriMeshGrid;

  TriMeshGrid grid_;
  TempMesh tri_mesh_;
  TriScalarType max_d_;

 public:
  void Init(PolyMesh &tri_mesh) {
    //vcg::tri::PolygonSupport<TempMesh,PolyMeshType>:(GuideSurf,poly_m);
    vcg::PolygonalAlgorithm<PolyMesh>::TriangulateToTriMesh(tri_mesh,tri_mesh_);
    vcg::tri::UpdateBounding<TempMesh>::Box(tri_mesh_);
    vcg::tri::UpdateNormal<TempMesh>::PerVertexNormalizedPerFace(tri_mesh_);
    vcg::tri::UpdateTopology<TempMesh>::FaceFace(tri_mesh_);
    vcg::tri::UpdateFlags<TempMesh>::FaceBorderFromFF(tri_mesh_);
    vcg::tri::MeshAssert<TempMesh>::VertexNormalNormalized(tri_mesh_);
    grid_.Set(tri_mesh_.face.begin(), tri_mesh_.face.end());
    max_d_ = tri_mesh_.bbox.Diag();
  }

  std::array<double, 3> Project(std::array<double, 3> query, size_t patch_id) override {
    TriCoordType test_pos = {query[0], query[1], query[2]};
    TriCoordType closest_pt;
    TriScalarType min_dist;
    TriCoordType norm, ip;
    vcg::tri::GetClosestFaceBase(tri_mesh_, grid_, test_pos, max_d_, min_dist, closest_pt, norm, ip);

    return {closest_pt[0], closest_pt[1], closest_pt[2]};
  }
};

}

#endif //QUADRETOPOLOGY__QR_PROJECTOR_H_
