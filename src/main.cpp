#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/quadric_error_simplification.h"
#include "geometrycentral/surface/rich_surface_mesh_data.h"
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/surface/boundary_first_flattening.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/transfer_functions.h"
#include "geometrycentral/utilities/utilities.h"

#include "happly.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

#include "args/args.hxx"
#include "imgui.h"

#include <Eigen/Dense>
#include <fmt/base.h>
#include <fmt/format.h>

#include "int_coarsen.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<SignpostIntrinsicTriangulation> intrinsic;

std::unique_ptr<ManifoldSurfaceMesh> n_mesh;
std::unique_ptr<VertexPositionGeometry> n_geometry;
std::unique_ptr<SignpostIntrinsicTriangulation> n_intrinsic;

std::unique_ptr<RichSurfaceMeshData> rich_data;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;
polyscope::CurveNetwork *psCurves;
polyscope::SurfaceMesh *psExtMesh;
polyscope::SurfaceMesh *psIntMesh;

// Structs
struct MeshIteration {
  int iter;
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  bool operator>(const MeshIteration& other) const { return iter > other.iter; }
};

// Some algorithm parameters
float param1 = 1.0;
double tolerance = 1.0;
int target;

// Meta information
int qed_iter_count = 0;
bool recordOutput = true;
bool saveGeoDist = false;
bool useIDT = false;
std::string outPath = "MeshAnalysis/QEDIterations/";
SurfacePoint geoDistSrcPt; 
SurfacePoint flattenTarget;
std::string meshName;
double minAngleDeg;
double inf = std::numeric_limits<double>::max();

std::priority_queue<MeshIteration, std::vector<MeshIteration>, std::greater<MeshIteration>> qemIters;

// Save the current resolution mesh information for use later
void saveMeshRes(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo) {
  qemIters.push(MeshIteration{qed_iter_count, mesh.copy(), geo.copy()});
}

void visualizeIntrinsicTriangulation() {
  
  EdgeData<std::vector<SurfacePoint>> edgePaths = intrinsic->traceAllIntrinsicEdgesAlongInput();

  std::vector<glm::vec3> nodes;
  std::vector<std::array<size_t, 2>> edges;

  for (Edge e : intrinsic->mesh.edges()) {
      const std::vector<SurfacePoint>& path = edgePaths[e];
      if (path.empty()) continue;

      size_t startNodeIdx = nodes.size();
      
      for (const SurfacePoint& p : path) {
          Vector3 pos3D = p.interpolate(geometry->vertexPositions);
          nodes.push_back(glm::vec3(pos3D.x, pos3D.y, pos3D.z));
      }
      
      for (size_t i = 0; i < path.size() - 1; i++) {
          edges.push_back({startNodeIdx + i, startNodeIdx + i + 1});
      }
  }
  psCurves = polyscope::registerCurveNetwork("idt", nodes, edges);
  psCurves->setRadius(0.002);
  psCurves->setColor(glm::vec3(1.0, 0.8, 0.0));
  psCurves->setEnabled(true);
}

// Select a point on the mesh
void selectPoint(SurfacePoint& pt, std::string pointCloudVisName) {
  if (polyscope::hasPointCloud(pointCloudVisName)) {
    polyscope::getPointCloud(pointCloudVisName)->setEnabled(false);
  }
  if (polyscope::hasCurveNetwork(pointCloudVisName)) {
    polyscope::getCurveNetwork(pointCloudVisName)->setEnabled(false);
  }

  long long int pick = psMesh->selectVertex();
  if (pick >= 0) {
    pt = mesh->vertex(pick);
  }

  polyscope::registerPointCloud(pointCloudVisName, std::vector<Vector3>{pt.interpolate(geometry->inputVertexPositions)})
  ->setPointRadius(0.005)
  ->setPointColor(glm::vec3{255, 128, 0});

  if (polyscope::hasPointCloud(pointCloudVisName)) {
    polyscope::getPointCloud(pointCloudVisName)->setEnabled(true);
  }
  if (polyscope::hasCurveNetwork(pointCloudVisName)) {
    polyscope::getCurveNetwork(pointCloudVisName)->setEnabled(true);
  }
}



// Compute Geodesic Distance via Heat Method from source point
void computeGeoDistFromSrc() {
  if (geoDistSrcPt == SurfacePoint()) {
    return;
  }

  intrinsic->flipToDelaunay();
  VertexData<double> dist;
  if (useIDT) {
    HeatMethodDistanceSolver heatSolver(*intrinsic);
    dist = heatSolver.computeDistance(geoDistSrcPt);
  } else {
    HeatMethodDistanceSolver heatSolver(*geometry);
    dist = heatSolver.computeDistance(geoDistSrcPt);
  }

  psMesh->addVertexScalarQuantity("Geodesic Distance", dist)
    ->setIsolinesEnabled(true)
    ->setColorMap("reds")
    ->setOnscreenColorbarEnabled(true)
    ->setEnabled(true);
    if (saveGeoDist) {
      psMesh->addVertexScalarQuantity("Geodesic Distance", dist)->exportColorbarToSVG("colorbar.svg");
    }
}


// Override GC's QEM implementation to have a target vertex count
void quadricErrorSimplify(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, size_t count) {
  double tol = 5.0;
  MutationManager mm(mesh, geo);
  if (mesh.nVertices() < count) return;
  
  auto toEigen = [](Vector3 v) -> Eigen::Vector3d {
    Eigen::Vector3d ret;
    ret << v.x, v.y, v.z;
    return ret;
  };
  auto fromEigen = [](Eigen::Vector3d v) -> Vector3 { return Vector3{v(0), v(1), v(2)}; };

  VertexData<Quadric> Q(mesh, Quadric());

  geo.requireFaceNormals();

  for (Face f : mesh.faces()) {
    Eigen::Vector3d n = toEigen(geo.faceNormals[f]);
    Eigen::Matrix3d M = n * n.transpose();
    for (Vertex v : f.adjacentVertices()) {
      Eigen::Vector3d q = toEigen(geo.inputVertexPositions[v]);
      double d = -n.dot(q);

      Q[v] += Quadric(M, d * n, d * d);
    }
  }

  using PotentialEdge = std::tuple<double, Edge>;

  auto cmp = [](const PotentialEdge& a, const PotentialEdge& b) -> bool { return std::get<0>(a) > std::get<0>(b); };

  std::priority_queue<PotentialEdge, std::vector<PotentialEdge>, decltype(cmp)> edgesToCheck(cmp);

  for (Edge e : mesh.edges()) {
    Quadric Qe = Q[e.halfedge().tailVertex()] + Q[e.halfedge().tipVertex()];
    Eigen::Vector3d q = Qe.optimalPoint();
    double cost = Qe.cost(q);
    edgesToCheck.push(std::make_tuple(cost, e));
  }

  while (!edgesToCheck.empty() && mesh.nVertices()>target) {
    PotentialEdge best = edgesToCheck.top();
    edgesToCheck.pop();

    // Stop when collapse becomes too expensive
    double cost = std::get<0>(best);
    if (cost > tol) break;
    if (!std::isfinite(cost)) continue; // numerical safety

    Edge e = std::get<1>(best);
    if (e.isDead()) continue; // edge no longer exists

    Vertex v1 = e.halfedge().tailVertex();
    Vertex v2 = e.halfedge().tipVertex();

    // Get edge quadric
    Quadric Qe(Q[v1], Q[v2]);
    Eigen::Vector3d q = Qe.optimalPoint();

    // If either vertex has been collapsed since the edge was pushed
    // onto the queue, the old cost was wrong. In that case, give up
    if (abs(cost - Qe.cost(q)) > 1e-8) continue;
    if (!q.array().isFinite().all()) continue; // numerical safety


    // if (geoDistSrcPt != SurfacePoint()) {
    //   for (Face f : e.adjacentFaces()) {
    //     if (geoDistSrcPt.face == f) {
          
    //     }
    //   }
    // }

    Vertex v = mm.collapseEdge(e, fromEigen(q));
    if (v == Vertex()) continue;
    Q[v] = Qe;

    for (Edge f : v.adjacentEdges()) {
      Quadric Qf(Q[f.halfedge().tailVertex()], Q[f.halfedge().tipVertex()]);
      Eigen::Vector3d q = Qf.optimalPoint();
      double cost = Qf.cost(q);
      edgesToCheck.push(std::make_tuple(cost, f));
    }
  }

  mesh.compress();
  geo.refreshQuantities();
  intrinsic.reset(new SignpostIntrinsicTriangulation(mesh, geo));
  intrinsic->flipToDelaunay();

  // if (geoDistSrcPt != SurfacePoint()) { geoDistSrcPt = geoDistSrcPt.nearestVertex(); }
  return;
}

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
extractVertNbhd(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, Vertex i) {
  std::vector<std::vector<size_t>> faces;
  std::vector<Vector3> positions;
  std::unordered_map<size_t, size_t> local;
  auto vInd = mesh.getVertexIndices();

  positions.push_back(geo.inputVertexPositions[i]);
  local[vInd[i]] = 0;

  size_t localIdx = 1;
  for (Vertex v : i.adjacentVertices()) {
    positions.push_back(geo.inputVertexPositions[v]);
    local[vInd[v]] = localIdx++;
  }

  for (Halfedge he : i.outgoingHalfedges()) {
    if (!he.isInterior()) {
      continue;
    }
    std::vector<size_t> faceList;
    for (Vertex v : he.face().adjacentVertices()) {
      faceList.push_back(local[vInd[v]]);
    }
    faces.push_back(faceList);
  }

  std::unique_ptr<ManifoldSurfaceMesh> newMesh;
  std::unique_ptr<VertexPositionGeometry> newGeo;

  std::tie(newMesh, newGeo) = makeManifoldSurfaceMeshAndGeometry(faces, positions);

  return std::make_tuple(std::move(newMesh), std::move(newGeo));
}

// Save current mesh data
void saveMesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, std::string name) {
  // Save obj data;
  std::filesystem::create_directories(outPath);
  writeSurfaceMesh(mesh, geo, fmt::format("{0}{1}_qed_iteration_{2}.obj", outPath, name, target));
}

void flattenAndVis() {
  auto psNbhd = polyscope::getSurfaceMesh("nbhd");
  VertexData<Vector2> param = parameterizeBFF(*n_mesh, *n_geometry);
  std::unordered_map<Vertex, size_t> vertMap;
}


// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
  
  psMesh->setSelectionMode(polyscope::MeshSelectionMode::Auto);

  ImGui::SliderInt("Target Vertex Count", &target, 5, mesh->nVertices());

  if (ImGui::Button("Run QEM Simplification")) {
    if (recordOutput) { saveMesh(*mesh, *geometry, psMesh->getName()); }
    if (polyscope::hasCurveNetwork("idt")) {
      polyscope::removeCurveNetwork("idt");
    }

    Vector3 oldSrcPos;
    if (geoDistSrcPt != SurfacePoint()) {
      oldSrcPos = geoDistSrcPt.interpolate(geometry->inputVertexPositions);
    }

    quadricErrorSimplify(*mesh, *geometry, static_cast<size_t>(target));
    psMesh = polyscope::registerSurfaceMesh(
      meshName,
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh)
    );
    
    if (polyscope::hasPointCloud("src") && geoDistSrcPt != SurfacePoint()) {
      polyscope::registerPointCloud("src", std::vector<Vector3>{geoDistSrcPt.interpolate(geometry->inputVertexPositions)})
      ->setPointRadius(0.005)
      ->setPointColor(glm::vec3{255, 128, 0})
      ->setEnabled(true);
    }
  }

  if (ImGui::Button("Compute Geodesic Distance via Heat Method")) {
    computeGeoDistFromSrc();
  }

  if (ImGui::Button("Render Intrinsic Triangulation")) {
    // if (polyscope::hasCurveNetwork("idt")) { return; }
    visualizeIntrinsicTriangulation();
  }

  

  if (ImGui::Button("Run ICE Simplification")) {
    size_t n = std::round(intrinsic->mesh.nVertices() * 0.90);
    intrinsicallyCoarsen(*intrinsic, 200);
  }

  if (ImGui::Button("Create Neighborhood")) {
    if (flattenTarget == SurfacePoint()) { return; }
    std::tie(n_mesh, n_geometry) = extractVertNbhd(*mesh, *geometry, flattenTarget.nearestVertex());
    polyscope::registerSurfaceMesh("nbhd", n_geometry->inputVertexPositions, n_mesh->getFaceVertexList());
  }

  if (ImGui::Button("Flatten and Remove Int. Vertex")) {
    if (flattenTarget == SurfacePoint()) { return; }
    auto tgt_vert = flattenTarget.nearestVertex();
    tgt_vert = intrinsic->intrinsicMesh->vertex(tgt_vert.getIndex());

    auto rslt = flattenVertexCETM(*intrinsic, tgt_vert);
    fmt::print("Result: {}\n", rslt);
    if (rslt) { 
      if(flatVertexRemoval(*intrinsic, tgt_vert)) {
        flattenTarget = SurfacePoint();
      }
    }

  }

  if (ImGui::Button("Select Source Vertex")) {
    selectPoint(geoDistSrcPt, "src");
  }

  if (ImGui::Button("Select Flatten Target Vertex")) {
    selectPoint(flattenTarget, "tgt");
  }

  ImGui::Checkbox("Save Geodesic Distance SVG", &saveGeoDist);
  ImGui::Checkbox("Use iDT for Geodesic Computation", &useIDT); 

  if (ImGui::Button("Write Intrinsic Mesh to .PLY")) {
    intrinsic->requireVertexGaussianCurvatures();
    intrinsic->requireEdgeCotanWeights();
    intrinsic->requireEdgeLengths();
    RichSurfaceMeshData richData(intrinsic->mesh);
    richData.outputFormat = happly::DataFormat::ASCII;
    // richData.addMeshConnectivity();
    richData.addMeshConnectivity();
    richData.addVertexProperty("vertex_gaussian_curvature", intrinsic->vertexGaussianCurvatures);
    richData.addEdgeProperty("edge_cotan_weights", intrinsic->edgeCotanWeights);
    richData.write("mesh.ply");
  }

  geometry->requireEdgeCotanWeights();
  ImGui::Text("Min. Edge Cotan. Weight: {}");
}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");
  args::Positional<bool> recordOutput(parser, "record output", "Record output to \'~/MeshAnalysis/QEDIterations/\'.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();  

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));
  intrinsic.reset(new SignpostIntrinsicTriangulation(*mesh, *geometry));
  intrinsic->flipToDelaunay();


  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(args::get(inputFilename)),
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  meshName = psMesh->getName();
  target = psMesh->nVertices();


  // Set vertex tangent spaces
  geometry->requireVertexTangentBasis();
  VertexData<Vector3> vBasisX(*mesh);
  VertexData<Vector3> vBasisY(*mesh);
  for (Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
    vBasisY[v] = geometry->vertexTangentBasis[v][1];
  }
  auto vField =
      geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry);
  psMesh->addVertexTangentVectorQuantity("VF", vField, vBasisX, vBasisY);


  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
